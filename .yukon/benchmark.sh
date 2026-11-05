#!/usr/bin/env bash

# Trusted benchmark orchestration. The score is the median of five sequential
# runs, published only when every run passes verification and the proof-size limit.
#
# usage: benchmark.sh [--track x86|arm] [--prove-only DIR | --verify-only DIR]
# No track flag selects the original x86 contract.
#
# By default this proves and verifies on one host. The arm track's CI splits the two: the ARM
# job runs `--prove-only DIR`, which leaves each run's seed, timing, and proof in DIR together
# with a description of the proving host, and an x86 job runs `--verify-only DIR` with the
# pinned verifier and publishes the score.

set -euo pipefail
umask 077

script_dir=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=.yukon/lib/sandbox.sh
source "$script_dir/lib/sandbox.sh"

mode=full
runs_dir=""
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --track | --prove-only | --verify-only)
      [[ "$#" -ge 2 ]] || { echo "$1 needs a value" >&2; exit 2; }
      case "$1" in
        --track) LABINIUS_TRACK=$2 ;;
        --prove-only) mode=prove; runs_dir=$2 ;;
        --verify-only) mode=verify; runs_dir=$2 ;;
      esac
      shift 2
      ;;
    *)
      echo "usage: benchmark.sh [--track x86|arm] [--prove-only DIR | --verify-only DIR]" >&2
      exit 2
      ;;
  esac
done
export LABINIUS_TRACK
track=$(labinius_track)
score_file=$(labinius_score_file)

# Recover the successful local build when no CI workspace was supplied.
if [[ -z "${RUNNER_TEMP:-}" && -z "${GITHUB_ACTIONS:-}" ]]; then
  rm -f "$script_dir/$score_file"
  state="$script_dir/local-workspace-$track"
  if [[ ! -f "$state" || -L "$state" ]]; then
    echo "Run yukon setup --track $track first." >&2
    exit 1
  fi
  IFS= read -r RUNNER_TEMP < "$state"
  export RUNNER_TEMP
  if [[ "$RUNNER_TEMP" != /var/tmp/labinius-local.* || ! -d "$RUNNER_TEMP" || -L "$RUNNER_TEMP" || ! -O "$RUNNER_TEMP" ]]; then
    echo "Local workspace is unavailable. Rerun yukon setup --track $track." >&2
    exit 1
  fi
  root=$(git -C "$script_dir" rev-parse --show-toplevel)
  if [[ $(cat "$RUNNER_TEMP/source-checkout") != "$root" ||
        $(git -C "$RUNNER_TEMP/trusted-checkout" rev-parse HEAD) != $(git -C "$root" rev-parse HEAD) ]]; then
    echo "Checkout changed. Rerun yukon setup --track $track." >&2
    exit 1
  fi
  unset TRUSTED_DIR CONTESTANT_DIR
fi

: "${RUNNER_TEMP:?RUNNER_TEMP is required}"
trusted="${TRUSTED_DIR:-$RUNNER_TEMP/labinius-trusted}"
contestant="${CONTESTANT_DIR:-$RUNNER_TEMP/labinius-contestant}"
require_runner_temp_path "$trusted"
require_runner_temp_path "$contestant"

score="$script_dir/$score_file"
staged="$script_dir/labinius-fast-$track-score.tmp"
checked_json="$script_dir/labinius-fast-$track-runs.tmp"
rm -f "$score" "$staged" "$checked_json"
max_bytes=$((64 * 1024 * 1024))

# A trusted binary from setup, checked against the hash recorded beside it.
trusted_binary() {
  local path="$trusted/$1"
  if [[ ! -f "$path" || -L "$path" ]]; then
    echo "missing regular trusted binary $1" >&2
    return 1
  fi
  if [[ "$(sha256_file "$path")" != "$(awk '{print $1}' "$path.sha256")" ]]; then
    echo "trusted binary $1 changed since setup" >&2
    return 1
  fi
  printf '%s\n' "$path"
}

# Five sequential fresh-seed runs. Each run's proof, seed, and time go to $1, with the
# proving host's description; nothing here verifies.
prove_runs() {
  local out=$1 writer prover trial seed prove_ms proof_bytes
  local workdir table_dir table
  command -v bwrap >/dev/null || {
    echo "bubblewrap (bwrap) is required to run the untrusted prover" >&2
    return 1
  }
  prover="$contestant/solution/proof_generate/target/release/labinius-fast"
  if [[ ! -f "$prover" || -L "$prover" ]]; then
    echo "missing regular contestant prover binary" >&2
    return 1
  fi
  if [[ $track == arm ]]; then
    writer=("$(trusted_binary labinius-table-writer)")
  else
    writer=("$(trusted_binary labinius-fast-verify)" --write-table)
  fi

  local cpu_model avx512f
  cpu_model=$(awk -F: '/^model name/ { gsub(/^ +/, "", $2); print $2; exit }' /proc/cpuinfo)
  if [[ -z "$cpu_model" ]]; then
    # aarch64 /proc/cpuinfo has no model name line.
    cpu_model=$(lscpu 2>/dev/null | awk -F: '/^Model name/ { gsub(/^ +/, "", $2); print $2; exit }')
    cpu_model="${cpu_model:-$(uname -m)}"
  fi
  if grep -qw avx512f /proc/cpuinfo; then avx512f=true; else avx512f=false; fi
  python3 - "$out/host.json" "$cpu_model" "$(nproc)" "$avx512f" \
    "$(awk '/^MemTotal:/ { print $2; exit }' /proc/meminfo)" <<'PY'
import json, sys
path, cpu, nproc, avx512f, mem = sys.argv[1:6]
with open(path, "w", encoding="utf-8") as handle:
    json.dump({"cpuModel": cpu, "nproc": int(nproc), "avx512f": avx512f == "true",
               "memTotalKb": int(mem)}, handle)
PY
  printf '[]\n' > "$out/runs.json"

  workdir=$(mktemp -d "$RUNNER_TEMP/labinius-fast-proof.XXXXXX")
  table_dir=$(mktemp -d "$RUNNER_TEMP/labinius-fast-table.XXXXXX")
  table="$table_dir/witness-table.bin"
  for trial in 1 2 3 4 5; do
    # A fresh seed keeps a prover from reusing a cached proof. The witness table
    # is written before the timer; verification runs after it.
    seed=$(openssl rand -hex 32)
    rm -f "$table" "$workdir/proof.bin"
    env -i PATH=/usr/bin:/bin LABINIUS_SEED="$seed" "${writer[@]}" "$table"
    chmod 444 "$table"

    prove_ms=$(python3 - "$script_dir" "$prover" "$seed" "$workdir" "$table" <<'PY'
import subprocess
import sys
import time

script_dir, prover, seed, workdir, table = sys.argv[1:6]
cmd = [
    "/usr/bin/bash",
    f"{script_dir}/run-prover.sh",
    prover,
    workdir,
    seed,
    "/work/proof.bin",
    table,
]
start = time.perf_counter()
# Python's timeout wait polls at up to 50 ms intervals on these Linux hosts.
# Keep the hard deadline in GNU timeout and wait for process exit directly.
subprocess.run(["/usr/bin/timeout", "--signal=KILL", "900s", *cmd], check=True)
print((time.perf_counter() - start) * 1000)
PY
)

    proof_bytes=$(regular_file_size "$workdir/proof.bin")
    if (( proof_bytes <= 0 || proof_bytes > max_bytes )); then
      echo "empty or oversized proof on run ${trial}" >&2
      return 1
    fi
    mv "$workdir/proof.bin" "$out/proof-$trial.bin"
    python3 - "$out/runs.json" "$trial" "$seed" "$prove_ms" "$proof_bytes" <<'PY'
import json, sys
path, run, seed, prove_ms, proof_bytes = sys.argv[1:6]
with open(path, encoding="utf-8") as handle:
    runs = json.load(handle)
runs.append({"run": int(run), "seed": seed, "provingTimeMs": float(prove_ms),
             "proofBytes": int(proof_bytes)})
with open(path, "w", encoding="utf-8") as handle:
    json.dump(runs, handle)
PY
  done
  rm -rf "$workdir" "$table_dir"
}

# Check every proof in $1 against its seed with the trusted verifier, then publish the score
# from the recorded timings and the proving host's description.
verify_runs() {
  local dir=$1 verifier run seed proof proof_bytes verify_out
  verifier=$(trusted_binary labinius-fast-verify)
  python3 - "$dir/runs.json" "$dir/host.json" <<'PY'
import json, math, re, sys
runs = json.load(open(sys.argv[1], encoding="utf-8"))
host = json.load(open(sys.argv[2], encoding="utf-8"))
if not isinstance(runs, list) or [r.get("run") for r in runs] != [1, 2, 3, 4, 5]:
    raise SystemExit("runs.json does not describe runs 1 to 5 in order")
for r in runs:
    if not isinstance(r.get("seed"), str) or not re.fullmatch(r"[0-9a-f]{64}", r["seed"]):
        raise SystemExit("runs.json has an invalid seed")
    t = r.get("provingTimeMs")
    if isinstance(t, bool) or not isinstance(t, (int, float)) or not math.isfinite(t) or t < 0:
        raise SystemExit("runs.json has an invalid proving time")
if len({r["seed"] for r in runs}) != len(runs):
    raise SystemExit("runs.json reuses a seed")
if set(host) != {"cpuModel", "nproc", "avx512f", "memTotalKb"}:
    raise SystemExit("host.json has unexpected fields")
PY
  printf '[]\n' > "$checked_json"
  for run in 1 2 3 4 5; do
    seed=$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))[int(sys.argv[2]) - 1]["seed"])' \
      "$dir/runs.json" "$run")
    proof="$dir/proof-$run.bin"
    proof_bytes=$(regular_file_size "$proof")
    if (( proof_bytes <= 0 || proof_bytes > max_bytes )); then
      echo "empty or oversized proof on run ${run}" >&2
      return 1
    fi

    verify_out=$(env -i PATH=/usr/bin:/bin LABINIUS_SEED="$seed" python3 - "$verifier" "$proof" <<'PY'
import subprocess
import sys

checked = subprocess.run(
    [sys.argv[1], "--verify", sys.argv[2]],
    timeout=900,
    capture_output=True,
    text=True,
)
if checked.returncode != 0:
    sys.stderr.write(checked.stderr)
    sys.stderr.write(checked.stdout)
    raise SystemExit(checked.returncode)
sys.stdout.write(checked.stdout)
PY
)

    python3 - "$checked_json" "$dir/runs.json" "$run" "$proof_bytes" "$verify_out" <<'PY'
import json, sys
path, recorded, run, proof_bytes, verifier = sys.argv[1:6]
with open(path, encoding="utf-8") as handle:
    checked = json.load(handle)
with open(recorded, encoding="utf-8") as handle:
    record = json.load(handle)[int(run) - 1]
if record["proofBytes"] != int(proof_bytes):
    raise SystemExit(f"run {run}: proof size differs from the recorded one")
checked.append({
    "provingTimeMs": float(record["provingTimeMs"]),
    "proofBytes": int(proof_bytes),
    "verifier": json.loads(verifier),
})
with open(path, "w", encoding="utf-8") as handle:
    json.dump(checked, handle)
PY
  done

  python3 - "$script_dir/lib/score.py" "$checked_json" "$dir/host.json" "$staged" <<'PY'
import json, subprocess, sys
score_py, runs, host_path, out = sys.argv[1:5]
host = json.load(open(host_path, encoding="utf-8"))
subprocess.run([
    sys.executable, score_py, "--runs", runs,
    "--cpu-model", str(host["cpuModel"]), "--nproc", str(host["nproc"]),
    "--avx512f", "true" if host["avx512f"] is True else "false",
    "--mem-total-kb", str(host["memTotalKb"]), "--out", out,
], check=True)
PY
  mv "$staged" "$score"
}

cleanup() {
  rm -f "$staged" "$checked_json"
  if [[ -n "${scratch:-}" ]]; then
    rm -rf "$scratch"
  fi
}
trap cleanup EXIT

case "$mode" in
  full)
    scratch=$(mktemp -d "$RUNNER_TEMP/labinius-fast-runs.XXXXXX")
    prove_runs "$scratch"
    verify_runs "$scratch"
    ;;
  prove)
    require_runner_temp_path "$runs_dir"
    rm -rf "$runs_dir"
    mkdir -p "$runs_dir"
    prove_runs "$runs_dir"
    ;;
  verify)
    if [[ ! -d "$runs_dir" || -L "$runs_dir" ]]; then
      echo "proof directory is not a directory: $runs_dir" >&2
      exit 1
    fi
    verify_runs "$runs_dir"
    ;;
esac
