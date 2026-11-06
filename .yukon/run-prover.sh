#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=.yukon/lib/sandbox.sh
source "$script_dir/lib/sandbox.sh"

usage() {
  cat >&2 <<'EOF'
usage: run-prover.sh PROVER WORKDIR SEED SANDBOX_PROOF_PATH TABLE
         One sandboxed prover run, as benchmark.sh invokes it.
       run-prover.sh --track arm|x86 [--runs N] [--modulus Q] [--no-sandbox]
         Build the selected track's proof_generate crate and a checker from verifier/proof_verify,
         then prove and verify N fresh-seed instances (default 1). The track defaults to x86.
         The prover runs in the bubblewrap sandbox on Linux when bwrap is
         installed, and directly otherwise. LABINIUS_LOCAL_WORK sets the scratch directory.
EOF
  exit 2
}

if [[ "$#" -eq 5 && "$1" != --* ]]; then
  run_contestant_prover "$1" "$2" "$3" "$4" "$5"
  exit
fi

track=""
runs=1
modulus=""
sandbox=auto
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --track) [[ "$#" -ge 2 ]] || usage; track=$2; shift 2 ;;
    --runs) [[ "$#" -ge 2 ]] || usage; runs=$2; shift 2 ;;
    --modulus) [[ "$#" -ge 2 ]] || usage; modulus=$2; shift 2 ;;
    --no-sandbox) sandbox=off; shift ;;
    *) usage ;;
  esac
done
track=$(LABINIUS_TRACK=$track labinius_track)
solution_path=$(LABINIUS_TRACK=$track labinius_solution_path)
[[ "$runs" =~ ^[1-9][0-9]*$ ]] || usage
if [[ "$track" == x86 && "$(uname -m)" != x86_64 ]]; then
  echo "the x86 track builds AVX-512 kernels and needs an x86-64 host; use --track arm" >&2
  exit 1
fi
if [[ "$sandbox" == auto ]]; then
  if [[ "$(uname -s)" == Linux ]] && command -v bwrap >/dev/null 2>&1; then
    sandbox=on
  else
    sandbox=off
  fi
fi
if [[ "$sandbox" == on && -n "$modulus" ]]; then
  echo "--modulus is not passed through the sandbox; add --no-sandbox" >&2
  exit 2
fi

root=$(cd "$script_dir/.." && pwd)
tmp=${TMPDIR:-/tmp}
work=${LABINIUS_LOCAL_WORK:-${tmp%/}/labinius-local-$track}
mkdir -p "$work/input" "$work/out"

now_ms() { perl -MTime::HiRes=time -e 'printf "%d\n", time * 1000'; }

echo "==> building $solution_path/proof_generate" >&2
prover_dir="$root/$solution_path/proof_generate"
(cd "$prover_dir" && CARGO_TARGET_DIR="$prover_dir/target" cargo build --release --locked)
prover="$prover_dir/target/release/labinius-fast"

# The pinned verifier/labinius only builds on x86-64; the arm track checks against its portable
# port of the same library.
if [[ "$track" == arm ]]; then
  library="$root/solution-arm/labinius/crates/pcs"
else
  library="$root/verifier/labinius/crates/pcs"
fi
echo "==> building the checker against ${library#"$root"/}" >&2
checker_src="$work/proof_verify"
rm -rf "$checker_src"
mkdir -p "$checker_src"
cp -R "$root/verifier/proof_verify/." "$checker_src/"
python3 - "$checker_src/Cargo.toml" "$library" <<'PY'
import sys
from pathlib import Path
manifest, library = Path(sys.argv[1]), sys.argv[2]
old = 'labinius = { path = "../labinius/crates/pcs" }'
text = manifest.read_text()
if text.count(old) != 1:
    sys.exit(f"{manifest}: expected one '{old}' dependency line")
manifest.write_text(text.replace(old, f'labinius = {{ path = "{library}" }}'))
PY
(cd "$checker_src" && CARGO_TARGET_DIR="$work/checker-target" cargo build --release --locked)
checker="$work/checker-target/release/labinius-fast-verify"

table="$work/input/witness-table.bin"
proof="$work/out/proof.bin"
prover_args=()
if [[ -n "$modulus" ]]; then
  prover_args=(--modulus "$modulus")
fi
echo "==> track=$track sandbox=$sandbox runs=$runs" >&2
times=()
for run in $(seq 1 "$runs"); do
  seed=$(openssl rand -hex 32)
  rm -f "$table" "$proof"
  LABINIUS_SEED=$seed "$checker" --write-table "$table"
  chmod 0444 "$table"
  start=$(now_ms)
  if [[ "$sandbox" == on ]]; then
    run_contestant_prover "$prover" "$work/out" "$seed" /work/proof.bin "$table"
  else
    LABINIUS_SEED=$seed "$prover" --prove "$proof" --table "$table" ${prover_args[@]+"${prover_args[@]}"}
  fi
  elapsed=$(($(now_ms) - start))
  size=$(wc -c < "$proof" | tr -d ' ')
  verdict=$(LABINIUS_SEED=$seed "$checker" --verify "$proof")
  echo "run $run/$runs: seed=$seed prove_ms=$elapsed proof_bytes=$size verified=$verdict"
  times+=("$elapsed")
done
rm -f "$table"
echo "proof of the last run: $proof"
if ((runs > 1)); then
  median=$(printf '%s\n' "${times[@]}" | sort -n | awk '{a[NR]=$1} END {print (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2}')
  echo "median prove_ms=$median over $runs runs (${times[*]})"
fi
