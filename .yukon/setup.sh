#!/usr/bin/env bash
set -euo pipefail
umask 077

# usage: setup.sh --track x86|arm [--verifier-only]
#
# Every proof is checked by the pinned verifier, which builds and runs on x86-64 with AVX-512
# only. The x86 track builds it here next to the contestant prover. On an ARM host the arm track
# builds the contestant prover and the BLAKE3-only table writer instead; CI then checks its
# proofs in a separate x86 job, whose setup is `--track arm --verifier-only`. A local arm setup,
# which has no x86 host to hand proofs to, also builds a development verifier against the
# portable library in the trusted checkout's solution/arm/labinius.

script_dir=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=.yukon/lib/sandbox.sh
source "$script_dir/lib/sandbox.sh"

verifier_only=false
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --verifier-only) verifier_only=true; shift ;;
    --track)
      [[ "$#" -ge 2 ]] || { echo "--track needs a value: x86 or arm" >&2; exit 2; }
      LABINIUS_TRACK=$2
      shift 2
      ;;
    *) echo "usage: setup.sh --track x86|arm [--verifier-only]" >&2; exit 2 ;;
  esac
done
export LABINIUS_TRACK
track=$(labinius_track)
if $verifier_only; then
  host_arch=x86_64
else
  case "$track" in
    x86) host_arch=x86_64 ;;
    arm) host_arch=aarch64 ;;
  esac
fi

# GitHub supplies these paths. Prepare them here for local CLI use.
local_setup=false
if [[ -z "${RUNNER_TEMP:-}" && -z "${GITHUB_ACTIONS:-}" ]]; then
  local_setup=true
  rm -f "$script_dir/local-workspace-$track" "$script_dir/labinius-fast-$track-score.json"
  if [[ $(uname -s) != Linux || $(uname -m) != "$host_arch" ]]; then
    echo "Local setup of the $track track requires Ubuntu 24.04 on $host_arch, not $(uname -s) $(uname -m)." >&2
    exit 1
  fi
  source /etc/os-release
  if [[ ${ID:-} != ubuntu || ${VERSION_ID:-} != 24.04 ]]; then
    echo 'Local setup requires Ubuntu 24.04.' >&2
    exit 1
  fi
  if [[ $host_arch == x86_64 ]]; then
    grep -qw avx512f /proc/cpuinfo || { echo 'AVX-512F is required.' >&2; exit 1; }
  fi

  packages=(build-essential gcc-13 bubblewrap python3 openssl libssl-dev pkg-config curl ca-certificates git)
  missing=()
  for package in "${packages[@]}"; do
    if [[ $(dpkg-query -W -f='${Status}' "$package" 2>/dev/null || true) != 'install ok installed' ]]; then
      missing+=("$package")
    fi
  done
  if (( ${#missing[@]} )); then
    sudo apt-get update
    sudo apt-get install -y "${missing[@]}"
  fi

  if [[ ! -e /opt/cargo && ! -L /opt/cargo && ! -e /opt/rustup && ! -L /opt/rustup ]]; then
    installer=$(mktemp)
    trap 'rm -f "$installer"' EXIT
    curl --proto '=https' --tlsv1.2 -fsS https://sh.rustup.rs -o "$installer"
    sudo env CARGO_HOME=/opt/cargo RUSTUP_HOME=/opt/rustup sh "$installer" \
      -y --profile minimal --default-toolchain stable --no-modify-path
    rm -f "$installer"
    trap - EXIT
  fi
  env CARGO_HOME=/opt/cargo RUSTUP_HOME=/opt/rustup /opt/cargo/bin/rustc --version

  export UNTRUSTED_CHECKOUT="$(git -C "$script_dir" rev-parse --show-toplevel)"
  revision=$(git -C "$UNTRUSTED_CHECKOUT" rev-parse HEAD)
  export RUNNER_TEMP="$(mktemp -d /var/tmp/labinius-local.XXXXXX)"
  export TRUSTED_CHECKOUT="$RUNNER_TEMP/trusted-checkout"
  printf 'Local workspace: %s\n' "$RUNNER_TEMP"
  git clone --no-hardlinks --no-checkout "$UNTRUSTED_CHECKOUT" "$TRUSTED_CHECKOUT"
  git -C "$TRUSTED_CHECKOUT" checkout --detach "$revision"
  git -C "$TRUSTED_CHECKOUT" submodule update --init --recursive
  printf '%s\n' "$UNTRUSTED_CHECKOUT" > "$RUNNER_TEMP/source-checkout"
  unset TRUSTED_DIR TRUSTED_CARGO_HOME TRUSTED_TARGET_DIR CONTESTANT_DIR CONTESTANT_CARGO_HOME
fi

: "${RUNNER_TEMP:?RUNNER_TEMP is required}"
: "${TRUSTED_CHECKOUT:?TRUSTED_CHECKOUT is required}"
: "${UNTRUSTED_CHECKOUT:?UNTRUSTED_CHECKOUT is required}"

if [[ $(uname -m) != "$host_arch" ]]; then
  echo "this setup needs an $host_arch host, not $(uname -m)" >&2
  exit 1
fi

trusted="${TRUSTED_DIR:-$RUNNER_TEMP/labinius-trusted}"
trusted_cargo="${TRUSTED_CARGO_HOME:-$RUNNER_TEMP/labinius-trusted-cargo}"
trusted_target="${TRUSTED_TARGET_DIR:-$RUNNER_TEMP/labinius-trusted-target}"
contestant="${CONTESTANT_DIR:-$RUNNER_TEMP/labinius-contestant}"
contestant_cargo="${CONTESTANT_CARGO_HOME:-$RUNNER_TEMP/labinius-contestant-cargo}"

require_runner_temp_path "$trusted"
require_runner_temp_path "$trusted_cargo"
require_runner_temp_path "$trusted_target"
require_runner_temp_path "$contestant"
require_runner_temp_path "$contestant_cargo"
require_host_command /opt/cargo/bin/cargo
require_host_command /opt/cargo/bin/rustc
require_host_command /usr/bin/gcc-13

# Which trusted binaries this host needs. Each is installed with its hash beside it.
build_pinned_verifier=false
build_portable_verifier=false
build_table_writer=false
if $verifier_only || [[ $track == x86 ]]; then
  build_pinned_verifier=true
else
  build_table_writer=true
  build_portable_verifier=$local_setup
fi

if $build_pinned_verifier && [[ ! -f "$TRUSTED_CHECKOUT/verifier/labinius/crates/pcs/Cargo.toml" ]]; then
  echo "trusted checkout is missing its pinned verifier submodule" >&2
  exit 1
fi
if ! $verifier_only; then
  command -v bwrap >/dev/null || {
    echo "bubblewrap (bwrap) is required to build the untrusted solution" >&2
    exit 1
  }
  if [[ ! -d "$UNTRUSTED_CHECKOUT/solution" ]]; then
    echo "untrusted checkout has no solution directory" >&2
    exit 1
  fi
fi

mkdir -p "$trusted" "$trusted_cargo" "$trusted_target"
chmod 700 "$trusted" "$trusted_cargo" "$trusted_target"

# Build only protected trusted source outside the contestant sandbox. Its Cargo
# home and target are disposable and never share contestant state or host tools.
trusted_build() {
  local manifest=$1 bin=$2
  env -i \
    PATH=/opt/cargo/bin:/usr/bin:/bin \
    HOME="$trusted_cargo" \
    RUSTUP_HOME=/opt/rustup \
    CARGO_HOME="$trusted_cargo" \
    CARGO_TARGET_DIR="$trusted_target" \
    CC=/usr/bin/gcc-13 \
    CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/gcc-13 \
    CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/gcc-13 \
    RUSTFLAGS="-C target-cpu=native" \
    /opt/cargo/bin/cargo build \
      --locked \
      --release \
      --bin "$bin" \
      --manifest-path "$manifest"
  install -m 755 "$trusted_target/release/$bin" "$trusted/$bin"
  sha256_file "$trusted/$bin" > "$trusted/$bin.sha256"
  chmod 555 "$trusted/$bin"
  chmod 444 "$trusted/$bin.sha256"
}

if $build_pinned_verifier; then
  trusted_build "$TRUSTED_CHECKOUT/verifier/proof_verify/Cargo.toml" labinius-fast-verify
fi
if $build_table_writer; then
  trusted_build "$TRUSTED_CHECKOUT/verifier/table_writer/Cargo.toml" labinius-table-writer
fi
if $build_portable_verifier; then
  # Local development only: the same checker, linked against the trusted checkout's
  # portable port of the library because the pinned one does not build on ARM.
  verifier_src="$RUNNER_TEMP/labinius-trusted-verifier-src"
  portable_pcs="$TRUSTED_CHECKOUT/solution/arm/labinius/crates/pcs"
  require_runner_temp_path "$verifier_src"
  if [[ ! -f "$portable_pcs/Cargo.toml" ]]; then
    echo "trusted checkout is missing the portable labinius library" >&2
    exit 1
  fi
  rm -rf "$verifier_src"
  mkdir -p "$verifier_src"
  cp -R "$TRUSTED_CHECKOUT/verifier/proof_verify/." "$verifier_src/"
  rm -rf "$verifier_src/target"
  python3 - "$verifier_src/Cargo.toml" "$portable_pcs" <<'PY'
import sys
from pathlib import Path

manifest, pcs = Path(sys.argv[1]), sys.argv[2]
text = manifest.read_text(encoding="utf-8")
pinned = 'labinius = { path = "../labinius/crates/pcs" }'
if text.count(pinned) != 1:
    raise SystemExit("verifier manifest no longer names the pinned labinius path")
manifest.write_text(text.replace(pinned, f'labinius = {{ path = "{pcs}" }}'), encoding="utf-8")
PY
  trusted_build "$verifier_src/Cargo.toml" labinius-fast-verify
fi
chmod 555 "$trusted"

if ! $verifier_only; then
  # Copy data only: reject symlinks, special files, git control data, and paths
  # that cannot remain beneath solution/. This private RUNNER_TEMP copy is the
  # contestant's only writable source tree for native build-system outputs.
  python3 "$script_dir/stage-solution.py" \
    "$UNTRUSTED_CHECKOUT/solution" \
    "$contestant/solution"

  prover_dir="$contestant/solution/$track/proof_generate"
  if [[ ! -f "$prover_dir/Cargo.toml" || -L "$prover_dir/Cargo.toml" ]]; then
    echo "solution has no $track track prover at solution/$track/proof_generate" >&2
    exit 1
  fi
  mkdir -p "$contestant_cargo" "$prover_dir/target"
  chmod 700 "$contestant_cargo"

  run_contestant_build \
    "$contestant/solution" \
    "$prover_dir/target" \
    "$contestant_cargo" \
    fetch \
    /usr/bin/bash -eu -c '
      manifest=./Cargo.toml
      if [[ ! -x /workspace || ! -x /workspace/solution || ! -r "$manifest" ]]; then
        echo "staged contestant manifest is not readable inside the sandbox" >&2
        exit 1
      fi
      exec /opt/cargo/bin/cargo fetch --locked --manifest-path "$manifest"
    '

  run_contestant_build \
    "$contestant/solution" \
    "$prover_dir/target" \
    "$contestant_cargo" \
    offline \
    /usr/bin/bash -eu -c '
      manifest=./Cargo.toml
      if [[ ! -x /workspace || ! -x /workspace/solution || ! -r "$manifest" ]]; then
        echo "staged contestant manifest is not readable inside the sandbox" >&2
        exit 1
      fi
      exec /opt/cargo/bin/cargo build \
        --locked --offline --release --bin labinius-fast --manifest-path "$manifest"
    '

  prover="$prover_dir/target/release/labinius-fast"
  if [[ ! -f "$prover" || -L "$prover" ]]; then
    echo "contestant build did not produce a regular prover binary" >&2
    exit 1
  fi
fi

# Publish state only after every build succeeds. Each setup uses a fresh workspace.
if $local_setup; then
  state=$(mktemp "$script_dir/local-workspace-$track.XXXXXX")
  printf '%s\n' "$RUNNER_TEMP" > "$state"
  mv -f "$state" "$script_dir/local-workspace-$track"
  echo "Setup complete. Run yukon run --track $track; rerun yukon setup --track $track after solution edits."
fi
