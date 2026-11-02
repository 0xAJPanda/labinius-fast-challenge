#!/usr/bin/env bash
# Install what setup.sh expects on a hosted Ubuntu 24.04 runner: GCC 13, bubblewrap, and a
# read-only Rust toolchain at /opt/cargo and /opt/rustup. The dedicated x86 runner image ships
# these already.
set -euo pipefail

sudo apt-get update
sudo apt-get install -y build-essential gcc-13 bubblewrap python3 openssl libssl-dev \
  pkg-config curl ca-certificates git
if [[ ! -e /opt/cargo && ! -e /opt/rustup ]]; then
  installer=$(mktemp)
  curl --proto '=https' --tlsv1.2 -fsS https://sh.rustup.rs -o "$installer"
  sudo env CARGO_HOME=/opt/cargo RUSTUP_HOME=/opt/rustup sh "$installer" \
    -y --profile minimal --default-toolchain stable --no-modify-path
  rm -f "$installer"
fi
env CARGO_HOME=/opt/cargo RUSTUP_HOME=/opt/rustup /opt/cargo/bin/rustc --version

# Ubuntu 24.04 denies unprivileged user namespaces to bwrap by default.
if [[ -e /proc/sys/kernel/apparmor_restrict_unprivileged_userns ]]; then
  sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
fi
bwrap --ro-bind / / --unshare-pid --unshare-net -- true
