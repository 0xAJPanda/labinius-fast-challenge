#!/bin/sh
# Add the frozen execution profile without replacing the harness compiler flags.
set -eu
pgo_crate=false
for pgo_arg do
    if [ "$pgo_arg" = "--crate-name" ]; then pgo_crate=true; fi
done
if "$pgo_crate"; then
    pgo_root=$(CDPATH= cd "$(dirname "$0")/.." && pwd)
    exec "$@" -C "profile-use=$pgo_root/pgo.profdata"
fi
exec "$@"
