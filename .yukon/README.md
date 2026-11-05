# Dedicated runner security boundary

## Compatible x86 conversion

The x86 track preserves the standalone execution contract: `solution/`, no-flag
`setup.sh` and `benchmark.sh`, `.yukon/labinius-fast-score.json`, and the
`benchmark.yml` entrypoint. ARM lives in the sibling `solution-arm/` directory.
Yukon's in-place conversion (#1348) can therefore keep the x86 board ID and
history without rewriting submission archives.

Executor publication is a separate requirement. Commit the harness/executor
changes first, then pin `benchmark.yml`, `dispatch-x86.yml`, and `dispatch-arm.yml`
to that exact commit. Publish the referenced commit without squashing it away.
The runner-group owners must admit the new immutable executor references in both
groups before either dispatcher can use them. Retain the previous admission
entries while existing runs still need them. A passing local test does not prove
GitHub dispatch or runner admission.

Before production conversion, exercise the transition on dev with a standalone
x86 board, representative historical x86 submissions, and an ARM import. Confirm
that both workflows complete all five runs, verify every proof, and that Yukon
ingests the resulting scores. Recheck the current production x86 solution before
merging so a newer promotion is not overwritten. Pause production intake and
drain jobs before merging and converting; verify both tracks before reopening.

Local regression checks (Python 3 with PyYAML, Bash 4.3 or newer on `PATH`):

```bash
python3 -m unittest discover -s .github/tests -v
python3 -m unittest discover -s .yukon/tests -v
```

The track compatibility tests use a stub verifier and capture sandbox arguments.
They test path selection, score publication, and failure cleanup; they do not
replace the Linux sandbox smoke tests or real proof verification.

## Runner admission and isolation

The benchmark job targets organization runner group
`labinius-intel-8488c-avx512-c7i` and compatibility label
`labinius-avx512-c7i`. The group is an external admission boundary:
organization settings must keep it closed until review, then restrict it to the
`x86` track's workflow, `.github/workflows/benchmark-x86.yml`.
Repository labels alone do not prevent another workflow from requesting a
self-hosted runner.

The `arm` track's workflow, `benchmark-arm.yml`, has two jobs. The `prove` job
targets the admission-controlled group `labinius-graviton4-neoverse-v2-sve2-c8g` (labels
`ARM64`, `sve2`). That group's restriction must admit the `arm` workflow. The image ships the same
toolchain paths, `gcc-13`, and bubblewrap as the x86 host, and the job checks
the CPU with `/usr/local/bin/validate-graviton4-cpu`. The pinned `verifier/labinius` builds only on
x86-64, so this job writes input tables with the trusted
`verifier/table_writer`. That writer depends only on BLAKE3 and produces the
same bytes as `labinius-fast-verify --write-table`. The job proves five fresh seeds with
`benchmark.sh --prove-only` and uploads the proofs, `runs.json` (seeds, timings,
and sizes recorded by the harness), and `host.json` as a one-day artifact. Contestant code
writes only the proof files. The `verify` job runs on an x86 Blacksmith runner with AVX-512
(`blacksmith-8vcpu-ubuntu-2404`). It installs the toolchain with `.yukon/provision.sh`, then builds the pinned verifier with
`setup.sh --verifier-only`, validates the run records, and verifies every proof with
`benchmark.sh --verify-only`. It then writes the score from the ARM host's description and
deletes the proof artifact, which needs `actions: write`.

Every job checks out its trusted harness and verifier at the immutable executor commit
(`job.workflow_sha` in `job.workflow_repository`) and runs only that checkout's scripts. Yukon admits submissions that change
only the track's editable path, so that harness is the benchmark's own. The solution
repository and immutable commit are checked out separately: the optional
`solution_repository`/`solution_ref` inputs, defaulting to the dispatched commit.
The `arm` verify job never checks out the solution. Only regular files and real
directories beneath the selected `solution/` or `solution-arm/` are copied into disposable runner state at `solution/`;
symlinks, special files, `.git` control data, and path escapes are rejected. No
workflow or Yukon script from the solution checkout is executed.

The host must provide the preinstalled, read-only toolchain paths
`/opt/cargo/bin`, `/opt/rustup`, `/usr/bin/gcc-13`, plus `bwrap`, Python 3,
OpenSSL, and standard GNU userland. The trusted verifier and contestant use
separate disposable Cargo homes and target directories. Native build systems may
write only within the validated private staged solution under `RUNNER_TEMP`; the
original solution checkout, trusted checkout, and host toolchain remain outside
that writable tree. A restricted `cargo fetch --locked` phase has networking;
the subsequent offline compilation and prover execution have no network. The
prover sandbox sees only its binary,
runtime libraries, `/proc`, `/dev`, a fresh `/tmp`, its proof output
directory, and a read-only mount of the witness table at
`/input/witness-table.bin`.

## Score and witness input

The leaderboard score is the median prover-process time, in milliseconds, across five sequential runs. Each run draws a fresh 32-byte seed. Before the timer starts, the trusted verifier writes that seed's witness table with `--write-table`. On `arm`, `labinius-table-writer` writes the same bytes. The timed interval then covers prover startup, loading of that file, and the rest of the prover process. Table generation and trusted verification stay outside the timer.

A score file is written only when every run exits successfully, the proof is a regular file of at most 64 MiB, and the trusted verifier accepts the proof for that run's seed. Verification rebuilds the table from the seed with `trace` and does not read the file given to the prover. A failed run leaves the score unpublished: `.yukon/labinius-fast-score.json` for x86, `.yukon/labinius-fast-arm-score.json` for ARM.

`score` and `metrics.provingTimeMs` are that median. `metrics.proofBytes` is the proof size of the selected run. Runs are ordered by proving time and then by earlier execution index, and the middle position of that order is the selected run. `runs` lists the five executions in order, each with `run`, `provingTimeMs`, and `proofBytes`. The other `metrics` fields keep their previous meaning and come from the selected run's verifier output plus the host description.

The scored prover is invoked inside the sandbox as:

```text
labinius-fast --prove /work/proof.bin --table /input/witness-table.bin
```

`LABINIUS_SEED` is the hex seed for that run. `/input` is a read-only directory containing only `witness-table.bin`. The file is the BLAKE3 XOF keyed by the seed, updated with the context `labinius-fast/binary-table/v1`, and filled to `(1 << 31) / 8` bytes. Those bytes are the little-endian encoding of the 2^24 GHASH-field elements produced by `trace(seed)`. The baseline prover loads this file and does not compute the XOF during the timed run. Optional `--modulus Q` may follow `--table`.

The single-run numbers in the disposable-clone section below predate median-of-five scoring.

Bubblewrap reduces what the contestant process can see and write on the host;
it is not a VM reset or a claim of complete isolation. The runner image,
metadata firewall, restrictive bubblewrap/AppArmor policy, non-sudo runner
account, runner-group admission, and post-job disposal/cleanup remain required
layers. Do not place application keys or long-lived secrets in the host
repository or job environment.

On a disposable runner image clone, run the regression test:

```bash
export RUNNER_TEMP
RUNNER_TEMP=$(mktemp -d /var/tmp/labinius-runner-test.XXXXXX)
bash .yukon/tests/test-sandbox.sh
rm -rf "$RUNNER_TEMP"
```

It requires the same toolchain paths and packages as production. It proves that
solution symlinks are rejected, an inherited marker variable is cleared, an
arbitrary host path is not visible, proof FIFOs/symlinks are rejected,
`/opt/cargo/bin` is not writable, a native build can write only in the private
staged solution, and a legitimate dependency-free Rust crate still builds.

Because the current base includes the generic-modulus prover/verifier change,
run a genuine full benchmark smoke test rather than relying on an older proof:

```bash
export RUNNER_TEMP TRUSTED_CHECKOUT UNTRUSTED_CHECKOUT
export TRUSTED_DIR TRUSTED_CARGO_HOME TRUSTED_TARGET_DIR
export CONTESTANT_DIR CONTESTANT_CARGO_HOME
RUNNER_TEMP=$(mktemp -d /var/tmp/labinius-full-test.XXXXXX)
TRUSTED_CHECKOUT=/path/to/protected-harness-checkout
UNTRUSTED_CHECKOUT=/path/to/solution-input-checkout
git -C "$TRUSTED_CHECKOUT" submodule update --init --recursive
TRUSTED_DIR=$RUNNER_TEMP/labinius-trusted
TRUSTED_CARGO_HOME=$RUNNER_TEMP/labinius-trusted-cargo
TRUSTED_TARGET_DIR=$RUNNER_TEMP/labinius-trusted-target
CONTESTANT_DIR=$RUNNER_TEMP/labinius-contestant
CONTESTANT_CARGO_HOME=$RUNNER_TEMP/labinius-contestant-cargo
bash "$TRUSTED_CHECKOUT/.yukon/setup.sh" --track x86
bash "$TRUSTED_CHECKOUT/.yukon/benchmark.sh" --track x86
test -s "$TRUSTED_CHECKOUT/.yukon/labinius-fast-score.json"
rm -rf "$RUNNER_TEMP"
```

The `arm` track splits the same run across two hosts, as its workflow does:
`setup.sh --track arm` and `benchmark.sh --track arm --prove-only DIR` on aarch64, then
`setup.sh --track arm --verifier-only` and `benchmark.sh --track arm --verify-only DIR`
on an AVX-512 x86-64 host after copying `DIR`.

## Disposable-clone validation

The implementation through commit
`84644416716c75b5db9e8f95bcadb160f8ddcb9b` was validated on 2026-10-09 using
a disposable `c7i.2xlarge` clone of private AMI
`ami-05be5ff6eebeffbd7` in account `089601768650`, with an Intel Xeon Platinum
8488C, 8 CPUs, and 16 GiB RAM. The source archive SHA-256 was
`26b9cfb8d20d5e4c68ad3aa0cdbd4a130a5cd4689388b355d078b08f16f2b82d`.

Observed results from that one run:

- sandbox regression, setup, and benchmark all exited 0;
- score/proving time: `2972.2717100000295 ms`;
- proof size: `13,958,612` bytes; modulus: `55144219393` (36 bits);
- benchmark wall time: `5.05 s`; peak RSS: `540,496 KiB`;
- reported `avx512f: true`, `nproc: 8`, and `MemTotal: 16,076,996 KiB`;
- the temporary encrypted S3 transfer was deleted before untrusted execution,
  then the clone and its EBS volume were deleted with no cleanup errors.

This is a functional validation and one performance observation, not a stable
performance baseline or proof of complete isolation. It did not test GitHub
runner-group admission: the group must remain closed until its selected-workflow
restriction explicitly allows the protected reusable workflow reference. The
routing/docs commit after the tested implementation changes no sandbox, setup,
benchmark, verifier, or solution code.
