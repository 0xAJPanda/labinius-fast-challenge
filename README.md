# labinius.fast

Make a packed-binary Labinius polynomial commitment opening faster.

## The task

Each trial starts with a fresh 32-byte seed. It determines a table of **2³¹ Boolean evaluations**, packed into **2²⁴ elements of the GHASH field** GF(2¹²⁸), 128 evaluations per element. The trusted harness generates this table before timing begins and supplies it to the prover as a read-only file.

Your prover must commit to that table, reduce the GHASH evaluation claim to a single F162 claim (the **GHASH → F162 field switch**), and then open the Labinius commitment at the resulting F162 point. The trusted verifier independently reconstructs the table from the seed and checks the proof. You optimize the prover implementation, starting from your track's baseline prover; the instance shape and verifier contract stay fixed.

## Tracks

This challenge has two independently scored tracks. They share the instance, proof format, trusted verifier, and scoring rules, and differ only in hardware and in which directory you edit:

| Track | Editable path | Scored on | Baseline |
| --- | --- | --- | --- |
| `x86` | `solution/x86` | Intel Xeon Platinum 8488C, 8 CPUs, AVX-512 | AVX-512 Labinius kernels with a precomputed commitment key |
| `arm` | `solution/arm` | AWS Graviton4 (Neoverse V2, SVE2), 8 vCPUs, no AVX-512 (proofs verified on x86) | Labinius with portable fallbacks for every AVX-512 intrinsic, and a grouped-table mod-`Q` commitment |

A submission to one track may change only that track's directory. Each track's prover is built from its own `solution/<track>/proof_generate` crate, which must produce the `labinius-fast` executable. Scores from different tracks measure different hardware and are not comparable.

## Scoring

The score is the **median wall-clock time, in milliseconds, of five sequential prover runs**, each using a fresh seed. Lower is better. For sorted run times `t₁ ≤ t₂ ≤ t₃ ≤ t₄ ≤ t₅`, the score is `t₃`.

Timing includes:

- Prover process startup and loading the input table.
- Commitment generation and opening proof generation.
- Writing the serialized proof file.

Compilation, trusted table generation, and verification are excluded.

**All five runs must succeed and produce accepted proofs.** A failed run, timeout, invalid proof, empty proof, or proof larger than **64 MiB** prevents score publication. Each prover invocation has a **900-second timeout**. The harness records each run's time and proof size, together with the verified commitment modulus. Proof size and modulus are supporting metrics, not separate ranking objectives.

Compare results only when they use the same benchmark version and hardware. Preserve all five timings when reporting an experiment; a single fast trial is not the benchmark score.

## Commitment parameters

The commitment uses the conductor-1944 cyclotomic ring of degree 648:

```text
R_Q = (Z/QZ)[X] / (X^648 - X^324 + 1)
34,816,965,049 ≤ Q < 2^96
```

The lower bound comes from the challenge's module-BKZ security estimate. The verifier also requires `Q > 1,083,574`, twice the largest coefficient an accepted folded witness can have, so two different folded witnesses are never congruent modulo `Q`; the lower bound already implies this. The upper bound is a limit of the verifier's 128-bit arithmetic, not of the scheme.

**You are explicitly allowed to change the commitment modulus `Q` as part of your optimization.** The baseline prover's default `Q` is not fixed by the challenge. You may choose any integer `Q` within the bounds above that satisfies the trusted verifier's acceptance checks, including a power of two. The verifier does not require `Q` to be prime or have a particular factorization.

You may modify the modulus selection, arithmetic, and serialization implementation under your track's `solution/<track>` directory to support your chosen `Q`. The baseline prover's use of products of supported primes is an implementation choice you may replace, not a challenge restriction. The resulting proof must retain the verifier's required format, carry the chosen `Q`, and pass the unchanged trusted verifier.

The commitment modulus `Q` is separate from the GHASH field polynomial. Changing `Q` does not change the input field, table size, or opening being checked.

## What you may change

Only your track's directory is editable: `solution/x86` for the `x86` track, `solution/arm` for the `arm` track, each including its vendored Labinius crates. The other track's directory is read-only to you. You may improve algorithms, arithmetic, memory layout, allocation, vectorization, threading, and serialization, provided every output remains compatible with the trusted verifier.

The trusted harness and verifier are checked out separately from the submission. Changes to `verifier/`, `.yukon/`, `.github/`, or `benchmark.json` are not solver changes and must not be used to make a candidate pass. Do not change the instance shape, weaken verification, alter the score calculation, or rely on proofs cached from earlier instances.

Contestant builds and prover execution are sandboxed. Dependency fetching is a separate network-enabled phase; compilation and proving run without network access. The input table is read-only.

## Hardware and setup

Official `x86` measurements require an **Intel Xeon Platinum 8488C**, **eight available CPUs**, and **AVX-512**. Official `arm` measurements run on an **AWS Graviton4** (Neoverse V2 with SVE2) with **eight vCPUs** (Ubuntu 24.04); its proofs are then checked on x86 by the same pinned verifier as the `x86` track. Binaries are compiled with `-C target-cpu=native`. Measurements on different CPUs or with different CPU availability are useful for development but are not equivalent official results.

The runner environment supplies Rust/Cargo, GCC 13, bubblewrap, Python 3, OpenSSL, and standard Linux tools. The protected setup scripts define the toolchain locations and sandbox requirements. `yukon setup` needs Ubuntu 24.04 on a host matching the track: x86-64 with AVX-512 for `x86`, aarch64 for `arm`. An Apple Silicon Mac cannot run `yukon setup`, but it can build, prove, and verify the `arm` track outside the sandbox with `.yukon/run-prover.sh --track arm --runs 5`. That loop verifies with a development build of the verifier linked against `solution/arm/labinius`. It is useful for development but is not an official measurement.

Use the installation and authentication instructions from the challenge's **Participate** dialog. With Yukon installed:

```bash
yukon benchmark list
yukon benchmark show <benchmark-id-or-name>
yukon clone <benchmark-id-or-name>
```

Replace the placeholder with the registered challenge identifier. `yukon clone <setter>/<challenge>` selects the first track (`x86`); `yukon clone <setter>/<challenge>/arm` selects `arm`. Follow the CLI's printed `cd` instruction into the benchmark work directory, then read `benchmark.json` and the installed Yukon CLI skill before running:

```bash
yukon tracks
yukon setup --track <x86|arm>
yukon run --track <x86|arm>
yukon submissions --all
```

This challenge's schema-v2 manifest declares both tracks on one shared branch, with the non-overlapping editable paths `solution/x86` and `solution/arm`. Pass `--track` to `yukon setup`, `yukon run`, and `yukon submit`; `yukon submit --track <track>` archives only that track's directory. `yukon switch <track>` changes only Yukon's local track selection: it does not check out a branch or change your `HEAD`, index, or worktree, so your local branches and uncommitted work stay in place. Switch before editing the other track so traces and commands without `--track` use the right one. Each track writes its own score file, `.yukon/labinius-fast-<track>-score.json`. Use a supported runner environment. If local setup reports missing trusted-checkout paths or sandbox prerequisites, resolve the environment rather than editing the protected scripts or disabling isolation.

The upcoming harness defines how the prover receives the input file and where it writes its proof. Read the checked-out `.yukon/run-prover.sh` and prover argument parser for the exact interface; do not assume an older seed-only invocation still applies.

## Repository map

| Path | Purpose |
| --- | --- |
| `benchmark.json` | Schema-v2 manifest: both tracks' editable paths, commands, score paths, and workflows. |
| `solution/x86/proof_generate/` | `x86` prover crate producing the `labinius-fast` executable. |
| `solution/x86/proof_generate/src/parameters.rs` | `x86` baseline prover's instance constants and default modulus. |
| `solution/x86/proof_generate/src/modulus_choice.rs` | Moduli supported by the `x86` baseline prover. |
| `solution/x86/labinius/` | `x86` track's editable vendored implementation (AVX-512). |
| `solution/arm/proof_generate/` | `arm` prover crate producing the `labinius-fast` executable. |
| `solution/arm/labinius/` | `arm` track's editable vendored implementation, with portable intrinsics in `crates/pcs/src/arch/`. |
| `verifier/proof_verify/` | Trusted verifier and independent commitment-ring arithmetic. |
| `verifier/labinius/` | Pinned verifier dependency (x86-64 only). Both tracks' proofs are verified with it on x86. |
| `verifier/table_writer/` | Trusted, dependency-light writer of the input table, used on ARM, where the pinned verifier cannot run. |
| `.yukon/` | Trusted setup, input staging, sandboxing, execution, and scoring. |
| `.yukon/run-prover.sh` | One sandboxed prover run for the harness; with `--track`, a local build, prove, and verify loop. |
| `.github/workflows/benchmark-x86.yml` | `x86` track's dispatch workflow: one job that sets up, proves, verifies, and uploads the score. |
| `.github/workflows/benchmark-arm.yml` | `arm` track's dispatch workflow: proves on ARM, then verifies and scores on x86. |
| `.github/workflows/test-benchmark.yml` | Pull-request check that runs both track workflows. |

Each track has one self-contained workflow, named by its manifest `runner.workflow`.
Yukon dispatches it without inputs. Each job runs the trusted harness from its own checkout of
the dispatched commit and benchmarks the solution from a separate checkout. That checkout uses
the same commit by default, or the optional `solution_repository`/`solution_ref` inputs.
The pull-request check calls both workflows the same way, without inputs, so a pull request
is benchmarked exactly as a submission would be. Concurrency keys are prefixed with
the track, so a pull request's two calls do not cancel each other.

The `x86` workflow runs on the admission-controlled runner group
`labinius-intel-8488c-avx512-c7i`. The `arm` workflow's proving job runs on the
admission-controlled group `labinius-graviton4-neoverse-v2-sve2-c8g`: Graviton4,
Neoverse V2 with SVE2, 8 vCPUs. Each group's selected-workflow restriction must admit its
track's workflow. The ARM job uploads the five proofs with their seeds, timings,
and host description as a short-lived artifact. A second job, on a Blacksmith x86 runner with
AVX-512 (`blacksmith-8vcpu-ubuntu-2404`), verifies them with the pinned verifier,
writes the score, and deletes the proof artifact. Only the score artifact remains,
which keeps the run within Yukon's artifact size limit.

## Optimization workflow

1. Record the track, checkout revision, benchmark version, hardware, current promoted frontier, and baseline five-run result.
2. Read the prover, verifier, and harness. Distinguish implementation choices from verifier requirements.
3. Form a specific performance hypothesis and make a focused change under `solution/<track>`.
4. Run relevant implementation tests, then `yukon run --track <track>` in the supported environment. Confirm all five proofs are accepted.
5. Compare medians under identical conditions. Repeat measurements when the gain is small relative to run-to-run variation.
6. Record the change, commands, all five times, median, proof sizes, modulus, and any limitations.

Per-phase prover logs help identify bottlenecks, but only the harness's external measurement determines the score. Do not omit slow trials, substitute an internal timer, or compare an old single-run baseline with a new five-run result.

## Submission

Use `yukon submit --help` and the installed Yukon CLI skill for current submission requirements. Submit from the linked benchmark work directory with `--track <track>`, a public note file, the exact underlying model, and the agent or harness actually used. Declare reused submissions through the CLI's dependency mechanism; do not copy predecessor attribution.

A useful note includes the base revision, benchmark version, implementation changes, reproducible commands, measured results, failed experiments that inform the result, and known limitations. Keep credentials and private information out of public notes. Official scoring and promotion are determined by the trusted remote evaluation, not a local timing claim.
