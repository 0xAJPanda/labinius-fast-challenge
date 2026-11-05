#!/usr/bin/env python3
"""Exercise legacy x86 and sibling ARM paths without requiring runner hardware."""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SANDBOX = ROOT / ".yukon/lib/sandbox.sh"


def shell(code, *args, env=None):
    environment = dict(os.environ)
    environment.pop("LABINIUS_TRACK", None)
    if env:
        environment.update(env)
    return subprocess.run(
        ["bash", "-eu", "-c", 'source "$1"; shift; ' + code, "bash", str(SANDBOX), *args],
        capture_output=True, text=True, env=environment,
    )


class TrackCompatibilityTest(unittest.TestCase):
    def test_no_flag_keeps_x86_even_on_an_arm_host(self):
        result = shell(
            "uname() { echo aarch64; }; labinius_track; labinius_solution_path; labinius_score_file"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), ["x86", "solution", "labinius-fast-score.json"])

    def test_explicit_tracks_and_invalid_selection(self):
        for track, path, score in [
            ("x86", "solution", "labinius-fast-score.json"),
            ("arm", "solution-arm", "labinius-fast-arm-score.json"),
        ]:
            result = shell(
                'labinius_parse_track_args "$@"; labinius_solution_path; labinius_score_file',
                "--track", track,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.splitlines(), [path, score])
        self.assertNotEqual(shell("LABINIUS_TRACK=unknown; labinius_solution_path").returncode, 0)
        self.assertNotEqual(shell('labinius_parse_track_args --track').returncode, 0)

    def test_staging_keeps_the_sibling_out_of_the_build_mount(self):
        for track in ["x86", "arm"]:
            with self.subTest(track=track), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                checkout = root / "input"
                for path in ["solution", "solution-arm"]:
                    (checkout / path / "proof_generate").mkdir(parents=True)
                    (checkout / path / "proof_generate/Cargo.toml").write_text(path)
                # A malformed sibling must not break or enter the selected track's build.
                selected = "solution" if track == "x86" else "solution-arm"
                sibling = "solution-arm" if track == "x86" else "solution"
                (checkout / sibling / "escape").symlink_to("/etc/passwd")
                result = shell(
                    'LABINIUS_TRACK=$1; selected=$(labinius_solution_path); '
                    'python3 "$2/.yukon/stage-solution.py" "$3/$selected" "$RUNNER_TEMP/staged"; '
                    'mkdir -p "$RUNNER_TEMP/staged/proof_generate/target" "$RUNNER_TEMP/cargo"; '
                    # Capture bwrap arguments; the real isolation smoke test runs on Linux.
                    'bwrap() { printf "%s\\n" "$@"; }; '
                    'run_contestant_build "$RUNNER_TEMP/staged" '
                    '"$RUNNER_TEMP/staged/proof_generate/target" "$RUNNER_TEMP/cargo" offline true',
                    track, str(ROOT), str(checkout), env={"RUNNER_TEMP": tmp},
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual((root / "staged/proof_generate/Cargo.toml").read_text(), selected)
                self.assertFalse((root / "staged/escape").exists())
                argv = result.stdout.splitlines()
                self.assertEqual(argv[argv.index("--chdir") + 1], "/workspace/solution/proof_generate")
                self.assertNotIn(str(checkout / sibling), argv)
                self.assertIn("--unshare-net", argv)

    def test_verification_publishes_the_correct_score_and_clears_failed_scores(self):
        # A stub verifier tests orchestration, not cryptographic proof acceptance.
        for track, flags, filename in [
            ("x86", [], "labinius-fast-score.json"),
            ("x86", ["--track", "x86"], "labinius-fast-score.json"),
            ("arm", ["--track", "arm"], "labinius-fast-arm-score.json"),
        ]:
            with self.subTest(flags=flags), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                harness = root / "harness"
                shutil.copytree(ROOT / ".yukon", harness)
                trusted = root / "trusted"
                trusted.mkdir()
                verifier = trusted / "labinius-fast-verify"
                verifier.write_text(
                    "#!/bin/sh\nprintf '%s\\n' "
                    "'{\"field\":\"GHASH\",\"polynomialSize\":2147483648}'\n"
                )
                verifier.chmod(0o755)
                (trusted / "labinius-fast-verify.sha256").write_text(
                    hashlib.sha256(verifier.read_bytes()).hexdigest()
                )
                proofs = root / "proofs"
                proofs.mkdir()
                runs = []
                for n in range(1, 6):
                    (proofs / f"proof-{n}.bin").write_bytes(b"proof")
                    runs.append({"run": n, "seed": f"{n:064x}", "provingTimeMs": n * 10, "proofBytes": 5})
                (proofs / "runs.json").write_text(json.dumps(runs))
                (proofs / "host.json").write_text(json.dumps(
                    {"cpuModel": "test", "nproc": 8, "avx512f": track == "x86", "memTotalKb": 1}
                ))
                # macOS stat differs from the runner's GNU stat; shim only this OS boundary.
                tools = root / "tools"
                tools.mkdir()
                stat = tools / "stat"
                stat.write_text("#!/bin/sh\nwc -c < \"$4\" | tr -d ' '\n")
                stat.chmod(0o755)
                env = dict(os.environ, RUNNER_TEMP=tmp, TRUSTED_DIR=str(trusted),
                           PATH=f"{tools}:{os.environ['PATH']}")
                env.pop("LABINIUS_TRACK", None)
                command = ["bash", str(harness / "benchmark.sh"), *flags, "--verify-only", str(proofs)]
                result = subprocess.run(command, env=env, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                score = json.loads((harness / filename).read_text())
                self.assertEqual(score["score"], 30)
                other = "labinius-fast-arm-score.json" if track == "x86" else "labinius-fast-score.json"
                self.assertFalse((harness / other).exists())
                (proofs / "proof-3.bin").write_bytes(b"")
                result = subprocess.run(command, env=env, text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((harness / filename).exists())


if __name__ == "__main__":
    unittest.main()
