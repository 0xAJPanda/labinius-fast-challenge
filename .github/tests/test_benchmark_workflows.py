#!/usr/bin/env python3
"""Structural regression checks for the per-track benchmark workflows."""

import json
from pathlib import Path
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[2]
WORKFLOWS = ROOT / ".github" / "workflows"
CHECKOUT = "actions/checkout@93cb6efe18208431cddfb8368fd83d5badbf9bfd"
UPLOAD = "actions/upload-artifact@b7c566a772e6b6bfb58ed0dc250532a479d7789f"
DOWNLOAD = "actions/download-artifact@37930b1c2abaa49bbe596cd826c3c89aef350131"
TRUSTED_REPOSITORY = "${{ job.workflow_repository }}"
TRUSTED_REF = "${{ job.workflow_sha }}"
TRACKS = ("x86", "arm")
READ_ONLY = {"contents": "read", "actions": "read"}
SOLUTION_INPUTS = {"solution_repository", "solution_ref"}
PINNED_WORKFLOWS = {
    "x86": (
        "Layr-Labs/labinius-fast-challenge",
        "6ca90ac0fa5e1c6827141a97e2581ae38046a0c3",
    ),
    "arm": (
        "Layr-Labs/labinius-fast-challenge",
        "6ca90ac0fa5e1c6827141a97e2581ae38046a0c3",
    ),
}


def load_workflow(name: str) -> dict:
    # BaseLoader keeps `on` as a string instead of applying YAML 1.1's
    # surprising on/off boolean coercion.
    with (WORKFLOWS / name).open(encoding="utf-8") as workflow_file:
        return yaml.load(workflow_file, Loader=yaml.BaseLoader)


def step_named(workflow: dict, job: str, name: str) -> dict:
    steps = workflow["jobs"][job]["steps"]
    return next(step for step in steps if step["name"] == name)


def score_path(track: str) -> str:
    return f".yukon/labinius-fast-{track}-score.json"


class BenchmarkWorkflowTest(unittest.TestCase):
    def setUp(self) -> None:
        self.workflows = {track: load_workflow(f"benchmark-{track}.yml") for track in TRACKS}
        self.dispatchers = {
            track: load_workflow(f"dispatch-{track}.yml") for track in TRACKS
        }
        with (ROOT / "benchmark.json").open(encoding="utf-8") as manifest_file:
            self.manifest = json.load(manifest_file)

    def test_only_the_track_workflows_and_dispatchers_exist(self) -> None:
        self.assertEqual(
            sorted(path.name for path in WORKFLOWS.glob("*.yml")),
            [
                "benchmark-arm.yml",
                "benchmark-x86.yml",
                "dispatch-arm.yml",
                "dispatch-x86.yml",
            ],
        )

    def test_manifest_declares_one_track_per_solution_directory(self) -> None:
        self.assertEqual(self.manifest["schemaVersion"], 2)
        self.assertEqual(set(self.manifest), {"schemaVersion", "name", "tracks"})
        tracks = {track["name"]: track for track in self.manifest["tracks"]}
        self.assertEqual(list(tracks), list(TRACKS))
        for name, track in tracks.items():
            self.assertEqual(track["editablePaths"], [f"solution/{name}"])
            self.assertEqual(
                track["setupCommand"], ["bash", ".yukon/setup.sh", "--track", name]
            )
            self.assertEqual(
                track["benchmarkCommand"],
                ["bash", ".yukon/benchmark.sh", "--track", name],
            )
            self.assertEqual(track["scorePath"], score_path(name))
            self.assertEqual(track["direction"], "-")
            self.assertEqual(
                track["runner"],
                {
                    "provider": "github-actions",
                    "workflow": f"dispatch-{name}.yml",
                    "maxConcurrentWorkflows": 1,
                },
            )

    def test_workflows_are_dispatchable_and_callable(self) -> None:
        for track, workflow in self.workflows.items():
            self.assertEqual(set(workflow["on"]), {"workflow_dispatch", "workflow_call"})
            for trigger in workflow["on"].values():
                self.assertEqual(set(trigger["inputs"]), SOLUTION_INPUTS)
                for spec in trigger["inputs"].values():
                    self.assertEqual(spec["required"], "false")
            self.assertEqual(workflow["permissions"], READ_ONLY)
            self.assertEqual(
                workflow["concurrency"],
                {
                    "group": f"benchmark-{track}-${{{{ github.repository }}}}-"
                    "${{ github.workflow }}-${{ github.ref }}-"
                    "${{ inputs.solution_ref || github.sha }}",
                    "cancel-in-progress": "true",
                },
            )

    def test_harness_is_pinned_independently_of_the_caller(self) -> None:
        for workflow in self.workflows.values():
            for name in workflow["jobs"]:
                trusted = step_named(workflow, name, "Checkout trusted harness")
                self.assertEqual(trusted["uses"].split(" ")[0], CHECKOUT)
                self.assertEqual(trusted["with"]["repository"], TRUSTED_REPOSITORY)
                self.assertEqual(trusted["with"]["ref"], TRUSTED_REF)
                self.assertEqual(trusted["with"]["path"], ".trusted-src")
                self.assertEqual(trusted["with"]["persist-credentials"], "false")
                paths = step_named(workflow, name, "Resolve isolated paths")["run"]
                self.assertIn("TRUSTED_CHECKOUT=$GITHUB_WORKSPACE/.trusted-src", paths)
                self.assertIn("UNTRUSTED_CHECKOUT=$GITHUB_WORKSPACE/.untrusted-src", paths)
                for step in workflow["jobs"][name]["steps"]:
                    if ".yukon/" in step.get("run", ""):
                        self.assertIn('bash "$TRUSTED_CHECKOUT/.yukon/', step["run"])

    def test_solution_is_checked_out_separately(self) -> None:
        jobs = {"x86": ["benchmark"], "arm": ["prove"]}
        for track, names in jobs.items():
            workflow = self.workflows[track]
            for name in names:
                solution = step_named(workflow, name, "Checkout untrusted solution input")
                self.assertEqual(
                    solution["with"],
                    {
                        "repository": "${{ inputs.solution_repository || github.repository }}",
                        "ref": "${{ inputs.solution_ref || github.sha }}",
                        "path": ".untrusted-src",
                        "fetch-depth": "1",
                        "submodules": "false",
                        "persist-credentials": "false",
                    },
                )
        verify_steps = [step["name"] for step in self.workflows["arm"]["jobs"]["verify"]["steps"]]
        self.assertNotIn("Checkout untrusted solution input", verify_steps)

    def test_x86_proves_and_verifies_on_the_dedicated_host(self) -> None:
        workflow = self.workflows["x86"]
        self.assertEqual(set(workflow["jobs"]), {"benchmark"})
        job = workflow["jobs"]["benchmark"]
        self.assertEqual(
            job["runs-on"],
            {
                "group": "labinius-intel-8488c-avx512-c7i",
                "labels": "labinius-avx512-c7i",
            },
        )
        self.assertEqual(
            step_named(workflow, "benchmark", "Setup")["run"],
            'bash "$TRUSTED_CHECKOUT/.yukon/setup.sh" --track x86',
        )
        self.assertEqual(
            step_named(workflow, "benchmark", "Benchmark")["run"],
            'bash "$TRUSTED_CHECKOUT/.yukon/benchmark.sh" --track x86',
        )

    def test_arm_proves_on_arm_and_verifies_on_x86(self) -> None:
        workflow = self.workflows["arm"]
        self.assertEqual(set(workflow["jobs"]), {"prove", "verify"})
        prove, verify = workflow["jobs"]["prove"], workflow["jobs"]["verify"]
        self.assertEqual(
            prove["runs-on"],
            {"group": "labinius-graviton4-neoverse-v2-sve2-c8g", "labels": ["ARM64", "sve2"]},
        )
        self.assertIn(
            "/usr/local/bin/validate-graviton4-cpu",
            step_named(workflow, "prove", "Require dedicated Graviton4 host")["run"],
        )
        # The dedicated host ships the toolchain; only the hosted x86 verifier provisions one.
        prove_runs = [step.get("run", "") for step in prove["steps"]]
        self.assertFalse(any("provision.sh" in run for run in prove_runs))
        self.assertIn(
            "provision.sh", step_named(workflow, "verify", "Setup toolchain")["run"]
        )
        self.assertEqual(verify["runs-on"], "blacksmith-8vcpu-ubuntu-2404")
        self.assertEqual(verify["needs"], "prove")
        self.assertEqual(
            step_named(workflow, "prove", "Setup")["run"],
            'bash "$TRUSTED_CHECKOUT/.yukon/setup.sh" --track arm',
        )
        self.assertEqual(
            step_named(workflow, "prove", "Benchmark proofs")["run"],
            'bash "$TRUSTED_CHECKOUT/.yukon/benchmark.sh" --track arm --prove-only "$RUNNER_TEMP/labinius-arm-runs"',
        )
        self.assertEqual(
            step_named(workflow, "verify", "Setup verifier")["run"],
            'bash "$TRUSTED_CHECKOUT/.yukon/setup.sh" --track arm --verifier-only',
        )
        self.assertEqual(
            step_named(workflow, "verify", "Benchmark verification")["run"],
            'bash "$TRUSTED_CHECKOUT/.yukon/benchmark.sh" --track arm --verify-only "$RUNNER_TEMP/labinius-arm-runs"',
        )
        upload = step_named(workflow, "prove", "Upload proofs")
        download = step_named(workflow, "verify", "Download proofs")
        self.assertEqual(upload["uses"].split(" ")[0], UPLOAD)
        self.assertEqual(download["uses"].split(" ")[0], DOWNLOAD)
        self.assertEqual(upload["with"]["name"], download["with"]["name"])
        self.assertEqual(upload["with"]["path"], download["with"]["path"])
        self.assertEqual(prove["outputs"]["proofs"], "${{ steps.upload.outputs.artifact-id }}")
        self.assertEqual(verify["permissions"], {"contents": "read", "actions": "write"})
        delete = step_named(workflow, "verify", "Delete proofs")
        self.assertEqual(delete["if"], "always()")
        self.assertEqual(delete["env"]["ARTIFACT_ID"], "${{ needs.prove.outputs.proofs }}")

    def test_score_artifact_carries_the_manifest_score_path(self) -> None:
        jobs = {"x86": "benchmark", "arm": "verify"}
        for track, job in jobs.items():
            upload = step_named(self.workflows[track], job, "Upload score")
            self.assertEqual(upload["if"], "always()")
            self.assertEqual(upload["uses"].split(" ")[0], UPLOAD)
            self.assertEqual(upload["with"]["name"], f"benchmark-score-{track}")
            self.assertEqual(
                upload["with"]["path"],
                f"${{{{ github.workspace }}}}/.trusted-src/{score_path(track)}",
            )
            self.assertEqual(upload["with"]["if-no-files-found"], "error")

    def test_dispatchers_pin_reviewed_executors_and_forward_the_solution(self) -> None:
        for track in TRACKS:
            dispatcher = self.dispatchers[track]
            self.assertEqual(dispatcher["on"], {"workflow_dispatch": ""})
            self.assertEqual(dispatcher["permissions"], READ_ONLY)
            self.assertEqual(set(dispatcher["jobs"]), {"benchmark"})
            job = dispatcher["jobs"]["benchmark"]
            repository, sha = PINNED_WORKFLOWS[track]
            self.assertEqual(
                job["uses"],
                f"{repository}/.github/workflows/benchmark-{track}.yml@{sha}",
            )
            self.assertEqual(
                job["with"],
                {
                    "solution_repository": "${{ github.repository }}",
                    "solution_ref": "${{ github.sha }}",
                },
            )
        # The called arm workflow deletes its proof artifact, which the caller must allow.
        self.assertEqual(
            self.dispatchers["arm"]["jobs"]["benchmark"]["permissions"],
            {"contents": "read", "actions": "write"},
        )


if __name__ == "__main__":
    unittest.main()
