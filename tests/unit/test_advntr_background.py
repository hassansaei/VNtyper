"""Tests for the packaged adVNTR recipe-v1 background-fit adapter."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from tests.unit.test_advntr_calibration_policy import capabilities
from vntyper.modules.advntr import advntr_background as background
from vntyper.modules.advntr.advntr_calibration_policy import AdvntrToolPin, decode_advntr_capabilities

pytestmark = pytest.mark.unit


def _policy() -> dict[str, object]:
    return {
        "schema_version": "advntr-frameshift-policy-v1",
        "mode": "exact",
        "cutoff": 0.001,
        "minimum_read_support": 3,
    }


def _plan(tmp_path: Path) -> background.BackgroundFitPlan:
    policy = tmp_path / "diagnostic policy.json"
    policy.write_text(json.dumps(_policy()))
    return background.build_background_fit_plan(
        ("/env with spaces/python", "-m", "advntr"),
        capture_root=tmp_path / "capture root",
        labels_path=tmp_path / "training labels.json",
        diagnostic_policy_path=policy,
        output_directory=tmp_path / "fit output",
        partition="training",
        profile="candidate-a",
        folds=5,
        insert_lengths=8,
        source_cohort="SIMULATED",
        design="synthetic preregistered design",
    )


def _write_outputs(plan: background.BackgroundFitPlan) -> None:
    plan.output_directory.mkdir()
    documents: dict[str, object] = {
        "candidate-a.background.json": {
            "schema": "advntr.frameshift.background",
            "version": 1,
            "provenance": "SYNTHETIC",
            "default_probability": 0.1,
            "states": {"D1_1": 0.05},
        },
        "candidate-a.sidecar.json": {
            "schema": "advntr-bench.frameshift.background.sidecar",
            "version": 1,
            "profile_name": "candidate-a",
            "partition": "training",
            "background_recipe_id": "recipe-v1",
            "diagnostic_policy": _policy(),
            "preregistration_overrides": {},
            "capture_identity": {
                "producer": {
                    "package_version": "2.4.0",
                    "build_id": "a" * 64,
                    "source_revision": "b" * 40,
                }
            },
            "independent_control_groups": 15,
            "loader_proof": {"bad_key_probe": {"refused": True}},
        },
        "candidate-a.build-report.json": {"background_recipe_id": "recipe-v1", "diagnostic_policy": _policy()},
        "candidate-a.cv.json": {"fold_count": 5},
        "candidate-a.falsification.json": {"synthetic": True},
        "candidate-a.predictions.json": {"synthetic": True},
    }
    for name, document in documents.items():
        (plan.output_directory / name).write_text(json.dumps(document))
    (plan.output_directory / "candidate-a.states.tsv").write_text("state\tp\n")
    (plan.output_directory / "candidate-a.build-report.md").write_text("# synthetic\n")
    (plan.output_directory / "loader-refusal-probe.json").write_text("{}\n")


def test_builds_frozen_recipe_argv_without_shell_or_smoke_overrides(tmp_path: Path) -> None:
    plan = _plan(tmp_path)

    assert plan.argv[:5] == ("/env with spaces/python", "-m", "advntr", "fit-background", "--capture-root")
    assert plan.argv[plan.argv.index("--background-recipe") + 1] == "recipe-v1"
    assert plan.argv[plan.argv.index("--diagnostic-policy") + 1] == str(tmp_path / "diagnostic policy.json")
    assert "--screen-min-samples" not in plan.argv
    assert "--worktree" not in plan.argv
    assert "--note" not in plan.argv


@pytest.mark.parametrize(("field", "value"), [("folds", True), ("partition", "bad space"), ("profile", "../bad")])
def test_plan_rejects_ambiguous_or_wrong_typed_values(tmp_path: Path, field: str, value: object) -> None:
    kwargs: dict[str, object] = {
        "capture_root": tmp_path / "capture",
        "labels_path": tmp_path / "labels",
        "diagnostic_policy_path": tmp_path / "policy",
        "output_directory": tmp_path / "output",
        "partition": "training",
        "profile": "candidate",
        "folds": 5,
        "insert_lengths": 8,
        "source_cohort": "SIMULATED",
        "design": "synthetic",
    }
    kwargs[field] = value
    with pytest.raises(ValueError):
        background.build_background_fit_plan(("advntr",), **kwargs)  # type: ignore[arg-type]


def test_run_preflights_and_validates_exact_fit_inventory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    plan = _plan(tmp_path)
    plan.capture_root.mkdir()
    plan.labels_path.write_text("{}")
    events: list[str] = []

    def probe(*_args: object, **_kwargs: object) -> object:
        events.append("probe")
        return decode_advntr_capabilities(capabilities())

    monkeypatch.setattr(
        background,
        "probe_advntr_capabilities",
        probe,
    )

    def runner(argv: tuple[str, ...], **kwargs: object) -> subprocess.CompletedProcess[str]:
        events.append("fit")
        _write_outputs(plan)
        return subprocess.CompletedProcess(argv, 0, "{}\n", "")

    result = background.fit_background(
        plan,
        AdvntrToolPin("2.4.0", "a" * 64, "b" * 40),
        runner=runner,
    )

    assert events == ["probe", "fit"]
    assert set(result.artifact_sha256) == set(background.fit_output_names("candidate-a"))
    assert (
        result.background_sha256
        == hashlib.sha256((plan.output_directory / "candidate-a.background.json").read_bytes()).hexdigest()
    )


def test_fit_rejects_missing_extra_or_mismatched_recipe_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    plan = _plan(tmp_path)
    plan.capture_root.mkdir()
    plan.labels_path.write_text("{}")
    monkeypatch.setattr(
        background,
        "probe_advntr_capabilities",
        lambda *_args, **_kwargs: decode_advntr_capabilities(capabilities()),
    )

    def runner(argv: tuple[str, ...], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        _write_outputs(plan)
        sidecar = plan.output_directory / "candidate-a.sidecar.json"
        document = json.loads(sidecar.read_text())
        document["background_recipe_id"] = "changed"
        sidecar.write_text(json.dumps(document))
        return subprocess.CompletedProcess(argv, 0, "", "")

    with pytest.raises(ValueError, match="recipe"):
        background.fit_background(plan, AdvntrToolPin("2.4.0", "a" * 64, "b" * 40), runner=runner)


def test_failed_fit_is_not_decoded_as_background_evidence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    plan = _plan(tmp_path)
    plan.capture_root.mkdir()
    plan.labels_path.write_text("{}")
    monkeypatch.setattr(
        background,
        "probe_advntr_capabilities",
        lambda *_args, **_kwargs: decode_advntr_capabilities(capabilities()),
    )
    with pytest.raises(RuntimeError, match="background fit process failed"):
        background.fit_background(
            plan,
            AdvntrToolPin("2.4.0", "a" * 64, "b" * 40),
            runner=lambda argv, **_kwargs: subprocess.CompletedProcess(argv, 1, "", "private stderr"),
        )
