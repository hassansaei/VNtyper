"""Tests for the packaged adVNTR replay subprocess seam."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from tests.unit.test_advntr_calibration_policy import caller_values, capabilities, capture_policy
from vntyper.modules.advntr import advntr_replay as replay
from vntyper.modules.advntr.advntr_calibration_policy import (
    AdvntrToolPin,
    advntr_canonical_sha256,
    capture_policy_for_caller,
    decode_advntr_capabilities,
    decode_capture_policy,
)
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values

pytestmark = pytest.mark.unit


def _json_bytes(document: object) -> bytes:
    return json.dumps(document, sort_keys=True, separators=(",", ":")).encode() + b"\n"


def _plan(tmp_path: Path, *, mode: str = "legacy") -> replay.ReplayPlan:
    capture_root = tmp_path / "capture root"
    capture_root.mkdir(parents=True)
    record = {"schema_version": "advntr-frameshift-capture-v2", "locus": {"vntr_id": 17}}
    capture_raw = _json_bytes(record)
    capture_path = capture_root / "capture one.jsonl"
    capture_path.write_bytes(capture_raw)
    member = replay.ReplayCapture("synthetic-key", "capture one.jsonl", hashlib.sha256(capture_raw).hexdigest(), (17,))
    caller = decode_caller_policy_values(caller_values(mode=mode))
    raw_capture = capture_policy()
    raw_capture["parameters"]["caller_mode"] = mode  # type: ignore[index]
    projected = capture_policy_for_caller(decode_capture_policy(raw_capture), caller)
    background_path = None
    background_sha = None
    if mode == "exact":
        background_path = tmp_path / "background.json"
        background_path.write_bytes(b"{}\n")
        background_sha = hashlib.sha256(background_path.read_bytes()).hexdigest()
    plan = replay.build_replay_plan(
        ("/env with spaces/python", "-m", "advntr"),
        capture_root=capture_root,
        captures=(member,),
        manifest_path=tmp_path / "manifest file.json",
        policy_path=tmp_path / "policy file.json",
        output_directory=tmp_path / "replay output",
        capture_policy=projected,
        caller_policy=caller,
        background_path=background_path,
        background_sha256=background_sha,
    )
    plan.manifest_path.write_bytes(_json_bytes(replay.replay_manifest_document(plan.captures)))
    plan.policy_path.write_bytes(plan.policy_json + b"\n")
    return plan


def _output(plan: replay.ReplayPlan, *, audit: bool = False) -> dict[str, object]:
    manifest_raw = plan.manifest_path.read_bytes()
    policy_raw = plan.policy_path.read_bytes()
    capture = plan.captures[0]
    record = json.loads((plan.capture_root / capture.filename).read_text())
    result = {
        "schema_version": "advntr-frameshift-replay-result-v1",
        "vntr_id": 17,
        "capture_record_sha256": advntr_canonical_sha256(record),
        "policy_sha256": advntr_canonical_sha256(json.loads(policy_raw)),
        "capture_producer": {
            "package_version": "2.4.0",
            "build_id": "a" * 64,
            "source_revision": "b" * 40,
        },
        "capture_assets": {"opaque": "preserved"},
        "loaded_background_sha256": None if plan.background_sha256 is None else "d" * 64,
        "baseline_parity": True,
        "decision_visits": [],
        "calls": [],
        "warnings": ["synthetic warning"],
        "capture_audit": {
            "attribution_outside_trials": ["D1_1"] if audit else [],
            "calibrated_policy_domain_errors": [],
        },
    }
    return {
        "schema_version": "advntr-frameshift-replay-output-v1",
        "manifest_file_sha256": hashlib.sha256(manifest_raw).hexdigest(),
        "manifest_sha256": advntr_canonical_sha256(json.loads(manifest_raw)),
        "policy_file_sha256": hashlib.sha256(policy_raw).hexdigest(),
        "policy_sha256": advntr_canonical_sha256(json.loads(policy_raw)),
        "background_file_sha256": plan.background_sha256,
        "replay_producer": {
            "package_version": "2.4.0",
            "build_id": "a" * 64,
            "source_revision": "b" * 40,
        },
        "results": [
            {"key": "synthetic-key", "capture_sha256": capture.sha256, "vntrs": [{"vntr_id": 17, "result": result}]}
        ],
    }


def test_builds_exact_replay_argv_and_manifest_with_path_spaces(tmp_path: Path) -> None:
    plan = _plan(tmp_path)

    assert plan.argv[:5] == ("/env with spaces/python", "-m", "advntr", "replay-frameshift", "--capture-root")
    assert plan.argv[plan.argv.index("--manifest") + 1] == str(tmp_path / "manifest file.json")
    assert "--background" not in plan.argv
    assert replay.replay_manifest_document(plan.captures) == {
        "schema_version": "advntr-frameshift-replay-manifest-v1",
        "captures": [
            {
                "key": "synthetic-key",
                "filename": "capture one.jsonl",
                "sha256": plan.captures[0].sha256,
                "vntr_ids": [17],
            }
        ],
    }


def test_exact_replay_requires_and_legacy_replay_forbids_background(tmp_path: Path) -> None:
    exact = _plan(tmp_path / "exact", mode="exact")
    assert exact.argv[exact.argv.index("--background") + 1] == str(exact.background_path)
    legacy = _plan(tmp_path / "legacy")
    with pytest.raises(ValueError, match="legacy replay forbids"):
        replay.build_replay_plan(
            legacy.argv_prefix,
            capture_root=legacy.capture_root,
            captures=legacy.captures,
            manifest_path=legacy.manifest_path,
            policy_path=legacy.policy_path,
            output_directory=tmp_path / "other",
            capture_policy=legacy.capture_policy,
            caller_policy=legacy.caller_policy,
            background_path=tmp_path / "background",
            background_sha256="a" * 64,
        )


def test_run_validates_bindings_and_preserves_unassessable_audit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan = _plan(tmp_path)
    monkeypatch.setattr(
        replay,
        "probe_advntr_capabilities",
        lambda *_args, **_kwargs: decode_advntr_capabilities(capabilities()),
    )

    def runner(argv: tuple[str, ...], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        plan.output_directory.mkdir()
        (plan.output_directory / "replay.json").write_bytes(_json_bytes(_output(plan, audit=True)))
        return subprocess.CompletedProcess(argv, 0, "", "")

    result = replay.run_replay(plan, AdvntrToolPin("2.4.0", "a" * 64, "b" * 40), runner=runner)

    assert len(result.loci) == 1
    assert result.loci[0].assessable is False
    preserved = replay.replay_locus_document(result.loci[0])
    assert preserved["capture_audit"]["attribution_outside_trials"] == ["D1_1"]  # type: ignore[index]
    assert preserved["calls"] == []


def test_exact_replay_requires_loaded_background_binding(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    plan = _plan(tmp_path, mode="exact")
    monkeypatch.setattr(
        replay,
        "probe_advntr_capabilities",
        lambda *_args, **_kwargs: decode_advntr_capabilities(capabilities()),
    )

    def runner(argv: tuple[str, ...], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        document = _output(plan)
        document["results"][0]["vntrs"][0]["result"]["loaded_background_sha256"] = None  # type: ignore[index]
        plan.output_directory.mkdir()
        (plan.output_directory / "replay.json").write_bytes(_json_bytes(document))
        return subprocess.CompletedProcess(argv, 0, "", "")

    with pytest.raises(ValueError, match="loaded replay background"):
        replay.run_replay(plan, AdvntrToolPin("2.4.0", "a" * 64, "b" * 40), runner=runner)


@pytest.mark.parametrize("field", ["manifest_file_sha256", "policy_sha256", "replay_producer"])
def test_replay_rejects_wrong_output_bindings(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str) -> None:
    plan = _plan(tmp_path)
    monkeypatch.setattr(
        replay,
        "probe_advntr_capabilities",
        lambda *_args, **_kwargs: decode_advntr_capabilities(capabilities()),
    )

    def runner(argv: tuple[str, ...], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        document = _output(plan)
        document[field] = {} if field == "replay_producer" else "0" * 64
        plan.output_directory.mkdir()
        (plan.output_directory / "replay.json").write_bytes(_json_bytes(document))
        return subprocess.CompletedProcess(argv, 0, "", "")

    with pytest.raises(ValueError):
        replay.run_replay(plan, AdvntrToolPin("2.4.0", "a" * 64, "b" * 40), runner=runner)


def test_failed_replay_is_not_interpreted_as_negative_evidence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    plan = _plan(tmp_path)
    monkeypatch.setattr(
        replay,
        "probe_advntr_capabilities",
        lambda *_args, **_kwargs: decode_advntr_capabilities(capabilities()),
    )
    with pytest.raises(RuntimeError, match="replay process failed"):
        replay.run_replay(
            plan,
            AdvntrToolPin("2.4.0", "a" * 64, "b" * 40),
            runner=lambda argv, **_kwargs: subprocess.CompletedProcess(argv, 1, "", "private stderr"),
        )
