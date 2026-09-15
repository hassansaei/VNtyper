"""Caller controller binds policy evidence, payloads, and exposure ordering."""

import hashlib
import json
from argparse import Namespace
from dataclasses import replace
from importlib import import_module
from pathlib import Path
from types import MappingProxyType
from unittest.mock import patch

import pytest

from tests.unit.test_advntr_calibration_policy import capture_policy
from tests.unit.test_calibration_callers import _evidence, _protocol, _roster
from tests.unit.test_calibration_portable_background import background
from tests.unit.test_calibration_target_contract import study_document
from tests.unit.test_calibration_target_runs import run_document
from vntyper.scripts.calibration_advntr_runtime_policy import decode_advntr_runtime_policy
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_protocol import CallerProtocol
from vntyper.scripts.calibration_callers import CallerEvidenceRow
from vntyper.scripts.calibration_portable_background import project_portable_background
from vntyper.scripts.calibration_target_contract import CallerBaselinePlan, decode_target_study
from vntyper.scripts.calibration_target_runs import TargetRun, TargetRunAsset
from vntyper.scripts.canonical_json import canonical_json_bytes

pytestmark = pytest.mark.unit


def _caller_study():
    raw = study_document("callers")
    raw["baseline"]["producer"]["tool_versions"].update(advntr="2.4.1", advntr_build_id="b" * 64)
    return decode_target_study(raw)


def test_role_evidence_binding_recomputes_digest_for_actual_assets() -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    source = _evidence()

    bound = module.bind_caller_role_evidence(
        "validation", "1" * 64, source.eligible_roster_sha256, "2" * 64,
        "3" * 64, source.policies, source.replay_equivalence,
    )

    assert bound.phase == "validation"
    assert bound.run_manifest_sha256 == "2" * 64
    assert bound.sha256 != source.sha256


def test_policy_evidence_preserves_full_grid_and_proves_native_replay_parity() -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    protocol, roster = _protocol(), _roster()
    policy_ids = tuple(sorted({protocol.baseline_policy_sha256, *(item.candidate_id for item in protocol.candidates)}))
    selected = tuple(
        type(
            "Run",
            (),
            {
                "manifest_key": member.key,
                "policy_sha256": policy_id,
                "capture_policy_sha256": "a" * 64,
                "execution_kind": (
                    "baseline-rerun" if policy_id == protocol.baseline_policy_sha256 else "scalar-replay"
                ),
                "assets": {"kestrel_capture": object()},
            },
        )()
        for policy_id in policy_ids
        for member in roster.members
    )
    capture = type("Capture", (), {"provenance": type("P", (), {"capture_policy_sha256": "a" * 64})()})()

    def replay(member, *_args, **_kwargs):
        observation = CallerObservation(member.key, member.group_key, False, (), False, (), ())
        return CallerEvidenceRow(observation, "zero-candidate", "d" * 64)

    with (
        patch.object(module, "require_run_assets"),
        patch.object(module, "read_target_json", return_value={}),
        patch.object(module, "decode_kestrel_capture", return_value=capture),
        patch.object(module, "_advntr_baseline", return_value=((), True, b"capture")),
        patch.object(module, "_advntr_replay", return_value=((), True)),
        patch.object(module, "replayed_caller_observation", side_effect=replay),
        patch.object(module, "_native_row", side_effect=lambda _run, _member, _truth, row, _callers: row),
        patch.object(module, "_source_digest", return_value="d" * 64),
    ):
        policies, equivalence = module._policy_evidence(protocol, roster, object(), selected, policy_ids)

    assert tuple(policy.candidate_id for policy in policies) == policy_ids
    assert len(equivalence) == 1
    assert equivalence[0].baseline_rerun_sha256 == next(
        policy.sha256 for policy in policies if policy.execution_kind == "baseline-rerun"
    )


def test_payload_freezes_selected_profile_runtime_policy_and_portable_background(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    study = _caller_study()
    assert isinstance(study.protocol, CallerProtocol)
    assert isinstance(study.baseline, CallerBaselinePlan)
    policy = study.protocol.candidates[0].policy
    raw = canonical_json_bytes(project_portable_background(background()))
    path = tmp_path / "background.json"
    path.write_bytes(raw)
    capture_path = tmp_path / "capture.jsonl"
    capture_path.write_bytes(b"capture")
    run = TargetRun(
        "artifact-1",
        policy.sha256,
        "1" * 64,
        "2" * 64,
        study.baseline.assets_sha256,
        "3" * 64,
        "scalar-replay",
        0,
        MappingProxyType(
            {
                "advntr_model": TargetRunAsset(tmp_path / "model", "4" * 64, 1),
                "advntr_background": TargetRunAsset(path, hashlib.sha256(raw).hexdigest(), len(raw)),
                "advntr_capture": TargetRunAsset(
                    capture_path, hashlib.sha256(b"capture").hexdigest(), len(b"capture")
                ),
            }
        ),
        (17,),
    )
    source = type("Source", (), {"sha256": "5" * 64})()
    fitted = type("Background", (), {"background_bytes": raw, "background_sha256": hashlib.sha256(raw).hexdigest()})()

    with patch("vntyper.scripts.calibration_caller_profile_io._capture_records", return_value=({"capture_policy": capture_policy(), "producer": {"package_version": "2.4.1", "build_id": "b" * 64, "source_revision": "c" * 40}, "assets": {"model_sha256": "4" * 64}},)):
        manifest, descriptor, files, profile = module._payload(study, policy, source, (run,), fitted)

    assert descriptor.required_callers == ("advntr", "kestrel")
    assert set(files) == {"advntr-policy.json", "background.json", "caller-bundle.json", "decision-profile.json"}
    assert b"invented source label" not in files["background.json"]
    assert descriptor.decision_profile_sha256 == profile.digest
    assert {item.path for item in manifest.files} == set(files)
    runtime = decode_advntr_runtime_policy(json.loads(files["advntr-policy.json"]))
    assert runtime.capture_policy_sha256 != run.capture_policy_sha256
    assert runtime.advntr_revision == "c" * 40
    assert runtime.advntr_revision != study.baseline.producer.source_revision


def test_payload_rejects_unprojected_background_before_candidate_freeze(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    study = _caller_study()
    assert isinstance(study.protocol, CallerProtocol)
    assert isinstance(study.baseline, CallerBaselinePlan)
    policy = study.protocol.candidates[0].policy
    raw = canonical_json_bytes(background())
    path = tmp_path / "background.json"
    path.write_bytes(raw)
    capture_path = tmp_path / "capture.jsonl"
    capture_path.write_bytes(b"capture")
    run = TargetRun(
        "artifact-1", policy.sha256, "1" * 64, "2" * 64, study.baseline.assets_sha256, "3" * 64,
        "scalar-replay", 0,
        MappingProxyType({
            "advntr_model": TargetRunAsset(tmp_path / "model", "4" * 64, 1),
            "advntr_background": TargetRunAsset(path, hashlib.sha256(raw).hexdigest(), len(raw)),
            "advntr_capture": TargetRunAsset(
                capture_path, hashlib.sha256(b"capture").hexdigest(), len(b"capture")
            ),
        }),
        (17,),
    )
    with (
        patch("vntyper.scripts.calibration_caller_profile_io._capture_records", return_value=({"capture_policy": capture_policy(), "producer": {"package_version": "2.4.1", "build_id": "b" * 64, "source_revision": "c" * 40}, "assets": {"model_sha256": "4" * 64}},)),
        pytest.raises(ValueError, match="portable"),
    ):
        module._payload(
            study, policy, type("Source", (), {"sha256": "5" * 64})(), (run,),
            type("Background", (), {"background_bytes": raw, "background_sha256": hashlib.sha256(raw).hexdigest()})(),
        )


def test_wrong_fit_objective_fails_before_role_metadata_or_exposure(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    study_raw = study_document("callers")
    runs_raw = run_document(target="callers")
    args = Namespace(evidence=tmp_path / "evidence", exposure_ledger=tmp_path / "ledger", objective="length-total-v1")
    output = tmp_path / "output"
    output.mkdir()
    with (
        patch.object(module, "load_object", side_effect=[study_raw, runs_raw]),
        patch.object(module, "_source", side_effect=AssertionError("role opened")) as source,
        patch.object(module, "record_exposure", side_effect=AssertionError("exposed")) as exposure,
        pytest.raises(ValueError, match="objective"),
    ):
        module.fit_caller_bundle(args, output)
    source.assert_not_called()
    exposure.assert_not_called()
    assert not tuple(output.iterdir())


def test_failed_selection_is_a_complete_reported_outcome(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    study = _caller_study()
    runs = type("Runs", (), {"sha256": "6" * 64})()
    roster = type("Roster", (), {"sha256": "7" * 64})()
    training = type("Source", (), {"sha256": "8" * 64, "roster": roster})()
    selection = type("Source", (), {"sha256": "9" * 64, "roster": roster, "keys": ("artifact-1",)})()
    result = type(
        "Result",
        (),
        {"selection": type("Selection", (), {"status": "failed", "selected_candidate_id": None})()},
    )()
    output = tmp_path / "output"
    output.mkdir()
    args = Namespace(
        evidence=tmp_path / "evidence", exposure_ledger=tmp_path / "ledger", objective="caller-safety-v1",
        advntr_executable=tmp_path / "advntr",
    )
    with (
        patch.object(module, "load_object", side_effect=[{}, {}]),
        patch.object(module, "decode_target_study", return_value=study),
        patch.object(module, "decode_target_runs", return_value=runs),
        patch.object(module, "_source", side_effect=[training, selection]),
        patch.object(module, "_expose", return_value=object()),
        patch.object(module, "load_caller_source_evidence", return_value=object()),
        patch.object(
            module, "fit_caller_training_background",
            return_value=type("Background", (), {"training_evidence_sha256": "a" * 64})(),
        ),
        patch.object(module, "_require_training_background_arms"),
        patch.object(module, "evaluate_caller_grid", return_value=result),
        patch.object(module, "target_study_document", return_value={}),
        patch.object(module, "target_runs_document", return_value={}),
        patch.object(module, "role_source_document", return_value={}),
        patch.object(module, "caller_eligible_roster_document", return_value=[]),
        patch.object(module, "caller_role_evidence_document", return_value={}),
        patch.object(module, "caller_evaluation_document", return_value={}),
        patch.object(module, "render_caller_evaluation", return_value="failed report"),
    ):
        completed = module.fit_caller_bundle(args, output)

    assert completed is False
    assert (output / "report.html").read_text() == "failed report"
    assert (output / "checksums.json").is_file()
    assert not (output / "candidate.json").exists()


def test_public_commands_require_path_arguments(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    with pytest.raises(ValueError, match="Path evidence"):
        module.fit_caller_bundle(Namespace(evidence="bad", exposure_ledger=tmp_path, objective="caller-safety-v1"), tmp_path)
    with pytest.raises(ValueError, match="Path profile"):
        module.assess_caller_bundle(Namespace(profile="bad", intake=tmp_path, runs=tmp_path, exposure_ledger=tmp_path), tmp_path)


def test_profile_loader_rejects_incomplete_and_symlink_directories(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_profile_io")
    incomplete = tmp_path / "incomplete"
    incomplete.mkdir()
    with pytest.raises(ValueError, match="inventory"):
        module.load_caller_research_profile(incomplete)
    linked = tmp_path / "linked"
    linked.symlink_to(incomplete, target_is_directory=True)
    with pytest.raises(ValueError, match="nonsymlink"):
        module.load_caller_research_profile(linked)


def test_supplied_locked_truth_bytes_are_hash_checked_and_never_reopened(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    raw = canonical_json_bytes(
        {
            "schema_version": "calibration-caller-truth-v1",
            "rows": [{"key": "artifact-1", "genotype": "negative", "variants": []}],
        }
    )
    source = type(
        "Source",
        (),
        {
            "keys": ("artifact-1",),
            "truth_asset": TargetRunAsset(tmp_path / "sealed.json", hashlib.sha256(raw).hexdigest(), len(raw)),
        },
    )()
    with patch.object(module, "read_target_json", side_effect=AssertionError("truth reopened")):
        truth = module._opened_truth(source, raw)
    assert truth.by_key["artifact-1"].genotype == "negative"
    with pytest.raises(ValueError, match="sealed source"):
        module._opened_truth(source, raw + b" ")


def test_fixed_evaluation_rejects_native_assets_outside_selected_bundle(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_caller_controller")
    study = _caller_study()
    assert isinstance(study.protocol, CallerProtocol)
    policy = study.protocol.candidates[0].policy
    background_raw = canonical_json_bytes(project_portable_background(background()))
    background_path = tmp_path / "background.json"
    background_path.write_bytes(background_raw)
    capture_path = tmp_path / "capture.jsonl"
    capture_path.write_bytes(b"capture")
    run = TargetRun(
        "artifact-1", policy.sha256, "1" * 64, "2" * 64, study.baseline.assets_sha256, "3" * 64,
        "scalar-replay", 0,
        MappingProxyType({
            "advntr_model": TargetRunAsset(tmp_path / "model", "4" * 64, 1),
            "advntr_background": TargetRunAsset(
                background_path, hashlib.sha256(background_raw).hexdigest(), len(background_raw)
            ),
            "advntr_capture": TargetRunAsset(
                capture_path, hashlib.sha256(b"capture").hexdigest(), len(b"capture")
            ),
        }),
        (17,),
    )
    with patch("vntyper.scripts.calibration_caller_profile_io._capture_records", return_value=({"capture_policy": capture_policy(), "producer": {"package_version": "2.4.1", "build_id": "b" * 64, "source_revision": "c" * 40}, "assets": {"model_sha256": "4" * 64}},)):
        fitted = type(
            "Background", (),
            {"background_bytes": background_raw, "background_sha256": hashlib.sha256(background_raw).hexdigest()},
        )()
        _, _, files, _ = module._payload(
            study, policy, type("Source", (), {"sha256": "5" * 64})(), (run,), fitted
        )
    profile = type(
        "Profile", (),
        {"study": study, "selected_protocol_candidate_id": policy.sha256, "payload_files": files},
    )()
    changed = replace(run, assets=MappingProxyType({**run.assets, "advntr_model": replace(run.assets["advntr_model"], sha256="f" * 64)}))
    with (
        patch.object(module, "select_target_runs", return_value=(changed,)),
        pytest.raises(ValueError, match="run assets"),
    ):
        module._require_profile_run_bindings(profile, object(), ("artifact-1",))
