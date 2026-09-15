"""Confirmation evaluates only after durable exposure and one-use custody claims."""

import hashlib
from argparse import Namespace
from importlib import import_module

import pytest

from tests.unit.test_calibration_candidate import candidate_document, resign
from tests.unit.test_calibration_target_contract import study_document
from tests.unit.test_calibration_target_runs import run_document
from vntyper.scripts.calibration_artifact_io import load_object, write_json
from vntyper.scripts.calibration_candidate import decode_candidate
from vntyper.scripts.calibration_exposure_io import initialize_exposure_ledger
from vntyper.scripts.calibration_manifest import connected_leakage_groups
from vntyper.scripts.calibration_target_authority import encode_target_custodian_authority
from vntyper.scripts.calibration_target_contract import decode_target_study
from vntyper.scripts.calibration_target_custody import initialize_target_custody
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256

pytestmark = pytest.mark.unit


def setup_case(tmp_path, monkeypatch, *, passed=True):
    module = import_module("vntyper.scripts.calibration_confirmation_controller")
    ledger = tmp_path / "exposure.jsonl"
    ledger_id = initialize_exposure_ledger(ledger)
    raw_study = study_document()
    raw_study["exposure_ledger_id"] = ledger_id
    study = decode_target_study(raw_study)
    raw_candidate = candidate_document()
    raw_candidate.update(study_sha256=study.sha256, partition_sha256=study.partitions.sha256)
    candidate = decode_candidate(resign(raw_candidate))
    custody = tmp_path / "custody"
    initialize_target_custody(custody, exposure_ledger_id=ledger_id)
    calls = []

    def evaluate(runs, source, receipt, opened_truth):
        assert receipt.evidence_sha256 == source.sha256
        assert ledger.stat().st_size > 0
        assert list(custody.glob("*." + ("validation" if source.role == "validation" else "locked") + "-started.json"))
        if source.role == "locked-heldout":
            assert list(custody.glob("*.locked-consumption.json"))
            assert opened_truth == b'{"invented":true}'
        else:
            assert opened_truth is None
        calls.append(source.role)
        return module.ConfirmationOutcome({"scientific_passed": passed}, "<html>Invented result</html>", passed)

    backend = module.ConfirmationBackend(candidate, study, evaluate)
    monkeypatch.setattr(module, "load_confirmation_backend", lambda *_: backend)
    profile = tmp_path / "profile"
    profile.mkdir()
    return module, Namespace(target="length", profile=profile, exposure_ledger=ledger, custody=custody), backend, calls


def source_dir(tmp_path, backend, role):
    index = 2 if role == "validation" else 3
    key = f"artifact-{index}"
    root = tmp_path / role
    root.mkdir()
    truth_raw = b'{"invented":true}'
    truth_path = root / "truth.json"
    truth_path.write_bytes(truth_raw)
    runs = run_document()
    runs["runs"][0].update(manifest_key=key, input_sha256=str(index) * 64)
    write_json(root / "runs.json", runs)
    raw = {
        "schema_version": "calibration-role-source-v2",
        "study_sha256": backend.study.sha256,
        "role": role,
        "roster": [
            {"key": key, "group_key": connected_leakage_groups(backend.study.partitions)[key], "strata": ["all"]}
        ],
        "excluded": [],
        "evidence_domains": {key: "synthetic"},
        "identities_by_key": {
            key: [
                {"namespace": "physical-readset", "sha256": str(index) * 64},
                {"namespace": "specimen", "sha256": str(index + 2) * 64},
            ]
        },
        "truth_asset": {
            "path": str(truth_path),
            "sha256": hashlib.sha256(truth_raw).hexdigest(),
            "size_bytes": len(truth_raw),
        },
        "run_manifest_sha256": canonical_sha256(runs),
    }
    write_json(root / "source.json", raw)
    return root, raw


@pytest.mark.parametrize("passed", [True, False])
def test_validation_publishes_actual_scientific_outcome_and_cannot_retry(tmp_path, monkeypatch, passed):
    module, args, backend, calls = setup_case(tmp_path, monkeypatch, passed=passed)
    args.evidence, _ = source_dir(tmp_path, backend, "validation")
    output = tmp_path / "result"
    output.mkdir()
    assert module.confirm_calibration_bundle(args, output, role="validation") is passed
    attestation = load_object(output / "validation-attestation.json", "test")
    assert attestation["status"] == ("passed" if passed else "failed")
    assert attestation["metrics_sha256"] == canonical_sha256({"scientific_passed": passed})
    assert calls == ["validation"]
    assert not (output / "completion.json").exists()
    with pytest.raises(ValueError):
        module.confirm_calibration_bundle(args, output, role="validation")
    assert calls == ["validation"]


def test_wrong_source_refused_before_exposure_or_evaluation(tmp_path, monkeypatch):
    module, args, backend, calls = setup_case(tmp_path, monkeypatch)
    args.evidence, raw = source_dir(tmp_path, backend, "validation")
    before = args.exposure_ledger.read_bytes()
    raw["study_sha256"] = "0" * 64
    write_json(args.evidence / "source.json", raw)
    with pytest.raises(ValueError):
        module.confirm_calibration_bundle(args, tmp_path / "output", role="validation")
    assert args.exposure_ledger.read_bytes() == before
    assert not calls


@pytest.mark.parametrize("locked_passed", [True, False])
def test_locked_reads_once_after_consumption_and_completes(tmp_path, monkeypatch, locked_passed):
    from dataclasses import replace

    module, args, backend, calls = setup_case(tmp_path, monkeypatch)
    args.evidence, _ = source_dir(tmp_path, backend, "validation")
    validation_output = tmp_path / "validation-output"
    validation_output.mkdir()
    assert module.confirm_calibration_bundle(args, validation_output, role="validation")
    args.validation = validation_output / "validation-attestation.json"
    validation = load_object(args.validation, "test")
    args.evidence, source = source_dir(tmp_path, backend, "locked-heldout")
    args.authority = tmp_path / "authority.json"
    common = {
        key: validation[key]
        for key in (
            "target",
            "candidate_sha256",
            "candidate_id",
            "study_sha256",
            "protocol_sha256",
            "partition_sha256",
            "baseline_sha256",
            "exposure_ledger_id",
        )
    }
    authority = encode_target_custodian_authority(
        **common,
        validation_attestation_sha256=canonical_sha256(validation),
        validation_evidence_sha256=validation["evidence_sha256"],
        validation_run_manifest_sha256=validation["run_manifest_sha256"],
        validation_exposure_receipt_sha256=validation["exposure_receipt_sha256"],
        locked_heldout_evidence_sha256=canonical_sha256(source),
        locked_heldout_run_manifest_sha256=source["run_manifest_sha256"],
        locked_payload_sha256=source["truth_asset"]["sha256"],
        custodian_name="invented independent operator",
        attestation_id="invented approval",
    )
    write_json(args.authority, authority)
    original_evaluate = backend.evaluate

    def locked_evaluate(*values):
        result = original_evaluate(*values)
        return module.ConfirmationOutcome({"scientific_passed": locked_passed}, result.report, locked_passed)

    backend = replace(backend, evaluate=locked_evaluate)
    monkeypatch.setattr(module, "load_confirmation_backend", lambda *_: backend)
    output = tmp_path / "locked-output"
    output.mkdir()
    assert module.confirm_calibration_bundle(args, output, role="locked-heldout") is locked_passed
    assert calls == ["validation", "locked-heldout"]
    locked = load_object(output / "locked-attestation.json", "test")
    assert locked["status"] == ("passed" if locked_passed else "failed")
    if locked_passed:
        terminal = load_object(output / "completion.json", "test")
        assert terminal["locked_heldout_attestation_sha256"] == canonical_sha256(locked)
        assert terminal["custody_consumption_receipt_sha256"] == locked["custody_consumption_receipt_sha256"]
    else:
        assert not (output / "completion.json").exists()
    assert (output / "report.html").read_text() == "<html>Invented result</html>"
    assert canonical_json_bytes(load_object(output / "metrics.json", "test")) == (output / "metrics.json").read_bytes()


def test_confirmation_uses_real_fixed_length_model_and_its_training_baseline(tmp_path, monkeypatch):
    from tests.unit.test_calibration_length_controller import fit_fixture
    from vntyper.scripts.calibration_length_controller import fit_length_bundle
    from vntyper.scripts.calibration_role_source import decode_role_source
    from vntyper.scripts.calibration_target_runs import decode_target_runs

    module = import_module("vntyper.scripts.calibration_confirmation_controller")
    fit_args, study = fit_fixture(tmp_path)
    profile = tmp_path / "profile"
    profile.mkdir()
    assert fit_length_bundle(fit_args, profile)
    # Confirmation must not call the fitting entry point after the candidate freezes.
    monkeypatch.setattr(
        "vntyper.scripts.calibration_length_controller.fit_length_bundle", lambda *_: pytest.fail("refit")
    )
    backend = module.load_confirmation_backend("length", profile)
    evidence = tmp_path / "confirmation-inputs"
    evidence.mkdir()
    source = load_object(fit_args.evidence / "roles" / "validation" / "source.json", "test")
    write_json(evidence / "source.json", source)
    runs = load_object(fit_args.evidence / "runs.json", "test")
    write_json(evidence / "runs.json", runs)
    custody = tmp_path / "custody"
    initialize_target_custody(custody, exposure_ledger_id=study.exposure_ledger_id)
    args = Namespace(
        target="length", profile=profile, evidence=evidence, exposure_ledger=fit_args.exposure_ledger, custody=custody
    )
    output = tmp_path / "validation-output"
    output.mkdir()
    assert module.confirm_calibration_bundle(args, output, role="validation")
    result = load_object(output / "metrics.json", "test")
    assert result["phase"] == "validation"
    assert result["baseline_sha256"] == backend.candidate.baseline_sha256
    assert result["candidates"][0]["candidate_id"] == "affine-a"
    assert result["candidates"][0]["acceptance"]["status"] == "passed"
    assert len(result["candidates"]) == 1
    assert decode_role_source(source, study=study, runs=decode_target_runs(runs), expected_role="validation").keys == (
        "v1",
        "v2",
    )
    predictions = result["candidates"][0]["predictions"]
    assert sorted(row["prediction"] for row in predictions) == pytest.approx([60, 160])
    assert [row["baseline_prediction"] for row in predictions] == [110, 110]


def test_evaluator_failure_retires_claim_and_cannot_be_retried(tmp_path, monkeypatch):
    from dataclasses import replace

    module, args, backend, _ = setup_case(tmp_path, monkeypatch)
    args.evidence, _ = source_dir(tmp_path, backend, "validation")

    def broken(*_):
        raise ValueError("invented unreadable native result")

    backend = replace(backend, evaluate=broken)
    monkeypatch.setattr(module, "load_confirmation_backend", lambda *_: backend)
    output = tmp_path / "output"
    output.mkdir()
    with pytest.raises(ValueError, match="unreadable"):
        module.confirm_calibration_bundle(args, output, role="validation")
    assert list(args.custody.glob("*.retired.json"))
    assert not (output / "validation-attestation.json").exists()
    with pytest.raises(ValueError):
        module.confirm_calibration_bundle(args, output, role="validation")


@pytest.mark.parametrize("role", ["training", "development-assessment", "typo"])
def test_confirmation_rejects_non_confirmation_roles_before_profile_read(role, tmp_path, monkeypatch):
    module = import_module("vntyper.scripts.calibration_confirmation_controller")
    monkeypatch.setattr(module, "load_confirmation_backend", lambda *_: pytest.fail("profile read"))
    with pytest.raises(ValueError, match="role"):
        module.confirm_calibration_bundle(Namespace(), tmp_path, role=role)


@pytest.mark.parametrize("kind", ["missing-report", "nonboolean-pass", "missing-outcome"])
def test_incomplete_evaluation_cannot_create_attestation(tmp_path, monkeypatch, kind):
    from dataclasses import replace

    module, args, backend, _ = setup_case(tmp_path, monkeypatch)
    args.evidence, _ = source_dir(tmp_path, backend, "validation")

    def invalid(*_):
        if kind == "missing-outcome":
            return None
        return module.ConfirmationOutcome(
            {}, "" if kind == "missing-report" else "report", True if kind == "missing-report" else 1
        )

    backend = replace(backend, evaluate=invalid)
    monkeypatch.setattr(module, "load_confirmation_backend", lambda *_: backend)
    output = tmp_path / "output"
    output.mkdir()
    with pytest.raises(ValueError, match="incomplete scientific outcome"):
        module.confirm_calibration_bundle(args, output, role="validation")
    assert not (output / "validation-attestation.json").exists()
    assert list(args.custody.glob("*.retired.json"))


def test_unknown_target_cannot_open_a_research_profile(tmp_path):
    module = import_module("vntyper.scripts.calibration_confirmation_controller")
    with pytest.raises(ValueError, match="target"):
        module.load_confirmation_backend("unknown", tmp_path)
