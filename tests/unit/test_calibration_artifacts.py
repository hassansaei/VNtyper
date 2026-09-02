"""Filesystem codecs and end-to-end calibration operation artifacts."""

import fcntl
import hashlib
import logging
import os
import shutil
import threading
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.calibration_run_fixture import IDENTITY, OTHER_IDENTITY, write_schema_three_run
from vntyper import cli
from vntyper.scripts.calibration_artifact_io import thaw_json, verify_checksums, write_checksums, write_json
from vntyper.scripts.calibration_artifacts import (
    _observations,
    evaluate_artifact_bundle,
    extract_artifact_bundle,
    fit_artifact_bundle,
    validate_artifact_bundle,
)
from vntyper.scripts.calibration_contract import decode_attestation, decode_metrics
from vntyper.scripts.calibration_custody import claim_candidate, retire_candidate
from vntyper.scripts.calibration_features import decode_feature_artifact, decode_label_artifact
from vntyper.scripts.calibration_locked_artifacts import decode_locked_payload
from vntyper.scripts.calibration_manifest import decode_study_declaration
from vntyper.scripts.calibration_objective import calculate_metrics, select_candidate
from vntyper.scripts.calibration_profiles import build_generated_profile
from vntyper.scripts.calibration_role_inputs import build_role_inputs
from vntyper.scripts.calibration_run_extraction import decode_run_hashes
from vntyper.scripts.calibration_statistics import BootstrapInterval
from vntyper.scripts.calibration_study_binding import decode_study_binding
from vntyper.scripts.calibration_validation_attestation import decode_validation_attestation
from vntyper.scripts.calibration_workflow import extract_evidence
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.cli_calibrate import handle_calibrate
from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.decision_profile import resolve_decision_profile

pytestmark = pytest.mark.unit


def test_extract_fit_validate_and_one_use_evaluate_round_trip(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    _run_cli(
        [
            "calibrate",
            "extract",
            "--truth",
            str(truth),
            "--partitions",
            str(partitions),
            "--runs",
            str(runs),
            "--output",
            str(evidence),
        ]
    )

    assert (evidence / "study.json").is_file()
    assert (evidence / "groups.json").is_file()
    assert (evidence / "roles" / "policy-selection" / "features.json").is_file()
    locked_root = evidence / "roles" / "locked-heldout"
    assert {path.name for path in locked_root.iterdir()} == {
        "checksums.json",
        "member_declaration.json",
        "run_commitments.json",
    }
    assert load_strict_json_object((evidence / "checksums.json").read_bytes())["schema_version"] == (
        "calibration-checksums-v1"
    )
    declared_runs = load_strict_json_object(runs.read_bytes())["runs"]
    training_hashes = load_strict_json_object((evidence / "roles" / "training" / "run_hashes.json").read_bytes())
    assert training_hashes["run_hashes"]["train"] == declared_runs["train"]["artifacts"]

    candidate = tmp_path / "candidate"
    _run_cli(
        [
            "calibrate",
            "fit",
            "--evidence",
            str(evidence),
            "--objective",
            "lexicographic-safety-v1",
            "--output",
            str(candidate),
        ]
    )
    profile = candidate / "decision_profile.json"
    assert profile.is_file()
    for name in (
        "report.html",
        "grid.json",
        "intervals.json",
        "roc.tsv",
        "pr.tsv",
        "joint_surface.tsv",
        "abstentions.tsv",
    ):
        assert (candidate / name).is_file(), name
    assert "optional minimum k-mer depth" in (candidate / "report.html").read_text(encoding="utf-8")
    assert load_strict_json_object((candidate / "fit_attestation.json").read_bytes())["accessed_roles"] == [
        "training",
        "policy-selection",
    ]

    validation = tmp_path / "validation"
    _run_cli(
        [
            "calibrate",
            "validate",
            "--profile",
            str(profile),
            "--evidence",
            str(evidence),
            "--output",
            str(validation),
        ]
    )
    validation_attestation = decode_validation_attestation(
        load_strict_json_object((validation / "attestation.json").read_bytes())
    )
    assert validation_attestation.role == "validation"
    assert (validation / "report.html").is_file()
    assert (validation / "intervals.json").is_file()

    with pytest.raises(ValueError, match="custodian.*import|authority attestation"):
        evaluate_artifact_bundle(profile, evidence, tmp_path / "heldout")

    custodian_import = _custodian_import(tmp_path / "custodian", truth, partitions, runs, profile, validation)
    heldout = tmp_path / "heldout"
    assert evaluate_artifact_bundle(profile, custodian_import, heldout) is True
    assert decode_attestation(load_strict_json_object((heldout / "attestation.json").read_bytes())).role == (
        "locked-heldout"
    )
    assert (heldout / "custody_receipt.json").is_file()
    limitations = load_strict_json_object((heldout / "custody_limitations.json").read_bytes())
    assert limitations["local_verification_proves_external_independence"] is False

    with pytest.raises(ValueError, match="consumed|one-use|precommit"):
        evaluate_artifact_bundle(profile, custodian_import, tmp_path / "heldout-again")


def test_ordinary_extract_discloses_no_locked_values_or_custody_assertion(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()

    extract_artifact_bundle(truth, partitions, runs, evidence)

    locked_root = evidence / "roles" / "locked-heldout"
    assert {path.name for path in locked_root.iterdir()} == {
        "checksums.json",
        "member_declaration.json",
        "run_commitments.json",
    }
    locked_bytes = b"".join(path.read_bytes() for path in sorted(locked_root.iterdir()) if path.is_file())
    assert b"external-custodian" not in locked_bytes
    assert IDENTITY.encode() not in locked_bytes
    assert b"59dupC" not in locked_bytes
    for forbidden in ("features.json", "labels.json", "baseline.json", "locked_payload.json", "evidence_manifest.json"):
        assert not (locked_root / forbidden).exists()


def test_ordinary_extract_accepts_nonlocked_truth_without_reading_locked_run_roots(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    runs_document = load_strict_json_object(runs.read_bytes())
    run_values = runs_document["runs"]
    assert isinstance(run_values, dict) and isinstance(run_values["held"], dict)
    locked_root = Path(str(run_values["held"]["root"]))
    original_read = Path.read_bytes

    def refuse_locked_reads(path: Path) -> bytes:
        if path == locked_root or locked_root in path.parents:
            raise AssertionError(f"ordinary extraction opened locked run bytes: {path}")
        return original_read(path)

    with patch.object(Path, "read_bytes", refuse_locked_reads):
        assert extract_artifact_bundle(truth, partitions, runs, tmp_path / "evidence") is True


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("locked-row", "label rows|partition manifest"),
        ("unknown-root", "label.*fields|closed contract"),
        ("wrong-schema", "label.*schema"),
        ("unknown-row", "label rows|partition manifest"),
        ("duplicate-row", "label rows|partition manifest|unique"),
        ("out-of-order", "label rows|partition manifest|align|increasing"),
    ],
)
def test_ordinary_extract_requires_exact_canonical_nonlocked_truth(tmp_path: Path, mutation: str, message: str) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    document = load_strict_json_object(truth.read_bytes())
    labels = document["labels"]
    assert isinstance(labels, dict) and isinstance(labels["rows"], list)
    rows = labels["rows"]
    if mutation == "locked-row":
        rows.insert(0, _label_row("held"))
    elif mutation == "unknown-root":
        labels["attacker_authored"] = True
    elif mutation == "wrong-schema":
        labels["schema_version"] = "calibration-labels-v0"
    elif mutation == "unknown-row":
        rows.append(_label_row("unknown"))
    elif mutation == "duplicate-row":
        rows.insert(1, dict(rows[0]))
    else:
        rows[0], rows[1] = rows[1], rows[0]
    truth.write_bytes(canonical_json_bytes(document))

    with pytest.raises(ValueError, match=message):
        extract_artifact_bundle(truth, partitions, runs, tmp_path / "evidence")


def test_fit_refuses_an_objective_different_from_the_snapshotted_protocol(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)

    with pytest.raises(ValueError, match="objective"):
        fit_artifact_bundle(evidence, "f1", tmp_path / "candidate")


def test_fit_verifies_only_authorized_role_payloads(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    locked_payload = evidence / "roles" / "locked-heldout" / "locked_payload.json"
    original_read = Path.read_bytes

    def guarded_read(path: Path) -> bytes:
        if path == locked_payload:
            raise AssertionError("fit opened locked held-out evidence")
        return original_read(path)

    candidate = tmp_path / "candidate"
    candidate.mkdir()
    with patch.object(Path, "read_bytes", guarded_read):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    assert (candidate / "decision_profile.json").is_file()


def test_fit_refuses_tampered_policy_selection_artifacts(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    features = evidence / "roles" / "policy-selection" / "features.json"
    features.write_bytes(features.read_bytes() + b" ")

    with pytest.raises(ValueError, match="checksum"):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", tmp_path / "candidate")


def test_validation_refuses_a_self_checksummed_role_directory_swap(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation_role = evidence / "roles" / "validation"
    shutil.rmtree(validation_role)
    shutil.copytree(evidence / "roles" / "policy-selection", validation_role)

    with pytest.raises(ValueError, match="role bundle|validation"):
        validate_artifact_bundle(
            candidate / "decision_profile.json",
            evidence,
            tmp_path / "validation-attestation",
        )


def test_completed_failed_validation_writes_attestation_and_retirement_before_returning_false(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"

    with patch("vntyper.scripts.calibration_artifacts.select_candidate", return_value=None):
        assert validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation) is False

    attestation = load_strict_json_object((validation / "attestation.json").read_bytes())
    assert attestation["status"] == "failed"
    assert (validation / "retirement.json").is_file()
    with pytest.raises(ValueError, match="retired"):
        validate_artifact_bundle(candidate / "decision_profile.json", evidence, tmp_path / "validation-again")


def test_completed_failed_validation_retires_before_fallible_result_writes(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)

    with (
        patch("vntyper.scripts.calibration_artifacts.select_candidate", return_value=None),
        patch("vntyper.scripts.calibration_artifacts.write_json", side_effect=RuntimeError("disk failed")),
        pytest.raises(RuntimeError, match="disk failed"),
    ):
        validate_artifact_bundle(candidate / "decision_profile.json", evidence, tmp_path / "validation")

    custody = tmp_path / ".evidence.calibration-custody"
    assert len(tuple((custody / "retired").glob("*.json"))) == 1


def test_completed_failed_locked_evaluation_writes_attestation_and_retirement(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation)
    custodian_import = _custodian_import(
        tmp_path / "custodian",
        truth,
        partitions,
        runs,
        candidate / "decision_profile.json",
        validation,
    )
    heldout = tmp_path / "heldout"

    with patch("vntyper.scripts.calibration_artifacts.select_candidate", return_value=None):
        assert evaluate_artifact_bundle(candidate / "decision_profile.json", custodian_import, heldout) is False

    assert load_strict_json_object((heldout / "attestation.json").read_bytes())["status"] == "failed"
    assert (heldout / "retirement.json").is_file()


@pytest.mark.parametrize(("passed", "status"), [(True, "passed"), (False, "failed")])
def test_terminal_unlock_error_preserves_coherent_cli_outcome(tmp_path: Path, passed: bool, status: str) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    profile = candidate / "decision_profile.json"
    validation = tmp_path / "validation"
    validate_artifact_bundle(profile, evidence, validation)
    custodian = _custodian_import(tmp_path / "custodian", truth, partitions, runs, profile, validation)
    output = tmp_path / "heldout"
    real_flock = fcntl.flock

    def fail_unlock(descriptor: int, operation: int) -> None:
        if operation == fcntl.LOCK_UN:
            raise OSError("explicit unlock failed")
        real_flock(descriptor, operation)

    with (
        patch("vntyper.cli.setup_logging"),
        patch("vntyper.scripts.calibration_custody.fcntl.flock", side_effect=fail_unlock),
        patch(
            "vntyper.scripts.calibration_artifacts.select_candidate",
            side_effect=select_candidate if passed else lambda *_args: None,
        ),
    ):
        if passed:
            assert (
                cli.main(
                    [
                        "calibrate",
                        "evaluate",
                        "--profile",
                        str(profile),
                        "--evidence",
                        str(custodian),
                        "--output",
                        str(output),
                    ]
                )
                is None
            )
        else:
            with pytest.raises(SystemExit) as completed_failure:
                cli.main(
                    [
                        "calibrate",
                        "evaluate",
                        "--profile",
                        str(profile),
                        "--evidence",
                        str(custodian),
                        "--output",
                        str(output),
                    ]
                )
            assert completed_failure.value.code == 1

    verify_checksums(output)
    assert load_strict_json_object((output / "attestation.json").read_bytes())["status"] == status
    custody = tmp_path / ".custodian.calibration-custody"
    assert len(tuple((custody / ("completed" if passed else "retired")).glob("*.json"))) == 1
    assert not tuple((custody / ("retired" if passed else "completed")).glob("*.json"))


def test_post_consumption_interrupt_retires_before_reraising(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation)
    custodian = _custodian_import(
        tmp_path / "custodian", truth, partitions, runs, candidate / "decision_profile.json", validation
    )

    with (
        patch("vntyper.scripts.calibration_artifacts.decode_metrics", side_effect=KeyboardInterrupt()),
        pytest.raises(KeyboardInterrupt),
    ):
        evaluate_artifact_bundle(candidate / "decision_profile.json", custodian, tmp_path / "heldout")

    custody = tmp_path / ".custodian.calibration-custody"
    assert len(tuple((custody / "consumed").glob("*.json"))) == 1
    assert len(tuple((custody / "retired").glob("*.json"))) == 1


def test_claim_owner_excludes_retirement_through_locked_result_finalization(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    profile_path = candidate / "decision_profile.json"
    validation = tmp_path / "validation"
    validate_artifact_bundle(profile_path, evidence, validation)
    custodian = _custodian_import(tmp_path / "custodian", truth, partitions, runs, profile_path, validation)
    authority = load_strict_json_object((custodian / "authority_attestation.json").read_bytes())
    profile = resolve_decision_profile(profile_path)
    custody = tmp_path / ".custodian.calibration-custody"
    decoder_entered = threading.Event()
    release_decoder = threading.Event()
    evaluation_results: list[bool] = []
    evaluation_errors: list[BaseException] = []
    retirement_errors: list[BaseException] = []
    real_decoder = decode_metrics

    def blocking_decoder(value: object):
        decoder_entered.set()
        assert release_decoder.wait(timeout=5)
        return real_decoder(value)

    def evaluate() -> None:
        try:
            evaluation_results.append(evaluate_artifact_bundle(profile_path, custodian, tmp_path / "heldout"))
        except ValueError as error:
            evaluation_errors.append(error)

    def retire() -> None:
        try:
            retire_candidate(
                custody,
                profile.digest,
                str(authority["locked_payload_sha256"]),
                "concurrent-retirement",
            )
        except ValueError as error:
            retirement_errors.append(error)

    with patch("vntyper.scripts.calibration_artifacts.decode_metrics", side_effect=blocking_decoder):
        evaluation_thread = threading.Thread(target=evaluate)
        evaluation_thread.start()
        assert decoder_entered.wait(timeout=5)
        retirement_thread = threading.Thread(target=retire)
        retirement_thread.start()
        retirement_thread.join(timeout=0.2)
        assert retirement_thread.is_alive()
        release_decoder.set()
        evaluation_thread.join(timeout=5)
        retirement_thread.join(timeout=5)

    assert evaluation_results == [True] and not evaluation_errors
    assert len(retirement_errors) == 1 and "completed" in str(retirement_errors[0])
    assert len(tuple((custody / "completed").glob("*.json"))) == 1
    assert not tuple((custody / "retired").glob("*.json"))


@pytest.mark.parametrize("mutation", ["aggregate", "row-key", "row-order"])
def test_locked_evaluation_rejects_unreplayed_baseline_semantics(tmp_path: Path, mutation: str) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation)
    custodian = _custodian_import(
        tmp_path / "custodian", truth, partitions, runs, candidate / "decision_profile.json", validation
    )
    payload_path = custodian / "locked_payload.json"
    payload = load_strict_json_object(payload_path.read_bytes())
    baseline = payload["baseline"]
    assert isinstance(baseline, dict)
    expected = baseline["expected"]
    observed = baseline["observed"]
    assert isinstance(expected, dict) and isinstance(observed, dict)
    if mutation == "aggregate":
        assert isinstance(observed["aggregate"], dict)
        observed["aggregate"]["displayed"] = 999
    else:
        for projection in (expected, observed):
            rows = projection["rows"]
            assert isinstance(rows, list) and isinstance(rows[0], dict)
            rows[0]["manifest_key" if mutation == "row-key" else "order"] = "train" if mutation == "row-key" else 999
    payload_path.write_bytes(canonical_json_bytes(payload))
    _refresh_locked_authority(custodian)

    with pytest.raises(ValueError, match="baseline|row|aggregate|order"):
        evaluate_artifact_bundle(candidate / "decision_profile.json", custodian, tmp_path / "heldout")


def test_locked_payload_symlink_swap_after_precommit_is_rejected_and_retired(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation)
    custodian = _custodian_import(
        tmp_path / "custodian", truth, partitions, runs, candidate / "decision_profile.json", validation
    )
    payload = custodian / "locked_payload.json"
    replacement = tmp_path / "replacement-locked-payload.json"

    def swap_after_precommit(*args, **kwargs):
        precommit = claim_candidate(*args, **kwargs)
        replacement.write_bytes(payload.read_bytes())
        payload.unlink()
        os.symlink(replacement, payload)
        return precommit

    with (
        patch("vntyper.scripts.calibration_workflow.claim_candidate", side_effect=swap_after_precommit),
        pytest.raises(ValueError, match="symlink|regular|locked evidence"),
    ):
        evaluate_artifact_bundle(candidate / "decision_profile.json", custodian, tmp_path / "heldout")

    assert len(tuple((tmp_path / ".custodian.calibration-custody" / "retired").glob("*.json"))) == 1


def test_locked_evaluation_refuses_manifest_protocol_mismatched_to_payload(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation)
    custodian = _custodian_import(
        tmp_path / "custodian",
        truth,
        partitions,
        runs,
        candidate / "decision_profile.json",
        validation,
    )
    authority_path = custodian / "authority_attestation.json"
    authority = load_strict_json_object(authority_path.read_bytes())
    authority["protocol_sha256"] = "a" * 64
    authority_path.write_bytes(canonical_json_bytes(authority))
    write_checksums(custodian)

    with pytest.raises(ValueError, match="protocol|study bindings"):
        evaluate_artifact_bundle(
            candidate / "decision_profile.json",
            custodian,
            tmp_path / "heldout-attestation",
        )
    custody = tmp_path / ".custodian.calibration-custody"
    assert not tuple((custody / "precommits").glob("*.json"))
    assert not tuple((custody / "consumed").glob("*.json"))
    assert not tuple((custody / "retired").glob("*.json"))


def test_locked_evaluation_rejects_cross_study_reuse_with_identical_candidate_grid(tmp_path: Path) -> None:
    truth_a, partitions_a, runs_a = _inputs(tmp_path / "a")
    evidence_a = tmp_path / "evidence-a"
    evidence_a.mkdir()
    extract_artifact_bundle(truth_a, partitions_a, runs_a, evidence_a)
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fit_artifact_bundle(evidence_a, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence_a, validation)

    truth_b, partitions_b, runs_b = _inputs(tmp_path / "b")
    study_b = load_strict_json_object(partitions_b.read_bytes())
    protocol_b = study_b["protocol"]
    assert isinstance(protocol_b, dict)
    protocol_b["seed"] = 296
    partitions_b.write_bytes(canonical_json_bytes(study_b))
    custodian_b = _custodian_import(
        tmp_path / "custodian-b",
        truth_b,
        partitions_b,
        runs_b,
        candidate / "decision_profile.json",
        validation,
    )

    with pytest.raises(ValueError, match="validation|study|protocol|bindings"):
        evaluate_artifact_bundle(candidate / "decision_profile.json", custodian_b, tmp_path / "heldout")


def test_generated_profile_rejects_alternate_validation_run_commitments(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence_a = tmp_path / "evidence-a"
    extract_artifact_bundle(truth, partitions, runs, evidence_a)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence_a, "lexicographic-safety-v1", candidate)
    runs_document = load_strict_json_object(runs.read_bytes())
    run_values = runs_document["runs"]
    assert isinstance(run_values, dict)
    run_values["validate"] = write_schema_three_run(tmp_path / "alternate-validation", "validate", support=11)
    runs.write_bytes(canonical_json_bytes(runs_document))
    evidence_b = tmp_path / "evidence-b"
    extract_artifact_bundle(truth, partitions, runs, evidence_b)

    with pytest.raises(ValueError, match="study|dataset|binding|commitment"):
        validate_artifact_bundle(candidate / "decision_profile.json", evidence_b, tmp_path / "validation")


def test_locked_evaluation_rejects_alternate_held_run_commitments(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation)
    runs_document = load_strict_json_object(runs.read_bytes())
    run_values = runs_document["runs"]
    assert isinstance(run_values, dict)
    run_values["held"] = write_schema_three_run(tmp_path / "alternate-held", "held", support=11)
    runs.write_bytes(canonical_json_bytes(runs_document))
    custodian = _custodian_import(
        tmp_path / "custodian", truth, partitions, runs, candidate / "decision_profile.json", validation
    )

    with pytest.raises(ValueError, match="held|run|commitment|binding"):
        evaluate_artifact_bundle(candidate / "decision_profile.json", custodian, tmp_path / "heldout")


def test_locked_evaluation_rejects_arbitrary_validation_evidence_hash(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation = tmp_path / "validation"
    validate_artifact_bundle(candidate / "decision_profile.json", evidence, validation)
    custodian = _custodian_import(
        tmp_path / "custodian", truth, partitions, runs, candidate / "decision_profile.json", validation
    )
    validation_path = custodian / "validation_attestation.json"
    validation_document = load_strict_json_object(validation_path.read_bytes())
    validation_document["evidence_sha256"] = "9" * 64
    write_json(validation_path, validation_document)
    authority_path = custodian / "authority_attestation.json"
    authority = load_strict_json_object(authority_path.read_bytes())
    authority["validation_attestation_sha256"] = hashlib.sha256(validation_path.read_bytes()).hexdigest()
    write_json(authority_path, authority)
    write_checksums(custodian)

    with pytest.raises(ValueError, match="validation.*evidence|binding|commitment"):
        evaluate_artifact_bundle(candidate / "decision_profile.json", custodian, tmp_path / "heldout")


def test_fit_rejects_a_declared_but_unobserved_bootstrap_stratum(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    study = load_strict_json_object(partitions.read_bytes())
    protocol = study["protocol"]
    assert isinstance(protocol, dict)
    protocol["assay_classes"] = ["capture-short-read", "genome-short-read"]
    partitions.write_bytes(canonical_json_bytes(study))
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)

    with pytest.raises(ValueError, match="empty declared bootstrap strata.*genome-short-read:duplication"):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", tmp_path / "candidate")


def test_fit_does_not_infer_canonical_identity_correctness_from_display_name(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    truth_document = load_strict_json_object(truth.read_bytes())
    labels = truth_document["labels"]
    assert isinstance(labels, dict)
    rows = labels["rows"]
    assert isinstance(rows, list)
    for raw in rows:
        assert isinstance(raw, dict)
        raw["expected_identity"] = OTHER_IDENTITY
    truth.write_bytes(canonical_json_bytes(truth_document))
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)

    metrics = load_strict_json_object((candidate / "metrics.json").read_bytes())
    assert metrics["macro_exact_recovery"] == "0"
    assert metrics["wrong_displayed_names_all_tiers"] == 0


def test_candidate_a_to_b_is_counted_as_a_wrong_tier_a_displayed_name(tmp_path: Path) -> None:
    truth, partitions, runs_path = _inputs(tmp_path)
    runs_document = load_strict_json_object(runs_path.read_bytes())
    run_values = runs_document["runs"]
    assert isinstance(run_values, dict)
    for key in tuple(run_values):
        run_values[key] = write_schema_three_run(
            tmp_path / "wrong-selection-runs" / key,
            key,
            identity=OTHER_IDENTITY,
            name="58dupG",
            with_advntr=True,
            advntr_identity=IDENTITY,
            advntr_name="59dupC",
            selected_name="59dupC",
            reconciled_identity=IDENTITY,
        )
    study = decode_study_declaration(load_strict_json_object(partitions.read_bytes()))
    labels_document = load_strict_json_object(truth.read_bytes())["labels"]
    assert isinstance(labels_document, dict) and isinstance(labels_document["rows"], list)
    labels_document["rows"].insert(0, _label_row("held"))
    labels = decode_label_artifact(labels_document)
    evidence = extract_evidence(study, labels, run_values)

    profile = build_generated_profile(
        {
            "enabled": True,
            "minimum_record_count_margin": 1,
            "minimum_record_share": 0.5,
            "minimum_record_share_margin": 0.0,
            "xd_veto": "disabled",
            "abstain_on_inadmissible_advntr": False,
        },
        dataset_manifest_hash="b" * 64,
        partition_manifest_hash="c" * 64,
        seed=295,
        objective="lexicographic-safety-v1",
        generator_version="test",
    )
    observations, _ = _observations(profile, evidence)
    summary = calculate_metrics(
        observations,
        profile_sha256="a" * 64,
        free_parameter_count=0,
        required_strata=("capture-short-read:duplication",),
    )

    assert {row.displayed_name for row in observations} == {"58dupG"}
    assert {row.tier for row in observations} == {"A"}
    assert summary.metrics.wrong_displayed_names_all_tiers == 3
    assert summary.metrics.wrong_tier_a_displayed_names == 3


@pytest.mark.parametrize("forbidden", ["features", "baseline"])
def test_extract_truth_is_strictly_labels_only(tmp_path: Path, forbidden: str) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    document = load_strict_json_object(truth.read_bytes())
    document[forbidden] = {"attacker_authored": True}
    truth.write_bytes(canonical_json_bytes(document))

    with pytest.raises(ValueError, match="truth.*closed contract|labels.only"):
        extract_artifact_bundle(truth, partitions, runs, tmp_path / "evidence")


def test_candidate_p_value_is_the_conservative_maximum_of_both_required_endpoints(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    detection = BootstrapInterval(Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1, 100), 10_000, 1)
    exact = BootstrapInterval(Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1, 25), 10_000, 1)

    with patch("vntyper.scripts.calibration_artifacts.paired_group_bootstrap", side_effect=(detection, exact)):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)

    evaluation = load_strict_json_object((candidate / "evaluation.json").read_bytes())
    assert evaluation["holm_adjusted_p_value"] == "1/25"


def _inputs(root: Path) -> tuple[Path, Path, Path]:
    keys = ("held", "select", "train", "validate")
    roles = ("locked-heldout", "policy-selection", "training", "validation")
    labels = []
    members = []
    run_paths: dict[str, object] = {}
    for key, role in zip(keys, roles, strict=True):
        if role != "locked-heldout":
            labels.append(_label_row(key))
        members.append(
            {
                "key": key,
                "role": role,
                "provenance": "external-custodian" if role == "locked-heldout" else "development",
                "assay_class": "capture-short-read",
                "groups": {
                    namespace: [f"{namespace}:{key}"]
                    for namespace in (
                        "individual-family",
                        "simulated-pair",
                        "backbone-seed-lineage",
                        "replicate-rerun",
                        "depth-series-source",
                        "batch",
                        "repeat-context",
                    )
                },
            }
        )
        run = root / "runs" / key
        run_paths[key] = write_schema_three_run(run, key)

    truth_value = {
        "schema_version": "calibration-truth-v1",
        "labels": {"schema_version": "calibration-labels-v1", "rows": labels},
    }
    study_value = {
        "schema_version": "calibration-study-v1",
        "protocol": {
            "objective": "lexicographic-safety-v1",
            "bootstrap_iterations": 10000,
            "bootstrap_interval": "percentile",
            "multiplicity_method": "holm",
            "seed": 295,
            "maximum_free_parameters": 4,
            "minimum_stratum_count": 1,
            "maximum_abstention_fraction": 0.25,
            "assay_classes": ["capture-short-read"],
            "mutation_classes": ["duplication"],
            "candidate_grid": {
                "minimum_record_count_margin": [1],
                "minimum_record_share": [0.5],
                "minimum_record_share_margin": [0.0],
                "xd_veto": ["disabled"],
            },
        },
        "partitions": {"schema_version": "calibration-partitions-v1", "members": members},
    }
    runs_value = {"schema_version": "calibration-runs-v1", "runs": run_paths}
    truth_path = root / "truth.json"
    partitions_path = root / "partitions.json"
    runs_path = root / "runs.json"
    truth_path.write_bytes(canonical_json_bytes(truth_value))
    partitions_path.write_bytes(canonical_json_bytes(study_value))
    runs_path.write_bytes(canonical_json_bytes(runs_value))
    return truth_path, partitions_path, runs_path


def _label_row(key: str) -> dict[str, object]:
    return {
        "label_key": f"l-{key}",
        "manifest_key": key,
        "truth_status": "mutated",
        "expected_identity": IDENTITY,
        "expected_display_name": "59dupC",
        "mutation_class": "duplication",
    }


def _custodian_import(
    root: Path,
    truth_path: Path,
    study_path: Path,
    runs_path: Path,
    profile_path: Path,
    validation_root: Path,
) -> Path:
    """Create a separately supplied custodian fixture; production has no minting command."""
    root.mkdir()
    truth = load_strict_json_object(truth_path.read_bytes())
    study_raw = load_strict_json_object(study_path.read_bytes())
    runs_raw = load_strict_json_object(runs_path.read_bytes())
    study = decode_study_declaration(study_raw)
    ordinary_labels = truth["labels"]
    assert isinstance(ordinary_labels, dict) and isinstance(ordinary_labels["rows"], list)
    external_labels = {
        "schema_version": "calibration-labels-v1",
        "rows": [_label_row("held"), *ordinary_labels["rows"]],
    }
    labels = decode_label_artifact(external_labels)
    run_values = runs_raw["runs"]
    assert isinstance(run_values, dict)
    extracted = extract_evidence(study, labels, run_values)
    feature_rows = [
        {"feature_key": row.feature_key, "manifest_key": row.manifest_key, "features": dict(row.features)}
        for row in extracted.features.rows
    ]
    label_document = external_labels
    baseline = thaw_json(extracted.baseline)
    assert isinstance(baseline, dict)
    locked_keys = tuple(member.key for member in study.partitions.members if member.role == "locked-heldout")
    inputs = build_role_inputs(locked_keys, feature_rows, label_document["rows"], baseline, extracted.run_hashes)
    payload = {
        "schema_version": "calibration-locked-payload-v1",
        "study": study_raw,
        "features": inputs.features,
        "labels": inputs.labels,
        "baseline": inputs.baseline,
        "run_hashes": inputs.run_hashes,
    }
    payload_raw = canonical_json_bytes(payload)
    (root / "locked_payload.json").write_bytes(payload_raw)
    locked_evidence = decode_locked_payload(payload_raw)
    validation_attestation = load_strict_json_object((validation_root / "attestation.json").read_bytes())
    write_json(root / "validation_attestation.json", validation_attestation)
    study_binding_document = load_strict_json_object((validation_root / "study_binding.json").read_bytes())
    write_json(root / "study_binding.json", study_binding_document)
    binding = decode_study_binding(study_binding_document)
    profile = resolve_decision_profile(profile_path)
    metadata = profile.document["generated_metadata"]
    assert isinstance(metadata, dict) or hasattr(metadata, "get")
    authority = {
        "schema_version": "calibration-custodian-authority-v1",
        "authority_kind": "named-external-custodian-attestation",
        "custodian_name": "Independent Example Repository",
        "attestation_id": "IER-295-0001",
        "role": "locked-heldout",
        "status": "authorized",
        "study_sha256": study.sha256,
        "protocol_sha256": study.protocol.sha256,
        "partition_manifest_sha256": study.partitions.sha256,
        "profile_sha256": profile.digest,
        "profile_dataset_sha256": metadata.get("dataset_manifest_hash"),
        "locked_payload_sha256": hashlib.sha256(payload_raw).hexdigest(),
        "locked_dataset_sha256": locked_evidence.dataset_sha256,
        "locked_run_hashes_sha256": canonical_sha256(inputs.run_hashes),
        "validation_attestation_sha256": hashlib.sha256(
            (root / "validation_attestation.json").read_bytes()
        ).hexdigest(),
        "validation_evidence_sha256": validation_attestation["evidence_sha256"],
        "study_binding_sha256": hashlib.sha256((root / "study_binding.json").read_bytes()).hexdigest(),
        "run_commitments_sha256": binding.run_commitments_sha256,
        "validation_role_run_commitments_sha256": binding.role_run_commitments_sha256["validation"],
        "validation_role_run_artifacts_sha256": binding.role_run_artifacts_sha256["validation"],
        "locked_role_run_commitments_sha256": binding.role_run_commitments_sha256["locked-heldout"],
        "locked_role_run_artifacts_sha256": binding.role_run_artifacts_sha256["locked-heldout"],
    }
    write_json(root / "authority_attestation.json", authority)
    write_checksums(root)
    return root


def _refresh_locked_authority(root: Path) -> None:
    payload_raw = (root / "locked_payload.json").read_bytes()
    payload = load_strict_json_object(payload_raw)
    study = decode_study_declaration(payload["study"])
    features = decode_feature_artifact(payload["features"])
    labels = decode_label_artifact(payload["labels"])
    baseline = payload["baseline"]
    run_hashes_raw = payload["run_hashes"]
    assert isinstance(run_hashes_raw, dict)
    run_hashes = decode_run_hashes(run_hashes_raw)
    dataset_sha256 = canonical_sha256(
        {
            "study_sha256": study.sha256,
            "features_sha256": features.sha256,
            "labels_sha256": labels.sha256,
            "baseline_sha256": canonical_sha256(baseline),
            "run_artifact_sha256": run_hashes,
        }
    )
    authority_path = root / "authority_attestation.json"
    authority = load_strict_json_object(authority_path.read_bytes())
    authority["locked_payload_sha256"] = hashlib.sha256(payload_raw).hexdigest()
    authority["locked_dataset_sha256"] = dataset_sha256
    authority["locked_run_hashes_sha256"] = canonical_sha256(run_hashes)
    write_json(authority_path, authority)
    write_checksums(root)


def _run_cli(argv: list[str]) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    handle_calibrate(args, {}, parser, logging.INFO, None)
