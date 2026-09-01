"""Local precommit, one-use consumption, and candidate retirement guards."""

import hashlib
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.calibration_contract import decode_evidence_manifest
from vntyper.scripts.calibration_custody import (
    open_locked_payload,
    record_consumption,
    retire_candidate,
    write_precommit,
)
from vntyper.scripts.calibration_profiles import build_generated_profile
from vntyper.scripts.calibration_workflow import evaluate_locked_candidate

pytestmark = pytest.mark.unit


def _profile(margin: int = 1):
    return build_generated_profile(
        {
            "enabled": True,
            "minimum_record_count_margin": margin,
            "minimum_record_share": 0.5,
            "minimum_record_share_margin": 0.0,
            "xd_veto": "disabled",
            "abstain_on_inadmissible_advntr": False,
        },
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=295,
        objective="lexicographic-safety-v1",
        generator_version="2.0.26",
    )


def _evidence(payload_digest: str):
    return decode_evidence_manifest(
        {
            "schema_version": "calibration-evidence-v1",
            "role": "locked-heldout",
            "provenance": "external-custodian",
            "protocol_sha256": "c" * 64,
            "partition_manifest_sha256": "d" * 64,
            "features_sha256": payload_digest,
            "labels_sha256": "e" * 64,
            "baseline_sha256": "f" * 64,
        }
    )


def test_precommit_is_durable_and_hash_bound_before_locked_payload_opens(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked evidence\n")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    profile = _profile()
    precommit = write_precommit(tmp_path / "custody", profile.digest, "c" * 64, digest)

    assert precommit.path.is_file()
    assert precommit.profile_sha256 == profile.digest
    assert precommit.protocol_sha256 == "c" * 64
    assert precommit.evidence_sha256 == digest

    opened = open_locked_payload(payload, precommit, tmp_path / "custody")
    assert opened.payload == b"locked evidence\n"
    assert opened.receipt.path.is_file()


def test_second_use_is_refused_even_for_another_profile(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"one use")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    custody = tmp_path / "custody"
    first = write_precommit(custody, _profile(1).digest, "c" * 64, digest)
    open_locked_payload(payload, first, custody)
    second = write_precommit(custody, _profile(2).digest, "c" * 64, digest)

    with pytest.raises(ValueError, match="consumed|one-use"):
        open_locked_payload(payload, second, custody)


def test_tampered_precommit_is_refused_before_consumption(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    custody = tmp_path / "custody"
    precommit = write_precommit(custody, _profile().digest, "c" * 64, digest)
    precommit.path.write_text("{}\n", encoding="utf-8")

    with pytest.raises(ValueError, match="precommit"):
        open_locked_payload(payload, precommit, custody)
    assert not (custody / "consumed" / f"{digest}.json").exists()


def test_payload_hash_failure_still_consumes_the_evidence_hash(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"before")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    custody = tmp_path / "custody"
    precommit = write_precommit(custody, _profile().digest, "c" * 64, digest)
    payload.write_bytes(b"after")

    with pytest.raises(ValueError, match="hash"):
        open_locked_payload(payload, precommit, custody)
    assert (custody / "consumed" / f"{digest}.json").is_file()


def test_failed_open_consumes_evidence_and_retires_candidate(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    profile = _profile()
    evidence = _evidence(digest)

    with (
        patch("vntyper.scripts.calibration_workflow.open_locked_payload", side_effect=RuntimeError("custodian failed")),
        pytest.raises(RuntimeError, match="custodian failed"),
    ):
        evaluate_locked_candidate(
            profile,
            payload,
            evidence,
            tmp_path / "custody",
            evaluator=lambda raw: {"bytes": len(raw)},
        )

    retired = tuple((tmp_path / "custody" / "retired").glob("*.json"))
    assert len(retired) == 1


def test_failed_evaluator_cannot_reuse_locked_evidence(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    profile = _profile()
    evidence = _evidence(digest)

    with pytest.raises(RuntimeError, match="evaluation failed"):
        evaluate_locked_candidate(
            profile,
            payload,
            evidence,
            tmp_path / "custody",
            evaluator=lambda raw: (_ for _ in ()).throw(RuntimeError("evaluation failed")),
        )
    with pytest.raises(ValueError, match="consumed|one-use"):
        evaluate_locked_candidate(
            _profile(2),
            payload,
            evidence,
            tmp_path / "custody",
            evaluator=lambda raw: {"bytes": len(raw)},
        )


def test_successful_locked_evaluation_is_role_and_hash_bound(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    profile = _profile()

    result = evaluate_locked_candidate(
        profile,
        payload,
        _evidence(digest),
        tmp_path / "custody",
        evaluator=lambda raw: {"bytes": len(raw)},
    )

    assert result.result == {"bytes": 6}
    assert result.profile_sha256 == profile.digest
    assert result.evidence_sha256 == digest
    assert result.role == "locked-heldout"
    assert result.receipt.path.is_file()


def test_non_heldout_or_development_evidence_cannot_enter_evaluation(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    validation = decode_evidence_manifest(
        {
            "schema_version": "calibration-evidence-v1",
            "role": "validation",
            "provenance": "development",
            "protocol_sha256": "c" * 64,
            "partition_manifest_sha256": "d" * 64,
            "features_sha256": digest,
            "labels_sha256": "e" * 64,
            "baseline_sha256": "f" * 64,
        }
    )

    with pytest.raises(ValueError, match="held-out|custodian"):
        evaluate_locked_candidate(_profile(), payload, validation, tmp_path / "custody", evaluator=lambda raw: len(raw))


def test_consumption_and_retirement_records_are_append_only(tmp_path: Path) -> None:
    custody = tmp_path / "custody"
    receipt = record_consumption(custody, "a" * 64, "b" * 64, "c" * 64)
    retirement = retire_candidate(custody, "a" * 64, "b" * 64, "failed-validation")

    assert receipt.path.is_file() and retirement.is_file()
    with pytest.raises(ValueError):
        record_consumption(custody, "d" * 64, "b" * 64, "c" * 64)
    with pytest.raises(ValueError):
        retire_candidate(custody, "a" * 64, "b" * 64, "failed-again")


def test_custody_records_reject_malformed_hashes_and_empty_reasons(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="SHA-256"):
        write_precommit(tmp_path, "A" * 64, "b" * 64, "c" * 64)
    with pytest.raises(ValueError, match="reason"):
        retire_candidate(tmp_path, "a" * 64, "b" * 64, " ")
