"""Local precommit, one-use consumption, and candidate retirement guards."""

import hashlib
import threading
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts import calibration_custody
from vntyper.scripts.calibration_contract import decode_evidence_manifest
from vntyper.scripts.calibration_custody import (
    open_locked_payload,
    record_consumption,
    require_candidate_active,
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
            evidence.protocol_sha256,
            evidence.features_sha256,
            tmp_path / "custody",
            evaluator=lambda raw: {"bytes": len(raw)},
        )

    retired = tuple((tmp_path / "custody" / "retired").glob("*.json"))
    assert len(retired) == 1


@pytest.mark.parametrize("interruption", [KeyboardInterrupt(), SystemExit(9)])
def test_interrupted_locked_access_retires_the_pair_and_cannot_retry(
    tmp_path: Path, interruption: BaseException
) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    profile = _profile()
    evidence = _evidence(digest)
    custody = tmp_path / "custody"

    with (
        patch("vntyper.scripts.calibration_workflow.open_locked_payload", side_effect=interruption),
        pytest.raises(type(interruption)),
    ):
        evaluate_locked_candidate(
            profile,
            payload,
            evidence.protocol_sha256,
            evidence.features_sha256,
            custody,
            evaluator=lambda raw: len(raw),
        )

    assert len(tuple((custody / "retired").glob("*.json"))) == 1
    with pytest.raises(ValueError, match="precommit|retired|one-use"):
        evaluate_locked_candidate(
            profile,
            payload,
            evidence.protocol_sha256,
            evidence.features_sha256,
            custody,
            evaluator=lambda raw: len(raw),
        )


@pytest.mark.parametrize("interruption", [KeyboardInterrupt(), SystemExit(9)])
def test_interrupted_exclusive_write_removes_only_its_partial_target(
    tmp_path: Path, interruption: BaseException
) -> None:
    custody = tmp_path / "custody"
    neighbour = custody / "precommits" / "keep.json"
    neighbour.parent.mkdir(parents=True)
    neighbour.write_text("keep\n", encoding="utf-8")

    with (
        patch("vntyper.scripts.calibration_custody.os.fsync", side_effect=interruption),
        pytest.raises(type(interruption)),
    ):
        write_precommit(custody, "a" * 64, "b" * 64, "c" * 64)

    assert neighbour.read_text(encoding="utf-8") == "keep\n"
    assert not (custody / "precommits" / f"{'a' * 64}.{'c' * 64}.json").exists()


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
            evidence.protocol_sha256,
            evidence.features_sha256,
            tmp_path / "custody",
            evaluator=lambda raw: (_ for _ in ()).throw(RuntimeError("evaluation failed")),
        )
    with pytest.raises(ValueError, match="consumed|one-use"):
        evaluate_locked_candidate(
            _profile(2),
            payload,
            evidence.protocol_sha256,
            evidence.features_sha256,
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
        _evidence(digest).protocol_sha256,
        digest,
        tmp_path / "custody",
        evaluator=lambda raw: {"bytes": len(raw)},
    )

    assert result.result == {"bytes": 6}
    assert result.profile_sha256 == profile.digest
    assert result.evidence_sha256 == digest
    assert result.role == "locked-heldout"
    assert result.receipt.path.is_file()


def test_malformed_authority_hashes_cannot_enter_locked_evaluation(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="SHA-256"):
        evaluate_locked_candidate(
            _profile(), payload, "invalid", digest, tmp_path / "custody", evaluator=lambda raw: len(raw)
        )


def test_consumption_and_retirement_records_are_append_only(tmp_path: Path) -> None:
    custody = tmp_path / "custody"
    receipt = record_consumption(custody, "a" * 64, "b" * 64, "c" * 64)
    retirement = retire_candidate(custody, "a" * 64, "b" * 64, "failed-validation")

    assert receipt.path.is_file() and retirement.is_file()
    with pytest.raises(ValueError):
        record_consumption(custody, "d" * 64, "b" * 64, "c" * 64)
    with pytest.raises(ValueError):
        retire_candidate(custody, "a" * 64, "b" * 64, "failed-again")


def test_retired_profile_evidence_pair_is_refused_before_another_operation(tmp_path: Path) -> None:
    custody = tmp_path / "custody"
    require_candidate_active(custody, "a" * 64, "b" * 64)
    retire_candidate(custody, "a" * 64, "b" * 64, "failed-validation")

    with pytest.raises(ValueError, match="retired"):
        require_candidate_active(custody, "a" * 64, "b" * 64)
    require_candidate_active(custody, "c" * 64, "b" * 64)


def test_retirement_inserted_between_active_check_and_precommit_prevents_evaluation(tmp_path: Path) -> None:
    payload = tmp_path / "locked.json"
    payload.write_bytes(b"locked")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    profile = _profile()
    custody = tmp_path / "custody"
    evaluated = False
    retirement_entered = threading.Event()
    release_retirement = threading.Event()
    real_exclusive_write = calibration_custody._exclusive_write
    errors: list[ValueError] = []

    def blocking_write(path, raw, message):
        if path.parent.name == "retired":
            retirement_entered.set()
            assert release_retirement.wait(timeout=5)
        return real_exclusive_write(path, raw, message)

    def evaluator(_raw: bytes) -> int:
        nonlocal evaluated
        evaluated = True
        return 1

    def retire() -> None:
        retire_candidate(custody, profile.digest, digest, "concurrent-retirement")

    def evaluate() -> None:
        try:
            evaluate_locked_candidate(profile, payload, "c" * 64, digest, custody, evaluator=evaluator)
        except ValueError as error:
            errors.append(error)

    with patch("vntyper.scripts.calibration_custody._exclusive_write", side_effect=blocking_write):
        retirement_thread = threading.Thread(target=retire)
        retirement_thread.start()
        assert retirement_entered.wait(timeout=5)
        evaluation_thread = threading.Thread(target=evaluate)
        evaluation_thread.start()
        release_retirement.set()
        retirement_thread.join(timeout=5)
        evaluation_thread.join(timeout=5)

    assert not evaluated
    assert len(errors) == 1 and isinstance(errors[0], ValueError) and "retired" in str(errors[0])
    assert not tuple((custody / "consumed").glob("*.json"))


def test_custody_records_reject_malformed_hashes_and_empty_reasons(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="SHA-256"):
        write_precommit(tmp_path, "A" * 64, "b" * 64, "c" * 64)
    with pytest.raises(ValueError, match="reason"):
        retire_candidate(tmp_path, "a" * 64, "b" * 64, " ")
