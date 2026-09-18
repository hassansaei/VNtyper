"""Target custody records consumption before access and never reopen interrupted claims."""

import hashlib
from importlib import import_module
from pathlib import Path
from unittest.mock import Mock

import pytest

from tests.unit.test_calibration_target_access import access_inputs
from vntyper.scripts.calibration_exposure import decode_exposure_receipt
from vntyper.scripts.calibration_target_access import preflight_target_confirmation
from vntyper.scripts.calibration_target_attestation import (
    decode_target_locked_attestation,
    decode_target_validation_attestation,
    target_validation_attestation_document,
)
from vntyper.scripts.calibration_target_authority import (
    decode_target_custodian_authority,
    target_custodian_authority_document,
)

pytestmark = pytest.mark.unit


def _module():
    return import_module("vntyper.scripts.calibration_target_custody")


def _exposure(confirmation, sequence=1):
    return decode_exposure_receipt(
        {
            "schema_version": "calibration-target-exposure-receipt-v2",
            "target": confirmation.candidate.target,
            "role": confirmation.role,
            "study_sha256": confirmation.study.sha256,
            "partition_sha256": confirmation.study.partitions.sha256,
            "evidence_sha256": confirmation.evidence_sha256,
            "membership_sha256": "0" * 64,
            "exposure_ledger_id": confirmation.study.exposure_ledger_id,
            "sequence": sequence,
        }
    )


def _validation():
    candidate, study, validation, _ = access_inputs()
    confirmation = preflight_target_confirmation(
        candidate=candidate,
        study=study,
        role="validation",
        evidence_sha256=validation.evidence_sha256,
        run_manifest_sha256=validation.run_manifest_sha256,
    )
    exposure = _exposure(confirmation)
    validation = decode_target_validation_attestation(
        {
            **target_validation_attestation_document(validation),
            "exposure_receipt_sha256": exposure.sha256,
        }
    )
    return confirmation, exposure, validation


def _root(tmp_path, confirmation):
    root = tmp_path / "custody"
    _module().initialize_target_custody(root, exposure_ledger_id=confirmation.study.exposure_ledger_id)
    return root


def _locked(root, *, payload=b'{"synthetic":"truth"}\n'):
    confirmation, exposure, validation = _validation()
    with _module().claim_target_confirmation(root, confirmation, exposure) as claim:
        claim.finish(validation)
    _, _, _, previous = access_inputs()
    authority = decode_target_custodian_authority(
        {
            **target_custodian_authority_document(previous),
            "validation_attestation_sha256": validation.sha256,
            "validation_exposure_receipt_sha256": exposure.sha256,
            "locked_payload_sha256": hashlib.sha256(payload).hexdigest(),
        }
    )
    locked = preflight_target_confirmation(
        candidate=confirmation.candidate,
        study=confirmation.study,
        role="locked-heldout",
        evidence_sha256=authority.locked_heldout_evidence_sha256,
        run_manifest_sha256=authority.locked_heldout_run_manifest_sha256,
        validation=validation,
        authority=authority,
    )
    return locked, _exposure(locked, 2)


def _locked_result(confirmation, exposure, consumption_sha256):
    return decode_target_locked_attestation(
        {
            **target_validation_attestation_document(confirmation.validation),
            "schema_version": "calibration-target-locked-attestation-v2",
            "role": "locked-heldout",
            "evidence_sha256": confirmation.evidence_sha256,
            "run_manifest_sha256": confirmation.run_manifest_sha256,
            "exposure_receipt_sha256": exposure.sha256,
            "validation_attestation_sha256": confirmation.validation.sha256,
            "custodian_authority_sha256": confirmation.authority.sha256,
            "custody_consumption_receipt_sha256": consumption_sha256,
        }
    )


def test_validation_then_locked_consumes_before_read_and_finishes_once(tmp_path: Path):
    initial, _, _ = _validation()
    root = _root(tmp_path, initial)
    locked, exposure = _locked(root)
    payload = b'{"synthetic":"truth"}\n'
    with _module().claim_target_confirmation(root, locked, exposure) as claim:

        def reader():
            assert tuple(root.glob("*.locked-consumption.json"))
            return payload

        opened = claim.open_locked_payload(reader)
        assert opened.payload == payload and len(opened.receipt_sha256) == 64
        assert opened.receipt["exposure_receipt_sha256"] == exposure.sha256
        claim.finish(_locked_result(locked, exposure, opened.receipt_sha256))
        with pytest.raises(ValueError):
            claim.open_locked_payload(reader)
    with pytest.raises(ValueError, match="already|terminal|consumed"):
        _module().claim_target_confirmation(root, locked, exposure)
    assert not list(root.glob("*.retired.json"))


def test_interrupted_validation_retires_without_allowing_retry(tmp_path: Path):
    confirmation, exposure, _ = _validation()
    root = _root(tmp_path, confirmation)
    with pytest.raises(RuntimeError), _module().claim_target_confirmation(root, confirmation, exposure):
        raise RuntimeError("synthetic evaluator failure")
    assert len(list(root.glob("*.retired.json"))) == 1
    with pytest.raises(ValueError):
        _module().claim_target_confirmation(root, confirmation, exposure)


def test_crashed_claim_without_terminal_still_cannot_restart(tmp_path: Path):
    confirmation, exposure, _ = _validation()
    root = _root(tmp_path, confirmation)
    claim = _module().claim_target_confirmation(root, confirmation, exposure)
    claim.close()
    assert not list(root.glob("*.retired.json"))
    with pytest.raises(ValueError):
        _module().claim_target_confirmation(root, confirmation, exposure)


@pytest.mark.parametrize("failure", ["reader", "wrong-bytes", "nonbytes"])
def test_locked_access_failure_preserves_consumption(tmp_path: Path, failure):
    confirmation, _, _ = _validation()
    root = _root(tmp_path, confirmation)
    locked, exposure = _locked(root)
    reader = (
        Mock(side_effect=OSError("synthetic read failure"))
        if failure == "reader"
        else Mock(return_value=b"changed" if failure == "wrong-bytes" else "text")
    )
    with pytest.raises((OSError, ValueError)), _module().claim_target_confirmation(root, locked, exposure) as claim:
        claim.open_locked_payload(reader)
    reader.assert_called_once()
    assert len(list(root.glob("*.locked-consumption.json"))) == 1
    assert len(list(root.glob("*.retired.json"))) == 1
    with pytest.raises(ValueError):
        _module().claim_target_confirmation(root, locked, exposure)


def test_locked_requires_local_passed_validation_and_consumption_before_finish(tmp_path: Path):
    confirmation, _, _ = _validation()
    root = _root(tmp_path, confirmation)
    locked, exposure = _locked(root)
    with (
        _module().claim_target_confirmation(root, locked, exposure) as claim,
        pytest.raises(ValueError, match="consum"),
    ):
        claim.finish(_locked_result(locked, exposure, "0" * 64))
    fresh = tmp_path / "fresh"
    _module().initialize_target_custody(fresh, exposure_ledger_id=confirmation.study.exposure_ledger_id)
    with pytest.raises(ValueError, match="validation"):
        _module().claim_target_confirmation(fresh, locked, exposure)


def test_failed_validation_never_creates_passed_terminal(tmp_path: Path):
    confirmation, exposure, validation = _validation()
    root = _root(tmp_path, confirmation)
    failed = decode_target_validation_attestation(
        {**target_validation_attestation_document(validation), "status": "failed"}
    )
    with _module().claim_target_confirmation(root, confirmation, exposure) as claim:
        claim.finish(failed)
    assert len(list(root.glob("*.validation-result.json"))) == 1
    with pytest.raises(ValueError):
        _module().claim_target_confirmation(root, confirmation, exposure)


def test_external_private_initialization_and_wrong_ledger_refusal(tmp_path: Path):
    confirmation, exposure, _ = _validation()
    root = _root(tmp_path, confirmation)
    assert root.stat().st_mode & 0o777 == 0o700
    with pytest.raises(ValueError):
        _module().initialize_target_custody(root, exposure_ledger_id=confirmation.study.exposure_ledger_id)
    other = tmp_path / "other"
    _module().initialize_target_custody(other, exposure_ledger_id="0" * 64)
    with pytest.raises(ValueError, match="ledger"):
        _module().claim_target_confirmation(other, confirmation, exposure)
    repository = tmp_path / "repository"
    repository.mkdir()
    (repository / ".git").mkdir()
    with pytest.raises(ValueError, match="outside"):
        _module().initialize_target_custody(repository / "custody", exposure_ledger_id="0" * 64)
    with pytest.raises(ValueError):
        _module().initialize_target_custody(
            tmp_path / "inside", exposure_ledger_id="0" * 64, forbidden_roots=(tmp_path,)
        )


def test_caught_payload_failure_cannot_be_completed_as_a_success(tmp_path: Path):
    confirmation, _, _ = _validation()
    root = _root(tmp_path, confirmation)
    locked, exposure = _locked(root)
    with _module().claim_target_confirmation(root, locked, exposure) as claim:
        with pytest.raises(ValueError, match="payload bytes"):
            claim.open_locked_payload(lambda: b"wrong")
        with pytest.raises(ValueError, match="verified|consum"):
            claim.finish(_locked_result(locked, exposure, claim.consumption_sha256))
    assert len(list(root.glob("*.retired.json"))) == 1


def test_consumption_write_failure_never_opens_payload(tmp_path: Path):
    from unittest.mock import patch

    confirmation, _, _ = _validation()
    root = _root(tmp_path, confirmation)
    locked, exposure = _locked(root)
    reader = Mock(return_value=b'{"synthetic":"truth"}\n')
    with _module().claim_target_confirmation(root, locked, exposure) as claim:
        with patch.object(_module().os, "write", return_value=0), pytest.raises(ValueError):
            claim.open_locked_payload(reader)
        reader.assert_not_called()
    assert len(list(root.glob("*.locked-consumption.json"))) == 1
    with pytest.raises(ValueError):
        _module().claim_target_confirmation(root, locked, exposure)


def test_wrong_exposure_receipt_is_refused_before_any_claim_write(tmp_path: Path):
    from vntyper.scripts.calibration_exposure import exposure_receipt_document

    confirmation, exposure, _ = _validation()
    root = _root(tmp_path, confirmation)
    wrong = decode_exposure_receipt({**exposure_receipt_document(exposure), "role": "training"})
    before = set(root.iterdir())
    with pytest.raises(ValueError, match="exposure"):
        _module().claim_target_confirmation(root, confirmation, wrong)
    assert set(root.iterdir()) == before


def test_noncanonical_and_unknown_records_deny_access(tmp_path: Path):
    confirmation, exposure, _ = _validation()
    root = _root(tmp_path, confirmation)
    identity = root / "identity.json"
    raw = identity.read_bytes()
    identity.write_bytes(b" " + raw)
    with pytest.raises(ValueError, match="canonical"):
        _module().claim_target_confirmation(root, confirmation, exposure)
    identity.write_bytes(raw)
    (root / "unexpected.json").write_text("{}")
    with pytest.raises(ValueError, match="inventory"):
        _module().claim_target_confirmation(root, confirmation, exposure)


def test_root_and_lock_replacement_during_claim_refuse_writes(tmp_path: Path):
    confirmation, exposure, validation = _validation()
    root = _root(tmp_path, confirmation)
    claim = _module().claim_target_confirmation(root, confirmation, exposure)
    lock = root / f"{confirmation.candidate.sha256}.lock"
    lock.rename(root / "moved-lock")
    lock.write_bytes(b"")
    try:
        with pytest.raises(ValueError, match="lock changed"):
            claim.finish(validation)
    finally:
        claim.close()
    assert not list(root.glob("*.validation-result.json"))


def test_concurrent_validation_claims_have_only_one_winner(tmp_path: Path):
    from concurrent.futures import ThreadPoolExecutor

    confirmation, exposure, validation = _validation()
    root = _root(tmp_path, confirmation)

    def attempt(_):
        try:
            with _module().claim_target_confirmation(root, confirmation, exposure) as claim:
                claim.finish(validation)
            return True
        except ValueError:
            return False

    with ThreadPoolExecutor(max_workers=2) as pool:
        assert sorted(pool.map(attempt, (1, 2))) == [False, True]
    assert len(list(root.glob("*.validation-result.json"))) == 1


def test_preflight_target_custody_validation_and_lifecycle(tmp_path: Path):
    confirmation, exposure, validation = _validation()
    root = _root(tmp_path, confirmation)
    module = _module()

    # Clean custody passes validation preflight
    module.preflight_target_custody(root, confirmation)

    # Corrupt lock file fails preflight
    prefix = confirmation.candidate.sha256
    lock_file = root / f"{prefix}.lock"
    lock_file.write_bytes(b"non-empty")
    with pytest.raises(ValueError, match="lock is not an empty regular file"):
        module.preflight_target_custody(root, confirmation)
    lock_file.unlink()

    # Completed validation candidate passes locked preflight and fails validation preflight
    locked_confirmation, _ = _locked(root)
    module.preflight_target_custody(root, locked_confirmation)

    with pytest.raises(ValueError, match="validation was already claimed or completed"):
        module.preflight_target_custody(root, confirmation)
