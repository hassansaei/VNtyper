"""Confirmation preflight checks lineage before any evidence or custody access."""

from dataclasses import replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_candidate import candidate_document, resign
from tests.unit.test_calibration_portable_projection import _artifacts
from tests.unit.test_calibration_target_contract import study_document
from vntyper.scripts.calibration_candidate import decode_candidate
from vntyper.scripts.calibration_target_attestation import (
    decode_target_validation_attestation,
    target_validation_attestation_document,
)
from vntyper.scripts.calibration_target_authority import (
    decode_target_custodian_authority,
    target_custodian_authority_document,
)
from vntyper.scripts.calibration_target_contract import decode_target_study

pytestmark = pytest.mark.unit


def access_inputs(target="length"):
    study = decode_target_study(study_document(target))
    raw = candidate_document(target)
    raw.update(study_sha256=study.sha256, partition_sha256=study.partitions.sha256)
    candidate = decode_candidate(resign(raw))
    _, old_validation, _, old_authority, _ = _artifacts(target=target)
    common = {
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "study_sha256": study.sha256,
        "protocol_sha256": study.protocol.sha256,
        "partition_sha256": study.partitions.sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "exposure_ledger_id": study.exposure_ledger_id,
    }
    validation = decode_target_validation_attestation(
        {**target_validation_attestation_document(old_validation), **common}
    )
    authority = decode_target_custodian_authority(
        {
            **target_custodian_authority_document(old_authority),
            **common,
            "validation_attestation_sha256": validation.sha256,
        }
    )
    return candidate, study, validation, authority


def _preflight(target="length", **changes):
    candidate, study, validation, authority = access_inputs(target)
    arguments = {
        "candidate": candidate,
        "study": study,
        "role": "locked-heldout",
        "evidence_sha256": authority.locked_heldout_evidence_sha256,
        "run_manifest_sha256": authority.locked_heldout_run_manifest_sha256,
        "validation": validation,
        "authority": authority,
    }
    arguments.update(changes)
    return import_module("vntyper.scripts.calibration_target_access").preflight_target_confirmation(**arguments)


@pytest.mark.parametrize("target", ["length", "callers"])
def test_locked_preflight_binds_one_candidate_and_external_authority(target):
    module = import_module("vntyper.scripts.calibration_target_access")
    result = _preflight(target)
    document = module.confirmation_document(result)
    assert document["target"] == target and document["role"] == "locked-heldout"
    assert document["locked_payload_sha256"] == result.authority.locked_payload_sha256
    assert document["custodian_authority_sha256"] == result.authority.sha256
    assert "custodian_name" not in document
    with pytest.raises(ValueError):
        module.confirmation_document(replace(result, sha256="0" * 64))


def test_validation_preflight_has_no_locked_or_previous_validation_authority():
    module = import_module("vntyper.scripts.calibration_target_access")
    result = _preflight(role="validation", validation=None, authority=None)
    document = module.confirmation_document(result)
    assert document["locked_payload_sha256"] is None
    assert document["validation_attestation_sha256"] is None
    assert result.role == "validation"


@pytest.mark.parametrize(
    "changes",
    [
        {"role": "training"},
        {"role": []},
        {"validation": None},
        {"authority": None},
        {"evidence_sha256": "0" * 64},
        {"run_manifest_sha256": "0" * 64},
        {"role": "validation"},
    ],
)
def test_preflight_refuses_missing_or_wrong_authorization(changes):
    with pytest.raises(ValueError):
        _preflight(**changes)


@pytest.mark.parametrize(
    "field",
    [
        "candidate_sha256",
        "candidate_id",
        "study_sha256",
        "protocol_sha256",
        "partition_sha256",
        "baseline_sha256",
        "exposure_ledger_id",
        "validation_evidence_sha256",
        "validation_run_manifest_sha256",
        "validation_attestation_sha256",
        "validation_exposure_receipt_sha256",
    ],
)
def test_authority_cannot_rebind_any_prior_validation_field(field):
    _, _, _, authority = access_inputs()
    changed = decode_target_custodian_authority({**target_custodian_authority_document(authority), field: "0" * 64})
    with pytest.raises(ValueError):
        _preflight(authority=changed)


def test_failed_validation_never_authorizes_locked_access():
    _, _, validation, _ = access_inputs()
    changed = decode_target_validation_attestation(
        {**target_validation_attestation_document(validation), "status": "failed"}
    )
    with pytest.raises(ValueError, match="passed"):
        _preflight(validation=changed)


def test_candidate_must_match_the_opened_study_even_when_its_hash_is_valid():
    candidate, _, _, _ = access_inputs()
    from vntyper.scripts.calibration_candidate import candidate_document as project

    for field in ("study_sha256", "partition_sha256"):
        changed = decode_candidate(resign({**project(candidate), field: "0" * 64}))
        with pytest.raises(ValueError, match="study"):
            _preflight(candidate=changed)


def test_confirmation_projection_refuses_untyped_content():
    with pytest.raises(ValueError, match="typed"):
        import_module("vntyper.scripts.calibration_target_access").confirmation_document({})
