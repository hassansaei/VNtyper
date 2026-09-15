"""Aggregate-free portable projection of passed target calibration evidence."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NoReturn, cast

from vntyper.scripts.calibration_candidate import (
    CalibrationTarget,
    CandidateEnvelope,
    candidate_applicability_document,
    candidate_document,
    candidate_producer_document,
    decode_candidate,
)
from vntyper.scripts.calibration_target_attestation import (
    TargetLockedAttestation,
    TargetValidationAttestation,
    decode_target_locked_attestation,
    decode_target_validation_attestation,
    target_locked_attestation_document,
    target_validation_attestation_document,
)
from vntyper.scripts.calibration_target_authority import (
    TargetCustodianAuthority,
    decode_target_custodian_authority,
    target_custodian_authority_document,
)
from vntyper.scripts.calibration_target_completion import (
    TargetCompletion,
    decode_target_completion,
    target_completion_document,
)
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_DIGEST_FIELDS = (
    "candidate_sha256",
    "candidate_id",
    "payload_sha256",
    "study_sha256",
    "protocol_sha256",
    "partition_sha256",
    "baseline_sha256",
    "exposure_ledger_id",
    "applicability_sha256",
    "producer_sha256",
    "validation_attestation_sha256",
    "locked_heldout_attestation_sha256",
    "custodian_authority_sha256",
    "custody_completion_sha256",
)
_FIELDS = {"schema_version", "target", "disposition", *_DIGEST_FIELDS}


@dataclass(frozen=True)
class PortableApproval:
    """Runtime-safe approval identities without metrics, labels, or custodian names."""

    target: CalibrationTarget
    candidate_sha256: str
    candidate_id: str
    payload_sha256: str
    study_sha256: str
    protocol_sha256: str
    partition_sha256: str
    baseline_sha256: str
    exposure_ledger_id: str
    applicability_sha256: str
    producer_sha256: str
    validation_attestation_sha256: str
    locked_heldout_attestation_sha256: str
    custodian_authority_sha256: str
    custody_completion_sha256: str
    sha256: str

    @property
    def disposition(self) -> str:
        """Return the fixed evidence disposition represented by this projection."""
        return "passed-validation-and-locked-heldout"


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"portable calibration {label} must be lowercase SHA-256")
    return value


def decode_portable_approval(value: object) -> PortableApproval:
    """Decode one closed aggregate-free portable approval projection.

    The projection is checked under the repository's trusted-operator model. Its
    digest is an integrity binding, not cryptographic proof about its author.

    Args:
        value: Parsed JSON projection.

    Returns:
        Immutable approval identities and canonical digest.

    Raises:
        ValueError: If fields, target, disposition, or identities differ.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        _fail("portable calibration approval fields differ from the closed contract")
    if value["schema_version"] != "calibration-portable-approval-v2":
        _fail("portable calibration approval schema is unsupported")
    if value["disposition"] != "passed-validation-and-locked-heldout":
        _fail("portable calibration approval disposition differs")
    target = value["target"]
    if target not in ("callers", "length"):
        _fail("portable calibration approval target must be callers or length")
    return PortableApproval(
        target=cast(CalibrationTarget, target),
        candidate_sha256=_digest(value["candidate_sha256"], "candidate_sha256"),
        candidate_id=_digest(value["candidate_id"], "candidate_id"),
        payload_sha256=_digest(value["payload_sha256"], "payload_sha256"),
        study_sha256=_digest(value["study_sha256"], "study_sha256"),
        protocol_sha256=_digest(value["protocol_sha256"], "protocol_sha256"),
        partition_sha256=_digest(value["partition_sha256"], "partition_sha256"),
        baseline_sha256=_digest(value["baseline_sha256"], "baseline_sha256"),
        exposure_ledger_id=_digest(value["exposure_ledger_id"], "exposure_ledger_id"),
        applicability_sha256=_digest(value["applicability_sha256"], "applicability_sha256"),
        producer_sha256=_digest(value["producer_sha256"], "producer_sha256"),
        validation_attestation_sha256=_digest(value["validation_attestation_sha256"], "validation_attestation_sha256"),
        locked_heldout_attestation_sha256=_digest(
            value["locked_heldout_attestation_sha256"], "locked_heldout_attestation_sha256"
        ),
        custodian_authority_sha256=_digest(value["custodian_authority_sha256"], "custodian_authority_sha256"),
        custody_completion_sha256=_digest(value["custody_completion_sha256"], "custody_completion_sha256"),
        sha256=canonical_sha256(value),
    )


def _document(value: PortableApproval) -> dict[str, object]:
    return {
        "schema_version": "calibration-portable-approval-v2",
        "target": value.target,
        "disposition": value.disposition,
        **{field: getattr(value, field) for field in _DIGEST_FIELDS},
    }


def portable_approval_document(value: PortableApproval) -> dict[str, object]:
    """Project and revalidate one immutable portable approval."""
    if not isinstance(value, PortableApproval):
        _fail("portable calibration approval must be its decoded typed value")
    document = _document(value)
    if decode_portable_approval(document) != value:
        _fail("portable calibration approval differs from its canonical content")
    return document


def build_portable_approval(
    candidate: CandidateEnvelope,
    validation: TargetValidationAttestation,
    locked: TargetLockedAttestation,
    authority: TargetCustodianAuthority,
    completion: TargetCompletion,
) -> PortableApproval:
    """Verify passed target evidence and derive its aggregate-free projection.

    Args:
        candidate: Exact research-only candidate being authorized for runtime.
        validation: Passed validation result.
        locked: Passed locked-heldout result.
        authority: External custodian authorization checked during export.
        completion: Successful receipt-bound target custody terminal.

    Returns:
        Portable approval containing only runtime identities.

    Raises:
        ValueError: If any typed object or direct/transitive binding differs.
    """
    candidate = decode_candidate(candidate_document(candidate))
    validation = _validation(validation)
    authority = _authority(authority)
    locked = _locked(locked)
    completion = _completion(completion)
    if validation.status != "passed" or locked.status != "passed":
        _fail("portable calibration approval requires passed validation and locked-heldout attestations")
    _validate_common(candidate, validation, locked, authority, completion)
    _validate_lineage(validation, locked, authority, completion)
    document = {
        "schema_version": "calibration-portable-approval-v2",
        "target": candidate.target,
        "disposition": "passed-validation-and-locked-heldout",
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "payload_sha256": candidate.payload_sha256,
        "study_sha256": candidate.study_sha256,
        "protocol_sha256": validation.protocol_sha256,
        "partition_sha256": candidate.partition_sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "exposure_ledger_id": validation.exposure_ledger_id,
        "applicability_sha256": canonical_sha256(
            candidate_applicability_document(candidate.applicability, target=candidate.target)
        ),
        "producer_sha256": canonical_sha256(candidate_producer_document(candidate.producer)),
        "validation_attestation_sha256": validation.sha256,
        "locked_heldout_attestation_sha256": locked.sha256,
        "custodian_authority_sha256": authority.sha256,
        "custody_completion_sha256": completion.sha256,
    }
    return decode_portable_approval(document)


def validate_portable_approval_candidate(approval: PortableApproval, candidate: CandidateEnvelope) -> None:
    """Require a runtime candidate and metadata to match a portable approval.

    Args:
        approval: Decoded aggregate-free approval.
        candidate: Decoded candidate loaded from the same portable bundle.

    Raises:
        ValueError: If either typed object or any candidate binding differs.
    """
    approval = decode_portable_approval(portable_approval_document(approval))
    candidate = decode_candidate(candidate_document(candidate))
    expected = (
        candidate.target,
        candidate.sha256,
        candidate.candidate_id,
        candidate.payload_sha256,
        candidate.study_sha256,
        candidate.partition_sha256,
        candidate.baseline_sha256,
        canonical_sha256(candidate_applicability_document(candidate.applicability, target=candidate.target)),
        canonical_sha256(candidate_producer_document(candidate.producer)),
    )
    observed = (
        approval.target,
        approval.candidate_sha256,
        approval.candidate_id,
        approval.payload_sha256,
        approval.study_sha256,
        approval.partition_sha256,
        approval.baseline_sha256,
        approval.applicability_sha256,
        approval.producer_sha256,
    )
    if observed != expected:
        _fail("portable calibration approval and runtime candidate bindings differ")


def _validation(value: TargetValidationAttestation) -> TargetValidationAttestation:
    if not isinstance(value, TargetValidationAttestation):
        _fail("portable calibration validation must be a typed attestation")
    document = target_validation_attestation_document(value)
    return decode_target_validation_attestation(document)


def _locked(value: TargetLockedAttestation) -> TargetLockedAttestation:
    if not isinstance(value, TargetLockedAttestation):
        _fail("portable calibration locked result must be a typed attestation")
    document = target_locked_attestation_document(value)
    return decode_target_locked_attestation(document)


def _authority(value: TargetCustodianAuthority) -> TargetCustodianAuthority:
    if not isinstance(value, TargetCustodianAuthority):
        _fail("portable calibration authority must be its typed contract")
    document = target_custodian_authority_document(value)
    return decode_target_custodian_authority(document)


def _completion(value: TargetCompletion) -> TargetCompletion:
    if not isinstance(value, TargetCompletion):
        _fail("portable calibration completion must be its typed contract")
    document = target_completion_document(value)
    return decode_target_completion(document)


def _validate_common(
    candidate: CandidateEnvelope,
    validation: TargetValidationAttestation,
    locked: TargetLockedAttestation,
    authority: TargetCustodianAuthority,
    completion: TargetCompletion,
) -> None:
    expected = (
        candidate.target,
        candidate.sha256,
        candidate.candidate_id,
        candidate.study_sha256,
        candidate.partition_sha256,
        candidate.baseline_sha256,
    )
    for label, value in (
        ("validation", validation),
        ("locked-heldout", locked),
        ("authority", authority),
        ("completion", completion),
    ):
        observed = (
            value.target,
            value.candidate_sha256,
            value.candidate_id,
            value.study_sha256,
            value.partition_sha256,
            value.baseline_sha256,
        )
        if observed != expected:
            _fail(f"portable calibration {label} candidate or study bindings differ")
    if (
        len({validation.protocol_sha256, locked.protocol_sha256, authority.protocol_sha256, completion.protocol_sha256})
        != 1
    ):
        _fail("portable calibration protocol bindings differ")


def _validate_lineage(
    validation: TargetValidationAttestation,
    locked: TargetLockedAttestation,
    authority: TargetCustodianAuthority,
    completion: TargetCompletion,
) -> None:
    if (
        authority.validation_evidence_sha256 != validation.evidence_sha256
        or authority.validation_run_manifest_sha256 != validation.run_manifest_sha256
        or authority.validation_attestation_sha256 != validation.sha256
        or authority.exposure_ledger_id != validation.exposure_ledger_id
        or authority.validation_exposure_receipt_sha256 != validation.exposure_receipt_sha256
        or locked.evidence_sha256 != authority.locked_heldout_evidence_sha256
        or locked.run_manifest_sha256 != authority.locked_heldout_run_manifest_sha256
        or locked.validation_attestation_sha256 != validation.sha256
        or locked.custodian_authority_sha256 != authority.sha256
        or locked.exposure_ledger_id != validation.exposure_ledger_id
        or completion.validation_attestation_sha256 != validation.sha256
        or completion.locked_heldout_attestation_sha256 != locked.sha256
        or completion.custodian_authority_sha256 != authority.sha256
        or completion.locked_heldout_evidence_sha256 != locked.evidence_sha256
        or completion.locked_heldout_run_manifest_sha256 != locked.run_manifest_sha256
        or completion.exposure_ledger_id != locked.exposure_ledger_id
        or completion.validation_exposure_receipt_sha256 != validation.exposure_receipt_sha256
        or completion.locked_heldout_exposure_receipt_sha256 != locked.exposure_receipt_sha256
        or completion.custody_consumption_receipt_sha256 != locked.custody_consumption_receipt_sha256
    ):
        _fail("portable calibration attestation, authority, exposure, or completion lineage differs")
