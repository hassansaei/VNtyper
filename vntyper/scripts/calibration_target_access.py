"""Pure target-aware confirmation preflight, before outcome or custody I/O."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_candidate import CandidateEnvelope, candidate_document
from vntyper.scripts.calibration_exposure import require_digest
from vntyper.scripts.calibration_target_attestation import (
    TargetValidationAttestation,
    target_validation_attestation_document,
)
from vntyper.scripts.calibration_target_authority import TargetCustodianAuthority, target_custodian_authority_document
from vntyper.scripts.calibration_target_contract import TargetStudy, target_study_document
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)
ConfirmationRole = Literal["validation", "locked-heldout"]


@dataclass(frozen=True)
class TargetConfirmation:
    """Verified lineage; external custody and exposure recording are still required."""

    candidate: CandidateEnvelope
    study: TargetStudy
    role: ConfirmationRole
    evidence_sha256: str
    run_manifest_sha256: str
    validation: TargetValidationAttestation | None
    authority: TargetCustodianAuthority | None
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _common_bindings(candidate: CandidateEnvelope, study: TargetStudy) -> dict[str, str]:
    return {
        "target": candidate.target,
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "study_sha256": study.sha256,
        "protocol_sha256": study.protocol.sha256,
        "partition_sha256": study.partitions.sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "exposure_ledger_id": study.exposure_ledger_id,
    }


def _validate_locked(
    candidate: CandidateEnvelope,
    study: TargetStudy,
    evidence_sha256: str,
    run_manifest_sha256: str,
    validation: TargetValidationAttestation | None,
    authority: TargetCustodianAuthority | None,
) -> None:
    if not isinstance(validation, TargetValidationAttestation) or not isinstance(authority, TargetCustodianAuthority):
        _fail("locked confirmation requires validation and external authority")
    target_validation_attestation_document(validation)
    target_custodian_authority_document(authority)
    if validation.status != "passed":
        _fail("locked confirmation requires passed validation")
    for name, expected in _common_bindings(candidate, study).items():
        if getattr(validation, name) != expected or getattr(authority, name) != expected:
            _fail(f"locked confirmation {name} differs from the opened candidate/study")
    pairs = {
        "validation_evidence_sha256": validation.evidence_sha256,
        "validation_run_manifest_sha256": validation.run_manifest_sha256,
        "validation_attestation_sha256": validation.sha256,
        "validation_exposure_receipt_sha256": validation.exposure_receipt_sha256,
        "locked_heldout_evidence_sha256": evidence_sha256,
        "locked_heldout_run_manifest_sha256": run_manifest_sha256,
    }
    if any(getattr(authority, name) != expected for name, expected in pairs.items()):
        _fail("locked confirmation authority differs from prior validation or locked commitments")


def _document(
    candidate: CandidateEnvelope,
    study: TargetStudy,
    role: ConfirmationRole,
    evidence_sha256: str,
    run_manifest_sha256: str,
    validation: TargetValidationAttestation | None,
    authority: TargetCustodianAuthority | None,
) -> dict[str, object]:
    return {
        "schema_version": "calibration-target-confirmation-preflight-v2",
        **_common_bindings(candidate, study),
        "role": role,
        "evidence_sha256": evidence_sha256,
        "run_manifest_sha256": run_manifest_sha256,
        "validation_attestation_sha256": None if validation is None else validation.sha256,
        "custodian_authority_sha256": None if authority is None else authority.sha256,
        "locked_payload_sha256": None if authority is None else authority.locked_payload_sha256,
    }


def preflight_target_confirmation(
    *,
    candidate: CandidateEnvelope,
    study: TargetStudy,
    role: str,
    evidence_sha256: str,
    run_manifest_sha256: str,
    validation: TargetValidationAttestation | None = None,
    authority: TargetCustodianAuthority | None = None,
) -> TargetConfirmation:
    """Check target, study and authority lineage without reading outcomes.

    Args:
        candidate: Fixed research candidate, whose payload is verified separately.
        study: Opened immutable study declaration.
        role: Exactly validation or locked-heldout.
        evidence_sha256: Authorized evidence commitment from role metadata.
        run_manifest_sha256: Exact role run-manifest commitment.
        validation: Prior passed validation, required only for locked evaluation.
        authority: External locked authorization, required only for locked evaluation.

    Returns:
        Immutable verified lineage. This is not a consumption or exposure receipt.

    Raises:
        ValueError: If a binding, role, typed identity or prior authorization differs.
    """
    candidate_document(candidate)
    target_study_document(study)
    if (
        candidate.target != study.target
        or candidate.study_sha256 != study.sha256
        or candidate.partition_sha256 != study.partitions.sha256
        or candidate.applicability != study.applicability
        or candidate.producer != study.baseline.producer
    ):
        _fail("confirmation candidate differs from the opened study")
    if not isinstance(role, str) or role not in {"validation", "locked-heldout"}:
        _fail("confirmation role must be validation or locked-heldout")
    role_value = cast(ConfirmationRole, role)
    evidence = require_digest(evidence_sha256, "confirmation evidence")
    runs = require_digest(run_manifest_sha256, "confirmation run manifest")
    if role_value == "validation":
        if validation is not None or authority is not None:
            _fail("validation confirmation cannot accept previous or locked authority")
    else:
        _validate_locked(candidate, study, evidence, runs, validation, authority)
    document = _document(candidate, study, role_value, evidence, runs, validation, authority)
    return TargetConfirmation(
        candidate, study, role_value, evidence, runs, validation, authority, canonical_sha256(document)
    )


def confirmation_document(confirmation: TargetConfirmation) -> dict[str, object]:
    """Project preflight lineage after revalidating all immutable bindings.

    Args:
        confirmation: Previously checked confirmation lineage.

    Returns:
        Aggregate-free content, with no named external authority or paths.

    Raises:
        ValueError: If the typed lineage or its derived identity changed.
    """
    if not isinstance(confirmation, TargetConfirmation):
        _fail("confirmation projection requires typed lineage")
    expected = preflight_target_confirmation(
        candidate=confirmation.candidate,
        study=confirmation.study,
        role=confirmation.role,
        evidence_sha256=confirmation.evidence_sha256,
        run_manifest_sha256=confirmation.run_manifest_sha256,
        validation=confirmation.validation,
        authority=confirmation.authority,
    )
    if expected != confirmation:
        _fail("confirmation identity differs from its content")
    return _document(
        confirmation.candidate,
        confirmation.study,
        confirmation.role,
        confirmation.evidence_sha256,
        confirmation.run_manifest_sha256,
        confirmation.validation,
        confirmation.authority,
    )
