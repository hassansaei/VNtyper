"""Closed target-aware v2 calibration result attestations."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_candidate import CalibrationTarget
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_COMMON_FIELDS = {
    "schema_version",
    "target",
    "role",
    "status",
    "candidate_sha256",
    "candidate_id",
    "study_sha256",
    "protocol_sha256",
    "partition_sha256",
    "baseline_sha256",
    "evidence_sha256",
    "run_manifest_sha256",
    "metrics_sha256",
    "exposure_ledger_id",
    "exposure_receipt_sha256",
}
_LOCKED_FIELDS = _COMMON_FIELDS | {
    "validation_attestation_sha256",
    "custodian_authority_sha256",
    "custody_consumption_receipt_sha256",
}
_HASH_FIELDS = (
    "candidate_sha256",
    "candidate_id",
    "study_sha256",
    "protocol_sha256",
    "partition_sha256",
    "baseline_sha256",
    "evidence_sha256",
    "run_manifest_sha256",
    "metrics_sha256",
    "exposure_ledger_id",
    "exposure_receipt_sha256",
)


@dataclass(frozen=True)
class TargetValidationAttestation:
    """One target-aware validation outcome with ledger exposure identity."""

    target: CalibrationTarget
    status: Literal["passed", "failed"]
    candidate_sha256: str
    candidate_id: str
    study_sha256: str
    protocol_sha256: str
    partition_sha256: str
    baseline_sha256: str
    evidence_sha256: str
    run_manifest_sha256: str
    metrics_sha256: str
    exposure_ledger_id: str
    exposure_receipt_sha256: str
    sha256: str

    @property
    def role(self) -> str:
        """Return the fixed validation role."""
        return "validation"


@dataclass(frozen=True)
class TargetLockedAttestation:
    """One target-aware locked-heldout outcome and its prior authorizations."""

    target: CalibrationTarget
    status: Literal["passed", "failed"]
    candidate_sha256: str
    candidate_id: str
    study_sha256: str
    protocol_sha256: str
    partition_sha256: str
    baseline_sha256: str
    evidence_sha256: str
    run_manifest_sha256: str
    metrics_sha256: str
    exposure_ledger_id: str
    exposure_receipt_sha256: str
    validation_attestation_sha256: str
    custodian_authority_sha256: str
    custody_consumption_receipt_sha256: str
    sha256: str

    @property
    def role(self) -> str:
        """Return the fixed locked-heldout role."""
        return "locked-heldout"


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"target calibration {label} must be lowercase SHA-256")
    return value


def _decode_common(
    value: object, *, fields: set[str], schema: str, role: str
) -> tuple[Mapping[str, object], CalibrationTarget, Literal["passed", "failed"], dict[str, str]]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail("target calibration attestation fields differ from the closed contract")
    if value["schema_version"] != schema or value["role"] != role:
        _fail("target calibration attestation schema or role differs")
    target = value["target"]
    if target not in ("callers", "length"):
        _fail("target calibration attestation target must be callers or length")
    status = value["status"]
    if status not in ("passed", "failed"):
        _fail("target calibration attestation status must be passed or failed")
    hashes = {field: _digest(value[field], field) for field in _HASH_FIELDS}
    return value, cast(CalibrationTarget, target), cast(Literal["passed", "failed"], status), hashes


def decode_target_validation_attestation(value: object) -> TargetValidationAttestation:
    """Decode a closed target-aware validation attestation.

    Args:
        value: Parsed JSON object.

    Returns:
        Immutable validation result and exact canonical digest.

    Raises:
        ValueError: If fields, roles, targets, statuses, or digests differ.
    """
    raw, target, status, hashes = _decode_common(
        value,
        fields=_COMMON_FIELDS,
        schema="calibration-target-validation-attestation-v2",
        role="validation",
    )
    return TargetValidationAttestation(
        target=target,
        status=status,
        candidate_sha256=hashes["candidate_sha256"],
        candidate_id=hashes["candidate_id"],
        study_sha256=hashes["study_sha256"],
        protocol_sha256=hashes["protocol_sha256"],
        partition_sha256=hashes["partition_sha256"],
        baseline_sha256=hashes["baseline_sha256"],
        evidence_sha256=hashes["evidence_sha256"],
        run_manifest_sha256=hashes["run_manifest_sha256"],
        metrics_sha256=hashes["metrics_sha256"],
        exposure_ledger_id=hashes["exposure_ledger_id"],
        exposure_receipt_sha256=hashes["exposure_receipt_sha256"],
        sha256=canonical_sha256(raw),
    )


def decode_target_locked_attestation(value: object) -> TargetLockedAttestation:
    """Decode a closed target-aware locked-heldout attestation.

    Args:
        value: Parsed JSON object.

    Returns:
        Immutable locked result and prior authorization bindings.

    Raises:
        ValueError: If fields, roles, targets, statuses, or digests differ.
    """
    raw, target, status, hashes = _decode_common(
        value,
        fields=_LOCKED_FIELDS,
        schema="calibration-target-locked-attestation-v2",
        role="locked-heldout",
    )
    return TargetLockedAttestation(
        target=target,
        status=status,
        candidate_sha256=hashes["candidate_sha256"],
        candidate_id=hashes["candidate_id"],
        study_sha256=hashes["study_sha256"],
        protocol_sha256=hashes["protocol_sha256"],
        partition_sha256=hashes["partition_sha256"],
        baseline_sha256=hashes["baseline_sha256"],
        evidence_sha256=hashes["evidence_sha256"],
        run_manifest_sha256=hashes["run_manifest_sha256"],
        metrics_sha256=hashes["metrics_sha256"],
        exposure_ledger_id=hashes["exposure_ledger_id"],
        exposure_receipt_sha256=hashes["exposure_receipt_sha256"],
        validation_attestation_sha256=_digest(raw["validation_attestation_sha256"], "validation_attestation_sha256"),
        custodian_authority_sha256=_digest(raw["custodian_authority_sha256"], "custodian_authority_sha256"),
        custody_consumption_receipt_sha256=_digest(
            raw["custody_consumption_receipt_sha256"], "custody_consumption_receipt_sha256"
        ),
        sha256=canonical_sha256(raw),
    )


def _validation_document(value: TargetValidationAttestation) -> dict[str, object]:
    return {
        "schema_version": "calibration-target-validation-attestation-v2",
        "target": value.target,
        "role": value.role,
        "status": value.status,
        **{field: getattr(value, field) for field in _HASH_FIELDS},
    }


def _locked_document(value: TargetLockedAttestation) -> dict[str, object]:
    return {
        "schema_version": "calibration-target-locked-attestation-v2",
        "target": value.target,
        "role": value.role,
        "status": value.status,
        **{field: getattr(value, field) for field in _HASH_FIELDS},
        "validation_attestation_sha256": value.validation_attestation_sha256,
        "custodian_authority_sha256": value.custodian_authority_sha256,
        "custody_consumption_receipt_sha256": value.custody_consumption_receipt_sha256,
    }


def target_validation_attestation_document(value: TargetValidationAttestation) -> dict[str, object]:
    """Project and revalidate one immutable validation attestation."""
    if not isinstance(value, TargetValidationAttestation):
        _fail("target validation attestation must be its decoded typed value")
    document = _validation_document(value)
    if decode_target_validation_attestation(document) != value:
        _fail("target validation attestation differs from its canonical content")
    return document


def target_locked_attestation_document(value: TargetLockedAttestation) -> dict[str, object]:
    """Project and revalidate one immutable locked-heldout attestation."""
    if not isinstance(value, TargetLockedAttestation):
        _fail("target locked attestation must be its decoded typed value")
    document = _locked_document(value)
    if decode_target_locked_attestation(document) != value:
        _fail("target locked attestation differs from its canonical content")
    return document


def encode_target_validation_attestation(**values: object) -> dict[str, object]:
    """Encode and self-validate one target-aware validation attestation."""
    document = {"schema_version": "calibration-target-validation-attestation-v2", "role": "validation", **values}
    decode_target_validation_attestation(document)
    return document


def encode_target_locked_attestation(**values: object) -> dict[str, object]:
    """Encode and self-validate one target-aware locked-heldout attestation."""
    document = {
        "schema_version": "calibration-target-locked-attestation-v2",
        "role": "locked-heldout",
        **values,
    }
    decode_target_locked_attestation(document)
    return document
