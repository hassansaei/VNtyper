"""Closed target-aware v2 external-custodian authority contract."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NoReturn, cast

from vntyper.scripts.calibration_candidate import CalibrationTarget
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_FIELDS = {
    "schema_version",
    "authority_kind",
    "custodian_name",
    "attestation_id",
    "status",
    "role",
    "target",
    "candidate_sha256",
    "candidate_id",
    "study_sha256",
    "protocol_sha256",
    "partition_sha256",
    "baseline_sha256",
    "validation_evidence_sha256",
    "locked_heldout_evidence_sha256",
    "validation_run_manifest_sha256",
    "locked_heldout_run_manifest_sha256",
    "locked_payload_sha256",
    "validation_attestation_sha256",
    "exposure_ledger_id",
    "validation_exposure_receipt_sha256",
}
_HASH_FIELDS = (
    "candidate_sha256",
    "candidate_id",
    "study_sha256",
    "protocol_sha256",
    "partition_sha256",
    "baseline_sha256",
    "validation_evidence_sha256",
    "locked_heldout_evidence_sha256",
    "validation_run_manifest_sha256",
    "locked_heldout_run_manifest_sha256",
    "locked_payload_sha256",
    "validation_attestation_sha256",
    "exposure_ledger_id",
    "validation_exposure_receipt_sha256",
)


@dataclass(frozen=True)
class TargetCustodianAuthority:
    """External authorization of one exact target/candidate/locked payload."""

    custodian_name: str
    attestation_id: str
    target: CalibrationTarget
    candidate_sha256: str
    candidate_id: str
    study_sha256: str
    protocol_sha256: str
    partition_sha256: str
    baseline_sha256: str
    validation_evidence_sha256: str
    locked_heldout_evidence_sha256: str
    validation_run_manifest_sha256: str
    locked_heldout_run_manifest_sha256: str
    locked_payload_sha256: str
    validation_attestation_sha256: str
    exposure_ledger_id: str
    validation_exposure_receipt_sha256: str
    sha256: str

    @property
    def authority_kind(self) -> str:
        """Return the fixed external authority kind."""
        return "external-custodian"

    @property
    def status(self) -> str:
        """Return the fixed authorization state."""
        return "authorized"

    @property
    def role(self) -> str:
        """Return the fixed locked-heldout role."""
        return "locked-heldout"


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"target custodian {label} must be lowercase SHA-256")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        _fail(f"target custodian {label} must be non-empty text without surrounding whitespace")
    return value


def decode_target_custodian_authority(value: object) -> TargetCustodianAuthority:
    """Decode one complete target-aware external authority document.

    This validates content and internal identity. It does not establish that the
    named custodian is independent without an external trust boundary.

    Args:
        value: Parsed JSON authority object.

    Returns:
        Immutable authority and canonical content digest.

    Raises:
        ValueError: If the document is open, malformed, or not authorized.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        _fail("target custodian authority fields differ from the closed contract")
    if value["schema_version"] != "calibration-target-custodian-authority-v2":
        _fail("target custodian authority schema is unsupported")
    if (
        value["authority_kind"] != "external-custodian"
        or value["status"] != "authorized"
        or value["role"] != "locked-heldout"
    ):
        _fail("target custodian authority kind, status, or role differs")
    target = value["target"]
    if target not in ("callers", "length"):
        _fail("target custodian authority target must be callers or length")
    return TargetCustodianAuthority(
        custodian_name=_text(value["custodian_name"], "custodian_name"),
        attestation_id=_text(value["attestation_id"], "attestation_id"),
        target=cast(CalibrationTarget, target),
        candidate_sha256=_digest(value["candidate_sha256"], "candidate_sha256"),
        candidate_id=_digest(value["candidate_id"], "candidate_id"),
        study_sha256=_digest(value["study_sha256"], "study_sha256"),
        protocol_sha256=_digest(value["protocol_sha256"], "protocol_sha256"),
        partition_sha256=_digest(value["partition_sha256"], "partition_sha256"),
        baseline_sha256=_digest(value["baseline_sha256"], "baseline_sha256"),
        validation_evidence_sha256=_digest(value["validation_evidence_sha256"], "validation_evidence_sha256"),
        locked_heldout_evidence_sha256=_digest(
            value["locked_heldout_evidence_sha256"], "locked_heldout_evidence_sha256"
        ),
        validation_run_manifest_sha256=_digest(
            value["validation_run_manifest_sha256"], "validation_run_manifest_sha256"
        ),
        locked_heldout_run_manifest_sha256=_digest(
            value["locked_heldout_run_manifest_sha256"], "locked_heldout_run_manifest_sha256"
        ),
        locked_payload_sha256=_digest(value["locked_payload_sha256"], "locked_payload_sha256"),
        validation_attestation_sha256=_digest(value["validation_attestation_sha256"], "validation_attestation_sha256"),
        exposure_ledger_id=_digest(value["exposure_ledger_id"], "exposure_ledger_id"),
        validation_exposure_receipt_sha256=_digest(
            value["validation_exposure_receipt_sha256"], "validation_exposure_receipt_sha256"
        ),
        sha256=canonical_sha256(value),
    )


def _document(value: TargetCustodianAuthority) -> dict[str, object]:
    return {
        "schema_version": "calibration-target-custodian-authority-v2",
        "authority_kind": value.authority_kind,
        "custodian_name": value.custodian_name,
        "attestation_id": value.attestation_id,
        "status": value.status,
        "role": value.role,
        "target": value.target,
        **{field: getattr(value, field) for field in _HASH_FIELDS},
    }


def target_custodian_authority_document(value: TargetCustodianAuthority) -> dict[str, object]:
    """Project and revalidate one immutable target authority."""
    if not isinstance(value, TargetCustodianAuthority):
        _fail("target custodian authority must be its decoded typed value")
    document = _document(value)
    if decode_target_custodian_authority(document) != value:
        _fail("target custodian authority differs from its canonical content")
    return document


def encode_target_custodian_authority(**values: object) -> dict[str, object]:
    """Encode and self-validate one target-aware external authority."""
    document = {
        "schema_version": "calibration-target-custodian-authority-v2",
        "authority_kind": "external-custodian",
        "status": "authorized",
        "role": "locked-heldout",
        **values,
    }
    decode_target_custodian_authority(document)
    return document
