"""Closed target-aware v2 successful custody terminal."""

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
_HASH_FIELDS = (
    "candidate_sha256",
    "candidate_id",
    "study_sha256",
    "protocol_sha256",
    "partition_sha256",
    "baseline_sha256",
    "validation_attestation_sha256",
    "locked_heldout_attestation_sha256",
    "custodian_authority_sha256",
    "locked_heldout_evidence_sha256",
    "locked_heldout_run_manifest_sha256",
    "exposure_ledger_id",
    "validation_exposure_receipt_sha256",
    "locked_heldout_exposure_receipt_sha256",
    "custody_consumption_receipt_sha256",
)
_FIELDS = {"schema_version", "status", "target", *_HASH_FIELDS}


@dataclass(frozen=True)
class TargetCompletion:
    """Terminal success binding authorization, exposures, and consumption."""

    target: CalibrationTarget
    candidate_sha256: str
    candidate_id: str
    study_sha256: str
    protocol_sha256: str
    partition_sha256: str
    baseline_sha256: str
    validation_attestation_sha256: str
    locked_heldout_attestation_sha256: str
    custodian_authority_sha256: str
    locked_heldout_evidence_sha256: str
    locked_heldout_run_manifest_sha256: str
    exposure_ledger_id: str
    validation_exposure_receipt_sha256: str
    locked_heldout_exposure_receipt_sha256: str
    custody_consumption_receipt_sha256: str
    sha256: str

    @property
    def status(self) -> str:
        """Return the fixed successful terminal state."""
        return "completed"


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"target completion {label} must be lowercase SHA-256")
    return value


def decode_target_completion(value: object) -> TargetCompletion:
    """Decode one complete target-aware successful terminal.

    Args:
        value: Parsed JSON completion object.

    Returns:
        Immutable completion and canonical content digest.

    Raises:
        ValueError: If fields, status, target, or identities differ.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        _fail("target completion fields differ from the closed contract")
    if value["schema_version"] != "calibration-target-completion-v2" or value["status"] != "completed":
        _fail("target completion schema or status differs")
    target = value["target"]
    if target not in ("callers", "length"):
        _fail("target completion target must be callers or length")
    return TargetCompletion(
        target=cast(CalibrationTarget, target),
        candidate_sha256=_digest(value["candidate_sha256"], "candidate_sha256"),
        candidate_id=_digest(value["candidate_id"], "candidate_id"),
        study_sha256=_digest(value["study_sha256"], "study_sha256"),
        protocol_sha256=_digest(value["protocol_sha256"], "protocol_sha256"),
        partition_sha256=_digest(value["partition_sha256"], "partition_sha256"),
        baseline_sha256=_digest(value["baseline_sha256"], "baseline_sha256"),
        validation_attestation_sha256=_digest(value["validation_attestation_sha256"], "validation_attestation_sha256"),
        locked_heldout_attestation_sha256=_digest(
            value["locked_heldout_attestation_sha256"], "locked_heldout_attestation_sha256"
        ),
        custodian_authority_sha256=_digest(value["custodian_authority_sha256"], "custodian_authority_sha256"),
        locked_heldout_evidence_sha256=_digest(
            value["locked_heldout_evidence_sha256"], "locked_heldout_evidence_sha256"
        ),
        locked_heldout_run_manifest_sha256=_digest(
            value["locked_heldout_run_manifest_sha256"], "locked_heldout_run_manifest_sha256"
        ),
        exposure_ledger_id=_digest(value["exposure_ledger_id"], "exposure_ledger_id"),
        validation_exposure_receipt_sha256=_digest(
            value["validation_exposure_receipt_sha256"], "validation_exposure_receipt_sha256"
        ),
        locked_heldout_exposure_receipt_sha256=_digest(
            value["locked_heldout_exposure_receipt_sha256"], "locked_heldout_exposure_receipt_sha256"
        ),
        custody_consumption_receipt_sha256=_digest(
            value["custody_consumption_receipt_sha256"], "custody_consumption_receipt_sha256"
        ),
        sha256=canonical_sha256(value),
    )


def _document(value: TargetCompletion) -> dict[str, object]:
    return {
        "schema_version": "calibration-target-completion-v2",
        "status": value.status,
        "target": value.target,
        **{field: getattr(value, field) for field in _HASH_FIELDS},
    }


def target_completion_document(value: TargetCompletion) -> dict[str, object]:
    """Project and revalidate one immutable target completion."""
    if not isinstance(value, TargetCompletion):
        _fail("target completion must be its decoded typed value")
    document = _document(value)
    if decode_target_completion(document) != value:
        _fail("target completion differs from its canonical content")
    return document


def encode_target_completion(**values: object) -> dict[str, object]:
    """Encode and self-validate one target-aware successful terminal."""
    document = {"schema_version": "calibration-target-completion-v2", "status": "completed", **values}
    decode_target_completion(document)
    return document
