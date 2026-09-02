"""Closed lineage-rich validation attestations for custodian handoff."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass

from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_FIELDS = {
    "schema_version",
    "role",
    "status",
    "profile_sha256",
    "protocol_sha256",
    "evidence_sha256",
    "metrics_sha256",
    "study_sha256",
    "partition_manifest_sha256",
    "profile_dataset_sha256",
    "study_binding_sha256",
    "run_commitments_sha256",
    "role_run_commitments_sha256",
    "role_run_artifacts_sha256",
}


@dataclass(frozen=True)
class ValidationAttestation:
    """Passed or failed validation result with complete predeclared lineage."""

    status: str
    profile_sha256: str
    protocol_sha256: str
    evidence_sha256: str
    metrics_sha256: str
    study_sha256: str
    partition_manifest_sha256: str
    profile_dataset_sha256: str
    study_binding_sha256: str
    run_commitments_sha256: str
    role_run_commitments_sha256: str
    role_run_artifacts_sha256: str
    sha256: str

    @property
    def role(self) -> str:
        """Return the fixed role encoded by this closed schema."""
        return "validation"


def encode_validation_attestation(**values: object) -> dict[str, object]:
    """Encode and self-validate one complete validation attestation."""
    document = {"schema_version": "calibration-validation-attestation-v2", "role": "validation", **values}
    decode_validation_attestation(document)
    return document


def decode_validation_attestation(value: object) -> ValidationAttestation:
    """Decode the closed version-2 validation lineage contract."""
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        raise ValueError("calibration validation attestation fields differ")
    if value["schema_version"] != "calibration-validation-attestation-v2" or value["role"] != "validation":
        raise ValueError("calibration validation attestation schema or role differs")
    status = value["status"]
    if status not in {"passed", "failed"}:
        raise ValueError("calibration validation attestation status differs")
    hashes = {
        field: _digest(value[field], field.replace("_sha256", "").replace("_", " "))
        for field in _FIELDS
        if field.endswith("_sha256")
    }
    return ValidationAttestation(
        str(status),
        hashes["profile_sha256"],
        hashes["protocol_sha256"],
        hashes["evidence_sha256"],
        hashes["metrics_sha256"],
        hashes["study_sha256"],
        hashes["partition_manifest_sha256"],
        hashes["profile_dataset_sha256"],
        hashes["study_binding_sha256"],
        hashes["run_commitments_sha256"],
        hashes["role_run_commitments_sha256"],
        hashes["role_run_artifacts_sha256"],
        canonical_sha256(value),
    )


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"calibration validation {label} hash must be lowercase SHA-256")
    return value
