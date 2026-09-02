"""Strict bindings for a separately supplied external-custodian import."""

from __future__ import annotations

import hashlib
import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.calibration_secure_io import SecureDirectoryReader
from vntyper.scripts.calibration_study_binding import StudyBinding, decode_study_binding
from vntyper.scripts.calibration_validation_attestation import (
    ValidationAttestation,
    decode_validation_attestation,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_FIELDS = {
    "schema_version",
    "authority_kind",
    "custodian_name",
    "attestation_id",
    "role",
    "status",
    "study_sha256",
    "protocol_sha256",
    "partition_manifest_sha256",
    "profile_sha256",
    "profile_dataset_sha256",
    "locked_payload_sha256",
    "locked_dataset_sha256",
    "locked_run_hashes_sha256",
    "validation_attestation_sha256",
    "validation_evidence_sha256",
    "study_binding_sha256",
    "run_commitments_sha256",
    "validation_role_run_commitments_sha256",
    "validation_role_run_artifacts_sha256",
    "locked_role_run_commitments_sha256",
    "locked_role_run_artifacts_sha256",
}


@dataclass(frozen=True)
class CustodianAuthority:
    """Named external authority assertion and its complete artifact bindings."""

    custodian_name: str
    attestation_id: str
    study_sha256: str
    protocol_sha256: str
    partition_manifest_sha256: str
    profile_sha256: str
    profile_dataset_sha256: str
    locked_payload_sha256: str
    locked_dataset_sha256: str
    locked_run_hashes_sha256: str
    validation_attestation_sha256: str
    validation_evidence_sha256: str
    study_binding_sha256: str
    run_commitments_sha256: str
    validation_role_run_commitments_sha256: str
    validation_role_run_artifacts_sha256: str
    locked_role_run_commitments_sha256: str
    locked_role_run_artifacts_sha256: str


@dataclass
class CustodianImport:
    """Pinned custodian header and deferred no-follow locked payload reader."""

    authority: CustodianAuthority
    validation: ValidationAttestation
    study_binding: StudyBinding
    reader: SecureDirectoryReader

    def read_locked_payload(self) -> bytes:
        """Read locked bytes once from the pinned import directory after precommit."""
        return self.reader.read_file("locked_payload.json")

    def close(self) -> None:
        """Close the pinned import directory."""
        self.reader.close()

    def __enter__(self) -> CustodianImport:
        return self

    def __exit__(self, _type, _value, _traceback) -> None:
        self.close()


def load_custodian_import_header(evidence_path: Path) -> CustodianImport:
    """Load and bind only non-locked custodian header bytes before precommit.

    Args:
        evidence_path: Separately supplied custodian-import directory.

    Returns:
        Decoded authority and passed-validation candidates. Status is checked by
        the operation after it binds the explicit profile.

    Raises:
        ValueError: If inventory, file types, JSON, or checksums differ.
    """
    expected = {
        "authority_attestation.json",
        "validation_attestation.json",
        "locked_payload.json",
        "checksums.json",
        "study_binding.json",
    }
    try:
        reader = SecureDirectoryReader.open(evidence_path, expected)
    except ValueError as error:
        raise ValueError(f"custodian import is invalid: {error}") from error
    try:
        raw = reader.read_files(
            ("authority_attestation.json", "validation_attestation.json", "study_binding.json", "checksums.json")
        )
        authority = decode_custodian_authority(load_strict_json_object(raw["authority_attestation.json"]))
        validation_document = load_strict_json_object(raw["validation_attestation.json"])
        if canonical_json_bytes(validation_document) != raw["validation_attestation.json"]:
            raise ValueError("custodian import validation attestation must be canonical installed bytes")
        validation = decode_validation_attestation(validation_document)
        binding = decode_study_binding(load_strict_json_object(raw["study_binding.json"]))
        validation_sha256 = hashlib.sha256(raw["validation_attestation.json"]).hexdigest()
        if authority.validation_attestation_sha256 != validation_sha256:
            raise ValueError("custodian authority does not bind the exact validation-attestation bytes")
        if authority.study_binding_sha256 != hashlib.sha256(raw["study_binding.json"]).hexdigest():
            raise ValueError("custodian authority does not bind the exact study-binding bytes")
        checksums = load_strict_json_object(raw["checksums.json"])
        files = _mapping(checksums.get("files"), "custodian import checksum files")
        expected_files = expected - {"checksums.json"}
        if checksums.get("schema_version") != "calibration-checksums-v1" or set(files) != expected_files:
            raise ValueError("custodian import checksum fields differ")
        observed = {
            "authority_attestation.json": hashlib.sha256(raw["authority_attestation.json"]).hexdigest(),
            "validation_attestation.json": validation_sha256,
            "study_binding.json": hashlib.sha256(raw["study_binding.json"]).hexdigest(),
            "locked_payload.json": authority.locked_payload_sha256,
        }
        if files != observed:
            raise ValueError("custodian import checksum bindings differ")
        return CustodianImport(authority, validation, binding, reader)
    except BaseException:
        reader.close()
        raise


def decode_custodian_authority(value: object) -> CustodianAuthority:
    """Decode one closed named external-custodian authority attestation.

    Args:
        value: Candidate JSON object supplied separately by the custodian.

    Returns:
        The strict immutable authority assertion.

    Raises:
        ValueError: If fields, authority identity, role, status, or hashes differ.

    Notes:
        This decoder verifies exact local bytes and bindings. Without a separately
        configured cryptographic trust root it does not prove external independence.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"calibration custodian authority fields differ: expected {sorted(_FIELDS)}, got {actual}")
    if value["schema_version"] != "calibration-custodian-authority-v1":
        raise ValueError("calibration custodian authority schema version is unsupported")
    if value["authority_kind"] != "named-external-custodian-attestation":
        raise ValueError("calibration authority kind must be named-external-custodian-attestation")
    if value["role"] != "locked-heldout":
        raise ValueError("calibration custodian authority role must be locked-heldout")
    if value["status"] != "authorized":
        raise ValueError("calibration custodian authority status must be authorized")
    custodian = _text(value["custodian_name"], "custodian name")
    attestation_id = _text(value["attestation_id"], "attestation id")
    return CustodianAuthority(
        custodian,
        attestation_id,
        _digest(value["study_sha256"], "study"),
        _digest(value["protocol_sha256"], "protocol"),
        _digest(value["partition_manifest_sha256"], "partition manifest"),
        _digest(value["profile_sha256"], "profile"),
        _digest(value["profile_dataset_sha256"], "profile dataset"),
        _digest(value["locked_payload_sha256"], "locked payload"),
        _digest(value["locked_dataset_sha256"], "locked dataset"),
        _digest(value["locked_run_hashes_sha256"], "locked run hashes"),
        _digest(value["validation_attestation_sha256"], "validation attestation"),
        _digest(value["validation_evidence_sha256"], "validation evidence"),
        _digest(value["study_binding_sha256"], "study binding"),
        _digest(value["run_commitments_sha256"], "run commitments"),
        _digest(value["validation_role_run_commitments_sha256"], "validation role run commitments"),
        _digest(value["validation_role_run_artifacts_sha256"], "validation role run artifacts"),
        _digest(value["locked_role_run_commitments_sha256"], "locked role run commitments"),
        _digest(value["locked_role_run_artifacts_sha256"], "locked role run artifacts"),
    )


def validate_custodian_authority_bindings(
    authority: CustodianAuthority,
    *,
    study_sha256: str,
    protocol_sha256: str,
    partition_manifest_sha256: str,
    profile_sha256: str,
    profile_dataset_sha256: str,
    locked_payload_sha256: str,
    locked_dataset_sha256: str,
    locked_run_hashes_sha256: str,
    validation_attestation_sha256: str,
    validation_evidence_sha256: str,
    study_binding_sha256: str,
    run_commitments_sha256: str,
    validation_role_run_commitments_sha256: str,
    validation_role_run_artifacts_sha256: str,
    locked_role_run_commitments_sha256: str,
    locked_role_run_artifacts_sha256: str,
) -> None:
    """Require the authority assertion to bind every opened artifact exactly.

    Args:
        authority: Decoded external authority assertion.
        study_sha256: Opened complete study digest.
        protocol_sha256: Opened protocol digest.
        partition_manifest_sha256: Opened partition digest.
        profile_sha256: Explicit evaluated profile digest.
        profile_dataset_sha256: Fitted dataset commitment in profile metadata.
        locked_payload_sha256: Exact opened payload-byte digest.
        locked_dataset_sha256: Decoded locked dataset digest.
        locked_run_hashes_sha256: Decoded locked run-artifact commitment digest.
        validation_attestation_sha256: Exact passed validation-attestation digest.

    Raises:
        ValueError: If any binding differs or a supplied digest is malformed.
    """
    if not isinstance(authority, CustodianAuthority):
        raise ValueError("calibration custodian authority bindings require a decoded authority")
    expected = tuple(
        _digest(value, label)
        for value, label in (
            (study_sha256, "opened study"),
            (protocol_sha256, "opened protocol"),
            (partition_manifest_sha256, "opened partition manifest"),
            (profile_sha256, "opened profile"),
            (profile_dataset_sha256, "opened profile dataset"),
            (locked_payload_sha256, "opened locked payload"),
            (locked_dataset_sha256, "opened locked dataset"),
            (locked_run_hashes_sha256, "opened locked run hashes"),
            (validation_attestation_sha256, "opened validation attestation"),
            (validation_evidence_sha256, "opened validation evidence"),
            (study_binding_sha256, "opened study binding"),
            (run_commitments_sha256, "opened run commitments"),
            (validation_role_run_commitments_sha256, "opened validation role run commitments"),
            (validation_role_run_artifacts_sha256, "opened validation role run artifacts"),
            (locked_role_run_commitments_sha256, "opened locked role run commitments"),
            (locked_role_run_artifacts_sha256, "opened locked role run artifacts"),
        )
    )
    observed = (
        authority.study_sha256,
        authority.protocol_sha256,
        authority.partition_manifest_sha256,
        authority.profile_sha256,
        authority.profile_dataset_sha256,
        authority.locked_payload_sha256,
        authority.locked_dataset_sha256,
        authority.locked_run_hashes_sha256,
        authority.validation_attestation_sha256,
        authority.validation_evidence_sha256,
        authority.study_binding_sha256,
        authority.run_commitments_sha256,
        authority.validation_role_run_commitments_sha256,
        authority.validation_role_run_artifacts_sha256,
        authority.locked_role_run_commitments_sha256,
        authority.locked_role_run_artifacts_sha256,
    )
    if observed != expected:
        raise ValueError("calibration custodian authority protocol or artifact bindings differ from opened artifacts")


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"calibration {label} must be a non-empty string")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"calibration {label} hash must be lowercase SHA-256")
    return value


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be an object")
    return value
