"""Strict external-custodian authority assertions and complete bindings."""

import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from vntyper.scripts.calibration_artifact_io import write_checksums, write_json
from vntyper.scripts.calibration_custodian_import import (
    decode_custodian_authority,
    load_custodian_import_header,
    validate_custodian_authority_bindings,
)
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def _document() -> dict[str, object]:
    return {
        "schema_version": "calibration-custodian-authority-v1",
        "authority_kind": "named-external-custodian-attestation",
        "custodian_name": "Independent Example Repository",
        "attestation_id": "IER-295-0001",
        "role": "locked-heldout",
        "status": "authorized",
        "study_sha256": "a" * 64,
        "protocol_sha256": "b" * 64,
        "partition_manifest_sha256": "c" * 64,
        "profile_sha256": "d" * 64,
        "profile_dataset_sha256": "e" * 64,
        "locked_payload_sha256": "f" * 64,
        "locked_dataset_sha256": "1" * 64,
        "locked_run_hashes_sha256": "2" * 64,
        "validation_attestation_sha256": "3" * 64,
        "validation_evidence_sha256": "4" * 64,
        "study_binding_sha256": "5" * 64,
        "run_commitments_sha256": "6" * 64,
        "validation_role_run_commitments_sha256": "7" * 64,
        "validation_role_run_artifacts_sha256": "8" * 64,
        "locked_role_run_commitments_sha256": "9" * 64,
        "locked_role_run_artifacts_sha256": "0" * 64,
    }


def test_external_authority_round_trips_and_validates_every_binding() -> None:
    authority = decode_custodian_authority(_document())

    validate_custodian_authority_bindings(
        authority,
        study_sha256="a" * 64,
        protocol_sha256="b" * 64,
        partition_manifest_sha256="c" * 64,
        profile_sha256="d" * 64,
        profile_dataset_sha256="e" * 64,
        locked_payload_sha256="f" * 64,
        locked_dataset_sha256="1" * 64,
        locked_run_hashes_sha256="2" * 64,
        validation_attestation_sha256="3" * 64,
        validation_evidence_sha256="4" * 64,
        study_binding_sha256="5" * 64,
        run_commitments_sha256="6" * 64,
        validation_role_run_commitments_sha256="7" * 64,
        validation_role_run_artifacts_sha256="8" * 64,
        locked_role_run_commitments_sha256="9" * 64,
        locked_role_run_artifacts_sha256="0" * 64,
    )
    assert authority.custodian_name == "Independent Example Repository"


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("custodian_name", "", "custodian"),
        ("authority_kind", "local", "authority"),
        ("status", "passed", "authorized"),
        ("role", "validation", "locked-heldout"),
        ("profile_sha256", "A" * 64, "SHA-256"),
    ],
)
def test_external_authority_rejects_invalid_identity_fields(field: str, value: object, message: str) -> None:
    document = _document()
    document[field] = value
    with pytest.raises(ValueError, match=message):
        decode_custodian_authority(document)


def test_external_authority_rejects_unknown_fields_and_each_mismatched_binding() -> None:
    document = _document()
    document["local_custody_is_independent_proof"] = True
    with pytest.raises(ValueError, match="fields"):
        decode_custodian_authority(document)

    authority = decode_custodian_authority(_document())
    with pytest.raises(ValueError, match="bindings"):
        validate_custodian_authority_bindings(
            replace(authority, profile_sha256="9" * 64),
            study_sha256="a" * 64,
            protocol_sha256="b" * 64,
            partition_manifest_sha256="c" * 64,
            profile_sha256="d" * 64,
            profile_dataset_sha256="e" * 64,
            locked_payload_sha256="f" * 64,
            locked_dataset_sha256="1" * 64,
            locked_run_hashes_sha256="2" * 64,
            validation_attestation_sha256="3" * 64,
            validation_evidence_sha256="4" * 64,
            study_binding_sha256="5" * 64,
            run_commitments_sha256="6" * 64,
            validation_role_run_commitments_sha256="7" * 64,
            validation_role_run_artifacts_sha256="8" * 64,
            locked_role_run_commitments_sha256="9" * 64,
            locked_role_run_artifacts_sha256="0" * 64,
        )


def test_custodian_header_requires_authority_to_bind_exact_validation_bytes(tmp_path: Path) -> None:
    validation = {
        "schema_version": "calibration-attestation-v1",
        "role": "validation",
        "status": "passed",
        "profile_sha256": "d" * 64,
        "protocol_sha256": "b" * 64,
        "evidence_sha256": "4" * 64,
        "metrics_sha256": "5" * 64,
    }
    authority = _document()
    authority["locked_payload_sha256"] = hashlib.sha256(b"locked").hexdigest()
    authority["validation_attestation_sha256"] = canonical_sha256(validation)
    write_json(tmp_path / "authority_attestation.json", authority)
    write_json(tmp_path / "validation_attestation.json", validation)
    with (tmp_path / "validation_attestation.json").open("ab") as stream:
        stream.write(b" ")
    authority["validation_attestation_sha256"] = hashlib.sha256(
        (tmp_path / "validation_attestation.json").read_bytes()
    ).hexdigest()
    write_json(tmp_path / "authority_attestation.json", authority)
    (tmp_path / "locked_payload.json").write_bytes(b"locked")
    (tmp_path / "study_binding.json").write_bytes(b"{}\n")
    write_checksums(tmp_path)

    with pytest.raises(ValueError, match="canonical.*validation|validation.*canonical"):
        load_custodian_import_header(tmp_path)
