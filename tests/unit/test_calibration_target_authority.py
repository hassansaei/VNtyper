"""Target-aware external-custodian authority contract."""

from dataclasses import replace

import pytest

from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def authority_document() -> dict[str, object]:
    fields = {
        "candidate_sha256": "0" * 64,
        "candidate_id": "1" * 64,
        "study_sha256": "2" * 64,
        "protocol_sha256": "3" * 64,
        "partition_sha256": "4" * 64,
        "baseline_sha256": "5" * 64,
        "validation_evidence_sha256": "6" * 64,
        "locked_heldout_evidence_sha256": "7" * 64,
        "validation_run_manifest_sha256": "8" * 64,
        "locked_heldout_run_manifest_sha256": "9" * 64,
        "locked_payload_sha256": "a" * 64,
        "validation_attestation_sha256": "b" * 64,
        "exposure_ledger_id": "c" * 64,
        "validation_exposure_receipt_sha256": "d" * 64,
    }
    return {
        "schema_version": "calibration-target-custodian-authority-v2",
        "authority_kind": "external-custodian",
        "custodian_name": "Independent Example Repository",
        "attestation_id": "IER-332-0001",
        "status": "authorized",
        "role": "locked-heldout",
        "target": "length",
        **fields,
    }


def test_target_authority_round_trips_complete_lineage() -> None:
    from vntyper.scripts.calibration_target_authority import (
        decode_target_custodian_authority,
        encode_target_custodian_authority,
        target_custodian_authority_document,
    )

    authority = decode_target_custodian_authority(authority_document())
    assert authority.sha256 == canonical_sha256(authority_document())
    assert authority.locked_payload_sha256 == "a" * 64
    assert target_custodian_authority_document(authority) == authority_document()
    values = authority_document()
    for field in ("schema_version", "authority_kind", "status", "role"):
        values.pop(field)
    assert encode_target_custodian_authority(**values) == authority_document()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("authority_kind", "local"),
        ("status", "passed"),
        ("role", "validation"),
        ("target", "dominance"),
        ("custodian_name", ""),
        ("attestation_id", " IER "),
        ("exposure_ledger_id", "f" * 63),
        ("locked_payload_sha256", "F" * 64),
    ],
)
def test_target_authority_rejects_invalid_identity_or_binding(field: str, value: object) -> None:
    from vntyper.scripts.calibration_target_authority import decode_target_custodian_authority

    document = authority_document()
    document[field] = value
    with pytest.raises(ValueError):
        decode_target_custodian_authority(document)


def test_target_authority_is_closed_and_projection_revalidates_content() -> None:
    from vntyper.scripts.calibration_target_authority import (
        decode_target_custodian_authority,
        target_custodian_authority_document,
    )

    document = authority_document()
    document["local_custody_is_independent_proof"] = True
    with pytest.raises(ValueError, match="fields"):
        decode_target_custodian_authority(document)
    authority = decode_target_custodian_authority(authority_document())
    with pytest.raises(ValueError, match="canonical"):
        target_custodian_authority_document(replace(authority, validation_attestation_sha256="f" * 64))
    with pytest.raises(ValueError, match="typed"):
        target_custodian_authority_document({})  # type: ignore[arg-type]
