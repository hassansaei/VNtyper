"""Target-aware validation and locked-heldout attestation contracts."""

from dataclasses import replace

import pytest

from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def validation_document() -> dict[str, object]:
    return {
        "schema_version": "calibration-target-validation-attestation-v2",
        "target": "length",
        "role": "validation",
        "status": "passed",
        "candidate_sha256": "0" * 64,
        "candidate_id": "1" * 64,
        "study_sha256": "2" * 64,
        "protocol_sha256": "3" * 64,
        "partition_sha256": "4" * 64,
        "baseline_sha256": "5" * 64,
        "evidence_sha256": "6" * 64,
        "run_manifest_sha256": "7" * 64,
        "metrics_sha256": "8" * 64,
        "exposure_ledger_id": "9" * 64,
        "exposure_receipt_sha256": "a" * 64,
    }


def locked_document() -> dict[str, object]:
    raw = validation_document()
    raw.update(
        {
            "schema_version": "calibration-target-locked-attestation-v2",
            "role": "locked-heldout",
            "evidence_sha256": "b" * 64,
            "run_manifest_sha256": "c" * 64,
            "metrics_sha256": "d" * 64,
            "exposure_receipt_sha256": "e" * 64,
            "validation_attestation_sha256": "f" * 64,
            "custodian_authority_sha256": "0" * 64,
            "custody_consumption_receipt_sha256": "1" * 64,
        }
    )
    return raw


def test_target_attestations_round_trip_and_hash_exact_content() -> None:
    from vntyper.scripts.calibration_target_attestation import (
        decode_target_locked_attestation,
        decode_target_validation_attestation,
        encode_target_locked_attestation,
        encode_target_validation_attestation,
        target_locked_attestation_document,
        target_validation_attestation_document,
    )

    validation = decode_target_validation_attestation(validation_document())
    locked = decode_target_locked_attestation(locked_document())

    assert validation.sha256 == canonical_sha256(validation_document())
    assert locked.sha256 == canonical_sha256(locked_document())
    assert target_validation_attestation_document(validation) == validation_document()
    assert target_locked_attestation_document(locked) == locked_document()
    validation_values = validation_document()
    validation_values.pop("schema_version")
    validation_values.pop("role")
    locked_values = locked_document()
    locked_values.pop("schema_version")
    locked_values.pop("role")
    assert encode_target_validation_attestation(**validation_values) == validation_document()
    assert encode_target_locked_attestation(**locked_values) == locked_document()


@pytest.mark.parametrize("kind", ["validation", "locked"])
@pytest.mark.parametrize("mutation", ["missing", "extra", "wrong-target", "wrong-role", "bad-ledger", "bad-status"])
def test_target_attestations_reject_open_or_invalid_documents(kind: str, mutation: str) -> None:
    from vntyper.scripts.calibration_target_attestation import (
        decode_target_locked_attestation,
        decode_target_validation_attestation,
    )

    document = validation_document() if kind == "validation" else locked_document()
    if mutation == "missing":
        document.pop("baseline_sha256")
    elif mutation == "extra":
        document["approved"] = True
    elif mutation == "wrong-target":
        document["target"] = "dominance"
    elif mutation == "wrong-role":
        document["role"] = "training"
    elif mutation == "bad-ledger":
        document["exposure_ledger_id"] = "ledger/path"
    else:
        document["status"] = "authorized"
    decoder = decode_target_validation_attestation if kind == "validation" else decode_target_locked_attestation
    with pytest.raises(ValueError):
        decoder(document)


def test_target_attestation_projection_revalidates_direct_construction() -> None:
    from vntyper.scripts.calibration_target_attestation import (
        decode_target_locked_attestation,
        decode_target_validation_attestation,
        target_locked_attestation_document,
        target_validation_attestation_document,
    )

    validation = decode_target_validation_attestation(validation_document())
    locked = decode_target_locked_attestation(locked_document())
    with pytest.raises(ValueError, match="canonical"):
        target_validation_attestation_document(replace(validation, sha256="f" * 64))
    with pytest.raises(ValueError, match="canonical"):
        target_locked_attestation_document(replace(locked, candidate_id="f" * 64))
    with pytest.raises(ValueError, match="typed"):
        target_validation_attestation_document({})  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="typed"):
        target_locked_attestation_document({})  # type: ignore[arg-type]
