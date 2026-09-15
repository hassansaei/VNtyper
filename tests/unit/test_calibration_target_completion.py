"""Target-aware terminal custody completion contract."""

from dataclasses import replace

import pytest

from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def completion_document() -> dict[str, object]:
    names = (
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
    return {
        "schema_version": "calibration-target-completion-v2",
        "status": "completed",
        "target": "length",
        **{name: format(index, "x") * 64 for index, name in enumerate(names)},
    }


def test_target_completion_round_trips_and_hashes_exact_content() -> None:
    from vntyper.scripts.calibration_target_completion import (
        decode_target_completion,
        encode_target_completion,
        target_completion_document,
    )

    completion = decode_target_completion(completion_document())
    assert completion.sha256 == canonical_sha256(completion_document())
    assert target_completion_document(completion) == completion_document()
    values = completion_document()
    values.pop("schema_version")
    values.pop("status")
    assert encode_target_completion(**values) == completion_document()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("schema_version", "calibration-completion-v1"),
        ("status", "passed"),
        ("target", "dominance"),
        ("exposure_ledger_id", "ledger"),
        ("custody_consumption_receipt_sha256", True),
    ],
)
def test_target_completion_rejects_invalid_fields(field: str, value: object) -> None:
    from vntyper.scripts.calibration_target_completion import decode_target_completion

    document = completion_document()
    document[field] = value
    with pytest.raises(ValueError):
        decode_target_completion(document)


def test_target_completion_is_closed_and_projection_revalidates_content() -> None:
    from vntyper.scripts.calibration_target_completion import decode_target_completion, target_completion_document

    document = completion_document()
    document["retryable"] = True
    with pytest.raises(ValueError, match="fields"):
        decode_target_completion(document)
    completion = decode_target_completion(completion_document())
    with pytest.raises(ValueError, match="canonical"):
        target_completion_document(replace(completion, sha256="f" * 64))
    with pytest.raises(ValueError, match="typed"):
        target_completion_document({})  # type: ignore[arg-type]
