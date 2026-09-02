"""Closed lineage-rich validation attestation contracts."""

import pytest

from vntyper.scripts.calibration_validation_attestation import (
    decode_validation_attestation,
    encode_validation_attestation,
)

pytestmark = pytest.mark.unit


def _document() -> dict[str, object]:
    return encode_validation_attestation(
        status="passed",
        profile_sha256="1" * 64,
        protocol_sha256="2" * 64,
        evidence_sha256="3" * 64,
        metrics_sha256="4" * 64,
        study_sha256="5" * 64,
        partition_manifest_sha256="6" * 64,
        profile_dataset_sha256="7" * 64,
        study_binding_sha256="8" * 64,
        run_commitments_sha256="9" * 64,
        role_run_commitments_sha256="a" * 64,
        role_run_artifacts_sha256="b" * 64,
    )


def test_validation_attestation_round_trips_complete_lineage() -> None:
    decoded = decode_validation_attestation(_document())

    assert decoded.role == "validation"
    assert decoded.status == "passed"
    assert decoded.study_binding_sha256 == "8" * 64
    assert decoded.role_run_artifacts_sha256 == "b" * 64


@pytest.mark.parametrize("mutation", ["missing", "unknown", "bad-hash", "bad-status"])
def test_validation_attestation_rejects_open_or_malformed_documents(mutation: str) -> None:
    document = _document()
    if mutation == "missing":
        document.pop("study_sha256")
    elif mutation == "unknown":
        document["developer_assertion"] = True
    elif mutation == "bad-hash":
        document["evidence_sha256"] = "A" * 64
    else:
        document["status"] = "authorized"

    with pytest.raises(ValueError):
        decode_validation_attestation(document)
