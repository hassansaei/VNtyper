"""Exposure declarations validate identities without granting outcome access."""

from dataclasses import replace

import pytest

from vntyper.scripts.calibration_exposure import (
    decode_exposure_identities,
    decode_exposure_receipt,
    exposure_receipt_document,
    require_digest,
)
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def _receipt():
    return {
        "schema_version": "calibration-target-exposure-receipt-v2",
        "target": "length",
        "role": "validation",
        **dict.fromkeys(
            ("study_sha256", "partition_sha256", "evidence_sha256", "membership_sha256", "exposure_ledger_id"), "a" * 64
        ),
        "sequence": 1,
    }


@pytest.mark.parametrize("bad", [None, {}, "A" * 64, "a" * 63, "a" * 64 + "\n", True])
def test_digest_rejects_noncanonical_tokens(bad):
    with pytest.raises(ValueError):
        require_digest(bad, "example")


@pytest.mark.parametrize(
    "bad",
    [
        None,
        {},
        [],
        [None],
        [{"namespace": "specimen", "sha256": "a" * 64, "extra": 1}],
        [{"namespace": [], "sha256": "a" * 64}],
        [{"namespace": "backbone", "sha256": "a" * 64}],
        [{"namespace": "specimen", "sha256": "a" * 64}],
        [{"namespace": "specimen", "sha256": "a" * 64}, {"namespace": "backbone", "sha256": "b" * 64}],
        [{"namespace": "backbone", "sha256": "b" * 64}] * 2,
    ],
)
def test_identity_set_rejects_missing_unsorted_and_duplicate_tokens(bad):
    with pytest.raises(ValueError):
        decode_exposure_identities(bad)


def test_identity_set_preserves_all_audited_namespaces_and_multiple_families():
    rows = [
        {"namespace": name, "sha256": "a" * 64}
        for name in (
            "backbone",
            "family",
            "named-readset",
            "pair",
            "physical-readset",
            "seed-family",
            "specimen",
            "unnamed-readset",
        )
    ]
    assert decode_exposure_identities(rows) == tuple((row["namespace"], row["sha256"]) for row in rows)


@pytest.mark.parametrize(
    "changes",
    [
        {"schema_version": "v1"},
        {"extra": 1},
        {"target": []},
        {"role": {}},
        {"sequence": True},
        {"sequence": 0},
        {"sequence": 1.0},
        {"sequence": 2**53},
        {"membership_sha256": "bad"},
        {"exposure_ledger_id": None},
    ],
)
def test_receipt_rejects_wrong_schema_identity_and_sequence(changes):
    with pytest.raises(ValueError):
        decode_exposure_receipt({**_receipt(), **changes})


def test_receipt_projection_revalidates_derived_identity_and_type():
    document = _receipt()
    receipt = decode_exposure_receipt(document)
    assert receipt.sha256 == canonical_sha256(document)
    assert exposure_receipt_document(receipt) == document
    for changed in (replace(receipt, sequence=2), replace(receipt, sha256="0" * 64), document):
        with pytest.raises(ValueError):
            exposure_receipt_document(changed)
    with pytest.raises(ValueError):
        decode_exposure_receipt(None)
