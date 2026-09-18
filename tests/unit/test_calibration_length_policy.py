"""Unit tests for pure length calibration policy, measurement policy, and source truth validation."""

from __future__ import annotations

import pytest

from tests.unit.test_calibration_target_contract import study_document
from vntyper.scripts.calibration_length_metrics import EligibleLengthMember, LengthEligibleRoster
from vntyper.scripts.calibration_length_policy import (
    decode_length_source_truth,
    length_measurement_policy_sha256,
    length_protocol_for_roster,
)
from vntyper.scripts.calibration_target_contract import decode_target_study
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION

pytestmark = pytest.mark.unit


def _sample_study():
    raw_study = study_document()
    return decode_target_study(raw_study)


def test_length_measurement_policy_sha256_success() -> None:
    study = _sample_study()
    policy_hash = length_measurement_policy_sha256(study)
    expected_doc = {
        "schema_version": "length-measurement-policy-v1",
        "annotation_sha256": study.baseline.annotation_sha256,
        "counting_policy_sha256": study.baseline.counting_policy_sha256,
    }
    assert policy_hash == canonical_sha256(expected_doc)


def test_length_measurement_policy_sha256_requires_length_baseline() -> None:
    caller_study = decode_target_study(study_document("callers"))
    with pytest.raises(ValueError, match="length extraction requires a length study"):
        length_measurement_policy_sha256(caller_study)


def test_decode_length_source_truth_success() -> None:
    payload = {
        "schema_version": "calibration-length-truth-v1",
        "boundary_definition": TARGET_BOUNDARY_DEFINITION,
        "rows": [
            {"key": "specimen-1", "total_repeat_count": 42},
            {"key": "specimen-2", "total_repeat_count": 55.0},
        ],
    }
    keys = ("specimen-1", "specimen-2")
    truth = decode_length_source_truth(payload, keys)
    assert truth == {"specimen-1": 42.0, "specimen-2": 55.0}


def test_decode_length_source_truth_invalid_root() -> None:
    with pytest.raises(ValueError, match="length source truth fields differ"):
        decode_length_source_truth("not-a-dict", ())

    with pytest.raises(ValueError, match="length source truth fields differ"):
        decode_length_source_truth({"schema_version": "v1"}, ())


def test_decode_length_source_truth_invalid_schema_or_boundary() -> None:
    bad_schema = {
        "schema_version": "wrong-schema",
        "boundary_definition": TARGET_BOUNDARY_DEFINITION,
        "rows": [],
    }
    with pytest.raises(ValueError, match="exact repeat-count boundary"):
        decode_length_source_truth(bad_schema, ())

    bad_boundary = {
        "schema_version": "calibration-length-truth-v1",
        "boundary_definition": "wrong-boundary",
        "rows": [],
    }
    with pytest.raises(ValueError, match="exact repeat-count boundary"):
        decode_length_source_truth(bad_boundary, ())


def test_decode_length_source_truth_rows_must_be_list() -> None:
    doc = {
        "schema_version": "calibration-length-truth-v1",
        "boundary_definition": TARGET_BOUNDARY_DEFINITION,
        "rows": "not-a-list",
    }
    with pytest.raises(ValueError, match="rows must be a list"):
        decode_length_source_truth(doc, ())


def test_decode_length_source_truth_row_validation() -> None:
    base = {
        "schema_version": "calibration-length-truth-v1",
        "boundary_definition": TARGET_BOUNDARY_DEFINITION,
    }

    # Row is not a mapping
    with pytest.raises(ValueError, match="row fields differ"):
        decode_length_source_truth({**base, "rows": ["not-a-map"]}, ())

    # Row fields wrong
    with pytest.raises(ValueError, match="row fields differ"):
        decode_length_source_truth({**base, "rows": [{"key": "k"}]}, ())

    # Key is not a string
    with pytest.raises(ValueError, match="invalid or duplicate key"):
        decode_length_source_truth({**base, "rows": [{"key": 123, "total_repeat_count": 10}]}, ())

    # Duplicate key
    with pytest.raises(ValueError, match="invalid or duplicate key"):
        decode_length_source_truth(
            {**base, "rows": [{"key": "k", "total_repeat_count": 10}, {"key": "k", "total_repeat_count": 20}]},
            ("k",),
        )

    # Boolean count
    with pytest.raises(ValueError, match="positive integral counts"):
        decode_length_source_truth({**base, "rows": [{"key": "k", "total_repeat_count": True}]}, ("k",))

    # Non-number count
    with pytest.raises(ValueError, match="positive integral counts"):
        decode_length_source_truth({**base, "rows": [{"key": "k", "total_repeat_count": "42"}]}, ("k",))

    # Zero or negative count
    with pytest.raises(ValueError, match="positive integral counts"):
        decode_length_source_truth({**base, "rows": [{"key": "k", "total_repeat_count": 0}]}, ("k",))
    with pytest.raises(ValueError, match="positive integral counts"):
        decode_length_source_truth({**base, "rows": [{"key": "k", "total_repeat_count": -5}]}, ("k",))

    # Non-integer float
    with pytest.raises(ValueError, match="positive integral counts"):
        decode_length_source_truth({**base, "rows": [{"key": "k", "total_repeat_count": 42.5}]}, ("k",))

    # NaN / inf
    with pytest.raises(ValueError, match="positive integral counts"):
        decode_length_source_truth({**base, "rows": [{"key": "k", "total_repeat_count": float("nan")}]}, ("k",))

    # Overflow int to float
    huge_int = int("9" * 400)
    with pytest.raises(ValueError, match="count exceeds numeric range"):
        decode_length_source_truth({**base, "rows": [{"key": "k", "total_repeat_count": huge_int}]}, ("k",))

    # Mismatched keys
    with pytest.raises(ValueError, match="does not match the exact eligible artifact roster"):
        decode_length_source_truth({**base, "rows": [{"key": "k1", "total_repeat_count": 10}]}, ("k2",))


def test_length_protocol_for_roster_success() -> None:
    study = _sample_study()
    protocol = study.protocol
    member = EligibleLengthMember(
        key="specimen-1",
        group_key="group-1",
        strata=("stratum-a",),
    )
    roster = LengthEligibleRoster(
        members=(member,),
        sha256=canonical_sha256({"keys": ["specimen-1"]}),
    )
    derived = length_protocol_for_roster(protocol, roster)
    assert derived.eligible_roster_sha256 == roster.sha256
    assert derived.candidates == protocol.candidates
