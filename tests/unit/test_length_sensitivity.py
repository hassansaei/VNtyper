"""Two-tier length sensitivity policy: validation, tiers and summary projection."""

import json
import math
from pathlib import Path

import pytest

import vntyper
from vntyper.scripts import length_sensitivity as subject

pytestmark = pytest.mark.unit

POLICY = subject.LengthSensitivityPolicy(110.0, 150.0, 14.3)
RAW = {"caution_threshold": 110.0, "high_threshold": 150.0, "uncertainty_repeats": 14.3}


def test_resolve_reads_the_shipped_default() -> None:
    config = json.loads((Path(vntyper.__file__).parent / "config.json").read_text())
    assert subject.resolve_length_sensitivity_policy(config) == POLICY


def test_resolve_without_block_is_none() -> None:
    assert subject.resolve_length_sensitivity_policy({}) is None
    assert subject.resolve_length_sensitivity_policy({"length_estimation": {"enabled": False}}) is None


@pytest.mark.parametrize(
    "raw",
    [
        None,
        [],
        {**RAW, "extra": 1},
        {k: v for k, v in RAW.items() if k != "uncertainty_repeats"},
        {**RAW, "caution_threshold": True},
        {**RAW, "caution_threshold": "110"},
        {**RAW, "high_threshold": math.inf},
        {**RAW, "uncertainty_repeats": math.nan},
        {**RAW, "caution_threshold": 0},
        {**RAW, "uncertainty_repeats": -1},
        {**RAW, "caution_threshold": 150.0},
        {**RAW, "caution_threshold": 160.0},
    ],
)
def test_decode_rejects_malformed_policy(raw) -> None:
    with pytest.raises(ValueError, match="length sensitivity"):
        subject.decode_length_sensitivity_policy(raw)


def test_resolve_rejects_non_mapping_section() -> None:
    with pytest.raises(ValueError, match="length_estimation"):
        subject.resolve_length_sensitivity_policy({"length_estimation": []})


@pytest.mark.parametrize(
    ("estimate", "tier"),
    [(80.0, "below"), (110.0, "below"), (110.01, "caution"), (150.0, "caution"), (150.01, "high"), (175, "high")],
)
def test_tier_boundaries_are_strict(estimate, tier) -> None:
    assert subject.classify_length_sensitivity("estimated", estimate, POLICY) == tier


@pytest.mark.parametrize(
    ("status", "estimate"),
    [
        ("unavailable", None),
        ("disabled", None),
        ("measured-only", None),
        ("estimated", None),
        ("estimated", math.nan),
        ("estimated", True),
        ("estimated", "120"),
    ],
)
def test_not_estimated_is_not_assessed(status, estimate) -> None:
    assert subject.classify_length_sensitivity(status, estimate, POLICY) == "not-assessed"


def test_apply_records_tier_policy_and_codes_without_mutating() -> None:
    fields = {
        "length_estimation_status": "estimated",
        "estimated_total_repeat_count": 160.0,
        "length_estimation_warnings": ["feature_A_outside_training_range"],
    }
    result = subject.apply_length_sensitivity(fields, POLICY)
    assert result["length_sensitivity_tier"] == "high"
    assert result["length_sensitivity_policy"] == RAW
    assert result["length_estimation_warnings"] == [
        "feature_A_outside_training_range",
        subject.CAUTION_CODE,
        subject.HIGH_CODE,
    ]
    assert fields["length_estimation_warnings"] == ["feature_A_outside_training_range"]


def test_apply_creates_warnings_on_the_approved_path() -> None:
    result = subject.apply_length_sensitivity(
        {"length_estimation_status": "estimated", "estimated_total_repeat_count": 120.0}, POLICY
    )
    assert result["length_estimation_warnings"] == [subject.CAUTION_CODE]
    assert result["length_sensitivity_tier"] == "caution"


def test_apply_below_and_not_assessed_add_no_codes() -> None:
    below = subject.apply_length_sensitivity(
        {"length_estimation_status": "estimated", "estimated_total_repeat_count": 90.0}, POLICY
    )
    assert below["length_estimation_warnings"] == []
    assert below["length_sensitivity_tier"] == "below"
    unavailable = subject.apply_length_sensitivity({"length_estimation_status": "unavailable"}, POLICY)
    assert unavailable["length_sensitivity_tier"] == "not-assessed"


def test_apply_without_policy_or_disabled_is_identity() -> None:
    fields = {"length_estimation_status": "estimated", "estimated_total_repeat_count": 160.0}
    assert subject.apply_length_sensitivity(fields, None) == fields
    disabled = {"length_estimation_status": "disabled", "estimated_total_repeat_count": None}
    assert subject.apply_length_sensitivity(disabled, POLICY) == disabled


def test_apply_does_not_duplicate_existing_codes() -> None:
    fields = {
        "length_estimation_status": "estimated",
        "estimated_total_repeat_count": 120.0,
        "length_estimation_warnings": [subject.CAUTION_CODE],
    }
    assert subject.apply_length_sensitivity(fields, POLICY)["length_estimation_warnings"] == [subject.CAUTION_CODE]
