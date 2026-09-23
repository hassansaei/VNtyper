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


WORDS = {
    "length_sensitivity": {
        "caution": {"badge": "Above {threshold} repeats: caution", "help": "Caution help {threshold}."},
        "high": {
            "badge": "Above {threshold} repeats: high",
            "help": "High help {threshold}.",
            "notice_not_positive": "Estimate {estimate} ± {uncertainty} exceeds {threshold}.",
        },
    }
}


def _summary(estimate, tier):
    return {
        "length_estimation_status": "estimated",
        "estimated_total_repeat_count": estimate,
        "length_sensitivity_tier": tier,
        "length_sensitivity_policy": dict(RAW),
    }


def test_view_high_not_positive_has_notice() -> None:
    view = subject.build_sensitivity_view(_summary(160.456, "high"), WORDS, is_positive=False)
    assert view is not None
    assert view.tier == "high"
    assert view.badge == "Above 150 repeats: high"
    assert view.notice == "Estimate 160.46 ± 14 exceeds 150."
    assert view.help == "High help 150."
    assert view.uncertainty == "± 14"


def test_view_high_positive_has_no_notice() -> None:
    view = subject.build_sensitivity_view(_summary(160.0, "high"), WORDS, is_positive=True)
    assert view is not None and view.notice is None and view.badge is not None


def test_view_caution_badge_only() -> None:
    view = subject.build_sensitivity_view(_summary(120.0, "caution"), WORDS, is_positive=False)
    assert view is not None
    assert view.badge == "Above 110 repeats: caution"
    assert view.notice is None


def test_view_below_has_uncertainty_only() -> None:
    view = subject.build_sensitivity_view(_summary(90.0, "below"), WORDS, is_positive=False)
    assert view is not None
    assert (view.badge, view.notice, view.help, view.uncertainty) == (None, None, None, "± 14")


def test_view_not_assessed_is_silent() -> None:
    summary = {**_summary(None, "not-assessed"), "length_estimation_status": "unavailable"}
    view = subject.build_sensitivity_view(summary, WORDS, is_positive=False)
    assert view is not None
    assert (view.badge, view.notice, view.help, view.uncertainty) == (None, None, None, None)


def test_view_absent_without_recorded_tier_or_wording() -> None:
    assert subject.build_sensitivity_view({"length_estimation_status": "estimated"}, WORDS, is_positive=False) is None
    assert subject.build_sensitivity_view(_summary(160.0, "high"), {}, is_positive=False) is None


def test_view_rejects_tier_inconsistent_with_estimate() -> None:
    with pytest.raises(ValueError, match="tier differs"):
        subject.build_sensitivity_view(_summary(120.0, "high"), WORDS, is_positive=False)


def test_view_rejects_unknown_tier() -> None:
    with pytest.raises(ValueError, match="tier differs"):
        subject.build_sensitivity_view(_summary(120.0, "severe"), WORDS, is_positive=False)


@pytest.mark.parametrize(
    "words",
    [
        {"length_sensitivity": {"caution": WORDS["length_sensitivity"]["caution"]}},
        {"length_sensitivity": {**WORDS["length_sensitivity"], "extra": {}}},
        {"length_sensitivity": {**WORDS["length_sensitivity"], "caution": {"badge": "x"}}},
        {
            "length_sensitivity": {
                **WORDS["length_sensitivity"],
                "high": {"badge": "", "help": "h", "notice_not_positive": "n"},
            }
        },
    ],
)
def test_view_rejects_malformed_wording(words) -> None:
    with pytest.raises(ValueError, match="length sensitivity"):
        subject.build_sensitivity_view(_summary(160.0, "high"), words, is_positive=False)
