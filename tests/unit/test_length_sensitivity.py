"""Two-tier length sensitivity policy: validation, tiers and summary projection."""

import json
import math
from pathlib import Path

import pytest

import vntyper
from vntyper.scripts import length_sensitivity as subject

pytestmark = pytest.mark.unit

POLICY = subject.LengthSensitivityPolicy(110.0, 150.0, 14.3)
RAW = {"caution_threshold": 110.0, "high_threshold": 150.0, "typical_error_repeats": 14.3}


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
        {k: v for k, v in RAW.items() if k != "typical_error_repeats"},
        {**RAW, "uncertainty_repeats": 14.3},
        {**RAW, "caution_threshold": True},
        {**RAW, "caution_threshold": "110"},
        {**RAW, "high_threshold": math.inf},
        {**RAW, "typical_error_repeats": math.nan},
        {**RAW, "caution_threshold": 0},
        {**RAW, "typical_error_repeats": -1},
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
    [(80.0, "below"), (110.0, "below"), (110.1, "caution"), (150.0, "caution"), (150.1, "high"), (175, "high")],
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


LABELS = {
    "notice_prefix": "Caution:",
    "cohort_kpi_label": "Length Tier High",
    "cohort_kpi_detail": "{caution} caution, {assessed} assessed",
    "typical_error": "(typical error ≈{error} repeats)",
}
WORDS = {
    "length_sensitivity": {
        "labels": LABELS,
        "caution": {"badge": "Above {threshold} repeats: caution", "help": "Caution help {threshold}."},
        "high": {
            "badge": "Above {threshold} repeats: high",
            "help": "High help {threshold}.",
            "notice_not_positive": "Estimate {estimate} exceeds {threshold}.",
        },
    }
}


_CODES = {"caution": [subject.CAUTION_CODE], "high": [subject.CAUTION_CODE, subject.HIGH_CODE]}


def _summary(estimate, tier):
    return {
        "length_estimation_status": "estimated",
        "estimated_total_repeat_count": estimate,
        "length_estimation_warnings": list(_CODES.get(tier, [])),
        "length_sensitivity_tier": tier,
        "length_sensitivity_policy": dict(RAW),
        "length_model_source": "packaged-research",
        "length_count_convention": "complete",
    }


def test_view_high_not_positive_has_notice() -> None:
    view = subject.build_sensitivity_view(_summary(160.456, "high"), WORDS, is_positive=False)
    assert view is not None
    assert view.tier == "high"
    assert view.badge == "Above 150 repeats: high"
    assert view.notice == "Estimate 160.5 exceeds 150."
    assert view.help == "High help 150."
    assert view.typical_error == "(typical error ≈14 repeats)"


def test_view_high_positive_has_no_notice() -> None:
    view = subject.build_sensitivity_view(_summary(160.0, "high"), WORDS, is_positive=True)
    assert view is not None and view.notice is None and view.badge is not None


def test_view_caution_badge_only() -> None:
    view = subject.build_sensitivity_view(_summary(120.0, "caution"), WORDS, is_positive=False)
    assert view is not None
    assert view.badge == "Above 110 repeats: caution"
    assert view.notice is None


def test_view_below_has_typical_error_only() -> None:
    view = subject.build_sensitivity_view(_summary(90.0, "below"), WORDS, is_positive=False)
    assert view is not None
    assert (view.badge, view.notice, view.help, view.typical_error) == (
        None,
        None,
        None,
        "(typical error ≈14 repeats)",
    )


def test_view_not_assessed_is_silent() -> None:
    summary = {**_summary(None, "not-assessed"), "length_estimation_status": "unavailable"}
    view = subject.build_sensitivity_view(summary, WORDS, is_positive=False)
    assert view is not None
    assert (view.badge, view.notice, view.help, view.typical_error) == (None, None, None, None)


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
        {"length_sensitivity": {**WORDS["length_sensitivity"], "labels": {"notice_prefix": "Caution:"}}},
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


@pytest.mark.parametrize(
    ("convention", "estimate", "tier"),
    [
        (None, 120.0, "caution"),
        ("complete", 120.0, "caution"),
        ("canonical-only", 95.0, "caution"),
        ("canonical-only", 92.0, "below"),
        ("canonical-only", 140.0, "high"),
        ("source-reported", 160.0, "not-assessed"),
        ("unknown", 160.0, "not-assessed"),
    ],
)
def test_tier_is_compared_in_the_complete_count_frame(convention, estimate, tier) -> None:
    fields = {"length_estimation_status": "estimated", "estimated_total_repeat_count": estimate}
    if convention is not None:
        fields["length_count_convention"] = convention
    assert subject.apply_length_sensitivity(fields, POLICY)["length_sensitivity_tier"] == tier


def test_canonical_only_view_rechecks_in_the_complete_frame() -> None:
    summary = {**_summary(95.0, "caution"), "length_count_convention": "canonical-only"}
    view = subject.build_sensitivity_view(summary, WORDS, is_positive=False)
    assert view is not None and view.tier == "caution"


@pytest.mark.parametrize("source", ["local-research", None])
def test_typical_error_is_shown_only_for_the_packaged_model(source) -> None:
    summary = _summary(160.456, "high")
    if source is None:
        del summary["length_model_source"]
        del summary["length_count_convention"]
    else:
        summary["length_model_source"] = source
    view = subject.build_sensitivity_view(summary, WORDS, is_positive=False)
    assert view is not None
    assert view.typical_error is None
    assert view.notice == "Estimate 160.5 exceeds 150."


@pytest.mark.parametrize(
    ("estimate", "tier"),
    [(110.01, "below"), (110.04, "below"), (110.06, "caution"), (150.04, "caution"), (150.06, "high")],
)
def test_tier_agrees_with_the_displayed_tenth(estimate, tier) -> None:
    assert subject.classify_length_sensitivity("estimated", estimate, POLICY) == tier


def test_notice_never_shows_the_threshold_as_exceeding_itself() -> None:
    view = subject.build_sensitivity_view(_summary(150.06, "high"), WORDS, is_positive=False)
    assert view is not None
    assert view.notice == "Estimate 150.1 exceeds 150."


def test_canonical_only_notice_shows_the_complete_frame_value() -> None:
    summary = {**_summary(140.0, "high"), "length_count_convention": "canonical-only"}
    view = subject.build_sensitivity_view(summary, WORDS, is_positive=False)
    assert view is not None and view.notice == "Estimate 158 exceeds 150."


def test_value_drops_to_one_decimal_beside_the_typical_error() -> None:
    view = subject.build_sensitivity_view(_summary(123.456, "caution"), WORDS, is_positive=False)
    assert view is not None and view.value == "123.5"
    bare = {**_summary(123.456, "caution"), "length_model_source": "local-research"}
    bare_view = subject.build_sensitivity_view(bare, WORDS, is_positive=False)
    assert bare_view is not None and bare_view.value is None


def test_view_carries_the_configured_notice_prefix() -> None:
    view = subject.build_sensitivity_view(_summary(160.0, "high"), WORDS, is_positive=False)
    assert view is not None and view.notice_prefix == "Caution:"


@pytest.mark.parametrize(
    ("tier", "estimate", "codes"),
    [
        ("high", 160.0, [subject.CAUTION_CODE]),
        ("caution", 120.0, []),
        ("caution", 120.0, [subject.CAUTION_CODE, subject.HIGH_CODE]),
        ("below", 90.0, [subject.CAUTION_CODE]),
    ],
)
def test_view_rejects_codes_that_differ_from_the_tier(tier, estimate, codes) -> None:
    summary = {**_summary(estimate, tier), "length_estimation_warnings": codes}
    with pytest.raises(ValueError, match="warning codes differ"):
        subject.build_sensitivity_view(summary, WORDS, is_positive=False)


def test_apply_rejects_non_list_warnings() -> None:
    fields = {
        "length_estimation_status": "estimated",
        "estimated_total_repeat_count": 120.0,
        "length_estimation_warnings": "feature_A_outside_training_range",
    }
    with pytest.raises(ValueError, match="warnings must be a list"):
        subject.apply_length_sensitivity(fields, POLICY)


def test_cohort_kpi_text_comes_from_configuration() -> None:
    assert subject.cohort_kpi_text(WORDS, {"high": 2, "caution": 5, "assessed": 9}) == (
        "Length Tier High",
        "5 caution, 9 assessed",
    )
    assert subject.cohort_kpi_text({}, {"high": 2, "caution": 5, "assessed": 9}) is None
    assert subject.cohort_kpi_text(WORDS, None) is None


def test_a_summary_recorded_under_the_old_field_name_still_decodes() -> None:
    """Summaries written before the rename carry the same RMSE as ``uncertainty_repeats``."""
    legacy = {"caution_threshold": 110.0, "high_threshold": 150.0, "uncertainty_repeats": 14.3}
    assert subject.decode_length_sensitivity_policy(legacy) == POLICY
    summary = {**_summary(160.456, "high"), "length_sensitivity_policy": legacy}
    view = subject.build_sensitivity_view(summary, WORDS, is_positive=False)
    assert view is not None and view.typical_error == "(typical error ≈14 repeats)"
    assert (
        subject.apply_length_sensitivity(
            {"length_estimation_status": "estimated", "estimated_total_repeat_count": 90.0}, POLICY
        )["length_sensitivity_policy"]
        == RAW
    )


def test_the_typical_error_is_never_rendered_as_a_plus_minus_interval() -> None:
    """The RMSE is not a 95% interval (55/76 residuals within 14), so no view text says "±"."""
    for estimate, tier in ((90.0, "below"), (120.0, "caution"), (160.456, "high")):
        view = subject.build_sensitivity_view(_summary(estimate, tier), WORDS, is_positive=False)
        assert view is not None
        rendered = (view.badge, view.notice, view.help, view.typical_error, view.value)
        assert all("±" not in text for text in rendered if text is not None)


def test_the_shipped_wording_names_the_rmse_as_a_typical_error() -> None:
    config = json.loads((Path(vntyper.__file__).parent / "scripts" / "report_config.json").read_text())
    view = subject.build_sensitivity_view(_summary(90.0, "below"), config, is_positive=False)
    assert view is not None
    assert view.typical_error == "(typical error ≈14 repeats, leave-one-out RMSE; not an interval)"
