"""Research source and count convention remain explicit in report wording."""

from typing import Any

import pytest

from tests.unit.test_length_presentation import CONFIG
from vntyper.scripts.length_presentation import build_length_presentation

pytestmark = pytest.mark.unit

WORDS: dict[str, Any] = {
    "help": "Development model; variant calls are unchanged.",
    "extrapolation_help": "Some measurements lie outside the training range.",
    "sources": {"packaged-research": "Packaged research model.", "local-research": "Local research model."},
    "count_conventions": {
        "source-reported": "Assay-reported repeat counts; terminal-unit inclusion is unspecified.",
        "complete": "Complete repeat counts including invariant units.",
        "canonical-only": "Canonical repeat counts excluding invariant units.",
    },
}


def summary() -> dict[str, object]:
    return {
        "length_model_source": "packaged-research",
        "length_model_evidence_status": "research-development",
        "length_count_convention": "source-reported",
        "length_estimation_status": "unavailable",
        "estimated_total_repeat_count": None,
        "length_estimation_reasons": ["unsupported-assembly"],
        "length_model_sha256": "a" * 64,
        "length_calibration_id": None,
        "length_portable_approval_sha256": None,
        "length_standard_features": None,
        "length_standard_features_sha256": None,
    }


def test_unavailable_standard_model_has_no_fabricated_calibration_id() -> None:
    result = build_length_presentation(summary(), {**CONFIG, "standard_length_estimation": WORDS})
    assert result is not None
    assert result.status == "unavailable"
    assert result.calibration_id is None
    assert result.feature_values == ()
    assert result.help == " ".join(
        [WORDS["help"], WORDS["sources"]["packaged-research"], WORDS["count_conventions"]["source-reported"]]
    )


@pytest.mark.parametrize(
    "changes",
    [
        {"length_model_source": "approved"},
        {"length_calibration_id": "pretend-approved"},
        {"length_portable_approval_sha256": "b" * 64},
        {"length_count_convention": "unknown-convention"},
        {"estimated_total_repeat_count": 100},
        {"length_estimation_reasons": []},
        {"length_model_sha256": "invalid"},
        {"length_model_evidence_status": "validated"},
    ],
)
def test_invalid_research_summary_fails_closed(changes: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        build_length_presentation({**summary(), **changes}, {**CONFIG, "standard_length_estimation": WORDS})


def test_old_report_configuration_can_omit_standard_wording() -> None:
    assert build_length_presentation(summary(), CONFIG) is None


def test_measured_prediction_and_extrapolation_are_reported_with_configured_wording() -> None:
    from tests.unit.test_length_standard_model import _measurement
    from vntyper.scripts.length_standard_features import encode_standard_length_measurement

    measurement = _measurement()
    value = {
        **summary(),
        "length_estimation_status": "estimated",
        "estimated_total_repeat_count": 110.25,
        "length_estimation_reasons": [],
        "length_estimation_warnings": ["feature-range:A"],
        "length_standard_features": encode_standard_length_measurement(measurement),
        "length_standard_features_sha256": measurement.sha256,
    }
    result = build_length_presentation(value, {**CONFIG, "standard_length_estimation": WORDS})
    assert result is not None
    assert result.value == "110.25"
    assert result.feature_values[0][1] == "2"
    assert result.help.endswith(WORDS["extrapolation_help"])
    value["length_standard_features_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="digest differs"):
        build_length_presentation(value, {**CONFIG, "standard_length_estimation": WORDS})


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf"), True])
def test_invalid_estimated_count_is_rejected(value: object) -> None:
    with pytest.raises(ValueError):
        build_length_presentation(
            {
                **summary(),
                "length_estimation_status": "estimated",
                "estimated_total_repeat_count": value,
                "length_estimation_reasons": [],
            },
            {**CONFIG, "standard_length_estimation": WORDS},
        )


SENSITIVITY_WORDS: dict[str, Any] = {
    "caution": {"badge": "Above {threshold}: caution", "help": "Caution help."},
    "high": {"badge": "Above {threshold}: high", "help": "High help.", "notice_not_positive": "Notice {estimate}."},
}
POLICY = {"caution_threshold": 110.0, "high_threshold": 150.0, "uncertainty_repeats": 14.3}
FULL_CONFIG = {**CONFIG, "standard_length_estimation": WORDS, "length_sensitivity": SENSITIVITY_WORDS}


def estimated_summary(estimate: float, tier: str | None, warnings: list[str]) -> dict[str, object]:
    from tests.unit.test_length_standard_model import _measurement
    from vntyper.scripts.length_standard_features import encode_standard_length_measurement

    measurement = _measurement()
    value: dict[str, object] = {
        **summary(),
        "length_estimation_status": "estimated",
        "estimated_total_repeat_count": estimate,
        "length_count_convention": "complete",
        "length_estimation_reasons": [],
        "length_estimation_warnings": warnings,
        "length_standard_features": encode_standard_length_measurement(measurement),
        "length_standard_features_sha256": measurement.sha256,
    }
    if tier is not None:
        value["length_sensitivity_tier"] = tier
        value["length_sensitivity_policy"] = dict(POLICY)
    return value


def test_sensitivity_code_alone_does_not_claim_extrapolation() -> None:
    from vntyper.scripts.length_sensitivity import CAUTION_CODE

    result = build_length_presentation(estimated_summary(120.0, "caution", [CAUTION_CODE]), FULL_CONFIG)
    assert result is not None
    assert WORDS["extrapolation_help"] not in result.help
    assert result.help.endswith("Caution help.")


def test_high_tier_without_finding_carries_notice_badge_and_uncertainty() -> None:
    from vntyper.scripts.length_sensitivity import CAUTION_CODE, HIGH_CODE

    value = estimated_summary(160.5, "high", [CAUTION_CODE, HIGH_CODE])
    result = build_length_presentation(value, FULL_CONFIG, is_positive=False)
    assert result is not None
    assert result.sensitivity_tier == "high"
    assert result.warning_badge == "Above 150: high"
    assert result.notice_text == "Notice 160.5 ± 14."
    assert result.uncertainty_text == "± 14"
    positive = build_length_presentation(value, FULL_CONFIG, is_positive=True)
    assert positive is not None and positive.notice_text is None and positive.warning_badge == "Above 150: high"


def test_approved_path_receives_the_same_projection() -> None:
    from tests.unit.test_length_presentation import _summary as approved_summary

    value = {
        **approved_summary("estimated", 120.0),
        "length_sensitivity_tier": "caution",
        "length_sensitivity_policy": dict(POLICY),
    }
    result = build_length_presentation(value, {**CONFIG, "length_sensitivity": SENSITIVITY_WORDS}, is_positive=False)
    assert result is not None
    assert result.warning_badge == "Above 110: caution"
    assert result.notice_text is None


def test_summary_without_tier_has_no_sensitivity_fields() -> None:
    from vntyper.scripts.length_sensitivity import CAUTION_CODE

    value = {**estimated_summary(120.0, None, [CAUTION_CODE]), "length_warning_threshold": 110.0}
    result = build_length_presentation(value, FULL_CONFIG, is_positive=False)
    assert result is not None
    assert (result.sensitivity_tier, result.warning_badge, result.notice_text, result.uncertainty_text) == (
        None,
        None,
        None,
        None,
    )
