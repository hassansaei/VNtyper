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
