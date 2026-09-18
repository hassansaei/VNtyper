"""Pure presentation decisions for optional VNTR length results."""

from __future__ import annotations

from dataclasses import replace

import pytest

pytestmark = pytest.mark.unit


CONFIG = {
    "length_estimation": {
        "label": "Approximate total diploid VNTR repeat count",
        "unit": "repeat units",
        "help": "Research-only estimate across both alleles; it does not affect variant calls.",
        "status_labels": {
            "not-recorded": "Not recorded by this run",
            "disabled": "Not requested",
            "measured-only": "Features measured; no model applied",
            "estimated": "Estimated",
            "unavailable": "Unavailable",
        },
        "feature_labels": {
            "A": "Core/invariant depth ratio (A)",
            "F": "Array/flank depth ratio (F)",
        },
    }
}


def _summary(status: str, estimate: float | None) -> dict[str, object]:
    from tests.unit.test_length_estimation import features as make_features
    from vntyper.scripts.canonical_json import canonical_sha256
    from vntyper.scripts.length_features import encode_length_features

    features = encode_length_features(make_features(manifest_key="private-member-must-not-render"))
    has_model = status in {"estimated", "unavailable"}

    return {
        "length_estimation_status": status,
        "estimated_total_repeat_count": estimate,
        "length_calibration_id": "affine-a" if has_model else None,
        "length_estimation_reasons": [] if status != "unavailable" else ["low_denominator_depth"],
        "length_features": features,
        "length_features_sha256": canonical_sha256(features),
        "length_model_sha256": "d" * 64 if has_model else None,
    }


def test_legacy_disabled_measured_unavailable_and_estimated_states() -> None:
    from vntyper.scripts.length_presentation import build_length_presentation

    legacy = build_length_presentation({}, CONFIG)
    assert legacy is not None
    assert legacy.status == "not-recorded"
    assert legacy.value == "Not recorded by this run"

    disabled = build_length_presentation(
        {
            "length_estimation_status": "disabled",
            "estimated_total_repeat_count": None,
            "length_calibration_id": None,
            "length_estimation_reasons": [],
            "length_features": None,
            "length_features_sha256": None,
            "length_model_sha256": None,
        },
        CONFIG,
    )
    assert disabled is not None
    assert disabled.status_text == "Not requested"
    assert disabled.feature_values == ()

    measured = build_length_presentation(_summary("measured-only", None), CONFIG)
    assert measured is not None
    assert measured.value == "Not calculated"
    assert measured.feature_values == (("Core/invariant depth ratio (A)", "2"), ("Array/flank depth ratio (F)", "1.5"))

    unavailable = build_length_presentation(_summary("unavailable", None), CONFIG)
    assert unavailable is not None
    assert unavailable.reasons == ("low_denominator_depth",)

    estimated = build_length_presentation(_summary("estimated", 110.25), CONFIG)
    assert estimated is not None
    assert estimated.value == "110.25"
    assert estimated.calibration_id == "affine-a"
    assert estimated.model_sha256 == "d" * 64


def test_missing_old_config_hides_section_and_present_config_is_closed() -> None:
    from vntyper.scripts.length_presentation import build_length_presentation

    assert build_length_presentation({}, {}) is None
    with pytest.raises(ValueError, match="configuration"):
        build_length_presentation({}, {"length_estimation": {**CONFIG["length_estimation"], "extra": True}})
    with pytest.raises(ValueError, match="without status"):
        build_length_presentation({"estimated_total_repeat_count": 110}, CONFIG)


@pytest.mark.parametrize(
    "change",
    [
        {"length_estimation_status": "estimated", "estimated_total_repeat_count": None},
        {"length_estimation_status": "disabled", "length_features": _summary("estimated", 1)["length_features"]},
        {"length_estimation_status": "estimated", "estimated_total_repeat_count": float("inf")},
        {"length_estimation_status": "estimated", "length_features_sha256": "0" * 64},
    ],
)
def test_inconsistent_or_tampered_new_summary_fails_closed(change: dict[str, object]) -> None:
    from vntyper.scripts.length_presentation import build_length_presentation

    summary = {**_summary("estimated", 110.0), **change}
    with pytest.raises(ValueError):
        build_length_presentation(summary, CONFIG)


def test_presentation_is_frozen() -> None:
    from vntyper.scripts.length_presentation import build_length_presentation

    presentation = build_length_presentation(_summary("estimated", 110.0), CONFIG)
    assert presentation is not None
    with pytest.raises(AttributeError):
        replace(presentation, value="forged").feature_values.append(("x", "y"))  # type: ignore[attr-defined]
