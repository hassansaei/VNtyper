"""Pure total-diploid length inference and availability decisions."""

from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from importlib import import_module

import pytest

from tests.unit.test_length_features import _annotation, _annotation_raw, _context
from tests.unit.test_length_model import model_document
from vntyper.scripts.length_annotation import decode_length_annotation
from vntyper.scripts.length_features import DepthPosition, extract_length_features

pytestmark = pytest.mark.unit


def features(*, missing_region=None, depths=(100, 200, 200, 100, 100, 100), fragment_evidence=True):
    annotation_raw = _annotation_raw()
    if missing_region is None:
        annotation = _annotation()
    else:
        annotation_raw["regions"][missing_region] = None
        annotation_raw["array_boundary_geometry"] = {"array_only_bp": None, "target_only_bp": None}
        annotation = decode_length_annotation(annotation_raw)
    fragment_ids = tuple(f"fragment-{index}" for index in range(100))
    positions = tuple(
        DepthPosition(
            contig=annotation.contig,
            position_zero_based=position,
            depth=depth,
            supporting_fragment_ids=(fragment_ids if depth else ()) if fragment_evidence else None,
        )
        for position, depth in enumerate(depths)
    )
    return extract_length_features(positions, annotation, _context(annotation))


def model_for(measured, feature="A", **changes):
    raw = model_document(feature)
    context = measured.provenance.measurement_context
    raw["annotation_sha256"] = measured.annotation_sha256
    raw["counting_policy_sha256"] = measured.counting_policy_sha256
    raw["applicability"].update(
        assemblies=[measured.assembly],
        assay_classes=[measured.assay_class],
        input_scopes=[measured.input_scope],
        preprocessing_ids=[context.preprocessing_id],
        aligner_name=context.aligner.name,
        aligner_version=context.aligner.version,
        aligner_arguments_sha256=context.aligner.arguments_sha256,
        primary_secondary_marking=context.aligner.primary_secondary_marking,
        counting_policy_sha256=measured.counting_policy_sha256,
    )
    for section, values in changes.items():
        if section in {"applicability", "qc", "target", "feature_bounds"}:
            target = raw[section][feature] if section == "feature_bounds" else raw[section]
            target.update(values)
        else:
            raw[section] = values
            if section == "counting_policy_sha256":
                raw["applicability"]["counting_policy_sha256"] = values
    return import_module("vntyper.scripts.length_model").decode_length_model(raw)


def test_affine_a_returns_one_unrounded_total_diploid_prediction():
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    model = model_for(measured)
    result = m.estimate_total_repeats(measured, model, evidence_domain="synthetic")
    assert measured.a == 2
    assert result.status == "estimated"
    assert result.estimated_total_repeat_count == 110
    assert result.reasons == ()
    assert result.model_sha256 == model.sha256
    assert result.features_sha256 == measured.sha256
    assert not hasattr(result, "allele_1")
    assert not hasattr(result, "allele_2")
    with pytest.raises(FrozenInstanceError):
        result.status = "unavailable"


def test_affine_f_uses_only_f_and_does_not_round():
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    result = m.estimate_total_repeats(measured, model_for(measured, "F"), evidence_domain="synthetic")
    assert measured.f == 1.5
    assert result.estimated_total_repeat_count == 85.0


@pytest.mark.parametrize("feature,bounds", [("A", {"minimum": 2, "maximum": 3}), ("A", {"minimum": 1, "maximum": 2})])
def test_feature_boundaries_are_inclusive(feature, bounds):
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    result = m.estimate_total_repeats(
        measured, model_for(measured, feature, feature_bounds=bounds), evidence_domain="synthetic"
    )
    assert result.status == "estimated"


@pytest.mark.parametrize("feature,missing_region,reason", [("A", "CORE", "missing_A"), ("F", "ARRAY", "missing_F")])
def test_missing_required_feature_is_sample_unavailability(feature, missing_region, reason):
    m = import_module("vntyper.scripts.length_estimation")
    measured = features(missing_region=missing_region)
    result = m.estimate_total_repeats(measured, model_for(measured, feature), evidence_domain="synthetic")
    assert result.status == "unavailable"
    assert result.estimated_total_repeat_count is None
    assert result.reasons == (reason,)


@pytest.mark.parametrize(
    "change,reason",
    [
        ({"assemblies": ["other"]}, "unsupported_assembly"),
        ({"assay_classes": ["other"]}, "unsupported_assay_class"),
        ({"input_scopes": ["full"]}, "unsupported_input_scope"),
        ({"preprocessing_ids": ["other"]}, "unsupported_preprocessing"),
        ({"aligner_name": "other"}, "unsupported_aligner_name"),
        ({"aligner_version": "other"}, "unsupported_aligner_version"),
        ({"aligner_arguments_sha256": "9" * 64}, "unsupported_aligner_arguments"),
        ({"primary_secondary_marking": "other"}, "unsupported_primary_secondary_marking"),
    ],
)
def test_observable_applicability_mismatch_is_unavailable(change, reason):
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    result = m.estimate_total_repeats(measured, model_for(measured, applicability=change), evidence_domain="synthetic")
    assert result.reasons == (reason,)


def test_annotation_counting_policy_and_domain_are_explicit_applicability():
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    cases = (
        (model_for(measured, annotation_sha256="9" * 64), "synthetic", "unsupported_annotation"),
        (model_for(measured, counting_policy_sha256="9" * 64), "synthetic", "unsupported_counting_policy"),
        (model_for(measured), None, "evidence_domain_missing"),
        (model_for(measured), "external", "unsupported_evidence_domain"),
    )
    for model, domain, reason in cases:
        assert m.estimate_total_repeats(measured, model, evidence_domain=domain).reasons == (reason,)
    for invalid_domain in ("guessed", True):
        with pytest.raises(ValueError, match="evidence_domain"):
            m.estimate_total_repeats(measured, model_for(measured), evidence_domain=invalid_domain)


@pytest.mark.parametrize(
    "feature,qc_change,depths,reason",
    [
        (
            "A",
            {"minimum_denominator_mean_depth": 101},
            (100, 200, 200, 100, 100, 100),
            "low_invariant_mean_depth",
        ),
        (
            "A",
            {"minimum_denominator_covered_fraction": 1},
            (100, 200, 200, 0, 100, 100),
            "low_invariant_covered_fraction",
        ),
        (
            "A",
            {"minimum_denominator_supporting_fragments": 101},
            (100, 200, 200, 100, 100, 100),
            "low_invariant_support",
        ),
        (
            "F",
            {"minimum_denominator_mean_depth": 101},
            (100, 200, 200, 100, 100, 100),
            "low_combined_flanks_mean_depth",
        ),
        (
            "F",
            {"minimum_denominator_covered_fraction": 1},
            (0, 200, 200, 100, 100, 100),
            "low_combined_flanks_covered_fraction",
        ),
        (
            "F",
            {"minimum_denominator_supporting_fragments": 101},
            (100, 200, 200, 100, 100, 100),
            "low_combined_flanks_support",
        ),
    ],
)
def test_each_per_denominator_qc_failure_is_sample_unavailability(feature, qc_change, depths, reason):
    m = import_module("vntyper.scripts.length_estimation")
    measured = features(depths=depths)
    result = m.estimate_total_repeats(measured, model_for(measured, feature, qc=qc_change), evidence_domain="synthetic")
    assert result.status == "unavailable"
    assert reason in result.reasons


@pytest.mark.parametrize("feature", ["A", "F"])
def test_qc_threshold_equality_passes_for_mean_coverage_and_support(feature):
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    model = model_for(
        measured,
        feature,
        qc={
            "minimum_denominator_mean_depth": 100,
            "minimum_denominator_covered_fraction": 1,
            "minimum_denominator_supporting_fragments": 100,
        },
    )
    assert m.estimate_total_repeats(measured, model, evidence_domain="synthetic").status == "estimated"


def test_depth_only_features_cannot_invent_required_fragment_evidence():
    m = import_module("vntyper.scripts.length_estimation")
    measured = features(fragment_evidence=False)
    result = m.estimate_total_repeats(measured, model_for(measured), evidence_domain="synthetic")
    assert result.status == "unavailable"
    assert result.reasons == ("fragment_evidence_unavailable",)


@pytest.mark.parametrize(
    "feature,depths",
    [
        ("A", (200, 100, 100, 200, 200, 200)),
        ("F", (200, 100, 100, 100, 100, 200)),
    ],
)
def test_low_nonzero_numerator_does_not_change_denominator_qc_availability(feature, depths):
    m = import_module("vntyper.scripts.length_estimation")
    measured = features(depths=depths)
    model = model_for(
        measured,
        feature,
        qc={"minimum_denominator_mean_depth": 150},
        feature_bounds={"minimum": 0.25, "maximum": 3},
    )
    assert m.estimate_total_repeats(measured, model, evidence_domain="synthetic").status == "estimated"


def test_out_of_bounds_nonpositive_and_nonfinite_outputs_are_unavailable_not_clipped():
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    outside = model_for(measured, feature_bounds={"minimum": 2.1, "maximum": 3})
    assert m.estimate_total_repeats(measured, outside, evidence_domain="synthetic").reasons == (
        "feature_A_out_of_bounds",
    )
    negative = model_for(measured, intercept=-200)
    result = m.estimate_total_repeats(measured, negative, evidence_domain="synthetic")
    assert result.estimated_total_repeat_count is None
    assert result.reasons == ("nonpositive_prediction",)
    overflow = model_for(measured, coefficients=[1e308])
    assert m.estimate_total_repeats(measured, overflow, evidence_domain="synthetic").reasons == (
        "nonfinite_prediction",
    )


def test_dispatcher_has_explicit_disabled_measured_only_and_unavailable_states():
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    disabled = m.length_estimation_state(None, None, measurement_enabled=False)
    assert disabled.status == "disabled"
    with pytest.raises(ValueError, match="requires measurement_enabled"):
        m.length_estimation_state(measured, model_for(measured), measurement_enabled=False)
    measured_only = m.length_estimation_state(measured, None, measurement_enabled=True)
    assert measured_only.status == "measured-only"
    unavailable = m.length_estimation_state(
        None, model_for(measured), measurement_enabled=True, evidence_domain="synthetic"
    )
    assert unavailable.status == "unavailable"
    assert unavailable.reasons == ("length_features_missing",)
    with pytest.raises(ValueError, match="measurement_enabled"):
        m.length_estimation_state(measured, None, measurement_enabled=1)


def test_explicit_invalid_model_or_features_are_errors():
    m = import_module("vntyper.scripts.length_estimation")
    measured = features()
    model = model_for(measured)
    with pytest.raises(ValueError, match="LengthModel"):
        m.estimate_total_repeats(measured, {}, evidence_domain="synthetic")
    forged = replace(measured, sha256="0" * 64)
    with pytest.raises(ValueError, match="digest"):
        m.estimate_total_repeats(forged, model, evidence_domain="synthetic")
    model_raw = deepcopy(model_document())
    model_raw["coefficients"] = ["__import__('os')"]
    with pytest.raises(ValueError):
        import_module("vntyper.scripts.length_model").decode_length_model(model_raw)
