"""Pure applicability, QC, and affine total-length inference decisions."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Literal, cast

from vntyper.scripts.length_features import LengthFeatures, RegionFeatures, encode_length_features
from vntyper.scripts.length_model import LengthModel, encode_length_model

logger = logging.getLogger(__name__)

EvidenceDomain = Literal["synthetic", "external"]
LengthEstimationStatus = Literal["disabled", "measured-only", "estimated", "unavailable"]


@dataclass(frozen=True)
class LengthEstimate:
    """One unrounded total-diploid prediction or an explicit availability state."""

    status: LengthEstimationStatus
    estimated_total_repeat_count: float | None
    reasons: tuple[str, ...]
    model_sha256: str | None
    features_sha256: str | None


def _evidence_domain(value: object) -> EvidenceDomain | None:
    if value is None:
        return None
    if not isinstance(value, str) or value not in {"synthetic", "external"}:
        raise ValueError("evidence_domain must be synthetic, external, or null")
    return cast(EvidenceDomain, value)


def _result(
    status: LengthEstimationStatus,
    features: LengthFeatures | None,
    model: LengthModel | None,
    *,
    prediction: float | None = None,
    reasons: tuple[str, ...] = (),
) -> LengthEstimate:
    return LengthEstimate(
        status=status,
        estimated_total_repeat_count=prediction,
        reasons=reasons,
        model_sha256=None if model is None else model.sha256,
        features_sha256=None if features is None else features.sha256,
    )


def _validated_features(value: object) -> LengthFeatures:
    if not isinstance(value, LengthFeatures):
        raise ValueError("length inference features must be LengthFeatures")
    encode_length_features(value)
    return value


def _validated_model(value: object) -> LengthModel:
    if not isinstance(value, LengthModel):
        raise ValueError("length inference model must be LengthModel")
    encode_length_model(value)
    return value


def _applicability_reason(
    features: LengthFeatures, model: LengthModel, evidence_domain: EvidenceDomain | None
) -> str | None:
    if evidence_domain is None:
        return "evidence_domain_missing"
    applicable = model.applicability
    context = features.provenance.measurement_context
    comparisons = (
        (evidence_domain == applicable.domain, "unsupported_evidence_domain"),
        (features.annotation_sha256 == model.annotation_sha256, "unsupported_annotation"),
        (features.counting_policy_sha256 == model.counting_policy_sha256, "unsupported_counting_policy"),
        (features.assembly in applicable.assemblies, "unsupported_assembly"),
        (features.assay_class in applicable.assay_classes, "unsupported_assay_class"),
        (features.input_scope in applicable.input_scopes, "unsupported_input_scope"),
        (context.preprocessing_id in applicable.preprocessing_ids, "unsupported_preprocessing"),
        (context.aligner.name == applicable.aligner_name, "unsupported_aligner_name"),
        (context.aligner.version == applicable.aligner_version, "unsupported_aligner_version"),
        (context.aligner.arguments_sha256 == applicable.aligner_arguments_sha256, "unsupported_aligner_arguments"),
        (
            context.aligner.primary_secondary_marking == applicable.primary_secondary_marking,
            "unsupported_primary_secondary_marking",
        ),
    )
    return next((reason for matches, reason in comparisons if not matches), None)


def _combined_flanks(features: LengthFeatures) -> RegionFeatures | None:
    left = features.regions["LEFT_FLANK"]
    right = features.regions["RIGHT_FLANK"]
    if left is None or right is None:
        return None
    length = left.length_bp + right.length_bp
    depth_sum = left.depth_sum + right.depth_sum
    support = features.provenance.denominator_qc.supporting_fragment_counts["BOTH_FLANKS"]
    return RegionFeatures(
        length_bp=length,
        depth_sum=depth_sum,
        mean_depth=depth_sum / length,
        covered_fraction=(left.covered_fraction * left.length_bp + right.covered_fraction * right.length_bp) / length,
        supporting_fragment_count=support,
    )


def _required_denominators(features: LengthFeatures, feature: str) -> tuple[tuple[str, RegionFeatures | None], ...]:
    if feature == "A":
        return (("invariant", features.regions["INVARIANT"]),)
    return (("combined_flanks", _combined_flanks(features)),)


def _qc_reasons(features: LengthFeatures, model: LengthModel, feature: str) -> tuple[str, ...]:
    evidence = features.provenance.denominator_qc
    if evidence.evidence_kind == "unavailable":
        return ("fragment_evidence_unavailable",)
    if evidence.evidence_kind != model.qc["fragment_evidence_kind"]:
        return ("unsupported_fragment_evidence_kind",)
    minimum_mean = cast(float, model.qc["minimum_denominator_mean_depth"])
    minimum_covered = cast(float, model.qc["minimum_denominator_covered_fraction"])
    minimum_support = cast(int, model.qc["minimum_denominator_supporting_fragments"])
    reasons: list[str] = []
    for name, region in _required_denominators(features, feature):
        if region is None:
            continue
        if region.mean_depth < minimum_mean:
            reasons.append(f"low_{name}_mean_depth")
        if region.covered_fraction < minimum_covered:
            reasons.append(f"low_{name}_covered_fraction")
        if region.supporting_fragment_count is None:
            reasons.append(f"missing_{name}_support")
        elif region.supporting_fragment_count < minimum_support:
            reasons.append(f"low_{name}_support")
    return tuple(reasons)


def estimate_total_repeats(
    features: LengthFeatures,
    model: LengthModel,
    *,
    evidence_domain: EvidenceDomain | None = None,
) -> LengthEstimate:
    """Apply a validated affine model after exact applicability and denominator QC.

    Args:
        features: Immutable measured A/F feature artifact.
        model: Immutable affine research model.
        evidence_domain: Explicit provenance domain supplied by a bound upstream context.

    Returns:
        One unrounded total-diploid estimate or stable sample-unavailability reasons.

    Raises:
        ValueError: If an explicit typed artifact or evidence domain is invalid.
    """
    checked_model = _validated_model(model)
    checked_features = _validated_features(features)
    domain = _evidence_domain(evidence_domain)
    applicability_reason = _applicability_reason(checked_features, checked_model, domain)
    if applicability_reason is not None:
        return _result("unavailable", checked_features, checked_model, reasons=(applicability_reason,))

    feature_name = checked_model.feature_order[0]
    feature_value = checked_features.a if feature_name == "A" else checked_features.f
    if feature_value is None:
        return _result("unavailable", checked_features, checked_model, reasons=(f"missing_{feature_name}",))
    qc_reasons = _qc_reasons(checked_features, checked_model, feature_name)
    if qc_reasons:
        return _result("unavailable", checked_features, checked_model, reasons=qc_reasons)
    bounds = checked_model.feature_bounds[feature_name]
    if feature_value < bounds.minimum or feature_value > bounds.maximum:
        return _result(
            "unavailable", checked_features, checked_model, reasons=(f"feature_{feature_name}_out_of_bounds",)
        )

    prediction = checked_model.intercept + checked_model.coefficients[0] * feature_value
    if not math.isfinite(prediction):
        return _result("unavailable", checked_features, checked_model, reasons=("nonfinite_prediction",))
    if prediction <= 0:
        return _result("unavailable", checked_features, checked_model, reasons=("nonpositive_prediction",))
    return _result("estimated", checked_features, checked_model, prediction=prediction)


def length_estimation_state(
    features: LengthFeatures | None,
    model: LengthModel | None,
    *,
    measurement_enabled: bool,
    evidence_domain: EvidenceDomain | None = None,
) -> LengthEstimate:
    """Dispatch disabled, measured-only, unavailable, and estimated states.

    Args:
        features: Optional immutable measurement artifact.
        model: Optional immutable affine research model.
        measurement_enabled: Strict boolean enabling length measurement.
        evidence_domain: Explicit provenance domain supplied by a bound upstream context.

    Returns:
        The explicit length estimation state and any available artifact identities.

    Raises:
        ValueError: If a supplied artifact, flag, or evidence domain is invalid.
    """
    if not isinstance(measurement_enabled, bool):
        raise ValueError("measurement_enabled must be boolean")
    domain = _evidence_domain(evidence_domain)
    checked_features = None if features is None else _validated_features(features)
    checked_model = None if model is None else _validated_model(model)
    if not measurement_enabled:
        if checked_model is not None:
            raise ValueError("length model requires measurement_enabled")
        return _result("disabled", checked_features, checked_model)
    if checked_features is None:
        return _result("unavailable", None, checked_model, reasons=("length_features_missing",))
    if checked_model is None:
        return _result("measured-only", checked_features, None, reasons=checked_features.reasons)
    return estimate_total_repeats(checked_features, checked_model, evidence_domain=domain)
