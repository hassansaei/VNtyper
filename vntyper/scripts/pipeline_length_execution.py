"""Execute optional length measurement against a retained alignment plan."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.length_depth_io import read_length_depth
from vntyper.scripts.length_estimation import LengthEstimate, length_estimation_state
from vntyper.scripts.length_features import LengthFeatures, encode_length_features, extract_length_features
from vntyper.scripts.pipeline_length import (
    LengthPipelineConfiguration,
    encode_length_pipeline_configuration,
)


@dataclass(frozen=True)
class LengthPipelineMeasurement:
    """Measured features and their optional approved-model estimate."""

    features: LengthFeatures
    estimate: LengthEstimate
    configuration_sha256: str


def measure_pipeline_length(
    *,
    plan: AlignmentPlan,
    reference_path: Path | None,
    samtools_path: Path,
    configuration: LengthPipelineConfiguration,
) -> LengthPipelineMeasurement:
    """Measure and estimate length before the alignment plan is released.

    Args:
        plan: Retained preflight plan whose descriptor-backed view is still open.
        reference_path: Exact local FASTA consumed by depth and fragment readers.
        samtools_path: Exact samtools executable consumed by the frozen policy.
        configuration: Enabled, preflighted length configuration.

    Returns:
        Immutable feature and estimation result bound to the configuration.

    Raises:
        ValueError: If typed inputs, configuration, or reference selection differ.
        RuntimeError: If depth or alignment reading fails.
    """
    if not isinstance(plan, AlignmentPlan):
        raise ValueError("length measurement requires a retained AlignmentPlan")
    encode_length_pipeline_configuration(configuration)
    if (
        not configuration.measurement_enabled
        or configuration.annotation is None
        or configuration.measurement_context is None
    ):
        raise ValueError("length measurement requires an enabled configuration")
    if not isinstance(reference_path, Path) or not isinstance(samtools_path, Path):
        raise ValueError("length measurement requires explicit reference and samtools Paths")
    if plan.file_format == "cram" and plan.reference_path != str(reference_path):
        raise ValueError("length measurement CRAM reference differs from the proven alignment plan")
    depths = read_length_depth(
        Path(plan.view_path),
        reference_path,
        configuration.annotation,
        configuration.measurement_context,
        samtools_path,
    )
    features = extract_length_features(depths, configuration.annotation, configuration.measurement_context)
    estimate = length_estimation_state(
        features,
        configuration.model,
        measurement_enabled=True,
        evidence_domain=configuration.evidence_domain,
    )
    return LengthPipelineMeasurement(features, estimate, configuration.sha256)


def length_summary_fields(
    configuration: LengthPipelineConfiguration,
    measurement: LengthPipelineMeasurement | None,
) -> dict[str, object]:
    """Project stable optional length fields for the pipeline summary.

    Args:
        configuration: Resolved length configuration for this run.
        measurement: Completed result, required exactly when measurement is enabled.

    Returns:
        Summary fields including model, annotation, policy and feature identities.

    Raises:
        ValueError: If result presence, typed content, or configuration binding differs.
    """
    encode_length_pipeline_configuration(configuration)
    if not configuration.measurement_enabled:
        if measurement is not None:
            raise ValueError("disabled length configuration cannot have a measurement")
        return _summary(configuration, None, None)
    if not isinstance(measurement, LengthPipelineMeasurement):
        raise ValueError("enabled length configuration requires a completed measurement")
    if measurement.configuration_sha256 != configuration.sha256:
        raise ValueError("length measurement configuration digest differs")
    features_document = encode_length_features(measurement.features)
    if measurement.features.provenance.measurement_context != configuration.measurement_context:
        raise ValueError("length measurement features differ from the configured measurement context")
    expected = length_estimation_state(
        measurement.features,
        configuration.model,
        measurement_enabled=True,
        evidence_domain=configuration.evidence_domain,
    )
    if measurement.estimate != expected:
        raise ValueError("length measurement estimate differs from features and configuration")
    return _summary(configuration, measurement, features_document)


def _summary(
    configuration: LengthPipelineConfiguration,
    measurement: LengthPipelineMeasurement | None,
    features_document: dict[str, object] | None,
) -> dict[str, object]:
    estimate = measurement.estimate if measurement is not None else None
    features = measurement.features if measurement is not None else None
    return {
        "estimated_total_repeat_count": None if estimate is None else estimate.estimated_total_repeat_count,
        "length_estimation_status": "disabled" if estimate is None else estimate.status,
        "length_calibration_id": configuration.candidate_id,
        "length_estimation_reasons": [] if estimate is None else list(estimate.reasons),
        "length_features": features_document,
        "length_features_sha256": None if features is None else features.sha256,
        "length_configuration_sha256": configuration.sha256,
        "length_context_sha256": configuration.context_sha256,
        "length_annotation_sha256": None if configuration.annotation is None else configuration.annotation.sha256,
        "length_counting_policy_sha256": (
            None
            if configuration.measurement_context is None
            else configuration.measurement_context.counting_policy_sha256
        ),
        "length_model_sha256": None if configuration.model is None else configuration.model.sha256,
        "length_model_bundle_sha256": configuration.model_bundle_sha256,
        "length_portable_approval_sha256": configuration.portable_approval_sha256,
    }
