"""Resolve and route length measurement without growing pipeline orchestration."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from vntyper.scripts.length_sensitivity import LengthSensitivityPolicy, apply_length_sensitivity
from vntyper.scripts.pipeline_length import (
    LengthPipelineConfiguration,
    encode_length_pipeline_configuration,
    resolve_length_pipeline_configuration,
)
from vntyper.scripts.pipeline_length_execution import LengthMeasurementRunner, length_summary_fields
from vntyper.scripts.pipeline_standard_length import StandardLengthConfiguration, StandardLengthRunner


def validate_pipeline_length(
    configuration: LengthPipelineConfiguration | None, operator_paths: tuple[str | Path, ...]
) -> LengthPipelineConfiguration:
    """Resolve legacy optional length settings and validate input ownership paths.

    Args:
        configuration: Explicit preflighted configuration, or None for disabled legacy mode.
        operator_paths: Exact operator-owned length inputs protected by pipeline guards.

    Returns:
        Validated legacy configuration.

    Raises:
        ValueError: If configuration or protected paths violate the pipeline contract.
    """
    if configuration is None:
        configuration = resolve_length_pipeline_configuration(
            measurement_enabled=False, model_path=None, annotation_path=None, context_path=None
        )
    elif not isinstance(configuration, LengthPipelineConfiguration):
        raise ValueError("pipeline length_configuration must be resolved")
    encode_length_pipeline_configuration(configuration)
    if not isinstance(operator_paths, tuple) or any(not isinstance(path, (str, Path)) for path in operator_paths):
        raise ValueError("pipeline length_operator_paths must be a tuple of paths")
    return configuration


def build_length_consumer(
    configuration: LengthPipelineConfiguration,
    standard: StandardLengthConfiguration,
    *,
    assembly: str,
    reference: str | Path | None,
    project_root: str | Path,
    samtools: str | Path,
    explicit_bam_reference: str | Path | None = None,
    approved_factory: Callable[..., LengthMeasurementRunner] = LengthMeasurementRunner,
) -> LengthMeasurementRunner | StandardLengthRunner | None:
    """Choose exactly one measurement callback inside the proven alignment lifetime.

    Args:
        configuration: Legacy approved/measurement-only path selection.
        standard: Automatic or explicitly selected research path.
        assembly: Actual pipeline coordinate system.
        reference: Actual assembly reference selected by the pipeline.
        project_root: Original working directory for relative reference paths.
        samtools: Tool required by the strict approved measurement path.
        explicit_bam_reference: Optional explicit FASTA for standard BAM measurement.
        approved_factory: Existing injectable strict runner constructor.

    Returns:
        A single selected callback, or None when both modes are disabled.

    Raises:
        ValueError: If two length paths are enabled simultaneously.
    """
    if configuration.measurement_enabled and standard.enabled:
        raise ValueError("multiple length measurement paths are enabled")
    if configuration.measurement_enabled:
        return approved_factory(
            configuration=configuration, bwa_reference=reference, project_root=project_root, samtools_path=samtools
        )
    if standard.enabled:
        return StandardLengthRunner(standard, assembly, explicit_bam_reference or reference, Path(project_root))
    return None


def completed_length_summary(
    configuration: LengthPipelineConfiguration,
    runner: LengthMeasurementRunner | StandardLengthRunner,
    *,
    approved_projector: Callable[..., dict[str, object]] = length_summary_fields,
    sensitivity_policy: LengthSensitivityPolicy | None = None,
) -> dict[str, object]:
    """Return recorded results after the selected callback consumed its plan.

    Args:
        configuration: Existing approved/measurement-only configuration.
        runner: The callback passed to the coverage stage.
        approved_projector: Existing injectable strict summary projection.
        sensitivity_policy: Length sensitivity tier policy recorded with the estimate, or
            None when the run configuration has no policy.

    Returns:
        Completed summary fields for the selected model source.

    Raises:
        ValueError: If measurement was not completed.
    """
    if isinstance(runner, StandardLengthRunner):
        if not runner.summary:
            raise ValueError("standard length callback did not complete")
        fields = dict(runner.summary)
    else:
        fields = approved_projector(configuration, runner.result)
    return apply_length_sensitivity(fields, sensitivity_policy)
