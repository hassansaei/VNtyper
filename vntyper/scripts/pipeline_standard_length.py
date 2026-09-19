"""Automatic, locus-only research length estimation inside alignment lifetime."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.canonical_json import canonical_sha256

if TYPE_CHECKING:
    from vntyper.scripts.length_standard_features import StandardLengthMeasurement
    from vntyper.scripts.length_standard_model import ModelSource, StandardLengthModel, StandardLengthPrediction

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class StandardLengthConfiguration:
    """Research model selection; this object carries no approval authority."""

    enabled: bool
    source: str | None
    model: StandardLengthModel | None
    sha256: str
    warning_threshold: float = 110.0


def encode_standard_length_configuration(configuration: StandardLengthConfiguration) -> dict[str, object]:
    """Revalidate the model and its path-free resume identity.

    Args:
        configuration: Resolved immutable research configuration.

    Returns:
        Canonicalizable configuration identity.

    Raises:
        ValueError: If fields, model source, model content or digest were replaced.
    """
    if not isinstance(configuration, StandardLengthConfiguration) or not isinstance(configuration.enabled, bool):
        raise ValueError("standard length configuration must be resolved")
    if (
        isinstance(configuration.warning_threshold, bool)
        or not isinstance(configuration.warning_threshold, (int, float))
        or not math.isfinite(configuration.warning_threshold)
        or configuration.warning_threshold <= 0
    ):
        raise ValueError("standard length warning threshold must be a positive number")
    if configuration.enabled:
        from vntyper.scripts.length_standard_model import encode_standard_length_model

        if configuration.model is None:
            raise ValueError("enabled standard length configuration requires a model")
        encode_standard_length_model(configuration.model)
        if configuration.source != configuration.model.model_source:
            raise ValueError("standard length configuration source differs from model")
    elif configuration.model is not None or configuration.source is not None:
        raise ValueError("disabled standard length configuration cannot contain a model")
    document = {
        "schema_version": "standard-length-configuration-v1",
        "enabled": configuration.enabled,
        "source": configuration.source,
        "model_sha256": None if configuration.model is None else configuration.model.sha256,
        "warning_threshold": configuration.warning_threshold,
    }
    if configuration.sha256 != canonical_sha256(document):
        raise ValueError("standard length configuration digest differs from its model and settings")
    return document


def _load_model(path: Path, *, expected_source: ModelSource | None = None) -> StandardLengthModel:
    from vntyper.scripts.length_standard_model import load_standard_length_model

    return load_standard_length_model(path, expected_source=expected_source)


def resolve_standard_length_configuration(
    config: Mapping[str, object],
    *,
    enabled: bool | None,
    model_path: Path | None,
    approved_enabled: bool,
    warning_threshold: float | None = None,
) -> StandardLengthConfiguration:
    """Resolve configured default, explicit overrides and approved-path precedence.

    Args:
        config: Loaded pipeline configuration; absent section preserves older behavior.
        enabled: Explicit CLI override, or None to use configuration.
        model_path: Optional locally fitted research JSON model.
        approved_enabled: Whether explicit measurement/approved inference owns length output.
        warning_threshold: Optional sensitivity warning cutoff in repeat units (default: 110.0).

    Returns:
        Frozen research configuration with a path-free resume identity.

    Raises:
        ValueError: On malformed settings or conflicting explicit modes.
    """
    section = config.get("length_estimation", {})
    if not isinstance(section, Mapping):
        raise ValueError("length_estimation configuration must be an object")
    default = section.get("enabled", False)
    if not isinstance(default, bool) or (enabled is not None and not isinstance(enabled, bool)):
        raise ValueError("standard length enabled setting must be boolean")
    raw_threshold = section.get("warning_threshold", 110.0) if warning_threshold is None else warning_threshold
    if (
        isinstance(raw_threshold, bool)
        or not isinstance(raw_threshold, (int, float))
        or not math.isfinite(raw_threshold)
        or raw_threshold <= 0
    ):
        raise ValueError("standard length warning threshold must be a positive number")
    resolved_threshold = float(raw_threshold)
    if model_path is not None and not isinstance(model_path, Path):
        raise ValueError("standard length model path must be a Path")
    if (approved_enabled and (enabled is True or model_path is not None)) or (
        model_path is not None and enabled is False
    ):
        raise ValueError("explicit standard length options conflict with the selected length mode")
    reference_data = config.get("reference_data")
    reference_model = (
        reference_data.get("standard_length_model_grch38") if isinstance(reference_data, Mapping) else None
    )
    reference_path = Path(str(reference_model)) if reference_model else None
    resolved_path = model_path or reference_path
    active = not approved_enabled and (model_path is not None or (default if enabled is None else enabled))
    if active:
        if resolved_path is None:
            raise ValueError("standard length model reference is not configured in reference_data")
        expected_source: ModelSource | None = "packaged-research" if model_path is None else None
        model = _load_model(resolved_path, expected_source=expected_source)
        source = model.model_source
    else:
        model = None
        source = None
    identity = {
        "schema_version": "standard-length-configuration-v1",
        "enabled": active,
        "source": source,
        "model_sha256": None if model is None else model.sha256,
        "warning_threshold": resolved_threshold,
    }
    return StandardLengthConfiguration(active, source, model, canonical_sha256(identity), resolved_threshold)


def standard_length_summary(
    configuration: StandardLengthConfiguration,
    measurement: StandardLengthMeasurement | None,
    prediction: StandardLengthPrediction | None,
    reasons: tuple[str, ...] = (),
) -> dict[str, object]:
    """Project research output without granting approved calibration identity.

    Args:
        configuration: Enabled research model selection.
        measurement: Completed observation, or None when measurement is unavailable.
        prediction: Prediction bound to measurement/model, or None on measurement failure.
        reasons: Explicit failure codes when no prediction could be attempted.

    Returns:
        Summary fields consumed by the report and recorded with the pipeline run.

    Raises:
        ValueError: When a disabled model or inconsistent prediction is supplied.
    """
    if not configuration.enabled or configuration.model is None:
        raise ValueError("standard length summary requires an enabled model")
    encode_standard_length_configuration(configuration)
    if prediction is not None:
        if measurement is None or prediction.measurement_sha256 != measurement.sha256:
            raise ValueError("standard length prediction measurement differs")
        if prediction.model_sha256 != configuration.model.sha256:
            raise ValueError("standard length prediction model differs")
        if reasons:
            raise ValueError("standard length prediction cannot have separate failure reasons")
    elif not reasons:
        raise ValueError("unavailable standard length measurement requires reasons")
    feature_document = None
    if measurement is not None:
        from vntyper.scripts.length_standard_features import encode_standard_length_measurement

        feature_document = encode_standard_length_measurement(measurement)
    return {
        "length_estimation_status": "unavailable" if prediction is None else prediction.status,
        "estimated_total_repeat_count": None if prediction is None else prediction.estimated_repeat_count,
        "length_warning_threshold": configuration.warning_threshold,
        "length_estimation_reasons": list(reasons if prediction is None else prediction.reasons),
        "length_estimation_warnings": [] if prediction is None else list(prediction.warnings),
        "length_model_source": configuration.source,
        "length_model_evidence_status": "research-development",
        "length_count_convention": configuration.model.count_convention,
        "length_model_sha256": configuration.model.sha256,
        "length_standard_features": feature_document,
        "length_standard_features_sha256": None if measurement is None else measurement.sha256,
        "length_configuration_sha256": configuration.sha256,
        "length_calibration_id": None,
        "length_portable_approval_sha256": None,
        "length_model_bundle_sha256": None,
        "length_features": None,
        "length_features_sha256": None,
    }


@dataclass
class StandardLengthRunner:
    """Consume one retained alignment plan and retain its research summary."""

    configuration: StandardLengthConfiguration
    assembly: str
    bwa_reference: str | Path | None
    project_root: Path
    summary: dict[str, object] = field(default_factory=dict)
    _called: bool = False

    def __call__(self, plan: AlignmentPlan) -> None:
        """Extract once; unavailable optional measurement leaves genotyping intact."""
        encode_standard_length_configuration(self.configuration)
        if self._called:
            raise ValueError("standard length callback must run exactly once")
        self._called = True
        if self.assembly not in {"hg38", "GRCh38", "hg38_ensembl", "hg38_ncbi"}:
            self.summary = standard_length_summary(self.configuration, None, None, ("unsupported-assembly",))
            return
        from vntyper.scripts.pipeline_length_execution import length_reference_path

        try:
            reference = length_reference_path(plan, self.bwa_reference, self.project_root)
        except ValueError:
            self.summary = standard_length_summary(self.configuration, None, None, ("reference-unavailable",))
            return
        self._measure(plan, reference)

    def _measure(self, plan: AlignmentPlan, reference: Path) -> None:
        from vntyper.scripts.length_standard_io import read_standard_length_features
        from vntyper.scripts.length_standard_model import predict_standard_length

        if self.configuration.model is None:
            raise ValueError("standard length runner requires a model")
        try:
            measurement = read_standard_length_features(
                Path(plan.view_path).resolve(),
                reference.resolve(),
                assembly=self.assembly,
                index_path=None if plan.stable_index_path is None else Path(plan.stable_index_path).resolve(),
            )
            prediction = predict_standard_length(
                measurement, self.configuration.model, warning_threshold=self.configuration.warning_threshold
            )
        except (OSError, RuntimeError, ValueError) as error:
            logger.warning("Standard VNTR length measurement unavailable: %s", error)
            self.summary = standard_length_summary(self.configuration, None, None, ("measurement-unavailable",))
            return
        self.summary = standard_length_summary(self.configuration, measurement, prediction)
