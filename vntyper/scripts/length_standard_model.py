"""Strict portable linear model and inference for standard MUC1 length features."""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.length_standard_features import (
    STANDARD_ACCEPTED_CONTIGS,
    STANDARD_ANNOTATION_SHA256,
    STANDARD_ASSEMBLY,
    STANDARD_FEATURE_DEFINITION_SHA256,
    STANDARD_FEATURE_ORDER,
    STANDARD_REFERENCE_LOCUS_SHA256,
    StandardLengthMeasurement,
    encode_standard_length_measurement,
)

ModelSource = Literal["packaged-research", "local-research"]
PredictionStatus = Literal["estimated", "unavailable"]
COUNT_CONVENTIONS = ("source-reported", "complete", "canonical-only")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_ROOT_FIELDS = {
    "schema_version",
    "model_version",
    "target",
    "model_kind",
    "feature_order",
    "intercept",
    "coefficients",
    "feature_bounds",
    "qc",
    "assembly",
    "accepted_contigs",
    "annotation_sha256",
    "reference_locus_sha256",
    "feature_definition_sha256",
    "model_source",
}
_TARGET_FIELDS = {"name", "count_convention", "unit"}
_QC_FIELDS = {
    "minimum_denominator_mean_depth",
    "minimum_denominator_covered_fraction",
    "minimum_denominator_supporting_fragments",
    "minimum_eligible_reads",
}


@dataclass(frozen=True)
class StandardFeatureBounds:
    """Observed training range used only to flag extrapolation."""

    minimum: float
    maximum: float


@dataclass(frozen=True)
class StandardModelQc:
    """Frozen applicability thresholds for observable standard evidence."""

    minimum_denominator_mean_depth: float
    minimum_denominator_covered_fraction: float
    minimum_denominator_supporting_fragments: int
    minimum_eligible_reads: int


@dataclass(frozen=True)
class StandardLengthModel:
    """Immutable research linear model with no executable artifact content."""

    model_version: str
    target_name: str
    count_convention: str
    unit: str
    model_kind: str
    feature_order: tuple[str, ...]
    intercept: float
    coefficients: tuple[float, ...]
    feature_bounds: Mapping[str, StandardFeatureBounds]
    qc: StandardModelQc
    assembly: str
    accepted_contigs: tuple[str, ...]
    annotation_sha256: str
    reference_locus_sha256: str
    feature_definition_sha256: str
    model_source: ModelSource
    sha256: str


@dataclass(frozen=True)
class StandardLengthPrediction:
    """One unrounded standard estimate, stable failure reasons, and range warnings."""

    status: PredictionStatus
    estimated_repeat_count: float | None
    reasons: tuple[str, ...]
    warnings: tuple[str, ...]
    measurement_sha256: str
    model_sha256: str


def _number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"standard length model {label} must be finite numeric")
    try:
        result = float(value)
    except OverflowError:
        raise ValueError(f"standard length model {label} must be finite numeric") from None
    if not math.isfinite(result):
        raise ValueError(f"standard length model {label} must be finite numeric")
    return result


def _positive_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"standard length model {label} must be a positive integer")
    return value


def _model_payload(model: StandardLengthModel) -> dict[str, object]:
    return {
        "schema_version": "standard-length-linear-model-v1",
        "model_version": model.model_version,
        "target": {
            "name": model.target_name,
            "count_convention": model.count_convention,
            "unit": model.unit,
        },
        "model_kind": model.model_kind,
        "feature_order": list(model.feature_order),
        "intercept": model.intercept,
        "coefficients": list(model.coefficients),
        "feature_bounds": {
            name: {"minimum": model.feature_bounds[name].minimum, "maximum": model.feature_bounds[name].maximum}
            for name in STANDARD_FEATURE_ORDER
        },
        "qc": {
            "minimum_denominator_mean_depth": model.qc.minimum_denominator_mean_depth,
            "minimum_denominator_covered_fraction": model.qc.minimum_denominator_covered_fraction,
            "minimum_denominator_supporting_fragments": model.qc.minimum_denominator_supporting_fragments,
            "minimum_eligible_reads": model.qc.minimum_eligible_reads,
        },
        "assembly": model.assembly,
        "accepted_contigs": list(model.accepted_contigs),
        "annotation_sha256": model.annotation_sha256,
        "reference_locus_sha256": model.reference_locus_sha256,
        "feature_definition_sha256": model.feature_definition_sha256,
        "model_source": model.model_source,
    }


def decode_standard_length_model(value: object) -> StandardLengthModel:
    """Decode the exact standard linear model JSON contract."""
    if not isinstance(value, Mapping) or set(value) != _ROOT_FIELDS:
        raise ValueError("standard length model fields differ")
    if value["schema_version"] != "standard-length-linear-model-v1":
        raise ValueError("standard length model schema is unsupported")
    version = value["model_version"]
    if not isinstance(version, str) or not version or version.strip() != version:
        raise ValueError("standard length model version must be non-empty trimmed text")
    target = value["target"]
    if not isinstance(target, Mapping) or set(target) != _TARGET_FIELDS:
        raise ValueError("standard length model target fields differ")
    if target["name"] != "source-reported-diploid-repeat-count" or target["unit"] != "repeat_units":
        raise ValueError("standard length model target is unsupported")
    convention = target["count_convention"]
    if convention not in COUNT_CONVENTIONS:
        raise ValueError("standard length model count convention is unsupported")
    if value["model_kind"] != "linear-standard13" or value["feature_order"] != list(STANDARD_FEATURE_ORDER):
        raise ValueError("standard length model kind or feature order differs")
    coefficients = value["coefficients"]
    if not isinstance(coefficients, list) or len(coefficients) != len(STANDARD_FEATURE_ORDER):
        raise ValueError("standard length model coefficients differ from feature order")
    bounds = value["feature_bounds"]
    if not isinstance(bounds, Mapping) or set(bounds) != set(STANDARD_FEATURE_ORDER):
        raise ValueError("standard length model feature bounds differ from feature order")
    decoded_bounds: dict[str, StandardFeatureBounds] = {}
    for name in STANDARD_FEATURE_ORDER:
        item = bounds[name]
        if not isinstance(item, Mapping) or set(item) != {"minimum", "maximum"}:
            raise ValueError("standard length model feature bound fields differ")
        minimum = _number(item["minimum"], "feature bound minimum")
        maximum = _number(item["maximum"], "feature bound maximum")
        if minimum > maximum:
            raise ValueError("standard length model feature bound minimum exceeds maximum")
        decoded_bounds[name] = StandardFeatureBounds(minimum, maximum)
    qc = value["qc"]
    if not isinstance(qc, Mapping) or set(qc) != _QC_FIELDS:
        raise ValueError("standard length model QC fields differ")
    minimum_mean = _number(qc["minimum_denominator_mean_depth"], "minimum denominator mean depth")
    minimum_covered = _number(qc["minimum_denominator_covered_fraction"], "minimum denominator covered fraction")
    if minimum_mean <= 0 or not 0 < minimum_covered <= 1:
        raise ValueError("standard length model QC thresholds are out of range")
    if value["assembly"] != STANDARD_ASSEMBLY or value["accepted_contigs"] != list(STANDARD_ACCEPTED_CONTIGS):
        raise ValueError("standard length model assembly or contigs differ")
    for field, expected in (
        ("annotation_sha256", STANDARD_ANNOTATION_SHA256),
        ("reference_locus_sha256", STANDARD_REFERENCE_LOCUS_SHA256),
        ("feature_definition_sha256", STANDARD_FEATURE_DEFINITION_SHA256),
    ):
        if value[field] != expected:
            raise ValueError(f"standard length model {field} differs")
    source = value["model_source"]
    if not isinstance(source, str) or source not in {"packaged-research", "local-research"}:
        raise ValueError("standard length model source is unsupported")
    model = StandardLengthModel(
        cast(str, version),
        "source-reported-diploid-repeat-count",
        cast(str, convention),
        "repeat_units",
        "linear-standard13",
        STANDARD_FEATURE_ORDER,
        _number(value["intercept"], "intercept"),
        tuple(_number(item, "coefficient") for item in coefficients),
        MappingProxyType(decoded_bounds),
        StandardModelQc(
            minimum_mean,
            minimum_covered,
            _positive_int(qc["minimum_denominator_supporting_fragments"], "minimum denominator support"),
            _positive_int(qc["minimum_eligible_reads"], "minimum eligible reads"),
        ),
        STANDARD_ASSEMBLY,
        STANDARD_ACCEPTED_CONTIGS,
        STANDARD_ANNOTATION_SHA256,
        STANDARD_REFERENCE_LOCUS_SHA256,
        STANDARD_FEATURE_DEFINITION_SHA256,
        cast(ModelSource, source),
        canonical_sha256(value),
    )
    return model


def encode_standard_length_model(model: StandardLengthModel) -> dict[str, object]:
    """Project a typed standard model after canonical content revalidation."""
    if not isinstance(model, StandardLengthModel):
        raise ValueError("standard length model must be typed")
    payload = _model_payload(model)
    decoded = decode_standard_length_model(payload)
    if decoded != model:
        raise ValueError("standard length model typed content or digest differs")
    return payload


def load_standard_length_model(path: Path | None = None) -> StandardLengthModel:
    """Load the packaged standard research model or one explicit local model."""
    selected = (
        Path(__file__).resolve().parents[1] / "data" / "length" / "grch38-standard-length-model-v1.json"
        if path is None
        else path
    )
    if not isinstance(selected, Path) or not selected.is_file():
        raise ValueError("standard length model path must be an existing file")
    raw = read_regular_path(selected)
    try:
        value = load_strict_json_object(raw)
    except (UnicodeDecodeError, ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"standard length model JSON is invalid: {error}") from error
    if canonical_json_bytes(value) != raw:
        raise ValueError("standard length model must use canonical JSON bytes")
    model = decode_standard_length_model(value)
    expected_source = "packaged-research" if path is None else "local-research"
    if model.model_source != expected_source:
        raise ValueError(f"standard length model source must be {expected_source}")
    if path is None:
        companion_path = selected.with_suffix(selected.suffix + ".sha256")
        companion = read_regular_path(companion_path)
        expected_digest = hashlib.sha256(raw).hexdigest().encode("ascii") + b"\n"
        if companion != expected_digest or model.sha256 != expected_digest.decode("ascii").strip():
            raise ValueError("packaged standard length model digest companion differs")
    return model


def _qc_reasons(measurement: StandardLengthMeasurement, qc: StandardModelQc) -> tuple[str, ...]:
    observed = measurement.qc
    reasons: list[str] = []
    for name in ("invariant", "left_flank", "right_flank"):
        if getattr(observed, f"{name}_mean_depth") < qc.minimum_denominator_mean_depth:
            reasons.append(f"low_{name}_mean_depth")
        if getattr(observed, f"{name}_covered_fraction") < qc.minimum_denominator_covered_fraction:
            reasons.append(f"low_{name}_covered_fraction")
        support = getattr(observed, f"{name}_supporting_fragments")
        if support is None:
            reasons.append(f"missing_{name}_support")
        elif support < qc.minimum_denominator_supporting_fragments:
            reasons.append(f"low_{name}_support")
    if observed.eligible_read_count < qc.minimum_eligible_reads:
        reasons.append("low_eligible_read_count")
    return tuple(reasons)


def predict_standard_length(
    measurement: StandardLengthMeasurement, model: StandardLengthModel
) -> StandardLengthPrediction:
    """Apply one validated raw-coefficient standard model without rounding."""
    encode_standard_length_measurement(measurement)
    encode_standard_length_model(model)
    identity_reasons = []
    if measurement.assembly != model.assembly:
        identity_reasons.append("unsupported_assembly")
    if measurement.contig not in model.accepted_contigs:
        identity_reasons.append("unsupported_contig")
    if measurement.annotation_sha256 != model.annotation_sha256:
        identity_reasons.append("unsupported_annotation")
    if measurement.reference_locus_sha256 != model.reference_locus_sha256:
        identity_reasons.append("unsupported_reference_locus")
    if measurement.feature_definition_sha256 != model.feature_definition_sha256:
        identity_reasons.append("unsupported_feature_definition")
    missing = tuple(f"missing_{name}" for name in STANDARD_FEATURE_ORDER if measurement.values[name] is None)
    reasons = tuple(identity_reasons) + measurement.reasons + missing + _qc_reasons(measurement, model.qc)
    reasons = tuple(dict.fromkeys(reasons))
    if reasons:
        return StandardLengthPrediction("unavailable", None, reasons, (), measurement.sha256, model.sha256)
    values = tuple(cast(float, measurement.values[name]) for name in STANDARD_FEATURE_ORDER)
    warnings = tuple(
        f"feature_{name}_outside_training_range"
        for name, observed in zip(STANDARD_FEATURE_ORDER, values, strict=True)
        if observed < model.feature_bounds[name].minimum or observed > model.feature_bounds[name].maximum
    )
    prediction = math.fsum(
        (
            model.intercept,
            *(coefficient * value for coefficient, value in zip(model.coefficients, values, strict=True)),
        )
    )
    if not math.isfinite(prediction):
        return StandardLengthPrediction(
            "unavailable", None, ("nonfinite_prediction",), warnings, measurement.sha256, model.sha256
        )
    if prediction <= 0:
        return StandardLengthPrediction(
            "unavailable", None, ("nonpositive_prediction",), warnings, measurement.sha256, model.sha256
        )
    return StandardLengthPrediction("estimated", prediction, (), warnings, measurement.sha256, model.sha256)
