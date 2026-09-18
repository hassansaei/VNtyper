"""Training-only Bayesian linear fit for the closed standard length predictors."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import cast

import numpy as np

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
from vntyper.scripts.length_standard_model import (
    COUNT_CONVENTIONS,
    StandardLengthModel,
    decode_standard_length_model,
)


def collapse_scaled_linear_model(
    *,
    scaled_intercept: float,
    scaled_coefficients: Sequence[float],
    means: Sequence[float],
    scales: Sequence[float],
) -> tuple[float, tuple[float, ...]]:
    """Collapse standardization and a linear estimator into raw-space coefficients."""
    if not (len(scaled_coefficients) == len(means) == len(scales)) or not scales:
        raise ValueError("scaled standard length vectors must have the same non-zero length")
    values = (scaled_intercept, *scaled_coefficients, *means, *scales)
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in values):
        raise ValueError("scaled standard length parameters must be numeric")
    numeric = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in numeric) or any(float(value) <= 0 for value in scales):
        raise ValueError("scaled standard length parameters must be finite with positive scales")
    raw = tuple(
        float(coefficient) / float(scale) for coefficient, scale in zip(scaled_coefficients, scales, strict=True)
    )
    intercept = float(scaled_intercept) - math.fsum(
        coefficient * float(mean) for coefficient, mean in zip(raw, means, strict=True)
    )
    if not math.isfinite(intercept) or not all(math.isfinite(value) for value in raw):
        raise ValueError("collapsed standard length parameters must be finite")
    return intercept, raw


def _fit_raw_bayesian(matrix: Sequence[Sequence[float]], targets: Sequence[float]) -> tuple[float, tuple[float, ...]]:
    try:
        from sklearn.linear_model import BayesianRidge
        from sklearn.preprocessing import StandardScaler
    except ImportError as error:
        raise RuntimeError("standard length fitting requires the optional scikit-learn training dependency") from error
    values = np.asarray(matrix, dtype=float)
    truth = np.asarray(targets, dtype=float)
    scaler = StandardScaler()
    standardized = scaler.fit_transform(values)
    estimator = BayesianRidge(max_iter=1000)
    estimator.fit(standardized, truth)
    return collapse_scaled_linear_model(
        scaled_intercept=float(estimator.intercept_),
        scaled_coefficients=tuple(float(value) for value in estimator.coef_),
        means=tuple(float(value) for value in scaler.mean_),
        scales=tuple(float(value) for value in scaler.scale_),
    )


def standard_length_training_reasons(measurement: StandardLengthMeasurement) -> tuple[str, ...]:
    """Return all stable row-level reasons that make a measurement ineligible for fitting."""
    encode_standard_length_measurement(measurement)
    qc = measurement.qc
    reasons = list(measurement.reasons)
    for name in ("invariant", "left_flank", "right_flank"):
        if getattr(qc, f"{name}_mean_depth") < 10:
            reasons.append(f"low_{name}_mean_depth")
        if getattr(qc, f"{name}_covered_fraction") < 0.9:
            reasons.append(f"low_{name}_covered_fraction")
        support = getattr(qc, f"{name}_supporting_fragments")
        if support is None:
            reasons.append(f"missing_{name}_support")
        elif support < 100:
            reasons.append(f"low_{name}_support")
    if qc.eligible_read_count < 1:
        reasons.append("low_eligible_read_count")
    return tuple(dict.fromkeys(reasons))


def fit_standard_length_model(
    matrix: Sequence[StandardLengthMeasurement],
    targets: Sequence[float],
    *,
    count_convention: str,
) -> StandardLengthModel:
    """Fit one training-only Bayesian linear model and export raw coefficients."""
    if count_convention not in COUNT_CONVENTIONS:
        raise ValueError("standard length count convention is unsupported")
    if not isinstance(matrix, Sequence) or isinstance(matrix, (str, bytes)) or len(matrix) < 2:
        raise ValueError("standard length fitting requires at least two measurements")
    if not isinstance(targets, Sequence) or isinstance(targets, (str, bytes)) or len(targets) != len(matrix):
        raise ValueError("standard length targets must match the measurement count")
    rows = tuple(matrix)
    if any(not isinstance(row, StandardLengthMeasurement) for row in rows):
        raise ValueError("standard length matrix must contain typed measurements")
    for row in rows:
        encode_standard_length_measurement(row)
    if any(standard_length_training_reasons(row) for row in rows):
        raise ValueError("standard length training requires complete predictors passing frozen QC")
    values = tuple(tuple(row.values[name] for name in STANDARD_FEATURE_ORDER) for row in rows)
    if any(value is None for row in values for value in row):
        raise ValueError("standard length training predictors must be complete")
    numeric_values = tuple(tuple(float(cast(float, value)) for value in row) for row in values)
    numeric_targets = []
    for target in targets:
        if isinstance(target, bool) or not isinstance(target, (int, float)):
            raise ValueError("standard length targets must be positive integral values")
        try:
            numeric = float(target)
        except OverflowError:
            raise ValueError("standard length targets must be positive integral values") from None
        if not math.isfinite(numeric) or numeric <= 0 or not numeric.is_integer():
            raise ValueError("standard length targets must be positive integral values")
        numeric_targets.append(numeric)
    intercept, coefficients = _fit_raw_bayesian(numeric_values, numeric_targets)
    if len(coefficients) != len(STANDARD_FEATURE_ORDER):
        raise ValueError("standard length fitted coefficient count differs")
    bounds = {}
    for index, name in enumerate(STANDARD_FEATURE_ORDER):
        column = tuple(row[index] for row in numeric_values)
        minimum = min(column)
        maximum = max(column)
        expansion = (maximum - minimum) * 0.1
        bounds[name] = {"minimum": minimum - expansion, "maximum": maximum + expansion}
    return decode_standard_length_model(
        {
            "schema_version": "standard-length-linear-model-v1",
            "model_version": "standard13-bayesian-v1",
            "target": {
                "name": "source-reported-diploid-repeat-count",
                "count_convention": count_convention,
                "unit": "repeat_units",
            },
            "model_kind": "linear-standard13",
            "feature_order": list(STANDARD_FEATURE_ORDER),
            "intercept": intercept,
            "coefficients": list(coefficients),
            "feature_bounds": bounds,
            "qc": {
                "minimum_denominator_mean_depth": 10,
                "minimum_denominator_covered_fraction": 0.9,
                "minimum_denominator_supporting_fragments": 100,
                "minimum_eligible_reads": 1,
            },
            "assembly": STANDARD_ASSEMBLY,
            "accepted_contigs": list(STANDARD_ACCEPTED_CONTIGS),
            "annotation_sha256": STANDARD_ANNOTATION_SHA256,
            "reference_locus_sha256": STANDARD_REFERENCE_LOCUS_SHA256,
            "feature_definition_sha256": STANDARD_FEATURE_DEFINITION_SHA256,
            "model_source": "local-research",
        }
    )
