"""Strict non-executable research models for total diploid repeat count."""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_candidate import (
    CandidateApplicability,
    CandidateProducer,
    candidate_applicability_document,
    candidate_producer_document,
    decode_candidate_applicability,
    decode_candidate_producer,
)
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

LengthModelKind = Literal["affine-A", "affine-F"]
TARGET_BOUNDARY_DEFINITION = "complete-core-plus-invariant-units-v1"
_FIELDS = {
    "schema_version",
    "target",
    "unit",
    "model_kind",
    "feature_order",
    "intercept",
    "coefficients",
    "annotation_sha256",
    "counting_policy_sha256",
    "applicability",
    "qc",
    "feature_bounds",
    "study_sha256",
    "training_evidence_sha256",
    "producer",
}
_QC_FIELDS = {
    "minimum_denominator_mean_depth",
    "minimum_denominator_covered_fraction",
    "minimum_denominator_supporting_fragments",
    "fragment_evidence_kind",
}
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))


@dataclass(frozen=True)
class FeatureBounds:
    """Inclusive finite bounds learned from the training-only feature range."""

    minimum: float
    maximum: float


@dataclass(frozen=True)
class LengthModel:
    """Immutable affine research model; decoding grants no promotion authority."""

    target_name: str
    boundary_definition: str
    unit: str
    model_kind: LengthModelKind
    feature_order: tuple[str, ...]
    intercept: float
    coefficients: tuple[float, ...]
    annotation_sha256: str
    counting_policy_sha256: str
    applicability: CandidateApplicability
    qc: Mapping[str, int | float | str]
    feature_bounds: Mapping[str, FeatureBounds]
    study_sha256: str
    training_evidence_sha256: str
    producer: CandidateProducer
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"length model {label} fields differ from the closed contract")
    return value


def _digest(value: object, field: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"length model {field} must be a lowercase SHA256 digest")
    return value


def _number(value: object, field: str, *, positive: bool = False, maximum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(f"length model {field} must be finite numeric")
    try:
        result = float(value)
    except OverflowError:
        _fail(f"length model {field} must be finite numeric")
    if not math.isfinite(result):
        _fail(f"length model {field} must be finite numeric")
    if positive and result <= 0:
        _fail(f"length model {field} must be positive")
    if maximum is not None and result > maximum:
        _fail(f"length model {field} exceeds its maximum")
    return result


def _qc(value: object) -> Mapping[str, int | float | str]:
    raw = _object(value, _QC_FIELDS, "QC")
    support = raw["minimum_denominator_supporting_fragments"]
    if isinstance(support, bool) or not isinstance(support, int) or support <= 0:
        _fail("length model minimum_denominator_supporting_fragments must be a positive integer")
    if raw["fragment_evidence_kind"] != "read-pair-identity-qc-proxy":
        _fail("length model fragment_evidence_kind must be read-pair-identity-qc-proxy")
    return MappingProxyType(
        {
            "minimum_denominator_mean_depth": _number(
                raw["minimum_denominator_mean_depth"], "minimum_denominator_mean_depth", positive=True
            ),
            "minimum_denominator_covered_fraction": _number(
                raw["minimum_denominator_covered_fraction"],
                "minimum_denominator_covered_fraction",
                positive=True,
                maximum=1,
            ),
            "minimum_denominator_supporting_fragments": support,
            "fragment_evidence_kind": "read-pair-identity-qc-proxy",
        }
    )


def decode_length_model_qc(value: object) -> Mapping[str, int | float | str]:
    """Decode the shared closed denominator-QC contract.

    Args:
        value: Parsed JSON-compatible QC object.

    Returns:
        Immutable validated model QC.

    Raises:
        ValueError: If fields, thresholds, or evidence identity are invalid.
    """
    return _qc(value)


def length_model_qc_document(qc: Mapping[str, int | float | str]) -> dict[str, object]:
    """Project immutable denominator QC after canonical revalidation.

    Args:
        qc: QC returned by :func:`decode_length_model_qc`.

    Returns:
        Fresh closed JSON-compatible QC object.

    Raises:
        ValueError: If the mapping is mutable or its content is invalid.
    """
    if not isinstance(qc, _MAPPING_PROXY_TYPE):
        _fail("length model QC must use a decoded immutable mapping")
    raw: dict[str, object] = dict(qc)
    if decode_length_model_qc(raw) != qc:
        _fail("length model QC differs from its decoded contract")
    return raw


def _feature_bounds(value: object, feature_order: tuple[str, ...]) -> Mapping[str, FeatureBounds]:
    if not isinstance(value, Mapping) or set(value) != set(feature_order):
        _fail("length model feature_bounds fields must match feature_order exactly")
    result = {}
    for feature in feature_order:
        raw = _object(value[feature], {"minimum", "maximum"}, "feature bounds")
        minimum = _number(raw["minimum"], "feature bound minimum")
        maximum = _number(raw["maximum"], "feature bound maximum")
        if minimum >= maximum:
            _fail("length model feature bound minimum must be less than maximum")
        result[feature] = FeatureBounds(minimum, maximum)
    return MappingProxyType(result)


def decode_length_model(value: object) -> LengthModel:
    """Decode a closed affine model without executing artifact-supplied code.

    Args:
        value: Parsed JSON-compatible model object.

    Returns:
        Immutable research model and canonical digest.

    Raises:
        ValueError: If schema, fields, identities, values or model kind are invalid.
    """
    raw = _object(value, _FIELDS, "root")
    if raw["schema_version"] != "length-model-v1":
        _fail("length model schema_version must be length-model-v1")
    target = _object(raw["target"], {"name", "boundary_definition"}, "target")
    if target["name"] != "total_diploid_repeat_count":
        _fail("length model target name must be total_diploid_repeat_count")
    if target["boundary_definition"] != TARGET_BOUNDARY_DEFINITION:
        _fail(f"length model boundary_definition must be {TARGET_BOUNDARY_DEFINITION}")
    if raw["unit"] != "repeat_units":
        _fail("length model unit must be repeat_units")
    kind = raw["model_kind"]
    if not isinstance(kind, str):
        _fail("length model model_kind must be affine-A or affine-F")
    if kind in {"physical-A", "physical-F"}:
        _fail("physical_model_geometry_evidence_unsupported")
    if kind not in {"affine-A", "affine-F"}:
        _fail("length model model_kind must be affine-A or affine-F")
    model_kind = cast(LengthModelKind, kind)
    feature = model_kind[-1]
    if raw["feature_order"] != [feature]:
        _fail("length model feature_order must exactly match its affine model kind")
    coefficients = raw["coefficients"]
    if not isinstance(coefficients, list) or len(coefficients) != 1:
        _fail("length model coefficients must contain exactly one value")
    annotation_digest = _digest(raw["annotation_sha256"], "annotation_sha256")
    counting_digest = _digest(raw["counting_policy_sha256"], "counting_policy_sha256")
    applicability = decode_candidate_applicability(raw["applicability"], target="length")
    if applicability.counting_policy_sha256 != counting_digest:
        _fail("length model applicability counting_policy_sha256 differs from the model binding")
    return LengthModel(
        target_name="total_diploid_repeat_count",
        boundary_definition=TARGET_BOUNDARY_DEFINITION,
        unit="repeat_units",
        model_kind=model_kind,
        feature_order=(feature,),
        intercept=_number(raw["intercept"], "intercept"),
        coefficients=(_number(coefficients[0], "coefficients"),),
        annotation_sha256=annotation_digest,
        counting_policy_sha256=counting_digest,
        applicability=applicability,
        qc=_qc(raw["qc"]),
        feature_bounds=_feature_bounds(raw["feature_bounds"], (feature,)),
        study_sha256=_digest(raw["study_sha256"], "study_sha256"),
        training_evidence_sha256=_digest(raw["training_evidence_sha256"], "training_evidence_sha256"),
        producer=decode_candidate_producer(raw["producer"]),
        sha256=canonical_sha256(raw),
    )


def _model_document(model: LengthModel) -> dict[str, object]:
    return {
        "schema_version": "length-model-v1",
        "target": {"name": model.target_name, "boundary_definition": model.boundary_definition},
        "unit": model.unit,
        "model_kind": model.model_kind,
        "feature_order": list(model.feature_order),
        "intercept": model.intercept,
        "coefficients": list(model.coefficients),
        "annotation_sha256": model.annotation_sha256,
        "counting_policy_sha256": model.counting_policy_sha256,
        "applicability": candidate_applicability_document(model.applicability, target="length"),
        "qc": dict(model.qc),
        "feature_bounds": {
            name: {"minimum": bounds.minimum, "maximum": bounds.maximum}
            for name, bounds in model.feature_bounds.items()
        },
        "study_sha256": model.study_sha256,
        "training_evidence_sha256": model.training_evidence_sha256,
        "producer": candidate_producer_document(model.producer),
    }


def encode_length_model(model: LengthModel) -> dict[str, object]:
    """Project a validated immutable model into canonicalizable JSON content.

    Args:
        model: Model returned by :func:`decode_length_model`.

    Returns:
        Fresh closed model object.

    Raises:
        ValueError: If typed content, immutability or its digest was forged.
    """
    if not isinstance(model, LengthModel):
        _fail("length model must be a LengthModel")
    if (
        not isinstance(model.feature_order, tuple)
        or not isinstance(model.coefficients, tuple)
        or not isinstance(model.qc, _MAPPING_PROXY_TYPE)
        or not isinstance(model.feature_bounds, _MAPPING_PROXY_TYPE)
        or any(not isinstance(bounds, FeatureBounds) for bounds in model.feature_bounds.values())
    ):
        _fail("length model must use decoded immutable collections")
    raw = _model_document(model)
    if decode_length_model(raw) != model:
        _fail("length model content or digest does not match its decoded contract")
    return raw
