"""Training-only affine length fitting and immutable baseline contracts.

The explicit training-role tag and roster binding form a trusted typed boundary;
they do not claim cryptographic proof that labels were previously unseen.
"""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Literal, NoReturn, cast

import numpy as np

from vntyper.scripts.calibration_candidate import (
    CandidateApplicability,
    CandidateProducer,
    candidate_applicability_document,
    candidate_producer_document,
    decode_candidate_applicability,
    decode_candidate_producer,
)
from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    length_eligible_roster_document,
)
from vntyper.scripts.calibration_length_protocol import LengthProtocol, length_protocol_document
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_estimation import EvidenceDomain, assess_length_feature
from vntyper.scripts.length_features import LengthFeatures, encode_length_features
from vntyper.scripts.length_model import (
    TARGET_BOUNDARY_DEFINITION,
    LengthModel,
    LengthModelKind,
    decode_length_model,
    decode_length_model_qc,
)

logger = logging.getLogger(__name__)

TrainingRole = Literal["training"]
FitStatus = Literal["fitted", "ineligible"]
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_QC_FIELDS = {
    "minimum_denominator_mean_depth",
    "minimum_denominator_covered_fraction",
    "minimum_denominator_supporting_fragments",
}
_METADATA_FIELDS = {
    "schema_version",
    "study_sha256",
    "training_evidence_sha256",
    "training_roster_sha256",
    "annotation_sha256",
    "counting_policy_sha256",
    "applicability",
    "qc",
    "producer",
    "maximum_condition_number",
}
_BASELINE_FIELDS = {
    "schema_version",
    "target",
    "unit",
    "baseline_kind",
    "mean_total_repeat_count",
    "study_sha256",
    "training_evidence_sha256",
    "training_roster_sha256",
    "independent_group_count",
}
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))


@dataclass(frozen=True)
class LengthTrainingRow:
    """One exact eligible training representative and its matched total truth."""

    key: str
    group_key: str
    role: TrainingRole
    features: LengthFeatures
    truth_boundary_definition: str
    total_truth_repeat_count: float
    evidence_domain: EvidenceDomain


@dataclass(frozen=True)
class LengthTrainingMetadata:
    """Frozen training bindings and numerical stability policy."""

    study_sha256: str
    training_evidence_sha256: str
    training_roster_sha256: str
    annotation_sha256: str
    counting_policy_sha256: str
    applicability: CandidateApplicability
    qc: Mapping[str, int | float]
    producer: CandidateProducer
    maximum_condition_number: float
    sha256: str


@dataclass(frozen=True)
class LengthBaseline:
    """Hash-bound training-mean baseline, frozen before selection."""

    mean_total_repeat_count: float
    study_sha256: str
    training_evidence_sha256: str
    training_roster_sha256: str
    independent_group_count: int
    sha256: str


@dataclass(frozen=True)
class LengthFitOutcome:
    """One declared hypothesis fit or an explicit whole-candidate failure."""

    candidate_id: str
    model_kind: str
    status: FitStatus
    reasons: tuple[str, ...]
    model: LengthModel | None


@dataclass(frozen=True)
class LengthTrainingResult:
    """Training-only baseline and diagnostics for every declared hypothesis."""

    metadata_sha256: str
    protocol_sha256: str
    training_roster_sha256: str
    baseline: LengthBaseline
    outcomes: tuple[LengthFitOutcome, ...]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"length training {label} fields differ from the closed contract")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        _fail(f"length training {label} must be non-empty text")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"length training {label} must be a lowercase SHA256 digest")
    return value


def _positive_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(f"length training {label} must be finite and positive")
    try:
        result = float(value)
    except OverflowError:
        _fail(f"length training {label} must be finite and positive")
    if not math.isfinite(result) or result <= 0:
        _fail(f"length training {label} must be finite and positive")
    return result


def _metadata_qc(value: object) -> Mapping[str, int | float]:
    raw = _object(value, _QC_FIELDS, "QC")
    model_qc = decode_length_model_qc({**raw, "fragment_evidence_kind": "read-pair-identity-qc-proxy"})
    return MappingProxyType({name: cast(int | float, model_qc[name]) for name in _QC_FIELDS})


def decode_length_training_metadata(value: object, training_roster: LengthEligibleRoster) -> LengthTrainingMetadata:
    """Decode metadata and bind its exact, separately frozen training roster.

    Args:
        value: Parsed JSON-compatible metadata object.
        training_roster: Roster frozen for training, independent of evaluation.

    Returns:
        Immutable metadata with a canonical digest.

    Raises:
        ValueError: If fields, identities, policies, or roster binding are invalid.
    """
    length_eligible_roster_document(training_roster)
    raw = _object(value, _METADATA_FIELDS, "metadata")
    if raw["schema_version"] != "length-training-metadata-v1":
        _fail("length training metadata schema_version must be length-training-metadata-v1")
    roster_sha256 = _digest(raw["training_roster_sha256"], "training roster digest")
    if roster_sha256 != training_roster.sha256:
        _fail("length training metadata training roster digest differs from the supplied roster")
    counting_sha256 = _digest(raw["counting_policy_sha256"], "counting policy digest")
    applicability = decode_candidate_applicability(raw["applicability"], target="length")
    if applicability.counting_policy_sha256 != counting_sha256:
        _fail("length training applicability counting policy differs from metadata")
    maximum_condition = _positive_number(raw["maximum_condition_number"], "maximum condition number")
    if maximum_condition < 1:
        _fail("length training maximum condition number must be at least one")
    return LengthTrainingMetadata(
        study_sha256=_digest(raw["study_sha256"], "study digest"),
        training_evidence_sha256=_digest(raw["training_evidence_sha256"], "training evidence digest"),
        training_roster_sha256=roster_sha256,
        annotation_sha256=_digest(raw["annotation_sha256"], "annotation digest"),
        counting_policy_sha256=counting_sha256,
        applicability=applicability,
        qc=_metadata_qc(raw["qc"]),
        producer=decode_candidate_producer(raw["producer"]),
        maximum_condition_number=maximum_condition,
        sha256=canonical_sha256(raw),
    )


def _metadata_document(metadata: LengthTrainingMetadata) -> dict[str, object]:
    return {
        "schema_version": "length-training-metadata-v1",
        "study_sha256": metadata.study_sha256,
        "training_evidence_sha256": metadata.training_evidence_sha256,
        "training_roster_sha256": metadata.training_roster_sha256,
        "annotation_sha256": metadata.annotation_sha256,
        "counting_policy_sha256": metadata.counting_policy_sha256,
        "applicability": candidate_applicability_document(metadata.applicability, target="length"),
        "qc": dict(metadata.qc),
        "producer": candidate_producer_document(metadata.producer),
        "maximum_condition_number": metadata.maximum_condition_number,
    }


def length_training_metadata_document(metadata: LengthTrainingMetadata) -> dict[str, object]:
    """Project metadata after validating its immutable content and digest.

    Args:
        metadata: Decoded training metadata.

    Returns:
        Fresh closed JSON-compatible metadata object.

    Raises:
        ValueError: If typed content is mutable, invalid, or stale.
    """
    if not isinstance(metadata, LengthTrainingMetadata):
        _fail("length training metadata must be LengthTrainingMetadata")
    if not isinstance(metadata.qc, _MAPPING_PROXY_TYPE):
        _fail("length training metadata must use immutable QC")
    raw = _metadata_document(metadata)
    for name in (
        "study_sha256",
        "training_evidence_sha256",
        "training_roster_sha256",
        "annotation_sha256",
        "counting_policy_sha256",
    ):
        _digest(raw[name], name)
    if metadata.applicability.counting_policy_sha256 != metadata.counting_policy_sha256:
        _fail("length training applicability counting policy differs from metadata")
    if _metadata_qc(raw["qc"]) != metadata.qc:
        _fail("length training metadata QC differs from its decoded contract")
    maximum_condition = _positive_number(raw["maximum_condition_number"], "maximum condition number")
    if maximum_condition < 1:
        _fail("length training maximum condition number must be at least one")
    if canonical_sha256(raw) != metadata.sha256:
        _fail("length training metadata digest differs from its canonical content")
    return raw


def decode_length_baseline(value: object) -> LengthBaseline:
    """Decode the closed training-mean baseline artifact.

    Args:
        value: Parsed JSON-compatible baseline object.

    Returns:
        Immutable baseline with canonical digest.

    Raises:
        ValueError: If fields, identities, target, or values are invalid.
    """
    raw = _object(value, _BASELINE_FIELDS, "baseline")
    target = _object(raw["target"], {"name", "boundary_definition"}, "baseline target")
    if target != {
        "name": "total_diploid_repeat_count",
        "boundary_definition": TARGET_BOUNDARY_DEFINITION,
    }:
        _fail("length training baseline target differs from total diploid repeat count")
    if raw["schema_version"] != "length-training-baseline-v1":
        _fail("length training baseline schema_version must be length-training-baseline-v1")
    if raw["unit"] != "repeat_units" or raw["baseline_kind"] != "training-mean-v1":
        _fail("length training baseline unit or kind is invalid")
    count = raw["independent_group_count"]
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
        _fail("length training baseline independent group count must be positive integer")
    return LengthBaseline(
        mean_total_repeat_count=_positive_number(raw["mean_total_repeat_count"], "baseline mean"),
        study_sha256=_digest(raw["study_sha256"], "baseline study digest"),
        training_evidence_sha256=_digest(raw["training_evidence_sha256"], "baseline training evidence digest"),
        training_roster_sha256=_digest(raw["training_roster_sha256"], "baseline training roster digest"),
        independent_group_count=count,
        sha256=canonical_sha256(raw),
    )


def _baseline_document(baseline: LengthBaseline) -> dict[str, object]:
    return {
        "schema_version": "length-training-baseline-v1",
        "target": {
            "name": "total_diploid_repeat_count",
            "boundary_definition": TARGET_BOUNDARY_DEFINITION,
        },
        "unit": "repeat_units",
        "baseline_kind": "training-mean-v1",
        "mean_total_repeat_count": baseline.mean_total_repeat_count,
        "study_sha256": baseline.study_sha256,
        "training_evidence_sha256": baseline.training_evidence_sha256,
        "training_roster_sha256": baseline.training_roster_sha256,
        "independent_group_count": baseline.independent_group_count,
    }


def length_baseline_document(baseline: LengthBaseline) -> dict[str, object]:
    """Project and revalidate a frozen training-mean baseline.

    Args:
        baseline: Baseline returned by training or its decoder.

    Returns:
        Fresh closed JSON-compatible baseline object.

    Raises:
        ValueError: If typed content or digest is invalid.
    """
    if not isinstance(baseline, LengthBaseline):
        _fail("length training baseline must be LengthBaseline")
    raw = _baseline_document(baseline)
    if decode_length_baseline(raw) != baseline:
        _fail("length training baseline differs from its canonical content or digest")
    return raw


def _training_rows(
    rows: Sequence[LengthTrainingRow], training_roster: LengthEligibleRoster
) -> tuple[LengthTrainingRow, ...]:
    if not isinstance(rows, (tuple, list)) or not rows:
        _fail("length training rows must be a non-empty typed sequence")
    checked: list[LengthTrainingRow] = []
    keys: set[str] = set()
    groups: set[str] = set()
    for row in rows:
        if not isinstance(row, LengthTrainingRow):
            _fail("length training rows must contain LengthTrainingRow values")
        _text(row.key, "row key")
        _text(row.group_key, "row group key")
        if row.key in keys:
            _fail("length training rows contain a duplicate key")
        if row.group_key in groups:
            _fail("length training rows contain a duplicate group")
        keys.add(row.key)
        groups.add(row.group_key)
        if row.role != "training":
            _fail("length training rows must have the training role")
        if row.truth_boundary_definition != TARGET_BOUNDARY_DEFINITION:
            _fail("length training truth boundary definition is incompatible")
        _positive_number(row.total_truth_repeat_count, "total truth repeat count")
        if not isinstance(row.evidence_domain, str) or row.evidence_domain not in {"synthetic", "external"}:
            _fail("length training evidence domain must be synthetic or external")
        encode_length_features(row.features)
        if row.features.manifest_key != row.key:
            _fail("length training feature manifest key differs from its row key")
        checked.append(row)
    expected = {(member.key, member.group_key) for member in training_roster.members}
    observed = {(row.key, row.group_key) for row in checked}
    if observed != expected:
        _fail("length training rows must match the frozen training roster exactly")
    return tuple(sorted(checked, key=lambda row: row.group_key))


def _model_qc(metadata: LengthTrainingMetadata) -> Mapping[str, int | float | str]:
    return decode_length_model_qc({**metadata.qc, "fragment_evidence_kind": "read-pair-identity-qc-proxy"})


def _baseline(rows: tuple[LengthTrainingRow, ...], metadata: LengthTrainingMetadata) -> LengthBaseline:
    baseline = LengthBaseline(
        mean_total_repeat_count=math.fsum(row.total_truth_repeat_count / len(rows) for row in rows),
        study_sha256=metadata.study_sha256,
        training_evidence_sha256=metadata.training_evidence_sha256,
        training_roster_sha256=metadata.training_roster_sha256,
        independent_group_count=len(rows),
        sha256="",
    )
    return replace(baseline, sha256=canonical_sha256(_baseline_document(baseline)))


def _failure(candidate_id: str, model_kind: str, *reasons: str) -> LengthFitOutcome:
    return LengthFitOutcome(candidate_id, model_kind, "ineligible", tuple(reasons), None)


def _affine_fit(
    candidate_id: str,
    model_kind: LengthModelKind,
    rows: tuple[LengthTrainingRow, ...],
    metadata: LengthTrainingMetadata,
    qc: Mapping[str, int | float | str],
) -> LengthFitOutcome:
    feature_name = cast(Literal["A", "F"], model_kind[-1])
    values: list[float] = []
    failures: list[str] = []
    for row in rows:
        assessment = assess_length_feature(
            row.features,
            feature_name=feature_name,
            annotation_sha256=metadata.annotation_sha256,
            counting_policy_sha256=metadata.counting_policy_sha256,
            applicability=metadata.applicability,
            qc=qc,
            evidence_domain=row.evidence_domain,
        )
        if assessment.reasons:
            failures.extend(assessment.reasons)
        else:
            values.append(cast(float, assessment.feature_value))
    if failures:
        return _failure(candidate_id, model_kind, *dict.fromkeys(failures))
    if len(values) < 2:
        return _failure(candidate_id, model_kind, "rank_deficient_design")
    minimum = min(values)
    maximum = max(values)
    if minimum == maximum:
        return _failure(candidate_id, model_kind, "zero_training_feature_range")

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            design = np.column_stack((np.ones(len(values), dtype=float), np.asarray(values, dtype=float)))
            truth = np.asarray([row.total_truth_repeat_count for row in rows], dtype=float)
            left, singular, right = np.linalg.svd(design, full_matrices=False)
    except FloatingPointError:
        return _failure(candidate_id, model_kind, "nonfinite_design")
    except np.linalg.LinAlgError:
        return _failure(candidate_id, model_kind, "svd_fit_failed")
    tolerance = np.finfo(float).eps * max(design.shape) * singular[0]
    if sum(singular > tolerance) != 2:
        return _failure(candidate_id, model_kind, "rank_deficient_design")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        condition = singular[0] / singular[-1]
    if not math.isfinite(condition) or condition > metadata.maximum_condition_number:
        return _failure(candidate_id, model_kind, "ill_conditioned_design")
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            coefficients = right.T @ ((left.T @ truth) / singular)
    except FloatingPointError:
        return _failure(candidate_id, model_kind, "nonfinite_fit")
    if coefficients.shape != (2,) or not np.all(np.isfinite(coefficients)):
        return _failure(candidate_id, model_kind, "nonfinite_fit")

    expansion = (maximum - minimum) * 0.1
    lower_bound = minimum - expansion
    upper_bound = maximum + expansion
    if not math.isfinite(lower_bound) or not math.isfinite(upper_bound):
        return _failure(candidate_id, model_kind, "nonfinite_feature_bounds")
    raw: dict[str, object] = {
        "schema_version": "length-model-v1",
        "target": {
            "name": "total_diploid_repeat_count",
            "boundary_definition": TARGET_BOUNDARY_DEFINITION,
        },
        "unit": "repeat_units",
        "model_kind": model_kind,
        "feature_order": [feature_name],
        "intercept": float(coefficients[0]),
        "coefficients": [float(coefficients[1])],
        "annotation_sha256": metadata.annotation_sha256,
        "counting_policy_sha256": metadata.counting_policy_sha256,
        "applicability": candidate_applicability_document(metadata.applicability, target="length"),
        "qc": dict(qc),
        "feature_bounds": {feature_name: {"minimum": lower_bound, "maximum": upper_bound}},
        "study_sha256": metadata.study_sha256,
        "training_evidence_sha256": metadata.training_evidence_sha256,
        "producer": candidate_producer_document(metadata.producer),
    }
    return LengthFitOutcome(candidate_id, model_kind, "fitted", (), decode_length_model(raw))


def fit_length_hypotheses(
    rows: Sequence[LengthTrainingRow],
    training_roster: LengthEligibleRoster,
    protocol: LengthProtocol,
    metadata: LengthTrainingMetadata,
) -> LengthTrainingResult:
    """Fit every declared affine candidate using training rows only.

    Args:
        rows: Exact eligible training representatives and matched total truth.
        training_roster: Separately frozen roster for this training role.
        protocol: Outcome-independent finite candidate and QC declaration.
        metadata: Hash-bound training inputs and stability threshold.

    Returns:
        Frozen baseline plus one fitted or ineligible outcome per candidate.

    Raises:
        ValueError: If typed inputs, role isolation, roster binding, or QC identity are invalid.
    """
    roster_document = length_eligible_roster_document(training_roster)
    length_protocol_document(protocol)
    metadata_document = length_training_metadata_document(metadata)
    if training_roster.sha256 != metadata.training_roster_sha256:
        _fail("length training roster differs from metadata")
    if dict(protocol.qc) != dict(metadata.qc):
        _fail("length training protocol and metadata QC differ")
    checked_rows = _training_rows(rows, training_roster)
    qc = _model_qc(metadata)
    outcomes = []
    for hypothesis in protocol.candidates:
        if hypothesis.model_kind in {"physical-A", "physical-F"}:
            outcomes.append(
                _failure(hypothesis.candidate_id, hypothesis.model_kind, "physical_model_geometry_evidence_unsupported")
            )
        else:
            outcomes.append(
                _affine_fit(
                    hypothesis.candidate_id,
                    cast(LengthModelKind, hypothesis.model_kind),
                    checked_rows,
                    metadata,
                    qc,
                )
            )
    return LengthTrainingResult(
        metadata_sha256=canonical_sha256(metadata_document),
        protocol_sha256=protocol.sha256,
        training_roster_sha256=canonical_sha256(roster_document),
        baseline=_baseline(checked_rows, metadata),
        outcomes=tuple(outcomes),
    )
