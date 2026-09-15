"""Offline reports for hash-bound total-length evaluation evidence."""

from __future__ import annotations

import logging
import re
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import NoReturn

from jinja2 import Environment, FileSystemLoader, StrictUndefined, select_autoescape

import vntyper
from vntyper.scripts.calibration_length import LengthBaseline, length_baseline_document
from vntyper.scripts.calibration_length_acceptance import (
    LengthGateResult,
    decode_length_prediction_context,
    evaluate_length_acceptance,
)
from vntyper.scripts.calibration_length_evaluation import (
    LengthCandidateEvaluation,
    LengthEvaluationResult,
    LengthPredictionEvidence,
    length_evaluation_document,
    select_length_candidates,
)
from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    LengthObservation,
    length_eligible_roster_document,
)
from vntyper.scripts.calibration_length_protocol import LengthHypothesis, LengthProtocol, length_protocol_document

logger = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator} ({100 * float(value):.2f}%)"


def _number(value: float | None) -> str:
    return "undefined" if value is None else f"{value:.6g}"


def _metric_rows(gate: LengthGateResult) -> list[tuple[str, str]] | None:
    metrics = gate.metrics
    if metrics is None:
        return None
    interval = gate.paired_error_difference_interval
    paired = "not applicable"
    if interval is not None:
        paired = f"{interval[0]:.6g} to {interval[1]:.6g} repeat units"
    return [
        ("Eligible independent groups", str(metrics.eligible_count)),
        ("Assessable predictions", str(metrics.assessable_count)),
        ("Mean absolute error", _number(metrics.mae)),
        ("Median absolute error", _number(metrics.median_absolute_error)),
        ("Root mean squared error", _number(metrics.rmse)),
        ("Bias", _number(metrics.bias)),
        ("R²", _number(metrics.r_squared)),
        ("Baseline mean absolute error", _number(metrics.baseline_mae)),
        ("Baseline-relative MAE improvement", _number(metrics.relative_mae_improvement)),
        ("Within tolerance", _fraction(metrics.within_tolerance)),
        ("Within tolerance one-sided 95% lower bound", _fraction(metrics.tolerance_lower)),
        ("Availability", _fraction(metrics.availability)),
        ("Availability one-sided 95% lower bound", _fraction(metrics.availability_lower)),
        ("central 95% paired error-difference interval", paired),
    ]


def _population(name: str, gate: LengthGateResult) -> dict[str, object]:
    return {"name": name, "status": gate.status, "reasons": gate.reasons, "rows": _metric_rows(gate)}


def _scale(value: float, lower: float, upper: float, start: float, end: float) -> float:
    if lower == upper:
        return (start + end) / 2
    return start + (value - lower) * (end - start) / (upper - lower)


def _plots(predictions: tuple[LengthPredictionEvidence, ...]) -> list[dict[str, object]]:
    assessable = tuple(item for item in predictions if item.prediction is not None)
    if not assessable:
        return []
    truths = [item.truth for item in assessable]
    predicted = [item.prediction for item in assessable if item.prediction is not None]
    residuals = [prediction - truth for truth, prediction in zip(truths, predicted, strict=True)]
    value_min = min(*truths, *predicted)
    value_max = max(*truths, *predicted)
    truth_min, truth_max = min(truths), max(truths)
    residual_min, residual_max = min(0, *residuals), max(0, *residuals)
    scatter_points = [
        {
            "x": _scale(truth, value_min, value_max, 45, 365),
            "y": _scale(prediction, value_min, value_max, 345, 25),
            "label": f"Observation {ordinal}: truth {truth:.6g}; prediction {prediction:.6g}",
        }
        for ordinal, (truth, prediction) in enumerate(zip(truths, predicted, strict=True), start=1)
    ]
    residual_points = [
        {
            "x": _scale(truth, truth_min, truth_max, 45, 365),
            "y": _scale(residual, residual_min, residual_max, 345, 25),
            "label": f"Observation {ordinal}: truth {truth:.6g}; residual {residual:.6g}",
        }
        for ordinal, (truth, residual) in enumerate(zip(truths, residuals, strict=True), start=1)
    ]
    zero_y = _scale(0, residual_min, residual_max, 345, 25)
    return [
        {
            "title": "Predicted versus truth",
            "x_label": "Truth (repeat units)",
            "y_label": "Predicted (repeat units)",
            "points": scatter_points,
            "reference": {"x1": 45, "y1": 345, "x2": 365, "y2": 25},
        },
        {
            "title": "Residuals",
            "x_label": "Truth (repeat units)",
            "y_label": "Prediction − truth",
            "points": residual_points,
            "reference": {"x1": 45, "y1": zero_y, "x2": 365, "y2": zero_y},
        },
    ]


def _observations(
    candidate: LengthCandidateEvaluation,
    roster: LengthEligibleRoster,
    baseline: LengthBaseline,
    common_truth: dict[tuple[str, str], float] | None,
) -> tuple[tuple[LengthObservation, ...], dict[tuple[str, str], float]]:
    expected = {(member.key, member.group_key) for member in roster.members}
    observed: set[tuple[str, str]] = set()
    truth: dict[tuple[str, str], float] = {}
    rows: list[LengthObservation] = []
    for prediction in candidate.predictions:
        identity = (prediction.key, prediction.group_key)
        if identity in observed:
            _fail("length report prediction identities must be unique")
        observed.add(identity)
        if prediction.baseline_prediction != baseline.mean_total_repeat_count:
            _fail("length report prediction differs from the frozen baseline prediction")
        if (prediction.prediction is None) != bool(prediction.availability_reasons):
            _fail("length report prediction availability and reasons are inconsistent")
        if not isinstance(prediction.features_sha256, str) or _DIGEST.fullmatch(prediction.features_sha256) is None:
            _fail("length report prediction feature identity must be a lowercase SHA256 digest")
        truth[identity] = prediction.truth
        rows.append(
            LengthObservation(
                prediction.key,
                prediction.group_key,
                prediction.truth,
                prediction.prediction,
                prediction.baseline_prediction,
            )
        )
    if observed != expected:
        _fail("length report prediction set must match the frozen roster exactly")
    if common_truth is not None and truth != common_truth:
        _fail("length report evaluated candidates must share one common truth population")
    return tuple(rows), truth


def _validate_status(candidate: LengthCandidateEvaluation, hypothesis: LengthHypothesis) -> None:
    if (
        candidate.candidate_id != hypothesis.candidate_id
        or candidate.model_kind != hypothesis.model_kind
        or candidate.free_parameters != hypothesis.free_parameters
    ):
        _fail("length report candidate roster differs from the frozen protocol")
    if candidate.status == "evaluated":
        if candidate.reasons or candidate.model_sha256 is None or candidate.acceptance is None:
            _fail("length report evaluated candidate state is inconsistent")
        if _DIGEST.fullmatch(candidate.model_sha256) is None:
            _fail("length report evaluated model identity must be a lowercase SHA256 digest")
    elif candidate.status == "fit-ineligible":
        if (
            not candidate.reasons
            or candidate.model_sha256 is not None
            or candidate.predictions
            or candidate.acceptance is not None
        ):
            _fail("length report fit-ineligible candidate state is inconsistent")
    elif candidate.status == "not-evaluated":
        if (
            not candidate.reasons
            or candidate.model_sha256 is None
            or candidate.predictions
            or candidate.acceptance is not None
        ):
            _fail("length report non-evaluated candidate state is inconsistent")
    else:
        _fail("length report candidate status is invalid")


def _candidate_entry(candidate: LengthCandidateEvaluation, protocol: LengthProtocol) -> dict[str, object]:
    acceptance = candidate.acceptance
    unavailable = Counter(reason for prediction in candidate.predictions for reason in prediction.availability_reasons)
    return {
        "id": candidate.candidate_id,
        "kind": candidate.model_kind,
        "status": candidate.status,
        "reasons": candidate.reasons,
        "unavailable": sorted(unavailable.items()),
        "populations": (
            []
            if acceptance is None
            else [
                _population("Pooled", acceptance.pooled),
                *(
                    _population(f"Mandatory stratum: {name}", acceptance.strata[name])
                    for name in protocol.required_strata
                ),
            ]
        ),
        "plots": _plots(candidate.predictions) if candidate.status == "evaluated" else [],
    }


def render_length_evaluation(
    result: LengthEvaluationResult,
    roster: LengthEligibleRoster,
    protocol: LengthProtocol,
    baseline: LengthBaseline,
) -> str:
    """Render anonymous offline HTML after recomputing metrics and selection.

    Args:
        result: Frozen local evidence produced by length model evaluation.
        roster: Exact role-specific eligible population bound by the protocol.
        protocol: Frozen candidate, uncertainty, strata, QC and gate rules.
        baseline: Frozen training-mean comparator used for every prediction.

    Returns:
        Standalone escaped HTML containing descriptive metrics and static SVGs.
        Rendering neither revalidates upstream feature/model provenance nor
        authorizes evidence exposure or model promotion.

    Raises:
        ValueError: If canonical content, bindings, populations, metrics, phase,
            selection or immutable candidate states are inconsistent.
    """
    length_evaluation_document(result)
    length_eligible_roster_document(roster)
    length_protocol_document(protocol)
    try:
        length_baseline_document(baseline)
    except ValueError:
        _fail("length report baseline binding is invalid")
    if (
        result.eligible_roster_sha256 != roster.sha256
        or result.protocol_sha256 != protocol.sha256
        or result.baseline_sha256 != baseline.sha256
        or roster.sha256 != protocol.eligible_roster_sha256
    ):
        _fail("length report artifacts differ from their exact frozen bindings")
    context = decode_length_prediction_context(
        {
            "schema_version": "length-prediction-context-v1",
            "protocol_sha256": protocol.sha256,
            "qc_sha256": protocol.qc_sha256,
        }
    )
    if context.sha256 != result.prediction_context_sha256:
        _fail("length report prediction context differs from the protocol binding")
    if len(result.candidates) != len(protocol.candidates):
        _fail("length report candidate roster differs from the frozen protocol")
    evaluated = 0
    common_truth: dict[tuple[str, str], float] | None = None
    for candidate, hypothesis in zip(result.candidates, protocol.candidates, strict=True):
        _validate_status(candidate, hypothesis)
        if candidate.status != "evaluated":
            continue
        evaluated += 1
        observations, truth = _observations(candidate, roster, baseline, common_truth)
        common_truth = truth
        rebuilt = evaluate_length_acceptance(observations, roster, context, protocol)
        if rebuilt != candidate.acceptance:
            _fail("length report stored acceptance differs from recomputed evidence")
    if result.phase == "policy-selection":
        if any(candidate.status == "not-evaluated" for candidate in result.candidates):
            _fail("length report policy selection cannot omit a fitted candidate")
    elif evaluated != 1:
        _fail("length report nonselection evidence must contain exactly one fixed evaluated candidate")
    if select_length_candidates(result.candidates, result.phase) != result.selection:
        _fail("length report stored selection differs from recomputed evidence")
    entries = [_candidate_entry(candidate, protocol) for candidate in result.candidates]
    environment = Environment(
        loader=FileSystemLoader(str(Path(vntyper.__file__).resolve().parent / "templates")),
        autoescape=select_autoescape(["html"]),
        undefined=StrictUndefined,
    )
    return environment.get_template("calibration_length_report.html").render(
        phase=result.phase,
        protocol=protocol.sha256,
        roster=roster.sha256,
        baseline=baseline.sha256,
        evidence=result.sha256,
        selection=result.selection,
        entries=entries,
    )
