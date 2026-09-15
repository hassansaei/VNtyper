"""Pure evaluation and policy selection for frozen total-length models."""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
from types import MappingProxyType
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_length import (
    LengthBaseline,
    LengthFitOutcome,
    length_baseline_document,
)
from vntyper.scripts.calibration_length_acceptance import (
    LengthAcceptance,
    LengthGateResult,
    LengthPredictionContext,
    decode_length_prediction_context,
    evaluate_length_acceptance,
)
from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    LengthMetrics,
    LengthObservation,
    length_eligible_roster_document,
)
from vntyper.scripts.calibration_length_protocol import LengthHypothesis, LengthProtocol, length_protocol_document
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_estimation import EvidenceDomain, estimate_total_repeats
from vntyper.scripts.length_features import LengthFeatures, encode_length_features
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION, LengthModel, encode_length_model

logger = logging.getLogger(__name__)

EvaluationPhase = Literal["policy-selection", "validation", "locked-heldout", "development-assessment"]
CandidateEvaluationStatus = Literal["evaluated", "fit-ineligible", "not-evaluated"]
SelectionStatus = Literal["selected", "no-feasible-candidate", "not-applicable"]
_PHASES = {"policy-selection", "validation", "locked-heldout", "development-assessment"}
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_QC_FIELDS = (
    "minimum_denominator_mean_depth",
    "minimum_denominator_covered_fraction",
    "minimum_denominator_supporting_fragments",
)
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))


@dataclass(frozen=True)
class LengthEvaluationRow:
    """One exact eligible evaluation representative and its bound features."""

    key: str
    group_key: str
    phase: EvaluationPhase
    features: LengthFeatures
    truth_boundary_definition: str
    total_truth_repeat_count: float
    evidence_domain: EvidenceDomain


@dataclass(frozen=True)
class LengthPredictionEvidence:
    """One model-produced prediction or explicit unavailable result."""

    key: str
    group_key: str
    truth: float
    prediction: float | None
    baseline_prediction: float
    availability_reasons: tuple[str, ...]
    features_sha256: str


@dataclass(frozen=True)
class LengthCandidateEvaluation:
    """Frozen evidence for one candidate declared by the protocol."""

    candidate_id: str
    model_kind: str
    free_parameters: int
    status: CandidateEvaluationStatus
    reasons: tuple[str, ...]
    model_sha256: str | None
    predictions: tuple[LengthPredictionEvidence, ...]
    acceptance: LengthAcceptance | None


@dataclass(frozen=True)
class LengthSelection:
    """Policy-selection decision; it never grants promotion authorization."""

    status: SelectionStatus
    selected_candidate_id: str | None
    selected_model_sha256: str | None
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class LengthEvaluationResult:
    """Hash-bound local evaluation evidence over exact inputs and models."""

    phase: EvaluationPhase
    protocol_sha256: str
    eligible_roster_sha256: str
    baseline_sha256: str
    prediction_context_sha256: str
    evaluation_rows_sha256: str
    candidates: tuple[LengthCandidateEvaluation, ...]
    selection: LengthSelection
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        _fail(f"length evaluation {label} must be non-empty text")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"length evaluation {label} must be a lowercase SHA256 digest")
    return value


def _exact_truth(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail("length evaluation truth must be a positive integral exact repeat count")
    try:
        result = float(value)
    except OverflowError:
        _fail("length evaluation truth must be a positive integral exact repeat count")
    if not math.isfinite(result) or result <= 0 or not result.is_integer():
        _fail("length evaluation truth must be a positive integral exact repeat count")
    return result


def _evaluation_rows(
    rows: Sequence[LengthEvaluationRow], roster: LengthEligibleRoster
) -> tuple[tuple[LengthEvaluationRow, ...], str, EvaluationPhase]:
    if not isinstance(rows, (tuple, list)) or not rows:
        _fail("length evaluation rows must be a non-empty typed sequence")
    checked: list[LengthEvaluationRow] = []
    keys: set[str] = set()
    groups: set[str] = set()
    phases: set[str] = set()
    for row in rows:
        if not isinstance(row, LengthEvaluationRow):
            _fail("length evaluation rows must contain LengthEvaluationRow values")
        _text(row.key, "row key")
        _text(row.group_key, "row group key")
        if row.key in keys or row.group_key in groups:
            _fail("length evaluation rows must have unique keys and groups")
        keys.add(row.key)
        groups.add(row.group_key)
        if not isinstance(row.phase, str) or row.phase not in _PHASES:
            _fail("length evaluation phase is not allowed")
        phases.add(row.phase)
        if row.truth_boundary_definition != TARGET_BOUNDARY_DEFINITION:
            _fail("length evaluation truth boundary definition is incompatible")
        _exact_truth(row.total_truth_repeat_count)
        if not isinstance(row.evidence_domain, str) or row.evidence_domain not in {"synthetic", "external"}:
            _fail("length evaluation evidence domain must be synthetic or external")
        encode_length_features(row.features)
        if row.features.manifest_key != row.key:
            _fail("length evaluation feature manifest key differs from its row key")
        checked.append(row)
    if len(phases) != 1:
        _fail("length evaluation rows must use one common phase")
    expected = {(member.key, member.group_key) for member in roster.members}
    observed = {(row.key, row.group_key) for row in checked}
    if observed != expected:
        _fail("length evaluation rows must match the frozen roster exactly")
    ordered = tuple(sorted(checked, key=lambda row: row.group_key))
    document = [
        {
            "key": row.key,
            "group_key": row.group_key,
            "phase": row.phase,
            "features_sha256": row.features.sha256,
            "truth_boundary_definition": row.truth_boundary_definition,
            "total_truth_repeat_count": row.total_truth_repeat_count,
            "evidence_domain": row.evidence_domain,
        }
        for row in ordered
    ]
    return ordered, canonical_sha256(document), cast(EvaluationPhase, next(iter(phases)))


def _fit_outcomes(
    outcomes: Sequence[LengthFitOutcome], protocol: LengthProtocol, baseline: LengthBaseline
) -> tuple[tuple[LengthFitOutcome, LengthHypothesis], ...]:
    if not isinstance(outcomes, (tuple, list)) or len(outcomes) != len(protocol.candidates):
        _fail("length fit outcomes must exactly match the protocol candidate roster")
    paired: list[tuple[LengthFitOutcome, LengthHypothesis]] = []
    for outcome, hypothesis in zip(outcomes, protocol.candidates, strict=True):
        if not isinstance(outcome, LengthFitOutcome):
            _fail("length fit outcomes must contain LengthFitOutcome values")
        if outcome.candidate_id != hypothesis.candidate_id or outcome.model_kind != hypothesis.model_kind:
            _fail("length fit outcomes must exactly match protocol candidate IDs and kinds")
        if not isinstance(outcome.reasons, tuple) or any(
            not isinstance(reason, str) or not reason or reason != reason.strip() for reason in outcome.reasons
        ):
            _fail("length fit outcome reasons must be immutable non-empty reason codes")
        if outcome.status == "fitted":
            if outcome.model is None or outcome.reasons:
                _fail("length fitted outcome requires one model and empty reasons")
            _model_binding(outcome.model, hypothesis, protocol, baseline)
        elif outcome.status == "ineligible":
            if outcome.model is not None or not outcome.reasons:
                _fail("length ineligible outcome requires no model and nonempty reasons")
        else:
            _fail("length fit outcome status is invalid")
        paired.append((outcome, hypothesis))
    return tuple(paired)


def _model_binding(
    model: LengthModel, hypothesis: LengthHypothesis, protocol: LengthProtocol, baseline: LengthBaseline
) -> None:
    encode_length_model(model)
    if model.model_kind != hypothesis.model_kind:
        _fail("length fitted model kind differs from its protocol hypothesis")
    if (
        model.study_sha256 != baseline.study_sha256
        or model.training_evidence_sha256 != baseline.training_evidence_sha256
    ):
        _fail("length fitted model study or training evidence differs from the baseline")
    if any(model.qc[name] != protocol.qc[name] for name in _QC_FIELDS):
        _fail("length fitted model QC differs from the protocol QC")


def _prediction_context(protocol: LengthProtocol) -> LengthPredictionContext:
    return decode_length_prediction_context(
        {
            "schema_version": "length-prediction-context-v1",
            "protocol_sha256": protocol.sha256,
            "qc_sha256": protocol.qc_sha256,
        }
    )


def _evaluate(
    outcome: LengthFitOutcome,
    hypothesis: LengthHypothesis,
    rows: tuple[LengthEvaluationRow, ...],
    baseline: LengthBaseline,
    roster: LengthEligibleRoster,
    context: LengthPredictionContext,
    protocol: LengthProtocol,
) -> LengthCandidateEvaluation:
    model = cast(LengthModel, outcome.model)
    predictions: list[LengthPredictionEvidence] = []
    observations: list[LengthObservation] = []
    for row in rows:
        estimate = estimate_total_repeats(row.features, model, evidence_domain=row.evidence_domain)
        prediction = estimate.estimated_total_repeat_count
        predictions.append(
            LengthPredictionEvidence(
                row.key,
                row.group_key,
                row.total_truth_repeat_count,
                prediction,
                baseline.mean_total_repeat_count,
                estimate.reasons,
                row.features.sha256,
            )
        )
        observations.append(
            LengthObservation(
                row.key,
                row.group_key,
                row.total_truth_repeat_count,
                prediction,
                baseline.mean_total_repeat_count,
            )
        )
    acceptance = evaluate_length_acceptance(observations, roster, context, protocol)
    return LengthCandidateEvaluation(
        hypothesis.candidate_id,
        hypothesis.model_kind,
        hypothesis.free_parameters,
        "evaluated",
        (),
        model.sha256,
        tuple(predictions),
        acceptance,
    )


def _unevaluated(outcome: LengthFitOutcome, hypothesis: LengthHypothesis) -> LengthCandidateEvaluation:
    if outcome.status == "ineligible":
        return LengthCandidateEvaluation(
            hypothesis.candidate_id,
            hypothesis.model_kind,
            hypothesis.free_parameters,
            "fit-ineligible",
            outcome.reasons,
            None,
            (),
            None,
        )
    return LengthCandidateEvaluation(
        hypothesis.candidate_id,
        hypothesis.model_kind,
        hypothesis.free_parameters,
        "not-evaluated",
        ("not_fixed_candidate",),
        cast(LengthModel, outcome.model).sha256,
        (),
        None,
    )


def _selection(candidates: tuple[LengthCandidateEvaluation, ...], phase: EvaluationPhase) -> LengthSelection:
    if phase != "policy-selection":
        return LengthSelection("not-applicable", None, None, ("selection_not_allowed_for_phase",))
    scored: list[tuple[LengthCandidateEvaluation, float]] = []
    for candidate in candidates:
        acceptance = candidate.acceptance
        if acceptance is None or acceptance.status != "passed":
            continue
        metrics = acceptance.pooled.metrics
        if metrics is not None and metrics.mae is not None:
            scored.append((candidate, metrics.mae))
    if not scored:
        return LengthSelection("no-feasible-candidate", None, None, ("no_candidate_passed_acceptance",))
    best_mae = min(mae for _, mae in scored)
    shortlist = [candidate for candidate, mae in scored if mae <= best_mae + 1]
    selected = min(shortlist, key=lambda item: (item.free_parameters, item.candidate_id))
    return LengthSelection("selected", selected.candidate_id, selected.model_sha256, ())


def _fraction(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def _metrics_document(metrics: LengthMetrics | None) -> dict[str, object] | None:
    if metrics is None:
        return None
    return {
        "eligible_count": metrics.eligible_count,
        "assessable_count": metrics.assessable_count,
        "mae": metrics.mae,
        "median_absolute_error": metrics.median_absolute_error,
        "rmse": metrics.rmse,
        "bias": metrics.bias,
        "r_squared": metrics.r_squared,
        "baseline_mae": metrics.baseline_mae,
        "relative_mae_improvement": metrics.relative_mae_improvement,
        "within_tolerance": _fraction(metrics.within_tolerance),
        "availability": _fraction(metrics.availability),
        "tolerance_lower": _fraction(metrics.tolerance_lower),
        "availability_lower": _fraction(metrics.availability_lower),
    }


def _gate_document(gate: LengthGateResult) -> dict[str, object]:
    return {
        "status": gate.status,
        "reasons": list(gate.reasons),
        "metrics": _metrics_document(gate.metrics),
        "paired_error_difference_interval": (
            None if gate.paired_error_difference_interval is None else list(gate.paired_error_difference_interval)
        ),
    }


def _acceptance_document(acceptance: LengthAcceptance | None) -> dict[str, object] | None:
    if acceptance is None:
        return None
    return {
        "status": acceptance.status,
        "protocol_sha256": acceptance.protocol_sha256,
        "eligible_roster_sha256": acceptance.eligible_roster_sha256,
        "prediction_context_sha256": acceptance.prediction_context_sha256,
        "pooled": _gate_document(acceptance.pooled),
        "strata": {name: _gate_document(gate) for name, gate in sorted(acceptance.strata.items())},
    }


def _candidate_document(candidate: LengthCandidateEvaluation) -> dict[str, object]:
    return {
        "candidate_id": candidate.candidate_id,
        "model_kind": candidate.model_kind,
        "free_parameters": candidate.free_parameters,
        "status": candidate.status,
        "reasons": list(candidate.reasons),
        "model_sha256": candidate.model_sha256,
        "predictions": [
            {
                "key": prediction.key,
                "group_key": prediction.group_key,
                "truth": prediction.truth,
                "prediction": prediction.prediction,
                "baseline_prediction": prediction.baseline_prediction,
                "availability_reasons": list(prediction.availability_reasons),
                "features_sha256": prediction.features_sha256,
            }
            for prediction in candidate.predictions
        ],
        "acceptance": _acceptance_document(candidate.acceptance),
    }


def _result_payload(result: LengthEvaluationResult) -> dict[str, object]:
    return {
        "schema_version": "length-evaluation-v1",
        "phase": result.phase,
        "protocol_sha256": result.protocol_sha256,
        "eligible_roster_sha256": result.eligible_roster_sha256,
        "baseline_sha256": result.baseline_sha256,
        "prediction_context_sha256": result.prediction_context_sha256,
        "evaluation_rows_sha256": result.evaluation_rows_sha256,
        "candidates": [_candidate_document(candidate) for candidate in result.candidates],
        "selection": {
            "status": result.selection.status,
            "selected_candidate_id": result.selection.selected_candidate_id,
            "selected_model_sha256": result.selection.selected_model_sha256,
            "reasons": list(result.selection.reasons),
        },
    }


def length_evaluation_document(result: LengthEvaluationResult) -> dict[str, object]:
    """Project local evaluation evidence after validating immutable content.

    This local artifact records scientific evaluation and selection evidence; it
    is not a portable model, custody authorization, or promotion decision.

    Args:
        result: Result returned by :func:`evaluate_length_hypotheses`.

    Returns:
        Fresh JSON-compatible content including its canonical payload digest.

    Raises:
        ValueError: If typed collections, identities, or digest were forged.
    """
    if not isinstance(result, LengthEvaluationResult) or not isinstance(result.candidates, tuple):
        _fail("length evaluation result must use the immutable typed contract")
    if not isinstance(result.selection, LengthSelection) or not isinstance(result.selection.reasons, tuple):
        _fail("length evaluation selection must use the immutable typed contract")
    for digest_name in (
        "protocol_sha256",
        "eligible_roster_sha256",
        "baseline_sha256",
        "prediction_context_sha256",
        "evaluation_rows_sha256",
        "sha256",
    ):
        _digest(getattr(result, digest_name), digest_name)
    for candidate in result.candidates:
        if (
            not isinstance(candidate, LengthCandidateEvaluation)
            or not isinstance(candidate.predictions, tuple)
            or not isinstance(candidate.reasons, tuple)
            or any(
                not isinstance(prediction, LengthPredictionEvidence)
                or not isinstance(prediction.availability_reasons, tuple)
                for prediction in candidate.predictions
            )
        ):
            _fail("length evaluation candidates must use immutable typed content")
        acceptance = candidate.acceptance
        if acceptance is not None:
            if not isinstance(acceptance, LengthAcceptance) or not isinstance(acceptance.strata, _MAPPING_PROXY_TYPE):
                _fail("length evaluation acceptance must use immutable typed content")
            gates = (acceptance.pooled, *acceptance.strata.values())
            if any(
                not isinstance(gate, LengthGateResult)
                or not isinstance(gate.reasons, tuple)
                or (gate.metrics is not None and not isinstance(gate.metrics, LengthMetrics))
                or (
                    gate.paired_error_difference_interval is not None
                    and not isinstance(gate.paired_error_difference_interval, tuple)
                )
                for gate in gates
            ):
                _fail("length evaluation acceptance must use immutable typed content")
    payload = _result_payload(result)
    if canonical_sha256(payload) != result.sha256:
        _fail("length evaluation digest differs from its canonical content")
    return {**payload, "sha256": result.sha256}


def evaluate_length_hypotheses(
    fit_outcomes: Sequence[LengthFitOutcome],
    baseline: LengthBaseline,
    roster: LengthEligibleRoster,
    protocol: LengthProtocol,
    rows: Sequence[LengthEvaluationRow],
    *,
    fixed_candidate_id: str | None = None,
) -> LengthEvaluationResult:
    """Evaluate frozen fitted models and select only on policy-selection evidence.

    Args:
        fit_outcomes: Full ordered fit outcome roster declared by the protocol.
        baseline: Frozen training-mean baseline used for every eligible member.
        roster: Exact role-specific eligible roster already frozen in the protocol.
        protocol: Outcome-independent finite candidate and acceptance declaration.
        rows: Exact evaluation features and matched total truth for one phase.
        fixed_candidate_id: Required outside policy selection; identifies the one
            model already authorized by the enclosing study/custody adapter.

    Returns:
        Immutable local predictions, acceptance evidence, and optional selection.
        No model is refitted and no phase label grants custody or promotion.

    Raises:
        ValueError: If any typed artifact, exact roster, phase, or model binding fails.
    """
    length_protocol_document(protocol)
    length_eligible_roster_document(roster)
    length_baseline_document(baseline)
    if roster.sha256 != protocol.eligible_roster_sha256:
        _fail("length evaluation roster differs from the frozen protocol roster")
    checked_rows, rows_sha256, phase = _evaluation_rows(rows, roster)
    paired = _fit_outcomes(fit_outcomes, protocol, baseline)
    if phase == "policy-selection":
        if fixed_candidate_id is not None:
            _fail("length policy-selection must not receive fixed_candidate_id")
    else:
        if not isinstance(fixed_candidate_id, str) or fixed_candidate_id not in {
            hypothesis.candidate_id for _, hypothesis in paired
        }:
            _fail("length nonselection evaluation requires a valid fixed_candidate_id")
        fixed = next(outcome for outcome, hypothesis in paired if hypothesis.candidate_id == fixed_candidate_id)
        if fixed.status != "fitted":
            _fail("length fixed_candidate_id must identify a fitted model")
    context = _prediction_context(protocol)
    candidates = tuple(
        _evaluate(outcome, hypothesis, checked_rows, baseline, roster, context, protocol)
        if outcome.status == "fitted" and (phase == "policy-selection" or hypothesis.candidate_id == fixed_candidate_id)
        else _unevaluated(outcome, hypothesis)
        for outcome, hypothesis in paired
    )
    result = LengthEvaluationResult(
        phase,
        protocol.sha256,
        roster.sha256,
        baseline.sha256,
        context.sha256,
        rows_sha256,
        candidates,
        _selection(candidates, phase),
        "",
    )
    return replace(result, sha256=canonical_sha256(_result_payload(result)))
