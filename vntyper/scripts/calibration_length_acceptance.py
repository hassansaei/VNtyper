"""Fixed total-length acceptance gates on predeclared independent observations."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Literal, NoReturn

from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    LengthMetrics,
    LengthObservation,
    bind_length_observations,
    calculate_length_metrics,
    length_eligible_roster_document,
    paired_length_error_interval,
)
from vntyper.scripts.calibration_length_protocol import LengthProtocol, length_protocol_document
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

GateStatus = Literal["passed", "failed", "insufficient-evidence"]
_CONTEXT_FIELDS = {"schema_version", "protocol_sha256", "qc_sha256"}
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class LengthPredictionContext:
    """Protocol/QC identity binding; it does not attest that pending QC was applied."""

    protocol_sha256: str
    qc_sha256: str
    sha256: str


@dataclass(frozen=True)
class LengthGateResult:
    """Metrics and all failed criteria for one frozen population."""

    status: GateStatus
    reasons: tuple[str, ...]
    metrics: LengthMetrics | None
    paired_error_difference_interval: tuple[float, float] | None


@dataclass(frozen=True)
class LengthAcceptance:
    """Pooled improvement and mandatory stratum accuracy/availability verdicts."""

    status: GateStatus
    protocol_sha256: str
    eligible_roster_sha256: str
    prediction_context_sha256: str
    pooled: LengthGateResult
    strata: Mapping[str, LengthGateResult]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, field: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"length prediction context {field} must be a lowercase SHA256 digest")
    return value


def decode_length_prediction_context(value: object) -> LengthPredictionContext:
    """Decode the protocol and declared-QC identity attached to predictions.

    This is an integrity binding only. Model integration must later prove that
    the declared QC was actually applied before scientific completion is claimed.

    Args:
        value: Closed JSON prediction-context object.

    Returns:
        Immutable context and canonical content digest.

    Raises:
        ValueError: If fields, schema or digests are invalid.
    """
    if not isinstance(value, Mapping) or set(value) != _CONTEXT_FIELDS:
        _fail("length prediction context fields differ from the closed contract")
    if value["schema_version"] != "length-prediction-context-v1":
        _fail("length prediction context schema_version must be length-prediction-context-v1")
    return LengthPredictionContext(
        _digest(value["protocol_sha256"], "protocol_sha256"),
        _digest(value["qc_sha256"], "qc_sha256"),
        canonical_sha256(value),
    )


def _prediction_context_document(context: LengthPredictionContext) -> dict[str, object]:
    return {
        "schema_version": "length-prediction-context-v1",
        "protocol_sha256": context.protocol_sha256,
        "qc_sha256": context.qc_sha256,
    }


def length_prediction_context_document(context: LengthPredictionContext) -> dict[str, object]:
    """Project a prediction context after canonical integrity validation.

    Args:
        context: Decoded immutable prediction context.

    Returns:
        Fresh JSON-compatible context object.

    Raises:
        ValueError: If typed context content or its digest was forged.
    """
    if not isinstance(context, LengthPredictionContext):
        _fail("length prediction context must be a LengthPredictionContext")
    raw = _prediction_context_document(context)
    if decode_length_prediction_context(raw) != context:
        _fail("length prediction context differs from its canonical content or digest")
    return raw


def _population_gate(
    rows: Sequence[LengthObservation], protocol: LengthProtocol, *, require_improvement: bool
) -> LengthGateResult:
    limits = protocol.acceptance
    metrics = calculate_length_metrics(
        rows, tolerance_absolute=limits["tolerance_absolute"], tolerance_relative=limits["tolerance_relative"]
    )
    reasons: list[str] = []
    sparse = metrics.eligible_count < limits["minimum_independent_count"]
    if sparse:
        reasons.append("insufficient_independent_groups")
    if metrics.mae is None:
        reasons.append("no_assessable_predictions")
    elif metrics.mae > limits["maximum_mae"]:
        reasons.append("maximum_mae")
    if metrics.tolerance_lower < Fraction(str(limits["minimum_tolerance_lower_bound"])):
        reasons.append("tolerance_lower_bound")
    if metrics.availability_lower < Fraction(str(limits["minimum_availability_lower_bound"])):
        reasons.append("availability_lower_bound")
    interval = None
    if require_improvement:
        if metrics.relative_mae_improvement is None:
            reasons.append("undefined_relative_mae_improvement")
        elif metrics.relative_mae_improvement < limits["minimum_relative_mae_improvement"]:
            reasons.append("relative_mae_improvement")
        interval = paired_length_error_interval(rows, seed=protocol.seed)
        if interval is None:
            reasons.append("undefined_paired_error_interval")
        elif interval[1] >= limits["paired_error_difference_upper_limit"]:
            reasons.append("paired_error_improvement")
    status: GateStatus = "insufficient-evidence" if sparse else "failed" if reasons else "passed"
    return LengthGateResult(status, tuple(reasons), metrics, interval)


def evaluate_length_acceptance(
    rows: Sequence[LengthObservation],
    roster: LengthEligibleRoster,
    prediction_context: LengthPredictionContext,
    protocol: LengthProtocol,
) -> LengthAcceptance:
    """Evaluate frozen predictions without fitting, selection, or promotion writes.

    Args:
        rows: One primary observation per predeclared independent group. Its
            stratum memberships are frozen separately in the roster.
        roster: Exact pre-outcome eligible members and stratum memberships.
        prediction_context: Integrity binding to protocol and declared QC policy.
        protocol: Validated study rules containing every mandatory stratum.

    Returns:
        All pooled and stratum verdicts. Absolute-error and proportion gates
        apply to every stratum; paired improvement applies to the pooled fixed
        baseline comparison. Sparse metadata populations stay distinct from
        model-induced unavailability. This verdict and prediction-context
        binding do not prove that pending model-stage QC was applied and are not
        an attestation of evidence custody, candidate identity, or applicability.

    Raises:
        ValueError: If observations, stratum labels, or protocol integrity fail.
    """
    length_protocol_document(protocol)
    length_eligible_roster_document(roster)
    length_prediction_context_document(prediction_context)
    if roster.sha256 != protocol.eligible_roster_sha256:
        _fail("length eligible roster digest differs from the frozen protocol")
    if prediction_context.protocol_sha256 != protocol.sha256 or prediction_context.qc_sha256 != protocol.qc_sha256:
        _fail("length prediction context differs from the protocol or QC binding")
    observations = bind_length_observations(rows, roster)
    allowed = set(protocol.required_strata)
    if any(not set(member.strata) <= allowed for member in roster.members):
        _fail("length eligible roster contains an undeclared stratum")
    pooled = _population_gate(observations, protocol, require_improvement=True)
    by_group = {row.group_key: row for row in observations}
    grouped = {
        name: [by_group[member.group_key] for member in roster.members if name in member.strata]
        for name in protocol.required_strata
    }
    strata = {
        name: _population_gate(grouped[name], protocol, require_improvement=False)
        if grouped[name]
        else LengthGateResult("insufficient-evidence", ("missing_required_stratum",), None, None)
        for name in protocol.required_strata
    }
    statuses = {pooled.status, *(item.status for item in strata.values())}
    status: GateStatus = (
        "insufficient-evidence"
        if "insufficient-evidence" in statuses
        else "failed"
        if "failed" in statuses
        else "passed"
    )
    return LengthAcceptance(
        status,
        protocol.sha256,
        roster.sha256,
        prediction_context.sha256,
        pooled,
        MappingProxyType(strata),
    )
