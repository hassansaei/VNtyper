"""Fixed total-length acceptance gates on predeclared independent observations."""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Literal

from vntyper.scripts.calibration_length_metrics import (
    LengthMetrics,
    LengthObservation,
    calculate_length_metrics,
    paired_length_error_interval,
)
from vntyper.scripts.calibration_length_protocol import LengthProtocol, length_protocol_document

logger = logging.getLogger(__name__)

GateStatus = Literal["passed", "failed", "insufficient-evidence"]


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
    pooled: LengthGateResult
    strata: Mapping[str, LengthGateResult]


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


def evaluate_length_acceptance(rows: Sequence[LengthObservation], protocol: LengthProtocol) -> LengthAcceptance:
    """Evaluate frozen predictions without fitting, selection, or promotion writes.

    Args:
        rows: One primary observation per predeclared independent group. Its
            stratum is assigned before observing prediction availability/errors.
        protocol: Validated study rules containing every mandatory stratum.

    Returns:
        All pooled and stratum verdicts. Absolute-error and proportion gates
        apply to every stratum; paired improvement applies to the pooled fixed
        baseline comparison. Sparse metadata populations stay distinct from
        model-induced unavailability. This verdict alone is not an attestation
        of evidence custody, candidate identity, or external applicability.

    Raises:
        ValueError: If observations, stratum labels, or protocol integrity fail.
    """
    length_protocol_document(protocol)
    pooled = _population_gate(rows, protocol, require_improvement=True)
    grouped: dict[str, list[LengthObservation]] = {}
    for row in rows:
        if row.stratum not in protocol.required_strata:
            message = "length acceptance observations contain an undeclared stratum"
            logger.error(message)
            raise ValueError(message)
        grouped.setdefault(row.stratum, []).append(row)
    strata = {
        name: _population_gate(grouped[name], protocol, require_improvement=False)
        if name in grouped
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
    return LengthAcceptance(status, protocol.sha256, pooled, MappingProxyType(strata))
