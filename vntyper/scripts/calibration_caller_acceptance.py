"""Frozen caller safety gates and finite, benefit-requiring policy selection."""

from __future__ import annotations

import logging
import random
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Literal, NoReturn

from vntyper.scripts.calibration_caller_metrics import CallerMetrics, CallerObservation, calculate_caller_metrics
from vntyper.scripts.calibration_caller_roster import CallerEligibleRoster, bind_caller_observations
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)
GateStatus = Literal["passed", "failed", "insufficient-evidence"]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


@dataclass(frozen=True)
class CallerGateRules:
    """Rules frozen by the enclosing study protocol before outcome access."""

    required_strata: tuple[str, ...]
    minimum_positive_groups: int
    minimum_negative_groups: int
    seed: int
    maximum_fpr_upper: Fraction = Fraction(1, 100)
    minimum_sensitivity_delta_lower: Fraction = Fraction(-1, 50)
    maximum_no_call_increase: Fraction = Fraction(1, 50)
    minimum_exact_benefit: Fraction = Fraction(1, 100)
    minimum_fpr_benefit: Fraction = Fraction(1, 200)

    def __post_init__(self) -> None:
        if (
            not isinstance(self.required_strata, tuple)
            or not self.required_strata
            or any(not isinstance(value, str) or not value or value.strip() != value for value in self.required_strata)
            or self.required_strata != tuple(sorted(set(self.required_strata)))
        ):
            _fail("caller gate strata must be sorted unique non-empty immutable strings")
        for value, minimum in ((self.minimum_positive_groups, 1), (self.minimum_negative_groups, 1), (self.seed, 0)):
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                _fail("caller gate group minima and seed must be valid integers")
        for name in (
            "maximum_fpr_upper",
            "minimum_sensitivity_delta_lower",
            "maximum_no_call_increase",
            "minimum_exact_benefit",
            "minimum_fpr_benefit",
        ):
            value = getattr(self, name)
            minimum = -1 if name == "minimum_sensitivity_delta_lower" else 0
            if not isinstance(value, Fraction) or not minimum <= value <= 1:
                _fail("caller gate rates must be exact bounded Fractions")
        if self.minimum_exact_benefit <= 0 or self.minimum_fpr_benefit <= 0:
            _fail("caller selection benefit thresholds must be positive")


@dataclass(frozen=True)
class CallerPopulationGate:
    """Baseline comparison on one complete independent population."""

    status: GateStatus
    reasons: tuple[str, ...]
    candidate: CallerMetrics | None
    baseline: CallerMetrics | None
    sensitivity_interval: tuple[Fraction, Fraction] | None
    no_call_increase: Fraction | None


@dataclass(frozen=True)
class CallerAcceptance:
    """Descriptive gate results; custody and candidate-policy bindings are separate."""

    status: GateStatus
    phase: str
    rules_sha256: str
    eligible_roster_sha256: str
    baseline_predictions_sha256: str
    pooled: CallerPopulationGate
    strata: Mapping[str, CallerPopulationGate]
    selection_benefit: bool


@dataclass(frozen=True)
class CallerSelectionEntry:
    """One protocol-enumerated candidate and its fixed selection evaluation."""

    candidate_id: str
    free_parameters: int
    acceptance: CallerAcceptance


def _rules_digest(rules: CallerGateRules) -> str:
    if not isinstance(rules, CallerGateRules):
        _fail("caller acceptance requires CallerGateRules")
    rules.__post_init__()
    return canonical_sha256(
        {
            "schema_version": "caller-acceptance-rules-v1",
            "required_strata": list(rules.required_strata),
            "minimum_positive_groups": rules.minimum_positive_groups,
            "minimum_negative_groups": rules.minimum_negative_groups,
            "seed": rules.seed,
            "bootstrap": "paired-group-central-95-percentile-10000-v1",
            "rates": {
                name: [getattr(rules, name).numerator, getattr(rules, name).denominator]
                for name in (
                    "maximum_fpr_upper",
                    "minimum_sensitivity_delta_lower",
                    "maximum_no_call_increase",
                    "minimum_exact_benefit",
                    "minimum_fpr_benefit",
                )
            },
        }
    )


def _sensitivity_interval(
    candidate: tuple[CallerObservation, ...], baseline: tuple[CallerObservation, ...], seed: int
) -> tuple[Fraction, Fraction] | None:
    differences = tuple(
        int(row.called_positive is True) - int(old.called_positive is True)
        for row, old in zip(candidate, baseline, strict=True)
        if row.truth_positive is True
    )
    count = len(differences)
    if count < 2:
        return None
    if len(set(differences)) == 1:
        value = Fraction(differences[0])
        return value, value
    generator = random.Random(seed)
    values = sorted(sum(generator.choices(differences, k=count)) for _ in range(10_000))
    # Same prespecified floor-index central percentile rule as calibration_statistics.
    return Fraction(values[249], count), Fraction(values[9749], count)


def _population(
    candidate: tuple[CallerObservation, ...], baseline: tuple[CallerObservation, ...], rules: CallerGateRules
) -> CallerPopulationGate:
    if not candidate:
        return CallerPopulationGate("insufficient-evidence", ("missing_required_stratum",), None, None, None, None)
    current = calculate_caller_metrics(candidate)
    old = calculate_caller_metrics(baseline)
    interval = _sensitivity_interval(candidate, baseline, rules.seed)
    increase = Fraction(current.no_calls - old.no_calls, current.eligible_count)
    reasons: list[str] = []
    sparse = False
    for count, minimum, reason in (
        (current.sensitivity.total, rules.minimum_positive_groups, "insufficient_positive_groups"),
        (current.false_positive_rate.total, rules.minimum_negative_groups, "insufficient_negative_groups"),
    ):
        if count < minimum:
            sparse = True
            reasons.append(reason)
    if current.wrong_tier_a_identity_groups:
        reasons.append("wrong_tier_a_identity")
    if current.fpr_one_sided_upper is None or current.fpr_one_sided_upper > rules.maximum_fpr_upper:
        reasons.append("false_positive_upper_bound")
    if interval is None:
        sparse = True
        reasons.append("undefined_sensitivity_interval")
    elif interval[0] < rules.minimum_sensitivity_delta_lower:
        reasons.append("sensitivity_noninferiority")
    if increase > rules.maximum_no_call_increase:
        reasons.append("no_call_increase")
    status: GateStatus = "insufficient-evidence" if sparse else "failed" if reasons else "passed"
    return CallerPopulationGate(status, tuple(reasons), current, old, interval, increase)


def _benefit(gate: CallerPopulationGate, rules: CallerGateRules) -> bool:
    if gate.candidate is None or gate.baseline is None:
        return False
    exact = gate.candidate.exact_variant_recovery.estimate
    previous_exact = gate.baseline.exact_variant_recovery.estimate
    fpr = gate.candidate.false_positive_rate.estimate
    previous_fpr = gate.baseline.false_positive_rate.estimate
    return (
        exact is not None and previous_exact is not None and exact - previous_exact >= rules.minimum_exact_benefit
    ) or (fpr is not None and previous_fpr is not None and previous_fpr - fpr >= rules.minimum_fpr_benefit)


def evaluate_caller_acceptance(
    candidate: Sequence[CallerObservation],
    baseline: Sequence[CallerObservation],
    roster: CallerEligibleRoster,
    rules: CallerGateRules,
    *,
    phase: str,
) -> CallerAcceptance:
    """Apply every frozen safety gate to pooled and mandatory-stratum outcomes.

    Args:
        candidate: Complete fixed candidate outcomes from production replay/rerun.
        baseline: Complete fixed baseline outcomes on the identical truth population.
        roster: Pre-outcome independent primary representatives and strata.
        rules: Thresholds/counts/seed frozen in the enclosing study protocol.
        phase: Authorized evidence phase; validation results cannot be selected.

    Returns:
        Immutable gate results and descriptive benefit. This function performs no
        fitting, evidence access, custody changes or promotion authorization.

    Raises:
        ValueError: If rules, population completeness or paired truth differ.
    """
    if not isinstance(phase, str) or phase not in {
        "policy-selection",
        "validation",
        "locked-heldout",
        "development-assessment",
    }:
        _fail("unsupported caller acceptance phase")
    rules_sha256 = _rules_digest(rules)
    current = bind_caller_observations(candidate, roster)
    old = bind_caller_observations(baseline, roster)
    if any(
        (row.truth_positive, row.truth_variants) != (previous.truth_positive, previous.truth_variants)
        for row, previous in zip(current, old, strict=True)
    ):
        _fail("caller acceptance candidate and baseline truth differ")
    pooled = _population(current, old, rules)
    strata = {}
    for name in rules.required_strata:
        keys = {member.key for member in roster.members if name in member.strata}
        strata[name] = _population(
            tuple(row for row in current if row.key in keys), tuple(row for row in old if row.key in keys), rules
        )
    statuses = {pooled.status, *(gate.status for gate in strata.values())}
    status: GateStatus = (
        "insufficient-evidence"
        if "insufficient-evidence" in statuses
        else ("failed" if "failed" in statuses else "passed")
    )
    baseline_digest = canonical_sha256(
        [
            {
                "key": row.key,
                "group": row.group_key,
                "truth": row.truth_positive,
                "truth_variants": None if row.truth_variants is None else list(row.truth_variants),
                "called": row.called_positive,
                "variants": list(row.called_variants),
                "tier_a": list(row.tier_a_variants),
            }
            for row in old
        ]
    )
    return CallerAcceptance(
        status,
        phase,
        rules_sha256,
        roster.sha256,
        baseline_digest,
        pooled,
        MappingProxyType(strata),
        _benefit(pooled, rules),
    )


def select_caller_candidate(entries: Sequence[CallerSelectionEntry]) -> str | None:
    """Select one feasible beneficial candidate using the frozen lexicographic objective.

    Args:
        entries: Protocol-enumerated candidates evaluated on policy-selection data
            only. The workflow owns phase authorization and policy-value bindings.

    Returns:
        One candidate identity, or None when no candidate passes with benefit.
        Lower false calls/wrong tier-A counts precede higher exact recovery and
        sensitivity, fewer free parameters and lexicographically smaller identity.

    Raises:
        ValueError: If candidate identities, complexity or comparison bindings differ.
    """
    if not isinstance(entries, (tuple, list)) or not entries:
        _fail("caller selection requires a non-empty finite candidate list")
    ids: set[str] = set()
    bindings: set[tuple[str, str, str]] = set()
    eligible: list[CallerSelectionEntry] = []
    for entry in entries:
        if not isinstance(entry, CallerSelectionEntry) or not isinstance(entry.acceptance, CallerAcceptance):
            _fail("caller selection requires typed candidate evaluations")
        if not isinstance(entry.candidate_id, str) or re.fullmatch(r"[0-9a-f]{64}", entry.candidate_id) is None:
            _fail("caller selection candidate identities must be SHA256 digests")
        if entry.candidate_id in ids:
            _fail("caller selection candidate identities must be unique")
        ids.add(entry.candidate_id)
        if (
            isinstance(entry.free_parameters, bool)
            or not isinstance(entry.free_parameters, int)
            or entry.free_parameters < 0
        ):
            _fail("caller selection free-parameter count must be a non-negative integer")
        result = entry.acceptance
        if result.phase != "policy-selection":
            _fail("caller candidate selection requires policy-selection evidence")
        bindings.add((result.rules_sha256, result.eligible_roster_sha256, result.baseline_predictions_sha256))
        if result.status == "passed" and result.selection_benefit:
            eligible.append(entry)
    if len(bindings) != 1:
        _fail("caller selection candidates must share rules, population and fixed baseline")
    return min(eligible, key=_selection_key).candidate_id if eligible else None


def _selection_key(entry: CallerSelectionEntry) -> tuple[object, ...]:
    metrics = entry.acceptance.pooled.candidate
    if metrics is None:
        _fail("a passing caller evaluation requires population metrics")
    exact = metrics.exact_variant_recovery.estimate
    sensitivity = metrics.sensitivity.estimate
    return (
        metrics.false_positives,
        metrics.wrong_tier_a_identity_groups,
        -(exact if exact is not None else Fraction(-1)),
        -(sensitivity if sensitivity is not None else Fraction(-1)),
        entry.free_parameters,
        entry.candidate_id,
    )
