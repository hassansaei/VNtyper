"""Independent-group caller metrics with explicit no-call and identity denominators."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import NoReturn

from vntyper.scripts.calibration_statistics import clopper_pearson_interval

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CallerObservation:
    """One frozen independent-group representative at a fixed operating point.

    Unknown truth is None. Missing positive variant identity is None, while an
    independently established negative truth has the empty variant tuple.
    A no-call is None and has no called identities. Tier-A identities are a
    subset of called identities, using the existing production tier assignment.
    """

    key: str
    group_key: str
    truth_positive: bool | None
    truth_variants: tuple[str, ...] | None
    called_positive: bool | None
    called_variants: tuple[str, ...]
    tier_a_variants: tuple[str, ...]


@dataclass(frozen=True)
class CallerRate:
    """Observed proportion and central 95% exact-binomial confidence interval."""

    events: int
    total: int
    estimate: Fraction | None
    lower: Fraction | None
    upper: Fraction | None


@dataclass(frozen=True)
class CallerMetrics:
    """Counts and rates; no-calls are neither true negatives nor false negatives."""

    eligible_count: int
    known_truth_count: int
    unknown_truth_count: int
    true_positives: int
    false_positives: int
    true_negatives: int
    false_negatives: int
    positive_no_calls: int
    negative_no_calls: int
    no_calls: int
    sensitivity: CallerRate
    specificity: CallerRate
    false_positive_rate: CallerRate
    precision: CallerRate
    no_call_rate: CallerRate
    assessability: CallerRate
    study_prevalence: Fraction | None
    exact_variant_recovery: CallerRate
    positive_truth_missing_identity: int
    identity_assessable_count: int
    wrong_identity_groups: int
    wrong_tier_a_identity_groups: int
    fpr_one_sided_upper: Fraction | None


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _variants(values: object) -> None:
    if (
        not isinstance(values, tuple)
        or any(not isinstance(value, str) or not value or value.strip() != value for value in values)
        or values != tuple(sorted(set(values)))
    ):
        _fail("caller variant identities must be sorted unique non-empty strings in a tuple")


def _observations(rows: Sequence[CallerObservation]) -> tuple[CallerObservation, ...]:
    if not isinstance(rows, (tuple, list)) or not rows:
        _fail("caller metrics require non-empty typed observations")
    keys: set[str] = set()
    groups: set[str] = set()
    for row in rows:
        if not isinstance(row, CallerObservation):
            _fail("caller metrics require CallerObservation values")
        for value in (row.key, row.group_key):
            if not isinstance(value, str) or not value or value.strip() != value:
                _fail("caller observation identifiers must be non-empty trimmed strings")
        for decision in (row.truth_positive, row.called_positive):
            if decision is not None and not isinstance(decision, bool):
                _fail("caller truth and calls must be boolean or None")
        _variants(row.called_variants)
        _variants(row.tier_a_variants)
        if row.truth_variants is not None:
            _variants(row.truth_variants)
        if row.truth_positive is None and row.truth_variants is not None:
            _fail("unknown caller truth cannot supply variant identity")
        if row.truth_positive is False and row.truth_variants != ():
            _fail("negative caller truth requires an empty identity tuple")
        if row.truth_positive is True and row.truth_variants == ():
            _fail("positive caller truth requires known identities or explicit missing identity")
        if row.called_positive is not True and row.called_variants:
            _fail("negative calls and no-calls cannot carry variant identities")
        if not set(row.tier_a_variants).issubset(row.called_variants):
            _fail("tier-A identities must be a subset of called variant identities")
        if row.key in keys or row.group_key in groups:
            _fail("caller metrics require unique specimen and independent-group representatives")
        keys.add(row.key)
        groups.add(row.group_key)
    return tuple(sorted(rows, key=lambda row: row.group_key))


def _rate(events: int, total: int) -> CallerRate:
    if total == 0:
        return CallerRate(0, 0, None, None, None)
    interval = clopper_pearson_interval(events, total)
    return CallerRate(events, total, interval.estimate, interval.lower, interval.upper)


def one_sided_fpr_upper(events: int, groups: int, *, alpha: Fraction = Fraction(1, 20)) -> Fraction:
    """Compute a one-sided exact upper false-positive proportion bound.

    Args:
        events: Independent negative-truth groups with a positive call.
        groups: All eligible independent negative-truth groups, including no-calls.
        alpha: Exact upper-tail error probability, strictly between zero and 1/2.

    Returns:
        Clopper-Pearson upper bound at confidence 1-alpha. The central interval
        invoked internally uses confidence 1-2*alpha, giving the required tail.

    Raises:
        ValueError: If counts or the tail probability are invalid.
    """
    if not isinstance(alpha, Fraction) or not 0 < alpha < Fraction(1, 2):
        _fail("one-sided FPR alpha must be an exact fraction between zero and one half")
    return clopper_pearson_interval(events, groups, confidence=1 - 2 * alpha).upper


def calculate_caller_metrics(rows: Sequence[CallerObservation]) -> CallerMetrics:
    """Summarize fixed calls without fitting thresholds or inferring truth.

    Args:
        rows: One prespecified representative per independent group. The caller
            must bind this complete roster before extracting predictions.

    Returns:
        Binary and exact-identity counts with explicit denominators. Sensitivity
        and specificity retain no-calls in their truth denominators, so
        specificity plus FPR can be less than one. Precision is conditional on
        known truth and accompanied by the study prevalence. Exact recovery is
        full identity-set agreement on positives whose identity is known.

    Raises:
        ValueError: If observations are empty, inconsistent or duplicated.
    """
    observations = _observations(rows)
    positives = tuple(row for row in observations if row.truth_positive is True)
    negatives = tuple(row for row in observations if row.truth_positive is False)
    known = len(positives) + len(negatives)
    tp = sum(row.called_positive is True for row in positives)
    fn = sum(row.called_positive is False for row in positives)
    tn = sum(row.called_positive is False for row in negatives)
    fp = sum(row.called_positive is True for row in negatives)
    positive_no_calls = sum(row.called_positive is None for row in positives)
    negative_no_calls = sum(row.called_positive is None for row in negatives)
    no_calls = sum(row.called_positive is None for row in observations)
    identity_rows = tuple(row for row in observations if row.truth_variants is not None)
    exact_positives = tuple(row for row in positives if row.truth_variants is not None)
    recovered = sum(
        row.called_positive is True and row.called_variants == row.truth_variants for row in exact_positives
    )
    wrong = sum(bool(set(row.called_variants) - set(row.truth_variants or ())) for row in identity_rows)
    wrong_tier_a = sum(bool(set(row.tier_a_variants) - set(row.truth_variants or ())) for row in identity_rows)
    count = len(observations)
    return CallerMetrics(
        count,
        known,
        count - known,
        tp,
        fp,
        tn,
        fn,
        positive_no_calls,
        negative_no_calls,
        no_calls,
        _rate(tp, len(positives)),
        _rate(tn, len(negatives)),
        _rate(fp, len(negatives)),
        _rate(tp, tp + fp),
        _rate(no_calls, count),
        _rate(count - no_calls, count),
        Fraction(len(positives), known) if known else None,
        _rate(recovered, len(exact_positives)),
        len(positives) - len(exact_positives),
        len(identity_rows),
        wrong,
        wrong_tier_a,
        one_sided_fpr_upper(fp, len(negatives)) if negatives else None,
    )
