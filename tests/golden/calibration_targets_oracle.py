"""Independent arithmetic for small invented calibration cohorts.

Only standard-library arithmetic is used; no production schema, decisions, fitting,
statistics, or fixture imports are allowed. This is not a read-level caller oracle.
Binomial inversion deliberately caps cohort size because it uses exact rational
polynomials, not the production numerical beta implementation.
"""

from __future__ import annotations

import math
import operator
from collections import Counter
from collections.abc import Sequence
from fractions import Fraction

Number = int | float | Fraction
Call = tuple[bool | None, bool | None, tuple[str, ...] | None, tuple[str, ...], tuple[str, ...]]
Summary = dict[str, int | Fraction | None]


def affine_fit(pairs: Sequence[tuple[Number, Number]]) -> tuple[Fraction, Fraction, Fraction]:
    """Solve two normal equations exactly; return intercept, slope, training mean.

    Args:
        pairs: Finite invented feature/truth pairs.

    Returns:
        Exact rational ordinary-least-squares coefficients and baseline mean.

    Raises:
        ValueError: If the design has no identifiable slope.
    """
    n = len(pairs)
    if n < 2:
        raise ValueError("oracle affine design requires two observations")
    x, y = zip(*((Fraction(x), Fraction(y)) for x, y in pairs), strict=True)
    sx, sy = sum(x), sum(y)
    xx = sum(value * value for value in x)
    xy = sum(a * b for a, b in zip(x, y, strict=True))
    determinant = n * xx - sx * sx
    if determinant == 0:
        raise ValueError("oracle affine design is singular")
    return (sy * xx - sx * xy) / determinant, (n * xy - sx * sy) / determinant, sy / n


def length_summary(rows: Sequence[tuple[Number, Number | None, Number]]) -> dict[str, Number | None]:
    """Independently score truth/prediction/baseline tuples on the fixed population.

    Args:
        rows: One invented observation per independent representative.

    Returns:
        Exact arithmetic except the square root in RMSE; undefined errors are null.

    Raises:
        ValueError: If no eligible observations exist.
    """
    if not rows:
        raise ValueError("oracle requires an eligible length population")
    paired = [(Fraction(t), Fraction(p), Fraction(b)) for t, p, b in rows if p is not None]
    size = len(paired)
    result: dict[str, Number | None] = {
        "eligible_count": len(rows),
        "assessable_count": size,
        "availability": Fraction(size, len(rows)),
        "within_tolerance": Fraction(sum(abs(p - t) <= max(Fraction(10), t / 10) for t, p, _ in paired), len(rows)),
        "mae": None,
        "rmse": None,
        "bias": None,
        "baseline_mae": None,
        "relative_mae_improvement": None,
    }
    if size:
        absolute = sum(abs(p - t) for t, p, _ in paired) / size
        baseline = sum(abs(b - t) for t, _, b in paired) / size
        result.update(
            mae=absolute,
            rmse=math.sqrt(sum((p - t) ** 2 for t, p, _ in paired) / size),
            bias=sum(p - t for t, p, _ in paired) / size,
            baseline_mae=baseline,
            relative_mae_improvement=None if baseline == 0 else 1 - absolute / baseline,
        )
    return result


def _rate(events: int, total: int) -> Fraction | None:
    return Fraction(events, total) if total else None


def caller_summary(rows: Sequence[Call]) -> Summary:
    """Count a truth/call contingency table, retaining no-calls and missing identity.

    Args:
        rows: Truth, call, optional truth identity, called identity, Tier A identity.

    Returns:
        Counts and fixed-denominator rates without substituting unknown identities.

    Raises:
        ValueError: If no eligible observations exist.
    """
    if not rows:
        raise ValueError("oracle requires an eligible caller population")
    table = Counter((truth, call) for truth, call, *_ in rows)
    positive = sum(table[True, call] for call in (True, False, None))
    negative = sum(table[False, call] for call in (True, False, None))
    no_calls = sum(table[truth, None] for truth in (True, False, None))
    known_identity = [(truth, call, ids, called, tier) for truth, call, ids, called, tier in rows if ids is not None]
    exact_denominator = sum(truth is True for truth, *_ in known_identity)
    recovered = sum(
        truth is True and call is True and set(ids) == set(called) for truth, call, ids, called, _ in known_identity
    )
    return {
        "eligible_count": len(rows),
        "known_truth_count": positive + negative,
        "unknown_truth_count": len(rows) - positive - negative,
        "true_positives": table[True, True],
        "false_negatives": table[True, False],
        "true_negatives": table[False, False],
        "false_positives": table[False, True],
        "positive_no_calls": table[True, None],
        "negative_no_calls": table[False, None],
        "no_calls": no_calls,
        "sensitivity": _rate(table[True, True], positive),
        "specificity": _rate(table[False, False], negative),
        "false_positive_rate": _rate(table[False, True], negative),
        "precision": _rate(table[True, True], table[True, True] + table[False, True]),
        "no_call_rate": _rate(no_calls, len(rows)),
        "assessability": _rate(len(rows) - no_calls, len(rows)),
        "exact_variant_recovery": _rate(recovered, exact_denominator),
        "positive_truth_missing_identity": positive - exact_denominator,
        "identity_assessable_count": len(known_identity),
        "wrong_identity_groups": sum(not set(called).issubset(ids) for _, _, ids, called, _ in known_identity),
        "wrong_tier_a_identity_groups": sum(not set(tier).issubset(ids) for _, _, ids, _, tier in known_identity),
    }


def cutoff_points(
    truths: Sequence[bool],
    scores: Sequence[Fraction | None],
    thresholds: Sequence[Fraction],
    comparison: str,
) -> tuple[tuple[Fraction, tuple[bool | None, ...], Summary], ...]:
    """Create exact invented score decisions and their independent operating points.

    Args:
        truths: Fixed truth population.
        scores: Exact rational scores; null denotes a fixed no-call.
        thresholds: Only cutoffs actually requested, with no fabricated endpoints.
        comparison: One of the four scalar comparisons.

    Returns:
        Strict-to-loose cutoff, calls, and summary triples.

    Raises:
        ValueError: If the comparison or input roster is invalid.
    """
    comparisons = {"<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge}
    if comparison not in comparisons or len(truths) != len(scores) or not truths:
        raise ValueError("oracle cutoff comparison or roster is invalid")
    result = []
    for threshold in sorted(thresholds, reverse=comparison.startswith(">")):
        calls = tuple(None if score is None else comparisons[comparison](score, threshold) for score in scores)
        rows = tuple(
            (truth, call, ("v",) if truth else (), ("v",) if call else (), ())
            for truth, call in zip(truths, calls, strict=True)
        )
        result.append((threshold, calls, caller_summary(rows)))
    return tuple(result)


def binomial_interval(events: int, total: int, tail: Fraction) -> tuple[float, float]:
    """Invert direct exact binomial polynomials, with ``tail`` error per endpoint.

    Args:
        events: Observed successes.
        total: Independent trials, bounded to 100 for this exact small-data oracle.
        tail: Exact one-tail error probability; 1/40 yields a central95% interval.

    Returns:
        Lower and upper endpoints rounded only at the final conversion to float.

    Raises:
        ValueError: If counts, supported size, or tail are invalid.
    """
    if type(events) is not int or type(total) is not int or not 0 <= events <= total <= 100 or total == 0:
        raise ValueError("oracle requires valid counts with 1..100 independent trials")
    if not isinstance(tail, Fraction) or not 0 < tail < Fraction(1, 2):
        raise ValueError("oracle tail must be an exact probability below one half")
    endpoints = []
    for lower in (True, False):
        if (lower and events == 0) or (not lower and events == total):
            endpoints.append(0.0 if lower else 1.0)
            continue
        left, right = Fraction(0), Fraction(1)
        for _ in range(48):
            p = (left + right) / 2
            indices = range(events, total + 1) if lower else range(events + 1)
            probability = sum(math.comb(total, k) * p**k * (1 - p) ** (total - k) for k in indices)
            if (probability < tail) == lower:
                left = p
            else:
                right = p
        endpoints.append(float((left + right) / 2))
    return endpoints[0], endpoints[1]
