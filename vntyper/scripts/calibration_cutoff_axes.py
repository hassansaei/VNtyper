"""Decision breakpoints derived from labelled data, and the complete policies they imply.

A Kestrel cutoff is a threshold compared against a continuous per-row statistic with a
``>=``/``<=`` comparator. Such a comparison changes a decision only at values the data
actually take: between two adjacent observed values every threshold yields the identical
partition of the cohort. The complete set of distinct outcomes is therefore obtained by
testing exactly the observed values, and nothing is gained by a hand-written grid --
which can only ever be a coarser, arbitrarily placed subsample of this set, and which
silently omits the breakpoints that separate the samples it was tuned on.

This module enumerates those breakpoints and turns each one into a *complete* validated
policy, so a curve over an axis is a curve over real, replayable operating points.

Why the depth gates move together
---------------------------------
``vntyper/scripts/confidence_rules.py`` evaluates ``CONFIDENCE_RULES`` as an ordered
first-match table:

* Row 0 ``subthreshold_or_nan``: ``Depth_Score`` NaN or ``< confidence_assignment.reporting_floor`` -> ``Negative``.
* Row 1 ``midband_demotion``: ``Depth_Score`` in the **closed** interval
  ``[depth_score_thresholds.low, depth_score_thresholds.high]`` -> ``Low_Precision``.
* Rows 2-5 all additionally require ``Depth_Score >= high`` (row 2) or ``> high`` (rows 3-5).
* Row 6 ``fallback_negative`` -> ``Negative``.

``vntyper/scripts/confidence_assignment.py`` then sets
``depth_confidence_pass = Confidence != "Negative"``.

A candidate whose ``Depth_Score`` lies between a **lowered** ``reporting_floor`` and an
**unchanged** ``depth_score_thresholds.low`` therefore escapes row 0, matches no other
row -- rows 1-5 all need at least ``low`` -- and falls through to row 6, which labels it
``Negative`` again. **Lowering the floor alone changes nothing.** It moves a row from one
``Negative`` verdict to another.

``vntyper/scripts/variant_parsing.py::filter_by_alt_values_and_finalize`` adds a second,
independent gate: a row whose ``ALT == "GG"`` passes ``alt_filter_pass`` only when
``Depth_Score >= alt_filtering.gg_depth_score_threshold``. The MUC1 dupC candidate is
literally ``ALT == "GG"``, so for the variant this pipeline exists to find, that gate
holds the row even once rows 0 and 1 have been satisfied.

All three gates must move together for a detection change, which is what
:data:`DEPTH_FLOOR_LINKED` -- the default depth axis -- does. :data:`GG_GATE_INDEPENDENT`
exists to measure the GG gate on its own and deliberately breaks that link.

Bounds, and what gets dropped
-----------------------------
Breakpoints are screened by constructing the candidate policy and offering it to
``decode_caller_policy_values``; whatever that function rejects is dropped with its own
message rather than raised, and is reported in :attr:`AxisBreakpoints.rejected`. The
bounds it enforces (``vntyper/scripts/calibration_caller_policy.py``) are: the four
depth-score pointers are finite numbers in ``[0, 1]``; ``var_active_region_threshold``
and the three alternate-depth pointers are non-negative integers; ``low <= high``; and
``alt_depth_thresholds`` must satisfy ``mid_low == low + 1`` and ``mid_high >= mid_low + 1``.
Two consequences are worth stating because they are not obvious: the linked depth axis is
bounded above by the baseline ``depth_score_thresholds.high`` (it sets ``low``, and
``low <= high`` is an invariant), and the alternate-depth band is bounded above by
``mid_high - 2``.

The attainable endpoint
-----------------------
Production floors are inclusive: ``Depth_Score < floor`` rejects, and the GG gate passes
``Depth_Score >= threshold``. A threshold equal to the largest observed value therefore
still passes that value, so the observed values alone never contain the operating point
that rejects every row. Without it a specificity floor can be reported unreachable, or a
worse cutoff chosen, even though a threshold just above the data attains it. The same
holds mirrored for a ``<=`` axis: a threshold equal to the smallest value still has that
value on the ``<=`` side, so the partition in which every row lies above the threshold
is missing. :func:`derive_axis` therefore adds one endpoint sentinel on the side of the
observed range that the comparator leaves unrepresented: the next representable float
above the largest value for a ``>=`` axis, and the next one below the smallest for a
``<=`` axis (``+1``/``-1`` on integer axes). It is recorded in :attr:`AxisBreakpoints.sentinel`,
excluded from ``observed_count``, screened like any other value, and omitted when the
baseline already lies beyond the observed range (the baseline then realizes that
partition). Being an extreme, it survives rank-space capping like the minimum and maximum.
Declared axes never receive one.

Exactness
---------
Observed values are carried as :class:`fractions.Fraction` so that unioning, sorting,
deduplication and rank subsampling are exact and never depend on floating-point
associativity. The conversion to ``float`` happens once, at the moment a breakpoint
becomes a policy value, and :attr:`AxisBreakpoints.values` stores the converted result --
so the axis document and the policy it produces agree bit for bit. Integer axes convert
to ``int`` instead, and a non-integral breakpoint on such an axis is rejected rather
than rounded.

The two adVNTR axes
--------------------
:data:`ADVNTR_CUTOFF` and :data:`ADVNTR_MIN_SUPPORT` register the other caller's two
tunable pointers alongside the five Kestrel axes above, but their comparators are
adVNTR's own p-value gate: ``calibrated_calling.cutoff`` is a strict ``<`` (a call needs
a p-value strictly below the threshold), and ``calibrated_calling.minimum_read_support``
is an inclusive ``>=`` (a call needs read support at or above the threshold). Neither
axis's breakpoints are enumerated from a Kestrel-style prefilter frame -- there is no
per-row ``Depth_Score`` equivalent for adVNTR's statistics -- so they are built by the
public :func:`observed_axis` from breakpoints a native replay already computed; see
``calibration_cutoff_advntr_axes`` for how those breakpoints and their statistics are
derived. :func:`derive_axis` refuses both axes outright, so the Kestrel path can never
be asked to enumerate them from data it cannot read.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
from types import MappingProxyType
from typing import Final, NoReturn, cast

import pandas as pd

from vntyper.scripts.calibration_caller_policy import (
    ADVNTR_CALLER_POLICY_POINTERS,
    KESTREL_CALLER_POLICY_POINTERS,
    CallerPolicyScalar,
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate

logger = logging.getLogger(__name__)

#: Stable axis names. They are part of the artifact contract and never derived from a
#: pointer spelling, so a pointer rename cannot silently rename a published axis.
DEPTH_FLOOR_LINKED: Final[str] = "depth_floor_linked"
GG_GATE_INDEPENDENT: Final[str] = "gg_gate_independent"
DEPTH_SCORE_HIGH: Final[str] = "depth_score_high"
ALT_DEPTH_BAND: Final[str] = "alt_depth_band"
ACTIVE_REGION: Final[str] = "var_active_region"
ADVNTR_CUTOFF: Final[str] = "advntr_cutoff"
ADVNTR_MIN_SUPPORT: Final[str] = "advntr_min_support"

#: Every final Kestrel gate except the two depth-linked ones. ``depth_confidence_pass``
#: and ``alt_filter_pass`` are exactly what an axis moves, so a row may not be excluded
#: from the observed set for failing them -- that would hide the breakpoints the axis
#: exists to find. The other four encode a structural or artifact judgement ("this row is
#: not a candidate variant at all"), which no cutoff can overturn. See
#: ``vntyper/scripts/subthreshold.py`` for the same distinction drawn for #266.
STRUCTURAL_GATE_COLUMNS: Final[tuple[str, ...]] = (
    "is_frameshift",
    "is_valid_frameshift",
    "motif_filter_pass",
    "flag_filter_pass",
)

_SCHEMA: Final[str] = "calibration-cutoff-axis-v1"
_SOURCES: Final[frozenset[str]] = frozenset({"observed-breakpoints", "declared"})
_MINIMUM_CAP: Final[int] = 3

# Resolved out of the policy module's own inventory: a pointer renamed there raises here
# at import time instead of producing an axis that quietly moves nothing.
_KESTREL: Final[Mapping[str, str]] = MappingProxyType(
    {pointer.removeprefix("/components/kestrel/"): pointer for pointer in KESTREL_CALLER_POLICY_POINTERS}
)
_FLOOR: Final[str] = _KESTREL["confidence_assignment/reporting_floor"]
_LOW: Final[str] = _KESTREL["confidence_assignment/depth_score_thresholds/low"]
_HIGH: Final[str] = _KESTREL["confidence_assignment/depth_score_thresholds/high"]
_GG: Final[str] = _KESTREL["alt_filtering/gg_depth_score_threshold"]
_ALT_LOW: Final[str] = _KESTREL["confidence_assignment/alt_depth_thresholds/low"]
_ALT_MID_LOW: Final[str] = _KESTREL["confidence_assignment/alt_depth_thresholds/mid_low"]
_ACTIVE: Final[str] = _KESTREL["confidence_assignment/var_active_region_threshold"]

_ADVNTR: Final[Mapping[str, str]] = MappingProxyType(
    {pointer.removeprefix("/components/advntr/"): pointer for pointer in ADVNTR_CALLER_POLICY_POINTERS}
)
_ADV_CUTOFF: Final[str] = _ADVNTR["calibrated_calling/cutoff"]
_ADV_SUPPORT: Final[str] = _ADVNTR["calibrated_calling/minimum_read_support"]


@dataclass(frozen=True)
class AxisBreakpoints:
    """One axis of complete decision breakpoints anchored on the shipped operating point.

    Attributes:
        axis: Stable axis name, one of the module-level axis constants.
        pointers: Sorted decision-profile JSON pointers this axis moves together.
        values: Ascending, deduplicated candidate values, always containing the baseline
            value of the axis's primary pointer. ``float`` for the depth-score axes and
            ``int`` for the integer ones.
        source: ``"observed-breakpoints"`` or ``"declared"``.
        observed_count: Distinct statistic values seen before any cap -- the size of the
            observed union for a derived axis, or of the declared set for a declared one.
            It counts rejected values too, and excludes the injected baseline unless the
            data also observed it.
        capped: True when ``values`` is a rank-space subsample of the accepted set.
        rejected: Ascending ``(value, reason)`` pairs dropped because the resulting policy
            would not decode; the reason is the decoder's own message.
        sentinel: The endpoint value added beyond the observed range, when it was added
            and accepted; it is then also a member of ``values``. ``None`` otherwise.
        fold_only: How many of ``values`` only a fold's training-derived inventory
            produced (see :func:`merge_axis_values`); 0 for an axis derived once.
        fold_capped: Ascending outer folds whose training-derived inventory was capped.
            ``capped`` describes the full-data inventory only.
    """

    axis: str
    pointers: tuple[str, ...]
    values: tuple[float, ...]
    source: str
    observed_count: int
    capped: bool
    rejected: tuple[tuple[float, str], ...]
    sentinel: float | None = None
    fold_only: int = 0
    fold_capped: tuple[int, ...] = ()


@dataclass(frozen=True)
class _AxisSpec:
    """Static definition of one axis: what it moves and what it reads."""

    primary: str
    pointers: tuple[str, ...]
    statistic: str
    integer: bool
    #: The production comparator the axis sweeps, ``>=``, ``<=`` or ``<``.
    comparison: str
    #: The caller whose policy pointers this axis moves, ``"kestrel"`` or ``"advntr"``.
    #: Kept last so existing positional constructions of this dataclass keep working.
    caller: str


_SPECS: Final[Mapping[str, _AxisSpec]] = MappingProxyType(
    {
        DEPTH_FLOOR_LINKED: _AxisSpec(
            _FLOOR, tuple(sorted((_FLOOR, _LOW, _GG))), "Depth_Score", False, ">=", "kestrel"
        ),
        GG_GATE_INDEPENDENT: _AxisSpec(_GG, (_GG,), "Depth_Score", False, ">=", "kestrel"),
        DEPTH_SCORE_HIGH: _AxisSpec(_HIGH, (_HIGH,), "Depth_Score", False, "<=", "kestrel"),
        ALT_DEPTH_BAND: _AxisSpec(
            _ALT_LOW, tuple(sorted((_ALT_LOW, _ALT_MID_LOW))), "Estimated_Depth_AlternateVariant", True, "<=", "kestrel"
        ),
        ACTIVE_REGION: _AxisSpec(_ACTIVE, (_ACTIVE,), "Estimated_Depth_Variant_ActiveRegion", True, "<=", "kestrel"),
        ADVNTR_CUTOFF: _AxisSpec(
            _ADV_CUTOFF, (_ADV_CUTOFF,), "advntr_min_pvalue_at_baseline_support", False, "<", "advntr"
        ),
        ADVNTR_MIN_SUPPORT: _AxisSpec(
            _ADV_SUPPORT, (_ADV_SUPPORT,), "advntr_max_read_support_below_baseline_cutoff", True, ">=", "advntr"
        ),
    }
)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _spec(axis: object) -> _AxisSpec:
    if not isinstance(axis, str) or axis not in _SPECS:
        _fail(f"cutoff axis name must be one of {sorted(_SPECS)}")
    return _SPECS[axis]


def axis_comparison(axis: str) -> str:
    """The production comparator one axis sweeps, ``>=``, ``<=`` or ``<``.

    Args:
        axis: One of the module-level axis names.

    Returns:
        The comparator: a row passes a ``>=`` axis at or above the threshold, sits on
        the ``<=`` side of a ``<=`` axis at or below it, and passes :data:`ADVNTR_CUTOFF`
        only strictly below the threshold.

    Raises:
        ValueError: For an unknown axis name.
    """
    return _spec(axis).comparison


def axis_caller(axis: str) -> str:
    """The caller whose policy pointers one axis moves.

    Args:
        axis: One of the module-level axis names.

    Returns:
        ``"kestrel"`` for the five Kestrel axes, ``"advntr"`` for :data:`ADVNTR_CUTOFF`
        and :data:`ADVNTR_MIN_SUPPORT`.

    Raises:
        ValueError: For an unknown axis name.
    """
    return _spec(axis).caller


def _endpoint_sentinel(spec: _AxisSpec, observed: set[Fraction], anchor: Fraction) -> Fraction | None:
    """The value just beyond the observed range on the side the comparator leaves unrepresented.

    Returns None for an empty observation set, and when the baseline already lies strictly
    beyond the observed range, because the baseline then realizes the same partition.
    """
    if not observed:
        return None
    upward = spec.comparison == ">="
    edge = max(observed) if upward else min(observed)
    if (anchor > edge) if upward else (anchor < edge):
        return None
    if spec.integer:
        return Fraction(math.floor(edge) + 1) if upward else Fraction(math.ceil(edge) - 1)
    target = math.inf if upward else -math.inf
    value = math.nextafter(float(edge), target)
    # ``float(edge)`` rounds a Fraction that did not come from a float; step until the
    # sentinel is strictly beyond the exact edge.
    while (Fraction(value) <= edge) if upward else (Fraction(value) >= edge):
        value = math.nextafter(value, target)
    return Fraction(value)


def _passes(value: object) -> bool:
    """One structural gate verdict, where only an explicit True counts as passing."""
    item = getattr(value, "item", None)
    if callable(item):
        value = item()
    return value is True


def _number(value: object) -> Fraction | None:
    """The exact value of one statistic cell, or None when it holds no measurement."""
    if value is None or value is pd.NA:
        return None
    item = getattr(value, "item", None)
    if callable(item):
        value = item()
    if isinstance(value, bool):
        _fail("cutoff axis statistic cells must be numbers, never booleans")
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, float):
        return Fraction(value) if math.isfinite(value) else None
    _fail("cutoff axis statistic cells must be numeric")


def _scalar(spec: _AxisSpec, value: Fraction) -> int | float | None:
    """The policy value a breakpoint becomes, or None when the axis cannot express it."""
    if not spec.integer:
        return float(value)
    return int(value) if value.denominator == 1 else None


def _overrides(baseline: CallerPolicyValues, spec: _AxisSpec, scalar: int | float) -> dict[str, CallerPolicyScalar]:
    """Resolve every pointer one breakpoint moves, keeping the policy admissible.

    The linked depth axis is the detection floor. Two production rules constrain it.
    A score at or above the floor but below the mid-band edge matches no confidence
    rule and falls through to Negative, so the edge must never sit above the floor.
    The policy decoder separately requires ``low <= high``. Setting the edge to
    ``min(floor, high)`` satisfies both: below the shipped edge the three gates move
    together, and above it the edge stays at ``high`` so the floor can keep rising.
    Clamping rather than refusing is what keeps the tightening arm of the curve
    reachable; refusing would truncate every sweep at the shipped edge.

    Args:
        baseline: Complete baseline policy supplying the unmoved values.
        spec: Static definition of the axis being moved.
        scalar: One admissible breakpoint on that axis.

    Returns:
        Pointer-to-value overrides to apply on top of the baseline.
    """
    if spec.primary == _FLOOR:
        high = cast("int | float", baseline.values[_HIGH])
        return {_FLOOR: scalar, _LOW: min(scalar, high), _GG: scalar}
    if spec.primary == _ALT_LOW:
        return {_ALT_LOW: scalar, _ALT_MID_LOW: cast(int, scalar) + 1}
    return {spec.primary: scalar}


def _policy(baseline: CallerPolicyValues, spec: _AxisSpec, scalar: int | float) -> CallerPolicyValues:
    raw = caller_policy_values_document(baseline)
    values = dict(baseline.values)
    values.update(_overrides(baseline, spec, scalar))
    return decode_caller_policy_values({**raw, "values": values})


def _baseline_value(baseline: CallerPolicyValues, spec: _AxisSpec, axis: str) -> Fraction:
    caller_policy_values_document(baseline)
    if spec.primary not in baseline.values:
        _fail(f"cutoff axis {axis} requires a baseline policy that includes {spec.caller}")
    value = baseline.values[spec.primary]
    return Fraction(cast("int | float", value))


def _screen(
    baseline: CallerPolicyValues, spec: _AxisSpec, candidates: Sequence[Fraction], anchor: Fraction
) -> tuple[list[tuple[Fraction, int | float]], list[tuple[float, str]]]:
    accepted: list[tuple[Fraction, int | float]] = []
    rejected: list[tuple[float, str]] = []
    for value in candidates:
        scalar = _scalar(spec, value)
        if scalar is None:
            rejected.append((float(value), f"{spec.primary} requires an integral breakpoint value"))
            continue
        try:
            _policy(baseline, spec, scalar)
        except ValueError as error:
            if value == anchor:
                _fail(f"cutoff axis {spec.primary} cannot hold its own baseline value: {error}")
            rejected.append((float(value), str(error)))
            continue
        accepted.append((value, scalar))
    return accepted, rejected


def _ranks(count: int, maximum: int, anchor: int) -> tuple[int, ...]:
    """Evenly spaced ranks including the first, the last and the anchor, never rounded values."""
    ranks = {round(index * (count - 1) / (maximum - 1)) for index in range(maximum)}
    if anchor not in ranks:
        interior = ranks - {0, count - 1}
        victim = min(interior, key=lambda rank: (abs(rank - anchor), rank))
        ranks = (ranks - {victim}) | {anchor}
    return tuple(sorted(ranks))


def _build(
    axis: str,
    spec: _AxisSpec,
    baseline: CallerPolicyValues,
    candidates: Sequence[Fraction],
    *,
    anchor: Fraction,
    observed_count: int,
    source: str,
    max_values: int | None,
    sentinel: Fraction | None = None,
) -> AxisBreakpoints:
    accepted, rejected = _screen(baseline, spec, candidates, anchor)
    capped = max_values is not None and len(accepted) > max_values
    if capped:
        index = [value for value, _ in accepted].index(anchor)
        accepted = [accepted[rank] for rank in _ranks(len(accepted), cast(int, max_values), index)]
    kept = next((scalar for value, scalar in accepted if value == sentinel), None)
    return AxisBreakpoints(
        axis=axis,
        pointers=spec.pointers,
        values=tuple(scalar for _, scalar in accepted),
        source=source,
        observed_count=observed_count,
        capped=capped,
        rejected=tuple(rejected),
        sentinel=kept if sentinel is not None else None,
    )


def _require_axis(axis: object) -> _AxisSpec:
    if not isinstance(axis, AxisBreakpoints):
        _fail("cutoff axis operations require AxisBreakpoints")
    spec = _spec(axis.axis)
    if axis.pointers != spec.pointers or axis.source not in _SOURCES or not axis.values:
        _fail(f"cutoff axis {axis.axis} content differs from its definition")
    if axis.sentinel is not None and (axis.source != "observed-breakpoints" or axis.sentinel not in axis.values):
        _fail(f"cutoff axis {axis.axis} sentinel must be one of its observed-breakpoint values")
    return spec


def eligible_statistic_values(prefilter: pd.DataFrame, statistic: str) -> tuple[Fraction, ...]:
    """Distinct values of one statistic among rows passing every structural gate.

    Rows failing ``depth_confidence_pass`` or ``alt_filter_pass`` are deliberately kept:
    those are the gates an axis moves, and excluding them would hide exactly the
    breakpoints that change a call. Rows failing any of :data:`STRUCTURAL_GATE_COLUMNS`
    are excluded, because no cutoff can turn a non-candidate into a candidate. A gate cell
    passes only when it is explicitly ``True``; a string spelling is unknown, not passing.

    Args:
        prefilter: Kestrel pre-filter frame carrying the statistic and the four gates.
        statistic: Column name of the continuous or integer statistic to enumerate.

    Returns:
        Ascending deduplicated exact values; cells holding no measurement are skipped.

    Raises:
        ValueError: If the frame, the column name or a statistic cell is malformed.
    """
    if not isinstance(prefilter, pd.DataFrame):
        _fail("cutoff axis eligibility requires a prefilter DataFrame")
    if not isinstance(statistic, str) or not statistic:
        _fail("cutoff axis statistic must be a nonempty column name")
    columns = (statistic, *STRUCTURAL_GATE_COLUMNS)
    missing = [column for column in columns if column not in prefilter.columns]
    if missing:
        _fail(f"cutoff axis prefilter is missing columns {missing}")
    values: set[Fraction] = set()
    for row in prefilter[list(columns)].itertuples(index=False, name=None):
        if not all(_passes(gate) for gate in row[1:]):
            continue
        number = _number(row[0])
        if number is not None:
            values.add(number)
    return tuple(sorted(values))


def derive_axis(
    axis: str,
    values_by_sample: Mapping[str, Sequence[Fraction]],
    *,
    baseline: CallerPolicyValues,
    max_values: int | None = None,
) -> AxisBreakpoints:
    """Union the per-sample observed values into ascending decision breakpoints.

    The union is taken over samples because a threshold is a cohort-wide decision: a value
    observed in one sample is a breakpoint for every sample it is compared against. One
    endpoint sentinel just beyond the observed range is added, above the maximum for a
    ``>=`` axis and below the minimum for a ``<=`` axis, so the partition no observed value
    can represent is reachable (see the module docstring).

    Args:
        axis: One of the five Kestrel module-level axis names.
        values_by_sample: Per-sample exact values, normally from
            :func:`eligible_statistic_values`. An empty mapping is valid.
        baseline: The shipped complete policy; its primary-pointer value is always kept.
        max_values: Optional cap of at least three, applied by even rank-space
            subsampling that retains the minimum, the maximum and the baseline value. The
            sentinel is the minimum or the maximum, so it is retained too.

    Returns:
        The axis, with breakpoints the policy decoder refuses recorded in ``rejected``.

    Raises:
        ValueError: If the axis name, the mapping, the cap or the baseline is invalid, if
            the axis is an adVNTR axis, or if the baseline value itself cannot be encoded
            on this axis.
    """
    spec = _spec(axis)
    if spec.caller != "kestrel":
        _fail(f"cutoff axis {axis} is an adVNTR axis; derive it with derive_advntr_axis")
    if max_values is not None and (type(max_values) is not int or max_values < _MINIMUM_CAP):
        _fail(f"cutoff axis max_values must be an integer of at least {_MINIMUM_CAP}")
    if not isinstance(values_by_sample, Mapping) or any(
        not isinstance(key, str) or not key for key in values_by_sample
    ):
        _fail("cutoff axis values must be keyed by nonempty sample identifiers")
    observed: set[Fraction] = set()
    for sample in sorted(values_by_sample):
        items = values_by_sample[sample]
        if isinstance(items, str) or not isinstance(items, Sequence):
            _fail(f"cutoff axis values for sample {sample} must be a sequence")
        for item in items:
            if not isinstance(item, Fraction):
                _fail(f"cutoff axis values for sample {sample} must be exact Fractions")
            observed.add(item)
    anchor = _baseline_value(baseline, spec, axis)
    sentinel = _endpoint_sentinel(spec, observed, anchor)
    extra = {anchor} if sentinel is None else {anchor, sentinel}
    return _build(
        axis,
        spec,
        baseline,
        sorted(observed | extra),
        anchor=anchor,
        observed_count=len(observed),
        source="observed-breakpoints",
        max_values=max_values,
        sentinel=sentinel,
    )


def declared_axis(axis: str, values: Sequence[float], *, baseline: CallerPolicyValues) -> AxisBreakpoints:
    """Build an axis from explicitly declared values, still baseline-anchored.

    This is the escape hatch for replaying a published axis or for probing values the
    cohort never produced. It applies the identical screening and anchoring as
    :func:`derive_axis`, so a declared axis is not a way around the policy bounds.

    Args:
        axis: One of the module-level axis names.
        values: Nonempty finite declared values; duplicates are collapsed.
        baseline: The shipped complete policy; its primary-pointer value is always kept.

    Returns:
        The axis with ``source="declared"`` and ``observed_count`` set to the number of
        distinct declared values.

    Raises:
        ValueError: If the axis name, the declared values or the baseline is invalid, or
            if the baseline value itself cannot be encoded on this axis.
    """
    spec = _spec(axis)
    if isinstance(values, str) or not isinstance(values, Sequence) or not values:
        _fail("declared cutoff axis values must be a nonempty sequence")
    declared: set[Fraction] = set()
    for item in values:
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item):
            _fail("declared cutoff axis values must be finite numbers, never booleans")
        declared.add(Fraction(item))
    anchor = _baseline_value(baseline, spec, axis)
    return _build(
        axis,
        spec,
        baseline,
        sorted(declared | {anchor}),
        anchor=anchor,
        observed_count=len(declared),
        source="declared",
        max_values=None,
    )


def observed_axis(
    axis: str,
    candidates: Sequence[Fraction],
    *,
    baseline: CallerPolicyValues,
    observed_count: int,
    sentinel: Fraction | None,
    max_values: int | None,
) -> AxisBreakpoints:
    """Build an axis from breakpoints a caller already computed, still baseline-anchored.

    This is the public counterpart of :func:`derive_axis` for axes whose breakpoints are
    not enumerated from a Kestrel-style prefilter frame -- most importantly the two adVNTR
    axes, whose breakpoints come from replaying capture evidence natively (see
    ``calibration_cutoff_advntr_axes``). Unlike :func:`derive_axis`, it does not union
    per-sample values or compute the endpoint sentinel itself; the caller supplies both
    the already-unioned candidates and the sentinel, and this function only adds the
    baseline anchor, screens every value against the policy decoder, and applies the same
    rank-space capping.

    Args:
        axis: One of the module-level axis names.
        candidates: Exact breakpoints already observed, deduplicated or not.
        baseline: The shipped complete policy; its primary-pointer value is always kept.
        observed_count: Distinct statistic values seen before any cap, as defined on
            :attr:`AxisBreakpoints.observed_count`.
        sentinel: The endpoint sentinel to record, or ``None`` when the caller determined
            none is needed.
        max_values: Optional cap of at least three, applied by even rank-space
            subsampling that retains the minimum, the maximum and the baseline value.

    Returns:
        The axis, with breakpoints the policy decoder refuses recorded in ``rejected``.

    Raises:
        ValueError: If the axis name, the candidates, the cap or the baseline is invalid,
            or if the baseline value itself cannot be encoded on this axis.
    """
    spec = _spec(axis)
    if max_values is not None and (type(max_values) is not int or max_values < _MINIMUM_CAP):
        _fail(f"cutoff axis max_values must be an integer of at least {_MINIMUM_CAP}")
    if isinstance(candidates, str) or not isinstance(candidates, Sequence):
        _fail(f"cutoff axis {axis} candidates must be a sequence of Fractions")
    for item in candidates:
        if not isinstance(item, Fraction):
            _fail(f"cutoff axis {axis} candidates must be exact Fractions")
    anchor = _baseline_value(baseline, spec, axis)
    extra = {anchor} | ({sentinel} if sentinel is not None else set())
    return _build(
        axis,
        spec,
        baseline,
        sorted(set(candidates) | extra),
        anchor=anchor,
        observed_count=observed_count,
        source="observed-breakpoints",
        max_values=max_values,
        sentinel=sentinel,
    )


def merge_axis_values(axis: AxisBreakpoints, extra: Iterable[int | float]) -> AxisBreakpoints:
    """Union extra breakpoints, normally fold inventories, into an axis.

    Args:
        axis: The full-data axis.
        extra: Values other inventories of the same axis hold; members already on the axis
            are not counted again.

    Returns:
        The axis with ascending ``values`` covering both sets and ``fold_only`` counting
        the values only ``extra`` contributed. Every other field is unchanged.

    Raises:
        ValueError: If the axis is forged, or an extra value is not a finite number
            (never a boolean), or is not an integer on an integer axis.
    """
    spec = _require_axis(axis)
    added: set[int | float] = set()
    for value in extra:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            _fail(f"cutoff axis {axis.axis} merged values must be finite numbers, never booleans")
        if spec.integer and not isinstance(value, int):
            _fail(f"cutoff axis {axis.axis} merged values must be integers")
        added.add(value)
    return replace(axis, values=tuple(sorted(set(axis.values) | added)), fold_only=len(added - set(axis.values)))


def axis_candidates(baseline: CallerPolicyValues, axis: AxisBreakpoints) -> tuple[CutoffCandidate, ...]:
    """One complete validated policy per breakpoint, with stable ``<axis>-NNNN`` ids.

    Every candidate is a whole policy, not a delta: the pointers the axis does not name
    keep their baseline values, so a candidate is replayable on its own and a curve over
    the axis is a curve over real operating points. The candidate whose value equals the
    baseline reproduces the baseline policy exactly, with empty ``parameters``.

    Args:
        baseline: The shipped complete policy the candidates are built from.
        axis: Breakpoints produced by :func:`derive_axis`, :func:`declared_axis` or
            :func:`observed_axis`.

    Returns:
        One candidate per breakpoint in ascending value order, ids ordered with them.

    Raises:
        ValueError: If the axis is forged, the baseline is invalid, or a breakpoint the
            axis carries does not decode into a policy.
    """
    spec = _require_axis(axis)
    caller_policy_values_document(baseline)
    candidates: list[CutoffCandidate] = []
    for index, value in enumerate(axis.values):
        policy = _policy(baseline, spec, value)
        changed = {
            pointer: scalar
            for pointer, scalar in _overrides(baseline, spec, value).items()
            if baseline.values[pointer] != scalar
        }
        candidates.append(CutoffCandidate(f"{axis.axis}-{index:04d}", policy, MappingProxyType(changed)))
    return tuple(candidates)


def axis_document(axis: AxisBreakpoints) -> dict[str, object]:
    """Project one axis as fresh canonical ``calibration-cutoff-axis-v1`` JSON.

    Args:
        axis: Breakpoints produced by :func:`derive_axis`, :func:`declared_axis` or
            :func:`observed_axis`.

    Returns:
        A JSON-compatible object carrying the axis name, the caller it moves, the
        production comparator it sweeps, the statistic it reads, the pointers it moves,
        its values, provenance, cap state, every rejected value with the reason it was
        dropped, and the endpoint sentinel (or null). Completeness is stated per
        inventory: ``breakpoint_completeness`` is ``"complete"`` only when neither the
        full-data inventory nor any fold's training-derived inventory was capped, because a
        capped inventory may miss distinct operating points.

    Raises:
        ValueError: If the axis content differs from its module-level definition.
    """
    spec = _require_axis(axis)
    return {
        "schema_version": _SCHEMA,
        "axis": axis.axis,
        "caller": spec.caller,
        "comparator": spec.comparison,
        "statistic": spec.statistic,
        "pointers": list(axis.pointers),
        "values": list(axis.values),
        "source": axis.source,
        "observed_count": axis.observed_count,
        # The v2.0.39 key, kept for calibration-cutoff-axis-v1 readers; full_data_capped names it.
        "capped": axis.capped,
        "rejected": [{"value": value, "reason": reason} for value, reason in axis.rejected],
        "endpoint_sentinel": axis.sentinel,
        "fold_only_values": axis.fold_only,
        "full_data_capped": axis.capped,
        "fold_capped": list(axis.fold_capped),
        "breakpoint_completeness": "complete" if not axis.capped and not axis.fold_capped else "capped-subsample",
    }
