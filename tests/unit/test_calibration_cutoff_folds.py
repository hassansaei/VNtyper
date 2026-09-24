"""Fold-local axis inventories: every value a fold may select comes from its training samples.

The regression cases are the reviewer's counterexamples (spec §14.1). A full-data axis lets
a held-out sample choose its own fold's reject-everything sentinel, and a full-data cap lets
a held-out value displace a training breakpoint. Deriving each fold's axis from that fold's
training samples alone removes both.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from fractions import Fraction

import pytest

from tests.unit.test_calibration_cutoff_advntr import _policy
from tests.unit.test_calibration_cutoff_axes import baseline as kestrel_baseline
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_cutoff_advntr_axes import derive_advntr_axis
from vntyper.scripts.calibration_cutoff_axes import (
    ADVNTR_CUTOFF,
    DEPTH_FLOOR_LINKED,
    GG_GATE_INDEPENDENT,
    AxisBreakpoints,
    axis_candidates,
    axis_document,
    derive_axis,
)
from vntyper.scripts.calibration_cutoff_folds import fold_axis, fold_candidate_inventories

pytestmark = pytest.mark.unit


def _up(value: float) -> float:
    return math.nextafter(value, math.inf)


def _advntr(cutoff: float, max_values: int | None = None):
    policy: CallerPolicyValues = _policy(cutoff=cutoff, support=3)

    def derive(statistics: Mapping[str, Fraction]) -> AxisBreakpoints:
        return derive_advntr_axis(ADVNTR_CUTOFF, statistics, baseline=policy, max_values=max_values)

    return derive


def _kestrel(axis: str = DEPTH_FLOOR_LINKED, max_values: int | None = None):
    policy = kestrel_baseline()

    def derive(values: Mapping[str, Sequence[Fraction]]) -> AxisBreakpoints:
        return derive_axis(axis, values, baseline=policy, max_values=max_values)

    return derive


def test_a_held_out_sample_never_becomes_its_own_folds_advntr_sentinel() -> None:
    """Baseline 0.5; training negative q=0.2, training positive without a statistic, held-out negative q=0.1.

    A full-data sentinel ``min(q) = 0.1`` rejects the held-out negative at exactly its own
    value. Its fold's training data only knows ``q = 0.2``, so the fold's sentinel is 0.2.
    """
    statistics = {"neg-train": Fraction(0.2), "neg-held": Fraction(0.1)}
    assignments = {"neg-train": 1, "pos-train": 1, "neg-held": 0}

    merged, folds = fold_axis(_advntr(0.5), statistics, assignments)

    assert folds[0] == frozenset({0.2, _up(0.2), 0.5})
    assert 0.2 in folds[0] and 0.1 not in folds[0]
    assert folds[1] == frozenset({0.1, _up(0.1), 0.5})
    # The full-data axis has the leaky sentinel; the fold's training sentinel is fold-only.
    assert merged.sentinel == 0.1
    assert merged.values == (0.1, _up(0.1), 0.2, _up(0.2), 0.5)
    assert merged.fold_only == 1 and merged.capped is False and merged.fold_capped == ()


def test_a_held_out_value_cannot_displace_a_training_breakpoint_under_the_cap() -> None:
    """Training minima {0.1, 0.3, 0.5, 0.7}, baseline 0.2, cap 4, held-out q = 0.4.

    Training-only: [0.1 (sentinel), up(0.1), 0.2, up(0.3), up(0.5), up(0.7)] -> ranks
    {0, 2, 3, 5} -> {0.1, 0.2, up(0.3), up(0.7)}. Adding the held-out 0.4 shifts the ranks of
    the full-data axis to {0, 2, 4, 6} -> {0.1, 0.2, up(0.4), up(0.7)}: up(0.4) displaces up(0.3).
    """
    training = {name: Fraction(value) for name, value in zip("abcd", (0.1, 0.3, 0.5, 0.7), strict=True)}
    statistics = {**training, "held": Fraction(0.4)}
    assignments = {"a": 1, "b": 1, "c": 1, "d": 1, "held": 0}
    derive = _advntr(0.2, max_values=4)

    merged, folds = fold_axis(derive, statistics, assignments)

    assert folds[0] == frozenset(derive(training).values)
    assert folds[0] == frozenset({0.1, 0.2, _up(0.3), _up(0.7)})
    assert _up(0.4) not in folds[0]
    assert _up(0.4) in derive(statistics).values  # the leak the fold derivation removes
    assert merged.capped is True
    assert all(len(values) <= 4 for values in folds.values())


def test_a_held_out_maximum_never_becomes_its_own_folds_kestrel_sentinel() -> None:
    """``depth_floor_linked`` is ``>=``: its sentinel sits just above the largest training score."""
    values = {"low": (Fraction(0.004),), "mid": (Fraction(0.006),), "top": (Fraction(0.009),)}
    assignments = {"low": 0, "mid": 0, "top": 1}

    merged, folds = fold_axis(_kestrel(), values, assignments)

    assert _up(0.006) in folds[1] and _up(0.009) not in folds[1]
    assert 0.009 not in folds[1]
    assert folds[1] == frozenset({0.004, 0.00469, 0.006, _up(0.006)})
    assert merged.sentinel == _up(0.009)
    assert _up(0.006) in merged.values and merged.fold_only >= 1


def test_a_fold_capped_axis_is_reported_as_a_capped_subsample() -> None:
    """Scores {0.002, 0.004, 0.008, 1.1}, baseline 0.00469, cap 4.

    Full data: 1.1 and its sentinel exceed the gate's range and are rejected, leaving
    {0.002, 0.004, 0.00469, 0.008} -- four values, uncapped. The fold holding 1.1 out gains
    the sentinel up(0.008), five values, and is capped.
    """
    scores = {"a": 0.002, "b": 0.004, "c": 0.008, "d": 1.1}
    values = {key: (Fraction(score),) for key, score in scores.items()}
    assignments = {"a": 0, "b": 1, "c": 2, "d": 3}

    merged, folds = fold_axis(_kestrel(max_values=4), values, assignments)

    assert merged.capped is False
    assert merged.fold_capped == (3,)
    assert all(len(inventory) <= 4 for inventory in folds.values())
    document = axis_document(merged)
    assert document["full_data_capped"] is False and document["fold_capped"] == [3]
    assert document["breakpoint_completeness"] == "capped-subsample"


def test_without_folds_the_full_axis_stands_alone() -> None:
    derive = _kestrel()
    values = {"only": (Fraction(0.003),)}

    merged, folds = fold_axis(derive, values, {})

    assert merged == derive(values) and folds == {}
    assert fold_candidate_inventories([(merged, axis_candidates(kestrel_baseline(), merged))], {merged.axis: {}}) == {}


def test_a_sample_without_a_fold_is_refused() -> None:
    with pytest.raises(ValueError, match="fold"):
        fold_axis(_kestrel(), {"stray": (Fraction(0.003),)}, {"other": 0, "another": 1})


def test_fold_inventories_map_each_folds_values_to_candidate_ids_across_axes() -> None:
    policy = kestrel_baseline()
    values = {"low": (Fraction(0.004),), "mid": (Fraction(0.006),), "top": (Fraction(0.009),)}
    assignments = {"low": 0, "mid": 0, "top": 1}
    linked, linked_folds = fold_axis(_kestrel(), values, assignments)
    gate, gate_folds = fold_axis(_kestrel(GG_GATE_INDEPENDENT), values, assignments)
    derived = [(linked, axis_candidates(policy, linked)), (gate, axis_candidates(policy, gate))]

    inventories = fold_candidate_inventories(derived, {linked.axis: linked_folds, gate.axis: gate_folds})

    assert set(inventories) == {0, 1}
    for fold in (0, 1):
        expected = {
            candidate.candidate_id
            for axis, candidates in derived
            for value, candidate in zip(axis.values, candidates, strict=True)
            if value in (linked_folds if axis is linked else gate_folds)[fold]
        }
        assert inventories[fold] == frozenset(expected)
    by_value = dict(zip(linked.values, (c.candidate_id for c in derived[0][1]), strict=True))
    assert by_value[_up(0.006)] in inventories[1] and by_value[_up(0.009)] not in inventories[1]


def test_fold_inventories_refuse_values_without_a_candidate_and_mismatched_axes() -> None:
    policy = kestrel_baseline()
    values = {"low": (Fraction(0.004),), "top": (Fraction(0.009),)}
    axis, folds = fold_axis(_kestrel(), values, {"low": 0, "top": 1})
    derived = [(axis, axis_candidates(policy, axis))]

    with pytest.raises(ValueError, match="candidate"):
        fold_candidate_inventories(derived, {axis.axis: {**folds, 0: folds[0] | {0.5}}})
    with pytest.raises(ValueError, match="axes"):
        fold_candidate_inventories(derived, {GG_GATE_INDEPENDENT: folds})
