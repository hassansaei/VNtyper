"""A curve needs one moving axis; anything else is a labelled table, not a ROC curve."""

import json
from dataclasses import replace
from fractions import Fraction
from types import MappingProxyType

import pytest

from tests.unit.test_calibration_caller_policy import policy_document, policy_values
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, decode_caller_policy_values
from vntyper.scripts.calibration_cutoff_axes import GG_GATE_INDEPENDENT, AxisBreakpoints
from vntyper.scripts.calibration_cutoff_curves import (
    AxisCurve,
    JointOperatingPoint,
    axis_curve_document,
    build_axis_curve,
    build_joint_points,
    joint_points_document,
)
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate

pytestmark = pytest.mark.unit

GG = "/components/kestrel/alt_filtering/gg_depth_score_threshold"
FLOOR = "/components/kestrel/confidence_assignment/reporting_floor"
BASELINE_GG = 0.00469
VALUES = (0.002, BASELINE_GG, 0.01)
PAIR = (0.002, BASELINE_GG)
TRUTH: dict[str, bool | None] = {"p1": True, "p2": True, "n1": False, "n2": False}


def policy(overrides: dict[str, object] | None = None) -> CallerPolicyValues:
    """A decoded Kestrel-only policy with the shipped values unless overridden."""
    raw = policy_document(include_advntr=False)
    policy_values(raw).update(overrides or {})
    return decode_caller_policy_values(raw)


def axis(values: tuple[float, ...] = VALUES, rejected: tuple[tuple[float, str], ...] = ()) -> AxisBreakpoints:
    """The single-pointer GG axis, anchored on the shipped GG threshold."""
    return AxisBreakpoints(
        axis=GG_GATE_INDEPENDENT,
        pointers=(GG,),
        values=values,
        source="observed-breakpoints",
        observed_count=len(values),
        capped=False,
        rejected=rejected,
    )


def candidate(index: int, values: tuple[float, ...] = VALUES) -> CutoffCandidate:
    """The complete policy one breakpoint implies, with the axis id convention."""
    value = values[index]
    changed = {} if value == BASELINE_GG else {GG: value}
    return CutoffCandidate(f"{GG_GATE_INDEPENDENT}-{index:04d}", policy({GG: value}), MappingProxyType(dict(changed)))


def observations(
    calls: dict[str, bool | None], truth: dict[str, bool | None] | None = None
) -> tuple[CallerObservation, ...]:
    """One independent-group representative per sample at a fixed operating point."""
    labels = TRUTH if truth is None else truth
    return tuple(
        CallerObservation(
            key,
            f"group-{key}",
            labels[key],
            None if labels[key] is None else (("v",) if labels[key] else ()),
            calls[key],
            ("v",) if calls[key] else (),
            (),
        )
        for key in sorted(labels)
    )


def scenario(
    calls_by_index: dict[int, dict[str, bool | None]],
    *,
    values: tuple[float, ...] = VALUES,
    truth: dict[str, bool | None] | None = None,
) -> tuple[AxisBreakpoints, tuple[CutoffCandidate, ...], dict[str, tuple[CallerObservation, ...]]]:
    """An axis, its candidates and the replayed arm for each of them."""
    candidates = tuple(candidate(index, values) for index in sorted(calls_by_index))
    arms = {
        item.candidate_id: observations(calls_by_index[index], truth)
        for index, item in zip(sorted(calls_by_index), candidates, strict=True)
    }
    return axis(values), candidates, arms


def curve(
    calls_by_index: dict[int, dict[str, bool | None]],
    *,
    values: tuple[float, ...] = VALUES,
    truth: dict[str, bool | None] | None = None,
) -> AxisCurve:
    """Build the one-axis curve for a replayed scenario."""
    breakpoints, candidates, arms = scenario(calls_by_index, values=values, truth=truth)
    return build_axis_curve(breakpoints, candidates, arms, comparison=">=", phase="policy-selection")


SEPARABLE: dict[int, dict[str, bool | None]] = {
    2: {"p1": True, "p2": False, "n1": False, "n2": False},
    1: {"p1": True, "p2": True, "n1": False, "n2": False},
    0: {"p1": True, "p2": True, "n1": True, "n2": False},
}


def test_clean_separable_cohort_is_perfect_at_the_shipped_threshold() -> None:
    result = curve(SEPARABLE)

    assert result.axis == GG_GATE_INDEPENDENT
    assert result.statistic == "Depth_Score"
    assert result.pointers == (GG,)
    assert result.baseline_threshold == BASELINE_GG
    assert [float(point.threshold) for point in result.curves.points] == [0.01, BASELINE_GG, 0.002]
    assert [(point.sensitivity, point.false_positive_rate) for point in result.curves.points] == [
        (Fraction(1, 2), Fraction(0)),
        (Fraction(1), Fraction(0)),
        (Fraction(1), Fraction(1, 2)),
    ]
    shipped = result.interval_by_threshold[repr(BASELINE_GG)]
    assert shipped["sensitivity"]["estimate"] == Fraction(1)
    assert shipped["specificity"]["estimate"] == Fraction(1)
    assert result.boundary_support == {
        "positives_within_band": 1,
        "negatives_within_band": 1,
        "band_low": Fraction(0.002),
        "band_high": Fraction(0.01),
    }


def test_loosening_trades_one_false_positive_for_two_true_positives() -> None:
    truth: dict[str, bool | None] = {"p1": True, "p2": True, "p3": True, "n1": False, "n2": False}
    result = curve(
        {
            1: {"p1": True, "p2": False, "p3": False, "n1": False, "n2": False},
            0: {"p1": True, "p2": True, "p3": True, "n1": True, "n2": False},
        },
        values=PAIR,
        truth=truth,
    )

    strict, loose = result.curves.points
    assert (strict.true_positives, strict.false_positives) == (1, 0)
    assert (loose.true_positives, loose.false_positives) == (3, 1)
    assert loose.true_positives - strict.true_positives == 2
    assert loose.false_positives - strict.false_positives == 1
    assert result.boundary_support["positives_within_band"] == 2
    assert result.boundary_support["negatives_within_band"] == 1


def test_off_axis_variation_is_refused_by_the_offending_candidate() -> None:
    breakpoints, candidates, arms = scenario(SEPARABLE)
    forged = replace(candidates[0], policy=policy({GG: 0.002, FLOOR: 0.001}))
    with pytest.raises(ValueError, match=f"{forged.candidate_id}.*{FLOOR}"):
        build_axis_curve(breakpoints, (forged, *candidates[1:]), arms, comparison=">=", phase="policy-selection")


def test_a_candidate_with_no_replayed_arm_is_refused() -> None:
    breakpoints, candidates, arms = scenario(SEPARABLE)
    del arms[candidates[0].candidate_id]
    with pytest.raises(ValueError, match=f"{candidates[0].candidate_id}.*no replayed observations"):
        build_axis_curve(breakpoints, candidates, arms, comparison=">=", phase="policy-selection")


@pytest.mark.parametrize("mode", ["group", "truth", "missing"])
def test_roster_disagreement_between_arms_is_refused(mode: str) -> None:
    breakpoints, candidates, arms = scenario(SEPARABLE)
    rows = arms[candidates[0].candidate_id]
    if mode == "group":
        rows = (replace(rows[0], group_key="other"), *rows[1:])
    elif mode == "truth":
        rows = (replace(rows[0], truth_positive=True, truth_variants=("v",)), *rows[1:])
    else:
        rows = rows[1:]
    arms[candidates[0].candidate_id] = rows
    with pytest.raises(ValueError, match=f"{candidates[0].candidate_id} replays a different roster"):
        build_axis_curve(breakpoints, candidates, arms, comparison=">=", phase="policy-selection")


def test_unknown_truth_keeps_its_own_denominator_and_stays_out_of_sensitivity() -> None:
    truth: dict[str, bool | None] = {"p1": True, "p2": True, "n1": False, "n2": False, "u1": None}
    result = curve(
        {
            1: {"p1": True, "p2": True, "n1": False, "n2": False, "u1": False},
            0: {"p1": True, "p2": True, "n1": False, "n2": False, "u1": True},
        },
        values=PAIR,
        truth=truth,
    )

    assert result.curves.eligible_count == 5
    assert result.curves.unknown_truth_count == 1
    loose = result.interval_by_threshold[repr(0.002)]
    assert loose["sensitivity"]["events"] == 2
    assert loose["sensitivity"]["total"] == 2
    assert loose["precision"]["total"] == 2
    assert result.boundary_support["positives_within_band"] == 0
    assert result.boundary_support["negatives_within_band"] == 0


def test_no_calls_stay_inside_their_truth_class_denominator() -> None:
    result = curve(
        {
            1: {"p1": True, "p2": None, "n1": False, "n2": False},
            0: {"p1": True, "p2": None, "n1": True, "n2": False},
        },
        values=PAIR,
    )

    assert result.curves.no_call_count == 1
    assert result.curves.positive_no_calls == 1
    assert result.curves.negative_no_calls == 0
    shipped = result.interval_by_threshold[repr(BASELINE_GG)]
    assert (shipped["sensitivity"]["events"], shipped["sensitivity"]["total"]) == (1, 2)
    assert shipped["sensitivity"]["estimate"] == Fraction(1, 2)


def test_tiny_denominators_widen_the_exact_clopper_pearson_bounds() -> None:
    truth: dict[str, bool | None] = {"p1": True, "n1": False}
    tiny = curve(
        {
            1: {"p1": True, "n1": False},
            0: {"p1": True, "n1": True},
        },
        values=PAIR,
        truth=truth,
    )

    one = tiny.interval_by_threshold[repr(BASELINE_GG)]
    assert one["sensitivity"]["lower"] == Fraction(1, 40)
    assert one["sensitivity"]["upper"] == Fraction(1)
    assert one["specificity"]["lower"] == Fraction(1, 40)
    two = curve(SEPARABLE).interval_by_threshold[repr(BASELINE_GG)]
    assert two["sensitivity"]["lower"] > one["sensitivity"]["lower"]


def test_zero_observed_false_positives_still_admit_a_material_rate() -> None:
    result = curve(SEPARABLE)

    shipped = result.interval_by_threshold[repr(BASELINE_GG)]
    assert shipped["false_positive_rate_one_sided_upper"] == Fraction(684455749451, 881583902934)
    assert shipped["false_positive_rate_one_sided_upper"] > 0


def test_boundary_support_is_zero_when_no_sample_lies_in_the_band() -> None:
    calls: dict[str, bool | None] = {"p1": True, "p2": False, "n1": False, "n2": False}
    result = curve({1: dict(calls), 0: dict(calls)}, values=PAIR)

    assert result.boundary_support["positives_within_band"] == 0
    assert result.boundary_support["negatives_within_band"] == 0
    assert result.boundary_support["band_low"] == Fraction(0.002)
    assert result.boundary_support["band_high"] == Fraction(BASELINE_GG)


@pytest.mark.parametrize(
    ("mode", "message"),
    [
        ("axis", "require AxisBreakpoints"),
        ("empty", "nonempty candidate"),
        ("arms", "keyed by candidate"),
        ("type", "CutoffCandidate"),
        ("foreign", "not a breakpoint"),
        ("range", "outside the breakpoints"),
        ("duplicate", "repeats a breakpoint"),
        ("pointers", "different policy pointers"),
        ("unanchored", "exactly one baseline"),
        ("doubled", "exactly one baseline"),
    ],
)
def test_forged_or_unanchored_axis_inputs_are_refused(mode: str, message: str) -> None:
    breakpoints, candidates, arms = scenario(SEPARABLE)
    supplied: object = candidates
    forged: object = breakpoints
    if mode == "axis":
        forged = "gg"
    elif mode == "empty":
        supplied = ()
    elif mode == "arms":
        arms = []  # type: ignore[assignment]
    elif mode == "type":
        supplied = (candidates[0].candidate_id, *candidates[1:])
    elif mode == "foreign":
        supplied = (replace(candidates[0], candidate_id="depth_score_high-0000"), *candidates[1:])
    elif mode == "range":
        supplied = (replace(candidates[0], candidate_id=f"{GG_GATE_INDEPENDENT}-0009"), *candidates[1:])
    elif mode == "duplicate":
        supplied = (candidates[0], candidates[0], *candidates[1:])
    elif mode == "pointers":
        supplied = (replace(candidates[0], policy=policy_with_advntr()), *candidates[1:])
    elif mode == "unanchored":
        supplied = (candidates[0], candidates[2])
    else:
        supplied = (replace(candidates[0], parameters=MappingProxyType({})), *candidates[1:])
    with pytest.raises(ValueError, match=message):
        build_axis_curve(forged, supplied, arms, comparison=">=", phase="policy-selection")  # type: ignore[arg-type]


def policy_with_advntr() -> CallerPolicyValues:
    """A policy whose pointer inventory differs from the Kestrel-only baseline."""
    raw = policy_document(include_advntr=True)
    policy_values(raw)[GG] = 0.002
    return decode_caller_policy_values(raw)


def test_axis_curve_document_round_trips_with_counts_and_rejected_breakpoints() -> None:
    breakpoints, candidates, arms = scenario(SEPARABLE)
    breakpoints = replace(breakpoints, rejected=((0.9, "gg_depth_score_threshold must be between zero and one"),))
    result = build_axis_curve(breakpoints, candidates, arms, comparison=">=", phase="policy-selection")

    document = axis_curve_document(result)

    assert json.loads(json.dumps(document)) == document
    assert document["schema_version"] == "calibration-cutoff-curve-v1"
    assert document["statistic"] == "Depth_Score"
    assert document["comparison"] == ">="
    assert document["thresholds"] == [0.01, BASELINE_GG, 0.002]
    assert document["baseline_threshold"] == BASELINE_GG
    assert document["rejected"] == [{"value": 0.9, "reason": "gg_depth_score_threshold must be between zero and one"}]
    rows = document["points"]
    assert [row["threshold"] for row in rows] == [0.01, BASELINE_GG, 0.002]
    assert [
        (row["true_positives"], row["false_positives"], row["true_negatives"], row["false_negatives"]) for row in rows
    ] == [(1, 0, 2, 1), (2, 0, 2, 0), (2, 1, 1, 0)]
    assert all(row["positive_no_calls"] == row["negative_no_calls"] == row["unknown_truth_count"] == 0 for row in rows)
    assert rows[1]["intervals"]["sensitivity"]["estimate"] == 1.0
    assert rows[0]["precision"] == 1.0
    assert document["boundary_support"]["band_high"] == 0.01


def test_axis_curve_document_refuses_a_foreign_object() -> None:
    with pytest.raises(ValueError, match="AxisCurve"):
        axis_curve_document("curve")  # type: ignore[arg-type]


def joint(
    label_values: dict[str, dict[str, object]],
) -> dict[str, tuple[CutoffCandidate, tuple[CallerObservation, ...]]]:
    """Labelled multi-parameter operating points over the same replayed cohort."""
    calls: dict[str, bool | None] = {"p1": True, "p2": True, "n1": True, "n2": False}
    return {
        label: (
            CutoffCandidate(f"joint-{index:04d}", policy(values), MappingProxyType(dict(values))),
            observations(calls),
        )
        for index, (label, values) in enumerate(label_values.items())
    }


def test_joint_points_carry_every_differing_pointer_and_are_label_ordered() -> None:
    labelled = joint({"zebra": {GG: 0.002, FLOOR: 0.001}, "alpha": {GG: 0.01}})

    points = build_joint_points(labelled, candidate(1))

    assert [point.label for point in points] == ["alpha", "zebra"]
    assert dict(points[0].values) == {GG: 0.01}
    assert dict(points[1].values) == {FLOOR: 0.001, GG: 0.002}
    assert tuple(points[1].values) == (GG, FLOOR)
    assert points[0].metrics.true_positives == 2
    assert points[0].metrics.false_positives == 1
    with pytest.raises(TypeError):
        points[0].values[GG] = 0.5  # type: ignore[index]


def test_joint_points_document_round_trips() -> None:
    points = build_joint_points(joint({"zebra": {GG: 0.002, FLOOR: 0.001}, "alpha": {GG: 0.01}}), candidate(1))

    document = joint_points_document(points)

    assert json.loads(json.dumps(document)) == document
    assert document["schema_version"] == "calibration-cutoff-joint-v1"
    assert document["point_count"] == 2
    rows = document["points"]
    assert [row["label"] for row in rows] == ["alpha", "zebra"]
    assert rows[1]["values"] == {FLOOR: 0.001, GG: 0.002}
    assert rows[0]["metrics"]["true_positives"] == 2
    assert rows[0]["metrics"]["sensitivity"]["estimate"] == 1.0


@pytest.mark.parametrize(
    ("mode", "message"),
    [
        ("baseline", "CutoffCandidate baseline"),
        ("empty", "nonempty labelled mapping"),
        ("label", "nonempty trimmed strings"),
        ("pair", "candidate and its observations"),
        ("candidate", "requires a CutoffCandidate"),
        ("rows", "sequence of observations"),
        ("pointers", "different policy pointers"),
    ],
)
def test_forged_joint_inputs_are_refused(mode: str, message: str) -> None:
    labelled: dict[str, object] = dict(joint({"alpha": {GG: 0.01}}))
    baseline: object = candidate(1)
    if mode == "baseline":
        baseline = "baseline"
    elif mode == "empty":
        labelled = {}
    elif mode == "label":
        labelled[" alpha "] = labelled.pop("alpha")
    elif mode == "pair":
        labelled["alpha"] = (candidate(0),)
    elif mode == "candidate":
        labelled["alpha"] = ("candidate", observations({"p1": True, "p2": True, "n1": True, "n2": False}))
    elif mode == "rows":
        labelled["alpha"] = (candidate(0), "rows")
    else:
        labelled["alpha"] = (
            CutoffCandidate("joint-0000", policy_with_advntr(), MappingProxyType({})),
            observations({"p1": True, "p2": True, "n1": True, "n2": False}),
        )
    with pytest.raises(ValueError, match=message):
        build_joint_points(labelled, baseline)  # type: ignore[arg-type]


def test_joint_points_document_refuses_a_foreign_row() -> None:
    with pytest.raises(ValueError, match="JointOperatingPoint"):
        joint_points_document(("alpha",))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="sequence of joint"):
        joint_points_document("alpha")  # type: ignore[arg-type]


def test_joint_operating_point_is_frozen() -> None:
    point = build_joint_points(joint({"alpha": {GG: 0.01}}), candidate(1))[0]
    assert isinstance(point, JointOperatingPoint)
    with pytest.raises(AttributeError):
        point.label = "other"  # type: ignore[misc]
