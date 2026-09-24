"""Observed-value cutoff axes enumerate exactly the decisions labelled data can change."""

import json
import math
from dataclasses import replace
from fractions import Fraction

import pandas as pd
import pytest

from tests.unit.test_calibration_caller_policy import policy_document, policy_values
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
from vntyper.scripts.calibration_cutoff_axes import (
    ACTIVE_REGION,
    ALT_DEPTH_BAND,
    DEPTH_FLOOR_LINKED,
    DEPTH_SCORE_HIGH,
    GG_GATE_INDEPENDENT,
    STRUCTURAL_GATE_COLUMNS,
    AxisBreakpoints,
    axis_candidates,
    axis_comparison,
    axis_document,
    declared_axis,
    derive_axis,
    eligible_statistic_values,
)
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit

FLOOR = "/components/kestrel/confidence_assignment/reporting_floor"
LOW = "/components/kestrel/confidence_assignment/depth_score_thresholds/low"
HIGH = "/components/kestrel/confidence_assignment/depth_score_thresholds/high"
GG = "/components/kestrel/alt_filtering/gg_depth_score_threshold"
ALT_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/low"
ALT_MID_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low"
ACTIVE = "/components/kestrel/confidence_assignment/var_active_region_threshold"

SHIPPED = {FLOOR: 0.00469, LOW: 0.00469, HIGH: 0.00515, GG: 0.00469}


def baseline(**overrides: object):
    """A decoded Kestrel-only baseline policy, shipped values unless overridden."""
    raw = policy_document(include_advntr=False)
    values = policy_values(raw)
    values.update(SHIPPED)
    values.update(overrides)
    return decode_caller_policy_values(raw)


def frame(*rows: dict[str, object]) -> pd.DataFrame:
    """A prefilter frame whose structural gates default to passing."""
    defaults: dict[str, object] = {
        "Depth_Score": 0.002,
        "Estimated_Depth_AlternateVariant": 30,
        "Estimated_Depth_Variant_ActiveRegion": 400,
        **dict.fromkeys(STRUCTURAL_GATE_COLUMNS, True),
    }
    return pd.DataFrame([{**defaults, **row} for row in rows])


def fractions(*values: float) -> list[Fraction]:
    return [Fraction(value) for value in values]


def test_structural_gate_columns_exclude_the_two_depth_linked_gates():
    assert STRUCTURAL_GATE_COLUMNS == (
        "is_frameshift",
        "is_valid_frameshift",
        "motif_filter_pass",
        "flag_filter_pass",
    )
    assert "depth_confidence_pass" not in STRUCTURAL_GATE_COLUMNS
    assert "alt_filter_pass" not in STRUCTURAL_GATE_COLUMNS


def test_eligible_values_exclude_rows_failing_any_structural_gate():
    rows = [{"Depth_Score": 0.001}]
    rows.extend({"Depth_Score": 0.9, gate: False} for gate in STRUCTURAL_GATE_COLUMNS)
    rows.append({"Depth_Score": 0.8, "flag_filter_pass": "true"})
    rows.append({"Depth_Score": 0.001})
    assert eligible_statistic_values(frame(*rows), "Depth_Score") == (Fraction(0.001),)


def test_eligible_values_are_exact_ascending_and_skip_missing_measurements():
    rows = ({"Depth_Score": 0.003}, {"Depth_Score": 0.001}, {"Depth_Score": float("nan")}, {"Depth_Score": None})
    observed = eligible_statistic_values(frame(*rows), "Depth_Score")
    assert observed == (Fraction(0.001), Fraction(0.003))
    assert all(isinstance(value, Fraction) for value in observed)
    assert eligible_statistic_values(frame(*rows), "Estimated_Depth_AlternateVariant") == (Fraction(30),)


def test_nullable_extension_cells_are_unwrapped_and_unknown_gates_never_pass():
    """A TSV round trip yields numpy scalars and pd.NA, not Python scalars."""
    data = pd.DataFrame(
        {
            "Depth_Score": pd.array([0.004, None, 0.9], dtype="Float64"),
            "Estimated_Depth_AlternateVariant": pd.array([30, 40, 50], dtype="Int64"),
            **{gate: pd.array([True, True, True], dtype="boolean") for gate in STRUCTURAL_GATE_COLUMNS},
        }
    )
    data.loc[2, "flag_filter_pass"] = None
    assert eligible_statistic_values(data, "Depth_Score") == (Fraction(0.004),)
    assert eligible_statistic_values(data, "Estimated_Depth_AlternateVariant") == (Fraction(30), Fraction(40))


def test_eligible_values_reject_malformed_frames_and_cells():
    with pytest.raises(ValueError):
        eligible_statistic_values([], "Depth_Score")
    with pytest.raises(ValueError):
        eligible_statistic_values(frame({}), "")
    with pytest.raises(ValueError):
        eligible_statistic_values(frame({}).drop(columns=["flag_filter_pass"]), "Depth_Score")
    with pytest.raises(ValueError):
        eligible_statistic_values(frame({}), "Missing_Column")
    with pytest.raises(ValueError):
        eligible_statistic_values(frame({"Depth_Score": "0.002"}), "Depth_Score")
    with pytest.raises(ValueError):
        eligible_statistic_values(frame({"Depth_Score": True}), "Depth_Score")


def test_linked_axis_moves_floor_low_and_gg_to_the_same_value():
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001)}, baseline=base)
    assert axis.pointers == tuple(sorted((FLOOR, LOW, GG)))
    assert axis.values == (0.001, 0.00469)
    candidate = axis_candidates(base, axis)[0]
    assert candidate.policy.values[FLOOR] == candidate.policy.values[LOW] == candidate.policy.values[GG] == 0.001
    assert set(candidate.parameters) == {FLOOR, LOW, GG}
    assert candidate.policy.values[HIGH] == base.values[HIGH]


def test_lowering_only_the_reporting_floor_is_not_what_the_linked_axis_does():
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001, 0.002)}, baseline=base)
    for candidate in axis_candidates(base, axis):
        moved = {pointer for pointer in (FLOOR, LOW, GG) if candidate.policy.values[pointer] != base.values[pointer]}
        assert moved in ({FLOOR, LOW, GG}, set())
        floor = candidate.policy.values[FLOOR]
        assert candidate.policy.values[LOW] == floor and candidate.policy.values[GG] == floor


def test_gg_gate_axis_breaks_the_link_deliberately_and_moves_only_the_gg_gate():
    base = baseline()
    axis = derive_axis(GG_GATE_INDEPENDENT, {"s1": fractions(0.001)}, baseline=base)
    assert axis.pointers == (GG,)
    candidate = next(item for item in axis_candidates(base, axis) if item.policy.values[GG] == 0.001)
    assert set(candidate.parameters) == {GG}
    assert candidate.policy.values[FLOOR] == base.values[FLOOR]
    assert candidate.policy.values[LOW] == base.values[LOW]


def test_depth_score_high_axis_relabels_confidence_without_moving_detection_gates():
    base = baseline()
    axis = derive_axis(DEPTH_SCORE_HIGH, {"s1": fractions(0.01)}, baseline=base)
    assert axis.pointers == (HIGH,)
    candidate = next(item for item in axis_candidates(base, axis) if item.policy.values[HIGH] == 0.01)
    assert set(candidate.parameters) == {HIGH}
    assert candidate.policy.values[FLOOR] == base.values[FLOOR] == candidate.policy.values[GG]


def test_alt_depth_band_keeps_mid_low_one_above_low_and_stays_integral():
    base = baseline()
    axis = derive_axis(ALT_DEPTH_BAND, {"s1": [Fraction(10), Fraction(20)]}, baseline=base)
    assert axis.pointers == tuple(sorted((ALT_LOW, ALT_MID_LOW)))
    # 9 is the endpoint sentinel below the smallest observed value of this ``<=`` axis.
    assert axis.values == (9, 10, 20) and axis.sentinel == 9
    assert all(type(value) is int for value in axis.values)
    candidate = axis_candidates(base, axis)[1]
    assert candidate.policy.values[ALT_LOW] == 10 and candidate.policy.values[ALT_MID_LOW] == 11


def test_active_region_axis_is_integral_and_moves_only_its_own_threshold():
    base = baseline()
    axis = derive_axis(ACTIVE_REGION, {"s1": [Fraction(150)]}, baseline=base)
    assert axis.pointers == (ACTIVE,) and axis.values == (149, 150, 200) and axis.sentinel == 149
    assert all(type(value) is int for value in axis.values)
    assert set(axis_candidates(base, axis)[0].parameters) == {ACTIVE}


def test_baseline_value_is_always_present_and_values_are_ascending_and_deduplicated():
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"a": fractions(0.002, 0.001), "b": fractions(0.002)}, baseline=base)
    assert axis.values == (0.001, 0.002, 0.00469)
    assert list(axis.values) == sorted(set(axis.values))
    assert axis.source == "observed-breakpoints"
    assert axis.observed_count == 2 and axis.capped is False and axis.rejected == ()


def test_empty_input_yields_the_baseline_only_axis():
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {}, baseline=base)
    assert axis.values == (0.00469,) and axis.observed_count == 0
    assert axis.capped is False and axis.rejected == ()
    assert derive_axis(ACTIVE_REGION, {"s1": []}, baseline=base).values == (200,)


def test_max_values_caps_in_rank_space_keeping_minimum_maximum_and_baseline():
    base = baseline()
    observed = fractions(0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005)
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": observed}, baseline=base, max_values=5)
    sentinel = math.nextafter(0.005, math.inf)
    # The sentinel above the observed maximum is the new maximum, so capping keeps it.
    assert axis.values == (0.0005, 0.002, 0.0035, 0.00469, sentinel)
    assert axis.capped is True and axis.observed_count == 10 and axis.sentinel == sentinel
    assert min(axis.values) == 0.0005 and max(axis.values) == sentinel and 0.00469 in axis.values


def test_max_values_keeps_a_boundary_baseline_and_does_not_cap_a_short_axis():
    base = baseline()
    observed = fractions(0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004)
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": observed}, baseline=base, max_values=3)
    assert axis.values == (0.0005, 0.0025, 0.00469) and axis.capped is True
    short = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001)}, baseline=base, max_values=3)
    assert short.values == (0.001, 0.00469) and short.capped is False


def test_out_of_range_breakpoints_are_rejected_with_a_reason_rather_than_raising():
    base = baseline()
    # 0.9 is an extreme but admissible depth-score ratio and is kept; only a value
    # outside the decoder's [0, 1] range is refused, and it is refused with a reason.
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001, 0.9, 1.5)}, baseline=base)
    assert axis.values == (0.001, 0.00469, 0.9)
    # The sentinel beyond 1.5 is out of range too, so it is rejected the same way.
    assert [value for value, _ in axis.rejected] == [1.5, math.nextafter(1.5, math.inf)]
    assert axis.sentinel is None
    assert all(isinstance(reason, str) and reason for _, reason in axis.rejected)
    assert axis.observed_count == 3


def test_non_integral_breakpoints_are_rejected_on_integer_axes():
    base = baseline()
    axis = declared_axis(ALT_DEPTH_BAND, [10, 10.5, 99], baseline=base)
    assert axis.values == (10, 20)
    assert [value for value, _ in axis.rejected] == [10.5, 99.0]


def test_an_inadmissible_alternate_depth_baseline_is_refused_before_any_axis_exists():
    """The band partition is a policy invariant, so it fails at construction.

    ``mid_low`` must be ``low + 1`` and ``mid_high`` at least ``mid_low + 1``, so a
    baseline on that ceiling is never decoded in the first place and no axis can be
    derived from it. Asserting the refusal here, rather than around ``derive_axis``,
    keeps the test honest about which call raises.
    """
    with pytest.raises(ValueError, match="alternate-depth partition"):
        baseline(**{ALT_LOW: 99, ALT_MID_LOW: 100})


def test_every_axis_can_hold_its_own_baseline_value():
    """Re-encoding a decoded baseline at its own anchor is always admissible.

    This is what makes the "baseline value is always present" guarantee safe: the
    anchor reproduces an already-valid policy, so it can never be screened out. The
    anchor guard in the screening loop therefore protects a future axis whose
    coupling could invalidate its own anchor, not any axis shipped today.
    """
    base = baseline()
    for axis_name in (DEPTH_FLOOR_LINKED, GG_GATE_INDEPENDENT, DEPTH_SCORE_HIGH, ALT_DEPTH_BAND, ACTIVE_REGION):
        axis = derive_axis(axis_name, {}, baseline=base)
        assert axis.rejected == ()
        assert len(axis.values) == 1


def test_declared_axis_is_baseline_anchored_and_validated():
    base = baseline()
    axis = declared_axis(GG_GATE_INDEPENDENT, [0.002, 0.002, 0.01], baseline=base)
    assert axis.values == (0.002, 0.00469, 0.01)
    assert axis.source == "declared" and axis.observed_count == 2 and axis.capped is False
    for invalid in ([], "0.1", [True], [float("inf")], ["0.1"]):
        with pytest.raises(ValueError):
            declared_axis(GG_GATE_INDEPENDENT, invalid, baseline=base)


def test_axis_candidates_are_complete_policies_with_unique_stable_ids():
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001, 0.002)}, baseline=base)
    candidates = axis_candidates(base, axis)
    assert [item.candidate_id for item in candidates] == [
        "depth_floor_linked-0000",
        "depth_floor_linked-0001",
        "depth_floor_linked-0002",
    ]
    assert len({item.policy.sha256 for item in candidates}) == len(candidates)
    for item in candidates:
        assert set(item.policy.values) == set(base.values)
        with pytest.raises(TypeError):
            item.parameters[FLOOR] = 1.0
    identical = candidates[-1]
    assert identical.policy.sha256 == base.sha256 and dict(identical.parameters) == {}


def test_axis_candidates_refuse_forged_axes_and_unusable_values():
    base = baseline()
    axis = derive_axis(GG_GATE_INDEPENDENT, {"s1": fractions(0.002)}, baseline=base)
    with pytest.raises(ValueError):
        axis_candidates(base, replace(axis, axis="not_an_axis"))
    with pytest.raises(ValueError):
        axis_candidates(base, replace(axis, pointers=(FLOOR,)))
    with pytest.raises(ValueError):
        axis_candidates(base, replace(axis, source="invented"))
    with pytest.raises(ValueError):
        axis_candidates(base, replace(axis, values=()))
    with pytest.raises(ValueError):
        axis_candidates(base, replace(axis, values=(1.5,)))
    with pytest.raises(ValueError):
        axis_candidates(base, object())


def test_axis_document_round_trips():
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001, 1.5)}, baseline=base)
    document = axis_document(axis)
    assert document["schema_version"] == "calibration-cutoff-axis-v1"
    assert document["axis"] == DEPTH_FLOOR_LINKED and document["statistic"] == "Depth_Score"
    assert document["pointers"] == list(axis.pointers) and document["values"] == list(axis.values)
    assert [entry["value"] for entry in document["rejected"]] == [1.5, math.nextafter(1.5, math.inf)]
    assert document["endpoint_sentinel"] is None
    assert json.loads(json.dumps(document)) == document
    assert canonical_sha256(document) == canonical_sha256(axis_document(axis))
    rebuilt = declared_axis(document["axis"], document["values"], baseline=base)
    assert rebuilt.values == axis.values and rebuilt.pointers == axis.pointers
    with pytest.raises(ValueError):
        axis_document(replace(axis, axis="not_an_axis"))


def test_unknown_axis_names_and_malformed_inputs_are_refused():
    base = baseline()
    with pytest.raises(ValueError):
        derive_axis("not_an_axis", {}, baseline=base)
    with pytest.raises(ValueError):
        declared_axis("not_an_axis", [0.1], baseline=base)
    for invalid in ({"": fractions(0.1)}, {1: fractions(0.1)}, {"s1": "0.1"}, {"s1": [0.1]}, []):
        with pytest.raises(ValueError):
            derive_axis(DEPTH_FLOOR_LINKED, invalid, baseline=base)
    for cap in (0, 2, 1.0, True):
        with pytest.raises(ValueError):
            derive_axis(DEPTH_FLOOR_LINKED, {}, baseline=base, max_values=cap)
    with pytest.raises(ValueError):
        derive_axis(DEPTH_FLOOR_LINKED, {}, baseline=object())
    with pytest.raises(ValueError):
        declared_axis(GG_GATE_INDEPENDENT, [0.1], baseline=object())


def test_axis_breakpoints_are_frozen_and_immutable():
    axis = derive_axis(ACTIVE_REGION, {"s1": [Fraction(150)]}, baseline=baseline())
    assert isinstance(axis, AxisBreakpoints)
    with pytest.raises(AttributeError):
        axis.values = ()


def test_linked_axis_reaches_thresholds_above_the_band_edge_by_clamping_low():
    """The upper ROC arm must be reachable, or the curve is truncated.

    Raising the detection floor above the mid-band edge is a valid tightening: the
    floor rule fires before the mid-band rule, so detection narrows. Only ``low`` is
    constrained by ``low <= high``, so the axis clamps the band edge instead of
    refusing the breakpoint. Refusing it would truncate the sweep at 0.00515 and make
    any optimum found at the top of the tested range an artifact of the truncation.
    """
    base = baseline()

    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s": fractions(0.002, 0.01, 0.05)}, baseline=base)

    assert axis.rejected == ()
    assert float(Fraction(1, 100)) in axis.values
    assert float(Fraction(1, 20)) in axis.values
    candidates = {candidate.parameters.get(FLOOR): candidate for candidate in axis_candidates(base, axis)}
    raised = candidates[0.01]
    assert raised.policy.values[FLOOR] == 0.01
    assert raised.policy.values[GG] == 0.01
    assert raised.policy.values[LOW] == SHIPPED[HIGH]
    lowered = candidates[0.002]
    assert lowered.policy.values[LOW] == 0.002


def test_linked_axis_never_strands_a_score_in_the_negative_fallback():
    """Every admitted score must land in a labelled confidence band.

    A score at or above the floor but below ``low`` matches no confidence rule and
    falls through to Negative, which is the defect that made earlier single-gate
    grids change nothing. The axis must keep ``low <= floor`` at every breakpoint.
    """
    base = baseline()

    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s": fractions(0.0009, 0.003, 0.00469, 0.02)}, baseline=base)

    for candidate in axis_candidates(base, axis):
        floor = candidate.policy.values[FLOOR]
        assert candidate.policy.values[LOW] <= floor
        assert candidate.policy.values[GG] == floor


def test_a_threshold_just_above_the_data_is_reachable_on_a_floor_axis():
    """Inclusive floors pass the maximum itself, so rejecting everything needs a sentinel.

    A positive at 0.005 and a negative at 0.01: without the sentinel the axis is
    (0.00469, 0.005, 0.01) and every candidate passes the negative, so no candidate reaches
    specificity 1.0. The next float above 0.01 rejects both rows.
    """
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"positive": fractions(0.005), "negative": fractions(0.01)}, baseline=base)
    sentinel = math.nextafter(0.01, math.inf)

    assert axis.values == (0.00469, 0.005, 0.01, sentinel)
    assert axis.sentinel == sentinel and axis.observed_count == 2
    top = axis_candidates(base, axis)[-1]
    assert top.policy.values[FLOOR] == sentinel and top.policy.values[GG] == sentinel
    assert top.policy.values[LOW] == SHIPPED[HIGH]
    assert not top.policy.values[FLOOR] <= 0.01  # the production comparator now rejects the maximum


def test_the_sentinel_of_a_less_or_equal_axis_lies_below_the_minimum():
    """The minimum sits between the shipped low (0.00469) and high, so ``low <= high`` holds."""
    base = baseline()
    axis = derive_axis(DEPTH_SCORE_HIGH, {"s1": fractions(0.005, 0.02)}, baseline=base)
    sentinel = math.nextafter(0.005, -math.inf)

    assert axis.values == (sentinel, 0.005, 0.00515, 0.02) and axis.sentinel == sentinel


def test_no_sentinel_when_the_baseline_already_lies_beyond_the_observed_range():
    base = baseline()
    floor = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001, 0.002)}, baseline=base)
    high = derive_axis(DEPTH_SCORE_HIGH, {"s1": fractions(0.01)}, baseline=base)

    assert floor.values == (0.001, 0.002, 0.00469) and floor.sentinel is None
    assert high.values == (0.00515, 0.01) and high.sentinel is None
    assert derive_axis(DEPTH_FLOOR_LINKED, {}, baseline=base).sentinel is None


def test_the_sentinel_is_added_when_the_baseline_equals_the_observed_edge():
    base = baseline()
    axis = derive_axis(DEPTH_FLOOR_LINKED, {"s1": fractions(0.001, 0.00469)}, baseline=base)

    assert axis.values == (0.001, 0.00469, math.nextafter(0.00469, math.inf))


def test_a_sentinel_the_decoder_refuses_lands_in_rejected():
    """An integer ``<=`` axis at zero cannot go to -1; the refusal is reported, not raised."""
    base = baseline()
    axis = derive_axis(ACTIVE_REGION, {"s1": [Fraction(0), Fraction(150)]}, baseline=base)

    assert axis.values == (0, 150, 200) and axis.sentinel is None
    assert [value for value, _ in axis.rejected] == [-1.0]


def test_a_sentinel_is_strictly_beyond_an_edge_that_is_not_a_float():
    """A Fraction with no exact float still gets a sentinel strictly beyond it."""
    base = baseline()
    edge = Fraction(1, 30)
    upward = derive_axis(DEPTH_FLOOR_LINKED, {"s1": [edge]}, baseline=base)
    downward = derive_axis(DEPTH_SCORE_HIGH, {"s1": [Fraction(1, 199)]}, baseline=base)

    assert upward.sentinel is not None and Fraction(upward.sentinel) > edge
    assert downward.sentinel is not None and Fraction(downward.sentinel) < Fraction(1, 199)


def test_declared_axes_never_carry_a_sentinel_and_forged_sentinels_are_refused():
    base = baseline()
    declared = declared_axis(GG_GATE_INDEPENDENT, [0.002, 0.01], baseline=base)
    derived = derive_axis(GG_GATE_INDEPENDENT, {"s1": fractions(0.01)}, baseline=base)

    assert declared.sentinel is None and axis_document(declared)["endpoint_sentinel"] is None
    assert axis_document(derived)["endpoint_sentinel"] == math.nextafter(0.01, math.inf)
    with pytest.raises(ValueError, match="sentinel"):
        axis_candidates(base, replace(declared, sentinel=0.01))
    with pytest.raises(ValueError, match="sentinel"):
        axis_candidates(base, replace(derived, sentinel=0.5))


def test_every_axis_declares_its_production_comparator():
    for name in (DEPTH_FLOOR_LINKED, GG_GATE_INDEPENDENT):
        assert axis_comparison(name) == ">="
    for name in (DEPTH_SCORE_HIGH, ALT_DEPTH_BAND, ACTIVE_REGION):
        assert axis_comparison(name) == "<="
    with pytest.raises(ValueError):
        axis_comparison("not_an_axis")
