"""adVNTR cutoff axes: statistics observed in native replay, breakpoints, checks."""

from __future__ import annotations

import json
import math
import subprocess
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path

import pytest

from tests.unit.advntr_grid_fakes import FIRST_VNTR_ID, advntr_grid_result, parity_capture
from tests.unit.advntr_grid_fakes import visit_document as _visit
from tests.unit.test_calibration_cutoff_advntr import _capture, _json_bytes, _policy, _ReplayTool
from vntyper.scripts.calibration_caller_policy import (
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_cutoff_advntr import (
    AdvntrCutoffGridResult,
    evaluate_advntr_cutoff_grid,
)
from vntyper.scripts.calibration_cutoff_advntr_axes import (
    PROBE_LADDERS,
    AdvntrVisit,
    advntr_baseline_parity,
    advntr_probe_policy,
    check_replay_consistency,
    derive_advntr_axis,
    predicted_call,
    probe_visits,
    sample_statistics,
    unrejectable_samples,
)
from vntyper.scripts.calibration_cutoff_axes import ADVNTR_CUTOFF, ADVNTR_MIN_SUPPORT, axis_candidates
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate

pytestmark = pytest.mark.unit

_CUT = "/components/advntr/calibrated_calling/cutoff"
_SUP = "/components/advntr/calibrated_calling/minimum_read_support"
_MODE = "/components/advntr/calibrated_calling/mode"


def _v(support: int, p: float) -> AdvntrVisit:
    return AdvntrVisit(support, Fraction(p))


def _with(policy: CallerPolicyValues, pointer: str, value: object) -> CallerPolicyValues:
    document = caller_policy_values_document(policy)
    values: dict[str, object] = dict(policy.values)
    values[pointer] = value
    return decode_caller_policy_values({**document, "values": values})


def _kestrel_only(policy: CallerPolicyValues) -> CallerPolicyValues:
    document = caller_policy_values_document(policy)
    values = {pointer: value for pointer, value in policy.values.items() if "/advntr/" not in pointer}
    return decode_caller_policy_values({**document, "required_callers": ["kestrel"], "values": values})


VISITS = {
    "pos-a": (_v(5, 0.0004), _v(2, 1e-9)),  # min p at s>=3: 0.0004 ; max s at p<0.001: 5
    "pos-b": (_v(4, 0.004), _v(9, 0.02)),  # min p: 0.004 ; max s with p<0.001: none
    "neg-c": (_v(3, 0.003),),  # min p: 0.003 ; none
    "unassessable": None,
    "empty": (),
}


def test_probe_policy_loosens_only_the_axis_pointer() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    cutoff_probe = advntr_probe_policy(ADVNTR_CUTOFF, baseline)
    assert cutoff_probe.values[_CUT] == math.nextafter(1.0, 0.0)
    assert cutoff_probe.values[_SUP] == 3
    support_probe = advntr_probe_policy(ADVNTR_MIN_SUPPORT, baseline)
    assert support_probe.values[_SUP] == 1 and support_probe.values[_CUT] == 0.001
    changed = {p for p in baseline.values if baseline.values[p] != cutoff_probe.values[p]}
    assert changed == {_CUT}
    changed = {p for p in baseline.values if baseline.values[p] != support_probe.values[p]}
    assert changed == {_SUP}


def test_probe_ladders_loosen_first() -> None:
    assert PROBE_LADDERS[ADVNTR_CUTOFF] == (math.nextafter(1.0, 0.0), 0.5, 0.1, 0.05)
    assert PROBE_LADDERS[ADVNTR_MIN_SUPPORT] == (1, 2)


def test_probe_policy_falls_back_to_the_next_admissible_rung(monkeypatch: pytest.MonkeyPatch) -> None:
    import vntyper.scripts.calibration_cutoff_advntr_axes as module

    baseline = _policy(cutoff=0.001, support=3)
    monkeypatch.setattr(module, "PROBE_LADDERS", {ADVNTR_CUTOFF: (1.5, 0.0, 0.25), ADVNTR_MIN_SUPPORT: (0, -1)})
    assert advntr_probe_policy(ADVNTR_CUTOFF, baseline).values[_CUT] == 0.25
    with pytest.raises(ValueError, match="advntr_min_support has no admissible permissive projection"):
        advntr_probe_policy(ADVNTR_MIN_SUPPORT, baseline)


def test_probe_policy_refuses_a_non_legacy_baseline() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    exact = _with(baseline, _MODE, "exact")
    with pytest.raises(ValueError, match="legacy calibrated_calling mode; the baseline mode is 'exact'"):
        advntr_probe_policy(ADVNTR_CUTOFF, exact)


def test_every_entry_point_refuses_an_advntr_less_baseline_and_an_unknown_axis() -> None:
    baseline = _kestrel_only(_policy())
    with pytest.raises(ValueError, match="require a baseline policy that includes adVNTR"):
        advntr_probe_policy(ADVNTR_CUTOFF, baseline)
    with pytest.raises(ValueError, match="require a baseline policy that includes adVNTR"):
        sample_statistics(ADVNTR_CUTOFF, VISITS, baseline)
    with pytest.raises(ValueError, match="require a baseline policy that includes adVNTR"):
        derive_advntr_axis(ADVNTR_CUTOFF, {}, baseline=baseline, max_values=None)
    legacy = _policy()
    for call in (
        lambda: advntr_probe_policy("kestrel_floor", legacy),
        lambda: sample_statistics("kestrel_floor", VISITS, legacy),
        lambda: derive_advntr_axis("kestrel_floor", {}, baseline=legacy, max_values=None),
        lambda: unrejectable_samples("kestrel_floor", {}),
    ):
        with pytest.raises(ValueError, match=r"adVNTR cutoff axis must be one of \['advntr_cutoff', 'advntr_min"):
            call()


def test_sample_statistics_follow_the_legacy_rule() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    assert sample_statistics(ADVNTR_CUTOFF, VISITS, baseline) == {
        "pos-a": Fraction(0.0004),
        "pos-b": Fraction(0.004),
        "neg-c": Fraction(0.003),
    }
    assert sample_statistics(ADVNTR_MIN_SUPPORT, VISITS, baseline) == {"pos-a": Fraction(5)}


def test_sample_statistics_use_inclusive_support_and_strict_cutoff() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    visits = {"edge": (_v(3, 0.001), _v(2, 0.0001))}
    # support 3 is eligible (>=); p == cutoff is not below the cutoff (<).
    assert sample_statistics(ADVNTR_CUTOFF, visits, baseline) == {"edge": Fraction(0.001)}
    assert sample_statistics(ADVNTR_MIN_SUPPORT, visits, baseline) == {"edge": Fraction(2)}


def test_cutoff_breakpoints_are_the_next_float_above_each_minimum_plus_a_rejecting_sentinel() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    stats = {"a": Fraction(0.0004), "b": Fraction(0.004), "c": Fraction(0.003)}
    axis = derive_advntr_axis(ADVNTR_CUTOFF, stats, baseline=baseline, max_values=None)
    expected = sorted({0.0004, math.nextafter(0.0004, 1), 0.001, math.nextafter(0.003, 1), math.nextafter(0.004, 1)})
    assert list(axis.values) == expected
    assert axis.sentinel == 0.0004  # c == min(q) calls nobody under strict <
    assert axis.observed_count == 3
    assert axis.capped is False and axis.rejected == ()
    candidates = axis_candidates(baseline, axis)
    assert [dict(candidate.parameters) for candidate in candidates] == [
        {_CUT: value} if value != 0.001 else {} for value in expected
    ]
    assert candidates[0].policy.values[_SUP] == 3


def test_a_zero_pvalue_cannot_be_rejected_and_the_sentinel_is_the_smallest_positive_float() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    stats = {"a": Fraction(0), "b": Fraction(0.004)}
    axis = derive_advntr_axis(ADVNTR_CUTOFF, stats, baseline=baseline, max_values=None)
    assert axis.sentinel == math.nextafter(0.0, 1.0)
    assert list(axis.values) == [math.nextafter(0.0, 1.0), 0.001, math.nextafter(0.004, 1)]
    assert unrejectable_samples(ADVNTR_CUTOFF, stats) == 1
    assert unrejectable_samples(ADVNTR_MIN_SUPPORT, {"a": Fraction(0)}) == 0


def test_the_sentinel_is_omitted_when_the_baseline_already_rejects_everything() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    axis = derive_advntr_axis(ADVNTR_CUTOFF, {"a": Fraction(0.004)}, baseline=baseline, max_values=None)
    assert axis.sentinel is None
    assert list(axis.values) == [0.001, math.nextafter(0.004, 1)]
    at_anchor = derive_advntr_axis(ADVNTR_CUTOFF, {"a": Fraction(0.001)}, baseline=baseline, max_values=None)
    assert at_anchor.sentinel is None


def test_an_empty_observation_is_the_anchor_alone() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    for axis_name, anchor in ((ADVNTR_CUTOFF, 0.001), (ADVNTR_MIN_SUPPORT, 3)):
        axis = derive_advntr_axis(axis_name, {}, baseline=baseline, max_values=None)
        assert axis.values == (anchor,) and axis.sentinel is None and axis.observed_count == 0


def test_cutoff_candidates_the_decoder_refuses_are_recorded() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    axis = derive_advntr_axis(ADVNTR_CUTOFF, {"a": Fraction(1)}, baseline=baseline, max_values=None)
    assert [value for value, _ in axis.rejected] == [math.nextafter(1.0, 2.0)]
    assert list(axis.values) == [0.001]


def test_support_breakpoints_are_observed_maxima_plus_one_above_the_largest() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    axis = derive_advntr_axis(
        ADVNTR_MIN_SUPPORT, {"a": Fraction(5), "b": Fraction(2)}, baseline=baseline, max_values=None
    )
    assert axis.values == (2, 3, 5, 6) and axis.sentinel == 6
    beyond = derive_advntr_axis(ADVNTR_MIN_SUPPORT, {"a": Fraction(2)}, baseline=baseline, max_values=None)
    assert beyond.values == (2, 3) and beyond.sentinel is None
    at_top = derive_advntr_axis(ADVNTR_MIN_SUPPORT, {"a": Fraction(3)}, baseline=baseline, max_values=None)
    assert at_top.values == (3, 4) and at_top.sentinel == 4


def test_the_cap_keeps_the_extremes_and_the_anchor() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    # Statistics are exact float values (they round-trip through float as breakpoints).
    stats = {f"s{i}": Fraction(i / 1000) for i in range(2, 40)}
    stats["low"] = Fraction(0.0001)
    axis = derive_advntr_axis(ADVNTR_CUTOFF, stats, baseline=baseline, max_values=5)
    assert axis.capped and len(axis.values) == 5
    assert axis.sentinel == 0.0001
    assert 0.001 in axis.values and axis.sentinel in axis.values and max(axis.values) == math.nextafter(0.039, 1)
    assert axis.values == (
        0.0001,
        0.001,
        math.nextafter(0.019, 1),
        math.nextafter(0.029, 1),
        math.nextafter(0.039, 1),
    )
    assert axis.observed_count == 39


def test_derive_refuses_inexact_statistics() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    with pytest.raises(ValueError, match="statistics must be exact Fractions"):
        derive_advntr_axis(ADVNTR_CUTOFF, {"a": 0.004}, baseline=baseline, max_values=None)  # type: ignore[dict-item]


@pytest.mark.parametrize("value", [Fraction(1, 10000), Fraction(1, 3), Fraction(10**400)])
def test_derive_refuses_statistics_that_are_not_exact_floats(value: Fraction) -> None:
    """Breakpoints and sentinels round through float, so a statistic must be a float value exactly."""
    baseline = _policy(cutoff=0.001, support=3)
    for axis in (ADVNTR_CUTOFF, ADVNTR_MIN_SUPPORT):
        with pytest.raises(ValueError, match="statistics must be exactly representable as floats"):
            derive_advntr_axis(axis, {"a": value}, baseline=baseline, max_values=None)


# --- probe_visits: statistics read from the native replay receipt ---------------------------------


class _VisitReplayTool(_ReplayTool):
    """The replay fake, with the decision visits of each sample key injected into its receipt."""

    def __init__(self, visits: dict[str, list[object]], **kwargs: object) -> None:
        super().__init__(**kwargs)  # type: ignore[arg-type]
        self.visits = visits

    def __call__(self, argv: tuple[str, ...], **kwargs: object) -> subprocess.CompletedProcess[str]:
        completed = super().__call__(argv, **kwargs)
        if argv[-2:] == ("capabilities", "--json") or completed.returncode:
            return completed
        receipt = Path(argv[argv.index("--output") + 1]) / "replay.json"
        document = json.loads(receipt.read_bytes())
        for row in document["results"]:
            for record in row["vntrs"]:
                record["result"]["decision_visits"] = self.visits.get(row["key"], [])
        receipt.write_bytes(_json_bytes(document))
        return completed


def _grid(tmp_path: Path, visits: dict[str, list[object]], **tool: object) -> AdvntrCutoffGridResult:
    baseline = _policy()
    paths = {}
    for index, key in enumerate(sorted(visits), start=17):
        paths[key] = tmp_path / f"{key}.jsonl"
        _capture(paths[key], vntr_id=index, baseline=baseline)
    return evaluate_advntr_cutoff_grid(
        paths,
        {"baseline": baseline},
        baseline_policy_id="baseline",
        executable_path=tmp_path / "advntr",
        output=tmp_path / "result",
        runner=_VisitReplayTool(visits, **tool),
    )


def test_probe_visits_keep_only_scored_dispositions(tmp_path: Path) -> None:
    visits: dict[str, list[object]] = {
        "sample-a": [
            _visit("called", 5, 0.0004),
            _visit("cutoff", 2, 0.5),
            _visit("cutoff", 4, 0),
            _visit("called", 7, 1),
            _visit("legacy-nonfinite", 6, None),
            _visit("insufficient-read-support", 1, None, statistic=False),
            _visit("outside-boundary", 9, 0.001),
            {"disposition": "legacy-nonfinite", "plan": None, "statistic": "not-an-object"},
            _visit("called", 0, 0.1),
        ],
        "sample-b": [],
    }
    result = _grid(tmp_path, visits)
    assert probe_visits(result, "baseline") == {
        "sample-a": (
            AdvntrVisit(5, Fraction(0.0004)),
            AdvntrVisit(2, Fraction(0.5)),
            AdvntrVisit(4, Fraction(0)),
            AdvntrVisit(7, Fraction(1)),
            AdvntrVisit(0, Fraction(0.1)),
        ),
        "sample-b": (),
    }


def test_an_audit_failed_sample_has_no_visits(tmp_path: Path) -> None:
    visits: dict[str, list[object]] = {"sample-audit": [_visit("called", 5, 0.0004)], "sample-ok": []}
    result = _grid(tmp_path, visits, audit_key="sample-audit")
    assert probe_visits(result, "baseline") == {"sample-audit": None, "sample-ok": ()}


def test_probe_visits_refuse_an_absent_policy(tmp_path: Path) -> None:
    result = _grid(tmp_path, {"sample-a": []})
    with pytest.raises(ValueError, match="adVNTR probe policy missing is absent from the replay grid"):
        probe_visits(result, "missing")


@pytest.mark.parametrize(
    ("visit", "message"),
    [
        (_visit("called", 3, 0.2, statistic=False), "statistic and plan must be objects"),
        (_visit("cutoff", 3, 0.2, statistic=False), "statistic and plan must be objects"),
        (_visit("called", 5, None), "must carry a numeric p-value"),
        (_visit("cutoff", 5, None), "must carry a numeric p-value"),
        (_visit("called", 5, True), "must carry a numeric p-value"),
        (_visit("called", 5, "0.1"), "must carry a numeric p-value"),
        (_visit("called", 5, 1.5), r"must be finite and lie in \[0, 1\]"),
        (_visit("called", 5, -0.1), r"must be finite and lie in \[0, 1\]"),
        (_visit("called", True, 0.1), "read support must be a non-negative integer"),
        (_visit("called", -1, 0.1), "read support must be a non-negative integer"),
        (_visit("called", 2.0, 0.1), "read support must be a non-negative integer"),
        ({"disposition": "called", "plan": [], "statistic": {"pvalue": 0.1}}, "statistic and plan must be objects"),
        ({"disposition": "called", "plan": {"read_support": 1}, "statistic": 0.1}, "statistic and plan must be"),
        ("called", "decision visit must be an object"),
    ],
)
def test_probe_visits_refuse_a_malformed_scored_visit(tmp_path: Path, visit: object, message: str) -> None:
    result = _grid(tmp_path, {"sample-a": [visit]})
    with pytest.raises(ValueError, match=message):
        probe_visits(result, "baseline")


def test_probe_visits_refuse_a_receipt_without_a_visit_list(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import vntyper.scripts.calibration_cutoff_advntr_axes as module

    result = _grid(tmp_path, {"sample-a": []})
    monkeypatch.setattr(module, "replay_locus_document", lambda _locus: {"decision_visits": None})
    with pytest.raises(ValueError, match="adVNTR replay decision visits must be a list"):
        probe_visits(result, "baseline")


# --- the synthetic native grid (tests/unit/advntr_grid_fakes.py) ---------------------------------------


def test_the_synthetic_grid_is_readable_by_probe_visits_and_predicts_its_calls() -> None:
    sample_visits: dict[str, Sequence[AdvntrVisit] | None] = {
        "pos-a": (_v(5, 0.0004), _v(2, 1e-9)),
        "neg-c": (_v(3, 0.003),),
        "unassessable": None,
    }
    result = advntr_grid_result(
        {"probe": None, "fixed": {"pos-a": False, "neg-c": True, "unassessable": None}},
        visits={"probe": sample_visits},
        thresholds={"probe": (0.001, 3)},
    )
    assert result.baseline_policy_id == "probe"
    probe, fixed = result.policies
    assert {row.key: row.called_positive for row in probe.samples} == {
        "neg-c": False,
        "pos-a": True,
        "unassessable": None,
    }
    assert probe.execution_id == "exec-0.001-3" and fixed.execution_id == "exec-fixed"
    assert {row.key: row.called_positive for row in fixed.samples}["neg-c"] is True
    # The support-2 visit is below the probe's support 3, so it is not scored.
    assert probe_visits(result, "probe") == {
        "neg-c": (AdvntrVisit(3, Fraction(0.003)),),
        "pos-a": (AdvntrVisit(5, Fraction(0.0004)),),
        "unassessable": None,
    }
    assert probe_visits(result, "fixed") == {"neg-c": (), "pos-a": (), "unassessable": None}
    with pytest.raises(AssertionError, match="needs visits and thresholds"):
        advntr_grid_result({"probe": None})


# --- predicted_call and replay consistency ---------------------------------------------------------


def test_predicted_calls_follow_the_strict_cutoff_and_inclusive_support() -> None:
    visits = (_v(3, 0.001),)
    assert predicted_call(visits, 0.001, 3) is False  # strict <
    assert predicted_call(visits, math.nextafter(0.001, 1), 3) is True
    assert predicted_call(visits, 0.5, 4) is False  # support >= 4 fails
    assert predicted_call(None, 0.5, 1) is None
    assert predicted_call((), 0.5, 1) is False


_CONSISTENCY_VISITS: dict[str, tuple[AdvntrVisit, ...] | None] = {
    "neg-b": (_v(5, 0.004),),
    "pos-a": (_v(5, 0.0004),),
    "unassessable": None,
}


def _grid_with_calls(
    native: Mapping[str, Mapping[str, bool | None]],
) -> tuple[AdvntrCutoffGridResult, list[CutoffCandidate], dict[str, tuple[AdvntrVisit, ...] | None]]:
    """Cutoff-axis candidates replayed with the rule's own calls, except where ``native`` overrides."""
    baseline = _policy(cutoff=0.001, support=3)
    axis = derive_advntr_axis(
        ADVNTR_CUTOFF,
        sample_statistics(ADVNTR_CUTOFF, _CONSISTENCY_VISITS, baseline),
        baseline=baseline,
        max_values=None,
    )
    candidates = list(axis_candidates(baseline, axis))
    policies: dict[str, Mapping[str, bool | None] | None] = {}
    thresholds: dict[str, tuple[float, int]] = {}
    for candidate in candidates:
        cutoff = float(candidate.policy.values[_CUT])  # type: ignore[arg-type]
        thresholds[candidate.candidate_id] = (cutoff, 3)
        policies[candidate.candidate_id] = native.get(candidate.candidate_id)
        if policies[candidate.candidate_id] is not None:
            calls = {key: predicted_call(items, cutoff, 3) for key, items in _CONSISTENCY_VISITS.items()}
            policies[candidate.candidate_id] = {**calls, **native[candidate.candidate_id]}
    result = advntr_grid_result(
        policies,
        visits={candidate.candidate_id: _CONSISTENCY_VISITS for candidate in candidates},
        thresholds=thresholds,
    )
    return result, candidates, dict(_CONSISTENCY_VISITS)


def test_replay_consistency_passes_when_native_calls_match() -> None:
    result, candidates, visits = _grid_with_calls({})
    assert len(candidates) == 4  # sentinel, nextafter(0.0004), baseline, nextafter(0.004)
    assert check_replay_consistency(result, candidates, visits) == len(candidates)


def test_replay_consistency_ignores_an_unassessable_sample_that_stays_unassessable() -> None:
    result, candidates, visits = _grid_with_calls({})
    assert all(
        row.called_positive is None for policy in result.policies for row in policy.samples if row.key == "unassessable"
    )
    assert check_replay_consistency(result, candidates, visits) == 4


def test_replay_consistency_fails_and_names_the_axis_candidate_and_threshold() -> None:
    result, candidates, visits = _grid_with_calls({})
    target = candidates[1]  # the smallest cutoff that calls pos-a
    result, candidates, visits = _grid_with_calls({target.candidate_id: {"pos-a": False}})
    cutoff = math.nextafter(0.0004, 1)
    with pytest.raises(ValueError, match="replay consistency") as error:
        check_replay_consistency(result, candidates, visits)
    message = str(error.value)
    assert f"axis {ADVNTR_CUTOFF}" in message
    assert target.candidate_id in message
    assert f"cutoff={cutoff!r}, support=3" in message
    assert "on 1 samples" in message


def test_replay_consistency_labels_a_baseline_candidate_and_a_support_axis_candidate() -> None:
    baseline = _policy(cutoff=0.001, support=3)
    moved = _with(baseline, _SUP, 6)
    candidates = [
        CutoffCandidate("anchor", baseline, {}),
        CutoffCandidate("support-6", moved, {_SUP: 6}),
    ]
    visits: dict[str, tuple[AdvntrVisit, ...] | None] = {"pos-a": (_v(5, 0.0004),)}
    wrong_anchor = advntr_grid_result({"anchor": {"pos-a": False}, "support-6": {"pos-a": False}})
    with pytest.raises(ValueError, match=r"axis baseline for candidate anchor \(cutoff=0.001, support=3\)"):
        check_replay_consistency(wrong_anchor, candidates, visits)
    wrong_support = advntr_grid_result({"anchor": {"pos-a": True}, "support-6": {"pos-a": True}})
    with pytest.raises(
        ValueError, match=rf"axis {ADVNTR_MIN_SUPPORT} for candidate support-6 \(cutoff=0.001, support=6\)"
    ):
        check_replay_consistency(wrong_support, candidates, visits)


def test_replay_consistency_refuses_a_missing_candidate_and_a_different_roster() -> None:
    result, candidates, visits = _grid_with_calls({})
    missing = CutoffCandidate("not-replayed", candidates[0].policy, candidates[0].parameters)
    with pytest.raises(ValueError, match="replay consistency: candidate not-replayed was not replayed"):
        check_replay_consistency(result, [missing], visits)
    with pytest.raises(ValueError, match="replayed a different roster"):
        check_replay_consistency(result, candidates, {**visits, "extra": ()})


def test_replay_consistency_refuses_a_non_legacy_or_advntr_less_candidate_policy() -> None:
    result, candidates, visits = _grid_with_calls({})
    exact = _with(candidates[0].policy, _MODE, "exact")
    with pytest.raises(ValueError, match="legacy calibrated_calling mode; the candidate .* mode is 'exact'"):
        check_replay_consistency(result, [CutoffCandidate(candidates[0].candidate_id, exact, {})], visits)
    kestrel = _kestrel_only(candidates[0].policy)
    with pytest.raises(ValueError, match="require a candidate .* policy that includes adVNTR"):
        check_replay_consistency(result, [CutoffCandidate(candidates[0].candidate_id, kestrel, {})], visits)


# --- adVNTR baseline parity against the captures' own native decisions ------------------------------


def _parity_inputs(
    tmp_path: Path, native: Mapping[str, bool | None], replayed: Mapping[str, bool | None]
) -> tuple[AdvntrCutoffGridResult, dict[str, Path]]:
    paths = {
        key: parity_capture(tmp_path / f"{key}.jsonl", FIRST_VNTR_ID + index, native[key])
        for index, key in enumerate(sorted(native))
    }
    return advntr_grid_result({"baseline": replayed, "probe": replayed}), paths


_NATIVE = {"neg-b": False, "pos-a": True, "unassessable": None}


def test_baseline_parity_compares_the_anchor_replay_with_the_native_capture_calls(tmp_path: Path) -> None:
    result, paths = _parity_inputs(tmp_path, _NATIVE, _NATIVE)
    assert advntr_baseline_parity(result, "baseline", paths) == {"proven": True, "sample_count": 3, "mismatches": []}


@pytest.mark.parametrize("flip", [{"pos-a": False}, {"neg-b": True}, {"unassessable": False}, {"pos-a": None}])
def test_baseline_parity_fails_on_one_flipped_sample(tmp_path: Path, flip: dict[str, bool | None]) -> None:
    result, paths = _parity_inputs(tmp_path, _NATIVE, {**_NATIVE, **flip})
    with pytest.raises(ValueError, match="baseline parity failed: .* for 1 samples"):
        advntr_baseline_parity(result, "baseline", paths)


def test_baseline_parity_refuses_an_absent_anchor_and_a_different_roster(tmp_path: Path) -> None:
    result, paths = _parity_inputs(tmp_path, _NATIVE, _NATIVE)
    with pytest.raises(ValueError, match="baseline parity: anchor missing was not replayed"):
        advntr_baseline_parity(result, "missing", paths)
    with pytest.raises(ValueError, match="disagree about the roster"):
        advntr_baseline_parity(result, "baseline", {key: paths[key] for key in ("neg-b", "pos-a")})


def test_baseline_parity_refuses_a_capture_of_a_different_locus(tmp_path: Path) -> None:
    result, paths = _parity_inputs(tmp_path, _NATIVE, _NATIVE)
    parity_capture(paths["pos-a"], 99, True)
    with pytest.raises(ValueError, match="capture for pos-a covers different loci than its replay"):
        advntr_baseline_parity(result, "baseline", paths)


@pytest.mark.parametrize("locus", [None, {"vntr_id": 0}, {"vntr_id": True}, {"vntr_id": "17"}])
def test_baseline_parity_refuses_a_capture_without_a_positive_vntr_id(tmp_path: Path, locus: object) -> None:
    result, paths = _parity_inputs(tmp_path, _NATIVE, _NATIVE)
    document = json.loads(paths["pos-a"].read_bytes())
    document["locus"] = locus
    paths["pos-a"].write_bytes(_json_bytes(document))
    with pytest.raises(ValueError, match="lacks a positive VNTR identifier"):
        advntr_baseline_parity(result, "baseline", paths)
