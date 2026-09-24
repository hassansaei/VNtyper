"""``vntyper calibrate optimize`` over the adVNTR axes, against hand-computed oracles (#269).

The native adVNTR grid is replaced by the synthetic fake in ``cutoff_optimize_fakes``; its
calls follow the legacy rule over :data:`~tests.unit.cutoff_optimize_fakes.PROBE_VISITS`,
while adVNTR baseline parity runs for real against one-record capture files.
"""

from __future__ import annotations

import html
from collections.abc import Callable, Mapping
from dataclasses import replace
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from tests.unit.advntr_grid_fakes import synthetic_capabilities
from tests.unit.cutoff_optimize_fakes import (
    ADV_CUT,
    PROBE_VISITS,
    STANDARD_COHORT,
    advntr_baseline_policy,
    namespace,
    run_advntr,
    run_optimize,
    up,
    visit,
    write_manifests,
)
from vntyper.modules.advntr.advntr_calibration_policy import advntr_capabilities_document
from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_cutoff_advntr import advntr_signature
from vntyper.scripts.calibration_cutoff_advntr_axes import AdvntrVisit, derive_advntr_axis
from vntyper.scripts.calibration_cutoff_axes import ADVNTR_CUTOFF, ADVNTR_MIN_SUPPORT, DEPTH_FLOOR_LINKED

pytestmark = pytest.mark.unit


#: The hand-computed adVNTR cutoff oracle over the scored roster (``alpha-repeat`` is a
#: dropped duplicate): value -> (TP, FN, TN, FP). A sample is called when ``q < cutoff``.
ADVNTR_ORACLE: dict[float, tuple[int, int, int, int]] = {
    0.0004: (0, 3, 2, 0),  # the sentinel min(q): calls nobody
    up(0.0004): (1, 2, 2, 0),  # alpha
    0.001: (1, 2, 2, 0),  # the baseline: alpha
    up(0.004): (2, 1, 2, 0),  # alpha, bravo
    up(0.006): (2, 1, 1, 1),  # + delta
    up(0.02): (2, 1, 0, 2),  # + charlie
    up(0.2): (2, 1, 0, 2),  # + foxtrot (unknown truth)
}


def test_the_advntr_arm_is_replayed_natively_and_never_approximated(tmp_path: Path) -> None:
    """adVNTR statistics come from the installed evaluator; this only mocks the seam."""
    seen: list[dict[str, Any]] = []
    rescued = {**PROBE_VISITS, "specimen-echo": (visit(5, 0.0001),)}  # adVNTR calls echo at baseline
    _, document, _ = run_advntr(
        tmp_path, seen=seen, visits=rescued, caller="both", axes=[DEPTH_FLOOR_LINKED], min_specificity=1.0
    )

    assert document["caller"] == "both"
    assert len(seen) == 1  # no adVNTR axis, so no probe grid
    assert seen[0]["executable_path"] == tmp_path / "advntr"
    assert seen[0]["output"].name == "advntr"
    assert set(seen[0]["captures"]) == set(PROBE_VISITS)
    advntr = document["provenance"]["advntr"]
    assert set(advntr) == {"sha256", "probe_sha256", "tool_identity", "probe_seconds", "main_seconds"}
    assert advntr["sha256"] == "d" * 64
    assert advntr["probe_sha256"] is None and advntr["probe_seconds"] is None
    assert advntr["tool_identity"] == advntr_capabilities_document(synthetic_capabilities())
    assert isinstance(advntr["main_seconds"], float) and advntr["main_seconds"] >= 0.0
    best = max(row["counts"]["true_positives"] for row in document["cutoffs"])
    assert best == 3


def test_caller_both_with_only_kestrel_axes_holds_advntr_at_baseline(tmp_path: Path) -> None:
    """With only Kestrel axes, the adVNTR arm of ``--caller both`` is replayed at its baseline policy."""
    from vntyper.scripts.calibration_cutoff_report import render_cutoff_report_html

    _, document, _ = run_advntr(tmp_path, caller="both", axes=[DEPTH_FLOOR_LINKED], min_specificity=1.0)
    scope = document["search_scope"]

    assert scope["searched_callers"] == ["kestrel"]
    assert scope["searched_axes"] == [DEPTH_FLOOR_LINKED]
    assert scope["advntr_policy"] == "held-at-baseline"
    assert scope["advntr_distinct_executions"] == 1
    assert scope["advntr_probe_executions"] is None
    assert scope["note"] == (
        "Only Kestrel axes were requested. The adVNTR arm was replayed at its baseline policy for every candidate."
    )
    assert document["baseline_parity"]["advntr"]["proven"] is True
    assert document["replay_consistency"] is None
    page = render_cutoff_report_html(document)
    assert "adVNTR arm was replayed at its baseline policy" in page


def test_a_kestrel_run_records_that_advntr_was_not_evaluated(tmp_path: Path) -> None:
    """The scope section is present on every run, not only on the adVNTR ones."""
    _, document, _ = run_optimize(tmp_path)

    assert document["search_scope"]["searched_callers"] == ["kestrel"]
    assert document["search_scope"]["advntr_policy"] == "not-evaluated"
    assert document["search_scope"]["advntr_distinct_executions"] is None
    assert document["search_scope"]["advntr_probe_executions"] is None
    assert document["provenance"]["advntr"] is None
    assert document["baseline_parity"]["advntr"] is None
    assert document["replay_consistency"] is None
    assert all(curve["status"] == "available" for curve in document["curves"])


def test_an_unexpected_extra_advntr_execution_aborts_the_run(tmp_path: Path) -> None:
    """One adVNTR policy replayed by two executions is refused, never assumed away."""

    def extra(result: Any, policies: Mapping[str, CallerPolicyValues]) -> Any:
        rows = list(result.policies)
        kestrel_only = next(index for index, row in enumerate(rows) if row.policy_id != result.baseline_policy_id)
        rows[kestrel_only] = replace(rows[kestrel_only], execution_id="exec-unexpected")
        return replace(result, policies=tuple(rows))

    with pytest.raises(ValueError, match=r"one distinct adVNTR policy was replayed by 2 executions"):
        run_advntr(tmp_path, override=extra, caller="both", axes=[DEPTH_FLOOR_LINKED])


def _main_grid_only(edit: Callable[[list[Any], Mapping[str, CallerPolicyValues]], None]) -> Any:
    """An override that edits the candidate grid's rows and leaves the probe grid intact."""

    def override(result: Any, policies: Mapping[str, CallerPolicyValues]) -> Any:
        if "baseline" in policies:
            return result
        rows = list(result.policies)
        edit(rows, policies)
        return replace(result, policies=tuple(rows))

    return override


def _set_execution(rows: list[Any], policy_id: str, execution_id: str) -> None:
    index = next(index for index, row in enumerate(rows) if row.policy_id == policy_id)
    rows[index] = replace(rows[index], execution_id=execution_id)


def _advntr_axis_ids(policies: Mapping[str, CallerPolicyValues]) -> list[str]:
    """The adVNTR-axis candidates that move the cutoff, in id order."""
    return sorted(pid for pid, policy in policies.items() if policy.values[ADV_CUT] != 0.001)


def test_two_advntr_policies_sharing_one_execution_abort_the_run(tmp_path: Path) -> None:
    """Distinct adVNTR policies must not be served by one native execution."""

    def merge(rows: list[Any], policies: Mapping[str, CallerPolicyValues]) -> None:
        first, second = _advntr_axis_ids(policies)[:2]
        shared = next(row.execution_id for row in rows if row.policy_id == first)
        _set_execution(rows, second, shared)

    with pytest.raises(ValueError, match=r"execution exec-\S+ served 2 distinct adVNTR policies"):
        run_advntr(tmp_path, override=_main_grid_only(merge), caller="advntr", min_specificity=1.0)


def test_a_split_and_a_merge_that_keep_the_execution_count_still_abort_the_run(tmp_path: Path) -> None:
    """Equal counts are not enough: the guard binds every signature to exactly one execution.

    Under ``--caller both`` the Kestrel candidates share the baseline adVNTR signature. One of
    them gets an execution of its own (a split), and two adVNTR-axis candidates share one
    execution (a merge), so the number of distinct executions still equals the number of
    distinct signatures.
    """

    def split_and_merge(rows: list[Any], policies: Mapping[str, CallerPolicyValues]) -> None:
        kestrel = next(pid for pid in sorted(policies) if pid.startswith(DEPTH_FLOOR_LINKED))
        _set_execution(rows, kestrel, "exec-split")
        first, second = _advntr_axis_ids(policies)[:2]
        _set_execution(rows, second, next(row.execution_id for row in rows if row.policy_id == first))
        assert len({row.execution_id for row in rows}) == len({advntr_signature(p) for p in policies.values()})

    with pytest.raises(ValueError, match="adVNTR policy"):
        run_advntr(tmp_path, override=_main_grid_only(split_and_merge), caller="both", min_specificity=1.0)


def test_caller_advntr_searches_the_advntr_cutoff_axis_by_native_replay(tmp_path: Path) -> None:
    """The adVNTR cutoff axis is derived from probe replays and scored against the hand oracle.

    Outer folds (3 folds, seed 20260915): fold 0 holds out bravo and foxtrot, fold 1 alpha and
    delta, fold 2 charlie and echo. No fold derives a value the full-data axis lacks: fold 0
    trains on q {0.0004, 0.006, 0.02}, fold 1 on {0.004, 0.02, 0.2} (the baseline 0.001 already
    rejects all, so no sentinel), fold 2 on {0.0004, 0.004, 0.006, 0.2}.
    """
    seen: list[dict[str, Any]] = []
    successful, document, _ = run_advntr(tmp_path, seen=seen, caller="advntr", min_specificity=1.0)

    assert successful is True
    rows = {row["value"]: row for row in document["cutoffs"] if row["policy_id"] != "baseline"}
    assert sorted(rows) == sorted(ADVNTR_ORACLE)
    for value, (tp, fn, tn, fp) in ADVNTR_ORACLE.items():
        counts = rows[value]["counts"]
        assert (counts["true_positives"], counts["false_negatives"]) == (tp, fn), value
        assert (counts["true_negatives"], counts["false_positives"]) == (tn, fp), value
    assert document["selection"]["axis"] == ADVNTR_CUTOFF
    assert document["selection"]["value"] == up(0.004)
    (axis,) = document["axes"]
    assert axis["axis"] == ADVNTR_CUTOFF
    assert axis["caller"] == "advntr"
    assert axis["comparator"] == "<"
    assert axis["breakpoint_completeness"] == "complete"
    assert axis["fold_only_values"] == 0
    assert axis["unrejectable_samples"] == 0
    scope = document["search_scope"]
    assert scope["advntr_policy"] == "searched"
    assert scope["searched_callers"] == ["advntr"]
    assert scope["advntr_probe_executions"] == 2  # the baseline and one cutoff probe
    assert document["replay_consistency"] == {"checked_candidates": 7, "mismatches": []}
    assert document["baseline_parity"]["advntr"]["proven"] is True
    assert document["evaluation"]["fold_admissibility"] == "training-derived-inventories"
    changed = {row["pointer"]: row for row in document["old_versus_derived"] if row["changed"]}
    assert set(changed) == {ADV_CUT}
    assert changed[ADV_CUT]["baseline_value"] == 0.001
    assert changed[ADV_CUT]["derived_value"] == up(0.004)
    assert document["profile"]["round_trip_matches_selected_policy"] is True

    probe, main = seen
    assert sorted(probe["policies"]) == ["baseline", f"probe-{ADVNTR_CUTOFF}"]
    assert probe["output"].name == "advntr-probe" and main["output"].name == "advntr"
    # X4: one execution per distinct adVNTR signature passed to the main grid, the anchor once.
    assert len({advntr_signature(policy) for policy in main["policies"].values()}) == 7
    assert document["search_scope"]["advntr_distinct_executions"] == 7
    advntr = document["provenance"]["advntr"]
    assert set(advntr) == {"sha256", "probe_sha256", "tool_identity", "probe_seconds", "main_seconds"}
    assert advntr["tool_identity"] == advntr_capabilities_document(synthetic_capabilities())
    assert advntr["probe_sha256"] == "d" * 64
    assert isinstance(advntr["probe_seconds"], float) and isinstance(advntr["main_seconds"], float)


def test_caller_both_searches_kestrel_and_advntr_axes_on_the_union(tmp_path: Path) -> None:
    """``--caller both`` defaults to one Kestrel and one adVNTR axis, each holding the other caller at baseline."""
    seen: list[dict[str, Any]] = []
    _, document, _ = run_advntr(tmp_path, seen=seen, caller="both", min_specificity=1.0)

    assert [axis["axis"] for axis in document["axes"]] == [DEPTH_FLOOR_LINKED, ADVNTR_CUTOFF]
    assert [axis["caller"] for axis in document["axes"]] == ["kestrel", "advntr"]
    _, main = seen
    signatures = {advntr_signature(policy) for policy in main["policies"].values()}
    # Every Kestrel candidate shares the baseline adVNTR signature; the adVNTR axis adds 6 more.
    assert len(signatures) == 7
    assert document["search_scope"]["advntr_distinct_executions"] == 1 + 6
    assert document["search_scope"]["searched_callers"] == ["advntr", "kestrel"]
    assert document["search_scope"]["advntr_policy"] == "searched"
    # Only adVNTR-axis candidates are checked: the seven values of ADVNTR_ORACLE.
    assert document["replay_consistency"] == {"checked_candidates": 7, "mismatches": []}
    assert document["baseline_parity"]["advntr"]["proven"] is True
    assert [curve["status"] for curve in document["curves"]] == ["available", "available"]
    advntr_values = {row["value"] for row in document["cutoffs"] if row["policy_id"].startswith(ADVNTR_CUTOFF)}
    assert advntr_values == set(ADVNTR_ORACLE)


def test_a_union_curve_whose_no_call_set_changes_is_published_unavailable(tmp_path: Path) -> None:
    """Echo is a Kestrel no-call everywhere; adVNTR calls it at every cutoff except the sentinel.

    On the either-caller union echo is therefore called on every adVNTR candidate but the
    sentinel 0.0001, where it stays a no-call, so the adVNTR axis has no fixed-denominator
    curve (spec 14.6). The Kestrel axis holds adVNTR at baseline, which calls echo, so its
    curve is unaffected.
    """
    from vntyper.scripts.calibration_cutoff_document import _CURVE_UNAVAILABLE

    rescued = {**PROBE_VISITS, "specimen-echo": (visit(5, 0.0001),)}
    successful, document, output = run_advntr(tmp_path, visits=rescued, caller="both", min_specificity=1.0)

    assert successful is True
    kestrel, advntr = document["curves"]
    assert kestrel["axis"] == DEPTH_FLOOR_LINKED and kestrel["status"] == "available"
    assert advntr == {"axis": ADVNTR_CUTOFF, "status": "unavailable", "reason": _CURVE_UNAVAILABLE}
    assert 0.0001 in {row["value"] for row in document["cutoffs"]}
    assert document["selection"]["axis"] == ADVNTR_CUTOFF  # the plateau needs outcome vectors, not the curve
    assert document["selection"]["plateau"]["selected_value"] == up(0.004)
    assert document["boundary_support"] == {"status": "unavailable", "reason": _CURVE_UNAVAILABLE, "warnings": []}
    curve_axes = {line.split("\t")[0] for line in (output / "roc-pr-curves.tsv").read_text().splitlines()[1:]}
    assert curve_axes == {DEPTH_FLOOR_LINKED}
    assert html.escape(_CURVE_UNAVAILABLE) in (output / "report.html").read_text(encoding="utf-8")


def test_advntr_replay_inconsistency_aborts_the_run(tmp_path: Path) -> None:
    """A native call the probe statistics do not predict aborts the run, naming the candidate."""

    def contradict(result: Any, policies: Mapping[str, CallerPolicyValues]) -> Any:
        if "baseline" in policies:
            return result  # the probe grid is left intact
        target = next(pid for pid, policy in policies.items() if policy.values[ADV_CUT] == up(0.004))
        rows = []
        for row in result.policies:
            if row.policy_id == target:
                samples = tuple(
                    replace(sample, called_positive=False) if sample.key == "specimen-bravo" else sample
                    for sample in row.samples
                )
                row = replace(row, samples=samples)
            rows.append(row)
        return replace(result, policies=tuple(rows))

    with pytest.raises(ValueError, match=rf"replay consistency failed on axis {ADVNTR_CUTOFF}"):
        run_advntr(tmp_path, override=contradict, caller="advntr")


def test_advntr_probe_and_candidate_grids_must_share_their_evidence(tmp_path: Path) -> None:
    """A candidate grid replayed by a different adVNTR build is not bound to the probe's statistics."""

    def rebuilt(result: Any, policies: Mapping[str, CallerPolicyValues]) -> Any:
        return result if "baseline" in policies else replace(result, capabilities=synthetic_capabilities("2.4.1"))

    with pytest.raises(ValueError, match="not bound to the same evidence"):
        run_advntr(tmp_path, override=rebuilt, caller="advntr")


def test_advntr_fold_admissibility_uses_training_derived_inventories(tmp_path: Path) -> None:
    """Each fold may select only cutoffs its own training samples derived, plus the baseline.

    Hand derivation (folds as in the oracle test; ``q`` from :data:`PROBE_VISITS`):

    * fold 0 trains on alpha, charlie, delta, echo: q {0.0004, 0.02, 0.006} gives
      {0.0004 (sentinel), up(0.0004), 0.001, up(0.006), up(0.02)}: 5 candidates + baseline = 6.
    * fold 1 trains on bravo, charlie, echo, foxtrot: q {0.004, 0.02, 0.2}; the baseline 0.001
      already rejects all, so no sentinel: {0.001, up(0.004), up(0.02), up(0.2)}: 4 + 1 = 5.
    * fold 2 trains on alpha, bravo, delta, foxtrot: q {0.0004, 0.004, 0.006, 0.2} gives
      {0.0004, up(0.0004), 0.001, up(0.004), up(0.006), up(0.2)}: 6 + 1 = 7.
    """
    _, document, _ = run_advntr(tmp_path, caller="advntr", min_specificity=1.0)
    evaluation = document["evaluation"]
    ids = {row["value"]: row["policy_id"] for row in document["cutoffs"] if row["policy_id"] != "baseline"}
    baseline = advntr_baseline_policy()

    assert evaluation["fold_admissibility"] == "training-derived-inventories"
    assert [fold["admissible_candidates"] for fold in evaluation["folds"]] == [6, 5, 7]
    for fold in evaluation["folds"]:
        training = set(fold["training_keys"])
        statistics = {
            key: min(visit.pvalue for visit in items)
            for key, items in PROBE_VISITS.items()
            if key in training and items
        }
        derived = derive_advntr_axis(ADVNTR_CUTOFF, statistics, baseline=baseline, max_values=None)
        admissible = {"baseline", *(ids[float(value)] for value in derived.values)}
        assert fold["admissible_candidates"] == len(admissible)
        assert fold["used_policy"] in admissible


#: Mixed read supports for a search of both adVNTR axes (baseline cutoff 0.001, support 3).
#:
#: ====================  =========================  ==============  ==============  ========
#: Sample                visits (support, p)        q (support>=3)  m (p < 0.001)   truth
#: ====================  =========================  ==============  ==============  ========
#: ``specimen-alpha``    (5, 0.0004)                0.0004          5               positive
#: ``specimen-bravo``    (2, 0.0002), (5, 0.004)    0.004           2               positive
#: ``specimen-charlie``  (4, 0.0008)                0.0008          4               negative
#: ``specimen-delta``    (5, 0.006)                 0.006           --              negative
#: ``specimen-echo``     --                         --              --              positive
#: ``specimen-foxtrot``  (1, 0.0001), (3, 0.2)      0.2             1               unknown
#: ====================  =========================  ==============  ==============  ========
#:
#: Bravo's only visit below the baseline cutoff has support 2, so the baseline support 3
#: misses it and only a lowered support rescues it.
MIXED_VISITS: dict[str, tuple[AdvntrVisit, ...] | None] = {
    "specimen-alpha": (visit(5, 0.0004),),
    "specimen-bravo": (visit(2, 0.0002), visit(5, 0.004)),
    "specimen-charlie": (visit(4, 0.0008),),
    "specimen-delta": (visit(5, 0.006),),
    "specimen-echo": (),
    "specimen-foxtrot": (visit(1, 0.0001), visit(3, 0.2)),
}

#: Support axis, value -> (TP, FN, TN, FP); called when a visit has support >= value and
#: p < 0.001. Observed m {1, 2, 4, 5}, the anchor 3 and the sentinel max(m) + 1 = 6.
MIXED_SUPPORT_ORACLE: dict[int, tuple[int, int, int, int]] = {
    1: (2, 1, 1, 1),  # alpha, bravo; charlie (FP); foxtrot is unknown truth
    2: (2, 1, 1, 1),
    3: (1, 2, 1, 1),  # the baseline: bravo is lost
    4: (1, 2, 1, 1),
    5: (1, 2, 2, 0),  # charlie (support 4) is rejected
    6: (0, 3, 2, 0),  # the sentinel: calls nobody
}

#: Cutoff axis, value -> (TP, FN, TN, FP); called when a visit has support >= 3 and
#: p < value. Full-data q {0.0004, 0.0008, 0.004, 0.006, 0.2}: next-up values, the anchor
#: 0.001 and the sentinel min(q) = 0.0004. Fold 1 trains on bravo, charlie and foxtrot, so
#: its own sentinel min(q) = 0.0008 is a value only a fold inventory derives.
MIXED_CUTOFF_ORACLE: dict[float, tuple[int, int, int, int]] = {
    0.0004: (0, 3, 2, 0),
    up(0.0004): (1, 2, 2, 0),
    0.0008: (1, 2, 2, 0),  # fold 1's sentinel
    up(0.0008): (1, 2, 1, 1),
    0.001: (1, 2, 1, 1),
    up(0.004): (2, 1, 1, 1),
    up(0.006): (2, 1, 0, 2),
    up(0.2): (2, 1, 0, 2),
}


def test_both_advntr_axes_with_mixed_supports_match_the_hand_oracle(tmp_path: Path) -> None:
    """Searching the cutoff and support axes together replays and checks every candidate once.

    Replay consistency checks 8 cutoff + 6 support = 14 candidates. The two anchors carry
    the same baseline adVNTR policy, so the candidate grid executes 8 + 6 - 1 = 13 distinct
    policies; the probe grid executes the baseline and one probe per axis, 3 in all.
    """
    seen: list[dict[str, Any]] = []
    _, document, _ = run_advntr(
        tmp_path,
        seen=seen,
        visits=MIXED_VISITS,
        caller="advntr",
        axes=[ADVNTR_CUTOFF, ADVNTR_MIN_SUPPORT],
        min_specificity=1.0,
    )

    for axis, oracle in ((ADVNTR_CUTOFF, MIXED_CUTOFF_ORACLE), (ADVNTR_MIN_SUPPORT, MIXED_SUPPORT_ORACLE)):
        rows = {row["value"]: row["counts"] for row in document["cutoffs"] if row["policy_id"].startswith(axis)}
        assert sorted(rows) == sorted(oracle), axis
        for value, (tp, fn, tn, fp) in oracle.items():
            counts = rows[value]
            assert (counts["true_positives"], counts["false_negatives"]) == (tp, fn), (axis, value)
            assert (counts["true_negatives"], counts["false_positives"]) == (tn, fp), (axis, value)
    cutoff, support = document["axes"]
    assert support["axis"] == ADVNTR_MIN_SUPPORT and support["endpoint_sentinel"] == 6
    assert support["fold_only_values"] == 0
    assert cutoff["endpoint_sentinel"] == 0.0004 and cutoff["fold_only_values"] == 1
    assert document["replay_consistency"] == {"checked_candidates": 14, "mismatches": []}
    scope = document["search_scope"]
    assert scope["advntr_distinct_executions"] == 13
    assert scope["advntr_probe_executions"] == 3
    probe, main = seen
    assert sorted(probe["policies"]) == ["baseline", f"probe-{ADVNTR_CUTOFF}", f"probe-{ADVNTR_MIN_SUPPORT}"]
    assert len({advntr_signature(policy) for policy in main["policies"].values()}) == 13


def test_optimize_breaks_selection_ties_by_policy_content(tmp_path: Path) -> None:
    """Candidate IDs follow global inventory rank, so optimize hands selection the policy digests."""
    from vntyper.scripts import calibration_cutoff_optimize as module

    real = module.evaluate_cutoff_arms
    seen: list[dict[str, Any]] = []

    def record(*call_args: Any, **call_kwargs: Any) -> Any:
        seen.append(call_kwargs)
        return real(*call_args, **call_kwargs)

    with patch.object(module, "evaluate_cutoff_arms", record):
        _, document, _ = run_advntr(tmp_path, caller="both", min_specificity=1.0)

    digests = {row["policy_id"]: row["policy_sha256"] for row in document["cutoffs"] if row["policy_id"] != "baseline"}
    assert seen[0]["tie_keys"] == digests


def test_a_multi_axis_search_holding_only_baseline_anchors_publishes_no_joint_points(tmp_path: Path) -> None:
    """Every visit has p = 1, so neither adVNTR axis observes a breakpoint beyond its anchor.

    The cutoff axis sees q = 1 everywhere: its next-up value exceeds 1 and is refused, and the
    baseline 0.001 already rejects every sample, so there is no sentinel. The support axis
    sees no visit below the baseline cutoff at all. Both axes hold their anchor alone, so
    there is no non-baseline point to tabulate.
    """
    flat = {key: (visit(5, 1.0),) for key in PROBE_VISITS}
    successful, document, output = run_advntr(
        tmp_path, visits=flat, caller="advntr", axes=[ADVNTR_CUTOFF, ADVNTR_MIN_SUPPORT], min_specificity=1.0
    )

    assert [len(axis["values"]) for axis in document["axes"]] == [1, 1]
    assert document["joint_points"] is None
    assert not (output / "joint-points.tsv").exists()
    assert (output / "report.html").is_file()
    assert successful is True


def test_captures_the_research_runtime_cannot_reproduce_are_refused_before_probing(tmp_path: Path) -> None:
    """A profile exported from these captures would run adVNTR under different capture semantics."""
    seen: list[dict[str, Any]] = []
    with pytest.raises(
        ValueError, match=r"capture settings the research runtime cannot reproduce: use_reference_alignment"
    ):
        run_advntr(tmp_path, seen=seen, caller="advntr", capture_parameters={"use_reference_alignment": False})

    assert seen == []  # neither the probe grid nor the candidate grid ran
    assert not (tmp_path / "derived").exists()


def test_the_capture_thread_count_is_not_part_of_the_research_contract(tmp_path: Path) -> None:
    successful, _, _ = run_advntr(tmp_path, caller="advntr", min_specificity=1.0, capture_parameters={"threads": 16})

    assert successful is True


def test_a_kestrel_search_over_an_exact_advntr_baseline_exports_nothing(tmp_path: Path) -> None:
    """The research runtime refuses exact adVNTR, so the export is refused before it is written.

    ``--caller kestrel`` never replays adVNTR, yet the exported profile carries the baseline's
    adVNTR pointers, exact mode included; the atomic output then publishes nothing at all.
    """
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = write_manifests(tmp_path, STANDARD_COHORT, advntr_baseline="exact")
    args = namespace(cohort_path, captures_path, caller="kestrel")

    with pytest.raises(ValueError, match="cutoff optimize cannot export .* adVNTR mode 'exact'"):
        atomic_output(tmp_path / "derived", lambda staging: run_cutoff_optimization(args, staging))
    assert not (tmp_path / "derived").exists()


def test_the_exported_advntr_profile_is_rendered_through_the_research_argv_builder(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_cutoff_optimize as module

    real = module.research_policy_argv
    rendered: list[tuple[str, ...]] = []

    def record(*call_args: Any, **call_kwargs: Any) -> tuple[str, ...]:
        rendered.append(real(*call_args, **call_kwargs))
        return rendered[-1]

    with patch.object(module, "research_policy_argv", record):
        _, document, _ = run_advntr(tmp_path, caller="advntr", min_specificity=1.0)

    (argv,) = rendered
    assert argv[:4] == ("-t", "1", "--frameshift-pvalue-cutoff", repr(document["selection"]["value"]))


def test_expected_decoder_refusals_are_reported_at_info_and_never_as_errors(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Screening a breakpoint the decoder refuses is routine, and the axis document records it.

    Every visit at p = 1 makes the cutoff axis test ``nextafter(1, +inf)``, which the decoder
    refuses in the full-data and in every fold inventory, and the Kestrel probe ladders
    include rungs the decoder refuses. None of that is an error; each refusal is logged once.
    """
    import logging

    from vntyper.scripts.calibration_cutoff_optimize import AXIS_PROBE

    caplog.set_level(logging.INFO)
    flat = {key: (visit(5, 1.0),) for key in PROBE_VISITS}
    _, document, _ = run_advntr(
        tmp_path / "advntr", visits=flat, caller="advntr", axes=[ADVNTR_CUTOFF], min_specificity=1.0
    )
    run_optimize(tmp_path / "kestrel", axes=list(AXIS_PROBE), objective="youden-j", min_specificity=None)

    assert document["axes"][0]["rejected"], "the test needs a refused breakpoint"
    assert [record.getMessage() for record in caplog.records if record.levelno >= logging.ERROR] == []
    refusals = [record.getMessage() for record in caplog.records if " refused: " in record.getMessage()]
    assert refusals and len(refusals) == len(set(refusals))
    assert any(ADVNTR_CUTOFF in message for message in refusals)
