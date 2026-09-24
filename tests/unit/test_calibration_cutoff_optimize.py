"""``vntyper calibrate optimize`` end to end, against a hand-computed oracle.

Every capture here is synthetic and every sample name is invented, so nothing in this
file can carry cohort identity. The cohort is small enough that the confusion matrix at
each tested cutoff is written out by hand in :data:`ORACLE` and compared against what the
command reports -- if the command ever starts computing its own truth, this disagrees.

The depth scores are chosen so the linked depth axis has a documented shape:

====================  ===========  =========================================
Sample                Depth_Score  Called at a floor of ...
====================  ===========  =========================================
``specimen-alpha``    0.014        every tested floor up to and including 0.014
``specimen-bravo``    0.004        0.004 and below -- the rescue
``specimen-delta``    0.001        0.001 and below -- the false positive
``specimen-charlie``  0.0004       0.0004 only -- the second false positive
====================  ===========  =========================================
"""

from __future__ import annotations

import argparse
import json
import math
import stat
from collections.abc import Callable, Mapping
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from types import MappingProxyType
from typing import Any, cast
from unittest.mock import patch

import pandas as pd
import pytest

from tests.builders import kestrel_config
from tests.unit.advntr_grid_fakes import advntr_grid_result, parity_capture, synthetic_capabilities
from tests.unit.test_calibration_cutoff_kestrel import _native_negative, _native_production
from tests.unit.test_calibration_kestrel_replay import _candidate, _capture, _policy, _raw
from vntyper.modules.advntr.advntr_calibration_policy import advntr_capabilities_document
from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_caller_policy import (
    ADVNTR_CALLER_POLICY_POINTERS,
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_cutoff_advntr_axes import AdvntrVisit, derive_advntr_axis, predicted_call
from vntyper.scripts.calibration_cutoff_axes import (
    ADVNTR_CUTOFF,
    ADVNTR_MIN_SUPPORT,
    DEPTH_FLOOR_LINKED,
    GG_GATE_INDEPENDENT,
)
from vntyper.scripts.calibration_kestrel_capture import (
    decode_kestrel_capture,
    kestrel_capture_document,
)
from vntyper.scripts.calibration_kestrel_replay import replay_kestrel_capture
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

pytestmark = pytest.mark.unit

_FLOOR = "/components/kestrel/confidence_assignment/reporting_floor"
_LOW = "/components/kestrel/confidence_assignment/depth_score_thresholds/low"
_GG = "/components/kestrel/alt_filtering/gg_depth_score_threshold"

#: ``(alternate depth, active-region depth)`` pairs whose ratio is the Depth_Score.
DEPTHS: dict[str, tuple[int, int]] = {
    "0.014": (7, 500),
    "0.008": (4, 500),
    "0.006": (3, 500),
    "0.004": (2, 500),
    "0.001": (5, 5000),
    "0.0004": (2, 5000),
}

#: The endpoint sentinel just above the largest observed Depth_Score: the one tested floor
#: that rejects every row, since production floors pass a score equal to the floor.
SENTINEL = math.nextafter(0.014, math.inf)

#: The hand-computed confusion matrix of the standard cohort at every tested floor.
#: Three positives (alpha, bravo, echo), two negatives (charlie, delta), one unknown
#: (foxtrot); ``echo`` has no candidates at all and is therefore a positive no-call at
#: every floor, which keeps the positive denominator at three throughout.
ORACLE: dict[str, dict[str, int]] = {
    "0.0004": {"true_positives": 2, "false_negatives": 0, "false_positives": 2, "true_negatives": 0},
    "0.001": {"true_positives": 2, "false_negatives": 0, "false_positives": 1, "true_negatives": 1},
    "0.004": {"true_positives": 2, "false_negatives": 0, "false_positives": 0, "true_negatives": 2},
    "0.00469": {"true_positives": 1, "false_negatives": 1, "false_positives": 0, "true_negatives": 2},
    "0.014": {"true_positives": 1, "false_negatives": 1, "false_positives": 0, "true_negatives": 2},
    repr(SENTINEL): {"true_positives": 0, "false_negatives": 2, "false_positives": 0, "true_negatives": 2},
}


def _frame(scores: tuple[str, ...]) -> pd.DataFrame:
    """Build a raw Kestrel frame with one candidate row per requested Depth_Score."""
    if not scores:
        return _raw().iloc[0:0]
    rows = []
    for ordinal, score in enumerate(scores):
        alternate, region = DEPTHS[score]
        row = _raw(depth_alt=alternate, depth_region=region)
        row.loc[0, "POS"] = 67 + ordinal
        rows.append(row)
    return pd.concat(rows, ignore_index=True)


#: The adVNTR half of an adVNTR-bearing baseline: legacy calling at cutoff 0.001, support 3.
_ADVNTR_BASELINE: dict[str, object] = {
    "/components/advntr/calibrated_calling/mode": "legacy",
    "/components/advntr/calibrated_calling/cutoff": 0.001,
    "/components/advntr/calibrated_calling/minimum_read_support": 3,
    "/components/advntr/calibrated_calling/rare_unit_fraction": None,
    "/components/advntr/calibrated_calling/adapter_filter": False,
    "/components/advntr/calibrated_calling/minimum_read_match_ratio": 0.6,
    "/components/advntr/calibrated_calling/prune_reverse": False,
}
_ADV_CUT = "/components/advntr/calibrated_calling/cutoff"
_ADV_SUP = "/components/advntr/calibrated_calling/minimum_read_support"


def _advntr_baseline_policy() -> CallerPolicyValues:
    """The standard Kestrel baseline plus the seven adVNTR legacy pointers."""
    kestrel = _policy(kestrel_config())
    document = caller_policy_values_document(kestrel)
    return decode_caller_policy_values(
        {**document, "required_callers": ["advntr", "kestrel"], "values": {**kestrel.values, **_ADVNTR_BASELINE}}
    )


def _write_capture(path: Path, scores: tuple[str, ...], policy: CallerPolicyValues | None = None) -> CallerPolicyValues:
    """Serialize one complete synthetic capture and return its baseline policy."""
    capture = _capture(_frame(scores), kestrel_config(), policy=policy)
    path.write_bytes(canonical_json_bytes(kestrel_capture_document(capture)))
    return capture.baseline_policy


def _write_manifests(
    tmp_path: Path,
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]],
    *,
    natives: dict[str, str] | None = None,
    advntr: bool = False,
    advntr_baseline: bool = False,
    advntr_calls: Mapping[str, bool | None] | None = None,
) -> tuple[Path, Path]:
    """Write the cohort TSV, the capture TSV and every capture the pair references.

    Args:
        tmp_path: Private directory the synthetic inputs are written into.
        cohort: Sample -> (genotype, Depth_Score labels, group id or None).
        natives: Optional sample -> ``"production"``, ``"negative"`` or ``"corrupt"``.
        advntr: Whether to declare an adVNTR capture column. Without ``advntr_calls`` the
            referenced files are placeholders the mocked native grid never opens.
        advntr_baseline: Whether the Kestrel captures declare an adVNTR-bearing baseline.
        advntr_calls: Sample -> the native baseline adVNTR call its one-record capture
            records, so adVNTR baseline parity is proven against real capture files.

    Returns:
        The cohort manifest path and the capture association manifest path.
    """
    tmp_path.mkdir(mode=0o700, parents=True, exist_ok=True)
    cohort_path = tmp_path / "cohort.tsv"
    captures_path = tmp_path / "captures.tsv"
    cohort_lines = ["sample_id\tbam\tassembly\tgenotype\tgroup_id"]
    header = "sample_id\tkestrel_capture\tnative_kestrel" + ("\tadvntr_capture" if advntr else "")
    capture_lines = [header]
    for name, (genotype, scores, group) in cohort.items():
        cohort_lines.append(f"{name}\t{name}.bam\tGRCh38\t{genotype}\t{group or ''}")
        capture = tmp_path / f"{name}.capture.json"
        _write_capture(capture, scores, _advntr_baseline_policy() if advntr_baseline else None)
        native = ""
        request = (natives or {}).get(name)
        if request is not None:
            native_path = tmp_path / f"{name}.native.tsv"
            if request == "negative":
                _native_negative(native_path)
            else:
                _native_production(
                    native_path, capture, changes={"Depth_Score": "0.99"} if request == "corrupt" else None
                )
            native = native_path.name
        row = f"{name}\t{capture.name}\t{native}"
        if advntr:
            locus = tmp_path / f"{name}.advntr.jsonl"
            if advntr_calls is None:
                locus.write_text("{}\n", encoding="utf-8")
            else:
                parity_capture(locus, 17 + len(capture_lines), advntr_calls.get(name, False))
            row += f"\t{locus.name}"
        capture_lines.append(row)
    cohort_path.write_text("\n".join(cohort_lines) + "\n", encoding="utf-8")
    captures_path.write_text("\n".join(capture_lines) + "\n", encoding="utf-8")
    return cohort_path, captures_path


STANDARD_COHORT: dict[str, tuple[str, tuple[str, ...], str | None]] = {
    "specimen-alpha": ("positive", ("0.014",), "family-1"),
    "specimen-alpha-repeat": ("positive", ("0.014",), "family-1"),
    "specimen-bravo": ("positive", ("0.004",), None),
    "specimen-charlie": ("negative", ("0.0004",), None),
    "specimen-delta": ("negative", ("0.001",), None),
    "specimen-echo": ("positive", (), None),
    "specimen-foxtrot": ("unknown", ("0.014",), None),
}


def _namespace(cohort_path: Path, captures_path: Path, **overrides: Any) -> argparse.Namespace:
    """Build the parsed-argument namespace the entry point reads."""
    values: dict[str, Any] = {
        "manifest": cohort_path,
        "captures": captures_path,
        "objective": "max-sensitivity-at-specificity",
        "min_sensitivity": None,
        "min_specificity": 1.0,
        "caller": "kestrel",
        "axes": None,
        "max_breakpoints": None,
        "folds": 3,
        "seed": 20260915,
        "workers": 1,
        "advntr_executable": None,
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def _run(
    tmp_path: Path,
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] | None = None,
    *,
    natives: dict[str, str] | None = None,
    **overrides: Any,
) -> tuple[bool, dict[str, Any], Path]:
    """Run the command through the real atomic adapter and read its report back."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = _write_manifests(tmp_path, cohort or STANDARD_COHORT, natives=natives)
    output = tmp_path / "derived"
    args = _namespace(cohort_path, captures_path, **overrides)
    successful = atomic_output(output, lambda staging: run_cutoff_optimization(args, staging))
    document = json.loads((output / "report.json").read_bytes())
    assert isinstance(document, dict)
    return successful, document, output


def _by_value(document: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Index the reported cutoff rows by the axis value each one tested."""
    return {repr(row["value"]): row for row in document["cutoffs"] if row["policy_id"] != "baseline"}


def test_the_reported_confusion_matrix_matches_the_hand_computed_oracle(tmp_path: Path) -> None:
    """Every tested cutoff agrees with the matrix written out in this module."""
    successful, document, _ = _run(tmp_path)

    assert successful is True
    assert document["schema_version"] == "calibration-cutoff-report-v1"
    rows = _by_value(document)
    assert sorted(rows) == sorted(ORACLE)
    for value, expected in ORACLE.items():
        counts = rows[value]["counts"]
        for field, number in expected.items():
            assert counts[field] == number, f"{value} {field}"


def test_the_declared_objective_and_its_constraints_are_reported_verbatim(tmp_path: Path) -> None:
    """The published decision states which quantity was maximized, and under what floor."""
    _, document, _ = _run(tmp_path, objective="youden-j", min_specificity=None, min_sensitivity=0.5)

    assert document["objective"] == {
        "objective": "youden-j",
        "min_sensitivity": 0.5,
        "min_specificity": None,
    }


def test_the_rescue_needs_all_three_linked_gates_and_not_the_floor_alone(tmp_path: Path) -> None:
    """A lowered floor on its own moves a row from one Negative verdict to another."""
    _, document, _ = _run(tmp_path)
    selected = document["selection"]

    assert selected["value"] == 0.004
    moved = {row["pointer"]: row for row in document["old_versus_derived"] if row["changed"]}
    assert set(moved) == {_FLOOR, _LOW, _GG}
    assert moved[_FLOOR]["baseline_value"] == 0.00469
    assert moved[_FLOOR]["derived_value"] == 0.004

    # The same rescue, attempted with the floor alone, does not reach a called endpoint.
    capture_path = tmp_path / "specimen-bravo.capture.json"
    capture = decode_kestrel_capture(load_strict_json_object(capture_path.read_bytes()))
    floor_only = _candidate(capture.baseline_policy, **{_FLOOR: 0.004})
    linked = _candidate(capture.baseline_policy, **{_FLOOR: 0.004, _LOW: 0.004, _GG: 0.004})
    kwargs = {"capture_policy_sha256": capture.provenance.capture_policy_sha256}

    assert replay_kestrel_capture(capture, floor_only, **kwargs).disposition == "no-call"
    assert replay_kestrel_capture(capture, linked, **kwargs).disposition == "called"


def test_lowering_the_cutoff_far_enough_admits_a_negative_and_the_report_shows_it(tmp_path: Path) -> None:
    """The specificity the rescue costs is a reported number, not an omission."""
    _, document, _ = _run(tmp_path)
    rows = _by_value(document)

    assert rows["0.004"]["counts"]["false_positives"] == 0
    assert rows["0.004"]["counts"]["specificity"] == 1.0
    assert rows["0.001"]["counts"]["false_positives"] == 1
    assert rows["0.001"]["counts"]["specificity"] == 0.5
    assert rows["0.0004"]["counts"]["specificity"] == 0.0
    # The one-sided upper bound admits a materially non-zero rate even at zero observed FP.
    upper = rows["0.004"]["metrics"]["fpr_one_sided_upper"]
    assert 0.0 < upper < 1.0


def test_the_baseline_arm_reproduces_the_shipped_outcome_exactly(tmp_path: Path) -> None:
    """The anchor candidate and the capture-derived baseline must agree per sample."""
    _, document, _ = _run(
        tmp_path,
        {
            "specimen-alpha": ("positive", ("0.014",), None),
            "specimen-charlie": ("negative", ("0.0004",), None),
        },
        natives={"specimen-alpha": "production", "specimen-charlie": "negative"},
    )
    parity = document["baseline_parity"]

    assert parity["proven"] is True
    assert parity["mismatches"] == []
    assert parity["native_exact_count"] == 2
    assert parity["capture_replay_authoritative_count"] == 0


def test_a_corrupted_native_baseline_aborts_the_run(tmp_path: Path) -> None:
    """A curve on evidence that cannot reproduce the shipped result is worthless."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] = {
        "specimen-alpha": ("positive", ("0.014",), None),
        "specimen-charlie": ("negative", ("0.0004",), None),
    }
    cohort_path, captures_path = _write_manifests(
        tmp_path, cohort, natives={"specimen-alpha": "corrupt", "specimen-charlie": "negative"}
    )
    args = _namespace(cohort_path, captures_path)

    with pytest.raises(ValueError, match="baseline parity"):
        atomic_output(tmp_path / "derived", lambda staging: run_cutoff_optimization(args, staging))


def test_a_baseline_arm_that_disagrees_with_its_anchor_aborts_the_run(tmp_path: Path) -> None:
    """The parity proof is a real comparison, not a comment describing one."""
    from vntyper.scripts import calibration_cutoff_optimize as module

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    args = _namespace(cohort_path, captures_path)
    real = module.replay_kestrel_grid

    def doctored(*call_args: Any, **call_kwargs: Any) -> Any:
        result = real(*call_args, **call_kwargs)
        anchor = next(name for name in result.policy_ids if name.endswith("-0003"))
        broken = dict(result.observations)
        endpoints = dict(broken[anchor])
        key = "specimen-alpha"
        endpoints[key] = replace(endpoints[key], disposition="no-call", called_positive=False)
        broken[anchor] = MappingProxyType(endpoints)
        return replace(result, observations=MappingProxyType(broken))

    with patch.object(module, "replay_kestrel_grid", doctored), pytest.raises(ValueError, match="baseline parity"):
        atomic_output(tmp_path / "derived", lambda staging: module.run_cutoff_optimization(args, staging))


def test_an_infeasible_objective_still_writes_the_report_and_names_the_constraint(tmp_path: Path) -> None:
    """A failed selection is an explicit published outcome, not a silent fallback."""
    successful, document, output = _run(tmp_path, objective="sensitivity", min_specificity=None, min_sensitivity=0.99)

    assert successful is False
    assert document["status"] == "infeasible"
    infeasible = document["selection"]["infeasible"]
    assert infeasible["unsatisfiable_constraints"] == ["min_sensitivity"]
    assert infeasible["best_achievable"]["min_sensitivity"] == pytest.approx(2 / 3)
    assert document["selection"]["policy_id"] is None
    assert (output / "report.html").is_file()
    assert document["profile"]["status"] == "unavailable"


def test_constraints_that_are_reachable_apart_but_not_together_are_named_as_such(tmp_path: Path) -> None:
    """ "Neither floor is too high, but no one cutoff meets both" is its own finding."""
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] = {
        # Both samples carry the identical Depth_Score, so no threshold separates them.
        "specimen-alpha": ("positive", ("0.004",), None),
        "specimen-charlie": ("negative", ("0.004",), None),
    }
    successful, document, _ = _run(tmp_path, cohort, objective="youden-j", min_sensitivity=0.5, min_specificity=0.5)
    infeasible = document["selection"]["infeasible"]

    assert successful is False
    assert infeasible["unsatisfiable_constraints"] == []
    assert infeasible["best_achievable"] == {"min_sensitivity": 1.0, "min_specificity": 1.0}
    assert "no single tested cutoff satisfies them all" in infeasible["note"]


def test_unknown_truth_and_no_calls_keep_fixed_denominators_at_every_cutoff(tmp_path: Path) -> None:
    """A no-call stays inside its truth class and unknown truth keeps its own count."""
    _, document, _ = _run(tmp_path)

    for row in document["cutoffs"]:
        counts = row["counts"]
        assert counts["positive_count"] == 3
        assert counts["negative_count"] == 2
        assert counts["unknown_truth_count"] == 1
        assert counts["positive_no_calls"] == 1
        assert counts["negative_no_calls"] == 0
    assert document["truth_set"]["by_genotype"] == {"negative": 2, "positive": 3, "unknown": 1}
    assert document["truth_set"]["by_assembly"] == {"GRCh38": 6}


def test_only_the_group_representative_is_scored_and_no_group_straddles_a_fold(tmp_path: Path) -> None:
    """Biological duplicates are dropped with a stated reason before anything is scored."""
    _, document, _ = _run(tmp_path)
    dropped = document["truth_set"]["dropped_duplicates"]

    assert [row["sample_id"] for row in dropped] == ["specimen-alpha-repeat"]
    assert dropped[0]["group_id"] == "declared:family-1"
    assert dropped[0]["representative"] == "specimen-alpha"
    assert dropped[0]["reason"] == "not-the-first-seen-representative-of-its-declared-group"
    assert document["truth_set"]["primary_count"] == 6
    assert document["truth_set"]["dropped_count"] == 1

    scored = {row["key"] for row in document["evaluation"]["rows"]}
    assert "specimen-alpha-repeat" not in scored
    folds: dict[str, set[int | None]] = {}
    for row in document["evaluation"]["rows"]:
        folds.setdefault(row["group_key"], set()).add(row["fold"])
    assert all(len(seen) == 1 for seen in folds.values())


def test_the_held_out_evaluation_is_reported_separately_from_the_full_data_fit(tmp_path: Path) -> None:
    """Selection happens inside training folds; the pooled held-out counts stand apart."""
    _, document, output = _run(tmp_path)
    evaluation = document["evaluation"]
    page = (output / "report.html").read_text(encoding="utf-8")

    assert page.index("Held-out performance (cross-validated)") < page.index("Every tested cutoff (descriptive")
    assert evaluation["held_out"] is not None
    assert evaluation["full_data_operating_points"].keys() >= {"baseline"}
    assert "descriptive" in evaluation["full_data_operating_points_scope"]
    assert evaluation["folds"], "grouped cross-validation produced no folds"
    for fold in evaluation["folds"]:
        assert set(fold["training_keys"]).isdisjoint(fold["held_out_keys"])


def test_the_exported_profile_round_trips_and_equals_the_selected_policy(tmp_path: Path) -> None:
    """The runtime resolver must recover exactly the policy the report published."""
    from vntyper.scripts.decision_profile import resolve_research_decision_profile

    _, document, output = _run(tmp_path)
    profile_path = output / "research-decision-profile.json"

    assert document["profile"]["round_trip_matches_selected_policy"] is True
    assert document["profile"]["sha256"] == document["profile"]["round_trip_sha256"]
    assert document["usage_hint"].startswith("vntyper pipeline --research-decision-profile")
    resolved = resolve_research_decision_profile(profile_path)
    kestrel = cast(Mapping[str, Any], resolved.components["kestrel"])
    assert kestrel["confidence_assignment"]["reporting_floor"] == 0.004
    assert kestrel["confidence_assignment"]["depth_score_thresholds"]["low"] == 0.004
    assert kestrel["alt_filtering"]["gg_depth_score_threshold"] == 0.004
    assert kestrel["confidence_assignment"]["depth_score_thresholds"]["high"] == 0.00515


def test_the_plateau_interval_is_wider_than_the_single_selected_value(tmp_path: Path) -> None:
    """Quoting one value from inside a step-function plateau is false precision."""
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] = {
        "specimen-alpha": ("positive", ("0.014",), None),
        "specimen-golf": ("positive", ("0.008",), None),
        "specimen-charlie": ("negative", ("0.0004",), None),
    }
    _, document, _ = _run(tmp_path, cohort, objective="youden-j", min_specificity=None)
    plateau = document["selection"]["plateau"]

    assert plateau["equivalent_values"] == [0.00469, 0.008]
    assert plateau["interval_low"] == 0.00469
    assert plateau["interval_high"] == 0.008
    assert plateau["single_value"] is False
    assert plateau["open_below"] == 0.0004
    # The fold holding out alpha trains on a maximum of 0.008, so its own sentinel just above
    # 0.008 is replayed too; it is the nearest value at which the outcome changes (golf lost).
    assert plateau["open_above"] == math.nextafter(0.008, math.inf)


def test_the_boundary_support_warns_when_no_truth_sample_constrains_the_cutoff(tmp_path: Path) -> None:
    """A cutoff no labelled sample sits near is an unsupported cutoff, and says so."""
    _, document, _ = _run(tmp_path)
    support = document["boundary_support"]

    assert support["positives_within_band"] >= 1
    assert support["negatives_within_band"] >= 1
    assert support["warnings"] == []


def test_the_axis_curve_and_the_joint_points_are_reported_separately(tmp_path: Path) -> None:
    """A two-axis run publishes one curve per axis and a labelled table beside them."""
    _, document, _ = _run(tmp_path, axes=[DEPTH_FLOOR_LINKED, GG_GATE_INDEPENDENT])

    axes = [curve["axis"] for curve in document["curves"]]
    assert axes == [DEPTH_FLOOR_LINKED, GG_GATE_INDEPENDENT]
    for curve in document["curves"]:
        assert curve["schema_version"] == "calibration-cutoff-curve-v1"
        assert curve["points"], "a curve with no points is not a curve"
    joint = document["joint_points"]
    assert joint is not None
    assert joint["schema_version"] == "calibration-cutoff-joint-v1"
    assert all("threshold" not in row for row in joint["points"])


def test_the_rejected_candidate_values_are_published_with_their_reasons(tmp_path: Path) -> None:
    """A breakpoint the policy decoder refuses is reported, not silently dropped."""
    _, document, _ = _run(tmp_path, axes=[DEPTH_FLOOR_LINKED])

    assert isinstance(document["rejected_candidates"], list)
    for row in document["rejected_candidates"]:
        assert set(row) == {"axis", "value", "reason"}
        assert row["reason"]


def test_the_output_directory_and_every_file_it_holds_are_private(tmp_path: Path) -> None:
    """The output holds cohort results, so it is 0700 with 0600 files throughout."""
    _, _, output = _run(tmp_path)

    assert stat.S_IMODE(output.stat().st_mode) == 0o700
    written = sorted(path.name for path in output.iterdir())
    assert written == [
        "checksums.json",
        "cutoffs.tsv",
        "folds.tsv",
        "old-versus-derived.tsv",
        "rejected-candidates.tsv",
        "report.html",
        "report.json",
        "research-decision-profile.json",
        "roc-pr-curves.tsv",
    ]
    for path in output.iterdir():
        assert stat.S_IMODE(path.stat().st_mode) == 0o600, path.name


def test_the_worker_count_reaches_the_replay(tmp_path: Path) -> None:
    """``--workers`` is plumbed through rather than accepted and discarded."""
    from vntyper.scripts import calibration_cutoff_optimize as module

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    args = _namespace(cohort_path, captures_path, workers=3)
    real = module.replay_kestrel_grid
    seen: list[int] = []

    def record(*call_args: Any, **call_kwargs: Any) -> Any:
        seen.append(call_kwargs["workers"])
        return real(*call_args, **{**call_kwargs, "workers": 1})

    with patch.object(module, "replay_kestrel_grid", record):
        atomic_output(tmp_path / "derived", lambda staging: module.run_cutoff_optimization(args, staging))

    assert seen == [3]


@pytest.mark.parametrize(
    ("damage", "message"),
    [
        ("missing", "missing, unsafe, or unreadable"),
        ("malformed", "is not strict JSON"),
        ("empty-manifest", "at least one sample"),
    ],
)
def test_broken_capture_evidence_fails_with_a_clear_message(tmp_path: Path, damage: str, message: str) -> None:
    """Missing, malformed and empty evidence are refused before anything is scored."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    if damage == "missing":
        (tmp_path / "specimen-alpha.capture.json").unlink()
    elif damage == "malformed":
        (tmp_path / "specimen-alpha.capture.json").write_bytes(b"")
    else:
        captures_path.write_text("sample_id\tkestrel_capture\tnative_kestrel\n", encoding="utf-8")
    args = _namespace(cohort_path, captures_path)

    with pytest.raises(ValueError, match=message):
        atomic_output(tmp_path / "derived", lambda staging: run_cutoff_optimization(args, staging))


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"objective": "not-an-objective"}, "unsupported cutoff search objective"),
        ({"objective": None}, "unsupported cutoff search objective"),
        ({"axes": []}, "at least one axis"),
        ({"axes": ["not-an-axis"]}, "cutoff axis name must be"),
        ({"axes": "depth_floor_linked"}, "at least one axis"),
        ({"axes": [DEPTH_FLOOR_LINKED, DEPTH_FLOOR_LINKED]}, "at most once"),
        ({"caller": "sideways"}, "caller must be kestrel, advntr, or both"),
        ({"caller": "both", "advntr_executable": None}, "requires --advntr-executable"),
        ({"caller": "advntr", "advntr_executable": None}, "requires --advntr-executable"),
        ({"axes": [ADVNTR_CUTOFF]}, r"--caller kestrel cannot search the adVNTR axis \['advntr_cutoff'\]"),
        (
            {"caller": "advntr", "axes": [DEPTH_FLOOR_LINKED, ADVNTR_CUTOFF, GG_GATE_INDEPENDENT]},
            r"--caller advntr cannot search the Kestrel axis \['depth_floor_linked', 'gg_gate_independent'\]",
        ),
        ({"manifest": "cohort.tsv"}, "must be Paths"),
        ({"advntr_executable": "advntr"}, "must be Paths"),
        ({"workers": 0}, "workers must be a positive integer"),
        ({"max_breakpoints": 2}, "max_breakpoints must be an integer of at least 3"),
        ({"folds": 1}, "folds must be an integer"),
    ],
)
def test_invalid_arguments_are_refused_before_any_evidence_is_read(
    tmp_path: Path, overrides: dict[str, Any], message: str
) -> None:
    """Argument validation happens before dependent I/O, as the cohort path does."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    args = _namespace(cohort_path, captures_path, **overrides)
    staging = tmp_path / "staging"
    staging.mkdir(mode=0o700)

    with pytest.raises(ValueError, match=message):
        run_cutoff_optimization(args, staging)


def test_a_staging_directory_that_already_holds_artifacts_is_refused(tmp_path: Path) -> None:
    """The entry point owns an empty private directory, never one with prior content."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    staging = tmp_path / "staging"
    staging.mkdir(mode=0o700)
    (staging / "stale.json").write_text("{}\n", encoding="utf-8")

    with pytest.raises(ValueError, match="empty staged directory"):
        run_cutoff_optimization(_namespace(cohort_path, captures_path), staging)


def test_an_unset_axis_option_derives_the_linked_depth_axis(tmp_path: Path) -> None:
    """``--axis`` is repeatable, so argparse leaves it None and the run owns the default."""
    _, document, _ = _run(tmp_path, axes=None)

    assert [entry["axis"] for entry in document["axes"]] == [DEPTH_FLOOR_LINKED]
    assert document["selection"]["axis"] == DEPTH_FLOOR_LINKED


def test_the_breakpoint_cap_subsamples_the_axis_and_says_so(tmp_path: Path) -> None:
    """Every capped inventory keeps its ends and its baseline, and the report records the cap.

    With a cap of 3, each inventory keeps its smallest value, the anchor 0.00469 and its
    sentinel S. The full-data inventory and folds 0 and 1 (which train on charlie's 0.0004)
    keep {0.0004, 0.00469, S}; fold 2 (charlie held out) keeps {0.001, 0.00469, S}. The
    replayed union is therefore {0.0004, 0.001, 0.00469, S}, with 0.001 derived by a fold only.
    """
    _, uncapped, _ = _run(tmp_path / "all", axes=[DEPTH_FLOOR_LINKED])
    _, capped, _ = _run(tmp_path / "few", axes=[DEPTH_FLOOR_LINKED], max_breakpoints=3)

    assert uncapped["axes"][0]["capped"] is False
    assert uncapped["axes"][0]["values"] == [0.0004, 0.001, 0.004, 0.00469, 0.014, SENTINEL]
    assert uncapped["axes"][0]["endpoint_sentinel"] == SENTINEL
    assert uncapped["axes"][0]["breakpoint_completeness"] == "complete"
    assert capped["max_breakpoints"] == 3
    axis = capped["axes"][0]
    assert axis["capped"] is True and axis["full_data_capped"] is True
    assert axis["fold_capped"] == [0, 1, 2]
    assert axis["breakpoint_completeness"] == "capped-subsample"
    assert axis["values"] == [0.0004, 0.001, 0.00469, SENTINEL]
    assert axis["fold_only_values"] == 1
    assert len(axis["values"]) - axis["fold_only_values"] <= 3  # the full-data inventory respects the cap
    # Each fold admits at most its capped inventory plus the baseline.
    assert [fold["admissible_candidates"] for fold in capped["evaluation"]["folds"]] == [4, 4, 4]


def test_a_capture_manifest_that_omits_a_primary_sample_is_refused(tmp_path: Path) -> None:
    """Every scored sample must have evidence; a partial roster is a manifest defect."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    kept = [line for line in captures_path.read_text(encoding="utf-8").splitlines() if "specimen-bravo" not in line]
    captures_path.write_text("\n".join(kept) + "\n", encoding="utf-8")
    args = _namespace(cohort_path, captures_path)

    with pytest.raises(ValueError, match="declare every primary cohort sample"):
        atomic_output(tmp_path / "derived", lambda staging: run_cutoff_optimization(args, staging))


def _v(support: int, pvalue: float) -> AdvntrVisit:
    return AdvntrVisit(support, Fraction(pvalue))


#: Scored adVNTR visits every synthetic native replay reports (support 5, so the baseline
#: support 3 admits each). The decisive legacy p-value ``q`` per sample:
#:
#: ====================  =======  ==========
#: Sample                q        truth
#: ====================  =======  ==========
#: ``specimen-alpha``    0.0004   positive
#: ``specimen-bravo``    0.004    positive
#: ``specimen-delta``    0.006    negative
#: ``specimen-charlie``  0.02     negative
#: ``specimen-echo``     --       positive (no visit: never called)
#: ``specimen-foxtrot``  0.2      unknown
#: ====================  =======  ==========
PROBE_VISITS: dict[str, tuple[AdvntrVisit, ...] | None] = {
    "specimen-alpha": (_v(5, 0.0004),),
    "specimen-bravo": (_v(5, 0.004),),
    "specimen-charlie": (_v(5, 0.02),),
    "specimen-delta": (_v(5, 0.006),),
    "specimen-echo": (),
    "specimen-foxtrot": (_v(5, 0.2),),
}


def _up(value: float) -> float:
    return math.nextafter(value, math.inf)


#: The hand-computed adVNTR cutoff oracle over the scored roster (``alpha-repeat`` is a
#: dropped duplicate): value -> (TP, FN, TN, FP). A sample is called when ``q < cutoff``.
ADVNTR_ORACLE: dict[float, tuple[int, int, int, int]] = {
    0.0004: (0, 3, 2, 0),  # the sentinel min(q): calls nobody
    _up(0.0004): (1, 2, 2, 0),  # alpha
    0.001: (1, 2, 2, 0),  # the baseline: alpha
    _up(0.004): (2, 1, 2, 0),  # alpha, bravo
    _up(0.006): (2, 1, 1, 1),  # + delta
    _up(0.02): (2, 1, 0, 2),  # + charlie
    _up(0.2): (2, 1, 0, 2),  # + foxtrot (unknown truth)
}

GridOverride = Callable[[Any, Mapping[str, CallerPolicyValues]], Any]


def _native_grid(
    seen: list[dict[str, Any]],
    *,
    visits: Mapping[str, tuple[AdvntrVisit, ...] | None] = PROBE_VISITS,
    override: GridOverride | None = None,
) -> Callable[..., Any]:
    """A native adVNTR grid fake: every policy's calls follow the legacy rule over ``visits``."""

    def grid(capture_paths: Mapping[str, Path], policies: Mapping[str, CallerPolicyValues], **kwargs: Any) -> Any:
        seen.append({"captures": dict(capture_paths), "policies": dict(policies), **kwargs})
        thresholds = {pid: (float(p.values[_ADV_CUT]), int(p.values[_ADV_SUP])) for pid, p in policies.items()}  # type: ignore[arg-type]
        result = advntr_grid_result(
            dict.fromkeys(policies),
            visits=dict.fromkeys(policies, visits),
            thresholds=thresholds,
            baseline_policy_id=kwargs["baseline_policy_id"],
            capture_paths=capture_paths,
        )
        return result if override is None else override(result, policies)

    return grid


def _run_advntr(
    tmp_path: Path,
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] | None = None,
    *,
    seen: list[dict[str, Any]] | None = None,
    visits: Mapping[str, tuple[AdvntrVisit, ...] | None] = PROBE_VISITS,
    override: GridOverride | None = None,
    **overrides: Any,
) -> tuple[bool, dict[str, Any], Path]:
    """Run an adVNTR-bearing optimize with both native grid seams replaced by one fake.

    The captures carry an adVNTR-bearing baseline and one-record adVNTR captures whose
    recorded baseline call is the legacy rule's own, so adVNTR baseline parity runs for real.
    """
    from vntyper.scripts import calibration_cutoff_advntr_axes as axes_module
    from vntyper.scripts import calibration_cutoff_optimize as module

    calls = {key: predicted_call(items, 0.001, 3) for key, items in visits.items()}
    cohort_path, captures_path = _write_manifests(
        tmp_path, cohort or STANDARD_COHORT, advntr=True, advntr_baseline=True, advntr_calls=calls
    )
    args = _namespace(cohort_path, captures_path, advntr_executable=tmp_path / "advntr", **overrides)
    grid = _native_grid([] if seen is None else seen, visits=visits, override=override)
    output = tmp_path / "derived"
    with (
        patch.object(module, "evaluate_advntr_cutoff_grid", grid),
        patch.object(axes_module, "evaluate_advntr_cutoff_grid", grid),
    ):
        successful = atomic_output(output, lambda staging: module.run_cutoff_optimization(args, staging))
    document = json.loads((output / "report.json").read_bytes())
    assert isinstance(document, dict)
    return successful, document, output


def _signature(policy: CallerPolicyValues) -> tuple[object, ...]:
    return tuple(policy.values[pointer] for pointer in ADVNTR_CALLER_POLICY_POINTERS)


def test_the_advntr_arm_is_replayed_natively_and_never_approximated(tmp_path: Path) -> None:
    """adVNTR statistics come from the installed evaluator; this only mocks the seam."""
    seen: list[dict[str, Any]] = []
    rescued = {**PROBE_VISITS, "specimen-echo": (_v(5, 0.0001),)}  # adVNTR calls echo at baseline
    _, document, _ = _run_advntr(
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

    _, document, _ = _run_advntr(tmp_path, caller="both", axes=[DEPTH_FLOOR_LINKED], min_specificity=1.0)
    scope = document["search_scope"]

    assert scope["searched_caller"] == "kestrel"
    assert scope["searched_axes"] == [DEPTH_FLOOR_LINKED]
    assert scope["advntr_policy"] == "held-at-baseline"
    assert scope["advntr_distinct_executions"] == 1
    assert "#269" in scope["note"]
    page = render_cutoff_report_html(document)
    assert "adVNTR arm was replayed at its baseline policy" in page


def test_a_kestrel_run_records_that_advntr_was_not_evaluated(tmp_path: Path) -> None:
    """The scope section is present on every run, not only on the adVNTR ones."""
    _, document, _ = _run(tmp_path)

    assert document["search_scope"]["advntr_policy"] == "not-evaluated"
    assert document["search_scope"]["advntr_distinct_executions"] is None
    assert document["provenance"]["advntr"] is None


def test_an_unexpected_extra_advntr_execution_aborts_the_run(tmp_path: Path) -> None:
    """The execution count is checked against the distinct adVNTR policies, never assumed."""

    def extra(result: Any, policies: Mapping[str, CallerPolicyValues]) -> Any:
        rows = list(result.policies)
        kestrel_only = next(index for index, row in enumerate(rows) if row.policy_id != result.baseline_policy_id)
        rows[kestrel_only] = replace(rows[kestrel_only], execution_id="exec-unexpected")
        return replace(result, policies=tuple(rows))

    with pytest.raises(ValueError, match=r"adVNTR executions \(2\) differ from the distinct adVNTR policies"):
        _run_advntr(tmp_path, override=extra, caller="both", axes=[DEPTH_FLOOR_LINKED])


def test_caller_advntr_searches_the_advntr_cutoff_axis_by_native_replay(tmp_path: Path) -> None:
    """The adVNTR cutoff axis is derived from probe replays and scored against the hand oracle.

    Outer folds (3 folds, seed 20260915): fold 0 holds out bravo and foxtrot, fold 1 alpha and
    delta, fold 2 charlie and echo. No fold derives a value the full-data axis lacks: fold 0
    trains on q {0.0004, 0.006, 0.02}, fold 1 on {0.004, 0.02, 0.2} (the baseline 0.001 already
    rejects all, so no sentinel), fold 2 on {0.0004, 0.004, 0.006, 0.2}.
    """
    seen: list[dict[str, Any]] = []
    successful, document, _ = _run_advntr(tmp_path, seen=seen, caller="advntr", min_specificity=1.0)

    assert successful is True
    rows = {row["value"]: row for row in document["cutoffs"] if row["policy_id"] != "baseline"}
    assert sorted(rows) == sorted(ADVNTR_ORACLE)
    for value, (tp, fn, tn, fp) in ADVNTR_ORACLE.items():
        counts = rows[value]["counts"]
        assert (counts["true_positives"], counts["false_negatives"]) == (tp, fn), value
        assert (counts["true_negatives"], counts["false_positives"]) == (tn, fp), value
    assert document["selection"]["axis"] == ADVNTR_CUTOFF
    assert document["selection"]["value"] == _up(0.004)
    (axis,) = document["axes"]
    assert axis["axis"] == ADVNTR_CUTOFF
    assert axis["caller"] == "advntr"
    assert axis["comparator"] == "<"
    assert axis["breakpoint_completeness"] == "complete"
    assert axis["fold_only_values"] == 0
    assert document["evaluation"]["fold_admissibility"] == "training-derived-inventories"
    changed = {row["pointer"]: row for row in document["old_versus_derived"] if row["changed"]}
    assert set(changed) == {_ADV_CUT}
    assert changed[_ADV_CUT]["baseline_value"] == 0.001
    assert changed[_ADV_CUT]["derived_value"] == _up(0.004)
    assert document["profile"]["round_trip_matches_selected_policy"] is True

    probe, main = seen
    assert sorted(probe["policies"]) == ["baseline", f"probe-{ADVNTR_CUTOFF}"]
    assert probe["output"].name == "advntr-probe" and main["output"].name == "advntr"
    # X4: one execution per distinct adVNTR signature passed to the main grid, the anchor once.
    assert len({_signature(policy) for policy in main["policies"].values()}) == 7
    assert document["search_scope"]["advntr_distinct_executions"] == 7
    advntr = document["provenance"]["advntr"]
    assert set(advntr) == {"sha256", "probe_sha256", "tool_identity", "probe_seconds", "main_seconds"}
    assert advntr["tool_identity"] == advntr_capabilities_document(synthetic_capabilities())
    assert advntr["probe_sha256"] == "d" * 64
    assert isinstance(advntr["probe_seconds"], float) and isinstance(advntr["main_seconds"], float)


def test_caller_both_searches_kestrel_and_advntr_axes_on_the_union(tmp_path: Path) -> None:
    """``--caller both`` defaults to one Kestrel and one adVNTR axis, each holding the other caller at baseline."""
    seen: list[dict[str, Any]] = []
    _, document, _ = _run_advntr(tmp_path, seen=seen, caller="both", min_specificity=1.0)

    assert [axis["axis"] for axis in document["axes"]] == [DEPTH_FLOOR_LINKED, ADVNTR_CUTOFF]
    assert [axis["caller"] for axis in document["axes"]] == ["kestrel", "advntr"]
    _, main = seen
    signatures = {_signature(policy) for policy in main["policies"].values()}
    # Every Kestrel candidate shares the baseline adVNTR signature; the adVNTR axis adds 6 more.
    assert len(signatures) == 7
    assert document["search_scope"]["advntr_distinct_executions"] == 1 + 6
    advntr_values = {row["value"] for row in document["cutoffs"] if row["policy_id"].startswith(ADVNTR_CUTOFF)}
    assert advntr_values == set(ADVNTR_ORACLE)


def test_advntr_replay_inconsistency_aborts_the_run(tmp_path: Path) -> None:
    """A native call the probe statistics do not predict aborts the run, naming the candidate."""

    def contradict(result: Any, policies: Mapping[str, CallerPolicyValues]) -> Any:
        if "baseline" in policies:
            return result  # the probe grid is left intact
        target = next(pid for pid, policy in policies.items() if policy.values[_ADV_CUT] == _up(0.004))
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
        _run_advntr(tmp_path, override=contradict, caller="advntr")


def test_advntr_probe_and_candidate_grids_must_share_their_evidence(tmp_path: Path) -> None:
    """A candidate grid replayed by a different adVNTR build is not bound to the probe's statistics."""

    def rebuilt(result: Any, policies: Mapping[str, CallerPolicyValues]) -> Any:
        return result if "baseline" in policies else replace(result, capabilities=synthetic_capabilities("2.4.1"))

    with pytest.raises(ValueError, match="not bound to the same evidence"):
        _run_advntr(tmp_path, override=rebuilt, caller="advntr")


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
    _, document, _ = _run_advntr(tmp_path, caller="advntr", min_specificity=1.0)
    evaluation = document["evaluation"]
    ids = {row["value"]: row["policy_id"] for row in document["cutoffs"] if row["policy_id"] != "baseline"}
    baseline = _advntr_baseline_policy()

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


def test_the_provenance_digests_cover_the_manifests_and_every_capture(tmp_path: Path) -> None:
    """The report commits to the exact bytes the decision was derived from."""
    _, document, _ = _run(tmp_path)
    provenance = document["provenance"]

    assert len(provenance["cohort_manifest_sha256"]) == 64
    assert len(provenance["capture_manifest_sha256"]) == 64
    assert len(provenance["baseline_policy_sha256"]) == 64
    assert len(provenance["grid_replay_sha256"]) == 64
    assert set(provenance["capture_file_sha256"]) == {
        "specimen-alpha",
        "specimen-bravo",
        "specimen-charlie",
        "specimen-delta",
        "specimen-echo",
        "specimen-foxtrot",
    }
    assert provenance["generator_version"]


def test_a_specificity_floor_reachable_only_by_rejecting_everything_is_feasible(tmp_path: Path) -> None:
    """A negative scoring above every positive can only be excluded by the endpoint sentinel."""
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] = {
        "specimen-bravo": ("positive", ("0.006",), None),
        "specimen-charlie": ("negative", ("0.014",), None),
    }
    successful, document, _ = _run(tmp_path, cohort, objective="max-sensitivity-at-specificity", min_specificity=1.0)

    assert successful is True
    assert document["selection"]["value"] == SENTINEL
    assert _by_value(document)[repr(SENTINEL)]["counts"]["specificity"] == 1.0
    others = [row["counts"]["specificity"] for row in document["cutoffs"] if row["value"] != SENTINEL]
    assert max(others) == 0.0


def test_every_fold_selects_only_breakpoints_its_training_samples_observed(tmp_path: Path) -> None:
    """A held-out sample's Depth_Score must never become the cutoff its own fold uses.

    Hand derivation (3 folds, seed 20260915: fold 0 holds out bravo and foxtrot, fold 1 alpha
    and delta, fold 2 charlie and echo). Each fold's inventory is its training scores, the
    anchor 0.00469 and the sentinel just above the training maximum (0.014 in every fold):

    * fold 0 trains on alpha 0.014, charlie 0.0004, delta 0.001: {0.0004, 0.001, 0.00469, 0.014, S}
    * fold 1 trains on bravo 0.004, charlie 0.0004, foxtrot 0.014: {0.0004, 0.004, 0.00469, 0.014, S}
    * fold 2 trains on alpha 0.014, bravo 0.004, delta 0.001, foxtrot 0.014: {0.001, 0.004, 0.00469, 0.014, S}

    so each fold admits 5 candidates plus the baseline.
    """
    _, document, _ = _run(tmp_path)
    evaluation = document["evaluation"]
    scores = {name: {float(score) for score in values} for name, (_, values, _) in STANDARD_COHORT.items()}
    ids = {row["value"]: row["policy_id"] for row in document["cutoffs"] if row["policy_id"] != "baseline"}

    assert evaluation["fold_admissibility"] == "training-derived-inventories"
    assert [fold["admissible_candidates"] for fold in evaluation["folds"]] == [6, 6, 6]
    for fold in evaluation["folds"]:
        training = {value for key in fold["training_keys"] for value in scores[key]}
        inventory = {*training, 0.00469, SENTINEL}
        admissible = {"baseline", *(ids[value] for value in inventory)}
        assert fold["admissible_candidates"] == len(admissible)
        assert fold["used_policy"] in admissible


def test_fold_inventories_derived_on_another_allocation_abort_the_run(tmp_path: Path) -> None:
    """Same fold numbers are not enough: inventories must come from the evaluation's own folds."""
    from vntyper.scripts import calibration_cutoff_optimize as module

    real = module.outer_fold_assignments

    def rotated(*call_args: Any, **call_kwargs: Any) -> dict[str, int]:
        assignments = real(*call_args, **call_kwargs)
        return {key: (fold + 1) % 3 for key, fold in assignments.items()}

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    args = _namespace(cohort_path, captures_path)
    with patch.object(module, "outer_fold_assignments", rotated), pytest.raises(ValueError, match="different outer"):
        atomic_output(tmp_path / "derived", lambda staging: module.run_cutoff_optimization(args, staging))


def test_a_singleton_cohort_degrades_to_no_fold_inventories_and_unavailable_cross_validation(
    tmp_path: Path,
) -> None:
    """One independent group cannot be cross-validated: no fold inventories, no folds.

    A single group carries a single truth label, so the ROC/PR curve (which needs both
    classes) cannot be published and the run stops there, as it did before fold-local
    inventories existed. What this pins is that the fold machinery itself degrades cleanly:
    the evaluation receives ``fold_inventories == {}`` and reports unavailable cross-validation.
    """
    from vntyper.scripts import calibration_cutoff_optimize as module

    cohort: dict[str, tuple[str, tuple[str, ...], str | None]] = {
        "specimen-alpha": ("positive", ("0.014",), "family-1"),
        "specimen-alpha-repeat": ("positive", ("0.004",), "family-1"),
    }
    cohort_path, captures_path = _write_manifests(tmp_path, cohort)
    args = _namespace(cohort_path, captures_path, objective="sensitivity", min_specificity=None)
    real = module.evaluate_cutoff_arms
    seen: list[tuple[dict[str, Any], dict[str, Any]]] = []

    def record(*call_args: Any, **call_kwargs: Any) -> Any:
        result = real(*call_args, **call_kwargs)
        seen.append((call_kwargs, result))
        return result

    with patch.object(module, "evaluate_cutoff_arms", record), pytest.raises(ValueError, match="both positive and"):
        atomic_output(tmp_path / "derived", lambda staging: module.run_cutoff_optimization(args, staging))

    ((kwargs, evaluation),) = seen
    assert kwargs["fold_inventories"] == {}
    assert evaluation["cross_validation_available"] is False
    assert evaluation["status_reason"] == "insufficient-cross-validation-groups"
    assert evaluation["folds"] == []
    assert evaluation["fold_admissibility"] == "training-derived-inventories"


def test_the_report_comparators_match_the_axis_definitions() -> None:
    """The optimize probe table and the axis module must agree on every comparator."""
    from vntyper.scripts.calibration_cutoff_axes import axis_comparison
    from vntyper.scripts.calibration_cutoff_optimize import AXIS_COMPARISON

    assert {name: axis_comparison(name) for name in AXIS_COMPARISON} == dict(AXIS_COMPARISON)
    assert len(AXIS_COMPARISON) == 7
    assert AXIS_COMPARISON[ADVNTR_CUTOFF] == "<" and AXIS_COMPARISON[ADVNTR_MIN_SUPPORT] == ">="
