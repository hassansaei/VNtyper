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
import stat
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from types import MappingProxyType
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest

from tests.builders import kestrel_config
from tests.unit.test_calibration_cutoff_kestrel import _native_negative, _native_production
from tests.unit.test_calibration_kestrel_replay import _candidate, _capture, _raw
from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_cutoff_axes import DEPTH_FLOOR_LINKED, GG_GATE_INDEPENDENT
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


def _write_capture(path: Path, scores: tuple[str, ...]) -> CallerPolicyValues:
    """Serialize one complete synthetic capture and return its baseline policy."""
    capture = _capture(_frame(scores), kestrel_config())
    path.write_bytes(canonical_json_bytes(kestrel_capture_document(capture)))
    return capture.baseline_policy


def _write_manifests(
    tmp_path: Path,
    cohort: dict[str, tuple[str, tuple[str, ...], str | None]],
    *,
    natives: dict[str, str] | None = None,
    advntr: bool = False,
) -> tuple[Path, Path]:
    """Write the cohort TSV, the capture TSV and every capture the pair references.

    Args:
        tmp_path: Private directory the synthetic inputs are written into.
        cohort: Sample -> (genotype, Depth_Score labels, group id or None).
        natives: Optional sample -> ``"production"``, ``"negative"`` or ``"corrupt"``.
        advntr: Whether to declare an adVNTR capture column; the referenced files exist
            but are never opened, because the native adVNTR grid is mocked out.

    Returns:
        The cohort manifest path and the capture association manifest path.
    """
    cohort_path = tmp_path / "cohort.tsv"
    captures_path = tmp_path / "captures.tsv"
    cohort_lines = ["sample_id\tbam\tassembly\tgenotype\tgroup_id"]
    header = "sample_id\tkestrel_capture\tnative_kestrel" + ("\tadvntr_capture" if advntr else "")
    capture_lines = [header]
    for name, (genotype, scores, group) in cohort.items():
        cohort_lines.append(f"{name}\t{name}.bam\tGRCh38\t{genotype}\t{group or ''}")
        capture = tmp_path / f"{name}.capture.json"
        _write_capture(capture, scores)
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
            locus.write_text("{}\n", encoding="utf-8")
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
        "axes": [DEPTH_FLOOR_LINKED],
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

    cohort = {
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
    _, document, _ = _run(tmp_path)
    evaluation = document["evaluation"]

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
    kestrel = resolved.components["kestrel"]
    assert kestrel["confidence_assignment"]["reporting_floor"] == 0.004
    assert kestrel["confidence_assignment"]["depth_score_thresholds"]["low"] == 0.004
    assert kestrel["alt_filtering"]["gg_depth_score_threshold"] == 0.004
    assert kestrel["confidence_assignment"]["depth_score_thresholds"]["high"] == 0.00515


def test_the_plateau_interval_is_wider_than_the_single_selected_value(tmp_path: Path) -> None:
    """Quoting one value from inside a step-function plateau is false precision."""
    cohort = {
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
    assert plateau["open_above"] == 0.014


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
        ({"axes": []}, "at least one axis"),
        ({"axes": ["not-an-axis"]}, "cutoff axis name must be"),
        ({"caller": "sideways"}, "caller must be kestrel, advntr, or both"),
        ({"caller": "both", "advntr_executable": None}, "requires --advntr-executable"),
        ({"manifest": "cohort.tsv"}, "must be Paths"),
        ({"workers": 0}, "workers must be a positive integer"),
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


def test_a_capture_manifest_that_omits_a_primary_sample_is_refused(tmp_path: Path) -> None:
    """Every scored sample must have evidence; a partial roster is a manifest defect."""
    from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT)
    kept = [line for line in captures_path.read_text(encoding="utf-8").splitlines() if "specimen-bravo" not in line]
    captures_path.write_text("\n".join(kept) + "\n", encoding="utf-8")
    args = _namespace(cohort_path, captures_path)

    with pytest.raises(ValueError, match="declare every primary cohort sample"):
        atomic_output(tmp_path / "derived", lambda staging: run_cutoff_optimization(args, staging))


def _advntr_result(policy_ids: Sequence[str], keys: Sequence[str], positive: str) -> Any:
    """A synthetic native adVNTR grid result; the executable itself is never run."""
    from vntyper.modules.advntr.advntr_calibration_policy import AdvntrCapabilities
    from vntyper.scripts.calibration_cutoff_advntr import (
        AdvntrCutoffGridResult,
        AdvntrCutoffPolicyResult,
        AdvntrCutoffSample,
    )

    samples = tuple(AdvntrCutoffSample(key, True, key == positive, ()) for key in keys)
    capabilities = AdvntrCapabilities("2.3.0", "synthetic", None, (), (), (), (), "a" * 64)
    return AdvntrCutoffGridResult(
        Path("/nonexistent"),
        policy_ids[0],
        "b" * 64,
        capabilities,
        tuple(AdvntrCutoffPolicyResult(name, "c" * 64, "exec", samples) for name in policy_ids),
        "d" * 64,
    )


@pytest.mark.parametrize(("caller", "expected_true_positives"), [("both", 3), ("advntr", 1)])
def test_the_advntr_arm_is_replayed_natively_and_never_approximated(
    tmp_path: Path, caller: str, expected_true_positives: int
) -> None:
    """adVNTR statistics come from the installed evaluator; this only mocks the seam."""
    from vntyper.scripts import calibration_cutoff_optimize as module

    cohort_path, captures_path = _write_manifests(tmp_path, STANDARD_COHORT, advntr=True)
    args = _namespace(
        cohort_path, captures_path, caller=caller, advntr_executable=tmp_path / "advntr", min_specificity=1.0
    )
    seen: dict[str, Any] = {}

    def grid(capture_paths: Any, policies: Any, **kwargs: Any) -> Any:
        seen.update(kwargs)
        seen["captures"] = dict(capture_paths)
        return _advntr_result(sorted(policies), sorted(capture_paths), "specimen-echo")

    output = tmp_path / "derived"
    with patch.object(module, "evaluate_advntr_cutoff_grid", grid):
        atomic_output(output, lambda staging: module.run_cutoff_optimization(args, staging))
    document = json.loads((output / "report.json").read_bytes())

    assert document["caller"] == caller
    assert seen["executable_path"] == tmp_path / "advntr"
    assert set(seen["captures"]) == {
        "specimen-alpha",
        "specimen-bravo",
        "specimen-charlie",
        "specimen-delta",
        "specimen-echo",
        "specimen-foxtrot",
    }
    assert document["provenance"]["advntr"] == {"sha256": "d" * 64}
    best = max(row["counts"]["true_positives"] for row in document["cutoffs"])
    assert best == expected_true_positives


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
