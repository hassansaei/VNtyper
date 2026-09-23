"""The cutoff report writers: canonical JSON, joinable TSVs and an offline page.

The page has to survive being opened from a private directory on a machine with no
network, so every one of these assertions about ``report.html`` is load-bearing: a single
``<script src=...>`` or web font turns a published decision into a blank page.
"""

from __future__ import annotations

import json
import stat
from pathlib import Path
from typing import Any

import pytest

pytestmark = pytest.mark.unit


def _curve(axis: str) -> dict[str, Any]:
    """One structurally complete axis curve document with two ordered points."""
    return {
        "schema_version": "calibration-cutoff-curve-v1",
        "axis": axis,
        "statistic": "Depth_Score",
        "pointers": ["/components/kestrel/confidence_assignment/reporting_floor"],
        "comparison": ">=",
        "fixed_policy_sha256": "a" * 64,
        "baseline_threshold": 0.00469,
        "thresholds": [0.00469, 0.004],
        "eligible_count": 4,
        "unknown_truth_count": 1,
        "no_call_count": 1,
        "positive_no_calls": 1,
        "negative_no_calls": 0,
        "boundary_support": {
            "positives_within_band": 1,
            "negatives_within_band": 0,
            "band_low": 0.004,
            "band_high": 0.00469,
        },
        "rejected": [{"value": 2.0, "reason": "reporting_floor must lie in [0,1]"}],
        "points": [
            {
                "candidate_id": "b" * 64,
                "threshold": 0.00469,
                "threshold_key": "0.00469",
                "true_positives": 1,
                "false_positives": 0,
                "true_negatives": 2,
                "false_negatives": 1,
                "positive_no_calls": 1,
                "negative_no_calls": 0,
                "unknown_truth_count": 1,
                "false_positive_rate": 0.0,
                "sensitivity": 0.3333333333333333,
                "precision": 1.0,
                "intervals": {},
            },
            {
                "candidate_id": "c" * 64,
                "threshold": 0.004,
                "threshold_key": "0.004",
                "true_positives": 2,
                "false_positives": 0,
                "true_negatives": 2,
                "false_negatives": 0,
                "positive_no_calls": 1,
                "negative_no_calls": 0,
                "unknown_truth_count": 1,
                "false_positive_rate": 0.0,
                "sensitivity": 0.6666666666666666,
                "precision": 1.0,
                "intervals": {},
            },
        ],
    }


def _document(**overrides: Any) -> dict[str, Any]:
    """A complete report document; the writers must not need anything beyond this."""
    document: dict[str, Any] = {
        "schema_version": "calibration-cutoff-report-v1",
        "status": "selected",
        "successful": True,
        "usage_hint": "vntyper pipeline --research-decision-profile <output>/research-decision-profile.json",
        "caller": "kestrel",
        "search_scope": {
            "searched_caller": "kestrel",
            "searched_axes": ["depth_floor_linked"],
            "advntr_policy": "not-evaluated",
            "advntr_distinct_executions": None,
            "note": "Only Kestrel axes were searched; adVNTR was not evaluated.",
        },
        "objective": {"objective": "youden-j", "min_sensitivity": None, "min_specificity": 1.0},
        "folds_requested": 3,
        "seed": 20260915,
        "workers": 1,
        "max_breakpoints": None,
        "truth_set": {
            "sample_count": 5,
            "primary_count": 4,
            "dropped_count": 1,
            "dropped_duplicates": [
                {
                    "sample_id": "specimen-alpha-repeat",
                    "group_id": "declared:family-1",
                    "representative": "specimen-alpha",
                    "reason": "not-the-first-seen-representative-of-its-declared-group",
                }
            ],
            "by_genotype": {"negative": 1, "positive": 2, "unknown": 1},
            "by_assembly": {"GRCh38": 4},
            "unknown_truth_count": 1,
            "truth_variant_count": 0,
        },
        "provenance": {
            "cohort_manifest_sha256": "d" * 64,
            "capture_manifest_sha256": "e" * 64,
            "baseline_policy_sha256": "f" * 64,
            "grid_replay_sha256": "0" * 64,
            "grid_input_sha256": "1" * 64,
            "capture_file_sha256": {"specimen-alpha": "2" * 64},
            "native_file_sha256": {"specimen-alpha": None},
            "generator_version": "2.0.35",
            "advntr": None,
        },
        "baseline_parity": {
            "proven": True,
            "anchor_candidate_ids": {"depth_floor_linked": "depth_floor_linked-0001"},
            "native_exact_count": 0,
            "capture_replay_authoritative_count": 1,
            "mismatches": [],
        },
        "axes": [],
        "cutoffs": [
            {
                "policy_id": "depth_floor_linked-0001",
                "axis": "depth_floor_linked",
                "value": 0.004,
                "policy_sha256": "c" * 64,
                "parameters": {"/components/kestrel/confidence_assignment/reporting_floor": 0.004},
                "counts": {
                    "eligible_count": 4,
                    "positive_count": 2,
                    "negative_count": 1,
                    "unknown_truth_count": 1,
                    "true_positives": 2,
                    "false_negatives": 0,
                    "true_negatives": 1,
                    "false_positives": 0,
                    "positive_no_calls": 0,
                    "negative_no_calls": 0,
                    "no_calls": 0,
                    "sensitivity": 1.0,
                    "specificity": 1.0,
                    "false_positive_rate": 0.0,
                    "precision": 1.0,
                    "f1": 1.0,
                    "balanced_accuracy": 1.0,
                    "youden_j": 1.0,
                    "no_call_rate": 0.0,
                },
                "metrics": {"fpr_one_sided_upper": 0.95},
            }
        ],
        "curves": [_curve("depth_floor_linked")],
        "joint_points": None,
        "rejected_candidates": [
            {"axis": "depth_floor_linked", "value": 2.0, "reason": "reporting_floor must lie in [0,1]"}
        ],
        "evaluation": {
            "status": "available",
            "folds": [
                {
                    "fold": 0,
                    "training_keys": ["specimen-alpha"],
                    "held_out_keys": ["specimen-charlie"],
                    "used_policy": "depth_floor_linked-0001",
                    "fallback_reason": None,
                    "selection": {"policy_id": "depth_floor_linked-0001", "reason": "selected"},
                }
            ],
            "held_out": {"counts": {"true_positives": 2, "false_positives": 0}},
            "full_data_operating_points": {},
            "full_data_operating_points_scope": "descriptive searched-cohort points",
            "rows": [],
        },
        "selection": {
            "policy_id": "depth_floor_linked-0001",
            "reason": "selected",
            "axis": "depth_floor_linked",
            "value": 0.004,
            "policy_sha256": "c" * 64,
            "plateau": {
                "axis": "depth_floor_linked",
                "comparison": ">=",
                "selected_value": 0.004,
                "equivalent_values": [0.004, 0.00469],
                "interval_low": 0.004,
                "interval_high": 0.00469,
                "open_below": 0.001,
                "open_above": 0.014,
                "width": 0.00069,
                "single_value": False,
                "note": "a threshold sweep is a step function",
            },
            "infeasible": None,
        },
        "old_versus_derived": [
            {
                "pointer": "/components/kestrel/confidence_assignment/reporting_floor",
                "baseline_value": 0.00469,
                "derived_value": 0.004,
                "changed": True,
            },
            {
                "pointer": "/components/kestrel/confidence_assignment/depth_score_thresholds/high",
                "baseline_value": 0.00515,
                "derived_value": 0.00515,
                "changed": False,
            },
        ],
        "boundary_support": {
            "positives_within_band": 1,
            "negatives_within_band": 0,
            "band_low": 0.004,
            "band_high": 0.00469,
            "warnings": ["no negative-truth sample lies inside the tested band of axis depth_floor_linked"],
        },
        "profile": {
            "status": "available",
            "path": "research-decision-profile.json",
            "sha256": "3" * 64,
            "profile_id": "vntyper-caller-generated-0123456789abcdef",
            "round_trip_matches_selected_policy": True,
            "round_trip_sha256": "3" * 64,
        },
        "limitations": "Research use only; derived cutoffs carry no deployment approval.",
    }
    document.update(overrides)
    return document


def _write(tmp_path: Path, document: dict[str, Any]) -> Path:
    from vntyper.scripts.calibration_cutoff_report import write_cutoff_reports

    output = tmp_path / "derived"
    output.mkdir(mode=0o700, parents=True)
    write_cutoff_reports(output, document)
    return output


def test_the_report_json_is_canonical_and_round_trips(tmp_path: Path) -> None:
    """The published JSON is the canonical encoding of exactly what it was given."""
    from vntyper.scripts.canonical_json import canonical_json_bytes

    document = _document()
    output = _write(tmp_path, document)
    raw = (output / "report.json").read_bytes()

    assert raw == canonical_json_bytes(document)
    assert json.loads(raw) == document


def test_every_written_file_is_private_and_checksummed(tmp_path: Path) -> None:
    """The directory holds cohort results, so nothing in it is group or world readable."""
    output = _write(tmp_path, _document())

    names = sorted(path.name for path in output.iterdir())
    assert names == [
        "checksums.json",
        "cutoffs.tsv",
        "folds.tsv",
        "old-versus-derived.tsv",
        "rejected-candidates.tsv",
        "report.html",
        "report.json",
        "roc-pr-curves.tsv",
    ]
    for path in output.iterdir():
        assert stat.S_IMODE(path.stat().st_mode) == 0o600, path.name
    checksums = json.loads((output / "checksums.json").read_bytes())
    assert checksums["schema_version"] == "calibration-checksums-v1"
    assert set(checksums["files"]) == set(names) - {"checksums.json"}


def test_the_cutoff_tsv_carries_every_count_rate_and_bound(tmp_path: Path) -> None:
    """The TSV is the joinable form of the same rows, not a shortened summary."""
    output = _write(tmp_path, _document())
    lines = (output / "cutoffs.tsv").read_text(encoding="utf-8").splitlines()

    header = lines[0].split("\t")
    assert header[:4] == ["policy_id", "axis", "value", "policy_sha256"]
    for column in ("true_positives", "false_positives", "true_negatives", "false_negatives"):
        assert column in header
    for column in ("positive_no_calls", "negative_no_calls", "unknown_truth_count"):
        assert column in header
    for column in ("sensitivity", "specificity", "precision", "f1", "youden_j"):
        assert column in header
    assert "fpr_one_sided_upper" in header
    assert len(lines) == 2
    row = dict(zip(header, lines[1].split("\t"), strict=True))
    assert row["policy_id"] == "depth_floor_linked-0001"
    assert row["true_positives"] == "2"


def test_the_curve_tsv_keeps_every_axis_apart(tmp_path: Path) -> None:
    """One file, but an explicit axis column: rows from two axes are never one curve."""
    document = _document(curves=[_curve("depth_floor_linked"), _curve("gg_gate_independent")])
    output = _write(tmp_path, document)
    lines = (output / "roc-pr-curves.tsv").read_text(encoding="utf-8").splitlines()

    header = lines[0].split("\t")
    assert header[0] == "axis"
    assert {line.split("\t")[0] for line in lines[1:]} == {"depth_floor_linked", "gg_gate_independent"}
    assert len(lines) == 5


def test_the_joint_points_are_written_only_when_a_run_produced_them(tmp_path: Path) -> None:
    """A labelled multi-axis table is a separate file, and never a curve row."""
    plain = _write(tmp_path / "plain", _document())
    assert not (plain / "joint-points.tsv").exists()

    joint = {
        "schema_version": "calibration-cutoff-joint-v1",
        "point_count": 1,
        "points": [
            {
                "label": "depth_floor_linked=0.004+gg_gate_independent=0.002",
                "values": {"/components/kestrel/confidence_assignment/reporting_floor": 0.004},
                "metrics": {"true_positives": 2, "false_positives": 0},
            }
        ],
    }
    output = _write(tmp_path / "joint", _document(joint_points=joint))
    lines = (output / "joint-points.tsv").read_text(encoding="utf-8").splitlines()

    assert lines[0].split("\t") == ["label", "values", "true_positives", "false_positives"]
    assert lines[1].startswith("depth_floor_linked=0.004+gg_gate_independent=0.002\t")


def test_the_old_versus_derived_table_lists_every_pointer(tmp_path: Path) -> None:
    """Unchanged pointers stay in the table; an omission would read as a change."""
    output = _write(tmp_path, _document())
    lines = (output / "old-versus-derived.tsv").read_text(encoding="utf-8").splitlines()

    assert lines[0].split("\t") == ["pointer", "baseline_value", "derived_value", "changed"]
    assert len(lines) == 3
    assert lines[1].endswith("\ttrue")
    assert lines[2].endswith("\tfalse")


def test_the_html_is_fully_offline(tmp_path: Path) -> None:
    """No network reference of any kind, including a stylesheet or a web font."""
    output = _write(tmp_path, _document())
    html = (output / "report.html").read_text(encoding="utf-8")

    assert "http://" not in html
    assert "https://" not in html
    assert "<script" not in html
    assert "<link " not in html
    assert "@import" not in html
    assert "<style>" in html
    assert "<svg" in html


def test_the_html_states_the_objective_the_selection_and_the_plateau(tmp_path: Path) -> None:
    """A reader must not have to open the JSON to learn what was decided and why."""
    output = _write(tmp_path, _document())
    html = (output / "report.html").read_text(encoding="utf-8")

    assert "youden-j" in html
    assert "depth_floor_linked-0001" in html
    assert "0.00469" in html
    assert "step function" in html
    assert "no negative-truth sample lies inside the tested band" in html
    assert "research-decision-profile" in html


def test_the_html_states_which_caller_the_search_varied(tmp_path: Path) -> None:
    """A ``--caller both`` page must say the adVNTR arm was held at its baseline policy."""
    scope = {
        "searched_caller": "kestrel",
        "searched_axes": ["depth_floor_linked", "gg_gate_independent"],
        "advntr_policy": "held-at-baseline",
        "advntr_distinct_executions": 1,
        "note": "The adVNTR arm was replayed at its baseline policy for every candidate (issue #269).",
    }
    output = _write(tmp_path, _document(caller="both", search_scope=scope))
    html = (output / "report.html").read_text(encoding="utf-8")

    assert "depth_floor_linked, gg_gate_independent" in html
    assert "held-at-baseline" in html
    assert "replayed at its baseline policy for every candidate (issue #269)" in html


def test_the_html_escapes_untrusted_text(tmp_path: Path) -> None:
    """Report text reaches the page as text; a reason string is never markup."""
    document = _document()
    document["rejected_candidates"] = [{"axis": "depth_floor_linked", "value": 2.0, "reason": "<b>bad</b> & worse"}]
    output = _write(tmp_path, document)
    html = (output / "report.html").read_text(encoding="utf-8")

    assert "<b>bad</b>" not in html
    assert "&lt;b&gt;bad&lt;/b&gt; &amp; worse" in html


def test_an_infeasible_report_still_renders_and_names_the_constraint(tmp_path: Path) -> None:
    """The failed-selection page is a complete page, not an error placeholder."""
    document = _document(
        status="infeasible",
        successful=False,
        curves=[],
        selection={
            "policy_id": None,
            "reason": "no-candidate-satisfies-constraints",
            "axis": None,
            "value": None,
            "policy_sha256": None,
            "plateau": None,
            "infeasible": {
                "unsatisfiable_constraints": ["min_specificity"],
                "best_achievable": {"min_specificity": 0.5},
                "note": "no tested cutoff satisfies min_specificity",
            },
        },
        profile={"status": "unavailable", "reason": "no cutoff satisfied the declared constraints"},
    )
    output = _write(tmp_path, document)
    html = (output / "report.html").read_text(encoding="utf-8")

    assert "min_specificity" in html
    assert "0.5" in html
    assert "no tested cutoff satisfies min_specificity" in html
    assert (output / "roc-pr-curves.tsv").read_text(encoding="utf-8").splitlines()[0].startswith("axis")


def test_the_writers_refuse_a_document_that_is_not_a_report(tmp_path: Path) -> None:
    """A wrong schema version is a programming error, and fails before any file is written."""
    from vntyper.scripts.calibration_cutoff_report import write_cutoff_reports

    output = tmp_path / "derived"
    output.mkdir(mode=0o700)

    with pytest.raises(ValueError, match="calibration-cutoff-report-v1"):
        write_cutoff_reports(output, _document(schema_version="something-else"))
    assert not any(output.iterdir())
