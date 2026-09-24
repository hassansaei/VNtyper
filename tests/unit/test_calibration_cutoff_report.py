"""The cutoff report writers: canonical JSON, joinable TSVs and an offline page.

The page has to survive being opened from a private directory on a machine with no
network, so every one of these assertions about ``report.html`` is load-bearing: a single
``<script src=...>`` or web font turns a published decision into a blank page.
"""

from __future__ import annotations

import json
import stat
from html import escape
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


#: Pooled held-out metrics in the exact shape ``evaluate_cutoff_arms`` publishes: held-out
#: specificity 22/27 = 0.815 sits below a requested 0.85 floor the full data met.
_HELD_OUT: dict[str, Any] = {
    "counts": {
        "true_positives": 50,
        "false_negatives": 5,
        "true_negatives": 22,
        "false_positives": 5,
        "no_calls": 0,
        "sensitivity": 50 / 55,
        "specificity": 22 / 27,
    },
    "exact": {
        "sensitivity": {"events": 50, "total": 55, "estimate": 50 / 55, "lower": 0.8005, "upper": 0.9698},
        "specificity": {"events": 22, "total": 27, "estimate": 22 / 27, "lower": 0.6192, "upper": 0.9370},
    },
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
            "searched_callers": ["kestrel"],
            "searched_axes": ["depth_floor_linked"],
            "advntr_policy": "not-evaluated",
            "advntr_distinct_executions": None,
            "advntr_probe_executions": None,
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
            "advntr": None,
        },
        "replay_consistency": None,
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
                    "admissible_candidates": 7,
                    "fallback_reason": None,
                    "selection": {"policy_id": "depth_floor_linked-0001", "reason": "selected"},
                }
            ],
            "fold_admissibility": "training-observed-breakpoints",
            "held_out": _HELD_OUT,
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
        "searched_callers": ["kestrel"],
        "searched_axes": ["depth_floor_linked", "gg_gate_independent"],
        "advntr_policy": "held-at-baseline",
        "advntr_distinct_executions": 1,
        "advntr_probe_executions": None,
        "note": "The adVNTR arm was replayed at its baseline policy for every candidate.",
    }
    output = _write(tmp_path, _document(caller="both", search_scope=scope))
    html = (output / "report.html").read_text(encoding="utf-8")

    assert "<dt>searched callers</dt><dd>kestrel</dd>" in html
    assert "depth_floor_linked, gg_gate_independent" in html
    assert "held-at-baseline" in html
    assert "replayed at its baseline policy for every candidate." in html
    assert "<dt>adVNTR baseline parity</dt><dd>not replayed</dd>" in html
    assert "<dt>adVNTR replay consistency (candidates checked)</dt><dd>not checked</dd>" in html


#: A searched adVNTR cutoff axis as ``axis_document`` publishes it, trimmed to what the page reads.
_ADVNTR_AXIS: dict[str, Any] = {
    "axis": "advntr_cutoff",
    "caller": "advntr",
    "comparator": "<",
    "breakpoint_completeness": "complete",
    "unrejectable_samples": 0,
}


def _advntr_document(**overrides: Any) -> dict[str, Any]:
    """A ``--caller advntr`` report whose adVNTR cutoff axis was searched natively."""
    from vntyper.scripts.calibration_cutoff_document import _SCOPE_NOTES

    document = _document(caller="advntr")
    document["search_scope"] = {
        "searched_callers": ["advntr"],
        "searched_axes": ["advntr_cutoff"],
        "advntr_policy": "searched",
        "advntr_distinct_executions": 7,
        "advntr_probe_executions": 2,
        "note": _SCOPE_NOTES["searched"],
    }
    document["axes"] = [dict(_ADVNTR_AXIS)]
    document["baseline_parity"] = {
        **document["baseline_parity"],
        "advntr": {"proven": True, "sample_count": 6, "mismatches": []},
    }
    document["replay_consistency"] = {"checked_candidates": 7, "mismatches": []}
    document["provenance"] = {
        **document["provenance"],
        "advntr": {
            "sha256": "4" * 64,
            "probe_sha256": "5" * 64,
            "tool_identity": {"package_version": "2.4.0", "build_id": "build-synthetic-0001"},
            "probe_seconds": 1.25,
            "main_seconds": 12.5,
        },
    }
    document.update(overrides)
    return document


def test_the_html_states_a_searched_advntr_scope_and_its_parity(tmp_path: Path) -> None:
    """A searched adVNTR axis is named with its note, its parity and its replay-consistency count."""
    from vntyper.scripts.calibration_cutoff_document import _SCOPE_NOTES

    html = (_write(tmp_path, _advntr_document()) / "report.html").read_text(encoding="utf-8")

    assert "<dt>searched callers</dt><dd>advntr</dd>" in html
    assert "<dt>adVNTR policy</dt><dd>searched</dd>" in html
    assert _SCOPE_NOTES["searched"] in html
    assert "<dt>adVNTR baseline parity</dt><dd>True</dd>" in html
    assert "<dt>adVNTR replay consistency (candidates checked)</dt><dd>7</dd>" in html
    assert html.index("Held-out performance (cross-validated)") < html.index("Every tested cutoff")
    assert "unrejectable" not in html and "no admissible cutoff can reject" not in html
    assert "subsampled by --max-breakpoints" not in html


def test_the_html_shows_the_advntr_tool_identity_and_both_durations(tmp_path: Path) -> None:
    html = (_write(tmp_path, _advntr_document()) / "report.html").read_text(encoding="utf-8")

    assert "<dt>adVNTR package version</dt><dd>2.4.0</dd>" in html
    assert "<dt>adVNTR build id</dt><dd>build-synthetic-0001</dd>" in html
    assert "<dt>adVNTR probe grid seconds</dt><dd>1.25</dd>" in html
    assert "<dt>adVNTR candidate grid seconds</dt><dd>12.5</dd>" in html


def test_a_kestrel_page_carries_no_advntr_tool_facts(tmp_path: Path) -> None:
    html = (_write(tmp_path, _document()) / "report.html").read_text(encoding="utf-8")

    assert "adVNTR package version" not in html
    assert "adVNTR probe grid seconds" not in html


def test_the_html_names_samples_no_admissible_cutoff_can_reject(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_cutoff_report import _UNREJECTABLE_NOTE

    axis = {**_ADVNTR_AXIS, "unrejectable_samples": 2}
    html = (_write(tmp_path, _advntr_document(axes=[axis])) / "report.html").read_text(encoding="utf-8")

    assert _UNREJECTABLE_NOTE.format(n=2) in html
    assert "2 sample(s) have a legacy p-value of exactly 0" in html


def test_the_html_warns_when_a_searched_inventory_was_capped(tmp_path: Path) -> None:
    """A capped search may miss distinct operating points, and the page says so (spec 14.2)."""
    from vntyper.scripts.calibration_cutoff_report import _CAPPED_WARNING

    complete = (_write(tmp_path / "complete", _advntr_document()) / "report.html").read_text(encoding="utf-8")
    capped_axis = {**_ADVNTR_AXIS, "breakpoint_completeness": "capped-subsample"}
    capped = (_write(tmp_path / "capped", _advntr_document(axes=[capped_axis])) / "report.html").read_text(
        encoding="utf-8"
    )

    assert _CAPPED_WARNING not in complete
    assert f"<div class='warn'>{_CAPPED_WARNING}</div>" in capped


def test_an_unavailable_union_curve_renders_its_reason_instead_of_a_plot(tmp_path: Path) -> None:
    """A curve whose no-call set changes is published as a reason, never as an SVG or a TSV row (spec 14.6)."""
    from vntyper.scripts.calibration_cutoff_document import _CURVE_UNAVAILABLE

    unavailable = {"axis": "advntr_cutoff", "status": "unavailable", "reason": _CURVE_UNAVAILABLE}
    document = _advntr_document(
        caller="both",
        curves=[{**_curve("depth_floor_linked"), "status": "available"}, unavailable],
        boundary_support={"status": "unavailable", "reason": _CURVE_UNAVAILABLE, "warnings": []},
    )
    output = _write(tmp_path, document)
    html = (output / "report.html").read_text(encoding="utf-8")
    curve_rows = (output / "roc-pr-curves.tsv").read_text(encoding="utf-8").splitlines()

    assert html.count(escape(_CURVE_UNAVAILABLE)) == 2  # the axis section and the boundary support
    assert "<h2>Axis advntr_cutoff</h2>" in html
    assert "ROC - advntr_cutoff" not in html
    assert "ROC - depth_floor_linked" in html
    assert "positives inside the tested band" not in html
    assert {line.split("\t")[0] for line in curve_rows[1:]} == {"depth_floor_linked"}


def test_the_html_shows_held_out_performance_before_any_full_data_number(tmp_path: Path) -> None:
    """The validated-looking full-data table must not be the first performance a reader sees."""
    objective = {"objective": "max-sensitivity-at-specificity", "min_sensitivity": None, "min_specificity": 0.85}
    output = _write(tmp_path, _document(objective=objective))
    html = (output / "report.html").read_text(encoding="utf-8")

    held = html.index("Held-out performance (cross-validated)")
    assert held < html.index("<h2>Selection</h2>") < html.index("Every tested cutoff (descriptive")
    assert "50/55 = 0.909 (95% CI 0.800-0.970)" in html
    assert "22/27 = 0.815 (95% CI 0.619-0.937)" in html
    assert "exclude selection uncertainty" in html
    assert "Held-out specificity 0.815 is below the requested floor 0.85" in html
    assert "Descriptive searched-cohort points" in html
    assert "not validated performance" in html
    assert "training-observed-breakpoints" in html
    section = html[held : html.index("<h2>Selection</h2>")]
    assert "depth_floor_linked-0001" in section and "<td>7</td>" in section


def test_the_html_says_so_when_no_held_out_estimate_exists(tmp_path: Path) -> None:
    document = _document()
    document["evaluation"] = {**document["evaluation"], "held_out": None, "status_reason": "insufficient-groups"}
    html = (_write(tmp_path, document) / "report.html").read_text(encoding="utf-8")

    assert "No held-out estimate is available: insufficient-groups" in html


def test_a_held_out_rate_without_eligible_samples_is_shown_as_undefined(tmp_path: Path) -> None:
    document = _document()
    held = {**_HELD_OUT, "exact": {**_HELD_OUT["exact"], "specificity": {"estimate": None}}}
    document["evaluation"] = {**document["evaluation"], "held_out": held}
    html = (_write(tmp_path, document) / "report.html").read_text(encoding="utf-8")

    assert "undefined (no eligible samples)" in html


def test_the_fold_tsv_records_how_many_candidates_each_fold_could_select(tmp_path: Path) -> None:
    output = _write(tmp_path, _document())
    lines = (output / "folds.tsv").read_text(encoding="utf-8").splitlines()

    assert lines[0].split("\t")[-1] == "admissible_candidates"
    assert lines[1].split("\t")[-1] == "7"


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


def test_an_advntr_axis_without_its_search_record_is_refused() -> None:
    """``unrejectable_samples`` comes from the probe search; an adVNTR axis without one is a wiring defect."""
    from types import SimpleNamespace

    from tests.unit.test_calibration_cutoff_optimize import _advntr_baseline_policy
    from vntyper.scripts.calibration_cutoff_axes import ADVNTR_CUTOFF, declared_axis
    from vntyper.scripts.calibration_cutoff_document import _axis_documents

    axis = declared_axis(ADVNTR_CUTOFF, [0.001, 0.004], baseline=_advntr_baseline_policy())
    inputs = SimpleNamespace(derived=[(axis, ())], advntr_search=None)

    with pytest.raises(ValueError, match="adVNTR axis advntr_cutoff has no adVNTR search record"):
        _axis_documents(inputs)  # type: ignore[arg-type]
