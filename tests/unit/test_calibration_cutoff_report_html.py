"""The offline report page: what it states, in which order, and what it refuses to draw.

The page has to survive being opened from a private directory on a machine with no
network, so the offline assertions are load-bearing: a single ``<script src=...>`` or web
font turns a published decision into a blank page. The document fixtures are shared with
the writer tests in ``test_calibration_cutoff_report``.
"""

from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Any

import pytest

from tests.unit.test_calibration_cutoff_report import _HELD_OUT, _curve, _document, _write

pytestmark = pytest.mark.unit


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
    from vntyper.scripts.calibration_cutoff_report_html import _UNREJECTABLE_NOTE

    axis = {**_ADVNTR_AXIS, "unrejectable_samples": 2}
    html = (_write(tmp_path, _advntr_document(axes=[axis])) / "report.html").read_text(encoding="utf-8")

    assert _UNREJECTABLE_NOTE.format(n=2) in html
    assert "2 sample(s) have a legacy p-value of exactly 0" in html


def test_the_html_warns_when_a_searched_inventory_was_capped(tmp_path: Path) -> None:
    """A capped search may miss distinct operating points, and the page says so (spec 14.2)."""
    from vntyper.scripts.calibration_cutoff_report_html import _CAPPED_WARNING

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


def test_the_html_escapes_untrusted_text(tmp_path: Path) -> None:
    """Report text reaches the page as text; a reason string is never markup."""
    document = _document()
    document["rejected_candidates"] = [{"axis": "depth_floor_linked", "value": 2.0, "reason": "<b>bad</b> & worse"}]
    output = _write(tmp_path, document)
    html = (output / "report.html").read_text(encoding="utf-8")

    assert "<b>bad</b>" not in html
    assert "&lt;b&gt;bad&lt;/b&gt; &amp; worse" in html


@pytest.mark.parametrize("count", [True, 1.0, "2", -1])
def test_a_malformed_unrejectable_count_is_refused(count: object) -> None:
    """A count that is not a non-negative integer is a malformed record, not a silent omission."""
    from vntyper.scripts.calibration_cutoff_report_html import render_cutoff_report_html

    document = _advntr_document(axes=[{**_ADVNTR_AXIS, "unrejectable_samples": count}])

    with pytest.raises(ValueError, match="advntr_cutoff unrejectable_samples must be a non-negative integer"):
        render_cutoff_report_html(document)


def test_the_unavailable_prefixes_are_rendered_from_their_constants(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_cutoff_document import _CURVE_UNAVAILABLE
    from vntyper.scripts.calibration_cutoff_report_html import _NO_BOUNDARY_SUPPORT, _NO_CURVE

    unavailable = {"axis": "advntr_cutoff", "status": "unavailable", "reason": _CURVE_UNAVAILABLE}
    document = _advntr_document(
        curves=[unavailable],
        boundary_support={"status": "unavailable", "reason": _CURVE_UNAVAILABLE, "warnings": []},
    )
    html = (_write(tmp_path, document) / "report.html").read_text(encoding="utf-8")

    assert f"<div class='warn'>{escape(_NO_CURVE)}{escape(_CURVE_UNAVAILABLE)}.</div>" in html
    assert f"<div class='warn'>{escape(_NO_BOUNDARY_SUPPORT)}{escape(_CURVE_UNAVAILABLE)}.</div>" in html
