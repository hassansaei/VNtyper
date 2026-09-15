"""Reports use complete populations and show only frozen validation points."""

from dataclasses import replace
from fractions import Fraction

import pytest

from vntyper.scripts.calibration_caller_curves import CallerOperatingPoint, build_caller_curves
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_report import CallerReportCandidate, render_caller_comparison
from vntyper.scripts.calibration_caller_roster import decode_caller_eligible_roster

pytestmark = pytest.mark.unit


def _inputs():
    rows = (
        CallerObservation("private-key-a", "g1", True, ("v",), True, ("v",), ()),
        CallerObservation("private-key-b", "g2", False, (), False, (), ()),
    )
    roster = decode_caller_eligible_roster(
        [
            {"key": rows[0].key, "group_key": "g1", "strata": ["nominal"]},
            {"key": rows[1].key, "group_key": "g2", "strata": ["nominal"]},
        ]
    )
    candidates = (
        CallerReportCandidate("a" * 64, rows, ()),
        CallerReportCandidate(
            "b" * 64, (replace(rows[0], called_positive=None, called_variants=()), rows[1]), ("sensitivity_gate",)
        ),
    )
    return candidates, roster


def _render(candidates=None, roster=None, **changes):
    default_candidates, default_roster = _inputs()
    options = {
        "phase": "policy-selection",
        "baseline_id": "a" * 64,
        "selected_id": None,
        "protocol_sha256": "c" * 64,
        "evidence_sha256": "d" * 64,
        "required_strata": ("nominal", "sparse"),
        "curves": (),
    }
    options.update(changes)
    return render_caller_comparison(candidates or default_candidates, roster or default_roster, **options)


def test_report_is_offline_escaped_and_contains_counts_intervals_and_rejections():
    html = _render()
    assert "policy-selection" in html
    assert "sensitivity_gate" in html
    assert "Sensitivity" in html and "Specificity" in html and "Precision" in html
    assert "No-call rate" in html and "Wrong tier-A identity groups" in html
    assert "0/1" in html and "1/2" in html
    assert "central 95%" in html and "one-sided 95%" in html
    assert "insufficient-evidence" in html
    assert "No candidate selected" in html
    assert "private-key" not in html
    assert "<script" not in html and "https://" not in html
    assert "lexicographic-safety-v1" in html


def test_report_rejects_an_omitted_no_call_before_calculating_denominators():
    candidates, roster = _inputs()
    incomplete = replace(candidates[1], observations=candidates[1].observations[1:])
    with pytest.raises(ValueError, match="roster"):
        _render((candidates[0], incomplete), roster)


def test_report_rejects_truth_changes_between_baseline_and_candidate():
    candidates, roster = _inputs()
    altered = replace(
        candidates[1],
        observations=(
            replace(candidates[1].observations[0], truth_variants=("different",)),
            candidates[1].observations[1],
        ),
    )
    with pytest.raises(ValueError, match="truth"):
        _render((candidates[0], altered), roster)


def test_validation_shows_only_frozen_baseline_and_candidate_and_no_search():
    candidates, roster = _inputs()
    candidates = (candidates[0], replace(candidates[1], rejection_reasons=()))
    html = _render(candidates, roster, phase="validation", selected_id="b" * 64)
    assert "Frozen candidate" in html
    with pytest.raises(ValueError, match="exactly"):
        _render(
            (*candidates, replace(candidates[1], candidate_id="e" * 64)),
            roster,
            phase="validation",
            selected_id="b" * 64,
        )
    with pytest.raises(ValueError, match="frozen"):
        _render(candidates, roster, phase="validation")


def test_labels_and_reasons_are_html_escaped():
    candidates, roster = _inputs()
    candidates = (candidates[0], replace(candidates[1], rejection_reasons=("<b>unsafe</b>",)))
    assert "&lt;b&gt;unsafe&lt;/b&gt;" in _render(candidates, roster)


@pytest.mark.parametrize(
    "changes",
    [
        {"phase": "training"},
        {"baseline_id": "e" * 64},
        {"selected_id": "e" * 64},
        {"selected_id": "a" * 64},
        {"protocol_sha256": "bad"},
        {"evidence_sha256": "bad"},
        {"required_strata": ()},
        {"required_strata": ("nominal", "nominal")},
    ],
)
def test_invalid_comparison_context_fails(changes):
    with pytest.raises(ValueError):
        _render(**changes)


def test_selected_candidate_cannot_have_selection_rejection_reasons():
    with pytest.raises(ValueError, match="rejected"):
        _render(selected_id="b" * 64)


def _curves():
    candidates, roster = _inputs()
    strict_rows = tuple(replace(row, called_positive=False, called_variants=()) for row in candidates[0].observations)
    loose_rows = tuple(replace(row, called_positive=True, called_variants=("v",)) for row in strict_rows)
    candidates = (replace(candidates[0], observations=strict_rows), replace(candidates[1], observations=loose_rows))
    curve = build_caller_curves(
        (
            CallerOperatingPoint("a" * 64, Fraction(9, 10), "e" * 64, strict_rows),
            CallerOperatingPoint("b" * 64, Fraction(1, 10), "e" * 64, loose_rows),
        ),
        comparison=">=",
        phase="policy-selection",
    )
    return candidates, roster, curve


def test_report_draws_revalidated_offline_roc_and_precision_recall_curves():
    candidates, roster, curve = _curves()
    html = _render(candidates, roster, curves=(curve,))
    assert html.count("<svg ") == 2
    assert html.count("<circle ") == 3  # Precision at no predicted positives is undefined.
    assert "ROC operating curve" in html
    assert "Precision–recall operating curve" in html
    assert "9/10" in html and "1/10" in html


@pytest.mark.parametrize("mode", ["mutable", "untyped", "empty", "point", "absent", "forged", "validation"])
def test_report_refuses_invalid_or_unauthorized_curve_data(mode):
    candidates, roster, curve = _curves()
    curves = (curve,)
    options = {}
    if mode == "mutable":
        curves = [curve]
    elif mode == "untyped":
        curves = (object(),)
    elif mode == "empty":
        curves = (replace(curve, points=()),)
    elif mode == "point":
        curves = (replace(curve, points=(object(),)),)
    elif mode == "absent":
        curves = (replace(curve, points=(replace(curve.points[0], candidate_id="f" * 64),)),)
    elif mode == "forged":
        curves = (replace(curve, eligible_count=999),)
    else:
        options = {"phase": "validation", "selected_id": "b" * 64}
    with pytest.raises(ValueError):
        _render(candidates, roster, curves=curves, **options)


@pytest.mark.parametrize("mode", ["duplicate", "untyped", "observations", "reasons", "id"])
def test_report_refuses_invalid_candidate_contracts(mode):
    candidates, roster = _inputs()
    if mode == "duplicate":
        candidates = (candidates[0], candidates[0])
    elif mode == "untyped":
        candidates = (object(),)
    elif mode == "observations":
        candidates = (replace(candidates[0], observations=list(candidates[0].observations)),)
    elif mode == "reasons":
        candidates = (replace(candidates[0], rejection_reasons=["reason"]),)
    else:
        candidates = (replace(candidates[0], candidate_id="bad"),)
    with pytest.raises(ValueError):
        _render(candidates, roster)


def test_failed_validation_keeps_frozen_candidate_rejection_reasons_visible():
    html = _render(phase="validation", selected_id="b" * 64)
    assert "Frozen candidate" in html and "sensitivity_gate" in html


def test_selected_and_unselected_labels_are_distinct_in_selection_report():
    candidates, roster = _inputs()
    candidates = (candidates[0], replace(candidates[1], rejection_reasons=()))
    assert "Selected candidate" in _render(candidates, roster, selected_id="b" * 64)
    assert "Unselected" in _render(candidates, roster)
