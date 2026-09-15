"""Offline reports revalidate frozen length evidence without exposing identities."""

from dataclasses import replace
from importlib import import_module
from pathlib import Path

import pytest

from tests.unit.test_calibration_length_evaluation import (
    evaluation_protocol,
    evaluation_rows,
    fitted_training,
)

pytestmark = pytest.mark.unit


def _evaluation(
    *,
    phase="policy-selection",
    candidates=(("affine-a", "affine-A"),),
    values=(1, 3),
    truths=(60, 160),
):
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    rows = evaluation_rows(values=values, truths=truths, phase=phase)
    protocol, roster = evaluation_protocol(rows, *candidates)
    fitted = fitted_training(protocol)
    result = evaluation.evaluate_length_hypotheses(
        fitted.outcomes,
        fitted.baseline,
        roster,
        protocol,
        rows,
        **({} if phase == "policy-selection" else {"fixed_candidate_id": "affine-a"}),
    )
    return result, roster, protocol, fitted.baseline


def _render(**options):
    report = import_module("vntyper.scripts.calibration_length_report")
    result, roster, protocol, baseline = _evaluation(**options)
    return report.render_length_evaluation(result, roster, protocol, baseline)


def _redigest(result):
    canonical = import_module("vntyper.scripts.canonical_json")
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    return replace(result, sha256=canonical.canonical_sha256(evaluation._result_payload(result)))


def test_report_is_offline_anonymous_and_shows_metrics_bounds_failures_and_selection():
    html = _render()

    for label in (
        "Mean absolute error",
        "Median absolute error",
        "Root mean squared error",
        "Bias",
        "R²",
        "Baseline-relative MAE improvement",
        "Within tolerance",
        "Availability",
        "one-sided 95% lower bound",
        "central 95% paired error-difference interval",
        "Selected candidate",
    ):
        assert label in html
    assert "policy-selection" in html
    assert html.count("<svg ") == 2
    assert "Predicted versus truth" in html and "Residuals" in html
    assert "evaluation-" not in html and "evaluation-group-" not in html
    assert "<script" not in html and "https://" not in html
    assert "Rendering does not authorize model promotion" in html

    escaped = _render(candidates=(("<b>candidate</b>", "affine-A"),))
    assert "&lt;b&gt;candidate&lt;/b&gt;" in escaped
    assert "<b>candidate</b>" not in escaped


def test_residual_plot_keeps_the_zero_reference_when_all_errors_have_one_sign():
    html = _render(truths=(50, 150))
    assert html.count("<line ") == 2  # Scatter identity and residual zero references.


def test_plots_show_numeric_minimum_and_maximum_axis_ticks():
    html = _render(values=(1, 3), truths=(60, 160))
    compact = " ".join(html.split())

    assert '<text class="x-tick" x="45" y="360" text-anchor="middle">60</text>' in compact
    assert '<text class="x-tick" x="365" y="360" text-anchor="middle">160</text>' in compact
    assert '<text class="y-tick" x="40" y="349" text-anchor="end">60</text>' in compact
    assert '<text class="y-tick" x="40" y="29" text-anchor="end">160</text>' in compact


def test_degenerate_plot_range_has_one_centered_numeric_tick_per_axis():
    report = import_module("vntyper.scripts.calibration_length_report")

    assert report._axis_ticks(60.0, 60.0, 45, 365) == [{"position": 205.0, "label": "60"}]
    assert report._axis_ticks(0.0, 0.0, 345, 25) == [{"position": 185.0, "label": "0"}]


@pytest.mark.parametrize("phase", ["validation", "locked-heldout", "development-assessment"])
def test_nonselection_report_has_one_fixed_evaluated_candidate_and_no_selection_claim(phase):
    html = _render(phase=phase)

    assert phase in html
    assert "fixed evaluated candidate" in html
    assert "Fixed evaluated candidate" not in html
    assert "Selection is not applicable in this phase" in html
    assert html.count("<svg ") == 2


def test_nonselection_report_keeps_unselected_siblings_visible_without_plotting_them():
    html = _render(
        phase="validation",
        candidates=(("affine-a", "affine-A"), ("affine-f", "affine-F")),
    )
    assert "not-evaluated" in html
    assert "not_fixed_candidate" in html
    assert "Prediction availability" in html
    assert "<p>not evaluated.</p>" in html
    assert html.count("<svg ") == 2


def test_fit_ineligible_and_unavailable_predictions_remain_visible_without_identifiers():
    result, roster, protocol, baseline = _evaluation(
        candidates=(("affine-a", "affine-A"), ("physical-a", "physical-A")),
        values=(1, 4),
        truths=(60, 210),
    )
    report = import_module("vntyper.scripts.calibration_length_report")
    html = report.render_length_evaluation(result, roster, protocol, baseline)

    assert "fit-ineligible" in html
    assert "physical_model_geometry_evidence_unsupported" in html
    assert "Unavailable predictions" in html
    assert "<p>not evaluated.</p>" in html
    assert "feature_A_out_of_bounds" in html
    assert "evaluation-" not in html

    all_unavailable = _render(values=(4, 5), truths=(210, 260))
    assert "feature_A_out_of_bounds: 2" in all_unavailable
    assert "no_assessable_predictions" in all_unavailable
    assert "<svg " not in all_unavailable


def test_report_recomputes_acceptance_and_rejects_forged_metric_labels():
    result, roster, protocol, baseline = _evaluation()
    candidate = result.candidates[0]
    pooled = candidate.acceptance.pooled
    forged_metrics = replace(pooled.metrics, mae=999.0)
    forged_acceptance = replace(candidate.acceptance, pooled=replace(pooled, metrics=forged_metrics))
    forged_candidate = replace(candidate, acceptance=forged_acceptance)
    forged = _redigest(replace(result, candidates=(forged_candidate,), sha256=""))
    report = import_module("vntyper.scripts.calibration_length_report")

    with pytest.raises(ValueError, match="acceptance differs"):
        report.render_length_evaluation(forged, roster, protocol, baseline)


def test_report_recomputes_selection_and_rejects_a_forged_selected_label():
    result, roster, protocol, baseline = _evaluation()
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    forged_selection = evaluation.LengthSelection("no-feasible-candidate", None, None, ("forged",))
    forged = _redigest(replace(result, selection=forged_selection, sha256=""))
    report = import_module("vntyper.scripts.calibration_length_report")

    with pytest.raises(ValueError, match="selection differs"):
        report.render_length_evaluation(forged, roster, protocol, baseline)


def test_evaluated_candidates_must_share_truth_and_the_frozen_baseline_for_each_member():
    result, roster, protocol, baseline = _evaluation(candidates=(("affine-a", "affine-A"), ("affine-f", "affine-F")))
    report = import_module("vntyper.scripts.calibration_length_report")
    second = result.candidates[1]
    changed_truth = replace(second.predictions[0], truth=second.predictions[0].truth + 1)
    forged = replace(second, predictions=(changed_truth, *second.predictions[1:]))
    changed = _redigest(replace(result, candidates=(result.candidates[0], forged), sha256=""))
    with pytest.raises(ValueError, match="common truth"):
        report.render_length_evaluation(changed, roster, protocol, baseline)

    first = result.candidates[0]
    changed_baseline = replace(first.predictions[0], baseline_prediction=baseline.mean_total_repeat_count + 1)
    forged = replace(first, predictions=(changed_baseline, *first.predictions[1:]))
    changed = _redigest(replace(result, candidates=(forged, result.candidates[1]), sha256=""))
    with pytest.raises(ValueError, match="baseline prediction"):
        report.render_length_evaluation(changed, roster, protocol, baseline)


def test_report_requires_exact_artifact_bindings_and_candidate_roster():
    result, roster, protocol, baseline = _evaluation()
    report = import_module("vntyper.scripts.calibration_length_report")
    with pytest.raises(ValueError, match="baseline binding"):
        report.render_length_evaluation(result, roster, protocol, replace(baseline, sha256="0" * 64))
    with pytest.raises(ValueError, match="candidate roster"):
        changed = _redigest(replace(result, candidates=(), sha256=""))
        report.render_length_evaluation(changed, roster, protocol, baseline)
    with pytest.raises(ValueError, match="exact frozen bindings"):
        changed = _redigest(replace(result, baseline_sha256="0" * 64, sha256=""))
        report.render_length_evaluation(changed, roster, protocol, baseline)
    with pytest.raises(ValueError, match="prediction context"):
        changed = _redigest(replace(result, prediction_context_sha256="0" * 64, sha256=""))
        report.render_length_evaluation(changed, roster, protocol, baseline)


def test_report_shows_sparse_mandatory_strata_and_undefined_metrics_explicitly():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    protocol_contract = import_module("vntyper.scripts.calibration_length_protocol")
    rows = evaluation_rows(truths=(60, 60))
    original, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    raw = protocol_contract.length_protocol_document(original)
    raw["required_strata"] = ["all", "sparse"]
    protocol = protocol_contract.decode_length_protocol(raw)
    fitted = fitted_training(protocol)
    result = evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, rows)
    baseline = fitted.baseline
    report = import_module("vntyper.scripts.calibration_length_report")

    html = report.render_length_evaluation(result, roster, protocol, baseline)
    assert "Mandatory stratum: all" in html
    assert "Mandatory stratum: sparse" in html
    assert "missing_required_stratum" in html
    assert "undefined: no metrics for this frozen population" in html
    assert "R²</th>" in html and "<td>undefined</td>" in html


@pytest.mark.parametrize(
    "mode",
    ["duplicate", "missing", "availability", "feature", "roster", "model", "evaluated", "fit", "status"],
)
def test_report_rejects_inconsistent_prediction_and_candidate_states(mode):
    result, roster, protocol, baseline = _evaluation()
    candidate = result.candidates[0]
    if mode == "duplicate":
        predictions = (candidate.predictions[0], candidate.predictions[0])
        forged = replace(candidate, predictions=predictions)
    elif mode == "missing":
        forged = replace(candidate, predictions=candidate.predictions[:1])
    elif mode == "availability":
        prediction = replace(candidate.predictions[0], availability_reasons=("unavailable",))
        forged = replace(candidate, predictions=(prediction, *candidate.predictions[1:]))
    elif mode == "feature":
        prediction = replace(candidate.predictions[0], features_sha256="bad")
        forged = replace(candidate, predictions=(prediction, *candidate.predictions[1:]))
    elif mode == "roster":
        forged = replace(candidate, candidate_id="other")
    elif mode == "model":
        forged = replace(candidate, model_sha256="bad")
    elif mode == "evaluated":
        forged = replace(candidate, reasons=("bad",))
    elif mode == "fit":
        forged = replace(candidate, status="fit-ineligible")
    else:
        forged = replace(candidate, status="unknown")
    changed = _redigest(replace(result, candidates=(forged,), sha256=""))
    report = import_module("vntyper.scripts.calibration_length_report")
    with pytest.raises(ValueError):
        report.render_length_evaluation(changed, roster, protocol, baseline)


def test_length_report_template_is_explicitly_packaged():
    root = Path(__file__).resolve().parents[2]
    manifest = (root / "MANIFEST.in").read_text(encoding="utf-8")
    project = (root / "pyproject.toml").read_text(encoding="utf-8")
    template = (root / "vntyper/templates/calibration_length_report.html").read_text(encoding="utf-8")
    assert "include vntyper/templates/calibration_length_report.html" in manifest
    assert '"templates/calibration_length_report.html"' in project
    assert all(len(line.rstrip()) <= 120 for line in template.splitlines())
