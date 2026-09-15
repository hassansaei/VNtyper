"""Pure model evaluation and policy-selection tests for total length."""

from dataclasses import FrozenInstanceError, replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_length import metadata_for, roster_for, training_row
from tests.unit.test_calibration_length_protocol import protocol_document
from tests.unit.test_length_estimation import model_for

pytestmark = pytest.mark.unit


def evaluation_rows(values=(1, 3), truths=(60, 160), *, phase="policy-selection", feature="A"):
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    return tuple(
        evaluation.LengthEvaluationRow(
            key=f"evaluation-{index}",
            group_key=f"evaluation-group-{index}",
            phase=phase,
            features=training_row(f"evaluation-{index}", value, truth, feature=feature).features,
            truth_boundary_definition="complete-core-plus-invariant-units-v1",
            total_truth_repeat_count=truth,
            evidence_domain="synthetic",
        )
        for index, (value, truth) in enumerate(zip(values, truths, strict=True), start=1)
    )


def evaluation_roster(rows):
    metrics = import_module("vntyper.scripts.calibration_length_metrics")
    return metrics.decode_length_eligible_roster(
        [
            {"key": row.key, "group_key": row.group_key, "strata": ["all"]}
            for row in sorted(rows, key=lambda item: item.group_key)
        ]
    )


def evaluation_protocol(rows, *candidates):
    protocol = import_module("vntyper.scripts.calibration_length_protocol")
    roster = evaluation_roster(rows)
    raw = protocol_document()
    raw["candidate_grid"] = sorted(
        ({"candidate_id": candidate_id, "model_kind": model_kind} for candidate_id, model_kind in candidates),
        key=lambda item: item["candidate_id"],
    )
    raw["maximum_candidate_count"] = 4
    raw["required_strata"] = ["all"]
    raw["eligible_roster_sha256"] = roster.sha256
    raw["acceptance"].update(
        minimum_independent_count=2,
        maximum_mae=2,
        minimum_relative_mae_improvement=0.01,
        minimum_tolerance_lower_bound=0.01,
        minimum_availability_lower_bound=0.01,
    )
    return protocol.decode_length_protocol(raw), roster


def fitted_training(protocol):
    training = import_module("vntyper.scripts.calibration_length")
    rows = tuple(training_row(str(index), x, 10 + 50 * x) for index, x in enumerate((1, 2, 3), start=1))
    roster = roster_for(rows)
    metadata = metadata_for(rows, roster)
    return training.fit_length_hypotheses(rows, roster, protocol, metadata)


def test_synthetic_fit_predict_acceptance_and_selection_uses_the_frozen_model_without_refit(monkeypatch):
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    training = import_module("vntyper.scripts.calibration_length")
    rows = evaluation_rows()
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    fitted = fitted_training(protocol)
    model_sha256 = fitted.outcomes[0].model.sha256
    monkeypatch.setattr(training, "fit_length_hypotheses", lambda *_args, **_kwargs: pytest.fail("evaluation refit"))

    result = evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, rows)

    assert result.phase == "policy-selection"
    assert result.selection.status == "selected"
    assert result.selection.selected_candidate_id == "affine-a"
    assert result.selection.selected_model_sha256 == model_sha256
    assert result.candidates[0].model_sha256 == model_sha256
    assert [item.prediction for item in result.candidates[0].predictions] == pytest.approx([60, 160])
    assert {item.baseline_prediction for item in result.candidates[0].predictions} == {
        fitted.baseline.mean_total_repeat_count
    }
    assert result.candidates[0].acceptance.status == "passed"
    assert result.baseline_sha256 == fitted.baseline.sha256
    assert result.eligible_roster_sha256 == roster.sha256
    assert result.evaluation_rows_sha256
    assert evaluation.length_evaluation_document(result)["sha256"] == result.sha256
    reordered = evaluation.evaluate_length_hypotheses(
        fitted.outcomes, fitted.baseline, roster, protocol, tuple(reversed(rows))
    )
    assert reordered == result
    with pytest.raises(FrozenInstanceError):
        result.phase = "development"


def test_inclusive_one_repeat_tie_uses_stable_candidate_id_for_equal_complexity():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    training = import_module("vntyper.scripts.calibration_length")
    rows = evaluation_rows()
    protocol, roster = evaluation_protocol(rows, ("a-affine-f", "affine-F"), ("z-affine-a", "affine-A"))
    fitted = fitted_training(evaluation_protocol(rows, ("z-affine-a", "affine-A"))[0])
    baseline = fitted.baseline
    measured = rows[0].features
    model_a = model_for(
        measured,
        "A",
        intercept=10,
        coefficients=[50],
        study_sha256=baseline.study_sha256,
        training_evidence_sha256=baseline.training_evidence_sha256,
    )
    model_f = model_for(
        measured,
        "F",
        intercept=-39,
        coefficients=[100],
        study_sha256=baseline.study_sha256,
        training_evidence_sha256=baseline.training_evidence_sha256,
    )
    outcomes = (
        training.LengthFitOutcome("a-affine-f", "affine-F", "fitted", (), model_f),
        training.LengthFitOutcome("z-affine-a", "affine-A", "fitted", (), model_a),
    )

    result = evaluation.evaluate_length_hypotheses(outcomes, baseline, roster, protocol, rows)

    metrics = {item.candidate_id: item.acceptance.pooled.metrics for item in result.candidates}
    assert metrics["z-affine-a"].mae == 0
    assert metrics["a-affine-f"].mae == 1
    assert result.selection.selected_candidate_id == "a-affine-f"


@pytest.mark.parametrize("phase", ["validation", "locked-heldout", "development-assessment"])
def test_nonselection_phases_require_and_evaluate_only_the_fixed_candidate(monkeypatch, phase):
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    estimation = import_module("vntyper.scripts.length_estimation")
    rows = evaluation_rows(phase=phase)
    protocol, roster = evaluation_protocol(rows, ("a-affine-f", "affine-F"), ("z-affine-a", "affine-A"))
    fitted = fitted_training(evaluation_protocol(rows, ("z-affine-a", "affine-A"))[0])
    measured = rows[0].features
    baseline = fitted.baseline
    outcomes = (
        import_module("vntyper.scripts.calibration_length").LengthFitOutcome(
            "a-affine-f",
            "affine-F",
            "fitted",
            (),
            model_for(
                measured,
                "F",
                study_sha256=baseline.study_sha256,
                training_evidence_sha256=baseline.training_evidence_sha256,
            ),
        ),
        fitted.outcomes[0],
    )
    original = estimation.estimate_total_repeats
    called = []

    def recording_estimate(features, model, *, evidence_domain=None):
        called.append(model.sha256)
        return original(features, model, evidence_domain=evidence_domain)

    monkeypatch.setattr(evaluation, "estimate_total_repeats", recording_estimate)
    with pytest.raises(ValueError, match="fixed_candidate_id"):
        evaluation.evaluate_length_hypotheses(outcomes, baseline, roster, protocol, rows)
    result = evaluation.evaluate_length_hypotheses(
        outcomes, baseline, roster, protocol, rows, fixed_candidate_id="z-affine-a"
    )
    assert len(called) == len(rows)
    assert set(called) == {fitted.outcomes[0].model.sha256}
    assert result.candidates[0].status == "not-evaluated"
    assert result.candidates[1].status == "evaluated"
    assert result.selection.status == "not-applicable"
    assert result.selection.selected_candidate_id is None


def test_old_development_phase_name_is_rejected():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    rows = evaluation_rows(phase="development")
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    fitted = fitted_training(protocol)
    with pytest.raises(ValueError, match="phase"):
        evaluation.evaluate_length_hypotheses(
            fitted.outcomes,
            fitted.baseline,
            roster,
            protocol,
            rows,
            fixed_candidate_id="affine-a",
        )


def test_full_protocol_fit_outcome_roster_and_fit_state_are_strict():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    training = import_module("vntyper.scripts.calibration_length")
    rows = evaluation_rows()
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"), ("physical-a", "physical-A"))
    fitted = fitted_training(protocol)
    result = evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, rows)
    assert result.candidates[1].status == "fit-ineligible"
    assert result.candidates[1].reasons == ("physical_model_geometry_evidence_unsupported",)
    with pytest.raises(ValueError, match="exactly match"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes[:1], fitted.baseline, roster, protocol, rows)
    malformed = replace(fitted.outcomes[0], reasons=("unexpected",))
    with pytest.raises(ValueError, match="fitted outcome"):
        evaluation.evaluate_length_hypotheses((malformed, fitted.outcomes[1]), fitted.baseline, roster, protocol, rows)
    malformed = training.LengthFitOutcome("affine-a", "affine-A", "ineligible", (), None)
    with pytest.raises(ValueError, match="ineligible outcome"):
        evaluation.evaluate_length_hypotheses((malformed, fitted.outcomes[1]), fitted.baseline, roster, protocol, rows)


def test_rows_require_exact_roster_feature_identity_integral_truth_and_one_allowed_phase():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    rows = evaluation_rows()
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    fitted = fitted_training(protocol)
    with pytest.raises(ValueError, match="roster exactly"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, rows[:1])
    swapped = (replace(rows[0], features=rows[1].features), rows[1])
    with pytest.raises(ValueError, match="manifest key"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, swapped)
    fractional = (replace(rows[0], total_truth_repeat_count=60.5), rows[1])
    with pytest.raises(ValueError, match="integral exact repeat count"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, fractional)
    mixed = (rows[0], replace(rows[1], phase="development-assessment"))
    with pytest.raises(ValueError, match="one common phase"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, mixed)
    invalid = (replace(rows[0], phase="training"), replace(rows[1], phase="training"))
    with pytest.raises(ValueError, match="phase"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, invalid)


def test_model_baseline_protocol_and_roster_bindings_are_revalidated():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    model_contract = import_module("vntyper.scripts.length_model")
    rows = evaluation_rows()
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    fitted = fitted_training(protocol)
    with pytest.raises(ValueError, match="baseline"):
        evaluation.evaluate_length_hypotheses(
            fitted.outcomes, replace(fitted.baseline, sha256="0" * 64), roster, protocol, rows
        )
    other_protocol, _ = evaluation_protocol(rows, ("affine-a", "affine-A"))
    raw = import_module("vntyper.scripts.calibration_length_protocol").length_protocol_document(other_protocol)
    raw["qc"]["minimum_denominator_mean_depth"] = 11
    changed_protocol = import_module("vntyper.scripts.calibration_length_protocol").decode_length_protocol(raw)
    with pytest.raises(ValueError, match="QC"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, changed_protocol, rows)
    model_raw = model_contract.encode_length_model(fitted.outcomes[0].model)
    model_raw["study_sha256"] = "f" * 64
    mismatched_model = model_contract.decode_length_model(model_raw)
    mismatched_outcome = replace(fitted.outcomes[0], model=mismatched_model)
    with pytest.raises(ValueError, match="study or training evidence"):
        evaluation.evaluate_length_hypotheses((mismatched_outcome,), fitted.baseline, roster, protocol, rows)
    with pytest.raises(ValueError, match="protocol roster"):
        evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster_for(rows), protocol, rows)


def test_unavailable_predictions_remain_in_denominator_and_input_content_is_digest_bound():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    rows = evaluation_rows(values=(1, 4), truths=(60, 210))
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    fitted = fitted_training(protocol)
    result = evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, rows)
    candidate = result.candidates[0]
    assert candidate.predictions[1].prediction is None
    assert candidate.predictions[1].availability_reasons == ("feature_A_out_of_bounds",)
    assert candidate.acceptance.pooled.metrics.eligible_count == 2
    assert candidate.acceptance.pooled.metrics.assessable_count == 1
    changed = (rows[0], replace(rows[1], total_truth_repeat_count=211))
    changed_result = evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, changed)
    assert changed_result.evaluation_rows_sha256 != result.evaluation_rows_sha256
    assert changed_result.sha256 != result.sha256


def test_policy_selection_rejects_a_fixed_candidate_and_result_digest_cannot_be_forged():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    rows = evaluation_rows()
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    fitted = fitted_training(protocol)
    with pytest.raises(ValueError, match="must not receive"):
        evaluation.evaluate_length_hypotheses(
            fitted.outcomes, fitted.baseline, roster, protocol, rows, fixed_candidate_id="affine-a"
        )
    result = evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, rows)
    with pytest.raises(ValueError, match="digest"):
        evaluation.length_evaluation_document(replace(result, sha256="0" * 64))
    mutable_candidate = replace(result.candidates[0], predictions=list(result.candidates[0].predictions))
    with pytest.raises(ValueError, match="immutable"):
        evaluation.length_evaluation_document(replace(result, candidates=(mutable_candidate,)))
    mutable_acceptance = replace(result.candidates[0].acceptance, strata=dict(result.candidates[0].acceptance.strata))
    mutable_candidate = replace(result.candidates[0], acceptance=mutable_acceptance)
    with pytest.raises(ValueError, match="immutable"):
        evaluation.length_evaluation_document(replace(result, candidates=(mutable_candidate,)))


def test_no_feasible_candidate_is_explicit_and_cannot_select_fit_failure():
    evaluation = import_module("vntyper.scripts.calibration_length_evaluation")
    rows = evaluation_rows()
    protocol, roster = evaluation_protocol(rows, ("physical-a", "physical-A"))
    fitted = fitted_training(protocol)
    result = evaluation.evaluate_length_hypotheses(fitted.outcomes, fitted.baseline, roster, protocol, rows)
    assert result.selection.status == "no-feasible-candidate"
    assert result.selection.reasons == ("no_candidate_passed_acceptance",)
    validation_rows = tuple(replace(row, phase="validation") for row in rows)
    with pytest.raises(ValueError, match="fitted model"):
        evaluation.evaluate_length_hypotheses(
            fitted.outcomes,
            fitted.baseline,
            roster,
            protocol,
            validation_rows,
            fixed_candidate_id="physical-a",
        )
