"""Reference length acceptance cannot hide bias, missingness or sparse strata."""

from dataclasses import replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_length_protocol import protocol_document
from vntyper.scripts.calibration_length_metrics import LengthObservation, decode_length_eligible_roster
from vntyper.scripts.calibration_length_protocol import decode_length_protocol

pytestmark = pytest.mark.unit


def population(stratum="nominal", count=60, error=0, *, memberships=None):
    if memberships is None:
        memberships = [stratum]
    rows = tuple(
        LengthObservation(f"{stratum}-sample-{i}", f"{stratum}-group-{i}", 100 + i, 100 + i + error, 90)
        for i in range(count)
    )
    members = [{"key": row.key, "group_key": row.group_key, "strata": sorted(memberships)} for row in rows]
    return rows, members


def evidence(*populations):
    rows = tuple(row for population_rows, _ in populations for row in population_rows)
    members = [member for _, population_members in populations for member in population_members]
    return tuple(sorted(rows, key=lambda row: row.group_key)), decode_length_eligible_roster(
        sorted(members, key=lambda member: member["group_key"])
    )


def protocol_for(roster, *strata, changes=None):
    raw = protocol_document()
    raw["required_strata"] = sorted(strata)
    raw["eligible_roster_sha256"] = roster.sha256
    if changes:
        changes(raw)
    return decode_length_protocol(raw)


def context_for(protocol):
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    return m.decode_length_prediction_context(
        {
            "schema_version": "length-prediction-context-v1",
            "protocol_sha256": protocol.sha256,
            "qc_sha256": protocol.qc_sha256,
        }
    )


def evaluate(rows, roster, protocol):
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    return m.evaluate_length_acceptance(rows, roster, context_for(protocol), protocol)


def test_perfect_fixed_predictions_pass_and_bind_the_frozen_inputs():
    rows, roster = evidence(population(), population("long"))
    protocol = protocol_for(roster, "nominal", "long")
    context = context_for(protocol)
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    result = m.evaluate_length_acceptance(rows, roster, context, protocol)
    assert result.status == "passed"
    assert result.pooled.metrics.eligible_count == 120
    assert result.pooled.paired_error_difference_interval[1] < 0
    assert tuple(result.strata) == ("long", "nominal")
    assert all(item.status == "passed" for item in result.strata.values())
    assert result.eligible_roster_sha256 == roster.sha256
    assert result.prediction_context_sha256 == context.sha256
    assert not hasattr(result, "promotion_eligible")
    with pytest.raises(TypeError):
        result.strata["other"] = result.pooled


def test_one_pooled_group_can_belong_to_multiple_mandatory_strata():
    rows, roster = evidence(population(memberships=["long", "nominal"]))
    result = evaluate(rows, roster, protocol_for(roster, "long", "nominal"))
    assert result.status == "passed"
    assert result.pooled.metrics.eligible_count == 60
    assert result.strata["long"].metrics.eligible_count == 60
    assert result.strata["nominal"].metrics.eligible_count == 60


def test_same_protocol_rejects_dropped_extra_or_swapped_outcomes():
    rows, roster = evidence(population())
    protocol = protocol_for(roster, "nominal")
    with pytest.raises(ValueError, match="outcome set"):
        evaluate(rows[:-1], roster, protocol)
    extra = LengthObservation("extra", "extra-group", 100, 100, 90)
    with pytest.raises(ValueError, match="outcome set"):
        evaluate((*rows, extra), roster, protocol)
    swapped = (
        replace(rows[0], group_key=rows[1].group_key),
        replace(rows[1], group_key=rows[0].group_key),
        *rows[2:],
    )
    with pytest.raises(ValueError, match="outcome set"):
        evaluate(swapped, roster, protocol)


def test_sparse_mandatory_stratum_cannot_be_pooled_away():
    rows, roster = evidence(population(count=120), population("long", count=59))
    result = evaluate(rows, roster, protocol_for(roster, "nominal", "long"))
    assert result.status == "insufficient-evidence"
    assert result.pooled.status == "passed"
    assert result.strata["long"].status == "insufficient-evidence"
    assert "insufficient_independent_groups" in result.strata["long"].reasons


def test_missing_mandatory_stratum_is_explicit():
    rows, roster = evidence(population())
    result = evaluate(rows, roster, protocol_for(roster, "nominal", "long"))
    assert result.status == "insufficient-evidence"
    assert result.strata["long"].metrics is None
    assert result.strata["long"].reasons == ("missing_required_stratum",)


def test_model_missingness_fails_coverage_without_shrinking_the_eligible_count():
    rows, roster = evidence(population())
    rows = (replace(rows[0], prediction=None),) + rows[1:]
    result = evaluate(rows, roster, protocol_for(roster, "nominal"))
    assert result.status == "failed"
    assert result.strata["nominal"].metrics.eligible_count == 60
    assert result.strata["nominal"].metrics.assessable_count == 59
    assert "availability_lower_bound" in result.strata["nominal"].reasons
    assert "insufficient_independent_groups" not in result.strata["nominal"].reasons


def test_long_array_bias_fails_absolute_error_despite_relative_tolerance():
    rows, roster = evidence(population("long"))
    rows = tuple(replace(row, truth=row.truth + 400, prediction=row.truth + 420) for row in rows)
    result = evaluate(rows, roster, protocol_for(roster, "long"))
    assert result.status == "failed"
    assert result.strata["long"].metrics.within_tolerance == 1
    assert result.strata["long"].reasons == ("maximum_mae",)


def test_constant_or_zero_baseline_cannot_pass_improvement_claim():
    rows, roster = evidence(population())
    constant = tuple(
        replace(row, truth=100 + i % 2, prediction=100.5, baseline_prediction=100.5) for i, row in enumerate(rows)
    )
    result = evaluate(constant, roster, protocol_for(roster, "nominal"))
    assert result.status == "failed"
    assert "relative_mae_improvement" in result.pooled.reasons
    assert "paired_error_improvement" in result.pooled.reasons
    zero = tuple(replace(row, baseline_prediction=row.truth) for row in rows)
    assert (
        "undefined_relative_mae_improvement" in evaluate(zero, roster, protocol_for(roster, "nominal")).pooled.reasons
    )


def test_no_assessable_predictions_have_explicit_failed_gates():
    rows, roster = evidence(population())
    missing = tuple(replace(row, prediction=None) for row in rows)
    result = evaluate(missing, roster, protocol_for(roster, "nominal"))
    assert result.status == "failed"
    assert "no_assessable_predictions" in result.pooled.reasons
    assert "undefined_paired_error_interval" in result.pooled.reasons


@pytest.mark.parametrize("error,passed", [(9.999, True), (10, True), (10.001, False)])
def test_absolute_mae_limit_is_inclusive_without_loosening_it(error, passed):
    rows, roster = evidence(population(error=error))
    result = evaluate(rows, roster, protocol_for(roster, "nominal"))
    assert ("maximum_mae" not in result.strata["nominal"].reasons) is passed


def test_roster_membership_and_all_typed_bindings_are_integrity_checked():
    rows, roster = evidence(population())
    protocol = protocol_for(roster, "nominal")
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    undeclared_rows, undeclared = evidence(population(memberships=["undeclared"]))
    undeclared_protocol = protocol_for(undeclared, "nominal")
    with pytest.raises(ValueError, match="strat"):
        m.evaluate_length_acceptance(undeclared_rows, undeclared, context_for(undeclared_protocol), undeclared_protocol)
    with pytest.raises(ValueError, match="digest"):
        m.evaluate_length_acceptance(rows, roster, context_for(protocol), replace(protocol, sha256="0" * 64))
    other_rows, other_roster = evidence(population("other"))
    with pytest.raises(ValueError, match="roster"):
        m.evaluate_length_acceptance(other_rows, other_roster, context_for(protocol), protocol)
    context = context_for(protocol)
    with pytest.raises(ValueError, match="context"):
        m.evaluate_length_acceptance(rows, roster, replace(context, qc_sha256="0" * 64), protocol)


def test_prediction_context_is_closed_hash_bound_and_only_an_integrity_binding():
    rows, roster = evidence(population())
    protocol = protocol_for(roster, "nominal")
    context = context_for(protocol)
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    document = m.length_prediction_context_document(context)
    assert document["schema_version"] == "length-prediction-context-v1"
    assert "qc_applied" not in document
    assert "promotion_eligible" not in document
    with pytest.raises(ValueError, match="fields"):
        m.decode_length_prediction_context({**document, "qc_applied": True})
    with pytest.raises(ValueError, match="schema_version"):
        m.decode_length_prediction_context({**document, "schema_version": "length-prediction-context-v0"})
    with pytest.raises(ValueError, match="qc_sha256"):
        m.decode_length_prediction_context({**document, "qc_sha256": "bad"})
    with pytest.raises(ValueError, match="canonical"):
        m.length_prediction_context_document(replace(context, sha256="0" * 64))
    assert evaluate(rows, roster, protocol).status == "passed"


def test_study_specific_tolerance_and_stricter_paired_limit_are_consumed():
    rows, roster = evidence(population(error=2))

    def changes(raw):
        raw["acceptance"].update(
            tolerance_absolute=1,
            tolerance_relative=0.001,
            paired_error_difference_upper_limit=-100,
        )

    result = evaluate(rows, roster, protocol_for(roster, "nominal", changes=changes))
    assert result.status == "failed"
    assert "tolerance_lower_bound" in result.pooled.reasons
    assert "tolerance_lower_bound" in result.strata["nominal"].reasons
    assert "paired_error_improvement" in result.pooled.reasons
