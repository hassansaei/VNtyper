"""Reference length acceptance cannot hide bias, missingness or sparse strata."""

from dataclasses import replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_length_protocol import protocol_document
from vntyper.scripts.calibration_length_metrics import LengthObservation
from vntyper.scripts.calibration_length_protocol import decode_length_protocol

pytestmark = pytest.mark.unit


def rows_for(stratum="nominal", count=60, error=0):
    return tuple(
        LengthObservation(f"{stratum}-sample-{i}", f"{stratum}-group-{i}", stratum, 100 + i, 100 + i + error, 90)
        for i in range(count)
    )


def protocol_for(*strata):
    raw = protocol_document()
    raw["required_strata"] = sorted(strata)
    return decode_length_protocol(raw)


def test_perfect_fixed_predictions_pass_with_enough_independent_groups_in_every_stratum():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    result = m.evaluate_length_acceptance(rows_for() + rows_for("long"), protocol_for("nominal", "long"))
    assert result.status == "passed"
    assert result.pooled.metrics.eligible_count == 120
    assert result.pooled.paired_error_difference_interval[1] < 0
    assert tuple(result.strata) == ("long", "nominal")
    assert all(item.status == "passed" for item in result.strata.values())
    with pytest.raises(TypeError):
        result.strata["other"] = result.pooled


def test_sparse_mandatory_stratum_cannot_be_pooled_away():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    result = m.evaluate_length_acceptance(
        rows_for(count=120) + rows_for("long", count=59), protocol_for("nominal", "long")
    )
    assert result.status == "insufficient-evidence"
    assert result.pooled.status == "passed"
    assert result.strata["long"].status == "insufficient-evidence"
    assert "insufficient_independent_groups" in result.strata["long"].reasons


def test_missing_mandatory_stratum_is_explicit():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    result = m.evaluate_length_acceptance(rows_for(), protocol_for("nominal", "long"))
    assert result.status == "insufficient-evidence"
    assert result.strata["long"].metrics is None
    assert result.strata["long"].reasons == ("missing_required_stratum",)


def test_model_missingness_fails_coverage_without_shrinking_the_eligible_count():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    rows = rows_for()
    rows = (replace(rows[0], prediction=None),) + rows[1:]
    result = m.evaluate_length_acceptance(rows, protocol_for("nominal"))
    assert result.status == "failed"
    assert result.strata["nominal"].metrics.eligible_count == 60
    assert result.strata["nominal"].metrics.assessable_count == 59
    assert "availability_lower_bound" in result.strata["nominal"].reasons
    assert "insufficient_independent_groups" not in result.strata["nominal"].reasons


def test_long_array_bias_fails_absolute_error_despite_perfect_correlation_and_relative_tolerance():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    rows = tuple(replace(row, truth=row.truth + 400, prediction=row.truth + 420) for row in rows_for("long"))
    result = m.evaluate_length_acceptance(rows, protocol_for("long"))
    assert result.status == "failed"
    assert result.strata["long"].metrics.within_tolerance == 1
    assert result.strata["long"].reasons == ("maximum_mae",)


def test_constant_baseline_predictions_cannot_pass_an_improvement_claim():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    rows = tuple(
        replace(row, truth=100 + i % 2, prediction=100.5, baseline_prediction=100.5) for i, row in enumerate(rows_for())
    )
    result = m.evaluate_length_acceptance(rows, protocol_for("nominal"))
    assert result.status == "failed"
    assert "relative_mae_improvement" in result.pooled.reasons
    assert "paired_error_improvement" in result.pooled.reasons
    assert result.pooled.metrics.mae == 0.5


def test_zero_baseline_error_and_no_assessable_predictions_have_explicit_failed_gates():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    zero = tuple(replace(row, baseline_prediction=row.truth) for row in rows_for())
    result = m.evaluate_length_acceptance(zero, protocol_for("nominal"))
    assert result.status == "failed"
    assert "undefined_relative_mae_improvement" in result.pooled.reasons
    missing = tuple(replace(row, prediction=None) for row in zero)
    result = m.evaluate_length_acceptance(missing, protocol_for("nominal"))
    assert result.status == "failed"
    assert "no_assessable_predictions" in result.pooled.reasons
    assert "undefined_paired_error_interval" in result.pooled.reasons


@pytest.mark.parametrize("error,passed", [(9.999, True), (10, True), (10.001, False)])
def test_absolute_mae_limit_is_inclusive_without_loosening_it(error, passed):
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    result = m.evaluate_length_acceptance(rows_for(error=error), protocol_for("nominal"))
    assert ("maximum_mae" not in result.strata["nominal"].reasons) is passed


def test_unexpected_strata_and_forged_protocol_are_integrity_errors():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    protocol = protocol_for("nominal")
    with pytest.raises(ValueError, match="strat"):
        m.evaluate_length_acceptance(rows_for("undeclared"), protocol)
    with pytest.raises(ValueError, match="digest"):
        m.evaluate_length_acceptance(rows_for(), replace(protocol, sha256="0" * 64))


def test_study_specific_tolerance_is_used_by_every_stratum():
    m = import_module("vntyper.scripts.calibration_length_acceptance")
    raw = protocol_document()
    raw["required_strata"] = ["nominal"]
    raw["acceptance"].update(tolerance_absolute=1, tolerance_relative=0.001)
    result = m.evaluate_length_acceptance(rows_for(error=2), decode_length_protocol(raw))
    assert result.status == "failed"
    assert "tolerance_lower_bound" in result.pooled.reasons
    assert "tolerance_lower_bound" in result.strata["nominal"].reasons
