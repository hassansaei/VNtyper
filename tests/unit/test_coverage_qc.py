"""
Unit tests for the coverage QC verdict (#172).

The verdict is the one thing that decides `quality_metrics_pass`, and it is also written
into the coverage TSV, so a disagreement between the two would put a FAIL column beside a
PASS sentence. These tests pin the rule, its boundaries and its rounding contract.
"""

import pytest

from vntyper.scripts.coverage_qc import (
    COVERAGE_QC_FAIL,
    COVERAGE_QC_PASS,
    REASON_MEAN,
    REASON_UNCOVERED,
    evaluate_coverage_qc,
)

pytestmark = pytest.mark.unit


def test_a_sample_above_both_thresholds_passes():
    qc = evaluate_coverage_qc(250.0, 5.0, 100, 50.0)

    assert qc.passed is True
    assert qc.status == COVERAGE_QC_PASS
    assert qc.reasons == ()


def test_a_low_mean_fails_and_names_the_mean():
    qc = evaluate_coverage_qc(99.0, 5.0, 100, 50.0)

    assert qc.passed is False
    assert qc.status == COVERAGE_QC_FAIL
    assert qc.reasons == (REASON_MEAN,)


def test_a_patchy_vntr_fails_even_with_acceptable_mean():
    """#172's headline case, and the reason the issue exists.

    Half the VNTR uncovered is precisely where a frameshift call can be missed, yet
    before this change the sample passed QC on its mean alone.
    """
    qc = evaluate_coverage_qc(250.0, 80.0, 100, 50.0)

    assert qc.passed is False
    assert qc.reasons == (REASON_UNCOVERED,)


def test_both_failures_are_reported_in_declaration_order():
    qc = evaluate_coverage_qc(10.0, 90.0, 100, 50.0)

    assert qc.reasons == (REASON_MEAN, REASON_UNCOVERED)


@pytest.mark.parametrize("mean", [None, 250.0])
@pytest.mark.parametrize("percent", [None, 5.0])
def test_a_metric_that_was_never_measured_does_not_fail_the_gate(mean, percent):
    """Preserved from the pre-#172 behaviour, deliberately.

    `test_screening_summary.py` pins that an unmeasured sample reports as passing. That
    is a displayed interpretation for every run with no Coverage Calculation step, so
    #172 does not change it - it only adds the second metric.
    """
    assert evaluate_coverage_qc(mean, percent, 100, 50.0).passed is True


def test_the_boundaries_pass_on_equality():
    """Mean fails strictly below; uncovered fails strictly above. Asymmetric on purpose:
    it matches `threshold_icon`'s `higher_better` argument in report_formatting."""
    assert evaluate_coverage_qc(100.0, 50.0, 100, 50.0).passed is True
    assert evaluate_coverage_qc(99.99, 50.0, 100, 50.0).passed is False
    assert evaluate_coverage_qc(100.0, 50.01, 100, 50.0).passed is False


def test_the_verdict_is_a_function_of_the_published_figures():
    """The rounding contract (#172, adversarial review A1).

    `format_coverage_summary` writes mean and percent with `:.2f`, and the report reads
    those strings back. If the writer evaluated the raw value and the reader the rounded
    one, a raw mean of 99.999 would emit FAIL into the TSV while the report recomputed
    PASS. Callers round before calling; this test pins that both sides then agree, and
    that the answer matches what the report prints.
    """
    raw, published = 99.999, round(99.999, 2)

    assert published == 100.0
    assert evaluate_coverage_qc(raw, 0.0, 100, 50.0).passed is False
    assert evaluate_coverage_qc(published, 0.0, 100, 50.0).passed is True


def test_the_verdict_is_frozen():
    """A verdict a consumer could mutate is not a verdict."""
    qc = evaluate_coverage_qc(250.0, 5.0, 100, 50.0)

    with pytest.raises(AttributeError):
        qc.passed = False  # type: ignore[misc]
