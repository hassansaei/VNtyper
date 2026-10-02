"""
Unit tests for the coverage QC verdict (#172).

The verdict is the one thing that decides `quality_metrics_pass`, and it is also written
into the coverage TSV, so a disagreement between the two would put a FAIL column beside a
PASS sentence. These tests pin the rule, its boundaries and its rounding contract.
"""

import pytest

from vntyper.scripts.coverage_qc import (
    COVERAGE_QC_FAIL,
    COVERAGE_QC_NOT_EVALUATED,
    COVERAGE_QC_PASS,
    COVERAGE_QC_REDUCED,
    REASON_MEAN,
    REASON_MEAN_LOW,
    REASON_NOT_MEASURED,
    REASON_UNCOVERED,
    evaluate_coverage_qc,
    resolve_low_mean_threshold,
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


def test_an_unmeasured_sample_is_not_reported_as_passing():
    """A verdict of PASS for a sample with no coverage data is a claim never checked.

    The report would print "Coverage QC: PASS" beside "Mean Coverage: Not calculated".
    `passed` stays True so the screening axis keeps the behaviour pinned by
    `test_coverage_that_was_never_measured_passes_the_quality_gate`; only the displayed
    status becomes honest. Found by adversarial review of the PR.
    """
    qc = evaluate_coverage_qc(None, None, 100, 50.0)

    assert qc.status == COVERAGE_QC_NOT_EVALUATED
    assert qc.reasons == (REASON_NOT_MEASURED,)
    assert qc.passed is True, "the screening axis must not change; only the status does"


def test_one_measured_metric_is_still_evaluated():
    """NOT_EVALUATED is for *nothing* measured, not for a partially populated summary."""
    assert evaluate_coverage_qc(250.0, None, 100, 50.0).status == COVERAGE_QC_PASS
    assert evaluate_coverage_qc(10.0, None, 100, 50.0).status == COVERAGE_QC_FAIL
    assert evaluate_coverage_qc(None, 90.0, 100, 50.0).status == COVERAGE_QC_FAIL


# ---------------------------------------------------------------------------
# Three levels: adequate, reduced, low
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("mean", "status", "passed", "reasons"),
    [
        (100.0, COVERAGE_QC_PASS, True, ()),
        (99.99, COVERAGE_QC_REDUCED, True, (REASON_MEAN,)),
        (50.0, COVERAGE_QC_REDUCED, True, (REASON_MEAN,)),
        (49.99, COVERAGE_QC_FAIL, False, (REASON_MEAN_LOW,)),
        (0.0, COVERAGE_QC_FAIL, False, (REASON_MEAN_LOW,)),
    ],
)
def test_a_mean_between_the_two_thresholds_is_reduced_and_still_passes(mean, status, passed, reasons):
    """A real exome at 96x with 168x flank depth was reported as a failure.

    Downsampling 28 confirmed positives still detected 67 of 69 variants between 50x and
    100x, so that band qualifies a negative without failing the sample.
    """
    qc = evaluate_coverage_qc(mean, 5.0, 100, 50.0, low_mean_threshold=50)

    assert (qc.status, qc.passed, qc.reasons) == (status, passed, reasons)


def test_a_patchy_vntr_fails_at_a_reduced_mean_and_names_only_what_failed():
    """Catch REDUCED masking the uncovered-fraction failure, or the mean being blamed for it."""
    qc = evaluate_coverage_qc(75.0, 80.0, 100, 50.0, low_mean_threshold=50)

    assert qc.status == COVERAGE_QC_FAIL
    assert qc.reasons == (REASON_UNCOVERED,)


def test_both_low_failures_are_reported_in_declaration_order():
    qc = evaluate_coverage_qc(10.0, 90.0, 100, 50.0, low_mean_threshold=50)

    assert qc.reasons == (REASON_MEAN_LOW, REASON_UNCOVERED)


def test_a_reduced_mean_with_no_uncovered_figure_is_still_reduced():
    assert evaluate_coverage_qc(75.0, None, 100, 50.0, low_mean_threshold=50).status == COVERAGE_QC_REDUCED


@pytest.mark.parametrize("assembly", ["hg38", "GRCh38", "hg38_ensembl", "hg38_ncbi"])
def test_the_reduced_band_applies_on_every_spelling_of_the_assembly_it_was_measured_on(assembly):
    assert resolve_low_mean_threshold({}, assembly) == 50.0
    assert resolve_low_mean_threshold({"mean_vntr_coverage_low": 60}, assembly) == 60.0


@pytest.mark.parametrize("assembly", ["hg19", "GRCh37", "hg19_ensembl", None, "", "not-an-assembly"])
def test_other_and_unknown_assemblies_keep_the_single_line(assembly):
    """The GRCh37 window reads about 2.7 times higher for the same sample, so 50x there
    is far less depth than the 50x that was measured. A GRCh37 negative at 75x failed
    before this change and must still fail."""
    assert resolve_low_mean_threshold({}, assembly) is None
    low = resolve_low_mean_threshold({}, assembly)
    assert evaluate_coverage_qc(75.0, 0.0, 100, 50.0, low_mean_threshold=low).status == COVERAGE_QC_FAIL


def test_the_measured_assemblies_are_configurable():
    thresholds = {"reduced_depth_assemblies": ["GRCh37", "GRCh38"]}

    assert resolve_low_mean_threshold(thresholds, "hg19") == 50.0
    assert resolve_low_mean_threshold({"reduced_depth_assemblies": []}, "hg38") is None


@pytest.mark.parametrize("thresholds", [{"mean_vntr_coverage": 30}, {"mean_vntr_coverage_low": 100}])
def test_a_low_threshold_not_below_the_adequate_one_leaves_no_reduced_band(thresholds):
    """A replaced config can set the adequate threshold below the shipped low default."""
    assert resolve_low_mean_threshold(thresholds, "hg38") is None


def test_the_shipped_config_declares_the_band_for_grch38_only():
    import json
    from pathlib import Path

    import vntyper

    thresholds = json.loads((Path(vntyper.__file__).parent / "config.json").read_text())["thresholds"]

    assert thresholds["mean_vntr_coverage_low"] == 50
    assert thresholds["reduced_depth_assemblies"] == ["GRCh38"]
    assert resolve_low_mean_threshold(thresholds, "hg38") == 50.0
    assert resolve_low_mean_threshold(thresholds, "hg19") is None


def test_without_a_low_threshold_the_single_line_is_unchanged():
    """A direct caller that passes two thresholds gets the pre-existing verdict and reason."""
    qc = evaluate_coverage_qc(75.0, 5.0, 100, 50.0)

    assert (qc.status, qc.reasons) == (COVERAGE_QC_FAIL, (REASON_MEAN,))


# ---------------------------------------------------------------------------
# Judging a summary written before the region-wide coverage change (#171)
# ---------------------------------------------------------------------------


def test_a_pre_2_0_8_mean_is_corrected_before_it_is_judged():
    """Adversarial review of the PR: an old mean is not comparable with the thresholds.

    Before 2.0.8 the mean excluded uncovered bases, so a stored 150.0 at 40% uncovered
    stands for a region-wide 90.0 and should fail a 100x threshold - judging the stored
    figure passes it. `percent_uncovered` was already correct, so the identity
    `mean_old * (1 - pct/100)` recovers the region-wide value exactly; the golden-cohort
    gate confirmed it on 61 of 61 cases.
    """
    stored_mean, pct = 150.0, 40.0
    corrected = round(stored_mean * (1 - pct / 100), 2)

    assert corrected == 90.0
    assert evaluate_coverage_qc(stored_mean, pct, 100, 50.0).status == COVERAGE_QC_PASS, "the stored figure is lenient"
    assert evaluate_coverage_qc(corrected, pct, 100, 50.0).status == COVERAGE_QC_FAIL, "the corrected figure is right"
