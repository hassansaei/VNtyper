"""
coverage_qc.py

The sample-level coverage quality verdict (#172).

Before this module, ``quality_metrics_pass`` considered mean VNTR coverage only.
``percent_vntr_uncovered`` was configured with a threshold, computed on every run and
compared to nothing - it drove a red or green icon and no decision. A sample with
acceptable mean coverage and half the VNTR uncovered passed QC, which is the opposite of
the desirable failure mode: a patchy VNTR is exactly where a frameshift call can be
missed.

**The verdict is a function of the *published* figures, not the raw ones.**
``coverage_stats.format_coverage_summary`` writes ``mean`` and ``percent_uncovered`` with
two decimal places, and the report reads those strings back out of
``pipeline_summary.json``. Evaluating the raw value in one place and the rounded value in
the other lets the two disagree at a threshold boundary - a raw mean of 99.999 is below a
threshold of 100, but serialises as ``100.00``. Callers therefore round before calling,
and the report prints no ``FAIL`` beside a displayed ``100.00``.

Functions:
    evaluate_coverage_qc: Two metrics and two thresholds to a verdict
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)

#: The verdict written into the coverage summary's ``coverage_qc`` column.
COVERAGE_QC_PASS = "PASS"
COVERAGE_QC_FAIL = "FAIL"

#: The verdict when there is nothing to judge. Distinct from ``PASS`` on purpose: a report
#: that prints "Coverage QC: PASS" beside "Mean Coverage: Not calculated" is asserting a
#: quality claim it never checked. ``passed`` stays True for this state, so the screening
#: axis keeps the behaviour pinned by
#: ``test_screening_summary.py::test_coverage_that_was_never_measured_passes_the_quality_gate``
#: - only the displayed status becomes honest.
COVERAGE_QC_NOT_EVALUATED = "NOT_EVALUATED"

#: Reason identifiers, named after the ``config.json`` threshold keys they come from so a
#: consumer can look up the number that was applied.
REASON_MEAN = "mean_vntr_coverage"
REASON_UNCOVERED = "percent_vntr_uncovered"

#: Why a verdict could not be reached: no ``Coverage Calculation`` step ran at all.
REASON_NOT_MEASURED = "coverage_not_measured"


@dataclass(frozen=True)
class CoverageQC:
    """The coverage quality verdict for one sample.

    Attributes:
        passed: Whether the sample met every configured coverage threshold.
        status: :data:`COVERAGE_QC_PASS` or :data:`COVERAGE_QC_FAIL`. This is the value
            written to the ``coverage_qc`` column.
        reasons: The threshold keys that failed, in declaration order. Empty when passed.
    """

    passed: bool
    status: str
    reasons: tuple[str, ...]


def evaluate_coverage_qc(
    mean_vntr_coverage: float | None,
    percent_vntr_uncovered: float | None,
    mean_threshold: float,
    percent_threshold: float,
) -> CoverageQC:
    """Decide whether a sample's VNTR coverage meets both configured thresholds.

    Args:
        mean_vntr_coverage (float | None): Mean depth over the VNTR region, as published
            (two decimal places). ``None`` when no coverage step ran.
        percent_vntr_uncovered (float | None): Percentage of the region at zero depth, as
            published. ``None`` when no coverage step ran.
        mean_threshold (float): ``config.json``'s ``thresholds.mean_vntr_coverage``.
        percent_threshold (float): ``config.json``'s ``thresholds.percent_vntr_uncovered``.

    Returns:
        CoverageQC: The verdict. A metric that is ``None`` never fails the gate - an
        unmeasured sample reported as failing would change a displayed interpretation for
        every run with no coverage step, which is out of scope for #172.

    Note:
        The comparisons are asymmetric on purpose, matching ``threshold_icon``'s
        ``higher_better`` argument: the mean fails strictly *below* its threshold, the
        uncovered fraction fails strictly *above* its own. A sample at exactly 100x and
        exactly 50.0% uncovered passes both.
    """
    # Nothing was measured, so there is nothing to judge. Reporting PASS here would state a
    # quality claim that was never checked - the report would print "Coverage QC: PASS"
    # beside "Mean Coverage: Not calculated". `passed` stays True so the screening axis is
    # unchanged; only the status tells the truth.
    if mean_vntr_coverage is None and percent_vntr_uncovered is None:
        logger.info("Coverage QC not evaluated: no coverage was measured for this sample.")
        return CoverageQC(passed=True, status=COVERAGE_QC_NOT_EVALUATED, reasons=(REASON_NOT_MEASURED,))

    reasons: list[str] = []

    if mean_vntr_coverage is not None and mean_vntr_coverage < mean_threshold:
        reasons.append(REASON_MEAN)
    if percent_vntr_uncovered is not None and percent_vntr_uncovered > percent_threshold:
        reasons.append(REASON_UNCOVERED)

    passed = not reasons
    if not passed:
        logger.info(f"Coverage QC failed on: {', '.join(reasons)}")

    return CoverageQC(
        passed=passed,
        status=COVERAGE_QC_PASS if passed else COVERAGE_QC_FAIL,
        reasons=tuple(reasons),
    )
