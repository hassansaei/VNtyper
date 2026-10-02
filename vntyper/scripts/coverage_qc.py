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

**Three levels of depth, not two.** A single 100x pass/fail line called one real exome
in seven a failure, while downsampling 28 confirmed positive hg38 exomes still detected
67 of 69 variants between 50x and 100x and 124 of 124 above it; below 50x it was 35 of
59. So on GRCh38 a mean between the two thresholds is ``REDUCED`` - recorded, shown, and
still passing - and only a mean below the lower one fails. Other assemblies keep the
single line: the band was not measured there.

Functions:
    resolve_low_mean_threshold: The lower threshold that applies to one run, if any
    evaluate_coverage_qc: Two metrics and their thresholds to a verdict
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

#: The verdict written into the coverage summary's ``coverage_qc`` column.
COVERAGE_QC_PASS = "PASS"
COVERAGE_QC_FAIL = "FAIL"

#: The mean is below the adequate threshold and at or above the low one. It passes the
#: gate: detection is slightly less sensitive at this depth, which qualifies a result
#: without a finding and is no reason to distrust a call.
COVERAGE_QC_REDUCED = "REDUCED"

#: The verdict when there is nothing to judge. Distinct from ``PASS`` on purpose: a report
#: that prints "Coverage QC: PASS" beside "Mean Coverage: Not calculated" is asserting a
#: quality claim it never checked. ``passed`` stays True for this state, so the screening
#: axis keeps the behaviour pinned by
#: ``test_screening_summary.py::test_coverage_that_was_never_measured_passes_the_quality_gate``
#: - only the displayed status becomes honest.
COVERAGE_QC_NOT_EVALUATED = "NOT_EVALUATED"

#: Reason identifiers, named after the ``config.json`` threshold keys they come from so a
#: consumer can look up the number that was applied.
#: Shipped defaults for a configuration that omits the keys.
DEFAULT_LOW_MEAN_THRESHOLD = 50
DEFAULT_REDUCED_DEPTH_ASSEMBLIES: tuple[str, ...] = ("GRCh38",)

REASON_MEAN = "mean_vntr_coverage"
REASON_MEAN_LOW = "mean_vntr_coverage_low"
REASON_UNCOVERED = "percent_vntr_uncovered"

#: Why a verdict could not be reached: no ``Coverage Calculation`` step ran at all.
REASON_NOT_MEASURED = "coverage_not_measured"


@dataclass(frozen=True)
class CoverageQC:
    """The coverage quality verdict for one sample.

    Attributes:
        passed: Whether the sample met every configured coverage threshold.
        status: :data:`COVERAGE_QC_PASS`, :data:`COVERAGE_QC_REDUCED`,
            :data:`COVERAGE_QC_FAIL` or :data:`COVERAGE_QC_NOT_EVALUATED`. This is the
            value written to the ``coverage_qc`` column.
        reasons: The threshold keys the sample fell short of, in declaration order. Empty
            for ``PASS``; a ``REDUCED`` verdict names the mean and still passes.
    """

    passed: bool
    status: str
    reasons: tuple[str, ...]


def resolve_low_mean_threshold(thresholds: Mapping[str, Any], reference_assembly: str | None) -> float | None:
    """Return the low mean threshold that applies to this run, or ``None`` for a single line.

    The reduced band was measured on the GRCh38 window only. The GRCh37 window holds less
    of the repeat array, so its mean reads about 2.7 times higher for the same sample, and
    50x there is far less depth than 50x on GRCh38. The band therefore applies only to the
    coordinate systems named in ``thresholds.reduced_depth_assemblies``; every other run,
    and one whose assembly is not recorded, keeps the single line at
    ``thresholds.mean_vntr_coverage``.

    Args:
        thresholds (Mapping[str, Any]): ``config.json``'s ``thresholds`` block. ``.get`` with
            the shipped defaults throughout: ``--config-path`` replaces the whole config.
        reference_assembly (str | None): The run's declared assembly, in any accepted
            spelling, or ``None`` when it is not known.

    Returns:
        float | None: The low threshold, or ``None`` when the band does not apply or is
        empty because the low threshold is not below the adequate one.
    """
    if not reference_assembly:
        return None
    from vntyper.scripts.reference_registry import get_coordinate_system

    try:
        coordinate_system = get_coordinate_system(reference_assembly)
    except ValueError:
        return None
    if coordinate_system not in thresholds.get("reduced_depth_assemblies", DEFAULT_REDUCED_DEPTH_ASSEMBLIES):
        return None
    mean_threshold = thresholds.get("mean_vntr_coverage", 100)
    low_mean_threshold = thresholds.get("mean_vntr_coverage_low", DEFAULT_LOW_MEAN_THRESHOLD)
    if low_mean_threshold >= mean_threshold:
        logger.info(
            f"thresholds.mean_vntr_coverage_low ({low_mean_threshold}) is not below "
            f"thresholds.mean_vntr_coverage ({mean_threshold}); no reduced band applies."
        )
        return None
    return float(low_mean_threshold)


def evaluate_coverage_qc(
    mean_vntr_coverage: float | None,
    percent_vntr_uncovered: float | None,
    mean_threshold: float,
    percent_threshold: float,
    *,
    low_mean_threshold: float | None = None,
) -> CoverageQC:
    """Decide whether a sample's VNTR coverage meets the configured thresholds.

    Args:
        mean_vntr_coverage (float | None): Mean depth over the VNTR region, as published
            (two decimal places). ``None`` when no coverage step ran.
        percent_vntr_uncovered (float | None): Percentage of the region at zero depth, as
            published. ``None`` when no coverage step ran.
        mean_threshold (float): ``config.json``'s ``thresholds.mean_vntr_coverage``.
        percent_threshold (float): ``config.json``'s ``thresholds.percent_vntr_uncovered``.
        low_mean_threshold (float | None): ``config.json``'s
            ``thresholds.mean_vntr_coverage_low``. A mean at or above it and below
            ``mean_threshold`` is ``REDUCED`` and passes. ``None`` keeps the single
            pass/fail line at ``mean_threshold``.

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

    fail_below = mean_threshold if low_mean_threshold is None else low_mean_threshold
    reasons: list[str] = []
    reduced = False

    if mean_vntr_coverage is not None and mean_vntr_coverage < fail_below:
        reasons.append(REASON_MEAN if low_mean_threshold is None else REASON_MEAN_LOW)
    elif mean_vntr_coverage is not None and mean_vntr_coverage < mean_threshold:
        reduced = True
    if percent_vntr_uncovered is not None and percent_vntr_uncovered > percent_threshold:
        reasons.append(REASON_UNCOVERED)

    if reasons:
        logger.info(f"Coverage QC failed on: {', '.join(reasons)}")
        return CoverageQC(passed=False, status=COVERAGE_QC_FAIL, reasons=tuple(reasons))
    if reduced:
        logger.info(f"Coverage QC reduced: mean {mean_vntr_coverage} is below {mean_threshold}.")
        return CoverageQC(passed=True, status=COVERAGE_QC_REDUCED, reasons=(REASON_MEAN,))
    return CoverageQC(passed=True, status=COVERAGE_QC_PASS, reasons=())
