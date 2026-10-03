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
in seven a failure. Downsampling 42 confirmed positive GRCh38 exomes (378 runs) detected
the variant in 160 of 161 runs at 100x or more, 95 of 108 between 50x and 100x, and 48
of 109 below 50x. So a mean between the two lines is ``REDUCED`` - recorded, shown, and
still passing - and only a mean below the lower one fails. The lines are set per
assembly, because the window mean is not comparable across them; see
:func:`resolve_mean_thresholds`.

Functions:
    resolve_mean_thresholds: The mean-depth lines that apply to one run
    evaluate_coverage_qc: Two metrics and their thresholds to a verdict
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)

#: The verdict written into the coverage summary's ``coverage_qc`` column.
COVERAGE_QC_PASS = "PASS"
COVERAGE_QC_FAIL = "FAIL"

#: The mean is below the adequate line and at or above the low one. It passes the gate:
#: a weak variant signal can be missed at this depth, which qualifies a result without a
#: finding and is no reason to distrust a call.
COVERAGE_QC_REDUCED = "REDUCED"

#: The verdict when there is nothing to judge. Distinct from ``PASS`` on purpose: a report
#: that prints "Coverage QC: PASS" beside "Mean Coverage: Not calculated" is asserting a
#: quality claim it never checked. ``passed`` stays True for this state, so the screening
#: axis keeps the behaviour pinned by
#: ``test_screening_summary.py::test_coverage_that_was_never_measured_passes_the_quality_gate``
#: - only the displayed status becomes honest.
COVERAGE_QC_NOT_EVALUATED = "NOT_EVALUATED"

#: Shipped defaults for a configuration that omits the keys. Keep in step with
#: ``config.json``.
DEFAULT_MEAN_THRESHOLD = 100
DEFAULT_MEAN_THRESHOLDS_BY_ASSEMBLY: Mapping[str, Mapping[str, float]] = {
    "GRCh38": {"adequate": 100, "low": 50},
    "GRCh37": {"adequate": 290, "low": 145},
}

#: Reason identifiers. ``REASON_MEAN`` is a mean below the adequate line (``REDUCED``, or
#: ``FAIL`` where no band applies); ``REASON_MEAN_LOW`` is a mean below the low line.
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


@dataclass(frozen=True)
class MeanThresholds:
    """The mean-depth lines that apply to one run.

    Attributes:
        adequate: A mean at or above this is adequate.
        low: A mean below ``adequate`` and at or above this is ``REDUCED``; below it fails.
            ``None`` means a single line: anything below ``adequate`` fails.
    """

    adequate: float
    low: float | None


def _line(value: object, key: str) -> float:
    """One configured depth line, or a ``ValueError`` naming its key.

    Validated rather than coerced: a ``NaN`` line compares false against every mean, so
    it would silently pass every sample.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        message = f"thresholds.{key} must be a finite, non-negative number; got {value!r}"
        logger.error(message)
        raise ValueError(message)
    return float(value)


def resolve_mean_thresholds(thresholds: Mapping[str, Any], reference_assembly: str | None) -> MeanThresholds:
    """Return the mean-depth lines for this run's assembly.

    The window mean is not comparable across assemblies: the GRCh37 window holds about
    13.5 repeat units against about 58 on GRCh38, so the same reads give a GRCh37 mean
    2.89 times higher (median over 214 exomes realigned to both builds; 5th-95th
    percentile 2.82-2.96). Kestrel called the same variants on both builds in 326 of 326
    paired runs, so detection depends on the reads, and the GRCh37 lines are the GRCh38
    lines scaled by that ratio.

    Three cases, in order:

    1. ``thresholds.mean_vntr_coverage_by_assembly`` lists the run's coordinate system:
       its ``adequate`` and optional ``low`` apply.
    2. The block is absent but ``thresholds.mean_vntr_coverage`` is set: an operator's own
       configuration written before 2.0.42 (``--config-path`` replaces the whole config).
       Its single line applies on every assembly, as it did through 2.0.41.
    3. Otherwise the shipped defaults apply. A run whose assembly is unknown, or not
       listed, gets the single line at ``mean_vntr_coverage``.

    Args:
        thresholds (Mapping[str, Any]): ``config.json``'s ``thresholds`` block.
        reference_assembly (str | None): The run's declared assembly, in any accepted
            spelling, or ``None`` when it is not known.

    Returns:
        MeanThresholds: The lines to apply. ``low`` is ``None`` when the band does not
        apply or is empty because it is not below ``adequate``.

    Raises:
        ValueError: If a line that applies to this run is not a finite, non-negative
            number, or the per-assembly block or its entry is not a mapping.
    """
    single = MeanThresholds(
        adequate=_line(thresholds.get("mean_vntr_coverage", DEFAULT_MEAN_THRESHOLD), "mean_vntr_coverage"),
        low=None,
    )
    if "mean_vntr_coverage_by_assembly" in thresholds:
        by_assembly = thresholds["mean_vntr_coverage_by_assembly"]
    elif "mean_vntr_coverage" in thresholds:
        return single
    else:
        by_assembly = DEFAULT_MEAN_THRESHOLDS_BY_ASSEMBLY
    if not reference_assembly:
        return single
    from vntyper.scripts.reference_registry import get_coordinate_system

    try:
        coordinate_system = get_coordinate_system(reference_assembly)
    except ValueError:
        return single
    if not isinstance(by_assembly, Mapping):
        message = "thresholds.mean_vntr_coverage_by_assembly must be a mapping of assembly to lines"
        logger.error(message)
        raise ValueError(message)
    lines = by_assembly.get(coordinate_system)
    if lines is None:
        return single
    key = f"mean_vntr_coverage_by_assembly.{coordinate_system}"
    if not isinstance(lines, Mapping) or "adequate" not in lines:
        message = f"thresholds.{key} must be a mapping with an 'adequate' line"
        logger.error(message)
        raise ValueError(message)
    adequate = _line(lines["adequate"], f"{key}.adequate")
    low = lines.get("low")
    if low is not None:
        low = _line(low, f"{key}.low")
    if low is not None and low >= adequate:
        logger.info(
            f"The low mean line for {coordinate_system} ({low}) is not below its adequate "
            f"line ({adequate}); no reduced band applies."
        )
        low = None
    return MeanThresholds(adequate=adequate, low=None if low is None else float(low))


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
        mean_threshold (float): The adequate line from :func:`resolve_mean_thresholds`.
        percent_threshold (float): ``config.json``'s ``thresholds.percent_vntr_uncovered``.
        low_mean_threshold (float | None): The low line from
            :func:`resolve_mean_thresholds`. A mean at or above it and below
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
