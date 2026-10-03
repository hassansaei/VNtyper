"""Pure presentation decisions for coverage quality-control verdicts.

This module keeps small status-label decisions out of the oversized report renderer
and screening-state module. Interpretive screening wording remains configuration-owned.

Functions:
    coverage_qc_word: Translate durable status tokens for the report chip.
    coverage_qc_tone: Select the chip tone for a durable status token.
    coverage_uncovered_exceeded: Whether the verdict failed on the uncovered fraction.
    coverage_notes: Select and fill the configured sentences for a measured verdict.
    coverage_not_measured_note: Select the opt-in note for an unevaluated gate.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any, Final

from vntyper.scripts.coverage_qc import (
    COVERAGE_QC_FAIL,
    COVERAGE_QC_NOT_EVALUATED,
    COVERAGE_QC_PASS,
    COVERAGE_QC_REDUCED,
    REASON_MEAN,
    REASON_MEAN_LOW,
    REASON_UNCOVERED,
    CoverageQC,
)

logger = logging.getLogger(__name__)

#: The config key naming the sentence rendered when the coverage gate had nothing to
#: judge. Config-driven and opt-in: older report configurations render as before.
COVERAGE_NOT_MEASURED_NOTE_KEY: Final[str] = "coverage_not_measured_note"

#: The coverage QC verdict as the chip prints it. Status tokens are durable record
#: vocabulary, while the chip row speaks in words like every other chip.
_COVERAGE_QC_WORDS: Final[dict[str, str]] = {
    COVERAGE_QC_PASS: "Adequate",
    COVERAGE_QC_REDUCED: "Reduced",
    COVERAGE_QC_FAIL: "Insufficient",
    COVERAGE_QC_NOT_EVALUATED: "Not evaluated",
}

#: The config key holding the sentences rendered for a measured coverage verdict, and the
#: four sentences it must declare. Opt-in as a block: a configuration without it renders
#: no depth sentence, and one that declares it incompletely is malformed.
COVERAGE_NOTES_KEY: Final[str] = "coverage_notes"
_COVERAGE_NOTE_FIELDS: Final[frozenset[str]] = frozenset(
    {"reduced_no_finding", "low_no_finding", "uncovered", "action_no_finding"}
)

#: The chip tone per status. The strings are the shared token layer's tone names.
_COVERAGE_QC_TONES: Final[dict[str, str]] = {
    COVERAGE_QC_PASS: "ok",
    COVERAGE_QC_FAIL: "caution",
}


def coverage_qc_word(status: str) -> str:
    """Return one coverage QC status as the chip row prints it.

    Args:
        status: ``CoverageQC.status`` or a future durable status token.

    Returns:
        str: The display word for a known token; unknown tokens pass through unchanged.
    """
    return _COVERAGE_QC_WORDS.get(status, status)


def coverage_qc_tone(status: str) -> str:
    """Return the chip tone for one coverage QC status.

    Only a failing depth is a caution. A reduced depth qualifies a result without a
    finding in words and carries no colour, like a run with nothing to judge.

    Args:
        status: ``CoverageQC.status`` or a future durable status token.

    Returns:
        str: ``"ok"``, ``"caution"`` or ``"none"``; an unknown token is ``"none"``.
    """
    return _COVERAGE_QC_TONES.get(status, "none")


def coverage_not_measured_note(report_config: dict[str, Any], coverage_qc: CoverageQC) -> str:
    """Return the configured note when the coverage gate was not evaluated.

    ``coverage_qc.passed`` deliberately stays true for ``NOT_EVALUATED`` so the
    screening axis is unchanged. This presentation decision therefore uses the
    explicit status and never infers measuredness from pass/fail.

    Args:
        report_config: The parsed ``report_config.json``.
        coverage_qc: The coverage QC verdict.

    Returns:
        str: The configured note when nothing was measured, otherwise ``""``.
    """
    if coverage_qc.status != COVERAGE_QC_NOT_EVALUATED:
        return ""
    note = str(report_config.get(COVERAGE_NOT_MEASURED_NOTE_KEY, "") or "")
    if note:
        logger.info("Coverage quality gate was not evaluated; rendering the configured note.")
    return note


def coverage_uncovered_exceeded(coverage_qc: CoverageQC) -> bool:
    """Return whether the verdict failed on the uncovered fraction.

    A region with more than the configured share at zero depth is a gross failure of the
    alignment over the VNTR, not a matter of depth, and it qualifies a call as well as a
    result without one. ``status`` cannot say this: ``FAIL`` also means a low mean.

    Args:
        coverage_qc: The coverage QC verdict.

    Returns:
        bool: True when the uncovered-fraction threshold is among the failed reasons.
    """
    return REASON_UNCOVERED in coverage_qc.reasons


def _shown(value: float) -> str:
    """One figure as the note prints it: no trailing zeros, two decimals at most."""
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _threshold_shown(value: float) -> str:
    """A configured line as the note prints it, never rounded.

    Rounding a line at any precision could print "50x (below 50x)" for a mean of 50.00
    judged against a line just above 50, so the line is printed as the shortest string
    that round-trips to it, without a trailing ``.0``.
    """
    shown = repr(float(value))
    return shown.removesuffix(".0")


def coverage_notes(
    report_config: Mapping[str, Any],
    coverage_qc: CoverageQC,
    *,
    is_positive: bool,
    mean_vntr_coverage: float | None,
    percent_vntr_uncovered: float | None,
    mean_threshold: float,
    low_mean_threshold: float | None,
    percent_threshold: float,
) -> tuple[str, ...]:
    """Return the configured sentences that qualify a result by its measured coverage.

    Depth limits what a result without a finding means and says nothing against a call:
    no false positive appeared in 1,200 simulated negative runs across six depth levels,
    and the precision tier already reflects the depth supporting a call. So the two mean
    sentences and the action are selected only when the sample has no finding. A region
    that is mostly uncovered is a different fact and is stated for every result.

    The sentences are configuration; this fills in the measured figure and the threshold
    that was applied, so the wording cannot drift from the numbers and never states a
    threshold the run did not use.

    Args:
        report_config: The parsed ``report_config.json``.
        coverage_qc: The coverage QC verdict.
        is_positive: Whether either algorithm reported a finding.
        mean_vntr_coverage: The mean the verdict was evaluated on, or ``None``.
        percent_vntr_uncovered: The uncovered percentage it was evaluated on, or ``None``.
        mean_threshold: The adequate mean threshold applied.
        low_mean_threshold: The low mean threshold applied, or ``None`` for a single line.
        percent_threshold: The uncovered-percentage threshold applied.

    Returns:
        tuple[str, ...]: The sentences in rendering order; empty for a passing or
        unmeasured verdict, and for a configuration that declares no ``coverage_notes``.

    Raises:
        ValueError: If ``coverage_notes`` is declared with fields other than the four
            required sentences, or one of them is not non-empty text.
    """
    raw = report_config.get(COVERAGE_NOTES_KEY)
    if raw is None:
        return ()
    if not isinstance(raw, Mapping) or set(raw) != _COVERAGE_NOTE_FIELDS:
        message = "coverage_notes wording differs from the closed contract"
        logger.error(message)
        raise ValueError(message)
    for name in sorted(_COVERAGE_NOTE_FIELDS):
        if not isinstance(raw[name], str) or not raw[name].strip():
            message = f"coverage_notes {name} must be non-empty text"
            logger.error(message)
            raise ValueError(message)

    notes: list[str] = []
    mean = "" if mean_vntr_coverage is None else _shown(mean_vntr_coverage)
    if coverage_uncovered_exceeded(coverage_qc) and percent_vntr_uncovered is not None:
        notes.append(
            raw["uncovered"].format(percent=_shown(percent_vntr_uncovered), limit=_threshold_shown(percent_threshold))
        )
    if not is_positive and mean_vntr_coverage is not None:
        if coverage_qc.status == COVERAGE_QC_REDUCED:
            notes.append(raw["reduced_no_finding"].format(mean=mean, threshold=_threshold_shown(mean_threshold)))
        elif REASON_MEAN_LOW in coverage_qc.reasons and low_mean_threshold is not None:
            notes.append(raw["low_no_finding"].format(mean=mean, threshold=_threshold_shown(low_mean_threshold)))
        elif coverage_qc.status == COVERAGE_QC_FAIL and REASON_MEAN in coverage_qc.reasons:
            notes.append(raw["low_no_finding"].format(mean=mean, threshold=_threshold_shown(mean_threshold)))
    if not is_positive and coverage_qc.status == COVERAGE_QC_FAIL:
        notes.append(raw["action_no_finding"])
    if notes:
        logger.info("Rendering %d configured coverage note(s) for status %s.", len(notes), coverage_qc.status)
    return tuple(notes)
