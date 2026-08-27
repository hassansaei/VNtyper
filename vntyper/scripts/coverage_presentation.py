"""Pure presentation decisions for coverage quality-control verdicts.

This module keeps small status-label decisions out of the oversized report renderer
and screening-state module. Interpretive screening wording remains configuration-owned.

Functions:
    coverage_qc_word: Translate durable status tokens for the report chip.
    coverage_not_measured_note: Select the opt-in note for an unevaluated gate.
"""

from __future__ import annotations

import logging
from typing import Any, Final

from vntyper.scripts.coverage_qc import (
    COVERAGE_QC_FAIL,
    COVERAGE_QC_NOT_EVALUATED,
    COVERAGE_QC_PASS,
    CoverageQC,
)

logger = logging.getLogger(__name__)

#: The config key naming the sentence rendered when the coverage gate had nothing to
#: judge. Config-driven and opt-in: older report configurations render as before.
COVERAGE_NOT_MEASURED_NOTE_KEY: Final[str] = "coverage_not_measured_note"

#: The coverage QC verdict as the chip prints it. Status tokens are durable record
#: vocabulary, while the chip row speaks in words like every other chip.
_COVERAGE_QC_WORDS: Final[dict[str, str]] = {
    COVERAGE_QC_PASS: "Pass",
    COVERAGE_QC_FAIL: "Fail",
    COVERAGE_QC_NOT_EVALUATED: "Not evaluated",
}


def coverage_qc_word(status: str) -> str:
    """Return one coverage QC status as the chip row prints it.

    Args:
        status: ``CoverageQC.status`` or a future durable status token.

    Returns:
        str: The display word for a known token; unknown tokens pass through unchanged.
    """
    return _COVERAGE_QC_WORDS.get(status, status)


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
