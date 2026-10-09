"""Pure report state for the Kestrel/adVNTR cross-match step."""

from __future__ import annotations

from typing import Any

from vntyper.scripts.nomenclature import FLAG_CALLER_DISAGREEMENT
from vntyper.scripts.summary_steps import (
    STEP_ABSENT,
    STEP_CROSS_MATCH,
    STEP_KESTREL,
    STEP_UNREADABLE,
    get_step_data,
    get_step_state,
)

MATCH_MESSAGE = "At least one match was found between Kestrel and adVNTR results."
NO_MATCH_MESSAGE = "No matches were found between Kestrel and adVNTR results."
RECONCILED_MESSAGE = (
    "Kestrel and adVNTR name one allele after reconciliation; their raw records are written differently."
)


def _reconciled_to_one_allele(kestrel_rows: list[dict[str, Any]]) -> bool:
    """Whether reconciliation named one allele from both callers on some Kestrel row.

    The cross-match compares each caller's raw inserted or deleted bases, written in
    its own frame, so it cannot see that two differently written records are one
    molecule. Reconciliation can; when it names adVNTR's allele on a Kestrel row
    without a disagreement, the two callers agree, whatever their records look like.
    """
    return any(
        row.get("Nomenclature")
        and row.get("Nomenclature_Kestrel")
        and row.get("Nomenclature") == row.get("Nomenclature_adVNTR")
        and row.get("Nomenclature_Tier") in ("A", "B")
        and FLAG_CALLER_DISAGREEMENT not in str(row.get("Nomenclature_Flags", "")).split(";")
        for row in kestrel_rows
    )


def build_cross_match_summary(
    pipeline_summary: dict[str, Any], report_config: dict[str, Any]
) -> tuple[str, bool, bool]:
    """Return the cross-match sentence, match state, and assessability.

    Args:
        pipeline_summary: Parsed ``pipeline_summary.json``.
        report_config: Parsed report configuration.

    Returns:
        The message, whether the callers agree, and whether a comparison result
        was readable. The callers agree when a raw row matched, or when
        reconciliation named one allele from both. An absent step has no message. An unreadable step retains
        the legacy no-match sentence when an older config provides no dedicated
        wording, but remains structurally not assessable.
    """
    state = get_step_state(pipeline_summary, STEP_CROSS_MATCH)
    if state == STEP_ABSENT:
        return "", False, False
    if state == STEP_UNREADABLE:
        cross_match_config = report_config.get("cross_match")
        if isinstance(cross_match_config, dict):
            configured_message = cross_match_config.get("not_assessable_message")
            if isinstance(configured_message, str) and configured_message:
                return configured_message, False, False
        return NO_MATCH_MESSAGE, False, False

    data = get_step_data(pipeline_summary, STEP_CROSS_MATCH)
    if any(item.get("Match") == "Yes" for item in data):
        return MATCH_MESSAGE, True, True
    if _reconciled_to_one_allele(get_step_data(pipeline_summary, STEP_KESTREL)):
        return RECONCILED_MESSAGE, True, True
    return NO_MATCH_MESSAGE, False, True
