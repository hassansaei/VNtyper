"""Pure report state for the Kestrel/adVNTR cross-match step."""

from __future__ import annotations

from typing import Any

from vntyper.scripts.summary_steps import (
    STEP_ABSENT,
    STEP_CROSS_MATCH,
    STEP_UNREADABLE,
    get_step_data,
    get_step_state,
)

MATCH_MESSAGE = "At least one match was found between Kestrel and adVNTR results."
NO_MATCH_MESSAGE = "No matches were found between Kestrel and adVNTR results."


def build_cross_match_summary(
    pipeline_summary: dict[str, Any], report_config: dict[str, Any]
) -> tuple[str, bool, bool]:
    """Return the cross-match sentence, match state, and assessability.

    Args:
        pipeline_summary: Parsed ``pipeline_summary.json``.
        report_config: Parsed report configuration.

    Returns:
        The message, whether any row matched, and whether a comparison result
        was readable. An absent step has no message. An unreadable step retains
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
    is_positive = any(item.get("Match") == "Yes" for item in data)
    return (MATCH_MESSAGE if is_positive else NO_MATCH_MESSAGE), is_positive, True
