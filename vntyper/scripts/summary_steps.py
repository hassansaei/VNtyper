"""
summary_steps.py

Module Purpose:
---------------
Single source of truth for the pipeline-summary step names and the accessors
used to read a step's parsed result.

These five names are matched by exact string comparison across pipeline.py,
generate_report.py, cohort_summary.py and cross_match.py. A typo in any one of
them silently drops a section from the report rather than failing, which is why
they live here rather than as literals at each site (AGENTS.md trap 5).

Note that ``parsed_result`` is not uniformly shaped: tsv/csv steps produce
``{"comments": [...], "data": [...]}`` while json steps (BAM Header Parsing)
produce the raw object. ``get_step_data`` returns ``[]`` for the latter; use
``get_step_result`` when you want the whole object.
"""

import logging
from typing import Any, Final

logger = logging.getLogger(__name__)

STEP_BAM_HEADER: Final[str] = "BAM Header Parsing"
STEP_COVERAGE: Final[str] = "Coverage Calculation"
STEP_KESTREL: Final[str] = "Kestrel Genotyping"
STEP_ADVNTR: Final[str] = "adVNTR Genotyping"
STEP_CROSS_MATCH: Final[str] = "Cross-Match Variant Comparison"

STEP_NAMES: Final[frozenset[str]] = frozenset(
    {
        STEP_BAM_HEADER,
        STEP_COVERAGE,
        STEP_KESTREL,
        STEP_ADVNTR,
        STEP_CROSS_MATCH,
    }
)


def get_step(summary: dict[str, Any], step_name: str) -> dict[str, Any] | None:
    """Return the step entry with the given name.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.
        step_name: One of the ``STEP_*`` constants.

    Returns:
        Optional[dict[str, Any]]: The step mapping, or None if absent.
    """
    for step in summary.get("steps", []):
        if step.get("step") == step_name:
            return step
    return None


def get_step_result(summary: dict[str, Any], step_name: str) -> dict[str, Any]:
    """Return a step's ``parsed_result``, or an empty mapping if absent.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.
        step_name: One of the ``STEP_*`` constants.

    Returns:
        dict[str, Any]: The parsed result, or ``{}``.
    """
    step = get_step(summary, step_name)
    if step is None:
        return {}
    result = step.get("parsed_result")
    return result if isinstance(result, dict) else {}


def get_step_data(summary: dict[str, Any], step_name: str) -> list[dict[str, Any]]:
    """Return a step's tabular rows, or an empty list.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.
        step_name: One of the ``STEP_*`` constants.

    Returns:
        list[dict[str, Any]]: The ``data`` rows. Empty when the step is absent
        or its parsed result is not tabular.
    """
    data = get_step_result(summary, step_name).get("data")
    return data if isinstance(data, list) else []
