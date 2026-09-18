"""
vntyper/scripts/cohort_summary_parsing.py

Module Purpose:
---------------
Parse pipeline summary records for cohort aggregation.

Extracts algorithm rows and per-sample statistics (runtime, version, assembly,
pipeline, coverage, and optional length estimation metrics) from
``pipeline_summary.json``.

Extracted from ``cohort_inputs.py`` to keep pure parsing logic decoupled from
filesystem discovery and under file-length guidelines.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any

from vntyper.scripts.molecular_identity_presentation import identity_compatible_result_row
from vntyper.scripts.report_formatting import is_empty_result_row
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_BAM_HEADER, STEP_COVERAGE, STEP_KESTREL

logger = logging.getLogger(__name__)

PIPELINE_SUMMARY_FILENAME = "pipeline_summary.json"


def parse_pipeline_summary(summary: dict[str, Any]) -> tuple[list[dict], list[dict], dict[str, Any]]:
    """Extract the cohort's three inputs from one parsed ``pipeline_summary.json``.

    Steps this cohort does not consume are ignored, and a step recorded more than once
    leaves the last occurrence in place.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.

    Returns:
        tuple[list[dict], list[dict], dict[str, Any]]: The Kestrel rows, the adVNTR
        rows, and the per-sample statistics (``runtime``, ``version``, ``assembly``,
        ``pipeline``, ``coverage``, plus length estimation stats if recorded).

    Raises:
        ValueError: If ``pipeline_start`` or ``pipeline_end`` is present but is not an
            ISO-8601 timestamp. :func:`load_pipeline_summary_for_sample` catches this
            and drops the sample.
    """
    kestrel_data: list[dict] = []
    advntr_data: list[dict] = []
    additional_stats: dict[str, Any] = {}

    # Compute runtime from top-level timestamps if available
    pipeline_start = summary.get("pipeline_start")
    pipeline_end = summary.get("pipeline_end")
    if pipeline_start and pipeline_end:
        start_dt = datetime.fromisoformat(pipeline_start)
        end_dt = datetime.fromisoformat(pipeline_end)
        runtime_sec = (end_dt - start_dt).total_seconds()
        additional_stats["runtime"] = f"{runtime_sec:.2f} seconds"
    else:
        additional_stats["runtime"] = "N/A"

    # Pipeline version from top-level field
    additional_stats["version"] = summary.get("version", "N/A")

    # Initialize defaults for assembly, pipeline and coverage
    additional_stats["assembly"] = "N/A"
    additional_stats["pipeline"] = "N/A"
    additional_stats["coverage"] = {}
    schema_version = summary.get("schema_version")

    for step in summary.get("steps", []):
        if step.get("step") == STEP_KESTREL:
            rows = step.get("parsed_result", {}).get("data", [])
            kestrel_data = [
                identity_compatible_result_row(
                    row,
                    schema_version=schema_version,
                    positive=not is_empty_result_row(row),
                )
                for row in rows
            ]
        elif step.get("step") == STEP_ADVNTR:
            rows = step.get("parsed_result", {}).get("data", [])
            advntr_data = [
                identity_compatible_result_row(
                    row,
                    schema_version=schema_version,
                    positive=row.get("VID") != "Negative",
                )
                for row in rows
            ]
        elif step.get("step") == STEP_BAM_HEADER:
            parsed = step.get("parsed_result", {})
            additional_stats["assembly"] = parsed.get("assembly_text", "N/A")
            additional_stats["pipeline"] = parsed.get("alignment_pipeline", "N/A")
        elif step.get("step") == STEP_COVERAGE:
            parsed = step.get("parsed_result", {})
            data_list = parsed.get("data", [])
            if data_list:
                additional_stats["coverage"] = data_list[0]

    # Length estimation metrics and warnings if present in summary
    if "estimated_total_repeat_count" in summary:
        additional_stats["estimated_total_repeat_count"] = summary.get("estimated_total_repeat_count")
    if "length_estimation_warnings" in summary:
        warnings = summary.get("length_estimation_warnings")
        if isinstance(warnings, list):
            additional_stats["length_warning"] = ";".join(warnings) if warnings else "none"
        else:
            additional_stats["length_warning"] = str(warnings)

    return kestrel_data, advntr_data, additional_stats
