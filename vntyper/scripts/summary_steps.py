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

A step has **three** states and not two, which is what :func:`get_step_state`
exists to say. "Recorded" is not "produced a readable result": ``record_step``
writes the record either way, and when the result file is absent it swallows the
``FileNotFoundError`` into ``md5sum=None`` and lets ``parse_tsv`` turn the failure
into a comment on an empty ``data`` list -- which is exactly what a stage that
legitimately found nothing produces (#212). Every consumer that asks only whether
the step is present therefore renders a failed stage as a negative result.
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


#: The step is not in the summary at all: the stage was never asked to run.
STEP_ABSENT: Final[str] = "absent"

#: The step is recorded and its result was read. An empty result is a result: a
#: ``kestrel_result.tsv`` carrying a header and no rows is what a run that genotyped
#: and called nothing writes, and it is read perfectly well.
STEP_READ: Final[str] = "read"

#: The step is recorded and its result could not be read. Distinct from both of the
#: above and reportable as neither: the stage ran, so it was not "not performed", and
#: nothing was read, so nothing about the sample was established.
STEP_UNREADABLE: Final[str] = "unreadable"


def get_step_state(summary: dict[str, Any], step_name: str) -> str:
    """Say which of the three states a step is in.

    The signals are structural, and each is a shape ``summary.record_step`` writes:

    * ``result_file_missing`` -- added when the path the stage was supposed to have
      written does not exist (#212). Until now nothing read it;
    * a ``parsed_result`` that is not a mapping -- ``record_step`` initialises the key
      to ``None`` and only parsing replaces it, so ``None`` (or the key's absence in a
      hand-built summary) means no result was parsed at all;
    * a ``parsed_result`` carrying an ``error`` key -- ``record_step``'s parse-failure
      and unsupported-file-type shapes.

    ``parse_tsv``'s own failure path is deliberately **not** matched: it records the
    failure as a *comment* on an otherwise ordinary result, and recognising it would
    mean string-matching "Error parsing TSV file" against text a legitimate comment
    could carry. The case that produces it is the missing file, which the first signal
    already covers.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.
        step_name: One of the ``STEP_*`` step-name constants.

    Returns:
        str: :data:`STEP_ABSENT`, :data:`STEP_UNREADABLE` or :data:`STEP_READ`.
    """
    step = get_step(summary, step_name)
    if step is None:
        return STEP_ABSENT
    if step.get("result_file_missing"):
        logger.warning("Step %r recorded a missing result file; its result is not a negative.", step_name)
        return STEP_UNREADABLE
    result = step.get("parsed_result")
    if not isinstance(result, dict) or "error" in result:
        logger.warning("Step %r recorded no readable result; its result is not a negative.", step_name)
        return STEP_UNREADABLE
    return STEP_READ


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


def get_step_comments(summary: dict[str, Any], step_name: str) -> list[str]:
    """Return a step's result-file comment lines, or an empty list.

    ``summary.parse_tsv`` records every ``#`` line of a step's TSV here, with the leading
    hashes stripped. #266's below-reporting-floor note travels on this channel: it is a
    banner line, so it reaches the report without ever appearing in :func:`get_step_data`'s
    rows, and no consumer reading the table can mistake it for a call.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.
        step_name: One of the ``STEP_*`` constants.

    Returns:
        list[str]: The comment lines. Empty when the step is absent, recorded none, or
        recorded something that is not a list.
    """
    comments = get_step_result(summary, step_name).get("comments")
    return comments if isinstance(comments, list) else []


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
