"""
vntyper/scripts/summary_flattening.py

Module Purpose:
---------------
Flatten one ``pipeline_summary.json`` mapping into the cells of the operator-facing
``pipeline_summary.<csv|tsv>`` and ``pipeline_summary_rows.<csv|tsv>`` tables.

Three pure functions and no I/O (``summary.py`` keeps the writers):

* :func:`run_columns` -- the run's own provenance (schema version, decision policy, the
  six ``decision_profile_*`` fields, version, inputs, sample, reference selection, the
  resolved region, the Kestrel counting mode, the adVNTR model, start and end) as
  ``run_``-prefixed text cells, so every exported row can say which run produced it.
* :func:`step_rows` -- one row per recorded step. A single-row result explodes into
  ``parsed_result_data_<field>`` columns; any other row count is recorded as
  ``parsed_result_n_rows``; comment lines join with ``" | "``; a json-typed result
  flattens with ``_``.
* :func:`long_rows` -- every result row of every step as ``step``, ``row_index``,
  ``field``, ``value`` records, so the multi-row adVNTR and cross-match tables pivot in
  one line.

Every cell is text and none of them is JSON: the previous flattener embedded a list of
row mappings as ``"; "``-joined ``json.dumps`` text, which is what #119 reported.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any

logger = logging.getLogger(__name__)

#: Prefix of every run-provenance column.
RUN_PREFIX = "run"

#: Prefix of every column derived from a step's ``parsed_result``.
PARSED_RESULT_PREFIX = "parsed_result"

#: Fixed columns of every :func:`step_rows` row, in this order.
STEP_COLUMNS: tuple[str, ...] = (
    "step",
    "start",
    "end",
    "command",
    "result_file",
    "file_type",
    "md5sum",
    "result_file_missing",
)

#: Columns of every :func:`long_rows` row, in this order.
LONG_COLUMNS: tuple[str, ...] = ("step", "row_index", "field", "value")

KEY_SEPARATOR = "_"
SCALAR_LIST_SEPARATOR = "; "
COMMENT_SEPARATOR = " | "

_STEPS_KEY = "steps"


def _cell(value: Any) -> str:
    """Render one scalar as a text cell; a JSON null is a blank cell, not ``"None"``."""
    return "" if value is None else str(value)


def _is_scalar(value: Any) -> bool:
    return value is None or isinstance(value, (str, int, float, bool))


def _flatten(mapping: Mapping[str, Any], prefix: str) -> dict[str, str]:
    """Flatten ``mapping`` under ``prefix`` with ``_``.

    A scalar becomes a text cell, a list of scalars joins with ``"; "``, a nested mapping
    recurses, and anything else (a list of mappings, a mixed list) is omitted: there is no
    flat cell for it, and embedding it as JSON text is the defect this module replaces.
    """
    cells: dict[str, str] = {}
    for key, value in mapping.items():
        column = f"{prefix}{KEY_SEPARATOR}{key}"
        if _is_scalar(value):
            cells[column] = _cell(value)
        elif isinstance(value, Mapping):
            cells.update(_flatten(value, column))
        elif isinstance(value, list) and all(_is_scalar(item) for item in value):
            cells[column] = SCALAR_LIST_SEPARATOR.join(_cell(item) for item in value)
        else:
            logger.debug(f"Omitting {column}: a {type(value).__name__} has no flat cell.")
    return cells


def run_columns(summary: Mapping[str, Any]) -> dict[str, str]:
    """Return the run's provenance as ``run_``-prefixed text cells.

    Every top-level key except ``steps`` is covered, in the summary's own order: the
    scalars ``schema_version``, ``decision_policy``, ``advntr_evidence_digest``, the six
    ``decision_profile_*`` fields, ``pipeline_start``, ``pipeline_end``, ``version``,
    ``sample_name``, ``sample_name_is_explicit``, the four ``reference_*`` fields,
    ``region_resolved`` and ``kestrel_counting_mode``; the mappings ``input_files`` and
    ``advntr_model`` flatten with ``_``.

    Args:
        summary: A ``pipeline_summary.json`` mapping.

    Returns:
        dict[str, str]: Column name to text cell. ``steps`` has its own tables and is
        never a run column, even when it is empty.
    """
    return _flatten({key: value for key, value in summary.items() if key != _STEPS_KEY}, RUN_PREFIX)
