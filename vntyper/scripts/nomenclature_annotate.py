"""Attach MUC1 nomenclature columns to caller result frames.

Kept separate from :mod:`vntyper.scripts.nomenclature` so the translator stays pure
and pandas-free: this is the only place the two worlds meet.

The five columns are added in one pass for both callers, so the result TSVs, the
summary JSON built from those frames, and the HTML report all inherit them from a
single edit.

Nullable by design: ``Ambiguity_Interval`` and ``Repeat_Form`` are written only when
they mean something -- an interval wider than 1 bp, a variant inside a detectable
tract -- and are the empty string otherwise, never a placeholder like ``N/A``.

Research use only.
"""

from __future__ import annotations

import logging
from typing import Any

import pandas as pd

from vntyper.scripts.nomenclature import (
    Nomenclature,
    from_advntr,
    from_kestrel,
    reconcile,
    render,
)

logger = logging.getLogger(__name__)

#: The five columns, in this order, on every surface that carries them.
NOMENCLATURE_COLUMNS: tuple[str, ...] = (
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
)

#: A negative Kestrel run writes a different, 10-column schema whose first column is
#: singular ``Motif``. It carries no variant, so there is nothing to name.
_NEGATIVE_MARKERS = ("None", "Negative", "Not applicable")


def _is_negative(row: pd.Series) -> bool:
    """Is this a placeholder row rather than a called variant?

    Args:
        row: One result row.

    Returns:
        bool: True for the negative-run placeholder of either caller.
    """
    if str(row.get("Confidence", "")) == "Negative":
        return True
    if str(row.get("VID", "")) == "Negative":
        return True
    return str(row.get("Motif", "")) == "None" and "Motifs" not in row.index


def _cells(call: Nomenclature) -> dict[str, Any]:
    """Project a call onto the five columns.

    The displayed name comes from :func:`render`, never from ``call.name``: the tier
    decides what may be shown, and reading the field directly would leak a bare
    number from a tier that is not entitled to state one.

    Args:
        call: The reconciled call.

    Returns:
        dict[str, Any]: The five column values.
    """
    interval = ""
    if call.ambiguity is not None:
        low, high = call.ambiguity
        interval = f"{low}_{high}"

    return {
        "Nomenclature": render(call),
        "Nomenclature_Tier": call.tier,
        "Nomenclature_Flags": ";".join(call.flags),
        "Ambiguity_Interval": interval,
        "Repeat_Form": call.repeat_form or "",
    }


def annotate_kestrel_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Add the nomenclature columns to a Kestrel result frame.

    Args:
        frame: The final Kestrel frame, as written to ``kestrel_result.tsv``.

    Returns:
        pd.DataFrame: A copy with the five columns appended. The negative-run frame
        is returned unchanged: it has no variant, and padding it with five empty
        columns would widen a schema that deliberately differs.
    """
    if frame.empty or "Motifs" not in frame.columns:
        return frame

    annotated = frame.copy()
    cells: list[dict[str, Any]] = []
    for _, row in annotated.iterrows():
        if _is_negative(row):
            cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
            continue
        try:
            call = from_kestrel(str(row["Motifs"]), int(row["POS"]), str(row["REF"]), str(row["ALT"]))
        except (TypeError, ValueError):
            logger.debug("Kestrel row not translatable to nomenclature: %s", dict(row))
            cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
            continue
        cells.append(_cells(reconcile(call)))

    for column in NOMENCLATURE_COLUMNS:
        annotated[column] = [cell[column] for cell in cells]
    return annotated


def annotate_advntr_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Add the nomenclature columns to an adVNTR result frame.

    Unlike Kestrel, adVNTR may report several events for one sample. Each row is
    named on its own here; cross-row reconciliation belongs to whatever combines the
    two callers, which is the only place that can see both.

    Args:
        frame: The processed adVNTR frame, as written to ``output_adVNTR_result.tsv``.

    Returns:
        pd.DataFrame: A copy with the five columns appended, empty on the negative
        placeholder row.
    """
    if frame.empty or "Variant" not in frame.columns:
        return frame

    annotated = frame.copy()
    cells: list[dict[str, Any]] = []
    for _, row in annotated.iterrows():
        if _is_negative(row):
            cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
            continue
        support: int | None
        try:
            support = int(float(row.get("NumberOfSupportingReads", "")))
        except (TypeError, ValueError):
            support = None
        calls = from_advntr(str(row["Variant"]))
        if not calls:
            cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
            continue
        cells.append(_cells(reconcile(*calls, support=support)))

    for column in NOMENCLATURE_COLUMNS:
        annotated[column] = [cell[column] for cell in cells]
    return annotated
