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
from pathlib import Path
from typing import Any

import pandas as pd

from vntyper.scripts.nomenclature import (
    Nomenclature,
    from_advntr,
    from_kestrel,
    reconcile,
    render,
)
from vntyper.scripts.nomenclature_bam import BamRescuer, from_bam, is_candidate, refine

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


def annotate_kestrel_frame(frame: pd.DataFrame, output_dir: str | Path | None = None) -> pd.DataFrame:
    """Add the nomenclature columns to a Kestrel result frame.

    Args:
        frame: The final Kestrel frame, as written to ``kestrel_result.tsv``.
        output_dir: The Kestrel output directory, when the reads are available. Given
            one, rows the VCF cannot resolve are refined against ``output.bam`` --
            the rescue path that recovers a delins. Omitted, the VCF result stands
            and no BAM is opened.

    Returns:
        pd.DataFrame: A copy with the five columns appended. The negative-run frame
        is returned unchanged: it has no variant, and padding it with five empty
        columns would widen a schema that deliberately differs.
    """
    if frame.empty or "Motifs" not in frame.columns:
        return frame

    annotated = frame.copy()
    rescuer = _open_rescuer(output_dir)
    locus = _read_locus(output_dir)

    try:
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

            merged = reconcile(call)
            # VCF primary, BAM refines. The BAM is consulted only for candidates, so
            # the common path never opens it.
            if rescuer is not None and locus is not None and is_candidate(merged):
                contig, position = locus
                consensus = rescuer.rescue(contig, position)
                if consensus is not None:
                    merged = refine(merged, from_bam(contig, consensus))
            cells.append(_cells(merged))
    finally:
        if rescuer is not None:
            rescuer.close()

    for column in NOMENCLATURE_COLUMNS:
        annotated[column] = [cell[column] for cell in cells]
    return annotated


def _open_rescuer(output_dir: str | Path | None) -> BamRescuer | None:
    """Build a rescuer for a sample, or ``None`` when the reads are unavailable.

    Args:
        output_dir: The Kestrel output directory.

    Returns:
        BamRescuer | None: A rescuer owning one handle for the whole frame. The
        handle is opened lazily on first use, so building this costs nothing when no
        row turns out to be a candidate.
    """
    if output_dir is None:
        return None
    bam = Path(output_dir) / "output.bam"
    if not bam.is_file():
        return None
    return BamRescuer(bam)


def _read_locus(output_dir: str | Path | None) -> tuple[str, int] | None:
    """Read the called locus from ``output.bed``.

    Args:
        output_dir: The Kestrel output directory.

    Returns:
        tuple[str, int] | None: ``(pair name, 1-based position)``, or ``None`` when
        the BED is absent or malformed.
    """
    if output_dir is None:
        return None
    bed = Path(output_dir) / "output.bed"
    if not bed.is_file():
        return None
    fields = bed.read_text().split()
    if len(fields) < 2:
        return None
    try:
        return fields[0], int(fields[1])
    except ValueError:
        logger.debug("output.bed did not carry a numeric position: %r", fields[:2])
        return None


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
