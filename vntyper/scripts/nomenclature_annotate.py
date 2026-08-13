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
    confidence_note,
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
    "Nomenclature_Note",
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
        "Nomenclature_Note": confidence_note(call),
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
            locus = _row_locus(row)
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


def _row_locus(row: pd.Series) -> tuple[str, int] | None:
    """The locus a row was called at, taken from the row itself.

    Deliberately not read from ``output.bed``. That file carries one record per
    frame row, so a single parse of it names the *first* row's locus for every row;
    and its start column is ``POS_fasta - 1`` because BED is half-open, which read
    as a 1-based position shifts the rescue window one base left. ``Motif_fasta``
    and ``POS_fasta`` are what the BED is written from, so taking them directly is
    both correct per row and immune to that off-by-one.

    Args:
        row: One Kestrel result row.

    Returns:
        tuple[str, int] | None: ``(pair name, 1-based position)``, or ``None`` when
        the row does not carry them.
    """
    contig = row.get("Motif_fasta")
    position = row.get("POS_fasta")
    if contig is None or position is None:
        return None
    try:
        return str(contig), int(position)
    except (TypeError, ValueError):
        logger.debug("row carried a non-numeric POS_fasta: %r", position)
        return None


def reconcile_caller_outputs(kestrel_tsv: str | Path, advntr_tsv: str | Path) -> bool:
    """Reconcile the two callers' written results and rewrite both in place.

    Each caller names its own rows as it writes them, which is what puts the columns
    on every surface. But tier A requires **two independent sources agreeing**, and
    no single caller's stage can see the other -- so without this step production
    could never emit a tier-A name, however well the two agreed.

    Both files are rewritten so the report, the summary and the cohort tables all
    read the same reconciled verdict rather than two independent opinions.

    Args:
        kestrel_tsv: Path to ``kestrel_result.tsv``.
        advntr_tsv: Path to ``output_adVNTR_result.tsv``.

    Returns:
        bool: True when both files were read and rewritten.
    """
    kestrel_path, advntr_path = Path(kestrel_tsv), Path(advntr_tsv)
    if not kestrel_path.is_file() or not advntr_path.is_file():
        logger.debug("Cross-caller reconciliation skipped; one of the result files is absent.")
        return False

    try:
        header = [line for line in kestrel_path.read_text().splitlines() if line.startswith("##")]
        kestrel = pd.read_csv(kestrel_path, sep="\t", comment="#", dtype=str)
        advntr = pd.read_csv(advntr_path, sep="\t", dtype=str)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        logger.warning("Cross-caller reconciliation skipped; could not read a result file: %s", error)
        return False

    if kestrel.empty or advntr.empty or "Motifs" not in kestrel.columns:
        return False

    calls: list[Nomenclature] = []
    supports: dict[str, int | None] = {}

    for _, row in kestrel.iterrows():
        if _is_negative(row):
            continue
        try:
            calls.append(from_kestrel(str(row["Motifs"]), int(row["POS"]), str(row["REF"]), str(row["ALT"])))
        except (TypeError, ValueError):
            continue
        supports.setdefault("kestrel_vcf", _as_int(row.get("Estimated_Depth_AlternateVariant")))

    # Every adVNTR event at the locus, not just the first: a sample reporting several
    # simultaneous events is not describing one simple allele, and hiding the rest
    # would let a wrong name be promoted.
    for _, row in advntr.iterrows():
        if _is_negative(row):
            continue
        parsed = from_advntr(str(row.get("Variant", "")))
        if not parsed:
            continue
        calls.extend(parsed)
        depth = _as_int(row.get("NumberOfSupportingReads"))
        existing = supports.get("advntr")
        supports["advntr"] = depth if existing is None else min(existing, depth or existing)

    if not calls:
        return False

    merged = reconcile(*calls, supports=supports)
    cells = _cells(merged)

    for frame, path, is_kestrel in ((kestrel, kestrel_path, True), (advntr, advntr_path, False)):
        updated = frame.copy()
        for column, value in cells.items():
            if column in updated.columns:
                updated.loc[~updated.apply(_is_negative, axis=1), column] = value
        _write_tsv(updated, path, header if is_kestrel else [])

    logger.info("Cross-caller nomenclature reconciled: tier %s.", merged.tier)
    return True


def _as_int(value: object) -> int | None:
    """Read a depth cell as an integer, or ``None`` when it is not one."""
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _write_tsv(frame: pd.DataFrame, path: Path, header: list[str]) -> None:
    """Rewrite a result TSV, preserving its ``##`` header lines.

    Args:
        frame: The frame to write.
        path: Destination.
        header: ``##`` lines to re-emit above the table.
    """
    with path.open("w") as handle:
        if header:
            handle.write("\n".join(header) + "\n")
        frame.to_csv(handle, sep="\t", index=False)


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
