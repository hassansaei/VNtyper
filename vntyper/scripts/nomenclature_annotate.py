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

#: The columns, in this order, on every surface that carries them.
#:
#: ``Nomenclature`` is the reconciled verdict. The two per-caller columns beside it
#: are what each caller reported on its own, kept so a disagreement stays legible:
#: collapsing both files to a single verdict destroyed the evidence a reader needs
#: in order to weigh that verdict at all. They are named absolutely rather than
#: relatively ("the other caller"), so a row means the same thing in either file and
#: in a cohort table that merges them.
NOMENCLATURE_COLUMNS: tuple[str, ...] = (
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
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


def _cells(call: Nomenclature, kestrel: str = "", advntr: str = "") -> dict[str, Any]:
    """Project a call onto the nomenclature columns.

    The displayed name comes from :func:`render`, never from ``call.name``: the tier
    decides what may be shown, and reading the field directly would leak a bare
    number from a tier that is not entitled to state one.

    Args:
        call: The reconciled call.
        kestrel: What Kestrel reported, or ``""`` when it reported nothing.
        advntr: What adVNTR reported, or ``""`` when it did not run or reported
            nothing. adVNTR is an optional module, so empty is the ordinary case.

    Returns:
        dict[str, Any]: The column values.
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
        "Nomenclature_Kestrel": kestrel,
        "Nomenclature_adVNTR": advntr,
    }


def _summarise(calls: list[Nomenclature]) -> str:
    """One caller's verdict as display text.

    Distinct renderings joined, so a caller reporting several events says so rather
    than being reduced to its first one. Empty means the caller contributed nothing,
    which is different from contributing a call it could not name -- that renders as
    a statement of what is known.

    Args:
        calls: Every call from one caller.

    Returns:
        str: The joined renderings, or ``""`` when there are none.
    """
    seen: list[str] = []
    for call in calls:
        text = render(call)
        if text and text not in seen:
            seen.append(text)
    return ";".join(seen)


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
            call = _kestrel_row_call(row)
            if call is None:
                cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
                continue

            merged = reconcile(call)
            # VCF primary, BAM refines. The BAM is consulted only for candidates, so
            # the common path never opens it.
            bam_call = None
            if rescuer is not None and is_candidate(merged):
                bam_call, _ = _row_read_call(row, rescuer)
            merged = refine(merged, bam_call)
            # adVNTR is optional and has not run at this point, so its column is
            # empty; the cross-caller stage fills it in when it does run.
            cells.append(_cells(merged, kestrel=render(merged)))
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


def _kestrel_row_call(row: pd.Series) -> Nomenclature | None:
    """Translate one Kestrel VCF row, or ``None`` when it cannot be translated.

    Args:
        row: One Kestrel result row.

    Returns:
        Nomenclature | None: The named record.
    """
    try:
        return from_kestrel(str(row["Motifs"]), int(row["POS"]), str(row["REF"]), str(row["ALT"]))
    except (TypeError, ValueError):
        logger.debug("Kestrel row not translatable to nomenclature: %s", dict(row))
        return None


def _row_read_call(row: pd.Series, rescuer: BamRescuer) -> tuple[Nomenclature | None, int | None]:
    """What the reads say at one Kestrel row's locus.

    Args:
        row: One Kestrel result row.
        rescuer: The open rescuer for this sample.

    Returns:
        tuple: ``(call, supporting reads)``; ``(None, None)`` when the locus is
        unreadable or the reads carry no length-changing edit. The read count is
        returned alongside because the tier gate must weigh an agreement by the reads
        actually under it, not by an unrelated source's depth.
    """
    locus = _row_locus(row)
    if locus is None:
        return None, None
    contig, position = locus
    consensus = rescuer.rescue(contig, position)
    if consensus is None:
        return None, None
    return from_bam(contig, consensus), consensus.support


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


def reconcile_caller_outputs(
    kestrel_tsv: str | Path,
    advntr_tsv: str | Path,
    kestrel_dir: str | Path | None = None,
) -> bool:
    """Reconcile the two callers' written results and rewrite both in place.

    Each caller names its own rows as it writes them, which is what puts the columns
    on every surface. But this is the only place that can see *both*, so it is the
    only place that can do three things no caller stage can:

    1. Promote to tier A, which requires two independent sources agreeing.
    2. Let two independent sources outvote a third -- the reason ``insG`` families
       are recovered rather than lost to a Kestrel placement one base 3' of truth.
    3. Record what each caller said, so a disagreement stays legible instead of
       being collapsed into one verdict.

    adVNTR is an optional module. When it has not run this step never fires, and the
    Kestrel result stands exactly as its own stage wrote it.

    Args:
        kestrel_tsv: Path to ``kestrel_result.tsv``.
        advntr_tsv: Path to ``output_adVNTR_result.tsv``.
        kestrel_dir: Directory holding ``output.bam``. Defaults to the directory the
            Kestrel TSV sits in, which is where the pipeline writes it.

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

    supports: dict[str, int | None] = {}
    kestrel_rows = [row for _, row in kestrel.iterrows() if not _is_negative(row)]

    vcf_calls: list[Nomenclature] = []
    for row in kestrel_rows:
        call = _kestrel_row_call(row)
        if call is None:
            continue
        vcf_calls.append(call)
        supports.setdefault("kestrel_vcf", _as_int(row.get("Estimated_Depth_AlternateVariant")))

    advntr_calls = _advntr_calls(advntr, supports)

    if not vcf_calls and not advntr_calls:
        return False

    bam_calls = _read_calls(kestrel_rows, kestrel_dir or kestrel_path.parent, vcf_calls, advntr_calls, supports)

    # Order matters: the Kestrel VCF is offered first so that whenever nothing
    # outvotes it, it is the call that stands. adVNTR is optional; Kestrel is not.
    merged = reconcile(*vcf_calls, *bam_calls, *advntr_calls, supports=supports)
    for bam_call in bam_calls:
        # Still applied after the vote, because it carries one rule the vote cannot:
        # a delins is unrepresentable in Kestrel's VCF, so a delins seen in the reads
        # is better evidence than the closest shape the VCF could write.
        merged = refine(merged, bam_call)

    cells = _cells(
        merged,
        kestrel=_summarise([refine(reconcile(call), bam) for call, bam in _pair_rows(vcf_calls, bam_calls)]),
        advntr=_summarise(advntr_calls),
    )

    for frame, path, is_kestrel in ((kestrel, kestrel_path, True), (advntr, advntr_path, False)):
        updated = frame.copy()
        keep = ~updated.apply(_is_negative, axis=1)
        for column, value in cells.items():
            if column not in updated.columns:
                # The per-caller columns are this step's own output. A file written
                # before they existed must gain them rather than silently drop the
                # information they carry.
                updated[column] = ""
            updated.loc[keep, column] = value
        _write_tsv(updated, path, header if is_kestrel else [])

    logger.info(
        "Cross-caller nomenclature reconciled: tier %s (kestrel=%r advntr=%r).",
        merged.tier,
        cells["Nomenclature_Kestrel"],
        cells["Nomenclature_adVNTR"],
    )
    return True


def _pair_rows(
    vcf_calls: list[Nomenclature],
    bam_calls: list[Nomenclature],
) -> list[tuple[Nomenclature, Nomenclature | None]]:
    """Pair each Kestrel VCF call with the read consensus at its own locus.

    Args:
        vcf_calls: One per translatable Kestrel row, in row order.
        bam_calls: The read consensus per row, same order, short when the reads were
            not consulted.

    Returns:
        list: ``(vcf call, read call or None)`` pairs.
    """
    return [(call, bam_calls[index] if index < len(bam_calls) else None) for index, call in enumerate(vcf_calls)]


def _advntr_calls(advntr: pd.DataFrame, supports: dict[str, int | None]) -> list[Nomenclature]:
    """Every adVNTR event at the locus, with the depth of the thinnest one.

    Not just the first: a sample reporting several simultaneous events is not
    describing one simple allele, and hiding the rest would let a wrong name be
    promoted.

    Args:
        advntr: The adVNTR result frame.
        supports: Per-source depths, updated in place.

    Returns:
        list[Nomenclature]: The named events.
    """
    calls: list[Nomenclature] = []
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
    return calls


def _read_calls(
    kestrel_rows: list[pd.Series],
    kestrel_dir: str | Path,
    vcf_calls: list[Nomenclature],
    advntr_calls: list[Nomenclature],
    supports: dict[str, int | None],
) -> list[Nomenclature]:
    """The reads, as a third source, where the two callers leave something open.

    The reads are what separate a Kestrel misplacement from a real disagreement, and
    they are the only evidence available when the two callers conflict. They are
    consulted only when the callers do not already settle the locus, so a sample
    where the two agree still never opens the BAM.

    Args:
        kestrel_rows: The non-negative Kestrel rows, in file order.
        kestrel_dir: Directory holding ``output.bam``.
        vcf_calls: The Kestrel VCF calls.
        advntr_calls: The adVNTR calls.
        supports: Per-source depths, updated in place with the read count.

    Returns:
        list[Nomenclature]: One call per row the reads could speak to.
    """
    if not is_candidate(reconcile(*vcf_calls, *advntr_calls, supports=supports)):
        return []

    rescuer = _open_rescuer(kestrel_dir)
    if rescuer is None:
        return []

    calls: list[Nomenclature] = []
    try:
        for row in kestrel_rows:
            call, support = _row_read_call(row, rescuer)
            if call is None:
                continue
            calls.append(call)
            existing = supports.get("kestrel_bam")
            # The weakest read consensus sets the depth, as for every other source:
            # an agreement is only as strong as its thinnest contributing evidence.
            supports["kestrel_bam"] = support if existing is None else min(existing, support or existing)
    finally:
        rescuer.close()
    return calls


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
        merged = reconcile(*calls, support=support)
        # Kestrel's column is filled by the cross-caller stage, which is the only
        # place that can see it.
        cells.append(_cells(merged, advntr=render(merged)))

    for column in NOMENCLATURE_COLUMNS:
        annotated[column] = [cell[column] for cell in cells]
    return annotated
