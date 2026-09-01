"""Pandas-to-BAM adapter for nomenclature rescue and identity evidence."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from vntyper.scripts.molecular_identity import IdentityTranslation
from vntyper.scripts.nomenclature import Nomenclature
from vntyper.scripts.nomenclature_bam import BamRescuer, from_bam
from vntyper.scripts.nomenclature_bam_evidence import (
    BamCandidateBinding,
    BamIdentityLocus,
    BamLocusEvidence,
)

logger = logging.getLogger(__name__)


def open_rescuer(output_dir: str | Path | None) -> BamRescuer | None:
    """Build a rescuer for a sample, or ``None`` without haplotype records.

    Args:
        output_dir: Kestrel output directory.

    Returns:
        A lazily opened sample rescuer, or ``None`` when ``output.bam`` is absent.
    """
    if output_dir is None:
        return None
    bam = Path(output_dir) / "output.bam"
    if not bam.is_file():
        return None
    return BamRescuer(bam)


def row_haplotype_call(row: pd.Series, rescuer: BamRescuer) -> tuple[Nomenclature | None, int | None]:
    """Return the legacy BAM call and record support for one Kestrel row.

    Args:
        row: One Kestrel result row.
        rescuer: Open sample rescuer.

    Returns:
        ``(call, supporting records)`` or ``(None, None)`` when unavailable.
    """
    locus = row_locus(row)
    if locus is None:
        return None, None
    contig, position = locus
    consensus = rescuer.rescue(contig, position)
    if consensus is None:
        return None, None
    return from_bam(contig, consensus), consensus.supporting_haplotype_records


def row_haplotype_identity_call(
    row: pd.Series,
    rescuer: BamRescuer,
    component: Any,
) -> tuple[Nomenclature | None, int | None, IdentityTranslation | None, BamLocusEvidence | None]:
    """Resolve one row and retain its complete candidate-bound BAM evidence.

    Args:
        row: Current-run Kestrel row carrying A3 internal metadata.
        rescuer: Open sample rescuer.
        component: Stage-bound complete-context translation authority.

    Returns:
        ``(call, support, translation, evidence)``. Translation is present only
        for an exact winner proven by all complete supporting records; evidence is
        retained even when the compatibility vote abstains.
    """
    from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells
    from vntyper.scripts.molecular_identity_translation import bind_bam_translation

    locus = row_locus(row)
    if locus is None:
        return None, None, None, None
    contig, position = locus
    persisted = parse_selected_candidate_cells(row)
    motifs = str(row.get("Motifs", ""))
    pair = str(row.get("Motif_sequence", ""))
    if len(pair) != 120 or not set(pair) <= frozenset("ACGT"):
        consensus = rescuer.rescue(contig, position)
        if consensus is None:
            return None, None, None, None
        return from_bam(contig, consensus), consensus.supporting_haplotype_records, None, None
    identity_locus = BamIdentityLocus(
        motifs,
        pair,
        (BamCandidateBinding(persisted.selected_observation_ordinal, persisted.translation),),
    )
    consensus, evidence = rescuer.rescue_with_identity_evidence(contig, position, identity_locus, component)
    if consensus is None:
        return None, None, None, evidence
    bound_translation = None
    if consensus.bound_identity is not None:
        bound_translation = bind_bam_translation(persisted.translation, consensus.bound_identity)
    return from_bam(contig, consensus), consensus.supporting_haplotype_records, bound_translation, evidence


def row_locus(row: pd.Series) -> tuple[str, int] | None:
    """Return the row-owned pair name and one-based position.

    Args:
        row: One Kestrel result row.

    Returns:
        ``(pair name, position)`` or ``None`` when absent or malformed.
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
