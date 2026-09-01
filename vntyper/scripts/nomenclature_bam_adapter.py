"""Pandas-to-BAM adapter for nomenclature rescue and identity evidence."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from vntyper.scripts.molecular_identity import IdentityTranslation
from vntyper.scripts.nomenclature import Nomenclature, load_nomenclature_config, reconcile
from vntyper.scripts.nomenclature_bam import BamRescuer, from_bam, is_candidate
from vntyper.scripts.nomenclature_bam_evidence import (
    BamCandidateBinding,
    BamIdentityLocus,
)
from vntyper.scripts.nomenclature_bam_replay import BamReplayArtifact, BamReplayLocus

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


@dataclass(frozen=True)
class BamIdentityRowCandidate:
    """One persisted selected candidate and its source row position."""

    row_index: int
    observation_ordinal: int
    translation: IdentityTranslation


@dataclass(frozen=True)
class BamIdentityCandidateGroup:
    """Candidates sharing one exact BAM locus and complete pair context."""

    candidates: tuple[BamIdentityRowCandidate, ...]
    locus: tuple[str, int] | None
    identity_locus: BamIdentityLocus | None

    @property
    def candidate_observation_ordinals(self) -> tuple[int, ...]:
        """Return the stable sorted candidate ordinals owned by this group."""
        return tuple(candidate.observation_ordinal for candidate in self.candidates)


@dataclass(frozen=True)
class BamIdentityGroupProjection:
    """One grouped BAM pass projected back onto source-row alignment."""

    calls: tuple[Nomenclature | None, ...]
    supports: tuple[int | None, ...]
    translations: tuple[IdentityTranslation | None, ...]
    artifact: BamReplayArtifact


def group_identity_rows(rows: list[pd.Series]) -> tuple[BamIdentityCandidateGroup, ...]:
    """Group persisted candidates by exact BAM locus and complete pair context.

    Args:
        rows: Non-negative current-run Kestrel rows in public file order.

    Returns:
        Deterministically ordinal-sorted groups. Rows without complete locus/context
        remain singleton unavailable groups so replay coverage stays complete.

    Raises:
        ValueError: If persisted candidate metadata or group ordinals are invalid.
    """
    from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells

    grouped: dict[tuple[str, int, str, str], list[BamIdentityRowCandidate]] = {}
    incomplete: list[BamIdentityCandidateGroup] = []
    for row_index, row in enumerate(rows):
        persisted = parse_selected_candidate_cells(row)
        candidate = BamIdentityRowCandidate(
            row_index,
            persisted.selected_observation_ordinal,
            persisted.translation,
        )
        locus = row_locus(row)
        motifs = str(row.get("Motifs", ""))
        pair = str(row.get("Motif_sequence", ""))
        if locus is None or not motifs or len(pair) != 120 or not set(pair) <= frozenset("ACGT"):
            incomplete.append(BamIdentityCandidateGroup((candidate,), locus, None))
            continue
        grouped.setdefault((*locus, motifs, pair), []).append(candidate)

    complete: list[BamIdentityCandidateGroup] = []
    for (contig, position, motifs, pair), candidates in grouped.items():
        ordered = tuple(sorted(candidates, key=lambda candidate: candidate.observation_ordinal))
        identity_locus = BamIdentityLocus(
            motifs,
            pair,
            tuple(BamCandidateBinding(candidate.observation_ordinal, candidate.translation) for candidate in ordered),
        )
        complete.append(BamIdentityCandidateGroup(ordered, (contig, position), identity_locus))
    return tuple(sorted((*complete, *incomplete), key=lambda group: group.candidate_observation_ordinals))


def observe_identity_groups(
    groups: tuple[BamIdentityCandidateGroup, ...],
    row_count: int,
    eligible_row_indices: frozenset[int],
    rescuer: BamRescuer | None,
    component: Any,
) -> BamIdentityGroupProjection:
    """Consult each eligible exact group once and project results to source rows.

    Args:
        groups: Validated exact candidate groups.
        row_count: Number of source Kestrel rows.
        eligible_row_indices: Rows the legacy path would have consulted/refined.
        rescuer: Shared sample rescuer, or ``None`` when BAM is unavailable.
        component: Stage-bound complete-context translation authority.

    Returns:
        Aligned legacy calls/support/translations and one joint replay locus per
        candidate group.

    Raises:
        ValueError: If row alignment or group coverage is inconsistent.
    """
    from vntyper.scripts.molecular_identity_translation import bind_bam_translation
    from vntyper.scripts.nomenclature_bam_replay import validate_bam_replay_artifact_ordinals

    if isinstance(row_count, bool) or not isinstance(row_count, int) or row_count < 0:
        raise ValueError("BAM identity projection row count must be a non-negative integer")
    if not isinstance(eligible_row_indices, frozenset) or any(
        isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < row_count
        for index in eligible_row_indices
    ):
        raise ValueError("BAM identity eligible row indices must be valid source-row indices")
    calls: list[Nomenclature | None] = [None] * row_count
    supports: list[int | None] = [None] * row_count
    translations: list[IdentityTranslation | None] = [None] * row_count
    replay: list[BamReplayLocus] = []
    all_ordinals: list[int] = []
    all_row_indices: list[int] = []
    for group in groups:
        ordinals = group.candidate_observation_ordinals
        all_ordinals.extend(ordinals)
        all_row_indices.extend(candidate.row_index for candidate in group.candidates)
        wanted = tuple(candidate for candidate in group.candidates if candidate.row_index in eligible_row_indices)
        if not wanted:
            replay.append(BamReplayLocus(ordinals, "not-consulted", None))
            continue
        if rescuer is None or group.locus is None:
            replay.append(BamReplayLocus(ordinals, "unavailable", None))
            continue
        contig, position = group.locus
        if group.identity_locus is None:
            consensus = rescuer.rescue(contig, position)
            replay.append(BamReplayLocus(ordinals, "unavailable", None))
            if consensus is not None:
                call = from_bam(contig, consensus)
                for candidate in wanted:
                    calls[candidate.row_index] = call
                    supports[candidate.row_index] = consensus.supporting_haplotype_records
            continue
        consensus, evidence = rescuer.rescue_with_identity_evidence(contig, position, group.identity_locus, component)
        if evidence is None:
            replay.append(BamReplayLocus(ordinals, "unavailable", None))
            continue
        replay.append(BamReplayLocus(ordinals, "observed", evidence))
        if consensus is None:
            continue
        call = from_bam(contig, consensus)
        support = consensus.supporting_haplotype_records
        for candidate in wanted:
            calls[candidate.row_index] = call
            supports[candidate.row_index] = support
            if consensus.bound_identity is not None:
                translations[candidate.row_index] = bind_bam_translation(
                    candidate.translation,
                    consensus.bound_identity,
                )
    if tuple(sorted(all_row_indices)) != tuple(range(row_count)):
        raise ValueError("BAM identity groups must cover every source row exactly once")
    artifact = BamReplayArtifact(tuple(replay))
    validate_bam_replay_artifact_ordinals(artifact, tuple(sorted(all_ordinals)))
    return BamIdentityGroupProjection(tuple(calls), tuple(supports), tuple(translations), artifact)


def haplotype_calls(
    kestrel_rows: list[pd.Series],
    kestrel_dir: str | Path,
    vcf_calls: list[Nomenclature],
    advntr_calls: list[Nomenclature],
    supports: dict[str, int | None],
    *,
    bam_translations: list[IdentityTranslation | None] | None = None,
    identity_component: Any = None,
    bam_replay_loci: list[BamReplayLocus] | None = None,
) -> list[Nomenclature | None]:
    """Resolve grouped BAM evidence only where caller evidence leaves a candidate.

    Args:
        kestrel_rows: Non-negative Kestrel rows in file order.
        kestrel_dir: Directory holding ``output.bam``.
        vcf_calls: Named Kestrel VCF calls.
        advntr_calls: Named adVNTR calls.
        supports: Per-source quantities, updated with BAM record support.
        bam_translations: Optional aligned output for complete-candidate bindings.
        identity_component: Optional stage-bound complete-context translator.
        bam_replay_loci: Optional output with one replay state per exact group.

    Returns:
        One call per Kestrel row for identity-aware input; otherwise the legacy
        aligned calls or an empty list when BAM consultation is unnecessary.

    Raises:
        ValueError: If replay output is requested without aligned translations.
    """
    if bam_replay_loci is not None and bam_translations is None:
        raise ValueError("BAM replay output requires aligned identity translations")
    consultation_needed = is_candidate(reconcile(*vcf_calls, *advntr_calls, supports=supports))
    if bam_translations is not None:
        groups = group_identity_rows(kestrel_rows)
        if identity_component is None:
            from vntyper.scripts.identity_candidates import translation_component_from_config

            identity_component = translation_component_from_config(load_nomenclature_config())
        rescuer = open_rescuer(kestrel_dir) if consultation_needed else None
        try:
            projection = observe_identity_groups(
                groups,
                len(kestrel_rows),
                frozenset(range(len(kestrel_rows))) if consultation_needed else frozenset(),
                rescuer,
                identity_component,
            )
        finally:
            if rescuer is not None:
                rescuer.close()
        bam_translations.extend(projection.translations)
        if bam_replay_loci is not None:
            bam_replay_loci.extend(projection.artifact.loci)
        projected_calls = list(projection.calls)
        for call, support in zip(projected_calls, projection.supports, strict=True):
            if call is None:
                continue
            existing = supports.get("kestrel_bam")
            supports["kestrel_bam"] = support if existing is None else min(existing, support or existing)
        return projected_calls

    if not consultation_needed:
        return []
    rescuer = open_rescuer(kestrel_dir)
    if rescuer is None:
        return []

    calls: list[Nomenclature | None] = []
    try:
        for row in kestrel_rows:
            call, support = row_haplotype_call(row, rescuer)
            calls.append(call)
            if call is None:
                continue
            existing = supports.get("kestrel_bam")
            supports["kestrel_bam"] = support if existing is None else min(existing, support or existing)
    finally:
        rescuer.close()
    return calls


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
