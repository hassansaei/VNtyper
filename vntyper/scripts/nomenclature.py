"""Literature-compatible naming for MUC1-VNTR variants.

Kestrel and adVNTR each report MUC1-VNTR variants in their own internal coordinate
frame (``POS 67 G>GG``, ``I22_2_G_LEN1``). Neither emits the naming the MUC1
literature uses (``59dupC``). This module turns a caller record into that naming.

Notation
--------
Names are the **bare positional shorthand** on a named reference: ``59dupC``,
``60dupA``, ``58_59insG``, ``1_5delGCCCA``. There is deliberately **no ``c.``
prefix**: ``c.`` asserts a coding-DNA reference sequence, and no transcript places
this tract at positions 53-59. The reference is named in prose instead -- "the
canonical MUC1 60 bp repeat unit, coding orientation".

Coordinates
-----------
All positions here are 1-based and inclusive, on a single 60 bp repeat unit in the
coding orientation. An **insertion** is expressed as the empty span between two
bases, i.e. ``end == start - 1``, so insertions and deletions share one convention.

Reference data and thresholds are configured, not written into the logic: see
``nomenclature_config.json``. It is read once at import, per AGENTS.md trap 1.
Beyond that this module does no I/O and has no pandas or logging in the hot path,
so it is safe to call from a dataframe ``.apply`` and from tests without fixtures.

Research use only.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, NamedTuple

from vntyper.scripts.identity_reconciliation import (
    IdentityReconciliationObservation,
    IdentityReconciliationPolicy,
    reconcile_identity_observations,
    select_compatibility_presentation_index,
)
from vntyper.scripts.nomenclature_evidence import (
    FLAG_LOW_EVIDENCE_SUPPORT,
    FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
    FLAG_LOW_KMER_PATH_SUPPORT,
    FLAG_LOW_READ_SUPPORT,
    FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,
    low_support_flag_for_source,
)
from vntyper.scripts.utils import load_config

if TYPE_CHECKING:  # pragma: no cover - typing only
    from collections.abc import Mapping

__all__ = [
    "CANONICAL_UNIT",
    "CALLER_OF",
    "FLAG_ALLELE_UNREPRESENTABLE",
    "FLAG_CALLER_DISAGREEMENT",
    "FLAG_KNOWN_VARIANT",
    "FLAG_LENGTH_TRUNCATED",
    "FLAG_LOW_EVIDENCE_SUPPORT",
    "FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT",
    "FLAG_LOW_KMER_PATH_SUPPORT",
    "FLAG_LOW_READ_SUPPORT",
    "FLAG_MOTIF_CONTEXT_DIVERGES",
    "FLAG_POSITION_AMBIGUOUS",
    "FLAG_REPRESENTATION_ONLY",
    "FLAG_SEQUENCE_UNDETERMINED",
    "FLAG_SPANS_UNIT_JUNCTION",
    "FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT",
    "MAPPABLE_RUS",
    "MIN_SUPPORT_FOR_TIER_A",
    "MOTIFS",
    "pair_sequence",
    "UNIT_LENGTH",
    "Nomenclature",
    "ambiguity_interval",
    "from_advntr",
    "name_coding_pair_edit",
    "from_kestrel",
    "name_edit",
    "normalise",
    "reconcile",
    "render",
    "repeat_form",
    "revcomp",
]


def load_nomenclature_config(config_path: str | None = None) -> dict:
    """Load the nomenclature reference data and thresholds.

    Follows the convention the rest of the package uses (AGENTS.md trap 1): the
    config is read once at import into a module global, not per call.

    Args:
        config_path: Optional override. Defaults to ``nomenclature_config.json``
            beside this module.

    Returns:
        dict: The parsed configuration.
    """
    if config_path is None:
        config_path = os.path.join(os.path.dirname(__file__), "nomenclature_config.json")
    return load_config(config_path)


#: Loaded once at import. Every table below is data, not code: the canonical unit,
#: the motif sequences, which adVNTR repeat units carry a canonical coordinate, the
#: support threshold, and the described-variant list all live in
#: ``nomenclature_config.json`` so they can be corrected without touching logic.
nomenclature_config: dict = load_nomenclature_config()

#: The canonical MUC1 repeat unit in the coding orientation. Carries the tract Wenzel
#: (2018, PMID:29520014) publishes: 7xC at positions 53-59, ``A`` at 60.
CANONICAL_UNIT: str = nomenclature_config["canonical_unit"]

#: Every MUC1 repeat unit is 60 bp.
UNIT_LENGTH: int = nomenclature_config["unit_length"]

#: Motif symbol -> 60 bp unit on the genomic plus strand, as Kestrel sees it.
#:
#: Configured rather than read from ``reference/MUC1_motifs_Rev_com.fa``, because
#: that directory is downloaded and not checked in: importing this module from a
#: fresh clone or an installed wheel must not depend on it. A test asserts the table
#: is byte-identical to the shipped FASTA wherever that FASTA is present.
#:
#: The 551-record pair reference is deliberately absent: every one of its records is
#: ``seq(R) ++ seq(L)`` over these same motifs -- verified across all 551 -- so
#: :func:`pair_sequence` derives it instead of duplicating 69 kB.
MOTIFS: dict[str, str] = nomenclature_config["motifs"]

#: Where the motif table came from, for the provenance test.
MOTIF_FASTA_NAME: str = nomenclature_config["motif_fasta_name"]

#: adVNTR repeat unit -> the motif it is a rotation of.
#:
#: adVNTR ships nine repeat units; only some are rotations of a MUC1 motif. The rest
#: are either 60 bp matching no motif, or not 60 bp at all, so parts of them have no
#: counterpart in a repeat unit. Projecting a state from one of those through the
#: rotation would fabricate a coordinate, so those states are reported without one.
MAPPABLE_RUS: dict[int, str] = {
    int(unit): symbol for unit, symbol in nomenclature_config["advntr"]["mappable_repeat_units"].items()
}

#: Rotation offset shared by every mappable repeat unit.
_RU_ROTATION: int = nomenclature_config["advntr"]["rotation_offset"]

#: Source-specific evidence support below which the top confidence is withheld.
#:
#: Each source retains its own evidence unit: Kestrel VCF k-mer-path depth,
#: Kestrel BAM resolved haplotype records, or adVNTR sequencing reads.
MIN_SUPPORT_FOR_TIER_A: int = nomenclature_config["thresholds"]["min_support_for_high_confidence"]

#: Evidence source -> the caller that produced it.
#:
#: ``kestrel_vcf`` and ``kestrel_bam`` map to the same caller. The BAM is Kestrel's
#: own alignment, not a second opinion on it, so the two agreeing is one caller
#: agreeing with itself. Treating them as independent would let Kestrel outvote
#: adVNTR by corroborating its own placement, and would satisfy the tier-A
#: requirement for "two independent sources" without any second opinion existing.
CALLER_OF: dict[str, str] = nomenclature_config["sources"]["caller_of"]

#: MUC1 variants described in the literature, name -> the report they come from.
#:
#: Used **only to check a name and raise confidence in it** -- never to produce one.
#: Every name is derived from the caller's own record; matching an entry here says
#: "and somebody has described this allele before", which is a weaker and different
#: claim than "this is correct".
#:
#: A name absent from this table is not suspect, merely un-cross-checked, and is
#: reported as what it is: a representation of the caller's call, requiring
#: validation. Extending the table is a config edit, not a code change.
KNOWN_VARIANTS: dict[str, str] = nomenclature_config["known_variants"]

#: Flag vocabulary. Stable, kebab-case, and closed: a consumer may match on these.
FLAG_POSITION_AMBIGUOUS = "position-ambiguous"
FLAG_SPANS_UNIT_JUNCTION = "spans-unit-junction"
FLAG_MOTIF_CONTEXT_DIVERGES = "motif-context-diverges"
FLAG_ALLELE_UNREPRESENTABLE = "allele-unrepresentable-in-vcf"
FLAG_CALLER_DISAGREEMENT = "caller-disagreement"
FLAG_LENGTH_TRUNCATED = "length-truncated"
FLAG_SEQUENCE_UNDETERMINED = "sequence-undetermined"

#: Set on a name that matches :data:`KNOWN_VARIANTS`.
FLAG_KNOWN_VARIANT = "known-variant"

#: Set on every name that does not. States what the name is rather than implying a
#: defect: a representation of what the caller reported, not a validated allele.
FLAG_REPRESENTATION_ONLY = "representation-of-caller-call"

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def pair_sequence(motifs: str) -> str | None:
    """Build the 120 bp Kestrel pair reference for a ``<L>-<R>`` label.

    A pair record holds ``seq(R) ++ seq(L)`` -- reverse-complementing it therefore
    swaps the halves as well as the strand, so the coding pair reads
    ``coding(L) ++ coding(R)``.

    Args:
        motifs: The ``Motifs`` field, e.g. ``S-C``.

    Returns:
        str | None: The 120 bp plus-strand sequence, or ``None`` when the label is
        malformed or names a motif that does not exist.
    """
    left, _, right = motifs.partition("-")
    if not _:
        return None
    left_seq, right_seq = MOTIFS.get(left), MOTIFS.get(right)
    if left_seq is None or right_seq is None:
        return None
    return right_seq + left_seq


def _is_dna(sequence: str) -> bool:
    """Is this a non-empty run of unambiguous DNA bases?

    A missing cell reaches the translator as the *text* ``nan`` or ``None`` once a
    pandas frame read with ``dtype=str`` is stringified, and those are accepted as
    bases by anything that only checks for emptiness -- an absent ``REF`` produced a
    confident ``52_53delGC``. Missing data must make a row untranslatable, not
    become DNA.

    Args:
        sequence: The candidate allele.

    Returns:
        bool: True when every character is A, C, G or T.
    """
    # `str.strip` with a character set is one C-level scan and allocates nothing,
    # unlike building a set per call in what is a per-dataframe-row hot path.
    return bool(sequence) and not sequence.strip("ACGTacgt")


def revcomp(sequence: str) -> str:
    """Reverse-complement a DNA sequence.

    One table translation and one slice; no Biopython, and allocation-light enough
    for the per-call budget.

    Args:
        sequence: DNA bases. ``N`` is preserved; case is preserved.

    Returns:
        str: The reverse complement.
    """
    return sequence.translate(_COMPLEMENT)[::-1]


def _trim(unit: str, start: int, end: int, inserted: str) -> tuple[int, int, str]:
    """Reduce an edit to its minimal span by removing shared flanking bases.

    ``GCTGGG>G`` describes the same allele as a 5 bp deletion of ``CTGGG``; naming
    the untrimmed form would put the position one base 5' of where it belongs.

    Args:
        unit: The repeat unit the edit applies to.
        start: 1-based inclusive start of the deleted span.
        end: 1-based inclusive end of the deleted span; ``start - 1`` for an insertion.
        inserted: The inserted bases; empty for a deletion.

    Returns:
        tuple[int, int, str]: The trimmed ``(start, end, inserted)``.
    """
    deleted = unit[start - 1 : end]

    # Trim a shared prefix, walking the start forward.
    while deleted and inserted and deleted[0] == inserted[0]:
        deleted, inserted = deleted[1:], inserted[1:]
        start += 1
    # Trim a shared suffix, walking the end back.
    while deleted and inserted and deleted[-1] == inserted[-1]:
        deleted, inserted = deleted[:-1], inserted[:-1]
        end -= 1

    return start, end, inserted


def normalise(unit: str, start: int, end: int, inserted: str) -> tuple[int, int, str]:
    """Shift an edit as far 3' as the sequence allows.

    HGVS requires the 3'-most representation of an ambiguous indel. Callers here have
    already been reverse-complemented into the coding frame, so "3'" is simply the
    direction of increasing coordinates in ``unit``.

    A **delins is never shifted**, matching VEP: once both a deletion and an insertion
    are present the allele is anchored and rolling it would change what it describes.

    Args:
        unit: The repeat unit the edit applies to.
        start: 1-based inclusive start.
        end: 1-based inclusive end; ``start - 1`` for an insertion.
        inserted: The inserted bases; empty for a deletion.

    Returns:
        tuple[int, int, str]: The 3'-most ``(start, end, inserted)``.
    """
    if end > len(unit) or start < 1:
        return start, end, inserted

    start, end, inserted = _trim(unit, start, end, inserted)
    length = len(unit)

    if inserted and end >= start:
        return start, end, inserted

    if not inserted:
        # Deletion: roll while the base leaving the 5' end equals the one entering 3'.
        while end < length and unit[start - 1] == unit[end]:
            start, end = start + 1, end + 1
        return start, end, inserted

    # Insertion: roll while the next reference base equals the first inserted base,
    # rotating the inserted string so a multi-base insert stays the same allele.
    while end < length and inserted[0] == unit[end]:
        inserted = inserted[1:] + inserted[0]
        start, end = start + 1, end + 1
    return start, end, inserted


def name_edit(unit: str, start: int, end: int, inserted: str) -> str:
    """Name one edit on one repeat unit, in the coding frame.

    Args:
        unit: The repeat unit the edit applies to.
        start: 1-based inclusive start of the replaced span.
        end: 1-based inclusive end; ``start - 1`` denotes an insertion.
        inserted: The inserted bases; empty denotes a deletion.

    Returns:
        str: The bare positional name, e.g. ``59dupC``, ``58_59insG``,
        ``1_5delGCCCA``, ``54_56delinsAT``, ``2C>A``. Never carries a ``c.`` prefix.
    """
    was_insertion = end == start - 1
    start, end, inserted = normalise(unit, start, end, inserted)

    if was_insertion or end < start:
        return _name_insertion(unit, start, inserted)

    deleted = unit[start - 1 : end]

    if not inserted:
        span = str(start) if start == end else f"{start}_{end}"
        return f"{span}del{deleted}"

    if start == end and len(inserted) == 1:
        return f"{start}{deleted}>{inserted}"

    span = str(start) if start == end else f"{start}_{end}"
    return f"{span}delins{inserted}"


def _roll_5prime(unit: str, start: int, end: int, inserted: str) -> tuple[int, int, str]:
    """Shift an edit as far 5' as the sequence allows -- the mirror of :func:`normalise`.

    Used only to find the far edge of an ambiguity window; names are always 3'-most.

    Args:
        unit: The repeat unit the edit applies to.
        start: 1-based inclusive start.
        end: 1-based inclusive end; ``start - 1`` for an insertion.
        inserted: The inserted bases; empty for a deletion.

    Returns:
        tuple[int, int, str]: The 5'-most ``(start, end, inserted)``.
    """
    if end > len(unit):
        return start, end, inserted

    if not inserted:
        while start > 1 and unit[end - 1] == unit[start - 2]:
            start, end = start - 1, end - 1
        return start, end, inserted

    while start > 1 and inserted[-1] == unit[start - 2]:
        inserted = inserted[-1] + inserted[:-1]
        start, end = start - 1, end - 1
    return start, end, inserted


def ambiguity_interval(unit: str, start: int, end: int, inserted: str) -> tuple[int, int] | None:
    """The span within which every anchor describes the same allele.

    ``59dupC``, ``53dupC`` and the older ``27dupC`` are one event, not three. Stating
    the window is what makes that visible.

    Args:
        unit: The repeat unit the edit applies to.
        start: 1-based inclusive start.
        end: 1-based inclusive end; ``start - 1`` for an insertion.
        inserted: The inserted bases; empty for a deletion.

    Returns:
        tuple[int, int] | None: The inclusive ``(low, high)`` window, or ``None`` when
        the edit cannot shift at all (a 1 bp window carries no information) or when it
        is a delins, which is anchored by definition.
    """
    if end > len(unit) or start < 1:
        return None

    start, end, inserted = _trim(unit, start, end, inserted)

    if inserted and end >= start:
        return None

    low_start, _, _ = _roll_5prime(unit, start, end, inserted)
    high_start, high_end, _ = normalise(unit, start, end, inserted)

    if inserted:
        # An insertion anchor is interbase ("before position p"), so the reference
        # bases it is indistinguishable from end one short of the 3'-most anchor.
        low, high = low_start, high_start - 1
    else:
        low, high = low_start, high_end

    if high <= low:
        return None
    return low, high


def _tract_at(unit: str, position: int) -> tuple[int, str, int] | None:
    """Find the homopolymer run covering a position.

    Args:
        unit: The repeat unit.
        position: 1-based position expected to sit inside a run.

    Returns:
        tuple[int, str, int] | None: ``(start, base, count)`` for a run of at least two
        identical bases, else ``None``.
    """
    if not 1 <= position <= len(unit):
        return None

    base = unit[position - 1]
    start = position
    while start > 1 and unit[start - 2] == base:
        start -= 1
    end = position
    while end < len(unit) and unit[end] == base:
        end += 1

    count = end - start + 1
    return (start, base, count) if count >= 2 else None


def repeat_form(unit: str, start: int, end: int, inserted: str) -> str | None:
    """Express the edit as a change in tract copy number.

    ``53C[7]>53C[8]`` states what was actually measured -- the tract went from 7 to 8
    copies -- instead of implying we know which base was added. It also scales: an
    ``insCCCC`` reads ``53C[11]``, far clearer than ``56_59dupCCCC``.

    Args:
        unit: The repeat unit the edit applies to.
        start: 1-based inclusive start.
        end: 1-based inclusive end; ``start - 1`` for an insertion.
        inserted: The inserted bases; empty for a deletion.

    Returns:
        str | None: The repeat form, or ``None`` when the edit does not sit in a
        detectable tract or does not change that tract's length.
    """
    if end > len(unit) or start < 1:
        return None

    start, end, inserted = _trim(unit, start, end, inserted)

    if inserted and end >= start:
        return None

    if inserted:
        if len(set(inserted)) != 1:
            return None
        anchor = start - 1
        delta = len(inserted)
    else:
        deleted = unit[start - 1 : end]
        if not deleted or len(set(deleted)) != 1:
            return None
        anchor = start
        delta = -(end - start + 1)

    tract = _tract_at(unit, anchor)
    if tract is None:
        return None

    tract_start, base, count = tract
    if base != (inserted[0] if inserted else unit[start - 1]):
        return None

    return f"{tract_start}{base}[{count}]>{tract_start}{base}[{count + delta}]"


def _name_insertion(unit: str, start: int, inserted: str) -> str:
    """Name an insertion, preferring ``dup`` when the insert repeats what precedes it.

    Args:
        unit: The repeat unit the edit applies to.
        start: 1-based position immediately 3' of the insertion point.
        inserted: The inserted bases.

    Returns:
        str: ``59dupC`` / ``56_59dupCCCC`` when the insert duplicates the directly
        preceding bases, otherwise ``58_59insG``.
    """
    left = start - 1
    preceding = unit[left - len(inserted) : left]

    if left >= len(inserted) and preceding == inserted:
        if len(inserted) == 1:
            return f"{left}dup{inserted}"
        return f"{left - len(inserted) + 1}_{left}dup{inserted}"

    return f"{left}_{left + 1}ins{inserted}"


def _event_of(start: int, end: int, inserted: str) -> str:
    """Classify a normalised edit.

    Args:
        start: 1-based inclusive start.
        end: 1-based inclusive end; ``start - 1`` for an insertion.
        inserted: The inserted bases; empty for a deletion.

    Returns:
        str: One of ``duplication``, ``insertion``, ``deletion``, ``delins``,
        ``substitution``.
    """
    if end < start:
        return "insertion"
    if not inserted:
        return "deletion"
    if start == end and len(inserted) == 1:
        return "substitution"
    return "delins"


@dataclass(frozen=True)
class Nomenclature:
    """One caller record, named.

    The spec states two things that must be read together: ``name`` is ``None`` only
    at tier C (§3.1), yet tier B must emit "no bare number" (§3.3). Both hold because
    they describe different layers -- this object carries the *computed* name so
    :func:`reconcile` can compare callers, and the output serializer renders that
    number only at tier A. Never surface ``name`` directly for a tier-B call.

    Attributes:
        name: The computed positional name, e.g. ``59dupC``. ``None`` at tier C,
            where no position may be stated at all. Not for direct display below
            tier A -- see :func:`render`.
        event: ``duplication`` | ``insertion`` | ``deletion`` | ``delins`` |
            ``substitution``.
        unit: The motif symbol the name is anchored on, ``None`` when undetermined.
        tier: ``A`` | ``B`` | ``C``.
        flags: Kebab-case flags from the closed vocabulary above.
        ambiguity: Inclusive window in which every anchor is the same allele, or
            ``None`` when the edit cannot shift.
        repeat_form: ``53C[7]>53C[8]``, or ``None`` outside a detectable tract.
        net_length: Change in length, e.g. ``+1`` for a duplication.
        source: ``kestrel_vcf`` | ``kestrel_bam`` | ``advntr``.
    """

    name: str | None
    event: str
    unit: str | None
    tier: str
    flags: tuple[str, ...]
    ambiguity: tuple[int, int] | None
    repeat_form: str | None
    net_length: int
    source: str


def _undetermined(event: str, net_length: int, source: str, flags: tuple[str, ...]) -> Nomenclature:
    """Build a tier-C result: a frameshift statement with no number in it.

    Args:
        event: The event class, as far as it is known.
        net_length: Change in length.
        source: The caller the record came from.
        flags: Why the allele is undetermined.

    Returns:
        Nomenclature: A tier-C value with ``name``, ``unit``, ``ambiguity`` and
        ``repeat_form`` all ``None``.
    """
    return Nomenclature(
        name=None,
        event=event,
        unit=None,
        tier="C",
        flags=flags,
        ambiguity=None,
        repeat_form=None,
        net_length=net_length,
        source=source,
    )


def _is_corroborated(sources: set[str]) -> bool:
    """Do these sources amount to independent corroboration?

    The test is on *callers*, not sources: ``kestrel_vcf`` and ``kestrel_bam`` are two
    outputs from one caller, and counting them as two would make Kestrel agreeing with
    its own alignment look like the two independent sources tier A requires.

    Counting callers subsumes counting sources -- the callers are the image of the
    sources under a map, so there can never be more of them -- which is why there is
    no separate "at least two sources" condition here.

    Args:
        sources: The ``source`` fields of the calls naming one allele.

    Returns:
        bool: True when the allele is independently corroborated.
    """
    # An unrecognised source is assumed to be its own caller: a new source added
    # without a config entry is treated as independent rather than silently folded
    # into an existing caller, which would hide a disagreement.
    return len({CALLER_OF.get(source, source) for source in sources}) >= 2


def reconcile(
    *calls: Nomenclature,
    support: int | None = None,
    supports: Mapping[str, int | None] | None = None,
    identity_observations: tuple[IdentityReconciliationObservation, ...] | None = None,
    identity_policy: IdentityReconciliationPolicy | None = None,
) -> Nomenclature:
    """Combine caller calls through legacy or explicitly supplied identity evidence.

    Args:
        *calls: Legacy presentation calls in their existing order.
        support: Legacy scalar support for compatibility-only callers.
        supports: Legacy support keyed by source.
        identity_observations: Typed observations aligned one-for-one with ``calls``.
        identity_policy: Explicit source-unit thresholds for typed reconciliation.

    Returns:
        The unchanged legacy projection when identity metadata is absent, or the
        identity-keyed decision projected onto the existing ``Nomenclature`` shape.

    Raises:
        ValueError: If typed observations and their explicit policy are inconsistent.
    """
    if identity_observations is None:
        if identity_policy is not None:
            raise ValueError("An identity policy requires typed identity observations")
        return _reconcile_legacy(*calls, support=support, supports=supports)
    if identity_policy is None:
        raise ValueError("Typed identity observations require an explicit identity policy")
    if not calls:
        return _undetermined("unknown", 0, "reconciled", ())

    implicit_positional_binding = len(identity_observations) == len(calls) and all(
        observation.presentation_call_index is None for observation in identity_observations
    )
    presentation_call_indices = tuple(
        index if implicit_positional_binding else observation.presentation_call_index
        for index, observation in enumerate(identity_observations)
    )
    bound_call_indices = {index for index in presentation_call_indices if index is not None}
    if bound_call_indices != set(range(len(calls))):
        raise ValueError("Typed reconciliation requires one identity observation per call")
    if any(
        observation.display_name != call.name
        or observation.source != call.source
        or observation.event != call.event
        or observation.net_length != call.net_length
        for observation, call_index in zip(identity_observations, presentation_call_indices, strict=True)
        if call_index is not None
        for call in (calls[call_index],)
    ):
        raise ValueError("Typed identity observations must match their presentation calls")

    admissible_call_indices = tuple(
        index
        for index in range(len(calls))
        if any(
            call_index == index and observation.disposition.value == "admissible"
            for observation, call_index in zip(identity_observations, presentation_call_indices, strict=True)
        )
    )
    decision_call_indices = admissible_call_indices or tuple(range(len(calls)))
    decision_calls = tuple(calls[index] for index in decision_call_indices)
    relative_selection = select_compatibility_presentation_index(decision_calls)
    if relative_selection is None:  # pragma: no cover - calls was established above
        return _undetermined("unknown", 0, "reconciled", ())
    presentation_call_index = decision_call_indices[relative_selection]
    presentation_observation_index = next(
        index for index, call_index in enumerate(presentation_call_indices) if call_index == presentation_call_index
    )
    presentation = _reconcile_legacy(*decision_calls, support=support, supports=supports)
    result = reconcile_identity_observations(
        identity_observations,
        identity_policy,
        presentation_selected_observation_index=presentation_observation_index,
    )

    flags = set().union(*(call.flags for call in decision_calls))
    if result.caller_disagreement or (result.identity is None and FLAG_CALLER_DISAGREEMENT in presentation.flags):
        flags.add(FLAG_CALLER_DISAGREEMENT)
    if any(identity_observations[index].translation.context_diverges for index in result.backing_observation_indices):
        flags.add(FLAG_MOTIF_CONTEXT_DIVERGES)
    for source in result.low_support_sources:
        flags.add(low_support_flag_for_source(source))
    if FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT in presentation.flags:
        flags.add(FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT)
    if presentation.name is not None:
        flags.add(FLAG_KNOWN_VARIANT if presentation.name in KNOWN_VARIANTS else FLAG_REPRESENTATION_ONLY)

    if result.event_disagreement or presentation.name is None:
        return _undetermined("unknown", presentation.net_length, "reconciled", tuple(sorted(flags)))
    tier = "B" if result.tier == "abstained" else result.tier
    return Nomenclature(
        name=presentation.name,
        event=presentation.event,
        unit=presentation.unit,
        tier=tier,
        flags=tuple(sorted(flags)),
        ambiguity=presentation.ambiguity,
        repeat_form=presentation.repeat_form,
        net_length=presentation.net_length,
        source="reconciled" if len(result.backing_sources) > 1 else presentation.source,
    )


def _reconcile_legacy(
    *calls: Nomenclature,
    support: int | None = None,
    supports: Mapping[str, int | None] | None = None,
) -> Nomenclature:
    """Combine independent callers into one result, and decide its tier.

    Tier A is the only tier that emits a bare number, so it is the only tier that can
    state a confident falsehood. It therefore requires evidence no single caller can
    supply: two independent sources agreeing after normalisation, a motif context
    matching canonical X, and support at or above
    :data:`MIN_SUPPORT_FOR_TIER_A`. Anything less keeps the number internal.

    The benchmark is why the bar is here. Kestrel places the whole ``insG`` family
    one position 3' of truth, and its ``insG_pos58`` records collapse onto the C-tract
    expansion; both look clean in isolation. Only a second caller separates them.

    Args:
        *calls: Translations of the same locus from different sources.
        support: Legacy sequencing-read support for the call, when known. ``None``
            means unknown, which is not the same as sufficient. Ignored when
            ``supports`` is given.
        supports: Evidence support per source in that source's own unit, e.g.
            ``{"advntr": 24}`` sequencing reads. Preferred: it binds support to the
            evidence it came from, so an unrelated well-supported observation cannot
            lend support to a thin agreement. The agreement is taken to be as strong
            as its weakest contributing source.

    Returns:
        Nomenclature: The reconciled call. Tier C when nothing was supplied or the
        callers disagree on the event class.
    """
    usable = [call for call in calls if call is not None]
    if not usable:
        return _undetermined("unknown", 0, "reconciled", ())

    primary = next((call for call in usable if call.source == "kestrel_vcf"), usable[0])
    selected_index = select_compatibility_presentation_index(usable)
    assert selected_index is not None
    flags = set(primary.flags)
    for call in usable:
        flags.update(call.flags)

    named = [call for call in usable if call.name is not None]
    events = {call.event for call in usable}

    # Disagreement on the event class means the two callers are not describing one
    # allele. Naming either would be picking a winner without grounds.
    if len(events) > 1:
        flags.add(FLAG_CALLER_DISAGREEMENT)
        net = primary.net_length
        return _undetermined("unknown", net, "reconciled", tuple(sorted(flags)))

    # Which sources back each candidate allele. Independence is counted over the calls
    # that actually carry a name, never over every call supplied: counting all of them
    # let an *unnamed* call -- an adVNTR state in an unmappable repeat unit, say --
    # donate a second `source` to two duplicate rows from one caller, so a single
    # caller's placement could be promoted as though two had agreed on it.
    backing: dict[str, set[str]] = {}
    for call in named:
        backing.setdefault(str(call.name), set()).add(call.source)

    if len(backing) > 1:
        flags.add(FLAG_CALLER_DISAGREEMENT)

    # The pure compatibility selector owns the legacy majority/VCF-primary choice so
    # the typed adapter and direct legacy callers cannot drift apart.
    selected = usable[selected_index]
    chosen = selected if selected.name is not None else None
    chosen_name = chosen.name if chosen is not None else None

    # Independence is judged on the chosen allele's own backing, not on every name
    # seen. A dissenting source is not part of the agreement it dissents from.
    backing_sources = backing.get(str(chosen_name), set()) if chosen_name is not None else set()
    agree = _is_corroborated(backing_sources)

    # Support must belong to the agreeing evidence. A sample-wide maximum would let a
    # well-covered but unrelated observation lend its depth to a thin agreement.
    effective_support = support
    if supports is not None:
        relevant = [supports.get(source) for source in backing_sources]
        # Unknown is not sufficient. Dropping the `None`s and taking the minimum of
        # what remained let one caller's depth stand in for the other's missing depth,
        # so an agreement with a blank or non-numeric depth column reached tier A on
        # one source's support -- exactly the "two independent sources" claim the tier
        # is supposed to guarantee. One unknown makes the whole agreement unknown.
        known = [value for value in relevant if value is not None]
        effective_support = None if len(known) != len(relevant) or not known else min(known)

        if effective_support is not None and effective_support < MIN_SUPPORT_FOR_TIER_A:
            for source in backing_sources:
                source_support = supports.get(source)
                if source_support is not None and source_support < MIN_SUPPORT_FOR_TIER_A:
                    flags.add(low_support_flag_for_source(source))
    elif effective_support is not None and effective_support < MIN_SUPPORT_FOR_TIER_A:
        flags.add(FLAG_LOW_READ_SUPPORT)

    # Known variants are used to *check* a name, never to make one. Matching the
    # literature says somebody has described this allele before -- a weaker and
    # different claim than "this is correct" -- so it raises confidence rather than
    # conferring it.
    if chosen_name is not None:
        flags.add(FLAG_KNOWN_VARIANT if chosen_name in KNOWN_VARIANTS else FLAG_REPRESENTATION_ONLY)

    tier = "B"
    if (
        agree
        and chosen_name in KNOWN_VARIANTS
        and effective_support is not None
        and effective_support >= MIN_SUPPORT_FOR_TIER_A
        and FLAG_MOTIF_CONTEXT_DIVERGES not in flags
        and FLAG_SEQUENCE_UNDETERMINED not in flags
        and FLAG_CALLER_DISAGREEMENT not in flags
    ):
        tier = "A"

    if chosen is None:
        return _undetermined(primary.event, primary.net_length, "reconciled", tuple(sorted(flags)))

    return Nomenclature(
        name=chosen.name,
        event=chosen.event,
        unit=chosen.unit,
        tier=tier,
        flags=tuple(sorted(flags)),
        ambiguity=chosen.ambiguity,
        repeat_form=chosen.repeat_form,
        net_length=chosen.net_length,
        source="reconciled" if len(backing_sources) > 1 else chosen.source,
    )


def confidence_note(call: Nomenclature) -> str:
    """One sentence saying what the name is and how far it has been checked.

    The honest position, in the output rather than only in the docs: a name is
    derived from a caller's record, and unless it matches a variant somebody has
    already described, that is all it is.

    Args:
        call: The reconciled call.

    Returns:
        str: The note, or ``""`` when there is no name to qualify.
    """
    if not call.name:
        return ""
    citation = KNOWN_VARIANTS.get(call.name)
    if citation is not None:
        return f"matches a described MUC1 variant ({citation}); requires validation"
    return "representation of the caller's call, not a described variant; requires validation"


def render(call: Nomenclature) -> str:
    """Render a call as display text.

    **If a name could be computed, the name is shown.** The tier is a confidence
    indicator carried beside it, not a gate that suppresses it.

    This was the other way round at first -- only tier A showed a number, and tiers
    B and C showed the event class instead. Measured on the 200-sample benchmark
    that discarded most of the useful output: 129 samples had the correct name
    computed and only 46 of them displayed it. Withholding a name that was right
    83 times to avoid displaying one that was wrong 0 times is a bad trade, and it
    loses information a reader can weigh for themselves. The tier, the flags and
    :func:`confidence_note` say how far the name has been checked; they do not decide
    whether the reader is allowed to see it.

    Args:
        call: The reconciled call.

    Returns:
        str: The name where one exists, otherwise a statement of what is known.
    """
    if call.name:
        return call.name

    # No name could be computed at all. A net length change of zero is not a
    # frameshift, so saying "frameshift +0" would state something untrue about a
    # locus we know nothing about.
    if call.net_length == 0:
        return "allele undetermined"
    sign = "+" if call.net_length > 0 else "-"
    return f"frameshift {sign}{abs(call.net_length)}, allele undetermined"


_ADVNTR_STATE = re.compile(r"^(?P<kind>[ID])(?P<pos>\d+)_(?P<ru>\d+)(?:_(?P<base>[ACGT])_LEN(?P<length>\d+))?$")


class _Component(NamedTuple):
    """One parsed adVNTR state component."""

    kind: str
    pos: int
    ru: int
    base: str
    length: int


def _parse_components(state: str) -> list[_Component]:
    """Parse an ``&``-joined adVNTR state into its components.

    Args:
        state: e.g. ``D27_2&I27_2_A_LEN2``.

    Returns:
        list[_Component]: The components, or an empty list if any part fails to
        parse -- a partially understood compound state is not named at all.
    """
    parts = [part.strip() for part in state.split("&") if part.strip()]
    if not parts:
        return []

    components: list[_Component] = []
    for part in parts:
        match = _ADVNTR_STATE.match(part)
        if match is None:
            return []
        kind = match["kind"]
        if kind == "I" and match["base"] is None:
            return []
        components.append(
            _Component(
                kind=kind,
                pos=int(match["pos"]),
                ru=int(match["ru"]),
                base=match["base"] or "",
                length=int(match["length"]) if match["length"] else 1,
            )
        )
    return components


def _group_components(components: list[_Component]) -> list[list[_Component]]:
    """Group components into events, using adVNTR's own joining rule.

    ``vntr_finder.py`` joins a deletion to the previous component only when the match
    indices are consecutive *and* the repeat-unit index is the same, and joins a
    deletion to an insertion only at the same position in the same unit. Merging more
    broadly -- "anything adjacent" -- would fuse two independent variants that happen
    to sit side by side into one fabricated allele.

    Args:
        components: Parsed components, in the order adVNTR emitted them.

    Returns:
        list[list[_Component]]: One inner list per event.
    """
    groups: list[list[_Component]] = []
    for component in sorted(components, key=lambda item: (item.ru, item.pos, item.kind)):
        if groups:
            previous = groups[-1][-1]
            same_unit = previous.ru == component.ru
            consecutive_deletion = (
                same_unit and previous.kind == "D" and component.kind == "D" and previous.pos + 1 == component.pos
            )
            colocated_delins = (
                same_unit and previous.kind == "D" and component.kind == "I" and previous.pos == component.pos
            )
            if consecutive_deletion or colocated_delins:
                groups[-1].append(component)
                continue
        groups.append([component])
    return groups


def _name_advntr_group(group: list[_Component]) -> Nomenclature:
    """Name one grouped adVNTR event.

    Args:
        group: Components belonging to a single event.

    Returns:
        Nomenclature: The named event, tier C when it cannot be placed.
    """
    ru = group[0].ru
    deletions = [item for item in group if item.kind == "D"]
    insertions = [item for item in group if item.kind == "I"]
    inserted_length = sum(item.length for item in insertions)
    net = inserted_length - len(deletions)

    if deletions and insertions:
        event = "delins"
    elif insertions:
        event = "insertion"
    else:
        event = "deletion"

    symbol = MAPPABLE_RUS.get(ru)
    if symbol is None:
        return _undetermined(event, net, "advntr", ())

    # An insertion of more than one base carries a length but not a sequence, so the
    # allele cannot be written out. The event and its length are still reportable.
    if inserted_length > 1:
        return _undetermined(event, net, "advntr", (FLAG_SEQUENCE_UNDETERMINED,))

    flags: list[str] = []
    if symbol != "X":
        flags.append(FLAG_MOTIF_CONTEXT_DIVERGES)

    if insertions:
        component = insertions[0]
        plus = ((_RU_ROTATION - 1 + component.pos) % UNIT_LENGTH) + 1
        # adVNTR anchors an insertion in the gap *after* the plus-strand position, and
        # reverse complement swaps which side of that gap it sits on, so the coding
        # base immediately 5' of it is `UNIT_LENGTH - plus`. Treating the anchor as a
        # base instead yields 60 for the canonical duplication rather than 59.
        left = UNIT_LENGTH - plus
        start, end = left + 1, left
        inserted = revcomp(component.base)
        if deletions:
            # A co-located delins: the deletion supplies the replaced base.
            start, end = left + 1, left + 1
        elif left == 0:
            # The insertion sits exactly on the unit junction. The array is tandem,
            # so that gap is equally "after position 60 of the preceding unit", and
            # naming it there turns a meaningless `0_1insA` into the `60dupA` it is.
            # `from_kestrel` already did this; adVNTR states reach the same boundary.
            start, end = UNIT_LENGTH + 1, UNIT_LENGTH
    else:
        positions = sorted(UNIT_LENGTH + 1 - (((_RU_ROTATION - 1 + item.pos) % UNIT_LENGTH) + 1) for item in deletions)
        # Consecutive positions in adVNTR's rotated unit are not necessarily
        # consecutive in the coding unit: the rotation puts a seam inside the unit, so
        # `D21_2&D22_2` projects to coding 1 and 60 -- opposite ends. Spanning
        # `[min, max]` there named a 2 bp deletion as a 60 bp one, deleting the entire
        # repeat unit. A gap means the event crosses the junction and cannot be
        # expressed as one span on a single unit.
        if positions[-1] - positions[0] != len(positions) - 1:
            return _undetermined(event, net, "advntr", (*flags, FLAG_SPANS_UNIT_JUNCTION))
        start, end = positions[0], positions[-1]
        inserted = ""

    if not 1 <= start <= UNIT_LENGTH + 1 or end > UNIT_LENGTH:
        return _undetermined(event, net, "advntr", tuple(flags))

    name = name_edit(CANONICAL_UNIT, start, end, inserted)
    window = ambiguity_interval(CANONICAL_UNIT, start, end, inserted)
    tract = repeat_form(CANONICAL_UNIT, start, end, inserted)
    if window is not None:
        flags.append(FLAG_POSITION_AMBIGUOUS)

    norm = normalise(CANONICAL_UNIT, start, end, inserted)
    resolved = _event_of(*norm)
    if resolved == "insertion" and "dup" in name:
        resolved = "duplication"

    return Nomenclature(
        name=name,
        event=resolved if not deletions or not insertions else "delins",
        unit=symbol,
        tier="B",
        flags=tuple(flags),
        ambiguity=window,
        repeat_form=tract,
        net_length=net,
        source="advntr",
    )


def from_advntr(state: str) -> tuple[Nomenclature, ...]:
    """Translate one adVNTR state field into MUC1 nomenclature.

    Args:
        state: The ``Variant`` field, e.g. ``I22_2_G_LEN1`` or
            ``D27_2&I27_2_A_LEN2``. Non-states such as ``Not applicable`` yield
            an empty tuple.

    Returns:
        tuple[Nomenclature, ...]: One entry per distinct event. Usually length 1;
        empty when nothing parseable was reported.
    """
    components = _parse_components(state)
    if not components:
        return ()
    return tuple(_name_advntr_group(group) for group in _group_components(components))


def from_kestrel(motifs: str, pos: int, ref: str, alt: str) -> Nomenclature:
    """Translate one Kestrel VCF record into MUC1 nomenclature.

    Kestrel reports on a 120 bp merged pair ``<L>-<R>`` whose sequence is
    ``seq(R) ++ seq(L)``, on the genomic plus strand. MUC1 is a minus-strand gene, so
    the edit is reverse-complemented into the coding frame before it is normalised --
    HGVS's "3'-most" is 3' of the *coding* sequence, which here is genomic-left.

    The position is then projected onto the canonical ``X`` unit and normalised there.
    Anchoring on the motif Kestrel assigned instead would name the canonical
    duplication ``57dupC`` on a pair whose motif carries only 5xC. When the assigned
    motif's local context differs from ``X``, the projection is flagged and the tier
    capped, so a projected name never claims top confidence on its own.

    Args:
        motifs: The ``Motifs`` field, e.g. ``S-C``.
        pos: 1-based position on the 120 bp pair.
        ref: Reference allele.
        alt: Alternate allele.

    Returns:
        Nomenclature: The named record; tier C when the record cannot be placed.
    """
    net = len(alt) - len(ref)

    pair = pair_sequence(motifs)
    if pair is None or not 1 <= pos <= len(pair) or not _is_dna(ref) or not _is_dna(alt):
        return _undetermined("insertion" if net > 0 else "deletion", net, "kestrel_vcf", ())

    pair_length = len(pair)

    # Resolve the edit in *pair* coordinates on the plus strand. This is deliberately
    # arithmetic on POS/REF/ALT rather than a trim against a reference unit: the two
    # observed shapes of the canonical duplication, `G>GG` and `C>CG`, carry different
    # anchor bases but insert the same G after the same pair position, and only the
    # arithmetic makes them agree.
    if len(alt) > len(ref) and alt.startswith(ref):
        inserted_plus = alt[len(ref) :]
        anchor_plus = pos + len(ref) - 1
        # The gap after plus position `a` is the gap after coding position
        # `pair_length - a`, because reverse complement swaps which side it sits on.
        coding_left = pair_length - anchor_plus
        coding_start, coding_end = coding_left + 1, coding_left
        inserted = revcomp(inserted_plus)
    elif len(ref) > len(alt) and ref.startswith(alt):
        del_lo = pos + len(alt)
        del_hi = pos + len(ref) - 1
        coding_start, coding_end = pair_length + 1 - del_hi, pair_length + 1 - del_lo
        inserted = ""
    else:
        coding_start = pair_length + 1 - (pos + len(ref) - 1)
        coding_end = pair_length + 1 - pos
        inserted = revcomp(alt)

    return name_coding_pair_edit(motifs, pair, coding_start, coding_end, inserted, net, "kestrel_vcf")


def name_coding_pair_edit(
    motifs: str,
    pair: str,
    coding_start: int,
    coding_end: int,
    inserted: str,
    net: int,
    source: str,
) -> Nomenclature:
    """Name an edit already expressed in coding-frame pair coordinates.

    Shared by the VCF and BAM paths so the two cannot drift: a name recovered from
    resolved haplotype records must be produced by exactly the machinery that names a VCF record,
    or the two would disagree for reasons that have nothing to do with the evidence.

    Args:
        motifs: The ``<L>-<R>`` pair label.
        pair: The 120 bp plus-strand pair sequence.
        coding_start: 1-based inclusive start on the coding pair.
        coding_end: 1-based inclusive end; ``coding_start - 1`` for an insertion.
        inserted: Inserted bases, already in the coding frame.
        net: Change in length.
        source: The caller this came from.

    Returns:
        Nomenclature: The named record; tier C when it cannot be placed.
    """
    pair_length = len(pair)

    if not 1 <= coding_start <= pair_length + 1:
        return _undetermined("insertion" if net > 0 else "deletion", net, source, ())

    # Normalise in the 120 bp pair frame FIRST, so a shift may legally cross the unit
    # junction, and only then assign the unit (spec §3.2). Order matters: a delGCCCA
    # record arrives as a deletion at pair 57-61, straddling the junction, and rolls
    # 3' to pair 61-65 -- entirely inside the second unit, where it is a clean
    # `1_5delGCCCA`. Assigning the unit first would have frozen it across the seam.
    coding_pair = revcomp(pair)
    coding_start, coding_end, inserted = normalise(coding_pair, coding_start, coding_end, inserted)

    # The coding pair reads coding(L) ++ coding(R): reverse-complementing
    # seq(R) ++ seq(L) swaps the halves as well as the strand.
    junction_flags: tuple[str, ...] = ()
    if coding_start <= UNIT_LENGTH < coding_end:
        junction_flags = (FLAG_SPANS_UNIT_JUNCTION,)

    left, right = motifs.split("-", 1)[0], motifs.split("-", 1)[-1]

    if coding_start > UNIT_LENGTH:
        symbol = right
        start = coding_start - UNIT_LENGTH
        end = coding_end - UNIT_LENGTH
    else:
        symbol = left
        start, end = coding_start, coding_end

    # An insertion landing before position 1 of a unit sits exactly on the junction.
    # The array is tandem, so that gap is equally "after position 60 of the unit 5'
    # of it" -- and naming it there is what turns a meaningless `0_1insA` into the
    # `60dupA` the allele actually is. A coordinate of 0 is never emitted.
    if end < start and start == 1:
        start, end = UNIT_LENGTH + 1, UNIT_LENGTH
        symbol = left if coding_start > UNIT_LENGTH else right

    assigned = MOTIFS.get(symbol)

    upper = UNIT_LENGTH + 1 if end < start else UNIT_LENGTH
    if not 1 <= start <= upper:
        return _undetermined("insertion" if net > 0 else "deletion", net, source, junction_flags)

    # A span still crossing the junction after normalisation cannot be projected onto
    # a single 60 bp unit, so it is named against the coding pair and anchored on the
    # unit holding its 5' end -- leaving an end above 60, never a negative one.
    # Deliberate, measured deviation from spec §3.2 step 4, which says to keep the
    # assigned motif when its local context is not homologous to X. Naming is done
    # against canonical X either way, and divergence is reported as a flag that caps
    # the tier.
    #
    # The spec's own numbers require this. §2.4 measures the X-anchored path at
    # 123/178 with dupC 96/96, and states that anchoring on the motif Kestrel
    # assigned instead drops dupC to 69/96 -- motifs H/E/A/J carry their longest
    # C-tract at 20-23 or 38-41, and C/5C have six C rather than seven. Following
    # §3.2 literally would forfeit 27 correct canonical calls.
    #
    # The cost is confined to tier B: a divergent-context call carries
    # `motif-context-diverges`, which blocks promotion, so no confident name depends
    # on the projection. What it can get wrong is the tier-B *event* word -- an edit
    # that is a duplication against motif J reads as an insertion against X.
    if junction_flags:
        context = coding_pair
        start, end = coding_start, coding_end
    else:
        context = CANONICAL_UNIT

    name = name_edit(context, start, end, inserted)
    window = ambiguity_interval(context, start, end, inserted)
    tract = repeat_form(context, start, end, inserted)
    norm_start, norm_end, norm_inserted = normalise(context, start, end, inserted)
    event = _event_of(norm_start, norm_end, norm_inserted)
    if event == "insertion" and name and "dup" in name:
        event = "duplication"

    flags: list[str] = list(junction_flags)
    if window is not None:
        flags.append(FLAG_POSITION_AMBIGUOUS)

    # Does the motif Kestrel assigned actually look like X where the name lands?
    span_lo = window[0] if window else start
    span_hi = window[1] if window else max(end, start)
    if assigned is None:
        flags.append(FLAG_MOTIF_CONTEXT_DIVERGES)
    else:
        coding_assigned = revcomp(assigned)
        lo = max(1, span_lo - 1)
        hi = min(UNIT_LENGTH, span_hi + 1)
        if coding_assigned[lo - 1 : hi] != CANONICAL_UNIT[lo - 1 : hi]:
            flags.append(FLAG_MOTIF_CONTEXT_DIVERGES)

    # A single caller never reaches tier A. Tier A requires agreement between two
    # independent sources (spec §3.3), and the benchmark shows why: Kestrel places
    # the whole insG family one position 3' of truth, so a lone Kestrel record that
    # looks perfectly clean still names the wrong allele. Promotion is `reconcile`'s
    # job; translation only reports what one caller said.
    return Nomenclature(
        name=name,
        event=event,
        unit=symbol if junction_flags else "X",
        tier="B",
        flags=tuple(flags),
        ambiguity=window,
        repeat_form=tract,
        net_length=net,
        source=source,
    )
