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

This module is pure: no I/O, no pandas, no logging. It is safe to call from a
dataframe ``.apply`` and from tests without fixtures.

Research use only.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:  # pragma: no cover - typing only
    from collections.abc import Mapping

__all__ = [
    "CANONICAL_UNIT",
    "MAPPABLE_RUS",
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

#: The canonical MUC1 repeat unit in the coding orientation. Carries the tract Wenzel
#: (2018, PMID:29520014) publishes: 7xC at positions 53-59, ``A`` at 60. Verified
#: byte-identical to the MUC1 simulator's ``X`` unit.
CANONICAL_UNIT = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"

#: Every MUC1 repeat unit is 60 bp.
UNIT_LENGTH = 60

#: Flag vocabulary. Stable, kebab-case, and closed: a consumer may match on these.
FLAG_POSITION_AMBIGUOUS = "position-ambiguous"
FLAG_SPANS_UNIT_JUNCTION = "spans-unit-junction"
FLAG_MOTIF_CONTEXT_DIVERGES = "motif-context-diverges"
FLAG_ALLELE_UNREPRESENTABLE = "allele-unrepresentable-in-vcf"
FLAG_LOW_READ_SUPPORT = "low-read-support"
FLAG_CALLER_DISAGREEMENT = "caller-disagreement"
FLAG_LENGTH_TRUNCATED = "length-truncated"

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

#: Motif symbol -> 60 bp unit on the genomic plus strand, as Kestrel sees it.
#:
#: Embedded rather than read from ``reference/MUC1_motifs_Rev_com.fa`` because that
#: directory is downloaded, not checked in: importing this module from a fresh clone
#: or an installed wheel must not depend on it. ``tests`` asserts this table is
#: byte-identical to the shipped FASTA whenever that FASTA is present, so the two
#: cannot drift.
#:
#: The 551-record pair reference is deliberately *not* embedded: every one of its
#: records is ``seq(R) ++ seq(L)`` over these same 34 motifs -- verified across all
#: 551 -- so :func:`pair_sequence` derives it and 69 kB of duplication is avoided.
MOTIFS: dict[str, str] = {
    "1": "CACAGCATTCTTCTCAGTAGAGCTGGGCACTGAACTTCTCTGGGTAGCCGAAGTCTCCTT",
    "2": "CTGAGTGGTGGAGGAGCCTGAACCGGGGCTGTGGCTGGAGAGTACGCTGCTGGTCATACT",
    "3": "CCAGGTGGCAGCTGAACCTGAAGCTGGTTCCGTGGCCGGGGCCAGAGTGACATCCTGTCC",
    "4": "TGGCGGGGTGGTGGAGCCCAGGGCTGGCCTGGTGACTGGGACCGAGGTGACATCCTGTCC",
    "4p": "TGGTGGGGTGGTGGAGCCCAGGGCTGGCCTGGTGACTGGGACCGAGGTGACATCCTGTCC",
    "5": "TGGGGGGGCGGTGGAGCCCGGGGCTGGCTTGTTGTCCGGGGCTGAGGTGACATCGTGGGC",
    "5C": "TTGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCTGAGGTGACATCGTGGGC",
    "X": "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "A": "TGCGGGCGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGGGCCGAGGTGACACCGTGGGC",
    "B": "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGGGCCGAGGTGACACCGTGGGC",
    "C": "TTGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "D": "TGGGGGGGCGGTGGAGCCCGGGGCGGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "E": "TGCGGGCGCGGTGGAGCCCGGGGCGGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "F": "TGTGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "G": "TGCGGGCGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "H": "TGGGGCGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "I": "TGGGGGCGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "J": "TGCGGGCGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGTGCCGAGGTGACACCGTGGGC",
    "K": "TGGGGGGGCGGTGGAGCCCAGGGCCGGCCTGCTCTCCGGGGCCGAGGTGACACCGTGGGC",
    "V": "TGGGGGTGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "W": "CGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "6": "CGGGGCCGGGGTGGAGCCCGGGGCCCGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "6p": "CGGGGCCGGGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "7": "CGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGCTGGGGGGGCGGTGGAGCC",
    "8": "CAAGGCGGGCCTGTTGTCCGGGGCCGAGGTGACACCATGGGCTGGGGGGGCGGTGGAGCC",
    "9": "TGAGCCTGATGCAGAGCCTGAGGCCGAGGTGACATTGTGGACTGGAGGGGCGGTGGAGCC",
    "L": "TGGGGCGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCTGAGGTGACACCGTGGGC",
    "M": "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGTGCCGAGGTGACACCGTGGGC",
    "N": "TGGGGGGGCGGTGGAGCCCGTGGCCGGCCTGCTCTCCGGGGCCGAGGTGACACCGTGGGC",
    "O": "TGGGGGGGCGGTGGAGCCTGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "P": "TGGGGGGGCGGTGGAGCCCGGGGCTGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC",
    "Q": "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACTGTGGGC",
    "R": "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGTGCCCGAGGTGACACCGTGGGC",
    "S": "TGCGGGGGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGGGCTGAGGTGACACCGTGGGC",
}

#: Where the embedded table came from, for the provenance test.
MOTIF_FASTA_NAME = "MUC1_motifs_Rev_com.fa"

#: adVNTR repeat unit -> the motif it is a rotation of.
#:
#: adVNTR ships nine repeat units. Only these three are rotations of a MUC1 motif
#: (each by 39 bases, which is what makes ``plus_from_ru`` a single formula). RU1,
#: RU4 and RU8 are 60 bp but match no motif; RU3, RU7 and RU9 are 78, 48 and 54 bp,
#: so parts of them have no counterpart in a 60 bp unit at all. Projecting a state
#: from any of the other six through this map would fabricate a coordinate, so those
#: states are reported without one.
MAPPABLE_RUS: dict[int, str] = {2: "X", 5: "V", 6: "C"}

#: Rotation offset shared by every mappable repeat unit.
_RU_ROTATION = 39

#: Emitted when adVNTR reports an insertion length without its sequence.
#:
#: An extension to the spec's flag vocabulary. It is needed because adVNTR records
#: only the first inserted base -- ``vntr_finder.py``: "If there are run of
#: insertions, the sequence might differ, but we just take the first base" -- so a
#: ``LEN4`` state constrains the length and nothing else. Without a distinct flag
#: this state is indistinguishable from a fully determined one.
FLAG_SEQUENCE_UNDETERMINED = "sequence-undetermined"


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


#: Reads supporting a call below which tier A is withheld. The benchmark's
#: `output.bam` splits reads across many pair records, so per-locus support is often
#: 1-3 even where `Estimated_Depth_AlternateVariant` aggregates to tens; a name is
#: only as good as the reads under it.
MIN_SUPPORT_FOR_TIER_A = 5


def reconcile(
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
        support: Reads supporting the call, when known. ``None`` means unknown, which
            is not the same as sufficient. Ignored when ``supports`` is given.
        supports: Reads per source, e.g. ``{"advntr": 24}``. Preferred: it binds the
            depth to the evidence it came from, so an unrelated well-covered
            observation cannot lend its depth to a thin agreement. The agreement is
            taken to be as strong as its weakest contributing source.

    Returns:
        Nomenclature: The reconciled call. Tier C when nothing was supplied or the
        callers disagree on the event class.
    """
    usable = [call for call in calls if call is not None]
    if not usable:
        return _undetermined("unknown", 0, "reconciled", ())

    primary = next((call for call in usable if call.source == "kestrel_vcf"), usable[0])
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

    # Independence is counted over the calls that actually carry the name, never over
    # every call supplied. Counting all of them let an *unnamed* call -- an adVNTR
    # state in an unmappable repeat unit, say -- donate a second `source` to two
    # duplicate rows from one caller, so a single caller's placement could be
    # promoted as though two had agreed on it.
    naming_sources = {call.source for call in named}
    agree = len({call.name for call in named}) == 1 and len(naming_sources) > 1
    if named and len({call.name for call in named}) > 1:
        flags.add(FLAG_CALLER_DISAGREEMENT)

    # Support must belong to the agreeing evidence. A sample-wide maximum would let a
    # well-covered but unrelated observation lend its depth to a 1-read agreement.
    effective_support = support
    if supports is not None:
        relevant = [supports.get(source) for source in naming_sources]
        present = [value for value in relevant if value is not None]
        effective_support = min(present) if present else None

    if effective_support is not None and effective_support < MIN_SUPPORT_FOR_TIER_A:
        flags.add(FLAG_LOW_READ_SUPPORT)

    tier = "B"
    if (
        agree
        and effective_support is not None
        and effective_support >= MIN_SUPPORT_FOR_TIER_A
        and FLAG_MOTIF_CONTEXT_DIVERGES not in flags
        and FLAG_SEQUENCE_UNDETERMINED not in flags
        and FLAG_CALLER_DISAGREEMENT not in flags
    ):
        tier = "A"

    if not named:
        return _undetermined(primary.event, primary.net_length, "reconciled", tuple(sorted(flags)))

    chosen = named[0]
    return Nomenclature(
        name=chosen.name,
        event=chosen.event,
        unit=chosen.unit,
        tier=tier,
        flags=tuple(sorted(flags)),
        ambiguity=chosen.ambiguity,
        repeat_form=chosen.repeat_form,
        net_length=chosen.net_length,
        source="reconciled" if len(naming_sources) > 1 else chosen.source,
    )


def render(call: Nomenclature) -> str:
    """Render a call as display text, showing only what its tier permits.

    The tier ladder is an emission rule, not just metadata: tier A shows the bare
    name, tier B shows the event and its ambiguity window but never a bare number,
    and tier C states the frameshift and stops. This function is the single place
    that rule is applied, so no surface can bypass it by reading ``name`` directly.

    Args:
        call: The reconciled call.

    Returns:
        str: Display text.
    """
    if call.tier == "A" and call.name:
        return call.name

    if call.tier == "C" or not call.name:
        sign = "+" if call.net_length >= 0 else "-"
        return f"frameshift {sign}{abs(call.net_length)}, allele undetermined"

    text = call.event
    if call.ambiguity:
        low, high = call.ambiguity
        text = f"{text} in {low}_{high}"
    return f"{text}, position-ambiguous"


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
    else:
        positions = sorted(UNIT_LENGTH + 1 - (((_RU_ROTATION - 1 + item.pos) % UNIT_LENGTH) + 1) for item in deletions)
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
    if pair is None or not 1 <= pos <= len(pair) or not ref or not alt:
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
    the reads must be produced by exactly the machinery that names a VCF record,
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
