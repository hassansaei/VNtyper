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

__all__ = [
    "CANONICAL_UNIT",
    "UNIT_LENGTH",
    "ambiguity_interval",
    "name_edit",
    "normalise",
    "repeat_form",
    "revcomp",
]

#: The canonical MUC1 repeat unit in the coding orientation. Carries the tract Wenzel
#: (2018, PMID:29520014) publishes: 7xC at positions 53-59, ``A`` at 60. Verified
#: byte-identical to the MUC1 simulator's ``X`` unit.
CANONICAL_UNIT = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"

#: Every MUC1 repeat unit is 60 bp.
UNIT_LENGTH = 60

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


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
