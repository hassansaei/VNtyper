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
    "name_edit",
    "normalise",
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
