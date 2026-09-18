"""Construct repeat sequences and their truth without production inference imports.

Repeat count describes the explicitly supplied complete reference units. Small
edits change their bases, not that count. This module does not infer a clinical
label, caller identity, or biological independence from an edit or sequence hash.
"""

from __future__ import annotations

import hashlib
import json
import logging
import sys
from dataclasses import dataclass
from typing import NoReturn

logger = logging.getLogger(__name__)
_BASES = frozenset("ACGT")


@dataclass(frozen=True)
class RepeatEdit:
    """Edit at a zero-based unit/offset in the original, unedited repeat array."""

    repeat_index: int
    offset: int
    deleted_bases: str
    inserted_bases: str


@dataclass(frozen=True)
class SimulatedHaplotype:
    """Exact generated bases, construction evidence and unrounded repeat truth."""

    repeat_units: tuple[str, ...]
    repeat_unit_bp: int
    left_flank: str
    right_flank: str
    edits: tuple[RepeatEdit, ...]
    sequence: str
    repeat_count: int
    reference_repeat_bp: int
    repeat_start: int
    repeat_end: int
    sha256: str


@dataclass(frozen=True)
class DiploidTruth:
    """Phased construction counts; these do not assert independence or pathogenicity."""

    allele_repeat_counts: tuple[int, int]
    total_repeat_count: int
    haplotype_sha256: tuple[str, str]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _integer(value: object, label: str, *, minimum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _fail(f"simulation {label} must be an integer >= {minimum}")


def _dna(value: object, label: str, *, allow_empty: bool = False) -> None:
    if not isinstance(value, str) or (not value and not allow_empty) or not set(value) <= _BASES:
        _fail(f"simulation {label} must contain exact uppercase A/C/G/T bases")


def _edit_positions(
    units: tuple[str, ...], unit_bp: int, edits: tuple[RepeatEdit, ...]
) -> tuple[tuple[int, int, str], ...]:
    if not isinstance(edits, tuple):
        _fail("simulation edits must be an immutable ordered tuple")
    result = []
    unit_length_changes: dict[int, int] = {}
    previous_start = -1
    previous_end = -1
    for edit in edits:
        if not isinstance(edit, RepeatEdit):
            _fail("simulation edits must contain RepeatEdit values")
        _integer(edit.repeat_index, "repeat index", minimum=0)
        _integer(edit.offset, "edit offset", minimum=0)
        _dna(edit.deleted_bases, "deleted bases", allow_empty=True)
        _dna(edit.inserted_bases, "inserted bases", allow_empty=True)
        if edit.repeat_index >= len(units) or edit.offset > unit_bp:
            _fail("simulation edit lies outside its repeat unit")
        if edit.deleted_bases == edit.inserted_bases:
            _fail("simulation edit must change the sequence")
        stop = edit.offset + len(edit.deleted_bases)
        if stop > unit_bp or units[edit.repeat_index][edit.offset : stop] != edit.deleted_bases:
            _fail("simulation deleted bases differ from the original repeat unit")
        start = edit.repeat_index * unit_bp + edit.offset
        end = start + len(edit.deleted_bases)
        if start <= previous_start or start < previous_end:
            _fail("simulation edits must have unique increasing non-overlapping original positions")
        previous_start, previous_end = start, end
        unit_length_changes[edit.repeat_index] = (
            unit_length_changes.get(edit.repeat_index, 0) + len(edit.inserted_bases) - len(edit.deleted_bases)
        )
        result.append((start, end, edit.inserted_bases))
    if (
        any(abs(change) >= unit_bp for change in unit_length_changes.values())
        or abs(sum(unit_length_changes.values())) >= unit_bp
    ):
        _fail("simulation edits cannot change a whole unit count; construct different source units")
    return tuple(result)


def build_haplotype(
    *,
    repeat_units: tuple[str, ...],
    repeat_unit_bp: int,
    left_flank: str,
    right_flank: str,
    maximum_haplotype_bp: int,
    edits: tuple[RepeatEdit, ...] = (),
) -> SimulatedHaplotype:
    """Build exact repeat bases and original-unit truth under an explicit size budget.

    Args:
        repeat_units: Complete source units in their original order.
        repeat_unit_bp: Exact common source-unit length, excluding edits.
        left_flank: Exact unique left-flank sequence supplied by the protocol.
        right_flank: Exact unique right-flank sequence supplied by the protocol.
        maximum_haplotype_bp: Protocol resource limit; never an inferred length cap.
        edits: Sorted, non-overlapping changes in original unit coordinates.
            Net changes per unit and across the array must be smaller than a
            complete unit. Whole unit gains/losses belong in the source
            construction and its count.

    Returns:
        Immutable construction record including the edited sequence and digest.

    Raises:
        ValueError: If bases, coordinates, editing evidence or size are invalid.
    """
    _integer(repeat_unit_bp, "repeat unit length", minimum=1)
    _integer(maximum_haplotype_bp, "haplotype size budget", minimum=1)
    _dna(left_flank, "left flank")
    _dna(right_flank, "right flank")
    if not isinstance(repeat_units, tuple) or not repeat_units:
        _fail("simulation repeat units must be a non-empty immutable tuple")
    original_length = len(repeat_units) * repeat_unit_bp
    if original_length + len(left_flank) + len(right_flank) > maximum_haplotype_bp:
        _fail("simulation source haplotype exceeds its declared size budget")
    for unit in repeat_units:
        _dna(unit, "repeat unit")
        if len(unit) != repeat_unit_bp:
            _fail("simulation repeat unit length differs from its declared complete geometry")
    positions = _edit_positions(repeat_units, repeat_unit_bp, edits)
    edited_length = original_length + sum(len(inserted) - (end - start) for start, end, inserted in positions)
    if edited_length + len(left_flank) + len(right_flank) > maximum_haplotype_bp:
        _fail("simulation edited haplotype exceeds its declared size budget")
    original = "".join(repeat_units)
    parts: list[str] = []
    cursor = 0
    for start, end, inserted in positions:
        parts.extend((original[cursor:start], inserted))
        cursor = end
    parts.append(original[cursor:])
    sequence = left_flank + "".join(parts) + right_flank
    document = {
        "schema_version": "simulation-haplotype-v1",
        "repeat_units": list(repeat_units),
        "repeat_unit_bp": repeat_unit_bp,
        "left_flank": left_flank,
        "right_flank": right_flank,
        "edits": [[edit.repeat_index, edit.offset, edit.deleted_bases, edit.inserted_bases] for edit in edits],
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
    }
    digest = hashlib.sha256(
        json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("ascii")
    ).hexdigest()
    return SimulatedHaplotype(
        repeat_units,
        repeat_unit_bp,
        left_flank,
        right_flank,
        edits,
        sequence,
        len(repeat_units),
        original_length,
        len(left_flank),
        len(left_flank) + edited_length,
        digest,
    )


def diploid_truth(first: SimulatedHaplotype, second: SimulatedHaplotype) -> DiploidTruth:
    """Derive paired truth from independently revalidated sequence constructions.

    Args:
        first: First phased haplotype and construction record.
        second: Second phased haplotype and construction record.

    Returns:
        Integer allele counts and their total, with both construction identities.

    Raises:
        ValueError: If either construction record has been altered.
    """
    for haplotype in (first, second):
        if not isinstance(haplotype, SimulatedHaplotype):
            _fail("simulation diploid truth requires constructed haplotypes")
        for value in (
            haplotype.repeat_count,
            haplotype.reference_repeat_bp,
            haplotype.repeat_start,
            haplotype.repeat_end,
        ):
            _integer(value, "derived haplotype count or coordinate", minimum=0)
        expected = build_haplotype(
            repeat_units=haplotype.repeat_units,
            repeat_unit_bp=haplotype.repeat_unit_bp,
            left_flank=haplotype.left_flank,
            right_flank=haplotype.right_flank,
            maximum_haplotype_bp=sys.maxsize,
            edits=haplotype.edits,
        )
        if expected != haplotype:
            _fail("simulation haplotype differs from its original construction evidence")
    if first.repeat_unit_bp != second.repeat_unit_bp:
        _fail("simulation diploid haplotypes must share one repeat unit length")
    return DiploidTruth(
        (first.repeat_count, second.repeat_count),
        first.repeat_count + second.repeat_count,
        (first.sha256, second.sha256),
    )
