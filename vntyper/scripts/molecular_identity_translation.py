"""Complete-context translation into canonical MUC1 molecular identities."""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final

from vntyper.modules.advntr.advntr_variant_annotations import derive_ru_and_pos
from vntyper.scripts.molecular_identity import (
    CANONICAL_MUC1_X_CODING_UNIT,
    AdvntrRepresentation,
    CodingEdit,
    IdentityTranslation,
    KestrelRepresentation,
    MolecularIdentity,
    TranslationFailure,
    make_coding_edit,
    make_molecular_identity,
)

_UNIT_LENGTH: Final[int] = len(CANONICAL_MUC1_X_CODING_UNIT)
_PAIR_LENGTH: Final[int] = _UNIT_LENGTH * 2
_DNA: Final[frozenset[str]] = frozenset("ACGT")
_COMPLEMENT: Final[dict[int, int]] = str.maketrans("ACGT", "TGCA")
_ADVNTR_INSERTION: Final[re.Pattern[str]] = re.compile(r"^I(\d+)_([0-9]+)_([ACGT])_LEN(\d+)$")
_ADVNTR_DELETION: Final[re.Pattern[str]] = re.compile(r"^D(\d+)_([0-9]+)$")


@dataclass(frozen=True)
class _RawEdit:
    """One unnormalised edit on a 60-base coding unit."""

    start: int
    end: int
    inserted: str


def translate_kestrel_representation(
    representation: KestrelRepresentation,
    motif_map: Mapping[str, str],
) -> IdentityTranslation:
    """Translate one complete Kestrel motif-pair representation independently.

    Args:
        representation: Raw VCF fields plus the observed 120-base reference pair.
        motif_map: Expected checked-in genomic-plus motif sequences, keyed by motif.

    Returns:
        A resolved canonical identity or a closed unresolved translation.
    """
    if not _is_dna(representation.reference_allele) or not _is_dna(representation.alternate_allele):
        return _unresolved("invalid-allele")

    expected_pair = _expected_pair(representation.motifs, motif_map)
    if expected_pair is None:
        return _unresolved("missing-motif-context")
    expected_halves = (expected_pair[:_UNIT_LENGTH], expected_pair[_UNIT_LENGTH:])

    start_index = representation.position - 1
    end_index = start_index + len(representation.reference_allele)
    if (
        start_index < 0
        or end_index > _PAIR_LENGTH
        or representation.pair_sequence[start_index:end_index] != representation.reference_allele
    ):
        return _unresolved("invalid-allele")

    changed_start, changed_end = _changed_reference_interval(
        representation.position,
        representation.reference_allele,
        representation.alternate_allele,
    )
    affected_half = _affected_half(changed_start, changed_end)
    if affected_half is None:
        return _unresolved("pair-boundary-edit")

    observed_halves = (representation.pair_sequence[:_UNIT_LENGTH], representation.pair_sequence[_UNIT_LENGTH:])
    neighbour_half = 1 - affected_half
    if observed_halves[neighbour_half] != expected_halves[neighbour_half]:
        return _unresolved("motif-anchor-mismatch")
    if observed_halves[affected_half] != expected_halves[affected_half]:
        return _unresolved("reconstruction-mismatch")

    alternate_pair = (
        representation.pair_sequence[:start_index]
        + representation.alternate_allele
        + representation.pair_sequence[end_index:]
    )
    length_change = len(representation.alternate_allele) - len(representation.reference_allele)
    alternate_halves = _split_alternate_pair(alternate_pair, affected_half, length_change)
    if alternate_halves is None or alternate_halves[neighbour_half] != expected_halves[neighbour_half]:
        return _unresolved("reconstruction-mismatch")

    source_unit = _reverse_complement(observed_halves[affected_half])
    alternate_source_unit = _reverse_complement(alternate_halves[affected_half])
    source_edit = _minimal_raw_edit(source_unit, alternate_source_unit)
    if source_edit is None:
        return _unresolved("reconstruction-mismatch")

    alternate_x = _apply_raw_edit(CANONICAL_MUC1_X_CODING_UNIT, source_edit)
    identity = resolve_coding_pair_edit(CANONICAL_MUC1_X_CODING_UNIT, alternate_x)
    if identity is None or _apply_identity(CANONICAL_MUC1_X_CODING_UNIT, identity) != alternate_x:
        return _unresolved("reconstruction-mismatch")

    return IdentityTranslation(
        identity=identity,
        status="resolved",
        failure=None,
        context_diverges=source_unit != CANONICAL_MUC1_X_CODING_UNIT,
    )


def translate_advntr_representation(
    representation: AdvntrRepresentation,
    repeat_unit_motifs: Mapping[str, str],
    rotation_offset: int,
) -> IdentityTranslation:
    """Translate complete parsed adVNTR State evidence on canonical repeat unit X.

    Args:
        representation: State text plus its independently parsed repeat-unit and
            position annotations.
        repeat_unit_motifs: Injected repeat-unit identifiers mapped to motif names.
        rotation_offset: Injected one-based rotation from adVNTR to genomic-plus
            motif orientation.

    Returns:
        A resolved canonical identity or a closed unresolved translation. A bare
        State label and first-base-only multibase insertions fail closed.
    """
    if representation.repeat_units is None or representation.positions is None:
        return _unresolved("missing-motif-context")
    if not _valid_advntr_context(repeat_unit_motifs, rotation_offset):
        return _unresolved("missing-motif-context")

    ru_annotations, position_annotations = derive_ru_and_pos((representation.state,))
    ru_tokens = tuple(ru_annotations[0].split(","))
    position_tokens = tuple(position_annotations[0].split(","))
    if "." in ru_tokens or "." in position_tokens:
        return _unresolved("invalid-allele")
    parsed_positions = tuple(int(position) for position in position_tokens)
    if representation.repeat_units != ru_tokens or representation.positions != parsed_positions:
        return _unresolved("reconstruction-mismatch")
    if any(repeat_unit_motifs.get(unit) != "X" for unit in ru_tokens):
        return _unresolved("non-x-unit")
    if any(position > _UNIT_LENGTH for position in parsed_positions):
        return _unresolved("non-x-unit")

    raw_edits: list[_RawEdit] = []
    deletion_positions: set[int] = set()
    for part in representation.state.split("&"):
        insertion = _ADVNTR_INSERTION.fullmatch(part.strip())
        deletion = _ADVNTR_DELETION.fullmatch(part.strip())
        if insertion is not None:
            length = int(insertion.group(4))
            if length != 1:
                return _unresolved("reconstruction-mismatch")
            position = int(insertion.group(1))
            if _advntr_plus_position(position, rotation_offset) == _UNIT_LENGTH:
                return _unresolved("pair-boundary-edit")
            raw_edits.append(_advntr_insertion(position, insertion.group(3), rotation_offset))
        elif deletion is not None:
            position = int(deletion.group(1))
            deletion_positions.add(position)
            raw_edits.append(_advntr_deletion(position, rotation_offset))
        else:
            return _unresolved("invalid-allele")
    if _deletions_cross_advntr_seam(deletion_positions, rotation_offset):
        return _unresolved("pair-boundary-edit")

    identity = _resolve_raw_edit_set(CANONICAL_MUC1_X_CODING_UNIT, tuple(raw_edits))
    if identity is None:
        return _unresolved("reconstruction-mismatch")
    return IdentityTranslation(identity=identity, status="resolved", failure=None, context_diverges=False)


def bind_bam_translation(
    kestrel_translation: IdentityTranslation,
    bam_identity: MolecularIdentity,
) -> IdentityTranslation | None:
    """Bind BAM evidence only to the complete existing Kestrel identity.

    Args:
        kestrel_translation: The independently translated Kestrel representation.
        bam_identity: The complete normalized edit set observed together in one BAM
            haplotype record.

    Returns:
        The unchanged Kestrel translation when the complete identities agree;
        otherwise ``None``. Binding creates neither a caller nor a representation.
    """
    if kestrel_translation.status != "resolved" or kestrel_translation.identity != bam_identity:
        return None
    return kestrel_translation


def resolve_coding_pair_edit(reference_unit: str, alternate_unit: str) -> MolecularIdentity | None:
    """Resolve one complete alternate X unit to its normalized canonical identity.

    Args:
        reference_unit: Complete 60-base coding-orientation reference unit.
        alternate_unit: Complete coding-orientation unit after the edit.

    Returns:
        The exact normalized identity, or ``None`` when reconstruction cannot prove it.
    """
    if reference_unit != CANONICAL_MUC1_X_CODING_UNIT or not _is_dna(alternate_unit):
        return None
    raw_edit = _minimal_raw_edit(reference_unit, alternate_unit)
    if raw_edit is None:
        return None

    normalized = _normalise_raw_edit(reference_unit, raw_edit)
    try:
        edit = make_coding_edit(
            normalized.start,
            normalized.end,
            reference_unit[normalized.start - 1 : normalized.end],
            normalized.inserted,
        )
        identity = make_molecular_identity((edit,))
    except ValueError:
        return None
    if _apply_identity(reference_unit, identity) != alternate_unit:
        return None
    return identity


def _unresolved(failure: TranslationFailure) -> IdentityTranslation:
    """Build one closed unresolved translation."""
    return IdentityTranslation(identity=None, status="unresolved", failure=failure, context_diverges=False)


def _is_dna(sequence: object) -> bool:
    """Return whether a value is a non-empty uppercase DNA allele."""
    return isinstance(sequence, str) and bool(sequence) and set(sequence) <= _DNA


def _valid_advntr_context(repeat_unit_motifs: Mapping[str, str], rotation_offset: int) -> bool:
    """Return whether injected adVNTR coordinate context is closed and usable."""
    return (
        isinstance(repeat_unit_motifs, Mapping)
        and bool(repeat_unit_motifs)
        and all(
            isinstance(unit, str) and bool(unit) and isinstance(motif, str) and bool(motif)
            for unit, motif in repeat_unit_motifs.items()
        )
        and not isinstance(rotation_offset, bool)
        and isinstance(rotation_offset, int)
        and 1 <= rotation_offset <= _UNIT_LENGTH
    )


def _reverse_complement(sequence: str) -> str:
    """Return uppercase DNA in reverse-complement orientation."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _expected_pair(motifs: str, motif_map: Mapping[str, str]) -> str | None:
    """Resolve the expected 120-base plus-strand pair, right motif then left."""
    parts = motifs.split("-")
    if len(parts) != 2 or not all(parts):
        return None
    left, right = parts
    left_sequence = motif_map.get(left)
    right_sequence = motif_map.get(right)
    if (
        not isinstance(left_sequence, str)
        or not isinstance(right_sequence, str)
        or not _is_dna(left_sequence)
        or not _is_dna(right_sequence)
        or len(left_sequence) != _UNIT_LENGTH
        or len(right_sequence) != _UNIT_LENGTH
    ):
        return None
    return right_sequence + left_sequence


def _changed_reference_interval(position: int, reference: str, alternate: str) -> tuple[int, int]:
    """Return the inclusive changed reference span, or an insertion gap."""
    prefix = 0
    while prefix < len(reference) and prefix < len(alternate) and reference[prefix] == alternate[prefix]:
        prefix += 1
    reference_end = len(reference)
    alternate_end = len(alternate)
    while (
        reference_end > prefix
        and alternate_end > prefix
        and reference[reference_end - 1] == alternate[alternate_end - 1]
    ):
        reference_end -= 1
        alternate_end -= 1
    start = position + prefix
    end = position + reference_end - 1
    return start, end


def _affected_half(start: int, end: int) -> int | None:
    """Return the affected plus-strand half, closing on the pair junction."""
    if end < start:
        gap = start - 1
        if gap == _UNIT_LENGTH:
            return None
        return 0 if gap < _UNIT_LENGTH else 1
    if start <= _UNIT_LENGTH < end:
        return None
    return 0 if end <= _UNIT_LENGTH else 1


def _split_alternate_pair(pair: str, affected_half: int, length_change: int) -> tuple[str, str] | None:
    """Split an alternate pair while retaining a 60-base neighbour anchor."""
    affected_length = _UNIT_LENGTH + length_change
    if affected_length < 1 or len(pair) != _PAIR_LENGTH + length_change:
        return None
    split = affected_length if affected_half == 0 else _UNIT_LENGTH
    halves = (pair[:split], pair[split:])
    if len(halves[1 - affected_half]) != _UNIT_LENGTH:
        return None
    return halves


def _minimal_raw_edit(reference: str, alternate: str) -> _RawEdit | None:
    """Trim common sequence to one minimal raw coding edit."""
    prefix = 0
    while prefix < len(reference) and prefix < len(alternate) and reference[prefix] == alternate[prefix]:
        prefix += 1
    reference_end = len(reference)
    alternate_end = len(alternate)
    while (
        reference_end > prefix
        and alternate_end > prefix
        and reference[reference_end - 1] == alternate[alternate_end - 1]
    ):
        reference_end -= 1
        alternate_end -= 1
    if prefix == reference_end and prefix == alternate_end:
        return None
    start = prefix + 1
    end = reference_end if reference_end > prefix else start - 1
    return _RawEdit(start, end, alternate[prefix:alternate_end])


def _normalise_raw_edit(reference: str, edit: _RawEdit) -> _RawEdit:
    """Shift an insertion 3-prime within the closed 1-through-60 reference."""
    if edit.end >= edit.start or not edit.inserted:
        return edit
    start = edit.start
    inserted = edit.inserted
    while start > len(reference) and reference[start - 2] == inserted[-1]:
        inserted = inserted[-1] + inserted[:-1]
        start -= 1
    while start < len(reference) and reference[start - 1] == inserted[0]:
        inserted = inserted[1:] + inserted[0]
        start += 1
    return _RawEdit(start, start - 1, inserted)


def _apply_raw_edit(reference: str, edit: _RawEdit) -> str:
    """Apply one raw edit using one-based inclusive coordinates."""
    return reference[: edit.start - 1] + edit.inserted + reference[edit.end :]


def _apply_raw_edits(reference: str, edits: tuple[_RawEdit, ...]) -> str:
    """Apply non-overlapping raw edits from the highest coordinate down."""
    alternate = reference
    for edit in sorted(edits, key=lambda item: (item.start, item.end), reverse=True):
        alternate = _apply_raw_edit(alternate, edit)
    return alternate


def _resolve_raw_edit_set(reference: str, raw_edits: tuple[_RawEdit, ...]) -> MolecularIdentity | None:
    """Normalize distinct co-occurring changes without fusing unchanged sequence."""
    expected_alternate = _apply_raw_edits(reference, raw_edits)
    normalized_edits: list[CodingEdit] = []
    for component in _raw_edit_components(raw_edits):
        single_alternate = _apply_raw_edits(reference, component)
        identity = resolve_coding_pair_edit(reference, single_alternate)
        if identity is None:
            return None
        normalized_edits.extend(identity.edits)
    try:
        complete_identity = make_molecular_identity(tuple(sorted(normalized_edits, key=lambda edit: edit.start)))
    except ValueError:
        return None
    if _apply_identity(reference, complete_identity) != expected_alternate:
        return None
    return complete_identity


def _raw_edit_components(raw_edits: tuple[_RawEdit, ...]) -> tuple[tuple[_RawEdit, ...], ...]:
    """Group overlapping or adjacent raw changes without joining separated edits."""
    components: list[list[_RawEdit]] = []
    for current in sorted(raw_edits, key=lambda edit: (edit.start, edit.end)):
        if components and any(_raw_edits_connect(previous, current) for previous in components[-1]):
            components[-1].append(current)
        else:
            components.append([current])
    return tuple(tuple(component) for component in components)


def _raw_edits_connect(left: _RawEdit, right: _RawEdit) -> bool:
    """Return whether two edits belong to one connected replacement component."""
    left_insertion = left.end < left.start
    right_insertion = right.end < right.start
    if left_insertion and right_insertion:
        return left.start == right.start
    if left_insertion:
        return right.start <= left.start <= right.end + 1
    if right_insertion:
        return left.start <= right.start <= left.end + 1
    return right.start <= left.end + 1 and left.start <= right.end + 1


def _apply_identity(reference: str, identity: MolecularIdentity) -> str:
    """Apply every normalized identity edit from the highest coordinate down."""
    alternate = reference
    for edit in reversed(identity.edits):
        alternate = alternate[: edit.start - 1] + edit.inserted + alternate[edit.end :]
    return alternate


def _advntr_plus_position(position: int, rotation_offset: int) -> int:
    """Rotate one adVNTR HMM position into the genomic-plus X unit."""
    return ((rotation_offset - 1 + position) % _UNIT_LENGTH) + 1


def _advntr_insertion(position: int, base: str, rotation_offset: int) -> _RawEdit:
    """Project one complete adVNTR single-base insertion to coding orientation."""
    coding_left = _UNIT_LENGTH - _advntr_plus_position(position, rotation_offset)
    return _RawEdit(coding_left + 1, coding_left, _reverse_complement(base))


def _advntr_deletion(position: int, rotation_offset: int) -> _RawEdit:
    """Project one adVNTR deletion to a coding-reference position."""
    coding_position = _UNIT_LENGTH + 1 - _advntr_plus_position(position, rotation_offset)
    return _RawEdit(coding_position, coding_position, "")


def _deletions_cross_advntr_seam(positions: set[int], rotation_offset: int) -> bool:
    """Detect adjacent HMM deletions whose projection wraps across coding 60/1."""
    for position in positions:
        next_position = (position % _UNIT_LENGTH) + 1
        if next_position not in positions:
            continue
        left = _advntr_plus_position(position, rotation_offset)
        right = _advntr_plus_position(next_position, rotation_offset)
        if {left, right} == {1, _UNIT_LENGTH}:
            return True
    return False
