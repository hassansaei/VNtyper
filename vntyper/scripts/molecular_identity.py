"""Validated value types for canonical MUC1 molecular identities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal, TypeAlias

IdentityReference: TypeAlias = Literal["MUC1-X-60-coding-v1"]
TranslationStatus: TypeAlias = Literal["resolved", "unresolved"]
TranslationFailure: TypeAlias = Literal[
    "invalid-allele",
    "missing-motif-context",
    "motif-anchor-mismatch",
    "pair-boundary-edit",
    "non-x-unit",
    "reconstruction-mismatch",
]
IdentitySource: TypeAlias = Literal["kestrel", "advntr"]
EvidenceDispositionValue: TypeAlias = Literal["admissible", "identity-insufficient"]

CANONICAL_MOLECULAR_IDENTITY_REFERENCE: Final[IdentityReference] = "MUC1-X-60-coding-v1"
CANONICAL_MUC1_X_CODING_UNIT: Final[str] = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"

_TRANSLATION_FAILURES: frozenset[str] = frozenset(
    {
        "invalid-allele",
        "missing-motif-context",
        "motif-anchor-mismatch",
        "pair-boundary-edit",
        "non-x-unit",
        "reconstruction-mismatch",
    }
)
_EVIDENCE_DISPOSITIONS: frozenset[str] = frozenset({"admissible", "identity-insufficient"})
_IDENTITY_TIERS: frozenset[str] = frozenset({"A", "B", "C", "abstained"})


@dataclass(frozen=True)
class CodingEdit:
    """A normalized edit in the canonical 60-base MUC1 coding repeat."""

    start: int
    end: int
    deleted: str
    inserted: str

    def __post_init__(self) -> None:
        """Validate canonical coding coordinates and alleles."""
        _validate_position(self.start, "start")
        _validate_allele(self.deleted, "deleted")
        _validate_allele(self.inserted, "inserted")

        if not self.deleted and not self.inserted:
            raise ValueError("Coding edit requires a deleted or inserted allele")
        if not self.deleted:
            if isinstance(self.end, bool) or not isinstance(self.end, int):
                raise ValueError("Insertion end must be an integer coordinate")
            if self.end != self.start - 1:
                raise ValueError("Insertion coordinates require end == start - 1")
        else:
            _validate_position(self.end, "end")
            if self.end < self.start:
                raise ValueError("Non-insertion coordinates require end >= start")
            if len(self.deleted) != self.end - self.start + 1:
                raise ValueError("Deleted allele length must match the coding interval")
        _validate_canonical_edit(self)


@dataclass(frozen=True)
class MolecularIdentity:
    """A versioned, normalized tuple of canonical coding edits."""

    reference: IdentityReference
    edits: tuple[CodingEdit, ...]

    def __post_init__(self) -> None:
        """Validate the reference version and canonical edit ordering."""
        if self.reference != CANONICAL_MOLECULAR_IDENTITY_REFERENCE:
            raise ValueError(f"Unsupported molecular identity reference: {self.reference}")
        if not isinstance(self.edits, tuple):
            raise ValueError("Molecular identity edits must be a tuple")
        if not self.edits:
            raise ValueError("Molecular identity requires at least one edit")
        if any(not isinstance(edit, CodingEdit) for edit in self.edits):
            raise ValueError("Molecular identity edits must be CodingEdit values")

        ordered_edits = tuple(sorted(self.edits, key=lambda edit: (edit.start, edit.end)))
        if self.edits != ordered_edits:
            raise ValueError("Molecular identity edits must be sorted by coordinate")
        for previous, current in zip(self.edits, self.edits[1:], strict=False):
            if _edits_collide(previous, current):
                raise ValueError("Molecular identity edits must not overlap or form a compound collision")


@dataclass(frozen=True)
class IdentityTranslation:
    """The resolved or closed-unresolved translation of one caller representation."""

    identity: MolecularIdentity | None
    status: TranslationStatus
    failure: TranslationFailure | None
    context_diverges: bool

    def __post_init__(self) -> None:
        """Require a consistent translation state."""
        if self.identity is not None and not isinstance(self.identity, MolecularIdentity):
            raise ValueError("Identity translation identity must be a MolecularIdentity or None")
        if not isinstance(self.context_diverges, bool):
            raise ValueError("Identity translation context divergence must be a boolean")
        if self.status not in {"resolved", "unresolved"}:
            raise ValueError(f"Unsupported molecular identity translation status: {self.status}")
        if self.failure is not None and self.failure not in _TRANSLATION_FAILURES:
            raise ValueError(f"Unsupported molecular identity translation failure: {self.failure}")
        if self.status == "resolved":
            if self.identity is None or self.failure is not None:
                raise ValueError("Resolved identity translations require an identity and no failure")
            return
        if self.identity is not None or self.failure is None:
            raise ValueError("Unresolved identity translations require no identity and a failure")
        if self.context_diverges:
            raise ValueError("Context divergence is meaningful only for resolved identity translations")


@dataclass(frozen=True)
class FrameConsequence:
    """The net length change and resulting mechanical frameshift state."""

    net_length_change: int
    is_frameshift: bool

    def __post_init__(self) -> None:
        """Validate that the mechanical frameshift state matches the length change."""
        if isinstance(self.net_length_change, bool) or not isinstance(self.net_length_change, int):
            raise ValueError("Net length change must be an integer")
        if not isinstance(self.is_frameshift, bool):
            raise ValueError("Frameshift state must be a boolean")
        if self.is_frameshift != (self.net_length_change % 3 != 0):
            raise ValueError("Frame consequence is inconsistent with net length change")


@dataclass(frozen=True)
class KestrelRepresentation:
    """A Kestrel VCF representation in its caller-specific reference context."""

    motifs: str
    position: int
    reference_allele: str
    alternate_allele: str

    def __post_init__(self) -> None:
        """Validate the complete raw Kestrel representation key."""
        if not isinstance(self.motifs, str) or not self.motifs:
            raise ValueError("Kestrel motifs must be a non-empty string")
        if isinstance(self.position, bool) or not isinstance(self.position, int) or self.position < 1:
            raise ValueError("Kestrel position must be a positive integer")
        _validate_allele(self.reference_allele, "Kestrel reference allele", allow_empty=False)
        _validate_allele(self.alternate_allele, "Kestrel alternate allele", allow_empty=False)


@dataclass(frozen=True)
class AdvntrRepresentation:
    """An adVNTR representation, retaining its caller-local State context."""

    state: str
    repeat_unit: str | None
    position: int | None

    def __post_init__(self) -> None:
        """Validate the available adVNTR state and optional repeat coordinates."""
        if not isinstance(self.state, str) or not self.state:
            raise ValueError("adVNTR State must be a non-empty string")
        if self.repeat_unit is not None and (not isinstance(self.repeat_unit, str) or not self.repeat_unit):
            raise ValueError("adVNTR repeat unit must be a non-empty string when supplied")
        if self.position is not None and (
            isinstance(self.position, bool) or not isinstance(self.position, int) or self.position < 1
        ):
            raise ValueError("adVNTR position must be a positive integer when supplied")


@dataclass(frozen=True)
class IdentityObservation:
    """One caller observation together with translation and frame-consequence evidence."""

    source: IdentitySource
    representation: KestrelRepresentation | AdvntrRepresentation
    translation: IdentityTranslation
    frame_consequence: FrameConsequence

    def __post_init__(self) -> None:
        """Keep a source bound to its own caller representation."""
        if self.source not in {"kestrel", "advntr"}:
            raise ValueError(f"Identity observation source is unsupported: {self.source}")
        if not isinstance(self.translation, IdentityTranslation):
            raise ValueError("Identity observation translation must be an IdentityTranslation")
        if not isinstance(self.frame_consequence, FrameConsequence):
            raise ValueError("Identity observation frame consequence must be a FrameConsequence")
        if self.source == "kestrel" and isinstance(self.representation, KestrelRepresentation):
            return
        if self.source == "advntr" and isinstance(self.representation, AdvntrRepresentation):
            return
        raise ValueError("Identity observation source does not match representation")


@dataclass(frozen=True)
class EvidenceDisposition:
    """Whether a caller observation may support molecular agreement."""

    value: EvidenceDispositionValue

    def __post_init__(self) -> None:
        """Validate the closed molecular-agreement evidence vocabulary."""
        if self.value not in _EVIDENCE_DISPOSITIONS:
            raise ValueError(f"Unsupported evidence disposition: {self.value}")


@dataclass(frozen=True)
class IdentityDecision:
    """A selected molecular identity or an explicit abstention."""

    identity: MolecularIdentity | None
    tier: str
    molecular_agreement: bool
    abstention_reason: str | None

    def __post_init__(self) -> None:
        """Keep selected-identity and abstention states distinct."""
        if self.identity is not None and not isinstance(self.identity, MolecularIdentity):
            raise ValueError("Identity decision identity must be a MolecularIdentity or None")
        if not isinstance(self.tier, str) or not self.tier:
            raise ValueError("Identity decision tier must be a non-empty string")
        if self.tier not in _IDENTITY_TIERS:
            raise ValueError(f"Identity decision tier is unsupported: {self.tier}")
        if not isinstance(self.molecular_agreement, bool):
            raise ValueError("Molecular agreement must be a boolean")
        if self.tier == "abstained":
            if self.identity is not None or self.molecular_agreement:
                raise ValueError("An abstained identity decision cannot select or agree on an identity")
            if not isinstance(self.abstention_reason, str) or not self.abstention_reason:
                raise ValueError("An abstained identity decision requires an abstention reason")
            return
        if self.identity is None:
            raise ValueError(f"Identity decision tier {self.tier} does not accept an identity-free decision")
        if self.abstention_reason is not None:
            raise ValueError("A selected identity decision cannot include an abstention reason")
        if self.tier == "A" and not self.molecular_agreement:
            raise ValueError("Tier A requires molecular agreement")


def make_coding_edit(start: int, end: int, deleted: str, inserted: str) -> CodingEdit:
    """Create one validated canonical coding edit.

    Args:
        start: One-based start position in the canonical 60-base repeat.
        end: Inclusive end position, or ``start - 1`` for an insertion.
        deleted: Uppercase canonical bases removed by the edit.
        inserted: Uppercase canonical bases introduced by the edit.

    Returns:
        A validated immutable coding edit.

    Raises:
        ValueError: If coordinates or alleles are not canonical.
    """
    return CodingEdit(start=start, end=end, deleted=deleted, inserted=inserted)


def make_molecular_identity(edits: tuple[CodingEdit, ...]) -> MolecularIdentity:
    """Create a validated identity on the canonical MUC1 X-repeat reference.

    Args:
        edits: Sorted, non-overlapping canonical coding edits.

    Returns:
        A validated immutable molecular identity.

    Raises:
        ValueError: If the edit tuple is empty, unordered, or overlapping.
    """
    return MolecularIdentity(reference=CANONICAL_MOLECULAR_IDENTITY_REFERENCE, edits=edits)


def parse_molecular_identity(value: str) -> MolecularIdentity:
    """Parse a strictly serialized canonical molecular identity.

    Args:
        value: The stable molecular-identity wire format.

    Returns:
        The validated molecular identity represented by ``value``.

    Raises:
        ValueError: If the reference version, edit fields, or canonical ordering is invalid.
    """
    if not isinstance(value, str):
        raise ValueError("Molecular identity serialization must be a string")
    reference, separator, encoded_edits = value.partition("|")
    if reference != CANONICAL_MOLECULAR_IDENTITY_REFERENCE:
        raise ValueError(f"Unsupported molecular identity reference: {reference}")
    if not separator or not encoded_edits:
        raise ValueError("Molecular identity serialization requires at least one edit")

    edits: list[CodingEdit] = []
    for encoded_edit in encoded_edits.split(";"):
        fields = encoded_edit.split("|")
        if len(fields) != 4:
            raise ValueError("Molecular identity edit serialization requires four fields")
        start_text, end_text, deleted, inserted = fields
        if not start_text or not end_text or not deleted or not inserted:
            raise ValueError("Molecular identity edit serialization has an empty field")
        start = _parse_coordinate(start_text)
        end = _parse_coordinate(end_text)
        edits.append(CodingEdit(start, end, _decode_allele(deleted), _decode_allele(inserted)))
    return make_molecular_identity(tuple(edits))


def serialize_molecular_identity(identity: MolecularIdentity) -> str:
    """Serialize a molecular identity in its stable canonical wire format.

    Args:
        identity: A validated canonical molecular identity.

    Returns:
        The versioned canonical serialization.
    """

    def allele(value: str) -> str:
        return value or "-"

    encoded = ";".join(
        f"{edit.start}|{edit.end}|{allele(edit.deleted)}|{allele(edit.inserted)}" for edit in identity.edits
    )
    return f"{identity.reference}|{encoded}"


def _validate_canonical_edit(edit: CodingEdit) -> None:
    """Reject edits that do not use the canonical reference or normalized representation."""
    expected_deleted = CANONICAL_MUC1_X_CODING_UNIT[edit.start - 1 : edit.end]
    if edit.deleted != expected_deleted:
        raise ValueError("Deleted allele does not match the canonical MUC1 X coding reference")

    canonical_edit = _canonicalize_edit(edit)
    if canonical_edit is None:
        raise ValueError("Coding edit must change the canonical MUC1 X coding reference")
    if canonical_edit != (edit.start, edit.end, edit.deleted, edit.inserted):
        raise ValueError("Coding edit must use the minimal 3-prime-most canonical representation")


def _canonicalize_edit(edit: CodingEdit) -> tuple[int, int, str, str] | None:
    """Return the minimal 3-prime-most representation of an edit's alternate unit."""
    alternate = (
        CANONICAL_MUC1_X_CODING_UNIT[: edit.start - 1] + edit.inserted + CANONICAL_MUC1_X_CODING_UNIT[edit.end :]
    )
    prefix_length = 0
    while (
        prefix_length < len(CANONICAL_MUC1_X_CODING_UNIT)
        and prefix_length < len(alternate)
        and CANONICAL_MUC1_X_CODING_UNIT[prefix_length] == alternate[prefix_length]
    ):
        prefix_length += 1

    reference_end = len(CANONICAL_MUC1_X_CODING_UNIT)
    alternate_end = len(alternate)
    while (
        reference_end > prefix_length
        and alternate_end > prefix_length
        and CANONICAL_MUC1_X_CODING_UNIT[reference_end - 1] == alternate[alternate_end - 1]
    ):
        reference_end -= 1
        alternate_end -= 1

    deleted = CANONICAL_MUC1_X_CODING_UNIT[prefix_length:reference_end]
    inserted = alternate[prefix_length:alternate_end]
    if not deleted and not inserted:
        return None
    start = prefix_length + 1
    end = reference_end if deleted else start - 1
    return start, end, deleted, inserted


def _edits_collide(previous: CodingEdit, current: CodingEdit) -> bool:
    """Return whether sorted edits overlap an insertion anchor or compound boundary."""
    previous_is_insertion = not previous.deleted
    current_is_insertion = not current.deleted
    if previous_is_insertion and current_is_insertion:
        return previous.start == current.start
    if previous_is_insertion:
        return current.start <= previous.start
    if current_is_insertion:
        return current.start <= previous.end + 1
    return current.start <= previous.end


def _parse_coordinate(value: str) -> int:
    """Parse one canonical decimal coordinate token from the identity wire format."""
    if not value.isascii() or not value.isdecimal() or (len(value) > 1 and value.startswith("0")):
        raise ValueError("Molecular identity edit coordinates must be canonical decimal integers")
    return int(value)


def _validate_position(value: int, name: str) -> None:
    """Validate a canonical one-based coding position."""
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 60:
        raise ValueError(f"Coding edit {name} must be an integer from 1 through 60")


def _validate_allele(value: str, name: str, *, allow_empty: bool = True) -> None:
    """Validate a canonical uppercase DNA allele."""
    if not isinstance(value, str) or (not allow_empty and not value):
        raise ValueError(f"{name} must be a {'non-empty ' if not allow_empty else ''}uppercase DNA allele")
    if any(base not in {"A", "C", "G", "T"} for base in value):
        raise ValueError(f"{name} must contain only uppercase A, C, G, and T")


def _decode_allele(value: str) -> str:
    """Decode the empty-allele sentinel from a serialized identity edit."""
    if value == "-":
        return ""
    return value
