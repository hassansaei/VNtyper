"""Pure identity-candidate capture and legacy-selection overlay decisions."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Literal, Protocol, TypeAlias, cast

from vntyper.scripts.molecular_identity import (
    AdvntrRepresentation,
    FrameConsequence,
    IdentityObservation,
    IdentityTranslation,
    KestrelRepresentation,
    MolecularIdentity,
)
from vntyper.scripts.molecular_identity_translation import (
    translate_advntr_representation,
    translate_kestrel_representation,
)

CandidateSource: TypeAlias = Literal["kestrel", "advntr"]
RawKeyValues: TypeAlias = tuple[str, int, str, str] | tuple[str, tuple[str, ...], tuple[int, ...]]
SupportValue: TypeAlias = int | float

LEGACY_GATE_COLUMNS: tuple[str, ...] = (
    "is_frameshift",
    "is_valid_frameshift",
    "depth_confidence_pass",
    "alt_filter_pass",
    "motif_filter_pass",
    "flag_filter_pass",
)


class IdentityTranslator(Protocol):
    """Explicit translation dependency consumed by caller capture."""

    def translate_kestrel(self, representation: KestrelRepresentation) -> IdentityTranslation:
        """Translate one complete Kestrel representation."""
        ...

    def translate_advntr(self, representation: AdvntrRepresentation) -> IdentityTranslation:
        """Translate one complete adVNTR representation."""
        ...


@dataclass(frozen=True)
class IdentityTranslationComponent:
    """Immutable configured adapter over the complete-context translators."""

    kestrel_motifs: Mapping[str, str]
    advntr_repeat_unit_motifs: Mapping[str, str]
    advntr_rotation_offset: int

    def __post_init__(self) -> None:
        """Copy and validate every injected translation authority."""
        motif_map = _validated_string_map(self.kestrel_motifs, "Kestrel motif map")
        repeat_unit_map = _validated_string_map(self.advntr_repeat_unit_motifs, "adVNTR repeat-unit map")
        if (
            isinstance(self.advntr_rotation_offset, bool)
            or not isinstance(self.advntr_rotation_offset, int)
            or not 1 <= self.advntr_rotation_offset <= 60
        ):
            raise ValueError("adVNTR rotation offset must be an integer from 1 through 60")
        object.__setattr__(self, "kestrel_motifs", MappingProxyType(motif_map))
        object.__setattr__(self, "advntr_repeat_unit_motifs", MappingProxyType(repeat_unit_map))

    def translate_kestrel(self, representation: KestrelRepresentation) -> IdentityTranslation:
        """Translate with the injected checked-in motif authority.

        Args:
            representation: One complete raw Kestrel representation.

        Returns:
            Its resolved or closed-unresolved translation.
        """
        return translate_kestrel_representation(representation, self.kestrel_motifs)

    def translate_advntr(self, representation: AdvntrRepresentation) -> IdentityTranslation:
        """Translate with the injected RU mapping and rotation.

        Args:
            representation: One complete parsed adVNTR State representation.

        Returns:
            Its resolved or closed-unresolved translation.
        """
        return translate_advntr_representation(
            representation,
            self.advntr_repeat_unit_motifs,
            self.advntr_rotation_offset,
        )


@dataclass(frozen=True)
class RawRepresentationKey:
    """A caller-local raw representation key, never a molecular identity."""

    source: CandidateSource
    values: RawKeyValues

    def __post_init__(self) -> None:
        """Validate the exact closed key shape for each caller."""
        if self.source == "kestrel":
            if not _is_kestrel_key(self.values):
                raise ValueError("Kestrel raw representation key must be (Motifs, POS, REF, ALT)")
            return
        if self.source == "advntr":
            if not _is_advntr_key(self.values):
                raise ValueError("adVNTR raw representation key must contain State and complete RU/POS tuples")
            return
        raise ValueError(f"Unsupported candidate source: {self.source}")


@dataclass(frozen=True)
class IdentityCandidate:
    """One independently translated caller representation and its evidence."""

    row_key: RawRepresentationKey
    observation: IdentityObservation
    support: SupportValue
    depth: SupportValue
    blocking_gates: frozenset[str] = frozenset()
    flags: frozenset[str] = frozenset()
    eligible: bool = False

    def __post_init__(self) -> None:
        """Keep source, scalar evidence, blockers, and eligibility typed."""
        if self.row_key.source != self.observation.source:
            raise ValueError("Candidate row-key source does not match its observation")
        _validate_support(self.support, "Candidate support")
        _validate_support(self.depth, "Candidate depth")
        if not isinstance(self.blocking_gates, frozenset) or any(
            gate not in LEGACY_GATE_COLUMNS for gate in self.blocking_gates
        ):
            raise ValueError("Candidate blocking gates must use the six legacy gate names")
        if not isinstance(self.flags, frozenset) or any(not isinstance(flag, str) or not flag for flag in self.flags):
            raise ValueError("Candidate flags must be non-empty strings")
        if not isinstance(self.eligible, bool):
            raise ValueError("Candidate eligibility must be a boolean")

    @property
    def identity(self) -> MolecularIdentity | None:
        """Return the resolved molecular identity, if translation succeeded."""
        return self.observation.translation.identity


@dataclass(frozen=True)
class IdentityCandidateGroup:
    """Equivalent representations with conservative unioned blockers."""

    identity: MolecularIdentity
    candidates: tuple[IdentityCandidate, ...]
    blocking_gates: frozenset[str]
    flags: frozenset[str]
    context_diverges: bool


@dataclass(frozen=True)
class IdentityCandidateSet:
    """A caller's complete captured candidates plus legacy selected-row overlay."""

    source: CandidateSource
    candidates: tuple[IdentityCandidate, ...]
    selected_row_key: RawRepresentationKey | None = None

    def __post_init__(self) -> None:
        """Require one source and a selected key drawn from the candidate set."""
        if self.source not in {"kestrel", "advntr"}:
            raise ValueError(f"Unsupported candidate-set source: {self.source}")
        if not isinstance(self.candidates, tuple):
            raise ValueError("Identity candidates must be a tuple")
        if any(candidate.row_key.source != self.source for candidate in self.candidates):
            raise ValueError("Identity candidate set cannot mix caller sources")
        if self.selected_row_key is not None:
            if self.selected_row_key.source != self.source:
                raise ValueError("Selected raw key source does not match its candidate set")
            if self.selected_row_key not in {candidate.row_key for candidate in self.candidates}:
                raise ValueError("Selected raw key is absent from its candidate set")

    @property
    def identity_hypothesis_count(self) -> int:
        """Count distinct resolved identities surviving caller eligibility."""
        return len({candidate.identity for candidate in self.candidates if candidate.eligible and candidate.identity})

    @property
    def selected_candidate(self) -> IdentityCandidate:
        """Return the exact legacy-selected row without re-ranking.

        Raises:
            ValueError: If selection has not been overlaid or the key is ambiguous.
        """
        if self.selected_row_key is None:
            raise ValueError("Identity candidate set has no selected row")
        matches = [candidate for candidate in self.candidates if candidate.row_key == self.selected_row_key]
        if len(matches) != 1:
            raise ValueError("Selected raw representation key is not unique")
        return matches[0]

    @property
    def selected_group(self) -> IdentityCandidateGroup:
        """Return the selected identity's equivalence group and unioned blockers.

        Raises:
            ValueError: If the selected row is unresolved.
        """
        selected = self.selected_candidate
        if selected.identity is None:
            raise ValueError("An unresolved selected row has no identity group")
        equivalents = tuple(candidate for candidate in self.candidates if candidate.identity == selected.identity)
        return IdentityCandidateGroup(
            identity=selected.identity,
            candidates=equivalents,
            blocking_gates=frozenset().union(*(candidate.blocking_gates for candidate in equivalents)),
            flags=frozenset().union(*(candidate.flags for candidate in equivalents)),
            context_diverges=any(candidate.observation.translation.context_diverges for candidate in equivalents),
        )

    @property
    def selected_support(self) -> SupportValue:
        """Return support from the exact legacy-selected row only."""
        return self.selected_candidate.support

    @property
    def selected_depth(self) -> SupportValue:
        """Return depth from the exact legacy-selected row only."""
        return self.selected_candidate.depth

    @property
    def equivalent_representation_count(self) -> int:
        """Count separately retained representations equivalent to the selected row."""
        return len(self.selected_group.candidates)


def translation_component_from_config(config: Mapping[str, object]) -> IdentityTranslationComponent:
    """Build one explicit translation component from nomenclature configuration.

    Args:
        config: Complete nomenclature configuration containing motif and adVNTR authorities.

    Returns:
        An immutable translation component owning copied configuration values.

    Raises:
        KeyError: If a required configuration field is absent.
        TypeError: If a required configuration section is not a mapping.
        ValueError: If a configured authority is malformed.
    """
    motifs = config["motifs"]
    advntr = config["advntr"]
    if not isinstance(motifs, Mapping) or not isinstance(advntr, Mapping):
        raise TypeError("Identity translation configuration sections must be mappings")
    repeat_unit_motifs = advntr["mappable_repeat_units"]
    rotation_offset = advntr["rotation_offset"]
    if not isinstance(repeat_unit_motifs, Mapping):
        raise TypeError("adVNTR repeat-unit mapping must be a mapping")
    return IdentityTranslationComponent(
        cast(Mapping[str, str], motifs),
        cast(Mapping[str, str], repeat_unit_motifs),
        cast(int, rotation_offset),
    )


def capture_kestrel_observations(
    rows: Iterable[Mapping[str, object]],
    component: IdentityTranslator,
) -> IdentityCandidateSet:
    """Translate every complete Kestrel row before legacy motif projection.

    Args:
        rows: Confidence-scored rows retaining raw pair context and support.
        component: Explicit translation dependency for this run.

    Returns:
        An unselected candidate set retaining one observation per input row.

    Raises:
        KeyError: If a required raw or evidence field is absent.
        ValueError: If a field cannot form the strict typed representation.
    """
    candidates: list[IdentityCandidate] = []
    for row in rows:
        representation = KestrelRepresentation(
            motifs=_required_string(row, "Motifs"),
            position=_required_int(row, "POS"),
            reference_allele=_required_string(row, "REF"),
            alternate_allele=_required_string(row, "ALT"),
            pair_sequence=_required_string(row, "Motif_sequence"),
        )
        net_change = len(representation.alternate_allele) - len(representation.reference_allele)
        observation = IdentityObservation(
            source="kestrel",
            representation=representation,
            translation=component.translate_kestrel(representation),
            frame_consequence=FrameConsequence(net_change, net_change % 3 != 0),
        )
        candidates.append(
            IdentityCandidate(
                row_key=_kestrel_key(representation),
                observation=observation,
                support=_required_support(row, "Estimated_Depth_AlternateVariant"),
                depth=_required_support(row, "Estimated_Depth_Variant_ActiveRegion"),
                flags=_flags(row.get("Flag")),
            )
        )
    return IdentityCandidateSet("kestrel", tuple(candidates))


def capture_advntr_observations(
    rows: Iterable[Mapping[str, object]],
    component: IdentityTranslator,
) -> IdentityCandidateSet:
    """Translate complete positive adVNTR State rows before reconciliation.

    Args:
        rows: Positive rows with State/Variant, complete RU/POS tuples, and support.
        component: Explicit translation dependency for this run.

    Returns:
        A caller-local candidate set whose positive rows are eligible hypotheses.

    Raises:
        KeyError: If required State context or evidence is absent.
        ValueError: If typed RU/POS or numeric evidence is malformed.
    """
    candidates: list[IdentityCandidate] = []
    for row in rows:
        state = _state(row)
        repeat_units = _string_tuple(row, "RU")
        positions = _integer_tuple(row, "POS")
        representation = AdvntrRepresentation(state, repeat_units, positions)
        net_change = _required_int(row, "Net_indel_length")
        observation = IdentityObservation(
            source="advntr",
            representation=representation,
            translation=component.translate_advntr(representation),
            frame_consequence=FrameConsequence(net_change, net_change % 3 != 0),
        )
        candidates.append(
            IdentityCandidate(
                row_key=RawRepresentationKey("advntr", (state, repeat_units, positions)),
                observation=observation,
                support=_required_support(row, "NumberOfSupportingReads"),
                depth=_required_support(row, "MeanCoverage"),
                flags=_flags(row.get("Flag")),
                eligible=True,
            )
        )
    return IdentityCandidateSet("advntr", tuple(candidates))


def with_candidate_evidence(
    captured: IdentityCandidateSet,
    rows: Iterable[Mapping[str, object]],
) -> IdentityCandidateSet:
    """Attach all six Kestrel gates and flags after legacy stages have marked rows.

    Args:
        captured: Pre-projection Kestrel observations.
        rows: Post-projection rows retaining the four-field raw key and all gates.

    Returns:
        The same observations with conservative raw-key-level blocker/flag unions.

    Raises:
        ValueError: If called for another source, a gate is absent, or evidence is incomplete.
    """
    if captured.source != "kestrel":
        raise ValueError("Six-gate evidence overlay is defined only for Kestrel candidates")
    blockers_by_key: dict[RawRepresentationKey, set[str]] = defaultdict(set)
    flags_by_key: dict[RawRepresentationKey, set[str]] = defaultdict(set)
    seen: set[RawRepresentationKey] = set()
    for row in rows:
        key = _kestrel_key_from_row(row)
        for gate in LEGACY_GATE_COLUMNS:
            if gate not in row:
                raise ValueError(f"Required legacy gate column is missing: {gate}")
            value = row[gate]
            if not isinstance(value, bool):
                raise ValueError(f"Legacy gate {gate} must be a boolean")
            if not value:
                blockers_by_key[key].add(gate)
        flags_by_key[key].update(_flags(row.get("Flag")))
        seen.add(key)
    missing = {candidate.row_key for candidate in captured.candidates} - seen
    if missing:
        raise ValueError("Post-projection evidence is missing a captured Kestrel raw key")
    return IdentityCandidateSet(
        source="kestrel",
        candidates=tuple(
            replace(
                candidate,
                blocking_gates=frozenset(blockers_by_key[candidate.row_key]),
                flags=candidate.flags | frozenset(flags_by_key[candidate.row_key]),
            )
            for candidate in captured.candidates
        ),
        selected_row_key=captured.selected_row_key,
    )


def overlay_legacy_projection(
    captured: IdentityCandidateSet,
    legacy_passing_keys: Iterable[RawRepresentationKey],
    legacy_selected: RawRepresentationKey,
) -> IdentityCandidateSet:
    """Overlay unchanged legacy eligibility and selection without ranking candidates.

    Args:
        captured: Independently translated caller observations.
        legacy_passing_keys: Raw keys surviving the existing complete gate policy.
        legacy_selected: Exact raw key returned by the unchanged legacy selector.

    Returns:
        A candidate set with eligibility and selected-row metadata overlaid.

    Raises:
        ValueError: If keys are unknown or the selected row did not pass.
    """
    passing = frozenset(legacy_passing_keys)
    known = {candidate.row_key for candidate in captured.candidates}
    if not passing <= known:
        raise ValueError("Legacy projection contains an unknown raw representation key")
    if legacy_selected not in passing:
        raise ValueError("The selected Kestrel row did not pass the complete legacy projection")
    candidates = tuple(replace(candidate, eligible=candidate.row_key in passing) for candidate in captured.candidates)
    return IdentityCandidateSet(captured.source, candidates, legacy_selected)


def _validated_string_map(value: Mapping[str, str], name: str) -> dict[str, str]:
    """Copy one non-empty string mapping or fail closed."""
    if not isinstance(value, Mapping) or not value:
        raise ValueError(f"{name} must be a non-empty mapping")
    copied = dict(value)
    if any(
        not isinstance(key, str) or not key or not isinstance(item, str) or not item for key, item in copied.items()
    ):
        raise ValueError(f"{name} must contain non-empty string keys and values")
    return copied


def _validate_support(value: object, name: str) -> None:
    """Require finite non-negative numeric evidence."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be a finite non-negative number")


def _required_support(row: Mapping[str, object], column: str) -> SupportValue:
    """Read and validate one numeric evidence cell."""
    value = row[column]
    _validate_support(value, column)
    return cast(SupportValue, value)


def _required_string(row: Mapping[str, object], column: str, *, allow_empty: bool = False) -> str:
    """Read one strict string cell."""
    value = row[column]
    if not isinstance(value, str) or (not allow_empty and not value):
        raise ValueError(f"{column} must be a {'possibly empty ' if allow_empty else 'non-empty '}string")
    return value


def _required_int(row: Mapping[str, object], column: str) -> int:
    """Read one strict integer cell."""
    value = row[column]
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{column} must be an integer")
    return value


def _state(row: Mapping[str, object]) -> str:
    """Read the sole State/Variant value without ambiguous fallback."""
    values = [row[column] for column in ("State", "Variant") if column in row]
    if len(values) != 1 or not isinstance(values[0], str) or not values[0]:
        raise ValueError("adVNTR row must contain exactly one non-empty State or Variant field")
    return values[0]


def _string_tuple(row: Mapping[str, object], column: str) -> tuple[str, ...]:
    """Parse a complete comma-separated string tuple."""
    value = _required_string(row, column)
    values = tuple(value.split(","))
    if any(not item for item in values):
        raise ValueError(f"{column} must contain complete comma-separated values")
    return values


def _integer_tuple(row: Mapping[str, object], column: str) -> tuple[int, ...]:
    """Parse a complete comma-separated positive-integer tuple."""
    values = _string_tuple(row, column)
    if any(not value.isascii() or not value.isdecimal() for value in values):
        raise ValueError(f"{column} must contain comma-separated positive integers")
    parsed = tuple(int(value) for value in values)
    if any(value < 1 for value in parsed):
        raise ValueError(f"{column} must contain comma-separated positive integers")
    return parsed


def _kestrel_key(representation: KestrelRepresentation) -> RawRepresentationKey:
    """Build the exact four-field Kestrel representation key."""
    return RawRepresentationKey(
        "kestrel",
        (
            representation.motifs,
            representation.position,
            representation.reference_allele,
            representation.alternate_allele,
        ),
    )


def _kestrel_key_from_row(row: Mapping[str, object]) -> RawRepresentationKey:
    """Build the exact four-field key from a later Kestrel row."""
    return RawRepresentationKey(
        "kestrel",
        (
            _required_string(row, "Motifs"),
            _required_int(row, "POS"),
            _required_string(row, "REF"),
            _required_string(row, "ALT"),
        ),
    )


def _is_kestrel_key(values: object) -> bool:
    """Return whether values have the exact Kestrel key shape."""
    return (
        isinstance(values, tuple)
        and len(values) == 4
        and isinstance(values[0], str)
        and bool(values[0])
        and isinstance(values[1], int)
        and not isinstance(values[1], bool)
        and values[1] >= 1
        and isinstance(values[2], str)
        and bool(values[2])
        and isinstance(values[3], str)
        and bool(values[3])
    )


def _is_advntr_key(values: object) -> bool:
    """Return whether values have the exact complete adVNTR key shape."""
    return (
        isinstance(values, tuple)
        and len(values) == 3
        and isinstance(values[0], str)
        and bool(values[0])
        and isinstance(values[1], tuple)
        and bool(values[1])
        and all(isinstance(unit, str) and unit for unit in values[1])
        and isinstance(values[2], tuple)
        and len(values[1]) == len(values[2])
        and all(
            isinstance(position, int) and not isinstance(position, bool) and position >= 1 for position in values[2]
        )
    )


def _flags(value: object) -> frozenset[str]:
    """Parse non-placeholder semicolon-delimited flags."""
    if value is None or value in {"", "Not flagged", "Not applicable"}:
        return frozenset()
    if not isinstance(value, str):
        raise ValueError("Candidate Flag must be a string when present")
    return frozenset(flag for flag in value.split(";") if flag)
