"""Pure canonical-identity evidence retained from Kestrel BAM records."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Protocol

from vntyper.scripts.molecular_identity import CodingEdit, IdentityTranslation, KestrelRepresentation, MolecularIdentity

_MAXIMUM_MINIMUM_KMER_DEPTH = 2_147_483_647
_DNA = frozenset("ACGT")


class BamIdentityTranslator(Protocol):
    """Complete-context Kestrel translation dependency used for BAM edits."""

    def translate_kestrel(self, representation: KestrelRepresentation) -> IdentityTranslation:
        """Translate one complete pair representation."""
        ...


@dataclass(frozen=True)
class BamEditObservation:
    """One CIGAR-derived non-matching edit on the plus-strand pair reference.

    Attributes:
        start: Zero-based reference position.
        ref_span: Number of reference bases replaced.
        inserted: Number of haplotype bases inserted.
        bases: Complete inserted haplotype sequence in plus orientation.
    """

    start: int
    ref_span: int
    inserted: int
    bases: str

    def __post_init__(self) -> None:
        """Validate a complete, sequence-bearing edit projection."""
        _nonnegative_integer(self.start, "BAM edit start")
        _nonnegative_integer(self.ref_span, "BAM edit reference span")
        _nonnegative_integer(self.inserted, "BAM edit inserted span")
        if self.ref_span == 0 and self.inserted == 0:
            raise ValueError("BAM edit must change at least one reference or haplotype base")
        if not isinstance(self.bases, str) or len(self.bases) != self.inserted or not set(self.bases) <= _DNA:
            raise ValueError("BAM edit bases must be complete uppercase DNA matching the inserted span")


@dataclass(frozen=True)
class BamRecordObservation:
    """One eligible resolved haplotype record before identity binding.

    Attributes:
        edits: Complete non-matching edit observations consulted at this locus.
        minimum_kmer_depth: Optional typed XD value for this record.
    """

    edits: tuple[BamEditObservation, ...]
    minimum_kmer_depth: int | None

    def __post_init__(self) -> None:
        """Keep record edits immutable and XD typed without interpreting either."""
        if not isinstance(self.edits, tuple) or any(not isinstance(edit, BamEditObservation) for edit in self.edits):
            raise ValueError("BAM record edits must be a tuple of BamEditObservation values")
        _optional_minimum_kmer_depth(self.minimum_kmer_depth)


@dataclass(frozen=True)
class BamCandidateBinding:
    """One already captured Kestrel representation eligible for BAM binding.

    Attributes:
        observation_ordinal: Stable A3 Kestrel observation ordinal.
        translation: Complete-context translation already owned by that representation.
    """

    observation_ordinal: int
    translation: IdentityTranslation

    def __post_init__(self) -> None:
        """Validate representation identity without resolving an unresolved candidate."""
        _nonnegative_integer(self.observation_ordinal, "Kestrel observation ordinal")
        if not isinstance(self.translation, IdentityTranslation):
            raise ValueError("BAM candidate translation must be an IdentityTranslation")


@dataclass(frozen=True)
class BamIdentityLocus:
    """Complete pair context and existing Kestrel candidates for one BAM locus.

    Attributes:
        motifs: Kestrel pair label.
        pair_sequence: Complete observed 120-base pair sequence.
        candidates: Existing Kestrel representations; BAM cannot add to this tuple.
    """

    motifs: str
    pair_sequence: str
    candidates: tuple[BamCandidateBinding, ...]

    def __post_init__(self) -> None:
        """Validate complete pair context and stable unique candidate ordinals."""
        if not isinstance(self.motifs, str) or not self.motifs:
            raise ValueError("BAM identity locus motifs must be a non-empty string")
        if (
            not isinstance(self.pair_sequence, str)
            or len(self.pair_sequence) != 120
            or not set(self.pair_sequence) <= _DNA
        ):
            raise ValueError("BAM identity locus pair sequence must be exactly 120 uppercase DNA bases")
        if not isinstance(self.candidates, tuple) or any(
            not isinstance(candidate, BamCandidateBinding) for candidate in self.candidates
        ):
            raise ValueError("BAM identity locus candidates must be a tuple of BamCandidateBinding values")
        ordinals = tuple(candidate.observation_ordinal for candidate in self.candidates)
        if len(ordinals) != len(set(ordinals)):
            raise ValueError("BAM identity locus candidate ordinals must be unique")


@dataclass(frozen=True)
class BamIdentityEvidence:
    """Bound canonical identities retained for one eligible haplotype record.

    Attributes:
        identities: Existing Kestrel candidate identities supported by this record,
            deduplicated in candidate encounter order.
        candidate_observation_ordinals: Existing representation ordinals bound to
            each aligned identity. Equivalent representations remain metadata and
            never expand the record's vote.
        minimum_kmer_depth: Optional typed XD value retained observationally.
    """

    identities: tuple[MolecularIdentity, ...]
    candidate_observation_ordinals: tuple[tuple[int, ...], ...]
    minimum_kmer_depth: int | None

    def __post_init__(self) -> None:
        """Validate identity uniqueness, aligned bindings, and optional XD."""
        if not isinstance(self.identities, tuple) or any(
            not isinstance(identity, MolecularIdentity) for identity in self.identities
        ):
            raise ValueError("BAM record identities must be a tuple of MolecularIdentity values")
        if len(self.identities) != len(set(self.identities)):
            raise ValueError("BAM record identities must be unique")
        if not isinstance(self.candidate_observation_ordinals, tuple) or len(
            self.candidate_observation_ordinals
        ) != len(self.identities):
            raise ValueError("BAM record identity bindings must be aligned with identities")
        all_ordinals: list[int] = []
        for ordinals in self.candidate_observation_ordinals:
            if not isinstance(ordinals, tuple) or not ordinals:
                raise ValueError("Each BAM record identity requires aligned candidate bindings")
            for ordinal in ordinals:
                _nonnegative_integer(ordinal, "Bound Kestrel observation ordinal")
            if len(ordinals) != len(set(ordinals)):
                raise ValueError("Bound Kestrel observation ordinals must be unique")
            all_ordinals.extend(ordinals)
        if len(all_ordinals) != len(set(all_ordinals)):
            raise ValueError("One Kestrel observation ordinal cannot bind more than one identity")
        _optional_minimum_kmer_depth(self.minimum_kmer_depth)


@dataclass(frozen=True)
class BamLocusEvidence:
    """Replayable unweighted identity evidence for every eligible locus record.

    Attributes:
        records: Per-record bound identities and XD in encounter order.
        eligible_record_count: Exact denominator consulted at this locus.
        counts: Per-identity unweighted record counts.
    """

    records: tuple[BamIdentityEvidence, ...]
    eligible_record_count: int
    counts: Mapping[MolecularIdentity, int]

    def __post_init__(self) -> None:
        """Own immutable values and reject counts inconsistent with record sets."""
        if not isinstance(self.records, tuple) or any(
            not isinstance(record, BamIdentityEvidence) for record in self.records
        ):
            raise ValueError("BAM locus records must be a tuple of BamIdentityEvidence values")
        _nonnegative_integer(self.eligible_record_count, "BAM eligible record count")
        if self.eligible_record_count != len(self.records):
            raise ValueError("BAM eligible record count must equal the retained record count")
        if not isinstance(self.counts, Mapping):
            raise ValueError("BAM identity counts must be a mapping")
        supplied: dict[MolecularIdentity, int] = {}
        for identity, count in self.counts.items():
            if not isinstance(identity, MolecularIdentity):
                raise ValueError("BAM identity count keys must be MolecularIdentity values")
            _positive_integer(count, "BAM identity record count")
            supplied[identity] = count
        derived = _record_counts(self.records)
        if supplied != derived:
            raise ValueError("BAM identity counts must equal the derived record counts")
        object.__setattr__(self, "counts", MappingProxyType(derived))

    @property
    def record_identity_sets(self) -> tuple[tuple[MolecularIdentity, ...], ...]:
        """Return every eligible record's bound identity set in encounter order."""
        return tuple(record.identities for record in self.records)

    @property
    def xd_by_record(self) -> tuple[int | None, ...]:
        """Return optional XD once per eligible record in encounter order."""
        return tuple(record.minimum_kmer_depth for record in self.records)

    @property
    def bindings_by_record(self) -> tuple[tuple[tuple[int, ...], ...], ...]:
        """Return existing representation bindings aligned to record identities."""
        return tuple(record.candidate_observation_ordinals for record in self.records)

    @property
    def winning_identity(self) -> MolecularIdentity | None:
        """Return the unique unweighted record-count winner, abstaining on a tie."""
        if not self.counts:
            return None
        highest = max(self.counts.values())
        winners = tuple(identity for identity, count in self.counts.items() if count == highest)
        return winners[0] if len(winners) == 1 else None


def collect_locus_evidence(
    records: Iterable[BamRecordObservation],
    locus: BamIdentityLocus,
    component: BamIdentityTranslator,
) -> BamLocusEvidence:
    """Bind eligible records only to existing complete Kestrel candidate identities.

    One record contributes at most one vote to each bound identity, but may back
    several distinct co-occurring identities. XD is retained without affecting
    binding, counts, ordering, or winner selection.

    Args:
        records: Eligible resolved haplotype records in BAM encounter order.
        locus: Complete pair context and already captured candidate representations.
        component: Injected complete-context identity translation authority.

    Returns:
        Immutable per-record evidence and deterministic unweighted counts.

    Raises:
        ValueError: If an input is not one of the validated evidence values.
    """
    if not isinstance(locus, BamIdentityLocus):
        raise ValueError("BAM identity locus must be a BamIdentityLocus")
    retained: list[BamIdentityEvidence] = []
    for record in records:
        if not isinstance(record, BamRecordObservation):
            raise ValueError("BAM records must contain BamRecordObservation values")
        retained.append(_bind_record(record, locus, component))
    retained_records = tuple(retained)
    return BamLocusEvidence(retained_records, len(retained_records), _record_counts(retained_records))


def project_record_observation(
    edits: Iterable[tuple[int, int, int, str]],
    minimum_kmer_depth: int | None,
    window_start: int,
    window_end: int,
) -> BamRecordObservation:
    """Project raw parsed edits into one complete fail-closed locus record.

    Args:
        edits: Parsed ``(start, ref span, inserted span, bases)`` tuples.
        minimum_kmer_depth: Optional typed XD retained for the denominator record.
        window_start: Inclusive zero-based locus-window start.
        window_end: Inclusive zero-based locus-window end.

    Returns:
        Complete in-window edit observations, or an empty identity context if any
        consulted edit is incomplete.

    Raises:
        ValueError: If window boundaries or XD are invalid.
    """
    _nonnegative_integer(window_start, "BAM record projection window start")
    _nonnegative_integer(window_end, "BAM record projection window end")
    if window_end < window_start:
        raise ValueError("BAM record projection window end must not precede its start")
    observations: list[BamEditObservation] = []
    for raw_edit in edits:
        if not isinstance(raw_edit, tuple) or len(raw_edit) != 4:
            return BamRecordObservation((), minimum_kmer_depth)
        try:
            observation = BamEditObservation(*raw_edit)
        except (TypeError, ValueError):
            return BamRecordObservation((), minimum_kmer_depth)
        if window_start <= observation.start <= window_end:
            observations.append(observation)
    return BamRecordObservation(tuple(observations), minimum_kmer_depth)


def bind_complete_winner_identity(
    records: tuple[BamRecordObservation, ...],
    evidence: BamLocusEvidence,
    winning_edit: BamEditObservation,
    locus: BamIdentityLocus,
    component: BamIdentityTranslator,
    supporting_record_count: int,
) -> MolecularIdentity | None:
    """Bind one exact winner only when every complete supporting record proves it.

    Args:
        records: Complete observed records aligned with ``evidence``.
        evidence: Candidate-bound identity sets for those same records.
        winning_edit: Exact legacy winning edit.
        locus: Complete pair context and existing Kestrel candidates.
        component: Injected complete-context identity translation authority.
        supporting_record_count: Legacy unweighted support for ``winning_edit``.

    Returns:
        The sole common exact candidate identity, or ``None`` on unresolved,
        ambiguous, inconsistent, or non-matching evidence.

    Raises:
        ValueError: If boundary values are not validated evidence types.
    """
    if not isinstance(records, tuple) or any(not isinstance(record, BamRecordObservation) for record in records):
        raise ValueError("BAM winner record observations must be a tuple of BamRecordObservation values")
    if not isinstance(evidence, BamLocusEvidence):
        raise ValueError("BAM winner locus evidence must be a BamLocusEvidence")
    if len(records) != len(evidence.records):
        raise ValueError("BAM winner record observations and locus evidence must be aligned")
    if not isinstance(winning_edit, BamEditObservation):
        raise ValueError("BAM winning edit must be a BamEditObservation")
    if not isinstance(locus, BamIdentityLocus):
        raise ValueError("BAM winner identity locus must be a BamIdentityLocus")
    _nonnegative_integer(supporting_record_count, "BAM supporting record count")

    supporting_sets = [
        set(bound.identities)
        for record, bound in zip(records, evidence.records, strict=True)
        if winning_edit in record.edits
    ]
    if len(supporting_sets) != supporting_record_count or not supporting_sets:
        return None
    common_identities = set.intersection(*supporting_sets)
    if len(common_identities) != 1:
        return None
    common_identity = next(iter(common_identities))
    isolated_evidence = collect_locus_evidence(
        (BamRecordObservation((winning_edit,), None),),
        locus,
        component,
    )
    exact_winner_identities = isolated_evidence.record_identity_sets[0]
    return common_identity if exact_winner_identities == (common_identity,) else None


def _bind_record(
    record: BamRecordObservation,
    locus: BamIdentityLocus,
    component: BamIdentityTranslator,
) -> BamIdentityEvidence:
    from vntyper.scripts.molecular_identity_translation import bind_bam_translation

    record_edits = _translated_record_edits(record, locus, component)
    identities: list[MolecularIdentity] = []
    bindings: list[list[int]] = []
    identity_indices: dict[MolecularIdentity, int] = {}
    if record_edits is not None:
        for candidate in locus.candidates:
            identity = candidate.translation.identity
            if identity is None or not set(identity.edits) <= record_edits:
                continue
            if bind_bam_translation(candidate.translation, identity) is None:  # pragma: no cover - defensive
                continue
            index = identity_indices.get(identity)
            if index is None:
                identity_indices[identity] = len(identities)
                identities.append(identity)
                bindings.append([candidate.observation_ordinal])
            else:
                bindings[index].append(candidate.observation_ordinal)
    return BamIdentityEvidence(
        tuple(identities),
        tuple(tuple(ordinals) for ordinals in bindings),
        record.minimum_kmer_depth,
    )


def _translated_record_edits(
    record: BamRecordObservation,
    locus: BamIdentityLocus,
    component: BamIdentityTranslator,
) -> set[CodingEdit] | None:
    translated: set[CodingEdit] = set()
    for edit in record.edits:
        representation = _edit_representation(edit, locus)
        if representation is None:
            return None
        translation = component.translate_kestrel(representation)
        if translation.identity is None:
            return None
        translated.update(translation.identity.edits)
    return translated


def _edit_representation(edit: BamEditObservation, locus: BamIdentityLocus) -> KestrelRepresentation | None:
    pair = locus.pair_sequence
    if edit.start > len(pair) or edit.start + edit.ref_span > len(pair):
        return None
    if edit.ref_span == 0:
        if edit.start == len(pair):
            anchor = len(pair) - 1
            reference = pair[anchor]
            alternate = reference + edit.bases
        elif edit.start == 0:
            anchor = 0
            reference = pair[anchor]
            alternate = edit.bases + reference
        else:
            anchor = edit.start - 1
            reference = pair[anchor]
            alternate = reference + edit.bases
        position = anchor + 1
    elif edit.inserted == 0:
        if edit.start == 0:
            if edit.ref_span == len(pair):
                return None
            position = 1
            reference = pair[: edit.ref_span + 1]
            alternate = pair[edit.ref_span]
        else:
            position = edit.start
            reference = pair[edit.start - 1 : edit.start + edit.ref_span]
            alternate = pair[edit.start - 1]
    else:
        position = edit.start + 1
        reference = pair[edit.start : edit.start + edit.ref_span]
        alternate = edit.bases
    return KestrelRepresentation(locus.motifs, position, reference, alternate, pair)


def _record_counts(records: tuple[BamIdentityEvidence, ...]) -> dict[MolecularIdentity, int]:
    counts: dict[MolecularIdentity, int] = {}
    for record in records:
        for identity in record.identities:
            counts[identity] = counts.get(identity, 0) + 1
    return counts


def _nonnegative_integer(value: object, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")


def _positive_integer(value: object, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _optional_minimum_kmer_depth(value: object) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= _MAXIMUM_MINIMUM_KMER_DEPTH:
        raise ValueError("BAM minimum k-mer depth must be None or an integer in the signed Java range")
