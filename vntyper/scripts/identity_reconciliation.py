"""Pure molecular-identity keyed reconciliation decisions."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

from vntyper.scripts.molecular_identity import (
    EvidenceDisposition,
    IdentityDecision,
    IdentityTranslation,
    MolecularIdentity,
)

if TYPE_CHECKING:
    from vntyper.scripts.identity_candidates import IdentityTranslator

_CALLER_OF = {"kestrel_vcf": "kestrel", "kestrel_bam": "kestrel", "advntr": "advntr"}
_TIER_A_BLOCKING_FLAGS = frozenset({"motif-context-diverges", "sequence-undetermined", "caller-disagreement"})
_LEGACY_GATE_NAMES = frozenset(
    {
        "is_frameshift",
        "is_valid_frameshift",
        "depth_confidence_pass",
        "alt_filter_pass",
        "motif_filter_pass",
        "flag_filter_pass",
    }
)


class ReconciliationPresentationCall(Protocol):
    """The legacy presentation fields required by the identity adapter."""

    @property
    def name(self) -> str | None:
        """Return the legacy display name."""
        ...

    @property
    def source(self) -> str:
        """Return the evidence source."""
        ...

    @property
    def event(self) -> str:
        """Return the legacy event word."""
        ...

    @property
    def net_length(self) -> int:
        """Return the event's signed net length."""
        ...

    @property
    def flags(self) -> tuple[str, ...]:
        """Return the legacy presentation flags."""
        ...


@dataclass(frozen=True)
class IdentityReconciliationPolicy:
    """The separately named source-unit thresholds for identity reconciliation."""

    kestrel_min_alternate_kmer_path_depth: int
    advntr_min_sequencing_read_support: int

    def __post_init__(self) -> None:
        """Require positive integer thresholds without conflating their units."""
        _positive_integer(
            self.kestrel_min_alternate_kmer_path_depth,
            "Kestrel minimum alternate k-mer-path depth",
        )
        _positive_integer(
            self.advntr_min_sequencing_read_support,
            "adVNTR minimum sequencing-read support",
        )


@dataclass(frozen=True)
class IdentityReconciliationObservation:
    """One display observation with explicit identity and source-bound evidence."""

    translation: IdentityTranslation
    display_name: str | None
    source: str
    event: str
    net_length: int
    flags: frozenset[str]
    disposition: EvidenceDisposition
    known_variant_match: bool
    kestrel_alternate_kmer_path_depth: int | None = None
    advntr_sequencing_read_support: int | None = None
    blocking_gates: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        """Reject malformed or cross-wired source evidence."""
        if not isinstance(self.translation, IdentityTranslation):
            raise ValueError("Reconciliation translation must be an IdentityTranslation")
        if self.display_name is not None and (not isinstance(self.display_name, str) or not self.display_name):
            raise ValueError("Reconciliation display name must be a non-empty string or None")
        if not isinstance(self.source, str) or not self.source:
            raise ValueError("Reconciliation source must be a non-empty string")
        if not isinstance(self.event, str) or not self.event:
            raise ValueError("Reconciliation event must be a non-empty string")
        if isinstance(self.net_length, bool) or not isinstance(self.net_length, int):
            raise ValueError("Reconciliation net length must be an integer")
        if not isinstance(self.flags, frozenset) or any(not isinstance(flag, str) or not flag for flag in self.flags):
            raise ValueError("Reconciliation flags must be a frozenset of non-empty strings")
        if not isinstance(self.disposition, EvidenceDisposition):
            raise ValueError("Reconciliation disposition must be an EvidenceDisposition")
        if not isinstance(self.known_variant_match, bool):
            raise ValueError("Known-variant match must be a boolean")
        if not isinstance(self.blocking_gates, frozenset) or not self.blocking_gates <= _LEGACY_GATE_NAMES:
            raise ValueError("Reconciliation blockers must use the six legacy gate names")
        _optional_nonnegative_integer(
            self.kestrel_alternate_kmer_path_depth,
            "Kestrel alternate k-mer-path depth",
        )
        _optional_nonnegative_integer(
            self.advntr_sequencing_read_support,
            "adVNTR sequencing-read support",
        )
        if self.source == "kestrel_vcf":
            if self.advntr_sequencing_read_support is not None:
                raise ValueError("Kestrel VCF observations cannot carry adVNTR sequencing-read support")
        elif self.source == "advntr":
            if self.kestrel_alternate_kmer_path_depth is not None:
                raise ValueError("adVNTR observations cannot carry Kestrel alternate k-mer-path depth")
        elif self.kestrel_alternate_kmer_path_depth is not None or self.advntr_sequencing_read_support is not None:
            raise ValueError("Only the owning caller source may carry a named Tier-A evidence quantity")

    @property
    def identity(self) -> MolecularIdentity | None:
        """Return the resolved identity, if any."""
        return self.translation.identity

    @property
    def backs_identity(self) -> bool:
        """Return whether this observation may back molecular agreement."""
        return self.identity is not None and self.disposition.value == "admissible"


@dataclass(frozen=True)
class IdentityReconciliationResult:
    """Identity decision plus the exact legacy presentation/backing projection."""

    decision: IdentityDecision
    selected_observation_index: int | None
    backing_observation_indices: tuple[int, ...]
    backing_sources: tuple[str, ...]
    backing_callers: tuple[str, ...]
    caller_disagreement: bool
    event_disagreement: bool
    low_support_sources: tuple[str, ...]

    def __post_init__(self) -> None:
        """Validate selection and deterministic reconciliation metadata."""
        if not isinstance(self.decision, IdentityDecision):
            raise ValueError("Reconciliation result decision must be an IdentityDecision")
        if self.selected_observation_index is not None and (
            isinstance(self.selected_observation_index, bool)
            or not isinstance(self.selected_observation_index, int)
            or self.selected_observation_index < 0
        ):
            raise ValueError("Selected reconciliation observation index must be non-negative")
        for values, name in (
            (self.backing_observation_indices, "Backing observation indices"),
            (self.backing_sources, "Backing sources"),
            (self.backing_callers, "Backing callers"),
            (self.low_support_sources, "Low-support sources"),
        ):
            if not isinstance(values, tuple):
                raise ValueError(f"{name} must be a tuple")
        if any(
            isinstance(index, bool) or not isinstance(index, int) or index < 0
            for index in self.backing_observation_indices
        ):
            raise ValueError("Backing observation indices must be non-negative integers")
        if any(
            not isinstance(value, str) or not value
            for value in (*self.backing_sources, *self.backing_callers, *self.low_support_sources)
        ):
            raise ValueError("Reconciliation source and caller metadata must contain non-empty strings")
        if not isinstance(self.caller_disagreement, bool):
            raise ValueError("Caller disagreement must be a boolean")
        if not isinstance(self.event_disagreement, bool):
            raise ValueError("Event disagreement must be a boolean")
        if (
            len(set(self.backing_observation_indices)) != len(self.backing_observation_indices)
            or tuple(sorted(self.backing_observation_indices)) != self.backing_observation_indices
        ):
            raise ValueError("Backing observation indices must be unique and ordered")
        for values, name in (
            (self.backing_sources, "Backing sources"),
            (self.backing_callers, "Backing callers"),
            (self.low_support_sources, "Low-support sources"),
        ):
            if len(set(values)) != len(values):
                raise ValueError(f"{name} must be unique")
        if self.decision.identity is None:
            if self.selected_observation_index is not None or any(
                (self.backing_observation_indices, self.backing_sources, self.backing_callers, self.low_support_sources)
            ):
                raise ValueError("An identity-free reconciliation cannot carry backing metadata")
            return
        if (
            self.selected_observation_index is None
            or self.selected_observation_index not in self.backing_observation_indices
            or not self.backing_sources
            or not self.backing_callers
        ):
            raise ValueError("A selected identity requires complete backing metadata")
        expected_callers = tuple(dict.fromkeys(_caller(source) for source in self.backing_sources))
        if self.backing_callers != expected_callers:
            raise ValueError("Backing callers must match the backing source ownership")
        if not set(self.low_support_sources) <= set(self.backing_sources):
            raise ValueError("Low-support sources must back the selected identity")
        if self.decision.molecular_agreement != (len(self.backing_callers) >= 2):
            raise ValueError("Molecular agreement must match independent backing callers")

    @property
    def identity(self) -> MolecularIdentity | None:
        """Return the selected molecular identity."""
        return self.decision.identity

    @property
    def tier(self) -> str:
        """Return the identity decision tier."""
        return self.decision.tier

    @property
    def molecular_agreement(self) -> bool:
        """Return whether independent callers back the selected identity."""
        return self.decision.molecular_agreement


def build_identity_reconciliation_observations(
    kestrel_rows: Sequence[Mapping[str, object]],
    vcf_calls: Sequence[ReconciliationPresentationCall | None],
    bam_calls: Sequence[ReconciliationPresentationCall | None],
    advntr_rows: Sequence[Mapping[str, object]],
    advntr_calls_by_row: Sequence[Sequence[ReconciliationPresentationCall]],
    translation_component: IdentityTranslator,
    known_variant_names: frozenset[str],
) -> tuple[IdentityReconciliationObservation, ...] | None:
    """Adapt current caller rows to typed identity observations without I/O.

    Args:
        kestrel_rows: Positive Kestrel rows retaining A3 internal metadata.
        vcf_calls: Legacy Kestrel VCF display calls aligned with those rows.
        bam_calls: Optional legacy BAM refinement calls, not identity-bound until A5.
        advntr_rows: Complete positive adVNTR rows with State/RU/POS context.
        advntr_calls_by_row: Legacy adVNTR display calls aligned with those rows.
        translation_component: Explicit checked-in translation authority resolved at
            the stage boundary.
        known_variant_names: Configured known display names used only for tier gating.

    Returns:
        Typed observations, or ``None`` for a deliberate pre-A3 legacy artifact.

    Raises:
        ValueError: If current-run internal metadata or caller context is malformed.
    """
    from vntyper.scripts.identity_candidate_persistence import (
        IDENTITY_CAPTURE_COLUMNS,
        IDENTITY_SELECTION_COLUMNS,
        parse_selected_candidate_cells,
    )
    from vntyper.scripts.identity_candidates import capture_advntr_observations

    identity_columns = frozenset((*IDENTITY_CAPTURE_COLUMNS, *IDENTITY_SELECTION_COLUMNS))
    observed_columns = frozenset().union(*(frozenset(row.keys()) for row in kestrel_rows))
    present = identity_columns & observed_columns
    if not present:
        return None
    if len(kestrel_rows) != len(vcf_calls) or any(not identity_columns.issubset(row.keys()) for row in kestrel_rows):
        raise ValueError("Kestrel identity metadata is incomplete or misaligned")
    if len(advntr_rows) != len(advntr_calls_by_row):
        raise ValueError("adVNTR identity context is misaligned with its display calls")

    observations: list[IdentityReconciliationObservation] = []
    for row, call in zip(kestrel_rows, vcf_calls, strict=True):
        if call is None:
            continue
        persisted = parse_selected_candidate_cells(row)
        translation = persisted.translation
        if translation.identity is not None and translation.context_diverges != persisted.group_context_diverges:
            translation = IdentityTranslation(translation.identity, "resolved", None, persisted.group_context_diverges)
        observations.append(
            _make_observation(
                call,
                translation,
                known_variant_names,
                kmer_depth=_optional_integer_cell(row.get("Estimated_Depth_AlternateVariant")),
                extra_flags=persisted.flags,
                blocking_gates=persisted.blocking_gates,
            )
        )

    observations.extend(
        _make_observation(
            call,
            IdentityTranslation(None, "unresolved", "missing-motif-context", False),
            known_variant_names,
        )
        for call in bam_calls
        if call is not None
    )

    capture_rows: list[dict[str, object]] = []
    for row, calls in zip(advntr_rows, advntr_calls_by_row, strict=True):
        support = _optional_integer_cell(row.get("NumberOfSupportingReads"))
        coverage = _optional_integer_cell(row.get("MeanCoverage"))
        if support is None or coverage is None:
            raise ValueError("Current adVNTR identity reconciliation requires numeric support and coverage")
        capture_rows.append(
            {
                "Variant": str(row.get("Variant", "")),
                "RU": str(row.get("RU", "")),
                "POS": str(row.get("POS", "")),
                "Net_indel_length": sum(call.net_length for call in calls),
                "NumberOfSupportingReads": support,
                "MeanCoverage": coverage,
                "Flag": str(row.get("Flag", "")),
            }
        )
    candidates = capture_advntr_observations(capture_rows, translation_component)
    for row, calls, candidate in zip(advntr_rows, advntr_calls_by_row, candidates.candidates, strict=True):
        observations.extend(
            _make_observation(
                call,
                candidate.observation.translation,
                known_variant_names,
                advntr_reads=_optional_integer_cell(row.get("NumberOfSupportingReads")),
                extra_flags=candidate.flags,
            )
            for call in calls
        )
    return tuple(observations)


def reconcile_identity_observations(
    observations: tuple[IdentityReconciliationObservation, ...],
    policy: IdentityReconciliationPolicy,
) -> IdentityReconciliationResult:
    """Reconcile caller evidence by resolved canonical identity.

    Args:
        observations: Legacy-ordered observations. Kestrel VCF remains the explicit
            primary presentation when no unique corroborated identity outvotes it.
        policy: Separately named thresholds in each caller's evidence unit.

    Returns:
        The typed identity decision and exact selection/backing metadata.

    Raises:
        ValueError: If the inputs are not validated immutable reconciliation values.
    """
    if not isinstance(observations, tuple):
        raise ValueError("Identity reconciliation observations must be a tuple")
    if not isinstance(policy, IdentityReconciliationPolicy):
        raise ValueError("Identity reconciliation requires an explicit policy")
    if any(not isinstance(observation, IdentityReconciliationObservation) for observation in observations):
        raise ValueError("Identity reconciliation inputs must be typed observations")

    event_disagreement = len({observation.event for observation in observations}) > 1
    backing: dict[MolecularIdentity, list[int]] = {}
    for index, observation in enumerate(observations):
        if observation.backs_identity:
            assert observation.identity is not None
            backing.setdefault(observation.identity, []).append(index)

    if not backing:
        return IdentityReconciliationResult(
            decision=IdentityDecision(None, "abstained", False, "identity-unresolved"),
            selected_observation_index=None,
            backing_observation_indices=(),
            backing_sources=(),
            backing_callers=(),
            caller_disagreement=event_disagreement,
            event_disagreement=event_disagreement,
            low_support_sources=(),
        )

    corroborated = [identity for identity, indices in backing.items() if len(_callers(observations, indices)) >= 2]
    primary_index = next(
        (
            index
            for index, observation in enumerate(observations)
            if observation.source == "kestrel_vcf" and observation.backs_identity
        ),
        next(iter(backing.values()))[0],
    )
    chosen_identity = observations[primary_index].identity
    assert chosen_identity is not None
    if len(corroborated) == 1:
        chosen_identity = corroborated[0]

    selected_index = backing[chosen_identity][0]
    backing_indices = tuple(backing[chosen_identity])
    backing_sources = tuple(dict.fromkeys(observations[index].source for index in backing_indices))
    backing_callers = tuple(dict.fromkeys(_caller(source) for source in backing_sources))
    molecular_agreement = len(backing_callers) >= 2
    caller_disagreement = len(backing) > 1 or event_disagreement
    low_support_sources, source_support_sufficient = _source_support(
        observations,
        backing_indices,
        policy,
    )

    chosen = observations[selected_index]
    flags = frozenset().union(*(observation.flags for observation in observations))
    if caller_disagreement:
        flags |= {"caller-disagreement"}
    context_diverges = any(observations[index].translation.context_diverges for index in backing_indices)
    blocking_gates = frozenset().union(*(observations[index].blocking_gates for index in backing_indices))
    tier = "C" if event_disagreement else "B"
    if (
        tier == "B"
        and molecular_agreement
        and chosen.known_variant_match
        and source_support_sufficient
        and not context_diverges
        and not blocking_gates
        and not (_TIER_A_BLOCKING_FLAGS & flags)
    ):
        tier = "A"

    return IdentityReconciliationResult(
        decision=IdentityDecision(chosen_identity, tier, molecular_agreement, None),
        selected_observation_index=selected_index,
        backing_observation_indices=backing_indices,
        backing_sources=backing_sources,
        backing_callers=backing_callers,
        caller_disagreement=caller_disagreement,
        event_disagreement=event_disagreement,
        low_support_sources=low_support_sources,
    )


def _source_support(
    observations: tuple[IdentityReconciliationObservation, ...],
    backing_indices: tuple[int, ...],
    policy: IdentityReconciliationPolicy,
) -> tuple[tuple[str, ...], bool]:
    """Evaluate only the selected identity's two provenance-preserving quantities."""
    kestrel_values = [
        observations[index].kestrel_alternate_kmer_path_depth
        for index in backing_indices
        if observations[index].source == "kestrel_vcf"
    ]
    advntr_values = [
        observations[index].advntr_sequencing_read_support
        for index in backing_indices
        if observations[index].source == "advntr"
    ]
    all_values = (*kestrel_values, *advntr_values)
    complete = bool(kestrel_values and advntr_values) and all(value is not None for value in all_values)
    if not complete:
        return (), False
    kestrel_depth = min(value for value in kestrel_values if value is not None)
    advntr_reads = min(value for value in advntr_values if value is not None)
    low_sources: list[str] = []
    if kestrel_depth < policy.kestrel_min_alternate_kmer_path_depth:
        low_sources.append("kestrel_vcf")
    if advntr_reads < policy.advntr_min_sequencing_read_support:
        low_sources.append("advntr")
    return tuple(low_sources), not low_sources


def _make_observation(
    call: ReconciliationPresentationCall,
    translation: IdentityTranslation,
    known_variant_names: frozenset[str],
    *,
    kmer_depth: int | None = None,
    advntr_reads: int | None = None,
    extra_flags: frozenset[str] = frozenset(),
    blocking_gates: frozenset[str] = frozenset(),
) -> IdentityReconciliationObservation:
    return IdentityReconciliationObservation(
        translation=translation,
        display_name=call.name,
        source=call.source,
        event=call.event,
        net_length=call.net_length,
        flags=frozenset(call.flags) | extra_flags,
        disposition=EvidenceDisposition("admissible"),
        known_variant_match=call.name in known_variant_names,
        kestrel_alternate_kmer_path_depth=kmer_depth,
        advntr_sequencing_read_support=advntr_reads,
        blocking_gates=blocking_gates,
    )


def _optional_integer_cell(value: object) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _callers(
    observations: tuple[IdentityReconciliationObservation, ...],
    indices: list[int],
) -> set[str]:
    return {_caller(observations[index].source) for index in indices}


def _caller(source: str) -> str:
    return _CALLER_OF.get(source, source)


def _positive_integer(value: object, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _optional_nonnegative_integer(value: object, name: str) -> None:
    if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 0):
        raise ValueError(f"{name} must be a non-negative integer or None")
