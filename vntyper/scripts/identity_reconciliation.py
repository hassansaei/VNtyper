"""Pure molecular-identity keyed reconciliation decisions."""

from __future__ import annotations

import math
import re
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
_CANONICAL_NONNEGATIVE_INTEGER = re.compile(r"(?:0|[1-9][0-9]*)\Z")
_CANONICAL_NONNEGATIVE_NUMBER = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?\Z")
_ADMISSIBLE_DISPOSITION = EvidenceDisposition("admissible")
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
    advntr_mean_coverage: int | float | None = None
    blocking_gates: frozenset[str] = frozenset()
    presentation_call_index: int | None = None

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
        _optional_nonnegative_number(self.advntr_mean_coverage, "adVNTR MeanCoverage")
        _optional_nonnegative_integer(self.presentation_call_index, "Presentation call index")
        if self.source == "kestrel_vcf":
            if self.advntr_sequencing_read_support is not None or self.advntr_mean_coverage is not None:
                raise ValueError("Kestrel VCF observations cannot carry adVNTR evidence")
        elif self.source == "advntr":
            if self.kestrel_alternate_kmer_path_depth is not None:
                raise ValueError("adVNTR observations cannot carry Kestrel alternate k-mer-path depth")
        elif (
            self.kestrel_alternate_kmer_path_depth is not None
            or self.advntr_sequencing_read_support is not None
            or self.advntr_mean_coverage is not None
        ):
            raise ValueError("Only the owning caller source may carry a named Tier-A evidence quantity")

    @property
    def identity(self) -> MolecularIdentity | None:
        """Return the resolved identity, if any."""
        return self.translation.identity

    @property
    def backs_identity(self) -> bool:
        """Return whether this observation may back molecular agreement."""
        return self.identity is not None and self.disposition.value == "admissible"

    @property
    def is_unbound_advntr_row(self) -> bool:
        """Return whether this is one complete compound adVNTR row observation."""
        return (
            self.presentation_call_index is None
            and self.source == "advntr"
            and self.display_name is None
            and self.identity is not None
            and self.disposition.value == "admissible"
            and self.known_variant_match is False
            and self.advntr_sequencing_read_support is not None
            and self.advntr_mean_coverage is not None
            and self.kestrel_alternate_kmer_path_depth is None
        )


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
    presentation_selected_observation_index: int | None = None

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
        _optional_nonnegative_integer(
            self.presentation_selected_observation_index,
            "Presentation-selected observation index",
        )
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
        if self.event_disagreement and not self.caller_disagreement:
            raise ValueError("Event disagreement requires caller disagreement")
        if self.decision.tier == "A" and (
            self.caller_disagreement or self.event_disagreement or self.low_support_sources
        ):
            raise ValueError("Tier A cannot carry conflict or low-support metadata")
        if self.decision.identity is not None and ((self.decision.tier == "C") != self.event_disagreement):
            raise ValueError("Tier C must match event disagreement for a selected identity")
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
        if (
            self.presentation_selected_observation_index is not None
            and self.selected_observation_index != self.presentation_selected_observation_index
        ):
            raise ValueError("Selected identity metadata must align with the presentation-selected observation")

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


def select_compatibility_presentation_index(
    calls: Sequence[ReconciliationPresentationCall],
) -> int | None:
    """Select the exact legacy display row without making a molecular claim.

    Args:
        calls: Presentation calls in the existing VCF, BAM, adVNTR order.

    Returns:
        The selected call index, or ``None`` when no call was supplied.
    """
    if not calls:
        return None
    primary_index = next((index for index, call in enumerate(calls) if call.source == "kestrel_vcf"), 0)
    if len({call.event for call in calls}) > 1:
        return primary_index

    backing: dict[str, list[int]] = {}
    for index, call in enumerate(calls):
        if call.name is not None:
            backing.setdefault(call.name, []).append(index)
    corroborated = [
        name for name, indices in backing.items() if len({_caller(calls[index].source) for index in indices}) >= 2
    ]
    if len(corroborated) == 1:
        return backing[corroborated[0]][0]
    return next((index for index, call in enumerate(calls) if call.name is not None), primary_index)


def build_identity_reconciliation_observations(
    kestrel_rows: Sequence[Mapping[str, object]],
    vcf_calls: Sequence[ReconciliationPresentationCall | None],
    bam_calls: Sequence[ReconciliationPresentationCall | None],
    advntr_rows: Sequence[Mapping[str, object]],
    advntr_calls_by_row: Sequence[Sequence[ReconciliationPresentationCall]],
    translation_component: IdentityTranslator,
    known_variant_names: frozenset[str],
    *,
    bam_translations: Sequence[IdentityTranslation | None] | None = None,
) -> tuple[IdentityReconciliationObservation, ...] | None:
    """Adapt current caller rows to typed identity observations without I/O.

    Args:
        kestrel_rows: Positive Kestrel rows retaining A3 internal metadata.
        vcf_calls: Legacy Kestrel VCF display calls aligned with those rows.
        bam_calls: Optional legacy BAM refinement calls aligned to Kestrel rows.
        advntr_rows: Complete positive adVNTR rows with State/RU/POS context.
        advntr_calls_by_row: Legacy adVNTR display calls aligned with those rows.
        translation_component: Explicit checked-in translation authority resolved at
            the stage boundary.
        known_variant_names: Configured known display names used only for tier gating.
        bam_translations: Optional complete-candidate bindings for BAM calls. A
            binding remains Kestrel-internal evidence and does not create another
            caller or representation.

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
    from vntyper.scripts.identity_candidates import RawRepresentationKey, capture_advntr_observations

    identity_columns = frozenset((*IDENTITY_CAPTURE_COLUMNS, *IDENTITY_SELECTION_COLUMNS))
    observed_columns = frozenset().union(*(frozenset(row.keys()) for row in kestrel_rows))
    present = identity_columns & observed_columns
    if not present:
        return None
    if len(kestrel_rows) != len(vcf_calls) or any(not identity_columns.issubset(row.keys()) for row in kestrel_rows):
        raise ValueError("Kestrel identity metadata is incomplete or misaligned")
    if len(advntr_rows) != len(advntr_calls_by_row):
        raise ValueError("adVNTR identity context is misaligned with its display calls")
    if bam_translations is not None and len(bam_translations) != len(bam_calls):
        raise ValueError("BAM identity bindings are misaligned with BAM display calls")

    observations: list[IdentityReconciliationObservation] = []
    presentation_call_index = 0
    for row, call in zip(kestrel_rows, vcf_calls, strict=True):
        persisted = parse_selected_candidate_cells(row)
        row_key = RawRepresentationKey(
            "kestrel",
            (
                _required_nonempty_string(row.get("Motifs"), "Kestrel Motifs"),
                _strict_positive_integer(row.get("POS"), "Kestrel POS"),
                _required_nonempty_string(row.get("REF"), "Kestrel REF"),
                _required_nonempty_string(row.get("ALT"), "Kestrel ALT"),
            ),
        )
        if persisted.selected_row_key != row_key:
            raise ValueError("Kestrel persisted selected key does not match the actual result row")
        kmer_depth = _strict_nonnegative_integer(
            row.get("Estimated_Depth_AlternateVariant"),
            "Kestrel alternate k-mer-path depth",
        )
        if call is None:
            continue
        translation = persisted.translation
        if translation.identity is not None and translation.context_diverges != persisted.group_context_diverges:
            translation = IdentityTranslation(translation.identity, "resolved", None, persisted.group_context_diverges)
        observations.append(
            _make_observation(
                call,
                translation,
                known_variant_names,
                kmer_depth=kmer_depth,
                extra_flags=persisted.flags,
                blocking_gates=persisted.blocking_gates,
                presentation_call_index=presentation_call_index,
            )
        )
        presentation_call_index += 1

    aligned_bam_translations = bam_translations or tuple(None for _ in bam_calls)
    for call, bound_translation in zip(bam_calls, aligned_bam_translations, strict=True):
        if call is not None:
            observations.append(
                _make_observation(
                    call,
                    bound_translation or IdentityTranslation(None, "unresolved", "missing-motif-context", False),
                    known_variant_names,
                    presentation_call_index=presentation_call_index,
                )
            )
            presentation_call_index += 1
        elif bound_translation is not None:
            raise ValueError("BAM identity binding requires an aligned BAM display call")

    capture_rows: list[dict[str, object]] = []
    for row, calls in zip(advntr_rows, advntr_calls_by_row, strict=True):
        support = _strict_nonnegative_integer(
            row.get("NumberOfSupportingReads"),
            "adVNTR sequencing-read support",
        )
        coverage = _strict_nonnegative_number(row.get("MeanCoverage"), "adVNTR MeanCoverage")
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
        support = _strict_nonnegative_integer(
            row.get("NumberOfSupportingReads"),
            "adVNTR sequencing-read support",
        )
        coverage = _strict_nonnegative_number(row.get("MeanCoverage"), "adVNTR MeanCoverage")
        if len(calls) == 1:
            observations.append(
                _make_observation(
                    calls[0],
                    candidate.observation.translation,
                    known_variant_names,
                    advntr_reads=support,
                    advntr_coverage=coverage,
                    extra_flags=candidate.flags,
                    presentation_call_index=presentation_call_index,
                )
            )
            presentation_call_index += 1
            continue
        for call in calls:
            observations.append(
                _make_observation(
                    call,
                    IdentityTranslation(None, "unresolved", "reconstruction-mismatch", False),
                    known_variant_names,
                    disposition=EvidenceDisposition("admissible"),
                    extra_flags=candidate.flags,
                    presentation_call_index=presentation_call_index,
                )
            )
            presentation_call_index += 1
        observations.append(
            IdentityReconciliationObservation(
                translation=candidate.observation.translation,
                display_name=None,
                source="advntr",
                event=_compound_event(calls),
                net_length=candidate.observation.frame_consequence.net_length_change,
                flags=candidate.flags,
                disposition=EvidenceDisposition("admissible"),
                known_variant_match=False,
                advntr_sequencing_read_support=support,
                advntr_mean_coverage=coverage,
            )
        )
    return tuple(observations)


def reconcile_identity_observations(
    observations: tuple[IdentityReconciliationObservation, ...],
    policy: IdentityReconciliationPolicy,
    *,
    presentation_selected_observation_index: int | None = None,
) -> IdentityReconciliationResult:
    """Reconcile caller evidence by resolved canonical identity.

    Args:
        observations: Legacy-ordered observations. Kestrel VCF remains the explicit
            primary presentation when no unique corroborated identity outvotes it.
        policy: Separately named thresholds in each caller's evidence unit.
        presentation_selected_observation_index: Explicit observation selected by
            the unchanged presentation vote. An unresolved or inadmissible selection
            remains visible but produces an abstained identity decision.

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
    if presentation_selected_observation_index is not None and (
        isinstance(presentation_selected_observation_index, bool)
        or not isinstance(presentation_selected_observation_index, int)
        or presentation_selected_observation_index < 0
        or presentation_selected_observation_index >= len(observations)
    ):
        raise ValueError("Presentation-selected observation index must identify an observation")

    decision_indices = tuple(
        index for index, observation in enumerate(observations) if observation.disposition.value == "admissible"
    )
    presentation_event_indices = tuple(
        index for index in decision_indices if observations[index].presentation_call_index is not None
    )
    event_indices = presentation_event_indices or decision_indices
    event_disagreement = len({observations[index].event for index in event_indices}) > 1
    backing: dict[MolecularIdentity, list[int]] = {}
    for index in decision_indices:
        observation = observations[index]
        if observation.backs_identity:
            assert observation.identity is not None
            backing.setdefault(observation.identity, []).append(index)

    corroborated = [identity for identity, indices in backing.items() if len(_callers(observations, indices)) >= 2]

    selected_observation = (
        observations[presentation_selected_observation_index]
        if presentation_selected_observation_index is not None
        else None
    )
    if len(corroborated) == 1:
        chosen_identity = corroborated[0]
        bound_indices = [
            index for index in backing[chosen_identity] if observations[index].presentation_call_index is not None
        ]
        selected_index = (
            min(bound_indices, key=lambda index: observations[index].presentation_call_index or 0)
            if bound_indices
            else backing[chosen_identity][0]
        )
        presentation_result_index = selected_index
    elif not backing or (selected_observation is not None and not selected_observation.backs_identity):
        return IdentityReconciliationResult(
            decision=IdentityDecision(None, "abstained", False, "identity-unresolved"),
            selected_observation_index=None,
            backing_observation_indices=(),
            backing_sources=(),
            backing_callers=(),
            caller_disagreement=len(backing) > 1 or event_disagreement,
            event_disagreement=event_disagreement,
            low_support_sources=(),
            presentation_selected_observation_index=presentation_selected_observation_index,
        )

    elif selected_observation is not None:
        selected_identity = selected_observation.identity
        assert selected_identity is not None
        chosen_identity = selected_identity
        assert presentation_selected_observation_index is not None
        selected_index = presentation_selected_observation_index
        presentation_result_index = selected_index
    else:
        primary_index = next(
            (
                index
                for index, observation in enumerate(observations)
                if observation.source == "kestrel_vcf" and observation.backs_identity
            ),
            next(iter(backing.values()))[0],
        )
        primary_identity = observations[primary_index].identity
        assert primary_identity is not None
        chosen_identity = primary_identity
        selected_index = backing[chosen_identity][0]
        presentation_result_index = selected_index

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
    flags = frozenset().union(*(observations[index].flags for index in decision_indices))
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
        presentation_selected_observation_index=presentation_result_index,
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
    advntr_coverage: int | float | None = None,
    disposition: EvidenceDisposition = _ADMISSIBLE_DISPOSITION,
    extra_flags: frozenset[str] = frozenset(),
    blocking_gates: frozenset[str] = frozenset(),
    presentation_call_index: int | None = None,
) -> IdentityReconciliationObservation:
    return IdentityReconciliationObservation(
        translation=translation,
        display_name=call.name,
        source=call.source,
        event=call.event,
        net_length=call.net_length,
        flags=frozenset(call.flags) | extra_flags,
        disposition=disposition,
        known_variant_match=call.name in known_variant_names,
        kestrel_alternate_kmer_path_depth=kmer_depth,
        advntr_sequencing_read_support=advntr_reads,
        advntr_mean_coverage=advntr_coverage,
        blocking_gates=blocking_gates,
        presentation_call_index=presentation_call_index,
    )


def _strict_nonnegative_integer(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a canonical finite non-negative integer")
    if isinstance(value, int):
        if value >= 0:
            return value
        raise ValueError(f"{name} must be a canonical finite non-negative integer")
    if isinstance(value, str) and _CANONICAL_NONNEGATIVE_INTEGER.fullmatch(value):
        return int(value)
    raise ValueError(f"{name} must be a canonical finite non-negative integer")


def _strict_positive_integer(value: object, name: str) -> int:
    parsed = _strict_nonnegative_integer(value, name)
    if parsed < 1:
        raise ValueError(f"{name} must be a canonical positive integer")
    return parsed


def _strict_nonnegative_number(value: object, name: str) -> int | float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite non-negative number")
    if isinstance(value, (int, float)):
        if math.isfinite(value) and value >= 0:
            return value
        raise ValueError(f"{name} must be a finite non-negative number")
    if isinstance(value, str) and _CANONICAL_NONNEGATIVE_NUMBER.fullmatch(value):
        parsed = float(value) if "." in value else int(value)
        if math.isfinite(parsed):
            return parsed
    raise ValueError(f"{name} must be a finite non-negative number")


def _required_nonempty_string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a non-empty string")
    return value


def _compound_event(calls: Sequence[ReconciliationPresentationCall]) -> str:
    events = {call.event for call in calls}
    return next(iter(events)) if len(events) == 1 else "compound"


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


def _optional_nonnegative_number(value: object, name: str) -> None:
    if value is not None and (
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0
    ):
        raise ValueError(f"{name} must be a finite non-negative number or None")
