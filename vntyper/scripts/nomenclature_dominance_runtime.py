"""Project one retained-evidence dominance decision onto fixed candidates."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any, Literal

from vntyper.scripts.identity_reconciliation import (
    IdentityReconciliationObservation,
    IdentityReconciliationPolicy,
    build_identity_reconciliation_observations,
)
from vntyper.scripts.molecular_identity import EvidenceDisposition, IdentityTranslation, MolecularIdentity
from vntyper.scripts.nomenclature import Nomenclature
from vntyper.scripts.nomenclature_bam_evidence import BamLocusEvidence
from vntyper.scripts.nomenclature_bam_replay import (
    BamReplayArtifact,
    BamReplayLocus,
    merge_bam_replay_artifacts,
    validate_bam_replay_artifact_ordinals,
)
from vntyper.scripts.nomenclature_bam_universe import validated_whole_locus_bam_evidence
from vntyper.scripts.nomenclature_dominance import (
    AbstentionReason,
    DominanceDecision,
    DominanceEvidence,
    evaluate_dominance,
)
from vntyper.scripts.nomenclature_identity_projection import IdentityAwareNomenclatureResult

if TYPE_CHECKING:
    from vntyper.modules.advntr.artifact_evidence import ArtifactEvidence
    from vntyper.scripts.nomenclature_decision_config import NomenclatureDecisionConfig

DOMINANCE_ABSTENTION_COLUMN = "__Dominance_Abstention_Reason"


@dataclass(frozen=True)
class DominanceSeamOutcome:
    """Explicit result returned only by an enabled production dominance seam."""

    evaluated: bool
    rewritten: bool
    dominance_outcome: Literal["selected", "abstained", "not-applicable"]

    def __post_init__(self) -> None:
        """Require enabled seams to report a complete actual evaluation."""
        if self.evaluated is not True:
            raise ValueError("enabled dominance seam outcome must record evaluation")
        if not isinstance(self.rewritten, bool):
            raise ValueError("enabled dominance seam rewrite status must be Boolean")
        if self.dominance_outcome not in {"selected", "abstained", "not-applicable"}:
            raise ValueError("enabled dominance seam must report a valid dominance outcome")


def _is_nomenclature(value: object) -> bool:
    """Recognize the frozen value across the suite's deliberate module reload."""
    value_type = type(value)
    return value_type.__module__ == Nomenclature.__module__ and value_type.__name__ == Nomenclature.__name__


def rows_carry_identity_metadata(rows: Sequence[Mapping[str, object]]) -> bool:
    """Return whether result rows carry any current-run identity metadata.

    Args:
        rows: Positive Kestrel result rows.

    Returns:
        Whether at least one internal capture or selection column is present.
    """
    from vntyper.scripts.identity_candidate_persistence import IDENTITY_CAPTURE_COLUMNS, IDENTITY_SELECTION_COLUMNS

    identity_columns = frozenset((*IDENTITY_CAPTURE_COLUMNS, *IDENTITY_SELECTION_COLUMNS))
    return any(column in row for row in rows for column in identity_columns)


def production_identity_observations(
    kestrel_rows: Sequence[Mapping[str, object]],
    vcf_calls: Sequence[Nomenclature | None],
    bam_calls: Sequence[Nomenclature | None],
    advntr_rows: Sequence[Mapping[str, object]],
    advntr_calls_by_row: Sequence[Sequence[Nomenclature]],
    bam_translations: Sequence[IdentityTranslation | None] | None = None,
    *,
    decision_config: NomenclatureDecisionConfig,
    translation_component: Any = None,
    artifact_evidence: ArtifactEvidence,
) -> tuple[tuple[IdentityReconciliationObservation, ...], IdentityReconciliationPolicy] | None:
    """Build typed production observations or retain the deliberate legacy path.

    Args:
        kestrel_rows: Positive Kestrel rows with optional identity metadata.
        vcf_calls: Kestrel VCF presentations aligned to result rows.
        bam_calls: Kestrel BAM refinements aligned to result rows.
        advntr_rows: Positive complete adVNTR rows.
        advntr_calls_by_row: adVNTR presentations aligned to positive rows.
        bam_translations: Optional typed bindings for BAM refinements.
        decision_config: Resolved immutable nomenclature policy.
        translation_component: Resolved identity translation authority.
        artifact_evidence: Verified governed adVNTR evidence.

    Returns:
        Typed observations and policy, or ``None`` for legacy rows.

    Raises:
        ValueError: If current-run identity metadata lacks translation authority
            or observation construction unexpectedly selects the legacy path.
    """
    if not rows_carry_identity_metadata(kestrel_rows):
        return None
    if translation_component is None:
        raise ValueError("Current-run identity metadata requires an explicit translation component")
    observations = build_identity_reconciliation_observations(
        kestrel_rows,
        vcf_calls,
        bam_calls,
        advntr_rows,
        advntr_calls_by_row,
        translation_component,
        frozenset(decision_config.known_variants),
        artifact_evidence=artifact_evidence,
        bam_translations=bam_translations,
    )
    if observations is None:  # pragma: no cover - identity columns were established above
        raise ValueError("Current-run identity metadata unexpectedly selected the legacy path")
    return observations, decision_config.identity_reconciliation


@dataclass(frozen=True)
class DominanceCandidate:
    """One existing canonical candidate and its already-fixed presentation."""

    identity: MolecularIdentity
    call: Nomenclature

    def __post_init__(self) -> None:
        """Reject untyped identities and presentation projections."""
        if not isinstance(self.identity, MolecularIdentity):
            raise ValueError("dominance candidate identity must be a MolecularIdentity")
        if not _is_nomenclature(self.call):
            raise ValueError("dominance candidate presentation must be a Nomenclature")


@dataclass(frozen=True)
class DominanceNomenclatureResult:
    """A fixed selected call, explicit abstention, or unchanged baseline."""

    call: Nomenclature | None
    selected_identity: MolecularIdentity | None
    abstention_reason: AbstentionReason | None
    decision: DominanceDecision

    def __post_init__(self) -> None:
        """Keep selected and abstained runtime projections disjoint."""
        if not isinstance(self.decision, DominanceDecision):
            raise ValueError("runtime dominance result requires a DominanceDecision")
        if self.call is not None and not _is_nomenclature(self.call):
            raise ValueError("runtime dominance call must be a Nomenclature or None")
        if self.selected_identity is not None and not isinstance(self.selected_identity, MolecularIdentity):
            raise ValueError("runtime dominance identity must be a MolecularIdentity or None")
        if self.abstention_reason is not None:
            if self.call is not None or self.selected_identity is not None:
                raise ValueError("runtime dominance abstention cannot retain a selected call or identity")
            if self.decision.outcome != "abstained" or self.decision.abstention_reason != self.abstention_reason:
                raise ValueError("runtime dominance abstention must match its decision")
        elif self.decision.outcome == "not-applicable":
            if self.call is None:
                raise ValueError("not-applicable runtime dominance must preserve the baseline call")
        elif self.call is None or self.selected_identity is None:
            raise ValueError("runtime dominance selection requires a call and canonical identity")


def fixed_kestrel_candidate_projections(
    identities: tuple[MolecularIdentity | None, ...],
    calls: tuple[Nomenclature | None, ...],
    whole_locus_call: Nomenclature,
) -> tuple[DominanceCandidate, ...]:
    """Bind canonical row identities to their fixed rendered name and locus tier.

    Args:
        identities: Persisted canonical identities aligned to Kestrel result rows.
        calls: Already-refined row presentation calls in the same order.
        whole_locus_call: Fixed reconciliation projection that owns tier and flags.

    Returns:
        Existing canonical candidates in row order. Missing identities or calls do
        not become candidates.

    Raises:
        ValueError: If inputs are untyped or row alignment differs.
    """
    if not isinstance(identities, tuple) or any(
        identity is not None and not isinstance(identity, MolecularIdentity) for identity in identities
    ):
        raise ValueError("fixed dominance identities must be a tuple of MolecularIdentity or None values")
    if not isinstance(calls, tuple) or any(call is not None and not _is_nomenclature(call) for call in calls):
        raise ValueError("fixed dominance calls must be a tuple of Nomenclature or None values")
    if len(identities) != len(calls):
        raise ValueError("fixed dominance identities and calls must stay row-aligned")
    if not _is_nomenclature(whole_locus_call):
        raise ValueError("fixed dominance whole-locus projection must be a Nomenclature")
    return tuple(
        DominanceCandidate(identity, replace(call, tier=whole_locus_call.tier, flags=whole_locus_call.flags))
        for identity, call in zip(identities, calls, strict=True)
        if identity is not None and call is not None
    )


def fixed_caller_candidate_projections(
    kestrel_identities: tuple[MolecularIdentity | None, ...],
    kestrel_calls: tuple[Nomenclature | None, ...],
    identity_observations: tuple[IdentityReconciliationObservation, ...],
    presentation_calls: tuple[Nomenclature, ...],
    whole_locus_call: Nomenclature,
) -> tuple[DominanceCandidate, ...]:
    """Close candidates over all typed caller observations and fixed displays.

    Kestrel row projections come first because they already incorporate any BAM
    presentation refinement. Remaining admissible typed VCF, BAM, and adVNTR observations
    can add a distinct canonical identity, but never replace a fixed projection
    and never derive identity from display text.

    Args:
        kestrel_identities: Persisted canonical identities aligned to Kestrel rows.
        kestrel_calls: Refined fixed Kestrel row presentation calls.
        identity_observations: Typed reconciliation observations for every caller.
        presentation_calls: Legacy calls bound by observation presentation indices.
        whole_locus_call: Fixed reconciliation projection that owns tier and flags.

    Returns:
        One deterministic fixed projection per existing canonical identity.

    Raises:
        ValueError: If observations or presentation bindings are malformed.
    """
    if not isinstance(identity_observations, tuple) or any(
        not isinstance(observation, IdentityReconciliationObservation) for observation in identity_observations
    ):
        raise ValueError("fixed dominance observations must be typed identity observations")
    if not isinstance(presentation_calls, tuple) or any(not _is_nomenclature(call) for call in presentation_calls):
        raise ValueError("fixed dominance presentation calls must be a tuple of Nomenclature values")
    candidates = list(fixed_kestrel_candidate_projections(kestrel_identities, kestrel_calls, whole_locus_call))
    existing = {candidate.identity for candidate in candidates}
    for observation in identity_observations:
        identity = observation.identity
        call_index = observation.presentation_call_index
        if not observation.backs_identity or identity is None or call_index is None or identity in existing:
            continue
        if call_index >= len(presentation_calls):
            raise ValueError("fixed dominance observation refers to an absent presentation call")
        call = presentation_calls[call_index]
        candidates.append(
            DominanceCandidate(identity, replace(call, tier=whole_locus_call.tier, flags=whole_locus_call.flags))
        )
        existing.add(identity)
    return tuple(candidates)


def retained_whole_locus_bam_evidence(
    artifact: BamReplayArtifact,
    authoritative_identities: Mapping[int, MolecularIdentity | None],
) -> BamLocusEvidence | None:
    """Merge complete observed replay loci without fabricating unavailable records.

    Args:
        artifact: Retained BAM replay artifact after lifecycle merging.
        authoritative_identities: Persisted candidate identities keyed by ordinal.

    Returns:
        All observed records in locus order, or ``None`` when BAM evidence is absent.

    Raises:
        ValueError: If ``artifact`` is not a validated replay artifact.
    """
    return validated_whole_locus_bam_evidence(artifact, authoritative_identities)


def evaluate_absent_bam_dominance(
    advntr_dispositions: tuple[EvidenceDisposition, ...],
    component: Mapping[str, object],
) -> DominanceDecision:
    """Evaluate one genuine no-BAM/no-locus seam through the real policy evaluator.

    Args:
        advntr_dispositions: Governed dispositions for visible positive adVNTR rows.
        component: Complete immutable dominance policy.

    Returns:
        The evaluator's required not-applicable decision.

    Raises:
        ValueError: If dispositions are malformed or absent BAM does not evaluate as
            not-applicable.
    """
    if not isinstance(advntr_dispositions, tuple) or any(
        not isinstance(disposition, EvidenceDisposition) for disposition in advntr_dispositions
    ):
        raise ValueError("absent-BAM dominance dispositions must be a tuple of EvidenceDisposition values")
    whole_locus_disposition = EvidenceDisposition(
        "identity-insufficient"
        if any(disposition.value == "identity-insufficient" for disposition in advntr_dispositions)
        else "admissible"
    )
    decision = evaluate_dominance(DominanceEvidence(None, whole_locus_disposition), component)
    if decision.outcome != "not-applicable":
        raise ValueError("absent BAM dominance evidence must evaluate as not-applicable")
    return decision


def retain_bam_replay(
    existing: BamReplayArtifact | None,
    current_loci: tuple[BamReplayLocus, ...],
    expected_ordinals: tuple[int, ...],
) -> BamReplayArtifact:
    """Validate current replay custody and retain prior observed evidence.

    Args:
        existing: Previously persisted replay, if present.
        current_loci: Current observation attempt for each selected locus.
        expected_ordinals: Selected Kestrel observation ordinals.

    Returns:
        Validated current replay merged with retained evidence when available.

    Raises:
        ValueError: If inputs are untyped or replay custody does not match the
            selected Kestrel candidates.
    """
    if existing is not None and not isinstance(existing, BamReplayArtifact):
        raise ValueError("existing dominance replay must be a BamReplayArtifact or None")
    if not isinstance(current_loci, tuple) or any(not isinstance(locus, BamReplayLocus) for locus in current_loci):
        raise ValueError("current dominance replay loci must be a tuple of BamReplayLocus values")
    if not isinstance(expected_ordinals, tuple) or any(not isinstance(ordinal, int) for ordinal in expected_ordinals):
        raise ValueError("dominance replay ordinals must be a tuple of integers")
    current = BamReplayArtifact(current_loci)
    validate_bam_replay_artifact_ordinals(current, expected_ordinals)
    return merge_bam_replay_artifacts(existing, current) if existing is not None else current


def reconcile_with_dominance(
    baseline: IdentityAwareNomenclatureResult,
    candidates: tuple[DominanceCandidate, ...],
    replay_artifact: BamReplayArtifact,
    advntr_disposition: EvidenceDisposition,
    component: Mapping[str, object],
    authoritative_identities: Mapping[int, MolecularIdentity | None],
) -> DominanceNomenclatureResult:
    """Evaluate dominance once and select only an existing candidate projection.

    Args:
        baseline: Existing fixed whole-locus reconciliation result.
        candidates: Canonical candidates paired with fixed name/tier projections.
        replay_artifact: Already-retained BAM evidence; no BAM is opened here.
        advntr_disposition: Governed whole-locus adVNTR disposition. Absence is
            represented by admissible, not by fabricated insufficiency.
        component: Complete immutable dominance profile component.
        authoritative_identities: Captured candidate identity keyed by ordinal.

    Returns:
        Unchanged baseline, one fixed candidate projection, or explicit abstention.

    Raises:
        ValueError: If inputs are untyped, candidate projections conflict, or a
            selected identity is absent from the existing candidate set.
    """
    if not isinstance(baseline, IdentityAwareNomenclatureResult):
        raise ValueError("runtime dominance requires an identity-aware baseline")
    if not isinstance(candidates, tuple) or any(
        not isinstance(candidate, DominanceCandidate) for candidate in candidates
    ):
        raise ValueError("runtime dominance candidates must be a tuple of DominanceCandidate values")
    if not isinstance(advntr_disposition, EvidenceDisposition):
        raise ValueError("runtime dominance requires a governed adVNTR disposition")
    projections: dict[MolecularIdentity, Nomenclature] = {}
    for candidate in candidates:
        previous = projections.setdefault(candidate.identity, candidate.call)
        if (previous.name, previous.tier) != (candidate.call.name, candidate.call.tier):
            raise ValueError("one canonical dominance candidate has conflicting fixed presentations")

    evidence = DominanceEvidence(
        retained_whole_locus_bam_evidence(replay_artifact, authoritative_identities),
        advntr_disposition,
    )
    decision = evaluate_dominance(evidence, component)
    if decision.outcome == "not-applicable":
        return DominanceNomenclatureResult(baseline.call, baseline.selected_identity, None, decision)
    if decision.outcome == "abstained":
        return DominanceNomenclatureResult(None, None, decision.abstention_reason, decision)
    assert decision.identity is not None
    selected = projections.get(decision.identity)
    if selected is None:
        raise ValueError("dominance-selected identity is absent from the existing canonical candidate projections")
    if decision.identity == baseline.selected_identity:
        selected = baseline.call
    return DominanceNomenclatureResult(selected, decision.identity, None, decision)


def reconcile_retained_dominance(
    baseline: IdentityAwareNomenclatureResult,
    identities: tuple[MolecularIdentity | None, ...],
    calls: tuple[Nomenclature | None, ...],
    identity_observations: tuple[IdentityReconciliationObservation, ...],
    presentation_calls: tuple[Nomenclature, ...],
    replay_artifact: BamReplayArtifact,
    advntr_dispositions: tuple[EvidenceDisposition, ...],
    component: Mapping[str, object],
    authoritative_identities: Mapping[int, MolecularIdentity | None],
) -> DominanceNomenclatureResult:
    """Evaluate retained whole-locus evidence against fixed Kestrel projections.

    Args:
        baseline: Existing identity-aware reconciliation result.
        identities: Persisted canonical Kestrel row identities.
        calls: Already-refined Kestrel row presentations.
        identity_observations: Typed observations over every caller candidate.
        presentation_calls: Fixed calls bound by those observations.
        replay_artifact: Retained BAM replay for those rows.
        advntr_dispositions: Governed dispositions for positive adVNTR rows.
        component: Complete immutable dominance profile component.
        authoritative_identities: Captured candidate identity keyed by ordinal.

    Returns:
        One dominance runtime projection from exactly one evaluator call.

    Raises:
        ValueError: If the dispositions are untyped or candidate inputs fail
            validation.
    """
    if not isinstance(advntr_dispositions, tuple) or any(
        not isinstance(disposition, EvidenceDisposition) for disposition in advntr_dispositions
    ):
        raise ValueError("runtime dominance adVNTR dispositions must be a tuple of EvidenceDisposition values")
    candidates = fixed_caller_candidate_projections(
        identities,
        calls,
        identity_observations,
        presentation_calls,
        baseline.call,
    )
    whole_locus_disposition = EvidenceDisposition(
        "identity-insufficient"
        if any(disposition.value == "identity-insufficient" for disposition in advntr_dispositions)
        else "admissible"
    )
    return reconcile_with_dominance(
        baseline,
        candidates,
        replay_artifact,
        whole_locus_disposition,
        component,
        authoritative_identities,
    )


def reconcile_complete_retained_dominance(
    baseline: IdentityAwareNomenclatureResult,
    candidate_rows: tuple[Mapping[str, object], ...],
    candidate_calls: tuple[Nomenclature | None, ...],
    identity_observations: tuple[IdentityReconciliationObservation, ...],
    presentation_calls: tuple[Nomenclature, ...],
    replay_artifact: BamReplayArtifact,
    advntr_dispositions: tuple[EvidenceDisposition, ...],
    component: Mapping[str, object],
) -> DominanceNomenclatureResult:
    """Evaluate complete persisted Kestrel candidates with their fixed names.

    Args:
        baseline: Legacy whole-locus identity-aware reconciliation.
        candidate_rows: Passing pre-result rows with authoritative identity metadata.
        candidate_calls: Fixed row calls aligned to ``candidate_rows``.
        identity_observations: Typed observations from all callers.
        presentation_calls: Fixed legacy caller presentations.
        replay_artifact: Producer-retained whole-locus BAM evidence.
        advntr_dispositions: Governed positive adVNTR dispositions.
        component: Complete immutable dominance policy.

    Returns:
        One dominance result selected only from the closed typed candidate set.

    Raises:
        ValueError: If candidate identity or fixed-presentation custody is incomplete.
    """
    from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells

    if len(candidate_rows) != len(candidate_calls) or not candidate_rows:
        raise ValueError("enabled dominance requires aligned complete Kestrel candidate projections")
    persisted = tuple(parse_selected_candidate_cells(row) for row in candidate_rows)
    authoritative = {candidate.selected_observation_ordinal: candidate.translation.identity for candidate in persisted}
    fixed_calls: list[Nomenclature | None] = []
    for row, call in zip(candidate_rows, candidate_calls, strict=True):
        fixed_name = str(row.get("Nomenclature", ""))
        if call is None or not fixed_name:
            raise ValueError("enabled dominance requires one fixed presentation per complete Kestrel candidate")
        fixed_calls.append(replace(call, name=fixed_name))
    return reconcile_retained_dominance(
        baseline,
        tuple(candidate.translation.identity for candidate in persisted),
        tuple(fixed_calls),
        identity_observations,
        presentation_calls,
        replay_artifact,
        advntr_dispositions,
        component,
        authoritative,
    )


def dominance_abstention_note(reason: AbstentionReason) -> str:
    """Return the technical output note for one closed abstention token.

    Args:
        reason: Closed dominance abstention reason.

    Returns:
        Stable human-readable technical note.
    """
    return f"Whole-locus dominance abstention: {reason}."
