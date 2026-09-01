"""Molecular-identity keyed reconciliation policy contracts."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from typing import Any

import pytest

from vntyper.scripts.identity_reconciliation import (
    IdentityReconciliationObservation,
    IdentityReconciliationPolicy,
    IdentityReconciliationResult,
    reconcile_identity_observations,
)
from vntyper.scripts.molecular_identity import (
    EvidenceDisposition,
    IdentityDecision,
    IdentityTranslation,
    make_coding_edit,
    make_molecular_identity,
)

pytestmark = pytest.mark.unit

POLICY = IdentityReconciliationPolicy(
    kestrel_min_alternate_kmer_path_depth=5,
    advntr_min_sequencing_read_support=5,
)
DUPC = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
DUPA = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))


def _observation(
    identity: Any,
    display_name: str | None,
    source: str,
    *,
    event: str = "duplication",
    disposition: str = "admissible",
    context_diverges: bool = False,
    kmer_depth: int | None = None,
    advntr_reads: int | None = None,
    known: bool = True,
    blocking_gates: frozenset[str] = frozenset(),
) -> IdentityReconciliationObservation:
    translation = (
        IdentityTranslation(identity, "resolved", None, context_diverges)
        if identity is not None
        else IdentityTranslation(None, "unresolved", "missing-motif-context", False)
    )
    return IdentityReconciliationObservation(
        translation=translation,
        display_name=display_name,
        source=source,
        event=event,
        net_length=1,
        flags=frozenset(),
        disposition=EvidenceDisposition(disposition),
        known_variant_match=known,
        kestrel_alternate_kmer_path_depth=kmer_depth,
        advntr_sequencing_read_support=advntr_reads,
        blocking_gates=blocking_gates,
    )


DISPLAY_EQUAL_IDENTITY_DIFFERENT = (
    _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
    _observation(DUPA, "59dupC", "advntr", advntr_reads=40),
)


def test_display_equal_identity_different_cannot_agree() -> None:
    decision = reconcile_identity_observations(DISPLAY_EQUAL_IDENTITY_DIFFERENT, POLICY)
    assert decision.molecular_agreement is False
    assert decision.tier != "A"
    assert decision.caller_disagreement is True


def test_display_different_identity_equal_agrees_without_replacing_legacy_display_row() -> None:
    observations = (
        _observation(DUPC, "legacy-vcf-display", "kestrel_vcf", kmer_depth=40),
        _observation(DUPC, "different-advntr-display", "advntr", advntr_reads=40),
    )

    result = reconcile_identity_observations(observations, POLICY)

    assert result.identity == DUPC
    assert result.molecular_agreement is True
    assert result.tier == "A"
    assert result.selected_observation_index == 0
    assert observations[result.selected_observation_index].display_name == "legacy-vcf-display"


def test_unresolved_observations_never_agree_through_display_event_or_position_words() -> None:
    observations = (
        _observation(None, "59dupC", "kestrel_vcf", kmer_depth=40),
        _observation(None, "59dupC", "advntr", advntr_reads=40),
    )

    result = reconcile_identity_observations(observations, POLICY)

    assert result.identity is None
    assert result.molecular_agreement is False
    assert result.tier == "abstained"
    assert result.backing_observation_indices == ()


def test_one_uniquely_corroborated_identity_retains_the_existing_vote() -> None:
    observations = (
        _observation(DUPA, "60dupA", "kestrel_vcf", kmer_depth=40),
        _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
        _observation(DUPC, "59dupC", "future_caller"),
    )

    result = reconcile_identity_observations(observations, POLICY)

    assert result.identity == DUPC
    assert result.selected_observation_index == 1
    assert result.molecular_agreement is True
    assert result.caller_disagreement is True


def test_two_corroborated_identities_are_conflict_with_vcf_primary_projection() -> None:
    observations = (
        _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
        _observation(DUPA, "60dupA", "future_one"),
        _observation(DUPA, "60dupA", "future_two"),
        _observation(DUPC, "different-display", "advntr", advntr_reads=40),
    )

    result = reconcile_identity_observations(observations, POLICY)

    assert result.identity == DUPC
    assert result.selected_observation_index == 0
    assert result.molecular_agreement is True
    assert result.caller_disagreement is True
    assert result.tier != "A"


def test_context_divergent_evidence_is_retained_but_blocks_tier_a() -> None:
    observations = (
        _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40, context_diverges=True),
        _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
    )

    result = reconcile_identity_observations(observations, POLICY)

    assert result.identity == DUPC
    assert result.molecular_agreement is True
    assert result.tier == "B"
    assert result.selected_observation_index == 0


def test_equivalent_group_blocking_gates_remain_tier_a_blockers() -> None:
    result = reconcile_identity_observations(
        (
            _observation(
                DUPC,
                "59dupC",
                "kestrel_vcf",
                kmer_depth=40,
                blocking_gates=frozenset({"motif_filter_pass"}),
            ),
            _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
        ),
        POLICY,
    )

    assert result.molecular_agreement is True
    assert result.tier == "B"


def test_kestrel_vcf_and_bam_are_one_caller_not_independent_corroboration() -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
            _observation(DUPC, "59dupC", "kestrel_bam"),
        ),
        POLICY,
    )

    assert result.molecular_agreement is False
    assert result.tier != "A"
    assert result.backing_callers == ("kestrel",)


@pytest.mark.parametrize(
    ("kmer_depth", "advntr_reads"),
    [(None, 40), (40, None), (4, 40), (40, 4)],
)
def test_tier_a_requires_known_source_bound_kmer_depth_and_advntr_reads(
    kmer_depth: int | None,
    advntr_reads: int | None,
) -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=kmer_depth),
            _observation(DUPC, "59dupC", "advntr", advntr_reads=advntr_reads),
        ),
        POLICY,
    )

    assert result.molecular_agreement is True
    assert result.tier == "B"


def test_non_backing_sources_cannot_lend_tier_a_support() -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=4),
            _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
            _observation(DUPA, "60dupA", "kestrel_vcf", kmer_depth=400),
        ),
        POLICY,
    )

    assert result.identity == DUPC
    assert result.tier == "B"
    assert result.low_support_sources == ("kestrel_vcf",)


def test_identity_insufficient_disposition_cannot_confer_agreement_or_tier_a() -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
            _observation(
                DUPC,
                "59dupC",
                "advntr",
                advntr_reads=40,
                disposition="identity-insufficient",
            ),
        ),
        POLICY,
    )

    assert result.identity == DUPC
    assert result.molecular_agreement is False
    assert result.tier != "A"
    assert result.backing_sources == ("kestrel_vcf",)


def test_identity_insufficient_observation_cannot_supply_the_selected_display_row() -> None:
    observations = (
        _observation(DUPC, "inadmissible-display", "advntr", advntr_reads=40, disposition="identity-insufficient"),
        _observation(DUPC, "vcf-display", "kestrel_vcf", kmer_depth=40),
        _observation(DUPC, "future-display", "future_caller"),
    )

    result = reconcile_identity_observations(observations, POLICY)

    assert result.selected_observation_index == 1


def test_unresolved_event_disagreement_is_retained_for_the_presentation_adapter() -> None:
    result = reconcile_identity_observations(
        (
            _observation(None, "59dupC", "kestrel_vcf", event="duplication", kmer_depth=40),
            _observation(None, "59dupC", "advntr", event="insertion", advntr_reads=40),
        ),
        POLICY,
    )

    assert result.event_disagreement is True
    assert result.caller_disagreement is True


def test_event_disagreement_keeps_identity_equality_separate_and_blocks_tier_a() -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", event="duplication", kmer_depth=40),
            _observation(DUPC, "59dupC", "advntr", event="insertion", advntr_reads=40),
        ),
        POLICY,
    )

    assert result.identity == DUPC
    assert result.molecular_agreement is True
    assert result.tier == "C"
    assert result.caller_disagreement is True


def test_policy_retains_two_separately_named_tier_a_thresholds() -> None:
    assert POLICY.kestrel_min_alternate_kmer_path_depth == 5
    assert POLICY.advntr_min_sequencing_read_support == 5


@pytest.mark.parametrize(
    "factory",
    [
        lambda: IdentityReconciliationPolicy(True, 5),
        lambda: IdentityReconciliationPolicy(5, 0),
        lambda: replace(_observation(DUPC, "59dupC", "kestrel_vcf"), source=""),
        lambda: replace(_observation(DUPC, "59dupC", "kestrel_vcf"), display_name=""),
        lambda: replace(_observation(DUPC, "59dupC", "kestrel_vcf"), net_length=True),
        lambda: replace(_observation(DUPC, "59dupC", "kestrel_vcf"), flags=("bad",)),
        lambda: replace(_observation(DUPC, "59dupC", "kestrel_vcf"), known_variant_match=1),
        lambda: replace(
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=5),
            advntr_sequencing_read_support=5,
        ),
        lambda: replace(_observation(DUPC, "59dupC", "advntr", advntr_reads=5), kestrel_alternate_kmer_path_depth=5),
        lambda: replace(_observation(DUPC, "59dupC", "kestrel_bam"), kestrel_alternate_kmer_path_depth=5),
    ],
)
def test_reconciliation_values_reject_invalid_states(factory: Callable[[], object]) -> None:
    with pytest.raises(ValueError):
        factory()


def test_reconciliation_requires_a_policy_and_typed_observation_tuple() -> None:
    with pytest.raises(ValueError):
        reconcile_identity_observations([], POLICY)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        reconcile_identity_observations((), None)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "changes",
    [
        {"selected_observation_index": None},
        {"backing_observation_indices": ()},
        {"backing_observation_indices": (0, 0)},
        {"backing_sources": ()},
        {"backing_callers": ()},
        {"low_support_sources": ("future_caller",)},
        {"decision": IdentityDecision(DUPC, "B", False, None)},
    ],
)
def test_result_rejects_inconsistent_selected_identity_metadata(changes: dict[str, object]) -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
            _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
        ),
        POLICY,
    )

    with pytest.raises(ValueError):
        replace(result, **changes)


def test_result_rejects_backing_metadata_on_an_abstained_decision() -> None:
    with pytest.raises(ValueError):
        IdentityReconciliationResult(
            decision=IdentityDecision(None, "abstained", False, "identity-unresolved"),
            selected_observation_index=None,
            backing_observation_indices=(),
            backing_sources=("advntr",),
            backing_callers=(),
            caller_disagreement=False,
            event_disagreement=False,
            low_support_sources=(),
        )
