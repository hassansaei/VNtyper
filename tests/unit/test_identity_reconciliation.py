"""Molecular-identity keyed reconciliation policy contracts."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from typing import Any

import pytest

from vntyper.scripts.identity_candidates import translation_component_from_config
from vntyper.scripts.identity_reconciliation import (
    IdentityReconciliationObservation,
    IdentityReconciliationPolicy,
    IdentityReconciliationResult,
    build_identity_reconciliation_observations,
    reconcile_identity_observations,
)
from vntyper.scripts.molecular_identity import (
    EvidenceDisposition,
    IdentityDecision,
    IdentityTranslation,
    make_coding_edit,
    make_molecular_identity,
    serialize_molecular_identity,
)
from vntyper.scripts.nomenclature import KNOWN_VARIANTS, Nomenclature, from_advntr, load_nomenclature_config, reconcile

pytestmark = pytest.mark.unit

POLICY = IdentityReconciliationPolicy(
    kestrel_min_alternate_kmer_path_depth=5,
    advntr_min_sequencing_read_support=5,
)
DUPC = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
DUPA = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))
DUPG = make_molecular_identity((make_coding_edit(60, 59, "", "G"),))
TRANSLATION_COMPONENT = translation_component_from_config(load_nomenclature_config())


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
    presentation_call_index: int | None = None,
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
        presentation_call_index=presentation_call_index,
    )


def _presentation_call(
    name: str,
    source: str,
    *,
    event: str = "duplication",
    net_length: int = 1,
) -> Nomenclature:
    return Nomenclature(name, event, "X", "B", (), None, None, net_length, source)


def _persisted_kestrel_row(
    *,
    identity: Any = DUPC,
    selected_raw_key: str = '{"source":"kestrel","values":["X-X",67,"G","GG"]}',
    alternate_depth: object = "40",
) -> dict[str, object]:
    return {
        "Motifs": "X-X",
        "POS": "67",
        "REF": "G",
        "ALT": "GG",
        "Estimated_Depth_AlternateVariant": alternate_depth,
        "__Identity_Raw_Representation_Key": selected_raw_key,
        "__Identity_Molecular_Identity": serialize_molecular_identity(identity),
        "__Identity_Translation_Status": "resolved",
        "__Identity_Translation_Failure": "absent",
        "__Identity_Context_Diverges": "false",
        "__Identity_Observation_Ordinal": "0",
        "__Identity_Selected_Raw_Representation_Key": selected_raw_key,
        "__Identity_Equivalent_Representation_Count": "1",
        "__Identity_Hypothesis_Count": "1",
        "__Identity_Group_Blocking_Gates": "[]",
        "__Identity_Group_Flags": "[]",
        "__Identity_Selected_Observation_Ordinal": "0",
        "__Identity_Group_Context_Diverges": "false",
    }


def _build_observations(
    kestrel_row: dict[str, object],
    kestrel_call: Nomenclature | None,
    *,
    advntr_row: dict[str, object] | None = None,
    advntr_calls: list[Nomenclature] | None = None,
) -> tuple[IdentityReconciliationObservation, ...] | None:
    return build_identity_reconciliation_observations(
        [kestrel_row],
        [kestrel_call],
        [],
        [] if advntr_row is None else [advntr_row],
        [] if advntr_calls is None else [advntr_calls],
        TRANSLATION_COMPONENT,
        frozenset(KNOWN_VARIANTS),
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


def test_unique_molecular_winner_precedes_a_competing_legacy_vcf_selection() -> None:
    observations = (
        _observation(DUPG, "competing-first", "future_a", presentation_call_index=0),
        _observation(DUPA, "vcf-a", "kestrel_vcf", kmer_depth=40, presentation_call_index=1),
        _observation(DUPC, "first-b-display", "future_b", presentation_call_index=2),
        _observation(DUPC, "second-b-display", "advntr", advntr_reads=40, presentation_call_index=3),
    )

    result = reconcile_identity_observations(
        observations,
        POLICY,
        presentation_selected_observation_index=1,
    )

    assert result.identity == DUPC
    assert result.selected_observation_index == 2
    assert result.presentation_selected_observation_index == 2
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


def test_explicit_compatibility_selection_keeps_later_vcf_in_a_corroborated_tie() -> None:
    observations = (
        _observation(DUPA, "60dupA", "future_one"),
        _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
        _observation(DUPA, "60dupA", "future_two"),
        _observation(DUPC, "different-display", "advntr", advntr_reads=40),
    )

    result = reconcile_identity_observations(
        observations,
        POLICY,
        presentation_selected_observation_index=1,
    )

    assert result.identity == DUPC
    assert result.selected_observation_index == 1
    assert result.caller_disagreement is True
    assert result.tier == "B"


def test_unbound_bam_presentation_cannot_hide_a_different_selected_identity() -> None:
    observations = (
        _observation(DUPC, "59_60insG", "kestrel_vcf", event="insertion", kmer_depth=40),
        _observation(None, "58_59insG", "kestrel_bam", event="insertion"),
        _observation(DUPA, "58_59insG", "advntr", event="insertion", advntr_reads=40),
    )

    result = reconcile_identity_observations(
        observations,
        POLICY,
        presentation_selected_observation_index=1,
    )

    assert result.identity is None
    assert result.tier == "abstained"
    assert result.selected_observation_index is None
    assert result.presentation_selected_observation_index == 1


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


def test_bound_bam_translation_stays_kestrel_internal_without_changing_candidate_counts() -> None:
    row = _persisted_kestrel_row()
    before = (
        row["__Identity_Equivalent_Representation_Count"],
        row["__Identity_Hypothesis_Count"],
    )
    vcf_call = _presentation_call("59dupC", "kestrel_vcf")
    bam_call = _presentation_call("59dupC", "kestrel_bam")

    observations = build_identity_reconciliation_observations(
        [row],
        [vcf_call],
        [bam_call],
        [],
        [],
        TRANSLATION_COMPONENT,
        frozenset(KNOWN_VARIANTS),
        bam_translations=[IdentityTranslation(DUPC, "resolved", None, False)],
    )

    assert observations is not None
    assert [observation.identity for observation in observations] == [DUPC, DUPC]
    decision = reconcile_identity_observations(observations, POLICY)
    assert decision.backing_callers == ("kestrel",)
    assert decision.molecular_agreement is False
    assert decision.tier != "A"
    assert (
        (
            row["__Identity_Equivalent_Representation_Count"],
            row["__Identity_Hypothesis_Count"],
        )
        == before
        == ("1", "1")
    )


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


def test_result_rejects_event_disagreement_without_caller_disagreement() -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", event="duplication", kmer_depth=40),
            _observation(DUPC, "59dupC", "advntr", event="insertion", advntr_reads=40),
        ),
        POLICY,
    )

    with pytest.raises(ValueError, match="Event disagreement requires caller disagreement"):
        replace(result, caller_disagreement=False)


@pytest.mark.parametrize(
    "changes",
    [
        {"decision": IdentityDecision(DUPC, "C", True, None)},
        {
            "decision": IdentityDecision(DUPC, "B", True, None),
            "event_disagreement": True,
            "caller_disagreement": True,
        },
    ],
)
def test_result_rejects_tier_c_event_disagreement_drift(changes: dict[str, object]) -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
            _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
        ),
        POLICY,
    )

    with pytest.raises(ValueError, match="Tier C must match event disagreement"):
        replace(result, **changes)


@pytest.mark.parametrize(
    "changes",
    [
        {"caller_disagreement": True},
        {"caller_disagreement": True, "event_disagreement": True},
        {"low_support_sources": ("kestrel_vcf",)},
    ],
)
def test_result_rejects_tier_a_conflict_or_low_support_metadata(changes: dict[str, object]) -> None:
    result = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
            _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
        ),
        POLICY,
    )
    assert result.tier == "A"

    with pytest.raises(ValueError, match="Tier A cannot carry conflict or low-support metadata"):
        replace(result, **changes)


def test_result_accepts_valid_tier_a_b_and_c_boundaries() -> None:
    tier_a = reconcile_identity_observations(
        (
            _observation(DUPC, "59dupC", "kestrel_vcf", kmer_depth=40),
            _observation(DUPC, "59dupC", "advntr", advntr_reads=40),
        ),
        POLICY,
    )
    tier_b = replace(tier_a, decision=IdentityDecision(DUPC, "B", True, None))
    tier_c = replace(
        tier_a,
        decision=IdentityDecision(DUPC, "C", True, None),
        caller_disagreement=True,
        event_disagreement=True,
    )

    assert (tier_a.tier, tier_b.tier, tier_c.tier) == ("A", "B", "C")


def test_builder_rejects_canonical_persisted_key_for_a_different_kestrel_row() -> None:
    mismatched = '{"source":"kestrel","values":["X-X",68,"G","GG"]}'

    with pytest.raises(ValueError, match="persisted selected key does not match"):
        _build_observations(
            _persisted_kestrel_row(selected_raw_key=mismatched),
            _presentation_call("59dupC", "kestrel_vcf"),
        )


def test_builder_parses_malformed_kestrel_metadata_even_without_a_presentation_call() -> None:
    malformed = _persisted_kestrel_row()
    malformed["__Identity_Selected_Raw_Representation_Key"] = '{"source":"kestrel","values":[]}'

    with pytest.raises(ValueError, match="Kestrel raw representation key"):
        _build_observations(malformed, None)


@pytest.mark.parametrize(
    "value",
    [True, 1.5, "1.0", "1e2", "01", " 1", "1 ", "NaN", "inf", "-1"],
)
def test_kestrel_alternate_depth_rejects_noncanonical_integer_evidence(value: object) -> None:
    with pytest.raises(ValueError, match="alternate k-mer-path depth"):
        _build_observations(
            _persisted_kestrel_row(alternate_depth=value),
            _presentation_call("59dupC", "kestrel_vcf"),
        )


@pytest.mark.parametrize("value", [0, 40, "0", "40"])
def test_kestrel_alternate_depth_accepts_typed_or_canonical_wire_integers(value: object) -> None:
    observations = _build_observations(
        _persisted_kestrel_row(alternate_depth=value),
        _presentation_call("59dupC", "kestrel_vcf"),
    )

    assert observations is not None
    assert observations[0].kestrel_alternate_kmer_path_depth == int(value)


@pytest.mark.parametrize(
    "value",
    [True, 1.5, "1.0", "1e2", "01", " 1", "1 ", "NaN", "inf", "-1"],
)
def test_advntr_read_support_rejects_noncanonical_integer_evidence(value: object) -> None:
    state = "I22_2_G_LEN1"
    with pytest.raises(ValueError, match="sequencing-read support"):
        _build_observations(
            _persisted_kestrel_row(),
            _presentation_call("59dupC", "kestrel_vcf"),
            advntr_row={
                "Variant": state,
                "RU": "2",
                "POS": "22",
                "NumberOfSupportingReads": value,
                "MeanCoverage": "153.98",
                "Flag": "Not flagged",
            },
            advntr_calls=from_advntr(state),
        )


@pytest.mark.parametrize("value", [True, -1, -1.0, float("nan"), float("inf"), "NaN", "inf", "-1"])
def test_advntr_mean_coverage_rejects_nonfinite_or_negative_numeric_evidence(value: object) -> None:
    state = "I22_2_G_LEN1"
    with pytest.raises(ValueError, match="MeanCoverage"):
        _build_observations(
            _persisted_kestrel_row(),
            _presentation_call("59dupC", "kestrel_vcf"),
            advntr_row={
                "Variant": state,
                "RU": "2",
                "POS": "22",
                "NumberOfSupportingReads": "40",
                "MeanCoverage": value,
                "Flag": "Not flagged",
            },
            advntr_calls=from_advntr(state),
        )


def test_fractional_advntr_mean_coverage_survives_the_capture_boundary() -> None:
    state = "I22_2_G_LEN1"
    observations = _build_observations(
        _persisted_kestrel_row(),
        _presentation_call("59dupC", "kestrel_vcf"),
        advntr_row={
            "Variant": state,
            "RU": "2",
            "POS": "22",
            "NumberOfSupportingReads": "40",
            "MeanCoverage": "153.98",
            "Flag": "Not flagged",
        },
        advntr_calls=from_advntr(state),
    )

    assert observations is not None
    advntr_observation = next(observation for observation in observations if observation.source == "advntr")
    assert advntr_observation.advntr_mean_coverage == 153.98


@pytest.mark.parametrize(
    ("state", "positions", "component_names"),
    [
        ("I22_2_G_LEN1&I23_2_A_LEN1", "22,23", ("59dupC", "58_59insT")),
        ("I22_2_G_LEN1&I24_2_A_LEN1", "22,24", ("59dupC", "57_58insT")),
    ],
)
def test_compound_advntr_row_has_one_row_identity_and_presentation_only_components(
    state: str,
    positions: str,
    component_names: tuple[str, str],
) -> None:
    calls = from_advntr(state)
    observations = _build_observations(
        _persisted_kestrel_row(),
        _presentation_call("59dupC", "kestrel_vcf"),
        advntr_row={
            "Variant": state,
            "RU": "2,2",
            "POS": positions,
            "NumberOfSupportingReads": "40",
            "MeanCoverage": "153.98",
            "Flag": "Not flagged",
        },
        advntr_calls=calls,
    )

    assert observations is not None
    advntr_observations = [observation for observation in observations if observation.source == "advntr"]
    row_observations = [
        observation for observation in advntr_observations if observation.presentation_call_index is None
    ]
    component_observations = [
        observation for observation in advntr_observations if observation.presentation_call_index is not None
    ]
    assert len(row_observations) == 1
    assert row_observations[0].identity is not None
    assert len(row_observations[0].identity.edits) == 2
    assert sum(len(edit.inserted) - len(edit.deleted) for edit in row_observations[0].identity.edits) == 2
    assert row_observations[0].disposition == EvidenceDisposition("admissible")
    assert row_observations[0].known_variant_match is False
    assert tuple(observation.display_name for observation in component_observations) == component_names
    assert all(observation.identity is None for observation in component_observations)
    assert all(observation.disposition == EvidenceDisposition("admissible") for observation in component_observations)
    assert all(observation.advntr_sequencing_read_support is None for observation in component_observations)
    assert all(observation.advntr_mean_coverage is None for observation in component_observations)
    advntr_only = reconcile_identity_observations(tuple(advntr_observations), POLICY)
    assert len(advntr_only.backing_observation_indices) == 1
    assert advntr_only.backing_sources == ("advntr",)
    assert advntr_only.backing_callers == ("advntr",)
    assert advntr_only.molecular_agreement is False

    kestrel_complete = replace(
        _observation(
            row_observations[0].identity,
            "complete-kestrel-projection",
            "kestrel_vcf",
            event="insertion",
            kmer_depth=40,
            known=False,
        ),
        presentation_call_index=0,
    )
    cross_caller = reconcile_identity_observations(
        (kestrel_complete, row_observations[0]),
        POLICY,
        presentation_selected_observation_index=0,
    )
    assert cross_caller.molecular_agreement is True
    assert cross_caller.event_disagreement is False
    assert cross_caller.tier == "B"


def test_noncontiguous_compound_keeps_exact_pr_a_legacy_projection_without_component_backing() -> None:
    state = "I22_2_G_LEN1&I24_2_A_LEN1"
    kestrel_call = _presentation_call("59dupC", "kestrel_vcf")
    component_calls = from_advntr(state)
    observations = _build_observations(
        _persisted_kestrel_row(),
        kestrel_call,
        advntr_row={
            "Variant": state,
            "RU": "2,2",
            "POS": "22,24",
            "NumberOfSupportingReads": "40",
            "MeanCoverage": "153.98",
            "Flag": "Not flagged",
        },
        advntr_calls=component_calls,
    )

    assert observations is not None
    typed = reconcile(
        kestrel_call,
        *component_calls,
        identity_observations=observations,
        identity_policy=POLICY,
    )
    legacy = reconcile(kestrel_call, *component_calls)
    advntr_observations = [observation for observation in observations if observation.source == "advntr"]

    assert typed == legacy
    assert typed.name is None
    assert typed.event == "unknown"
    assert typed.tier == "C"
    assert len([observation for observation in advntr_observations if observation.backs_identity]) == 1
    assert all(
        observation.advntr_sequencing_read_support is None
        for observation in advntr_observations
        if observation.presentation_call_index is not None
    )
