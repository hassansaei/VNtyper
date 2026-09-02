"""Runtime projection of generated dominance onto fixed caller candidates."""

from importlib import import_module

import pytest

from vntyper.scripts.identity_reconciliation import IdentityReconciliationObservation
from vntyper.scripts.molecular_identity import (
    EvidenceDisposition,
    IdentityTranslation,
    make_coding_edit,
    make_molecular_identity,
)
from vntyper.scripts.nomenclature import Nomenclature
from vntyper.scripts.nomenclature_bam_evidence import BamIdentityEvidence, BamLocusEvidence
from vntyper.scripts.nomenclature_bam_replay import BamReplayArtifact, BamReplayLocus
from vntyper.scripts.nomenclature_identity_projection import IdentityAwareNomenclatureResult

pytestmark = pytest.mark.unit

_DUPC = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
_DUPG = make_molecular_identity((make_coding_edit(59, 58, "", "G"),))


def _call(name: str, *, tier: str = "B", source: str = "kestrel_vcf") -> Nomenclature:
    return Nomenclature(name, "duplication", "X", tier, (), None, None, 1, source)


def _record(identity, xd: int | None = 10) -> BamIdentityEvidence:
    ordinal = 0 if identity == _DUPC else 1
    return BamIdentityEvidence((identity,), ((ordinal,),), xd)


def _artifact(*records: BamIdentityEvidence) -> BamReplayArtifact:
    counts = {
        identity: sum(identity in record.identities for record in records)
        for record in records
        for identity in record.identities
    }
    evidence = BamLocusEvidence(tuple(records), len(records), counts)
    return BamReplayArtifact((BamReplayLocus((0, 1), "observed", evidence),))


def _component(**changes: object) -> dict[str, object]:
    component: dict[str, object] = {
        "enabled": True,
        "minimum_record_count_margin": 1,
        "minimum_record_share": 0.5,
        "minimum_record_share_margin": 0.0,
        "xd_veto": "disabled",
        "abstain_on_inadmissible_advntr": False,
    }
    component.update(changes)
    return component


def _runtime():
    return import_module("vntyper.scripts.nomenclature_dominance_runtime")


def test_dominance_selects_only_an_existing_canonical_candidate_and_keeps_its_fixed_projection() -> None:
    """Mutation caught: the BAM display name is parsed into a fabricated delins identity."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC", tier="A"), _DUPC)
    candidates = (
        runtime.DominanceCandidate(_DUPC, _call("59dupC", tier="A")),
        runtime.DominanceCandidate(_DUPG, _call("59_60delinsAT", tier="A", source="kestrel_bam")),
    )

    result = runtime.reconcile_with_dominance(
        baseline,
        candidates,
        _artifact(_record(_DUPC), _record(_DUPG), _record(_DUPG)),
        EvidenceDisposition("admissible"),
        _component(),
        {0: _DUPC, 1: _DUPG},
    )

    assert result.selected_identity == _DUPG
    assert result.call == candidates[1].call
    assert result.call.name == "59_60delinsAT"
    assert result.call.tier == "A"
    assert result.abstention_reason is None


@pytest.mark.parametrize(
    ("component", "disposition", "reason"),
    [
        (_component(minimum_record_count_margin=3), EvidenceDisposition("admissible"), "insufficient-dominance"),
        (
            _component(abstain_on_inadmissible_advntr=True),
            EvidenceDisposition("identity-insufficient"),
            "inadmissible-advntr",
        ),
    ],
)
def test_dominance_abstention_suppresses_the_whole_locus_selection(
    component: dict[str, object], disposition: EvidenceDisposition, reason: str
) -> None:
    """Mutation caught: a veto falls through to the winner or preserves the baseline call."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC"), _DUPC)

    result = runtime.reconcile_with_dominance(
        baseline,
        (runtime.DominanceCandidate(_DUPC, baseline.call),),
        _artifact(_record(_DUPC), _record(_DUPC)),
        disposition,
        component,
        {0: _DUPC, 1: _DUPG},
    )

    assert result.call is None
    assert result.selected_identity is None
    assert result.abstention_reason == reason
    assert runtime.dominance_abstention_note(reason) == f"Whole-locus dominance abstention: {reason}."


def test_record_tie_abstains_even_when_raw_xd_favours_a_candidate() -> None:
    """Mutation caught: raw XD breaks the record tie or falls through to one identity."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC"), _DUPC)
    candidates = (
        runtime.DominanceCandidate(_DUPC, baseline.call),
        runtime.DominanceCandidate(_DUPG, _call("58dupG")),
    )

    result = runtime.reconcile_with_dominance(
        baseline,
        candidates,
        _artifact(_record(_DUPC, 1), _record(_DUPG, 10_000)),
        EvidenceDisposition("admissible"),
        _component(),
        {0: _DUPC, 1: _DUPG},
    )

    assert result.call is None
    assert result.selected_identity is None
    assert result.abstention_reason == "record-tie"


def test_selected_identity_absent_from_fixed_candidates_fails_closed() -> None:
    """Mutation caught: a dominance identity is rendered by guessing from display text."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC"), _DUPC)

    with pytest.raises(ValueError, match="absent from.*canonical candidate"):
        runtime.reconcile_with_dominance(
            baseline,
            (runtime.DominanceCandidate(_DUPC, baseline.call),),
            _artifact(_record(_DUPG), _record(_DUPG)),
            EvidenceDisposition("admissible"),
            _component(),
            {0: _DUPC, 1: _DUPG},
        )


def test_replay_identity_binding_must_match_the_authoritative_candidate_ordinal() -> None:
    """Mutation caught: runtime trusts a replay identity bound to another candidate."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC"), _DUPC)
    impossible = BamIdentityEvidence((_DUPG,), ((0,),), 10)

    with pytest.raises(ValueError, match="BAM.*binding.*candidate|candidate.*identity"):
        runtime.reconcile_with_dominance(
            baseline,
            (
                runtime.DominanceCandidate(_DUPC, baseline.call),
                runtime.DominanceCandidate(_DUPG, _call("58dupG")),
            ),
            _artifact(impossible, impossible),
            EvidenceDisposition("admissible"),
            _component(),
            {0: _DUPC, 1: _DUPG},
        )


def test_fixed_candidates_keep_row_rendering_but_use_the_fixed_whole_locus_tier() -> None:
    """Mutation caught: dominance itself assigns Tier A or drops a BAM-refined name."""
    runtime = _runtime()
    baseline = _call("59dupC", tier="A")
    bam_refined = _call("59_60delinsAT", tier="B", source="kestrel_bam")

    candidates = runtime.fixed_kestrel_candidate_projections((_DUPG,), (bam_refined,), baseline)

    assert candidates == (runtime.DominanceCandidate(_DUPG, candidates[0].call),)
    assert candidates[0].call.name == "59_60delinsAT"
    assert candidates[0].call.tier == "A"


def test_fixed_candidates_include_a_distinct_typed_advntr_projection() -> None:
    """Mutation caught: candidate closure includes only persisted Kestrel row identities."""
    runtime = _runtime()
    baseline = _call("59dupC", tier="B")
    advntr_call = _call("58dupG", tier="C", source="advntr")
    observation = IdentityReconciliationObservation(
        IdentityTranslation(_DUPG, "resolved", None, False),
        advntr_call.name,
        "advntr",
        advntr_call.event,
        advntr_call.net_length,
        frozenset(),
        EvidenceDisposition("admissible"),
        False,
        advntr_sequencing_read_support=10,
        advntr_mean_coverage=20,
        presentation_call_index=0,
    )

    candidates = runtime.fixed_caller_candidate_projections(
        (_DUPC,),
        (_call("59dupC", tier="C"),),
        (observation,),
        (advntr_call,),
        baseline,
    )

    assert tuple(candidate.identity for candidate in candidates) == (_DUPC, _DUPG)
    assert candidates[1].call.name == "58dupG"
    assert candidates[1].call.tier == "B"


def test_inadmissible_advntr_identity_is_masked_from_candidate_promotion() -> None:
    """Mutation caught: a governed adVNTR identity becomes a selectable candidate."""
    runtime = _runtime()
    advntr_call = _call("58dupG", source="advntr")
    observation = IdentityReconciliationObservation(
        IdentityTranslation(_DUPG, "resolved", None, False),
        advntr_call.name,
        "advntr",
        advntr_call.event,
        advntr_call.net_length,
        frozenset(),
        EvidenceDisposition("identity-insufficient"),
        False,
        advntr_sequencing_read_support=10,
        advntr_mean_coverage=20,
        presentation_call_index=0,
    )

    candidates = runtime.fixed_caller_candidate_projections(
        (_DUPC,),
        (_call("59dupC"),),
        (observation,),
        (advntr_call,),
        _call("59dupC"),
    )

    assert tuple(candidate.identity for candidate in candidates) == (_DUPC,)


def test_observed_loci_are_merged_as_unweighted_records_and_absence_stays_missing() -> None:
    """Mutation caught: only one locus is used, or unavailable BAM becomes zero support."""
    runtime = _runtime()
    first = BamLocusEvidence((_record(_DUPC),), 1, {_DUPC: 1})
    second_records = (
        BamIdentityEvidence((_DUPG,), ((1,),), 10),
        BamIdentityEvidence((_DUPG,), ((1,),), 10),
    )
    second = BamLocusEvidence(second_records, 2, {_DUPG: 2})
    replay = BamReplayArtifact(
        (
            BamReplayLocus((0,), "observed", first),
            BamReplayLocus((1,), "observed", second),
        )
    )

    merged = runtime.retained_whole_locus_bam_evidence(replay, {0: _DUPC, 1: _DUPG})

    assert merged is not None
    assert merged.eligible_record_count == 3
    assert merged.counts == {_DUPC: 1, _DUPG: 2}
    assert (
        runtime.retained_whole_locus_bam_evidence(
            BamReplayArtifact((BamReplayLocus((0,), "unavailable", None),)),
            {0: _DUPC},
        )
        is None
    )


def test_replay_retention_keeps_prior_observed_evidence_when_reobservation_is_unavailable() -> None:
    """Mutation caught: reconciliation discards the already-retained BAM evidence."""
    runtime = _runtime()
    existing = _artifact(_record(_DUPC), _record(_DUPC))

    retained = runtime.retain_bam_replay(
        existing,
        (BamReplayLocus((0, 1), "unavailable", None),),
        (0, 1),
    )

    assert retained == existing


def test_disabled_policy_is_byte_path_neutral_and_evaluates_once(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mutation caught: neutral dominance rewrites selection or evaluates more than once."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC", tier="A"), _DUPC)
    calls = 0
    real_evaluate = runtime.evaluate_dominance

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return real_evaluate(*args, **kwargs)

    monkeypatch.setattr(runtime, "evaluate_dominance", counted)
    component = _component(enabled=False)
    result = runtime.reconcile_with_dominance(
        baseline,
        (runtime.DominanceCandidate(_DUPC, baseline.call),),
        _artifact(_record(_DUPC), _record(_DUPC)),
        EvidenceDisposition("admissible"),
        component,
        {0: _DUPC, 1: _DUPG},
    )

    assert calls == 1
    assert result.call is baseline.call
    assert result.selected_identity is baseline.selected_identity
    assert result.abstention_reason is None


def test_disabled_policy_preserves_an_identity_abstained_legacy_baseline() -> None:
    """Mutation caught: result validation rejects a valid call with no selected identity."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC", tier="B"), None)

    result = runtime.reconcile_with_dominance(
        baseline,
        (),
        _artifact(_record(_DUPC), _record(_DUPC)),
        EvidenceDisposition("admissible"),
        _component(enabled=False),
        {0: _DUPC, 1: _DUPG},
    )

    assert result.call is baseline.call
    assert result.selected_identity is None
    assert result.decision.outcome == "not-applicable"
    assert result.abstention_reason is None


def test_absent_bam_evaluation_calls_the_real_evaluator_once(monkeypatch: pytest.MonkeyPatch) -> None:
    runtime = _runtime()
    calls = 0
    real_evaluate = runtime.evaluate_dominance

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return real_evaluate(*args, **kwargs)

    monkeypatch.setattr(runtime, "evaluate_dominance", counted)

    decision = runtime.evaluate_absent_bam_dominance(
        (EvidenceDisposition("admissible"),),
        _component(),
    )

    assert calls == 1
    assert decision.outcome == "not-applicable"


def test_absent_bam_evaluation_rejects_malformed_dispositions() -> None:
    runtime = _runtime()

    with pytest.raises(ValueError, match="dispositions"):
        runtime.evaluate_absent_bam_dominance([], _component())


def test_absent_bam_evaluation_fails_closed_on_an_impossible_evaluator_outcome(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runtime = _runtime()
    monkeypatch.setattr(
        runtime,
        "evaluate_dominance",
        lambda *_args, **_kwargs: runtime.DominanceDecision("abstained", None, "record-tie"),
    )

    with pytest.raises(ValueError, match="must evaluate as not-applicable"):
        runtime.evaluate_absent_bam_dominance((), _component())


def test_xd_missingness_veto_abstains_at_runtime_instead_of_promoting_the_covered_runner_up() -> None:
    """Mutation caught: a missing-XD veto on the winner falls through to the fully covered runner-up."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC"), _DUPC)
    candidates = (
        runtime.DominanceCandidate(_DUPC, baseline.call),
        runtime.DominanceCandidate(_DUPG, _call("58dupG", source="kestrel_bam")),
    )

    result = runtime.reconcile_with_dominance(
        baseline,
        candidates,
        _artifact(_record(_DUPC, xd=None), _record(_DUPC), _record(_DUPG, xd=99)),
        EvidenceDisposition("admissible"),
        _component(xd_veto="missingness"),
        {0: _DUPC, 1: _DUPG},
    )

    assert result.call is None
    assert result.selected_identity is None
    assert result.abstention_reason == "xd-missingness"


def test_xd_missingness_veto_never_promotes_a_candidate_record_counts_did_not_select() -> None:
    """Mutation caught: complete XD coverage is treated as a reason to prefer the runner-up."""
    runtime = _runtime()
    baseline = IdentityAwareNomenclatureResult(_call("59dupC"), _DUPC)
    candidates = (
        runtime.DominanceCandidate(_DUPC, baseline.call),
        runtime.DominanceCandidate(_DUPG, _call("58dupG", source="kestrel_bam")),
    )

    result = runtime.reconcile_with_dominance(
        baseline,
        candidates,
        _artifact(_record(_DUPC), _record(_DUPC), _record(_DUPG, xd=None)),
        EvidenceDisposition("admissible"),
        _component(xd_veto="missingness"),
        {0: _DUPC, 1: _DUPG},
    )

    assert result.selected_identity == _DUPC
    assert result.call == candidates[0].call
    assert result.abstention_reason is None
