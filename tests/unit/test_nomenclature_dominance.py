"""Closed dominance and whole-locus abstention decisions."""

from fractions import Fraction

import pytest

from vntyper.scripts.molecular_identity import EvidenceDisposition, make_coding_edit, make_molecular_identity
from vntyper.scripts.nomenclature_bam_evidence import BamIdentityEvidence, BamLocusEvidence
from vntyper.scripts.nomenclature_dominance import DominanceEvidence, evaluate_dominance

pytestmark = pytest.mark.unit

_DUPC = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
_DUPA = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))


def _record(*identities, xd: int | None = 5) -> BamIdentityEvidence:
    return BamIdentityEvidence(tuple(identities), tuple((index,) for index, _ in enumerate(identities)), xd)


def _evidence(
    records: tuple[BamIdentityEvidence, ...],
    *,
    disposition: str = "admissible",
) -> DominanceEvidence:
    counts = {
        identity: sum(identity in record.identities for record in records)
        for identity in (_DUPC, _DUPA)
        if any(identity in record.identities for record in records)
    }
    return DominanceEvidence(
        BamLocusEvidence(records, len(records), counts),
        EvidenceDisposition(disposition),  # type: ignore[arg-type]
    )


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


def test_disabled_or_missing_evidence_is_not_applicable() -> None:
    disabled = evaluate_dominance(_evidence((_record(_DUPC),)), _component(enabled=False))
    missing = evaluate_dominance(DominanceEvidence(None, EvidenceDisposition("admissible")), _component())
    empty = evaluate_dominance(_evidence(()), _component())

    assert (disabled.outcome, disabled.identity, disabled.abstention_reason) == ("not-applicable", None, None)
    assert (missing.outcome, missing.identity, missing.abstention_reason) == ("not-applicable", None, None)
    assert (empty.outcome, empty.identity, empty.abstention_reason) == ("not-applicable", None, None)


def test_record_tie_abstains_before_unequal_xd_can_choose_a_winner() -> None:
    evidence = _evidence(
        (
            _record(_DUPC, xd=1),
            _record(_DUPC, xd=1),
            _record(_DUPA, xd=10_000),
            _record(_DUPA, xd=10_000),
        )
    )

    decision = evaluate_dominance(evidence, _component(xd_veto="discordance"))

    assert decision.outcome == "abstained"
    assert decision.abstention_reason == "record-tie"
    assert decision.identity is None


def test_counts_and_shares_use_complete_record_denominator() -> None:
    evidence = _evidence(
        (
            _record(_DUPC, _DUPA),
            _record(_DUPC),
            _record(_DUPC),
            _record(_DUPA),
            _record(),
        )
    )

    decision = evaluate_dominance(evidence, _component())

    assert decision.outcome == "selected"
    assert decision.identity == _DUPC
    assert decision.top_record_count == 3
    assert decision.runner_up_record_count == 2
    assert decision.record_count_margin == 1
    assert decision.top_record_share == Fraction(3, 5)
    assert decision.runner_up_record_share == Fraction(2, 5)
    assert decision.record_share_margin == Fraction(1, 5)
    assert evidence.bam_evidence is not None
    assert sum((Fraction(count, 5) for count in evidence.bam_evidence.counts.values()), start=Fraction(0)) == Fraction(
        1, 1
    )


def test_cooccurring_identity_shares_may_sum_above_one_without_duplicate_votes() -> None:
    evidence = _evidence(
        (
            _record(_DUPC, _DUPA),
            _record(_DUPC, _DUPA),
            _record(_DUPC),
        )
    )

    decision = evaluate_dominance(evidence, _component())

    assert decision.top_record_count == 3
    assert decision.runner_up_record_count == 2
    assert decision.top_record_share + decision.runner_up_record_share == Fraction(5, 3)


@pytest.mark.parametrize(
    ("change", "records"),
    [
        ({"minimum_record_count_margin": 2}, (_record(_DUPC), _record(_DUPC), _record(_DUPA))),
        ({"minimum_record_share": 0.75}, (_record(_DUPC), _record(_DUPC), _record(_DUPA))),
        ({"minimum_record_share_margin": 0.5}, (_record(_DUPC), _record(_DUPC), _record(_DUPA))),
    ],
)
def test_each_dominance_floor_can_abstain(change: dict[str, object], records: tuple[BamIdentityEvidence, ...]) -> None:
    decision = evaluate_dominance(_evidence(records), _component(**change))

    assert (decision.outcome, decision.identity, decision.abstention_reason) == (
        "abstained",
        None,
        "insufficient-dominance",
    )


def test_governed_advntr_disposition_can_only_veto_the_count_winner() -> None:
    evidence = _evidence(
        (_record(_DUPC), _record(_DUPC), _record(_DUPA)),
        disposition="identity-insufficient",
    )

    ignored = evaluate_dominance(evidence, _component(abstain_on_inadmissible_advntr=False))
    vetoed = evaluate_dominance(evidence, _component(abstain_on_inadmissible_advntr=True))

    assert ignored.outcome == "selected" and ignored.identity == _DUPC
    assert (vetoed.outcome, vetoed.identity, vetoed.abstention_reason) == (
        "abstained",
        None,
        "inadmissible-advntr",
    )


def test_xd_missingness_veto_abstains_instead_of_selecting_the_runner_up() -> None:
    evidence = _evidence((_record(_DUPC, xd=None), _record(_DUPC, xd=5), _record(_DUPA, xd=99)))

    decision = evaluate_dominance(evidence, _component(xd_veto="missingness"))

    assert (decision.outcome, decision.identity, decision.abstention_reason) == (
        "abstained",
        None,
        "xd-missingness",
    )


def test_xd_concentration_veto_detects_one_record_dominating_winner_support() -> None:
    concentrated = _evidence(
        (_record(_DUPC, xd=100), _record(_DUPC, xd=1), _record(_DUPC, xd=1), _record(_DUPA, xd=30))
    )
    balanced = _evidence((_record(_DUPC, xd=10), _record(_DUPC, xd=10), _record(_DUPC, xd=10), _record(_DUPA, xd=30)))

    vetoed = evaluate_dominance(concentrated, _component(xd_veto="concentration"))
    selected = evaluate_dominance(balanced, _component(xd_veto="concentration"))

    assert vetoed.abstention_reason == "xd-concentration"
    assert selected.outcome == "selected" and selected.identity == _DUPC


def test_xd_discordance_veto_never_promotes_the_runner_up() -> None:
    evidence = _evidence((_record(_DUPC, xd=1), _record(_DUPC, xd=1), _record(_DUPA, xd=10)))

    decision = evaluate_dominance(evidence, _component(xd_veto="discordance"))

    assert (decision.outcome, decision.identity, decision.abstention_reason) == (
        "abstained",
        None,
        "xd-discordance",
    )


@pytest.mark.parametrize(
    "component",
    [
        {},
        _component(extra=True),
        _component(enabled=1),
        _component(minimum_record_count_margin=-1),
        _component(minimum_record_share=1.1),
        _component(minimum_record_share_margin=True),
        _component(xd_veto="winner"),
        _component(abstain_on_inadmissible_advntr=1),
    ],
)
def test_dominance_component_is_strict_and_typed(component: dict[str, object]) -> None:
    with pytest.raises(ValueError, match="dominance"):
        evaluate_dominance(_evidence((_record(_DUPC),)), component)


def test_dominance_evidence_and_decision_are_closed_typed_values() -> None:
    with pytest.raises(ValueError, match="BamLocusEvidence"):
        DominanceEvidence(object(), EvidenceDisposition("admissible"))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="EvidenceDisposition"):
        DominanceEvidence(None, object())  # type: ignore[arg-type]
