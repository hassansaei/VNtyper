"""Canonical-identity binding for resolved Kestrel BAM haplotype records."""

from __future__ import annotations

from types import MappingProxyType

import pytest

from vntyper.scripts.identity_candidates import IdentityTranslationComponent
from vntyper.scripts.molecular_identity import IdentityTranslation, make_coding_edit, make_molecular_identity
from vntyper.scripts.nomenclature_bam_evidence import (
    BamCandidateBinding,
    BamEditObservation,
    BamIdentityEvidence,
    BamIdentityLocus,
    BamLocusEvidence,
    BamRecordObservation,
    bind_complete_winner_identity,
    collect_locus_evidence,
    project_record_observation,
)

pytestmark = pytest.mark.unit

_X = "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC"
_PAIR = _X + _X
_COMPONENT = IdentityTranslationComponent({"X": _X}, {"2": "X"}, 39)
_DUPC = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
_DUPA = make_molecular_identity((make_coding_edit(60, 59, "", "A"),))
_DEL_C = make_molecular_identity((make_coding_edit(30, 30, "C", ""),))
_INS_T = make_molecular_identity((make_coding_edit(31, 30, "", "T"),))
_DELINS_AT = make_molecular_identity((make_coding_edit(30, 30, "C", "AT"),))


def _translation(identity=_DUPC) -> IdentityTranslation:
    return IdentityTranslation(identity, "resolved", None, False)


def _locus(*identities) -> BamIdentityLocus:
    return BamIdentityLocus(
        motifs="X-X",
        pair_sequence=_PAIR,
        candidates=tuple(
            BamCandidateBinding(index, _translation(identity)) for index, identity in enumerate(identities)
        ),
    )


def _record(*edits: BamEditObservation, xd: int | None = None) -> BamRecordObservation:
    return BamRecordObservation(edits=tuple(edits), minimum_kmer_depth=xd)


_DUPC_EDIT = BamEditObservation(start=67, ref_span=0, inserted=1, bases="G")
_DUPA_EDIT = BamEditObservation(start=61, ref_span=0, inserted=1, bases="T")


def test_record_votes_once_for_each_cooccurring_identity() -> None:
    evidence = collect_locus_evidence(
        [_record(_DUPC_EDIT, _DUPC_EDIT, _DUPA_EDIT)],
        _locus(_DUPC, _DUPA),
        _COMPONENT,
    )

    assert evidence.counts == {_DUPC: 1, _DUPA: 1}
    assert evidence.eligible_record_count == 1
    assert evidence.record_identity_sets == ((_DUPC, _DUPA),)
    assert evidence.xd_by_record == (None,)
    assert evidence.bindings_by_record == (((0,), (1,)),)


def test_denominator_keeps_records_without_a_supported_candidate_identity() -> None:
    evidence = collect_locus_evidence(
        [_record(_DUPC_EDIT, xd=5), _record(xd=181), _record(_DUPA_EDIT, xd=8_704)],
        _locus(_DUPC),
        _COMPONENT,
    )

    assert evidence.counts == {_DUPC: 1}
    assert evidence.eligible_record_count == 3
    assert evidence.record_identity_sets == ((_DUPC,), (), ())
    assert evidence.xd_by_record == (5, 181, 8_704)


def test_record_identity_sets_and_runner_up_order_are_retained_for_replay() -> None:
    records = [
        _record(_DUPC_EDIT, _DUPA_EDIT, xd=5),
        _record(_DUPC_EDIT, xd=None),
        _record(_DUPA_EDIT, xd=8_704),
        _record(xd=0),
    ]

    evidence = collect_locus_evidence(records, _locus(_DUPC, _DUPA), _COMPONENT)

    assert evidence.counts == {_DUPC: 2, _DUPA: 2}
    assert evidence.record_identity_sets == ((_DUPC, _DUPA), (_DUPC,), (_DUPA,), ())
    assert evidence.xd_by_record == (5, None, 8_704, 0)


def test_three_records_outvote_two_regardless_of_inverted_extreme_xd() -> None:
    evidence = collect_locus_evidence(
        [
            _record(_DUPC_EDIT, xd=0),
            _record(_DUPC_EDIT, xd=None),
            _record(_DUPC_EDIT, xd=5),
            _record(_DUPA_EDIT, xd=2_147_483_647),
            _record(_DUPA_EDIT, xd=2_147_483_647),
        ],
        _locus(_DUPC, _DUPA),
        _COMPONENT,
    )

    assert evidence.counts == {_DUPC: 3, _DUPA: 2}
    assert evidence.winning_identity is _DUPC


def test_unequal_xd_cannot_break_an_equal_record_count_tie() -> None:
    evidence = collect_locus_evidence(
        [
            _record(_DUPC_EDIT, xd=5),
            _record(_DUPC_EDIT, xd=5),
            _record(_DUPA_EDIT, xd=8_704),
            _record(_DUPA_EDIT, xd=2_147_483_647),
        ],
        _locus(_DUPC, _DUPA),
        _COMPONENT,
    )

    assert evidence.counts == {_DUPC: 2, _DUPA: 2}
    assert evidence.winning_identity is None


def test_equivalent_kestrel_representations_are_bindings_not_extra_votes() -> None:
    locus = BamIdentityLocus(
        motifs="X-X",
        pair_sequence=_PAIR,
        candidates=(BamCandidateBinding(7, _translation()), BamCandidateBinding(11, _translation())),
    )

    evidence = collect_locus_evidence([_record(_DUPC_EDIT)], locus, _COMPONENT)

    assert evidence.counts == {_DUPC: 1}
    assert evidence.records[0].identities == (_DUPC,)
    assert evidence.records[0].candidate_observation_ordinals == ((7, 11),)


def test_unresolved_or_absent_kestrel_candidates_cannot_be_created_by_bam() -> None:
    unresolved = IdentityTranslation(None, "unresolved", "missing-motif-context", False)
    locus = BamIdentityLocus("X-X", _PAIR, (BamCandidateBinding(3, unresolved),))

    evidence = collect_locus_evidence([_record(_DUPC_EDIT)], locus, _COMPONENT)

    assert evidence.counts == {}
    assert evidence.record_identity_sets == ((),)


def test_pair_4092_shape_cannot_fabricate_an_absent_delins_candidate() -> None:
    # One resolved record carries a contiguous 1X1I delins at plus position 61,
    # but the only existing Kestrel candidates are its separately representable
    # deletion and insertion shapes. The compound identity is deliberately absent.
    delins_record = _record(BamEditObservation(start=30, ref_span=1, inserted=2, bases="AT"), xd=8_704)

    positive_control = collect_locus_evidence((delins_record,), _locus(_DELINS_AT), _COMPONENT)
    evidence = collect_locus_evidence((delins_record for _ in range(2)), _locus(_DEL_C, _INS_T), _COMPONENT)

    assert positive_control.counts == {_DELINS_AT: 1}
    assert _DELINS_AT not in evidence.counts
    assert evidence.counts == {}
    assert evidence.record_identity_sets == ((), ())


def test_separate_records_cannot_be_combined_into_one_delins_identity() -> None:
    evidence = collect_locus_evidence(
        [
            _record(BamEditObservation(start=30, ref_span=1, inserted=0, bases="")),
            _record(BamEditObservation(start=30, ref_span=0, inserted=1, bases="A")),
        ],
        _locus(_DELINS_AT),
        _COMPONENT,
    )

    assert evidence.counts == {}
    assert evidence.record_identity_sets == ((), ())


def test_identity_evidence_is_frozen_and_validates_aligned_bindings() -> None:
    evidence = BamIdentityEvidence((_DUPC,), ((4,),), 5)
    assert evidence.identities == (_DUPC,)

    with pytest.raises(ValueError, match="aligned"):
        BamIdentityEvidence((_DUPC,), (), 5)
    with pytest.raises(ValueError, match="unique"):
        BamIdentityEvidence((_DUPC, _DUPC), ((4,), (5,)), 5)
    with pytest.raises(ValueError, match="one identity"):
        BamIdentityEvidence((_DUPC, _DUPA), ((4,), (4,)), 5)
    with pytest.raises(ValueError, match="minimum k-mer depth"):
        BamIdentityEvidence((), (), -1)


@pytest.mark.parametrize("value", [True, 1.5, 2_147_483_648])
def test_record_and_bound_evidence_reject_invalid_xd_values(value: object) -> None:
    with pytest.raises(ValueError, match="minimum k-mer depth"):
        BamRecordObservation((), value)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="minimum k-mer depth"):
        BamIdentityEvidence((), (), value)  # type: ignore[arg-type]


@pytest.mark.parametrize("value", [0, 2_147_483_647])
def test_record_and_bound_evidence_accept_xd_boundaries(value: int) -> None:
    assert BamRecordObservation((), value).minimum_kmer_depth == value
    assert BamIdentityEvidence((), (), value).minimum_kmer_depth == value


def test_locus_evidence_rejects_counts_inconsistent_with_records() -> None:
    record = BamIdentityEvidence((_DUPC,), ((4,),), None)

    with pytest.raises(ValueError, match="derived record counts"):
        BamLocusEvidence((record,), 1, {_DUPC: 2})
    with pytest.raises(ValueError, match="eligible record count"):
        BamLocusEvidence((record,), 2, {_DUPC: 1})


def test_locus_evidence_owns_a_read_only_deterministic_count_view() -> None:
    record = BamIdentityEvidence((_DUPC, _DUPA), ((4,), (5,)), None)
    supplied = {_DUPA: 1, _DUPC: 1}

    evidence = BamLocusEvidence((record,), 1, supplied)
    supplied[_DUPC] = 99

    assert isinstance(evidence.counts, MappingProxyType)
    assert tuple(evidence.counts) == (_DUPC, _DUPA)
    assert evidence.counts == {_DUPC: 1, _DUPA: 1}
    with pytest.raises(TypeError):
        evidence.counts[_DUPC] = 2  # type: ignore[index]


def test_empty_locus_has_no_winner_and_invalid_record_context_fails_closed() -> None:
    empty = collect_locus_evidence([], _locus(_DUPC), _COMPONENT)
    outside_pair = collect_locus_evidence(
        [_record(BamEditObservation(121, 0, 1, "G"))],
        _locus(_DUPC),
        _COMPONENT,
    )

    assert empty.winning_identity is None
    assert empty.eligible_record_count == 0
    assert outside_pair.record_identity_sets == ((),)


def test_unresolved_record_edit_closes_the_complete_record_binding() -> None:
    class UnresolvedComponent:
        def translate_kestrel(self, representation) -> IdentityTranslation:
            return IdentityTranslation(None, "unresolved", "motif-anchor-mismatch", False)

    evidence = collect_locus_evidence([_record(_DUPC_EDIT)], _locus(_DUPC), UnresolvedComponent())

    assert evidence.counts == {}
    assert evidence.record_identity_sets == ((),)


def test_complete_winner_binding_requires_one_common_exact_identity_from_every_supporter() -> None:
    complete_records = (
        _record(_DUPC_EDIT, _DUPA_EDIT),
        _record(_DUPC_EDIT, _DUPA_EDIT),
    )
    complete_evidence = collect_locus_evidence(complete_records, _locus(_DUPC, _DUPA), _COMPONENT)
    one_unresolved = (
        BamIdentityEvidence((_DUPC,), ((0,),), None),
        BamIdentityEvidence((), (), None),
    )

    assert (
        bind_complete_winner_identity(
            complete_records,
            BamLocusEvidence(one_unresolved, 2, {_DUPC: 1}),
            _DUPC_EDIT,
            _locus(_DUPC),
            _COMPONENT,
            2,
        )
        is None
    )
    assert (
        bind_complete_winner_identity(
            complete_records, complete_evidence, _DUPC_EDIT, _locus(_DUPC, _DUPA), _COMPONENT, 2
        )
        is None
    )

    sole_records = (_record(_DUPC_EDIT), _record(_DUPC_EDIT))
    sole_evidence = collect_locus_evidence(sole_records, _locus(_DUPC), _COMPONENT)
    assert bind_complete_winner_identity(sole_records, sole_evidence, _DUPC_EDIT, _locus(_DUPC), _COMPONENT, 2) == _DUPC
    assert bind_complete_winner_identity(sole_records, sole_evidence, _DUPC_EDIT, _locus(_DUPA), _COMPONENT, 2) is None
    assert bind_complete_winner_identity(sole_records, sole_evidence, _DUPC_EDIT, _locus(_DUPC), _COMPONENT, 3) is None


def test_complete_winner_binding_rejects_unvalidated_boundaries() -> None:
    evidence = collect_locus_evidence((), _locus(_DUPC), _COMPONENT)
    with pytest.raises(ValueError, match="record observations"):
        bind_complete_winner_identity([], evidence, _DUPC_EDIT, _locus(_DUPC), _COMPONENT, 0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="locus evidence"):
        bind_complete_winner_identity((), object(), _DUPC_EDIT, _locus(_DUPC), _COMPONENT, 0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="winning edit"):
        bind_complete_winner_identity((), evidence, object(), _locus(_DUPC), _COMPONENT, 0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="identity locus"):
        bind_complete_winner_identity((), evidence, _DUPC_EDIT, object(), _COMPONENT, 0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="supporting record count"):
        bind_complete_winner_identity((), evidence, _DUPC_EDIT, _locus(_DUPC), _COMPONENT, True)


def test_record_projection_retains_complete_window_edits_and_fails_closed() -> None:
    projected = project_record_observation(((30, 0, 1, "T"), (50, 0, 1, "")), 181, 20, 40)
    invalid_in_window = project_record_observation(((30, 0, 1, "T"), (35, 0, 1, "")), 181, 20, 40)
    unclassifiable = project_record_observation(((30, 0, 1, "T"), ("bad", 0, 1, "A")), 181, 20, 40)  # type: ignore[arg-type]

    assert projected == BamRecordObservation((BamEditObservation(30, 0, 1, "T"),), 181)
    assert invalid_in_window == BamRecordObservation((), 181)
    assert unclassifiable == BamRecordObservation((), 181)


def test_record_projection_rejects_invalid_window_boundaries() -> None:
    with pytest.raises(ValueError, match="window"):
        project_record_observation((), None, True, 40)
    with pytest.raises(ValueError, match="window"):
        project_record_observation((), None, 40, 20)


def test_collection_rejects_unvalidated_boundary_values() -> None:
    with pytest.raises(ValueError, match="BamIdentityLocus"):
        collect_locus_evidence([], object(), _COMPONENT)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="BamRecordObservation"):
        collect_locus_evidence([object()], _locus(_DUPC), _COMPONENT)  # type: ignore[list-item]


@pytest.mark.parametrize(
    "factory",
    [
        lambda: BamEditObservation(-1, 0, 1, "G"),
        lambda: BamEditObservation(1, 0, 0, ""),
        lambda: BamEditObservation(1, 0, 2, "G"),
        lambda: BamEditObservation(1, 0, 1, "g"),
        lambda: BamRecordObservation([], None),  # type: ignore[arg-type]
        lambda: BamCandidateBinding(True, _translation()),  # type: ignore[arg-type]
        lambda: BamCandidateBinding(0, object()),  # type: ignore[arg-type]
        lambda: BamIdentityLocus("", _PAIR, ()),
        lambda: BamIdentityLocus("X-X", _PAIR[:-1], ()),
        lambda: BamIdentityLocus("X-X", _PAIR, []),  # type: ignore[arg-type]
        lambda: BamIdentityLocus(
            "X-X",
            _PAIR,
            (BamCandidateBinding(0, _translation()), BamCandidateBinding(0, _translation(_DUPA))),
        ),
        lambda: BamIdentityEvidence([_DUPC], ((0,),), None),  # type: ignore[arg-type]
        lambda: BamIdentityEvidence((_DUPC,), ((),), None),
        lambda: BamIdentityEvidence((_DUPC,), ((0, 0),), None),
        lambda: BamLocusEvidence([], 0, {}),  # type: ignore[arg-type]
        lambda: BamLocusEvidence((), 0, []),  # type: ignore[arg-type]
        lambda: BamLocusEvidence((), 0, {"not-an-identity": 1}),  # type: ignore[dict-item]
        lambda: BamLocusEvidence(
            (BamIdentityEvidence((_DUPC,), ((0,),), None),),
            1,
            {_DUPC: 0},
        ),
    ],
)
def test_bam_evidence_inputs_reject_malformed_values(factory) -> None:
    with pytest.raises(ValueError):
        factory()
