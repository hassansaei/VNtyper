"""Complete-context translation into canonical MUC1 molecular identities."""

from typing import cast

import pytest

from vntyper.scripts.molecular_identity import (
    AdvntrRepresentation,
    CodingEdit,
    IdentityTranslation,
    KestrelRepresentation,
    make_coding_edit,
    make_molecular_identity,
    serialize_molecular_identity,
)
from vntyper.scripts.molecular_identity_translation import (
    bind_bam_translation,
    resolve_coding_pair_edit,
    translate_advntr_representation,
    translate_kestrel_representation,
)
from vntyper.scripts.nomenclature import MOTIFS, pair_sequence

pytestmark = pytest.mark.unit

TEST_MOTIF_MAP = {symbol: MOTIFS[symbol] for symbol in ("A", "C", "J", "S", "X")}
CANONICAL_DUPC = "MUC1-X-60-coding-v1|60|59|-|C"
CANONICAL_DUPA = "MUC1-X-60-coding-v1|60|59|-|A"


def kestrel_representation(
    observed_pair: str,
    *,
    motifs: str,
    position: int = 67,
    reference: str = "G",
    alternate: str = "GG",
) -> KestrelRepresentation:
    """Build a complete test representation with independent observed context."""
    return KestrelRepresentation(motifs, position, reference, alternate, observed_pair)


def required_pair(motifs: str) -> str:
    """Return one checked-in complete pair or fail the test fixture."""
    pair = pair_sequence(motifs)
    assert pair is not None
    return pair


def replace_neighbour(pair: str, base: str) -> str:
    """Replace a base in the half neighbouring the POS-67 affected unit."""
    assert base != pair[0]
    return base + pair[1:]


@pytest.mark.parametrize(
    ("motifs", "reference", "alternate"),
    [
        ("S-C", "G", "GG"),
        ("A-J", "C", "CG"),
    ],
)
def test_observed_dupc_shapes_resolve_to_one_identity(motifs: str, reference: str, alternate: str) -> None:
    """The two caller anchor shapes cannot split canonical dupC identity."""
    result = translate_kestrel_representation(
        kestrel_representation(required_pair(motifs), motifs=motifs, reference=reference, alternate=alternate),
        TEST_MOTIF_MAP,
    )

    assert result.status == "resolved"
    assert result.failure is None
    assert result.context_diverges is True
    assert result.identity is not None
    assert serialize_molecular_identity(result.identity) == CANONICAL_DUPC


@pytest.mark.parametrize("position", [62, 64, 67, 68])
def test_homopolymer_anchors_normalize_to_one_dupc_identity(position: int) -> None:
    """Every equivalent G anchor in the plus-strand run becomes coding 59dupC."""
    result = translate_kestrel_representation(
        kestrel_representation(required_pair("X-X"), motifs="X-X", position=position),
        TEST_MOTIF_MAP,
    )

    assert result.identity is not None
    assert serialize_molecular_identity(result.identity) == CANONICAL_DUPC


def test_dupa_remains_distinct_from_dupc() -> None:
    """A terminal A duplication must not collapse into the canonical C duplication."""
    result = translate_kestrel_representation(
        kestrel_representation(
            required_pair("X-X"),
            motifs="X-X",
            position=61,
            reference="T",
            alternate="TT",
        ),
        TEST_MOTIF_MAP,
    )

    assert result.identity is not None
    assert serialize_molecular_identity(result.identity) == CANONICAL_DUPA
    assert serialize_molecular_identity(result.identity) != CANONICAL_DUPC


def test_first_pair_half_is_oriented_without_swapping_the_identity() -> None:
    """An edit in the plus pair's first half still becomes coding terminal dupA."""
    result = translate_kestrel_representation(
        kestrel_representation(
            required_pair("X-X"),
            motifs="X-X",
            position=1,
            reference="T",
            alternate="TT",
        ),
        TEST_MOTIF_MAP,
    )

    assert result.identity is not None
    assert serialize_molecular_identity(result.identity) == CANONICAL_DUPA


def test_invalid_runtime_allele_is_unresolved() -> None:
    """A corrupted representation cannot bypass the typed allele boundary."""
    representation = kestrel_representation(required_pair("X-X"), motifs="X-X")
    object.__setattr__(representation, "alternate_allele", cast(str, "GN"))

    result = translate_kestrel_representation(representation, TEST_MOTIF_MAP)

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "invalid-allele"


def test_absent_motif_context_is_unresolved() -> None:
    """An observed pair cannot authenticate a motif key absent from the authority."""
    result = translate_kestrel_representation(
        kestrel_representation(required_pair("X-X"), motifs="missing-X"),
        TEST_MOTIF_MAP,
    )

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "missing-motif-context"


@pytest.mark.parametrize(
    ("motifs", "motif_map"),
    [
        ("X", TEST_MOTIF_MAP),
        ("X-X-X", TEST_MOTIF_MAP),
        ("X-X", {"X": "A" * 59}),
    ],
)
def test_malformed_expected_motif_context_is_absent(motifs: str, motif_map: dict[str, str]) -> None:
    """Malformed keys and non-60-base authority entries cannot define a pair."""
    result = translate_kestrel_representation(
        kestrel_representation(required_pair("X-X"), motifs=motifs),
        motif_map,
    )

    assert result.failure == "missing-motif-context"


def test_reference_allele_must_match_the_observed_pair() -> None:
    """A valid DNA REF that disagrees with the observed pair is still invalid."""
    result = translate_kestrel_representation(
        kestrel_representation(required_pair("X-X"), motifs="X-X", reference="A", alternate="AG"),
        TEST_MOTIF_MAP,
    )

    assert result.failure == "invalid-allele"


def test_changed_observed_source_unit_fails_reconstruction() -> None:
    """Only the neighbour has a distinct anchor failure; other pair drift is closed too."""
    pair = required_pair("X-X")
    changed_source = pair[:119] + ("A" if pair[-1] != "A" else "C")
    result = translate_kestrel_representation(
        kestrel_representation(changed_source, motifs="X-X"),
        TEST_MOTIF_MAP,
    )

    assert result.failure == "reconstruction-mismatch"


def test_wrong_neighbour_anchor_is_unresolved() -> None:
    """Edited-base agreement cannot hide divergence in the complete neighbour."""
    canonical_pair = required_pair("X-X")
    representation = kestrel_representation(replace_neighbour(canonical_pair, "A"), motifs="X-X")

    result = translate_kestrel_representation(representation, TEST_MOTIF_MAP)

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "motif-anchor-mismatch"


def test_pair_boundary_edit_is_unresolved() -> None:
    """One molecular identity cannot span the two repeat-unit halves."""
    pair = required_pair("X-X")
    representation = kestrel_representation(
        pair,
        motifs="X-X",
        position=60,
        reference=pair[59:61],
        alternate="A",
    )

    result = translate_kestrel_representation(representation, TEST_MOTIF_MAP)

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "pair-boundary-edit"


def test_pair_boundary_insertion_is_unresolved() -> None:
    """The gap between pair halves cannot be assigned to one X unit."""
    pair = required_pair("X-X")
    result = translate_kestrel_representation(
        kestrel_representation(
            pair,
            motifs="X-X",
            position=60,
            reference=pair[59],
            alternate=pair[59] + "A",
        ),
        TEST_MOTIF_MAP,
    )

    assert result.failure == "pair-boundary-edit"


def test_unchanged_allele_fails_reconstruction() -> None:
    """A representation that reconstructs no alternate unit has no molecular edit."""
    result = translate_kestrel_representation(
        kestrel_representation(required_pair("X-X"), motifs="X-X", alternate="G"),
        TEST_MOTIF_MAP,
    )

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "reconstruction-mismatch"


def test_advntr_uses_parsed_ru_and_position_context() -> None:
    """The complete single-base RU2 State resolves to canonical dupC."""
    result = translate_advntr_representation(AdvntrRepresentation(state="I22_2_G_LEN1", repeat_unit="2", position=22))

    assert result.identity is not None
    assert result.context_diverges is False
    assert serialize_molecular_identity(result.identity) == CANONICAL_DUPC


def test_advntr_bare_state_fails_closed() -> None:
    """State text alone is insufficient complete translation context."""
    result = translate_advntr_representation(
        AdvntrRepresentation(state="I22_2_G_LEN1", repeat_unit=None, position=None)
    )

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "missing-motif-context"


def test_advntr_parsed_position_mismatch_fails_reconstruction() -> None:
    """A stale parsed POS annotation cannot translate a different State position."""
    result = translate_advntr_representation(AdvntrRepresentation(state="I22_2_G_LEN1", repeat_unit="2", position=23))

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "reconstruction-mismatch"


def test_advntr_non_x_repeat_unit_is_unresolved() -> None:
    """A mapped but non-X repeat unit cannot inherit canonical-X identity."""
    result = translate_advntr_representation(AdvntrRepresentation(state="I23_6_G_LEN1", repeat_unit="6", position=23))

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "non-x-unit"


def test_advntr_position_outside_x_repeat_is_unresolved() -> None:
    """Modulo arithmetic cannot wrap an out-of-unit State position into canonical X."""
    result = translate_advntr_representation(AdvntrRepresentation(state="I61_2_G_LEN1", repeat_unit="2", position=61))

    assert result.identity is None
    assert result.failure == "non-x-unit"


def test_advntr_incomplete_multibase_insertion_fails_reconstruction() -> None:
    """adVNTR's first-base-only LEN4 State cannot invent the other three bases."""
    result = translate_advntr_representation(AdvntrRepresentation(state="I22_2_G_LEN4", repeat_unit="2", position=22))

    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "reconstruction-mismatch"


def test_advntr_unparseable_state_is_an_invalid_allele() -> None:
    """Parsed context cannot turn unrecognized State syntax into an identity."""
    result = translate_advntr_representation(AdvntrRepresentation(state="banana", repeat_unit="2", position=1))

    assert result.failure == "invalid-allele"


def test_advntr_repeat_unit_annotation_must_match_state() -> None:
    """A stale parsed RU annotation cannot authenticate a different State unit."""
    result = translate_advntr_representation(AdvntrRepresentation(state="I22_2_G_LEN1", repeat_unit="5", position=22))

    assert result.failure == "reconstruction-mismatch"


def test_advntr_complete_consecutive_deletion_reconstructs_the_unit() -> None:
    """Five complete deletion State parts become the exact canonical 1_5 deletion."""
    state = "D17_2&D18_2&D19_2&D20_2&D21_2"
    result = translate_advntr_representation(AdvntrRepresentation(state=state, repeat_unit="2", position=17))

    assert result.identity is not None
    assert serialize_molecular_identity(result.identity) == "MUC1-X-60-coding-v1|1|5|GCCCA|-"


def test_advntr_noncontiguous_complete_state_preserves_each_edit() -> None:
    """Separated complete State changes remain co-occurring edits, not one delins."""
    state = "I28_2_C_LEN1&I24_2_C_LEN1"
    result = translate_advntr_representation(AdvntrRepresentation(state=state, repeat_unit="2", position=28))

    assert result.identity is not None
    assert result.identity.edits == (
        CodingEdit(54, 53, "", "G"),
        CodingEdit(58, 57, "", "G"),
    )


def test_advntr_complete_colocated_delins_is_one_edit() -> None:
    """A complete deletion and insertion at one State position remain one replacement."""
    state = "D27_2&I27_2_A_LEN1"
    result = translate_advntr_representation(AdvntrRepresentation(state=state, repeat_unit="2", position=27))

    assert result.identity is not None
    assert result.identity.edits == (CodingEdit(55, 55, "C", "T"),)


def test_complete_unit_resolver_rejects_an_unchanged_or_foreign_reference() -> None:
    """The internal complete-unit seam cannot create identity from absent or foreign change."""
    canonical = required_pair("X-X")
    coding_x = canonical[:60][::-1].translate(str.maketrans("ACGT", "TGCA"))

    assert resolve_coding_pair_edit(coding_x, coding_x) is None
    assert resolve_coding_pair_edit("A" * 60, "C" + "A" * 59) is None


def test_bam_binding_preserves_complete_cooccurring_edit_set() -> None:
    """Binding evidence cannot keep only the first edit of a complete identity."""
    identity = make_molecular_identity(
        (
            make_coding_edit(10, 10, "G", "A"),
            make_coding_edit(60, 59, "", "C"),
        )
    )
    translation = IdentityTranslation(identity=identity, status="resolved", failure=None, context_diverges=True)

    bound = bind_bam_translation(translation, identity)

    assert bound is translation
    assert bound.identity is not None
    assert bound.identity.edits == (
        CodingEdit(10, 10, "G", "A"),
        CodingEdit(60, 59, "", "C"),
    )


def test_bam_binding_rejects_a_partial_edit_set() -> None:
    """A BAM subset is not complete evidence for a multi-edit Kestrel identity."""
    identity = make_molecular_identity(
        (
            make_coding_edit(10, 10, "G", "A"),
            make_coding_edit(60, 59, "", "C"),
        )
    )
    partial = make_molecular_identity((make_coding_edit(10, 10, "G", "A"),))
    translation = IdentityTranslation(identity=identity, status="resolved", failure=None, context_diverges=False)

    assert bind_bam_translation(translation, partial) is None


def test_bam_binding_does_not_resolve_an_unresolved_kestrel_call() -> None:
    """BAM evidence cannot create a caller representation or molecular identity."""
    unresolved = IdentityTranslation(
        identity=None,
        status="unresolved",
        failure="missing-motif-context",
        context_diverges=False,
    )
    identity = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))

    assert bind_bam_translation(unresolved, identity) is None
