"""Unit tests for canonical molecular identity values."""

from collections.abc import Callable

import pytest

from vntyper.scripts.molecular_identity import (
    AdvntrRepresentation,
    CodingEdit,
    EvidenceDisposition,
    FrameConsequence,
    IdentityDecision,
    IdentityObservation,
    IdentityTranslation,
    KestrelRepresentation,
    MolecularIdentity,
    make_coding_edit,
    make_molecular_identity,
    parse_molecular_identity,
    serialize_molecular_identity,
)

pytestmark = pytest.mark.unit


def test_canonical_dupc_has_stable_identity() -> None:
    """Canonical dupC retains its versioned serialization."""
    identity = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
    assert isinstance(identity, MolecularIdentity)
    assert serialize_molecular_identity(identity) == "MUC1-X-60-coding-v1|60|59|-|C"


def test_empty_identity_is_rejected() -> None:
    """An identity must include an edit."""
    with pytest.raises(ValueError, match="Molecular identity requires at least one edit"):
        make_molecular_identity(())


def test_direct_invalid_insertion_is_rejected() -> None:
    """Insertions use the canonical between-base coordinate convention."""
    with pytest.raises(ValueError, match="Insertion coordinates require end == start - 1"):
        CodingEdit(60, 60, "", "C")


@pytest.mark.parametrize(
    "factory",
    [
        lambda: CodingEdit(60, 59, "", "c"),
        lambda: CodingEdit(61, 60, "", "C"),
        lambda: CodingEdit(1, False, "", "C"),
    ],
)
def test_case_and_range_are_strict(factory: Callable[[], CodingEdit]) -> None:
    """Alleles and positions reject values outside the canonical coding model."""
    with pytest.raises(ValueError):
        factory()


def test_foreign_identity_version_is_rejected() -> None:
    """Serialized identities accept only the supported reference version."""
    with pytest.raises(ValueError, match="Unsupported molecular identity reference"):
        parse_molecular_identity("MUC1-X-60-coding-v2|60|59|-|C")


def test_translation_state_is_consistent() -> None:
    """Resolved translations require an identity and no failure reason."""
    with pytest.raises(ValueError):
        IdentityTranslation(identity=None, status="resolved", failure=None, context_diverges=False)


def test_parser_round_trips_multiple_sorted_edits() -> None:
    """The strict parser preserves each canonical edit in a stable identity."""
    serialized = "MUC1-X-60-coding-v1|10|10|A|G;60|59|-|C"
    assert serialize_molecular_identity(parse_molecular_identity(serialized)) == serialized


@pytest.mark.parametrize(
    "serialized",
    [
        "MUC1-X-60-coding-v1",
        "MUC1-X-60-coding-v1|60|59|-|-",
        "MUC1-X-60-coding-v1|60|59||C",
    ],
)
def test_parser_rejects_malformed_or_empty_edits(serialized: str) -> None:
    """Malformed wire values cannot create ambiguous identities."""
    with pytest.raises(ValueError):
        parse_molecular_identity(serialized)


def test_identity_rejects_unsorted_and_overlapping_edits() -> None:
    """Identity equality is only defined for sorted, non-overlapping edits."""
    first = make_coding_edit(10, 10, "A", "G")
    second = make_coding_edit(9, 9, "C", "T")
    overlap = make_coding_edit(10, 11, "AC", "")

    with pytest.raises(ValueError, match="sorted"):
        make_molecular_identity((first, second))
    with pytest.raises(ValueError, match="overlap"):
        make_molecular_identity((first, overlap))


def test_translation_rejects_unresolved_state_with_identity() -> None:
    """Unresolved translations cannot carry a comparable identity."""
    identity = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
    with pytest.raises(ValueError):
        IdentityTranslation(
            identity=identity,
            status="unresolved",
            failure="missing-motif-context",
            context_diverges=False,
        )


@pytest.mark.parametrize(
    "factory",
    [
        lambda: CodingEdit(10, 10, "", ""),
        lambda: CodingEdit(10, 9, "A", ""),
        lambda: CodingEdit(10, 11, "A", ""),
        lambda: MolecularIdentity(reference="MUC1-X-60-coding-v2", edits=(make_coding_edit(60, 59, "", "C"),)),  # type: ignore[arg-type]
        lambda: MolecularIdentity(reference="MUC1-X-60-coding-v1", edits=[]),  # type: ignore[arg-type]
    ],
)
def test_edit_and_identity_value_types_reject_noncanonical_states(factory: Callable[[], object]) -> None:
    """Direct value construction cannot bypass canonical edit normalization."""
    with pytest.raises(ValueError):
        factory()


@pytest.mark.parametrize(
    "translation",
    [
        lambda: IdentityTranslation(identity=None, status="unknown", failure=None, context_diverges=False),  # type: ignore[arg-type]
        lambda: IdentityTranslation(
            identity=None,
            status="unresolved",
            failure="unknown-failure",  # type: ignore[arg-type]
            context_diverges=False,
        ),
        lambda: IdentityTranslation(
            identity=None,
            status="unresolved",
            failure="missing-motif-context",
            context_diverges=True,
        ),
    ],
)
def test_translation_value_type_rejects_unknown_or_divergent_unresolved_states(
    translation: Callable[[], IdentityTranslation],
) -> None:
    """Closed translation states cannot invent a resolution or failure vocabulary."""
    with pytest.raises(ValueError):
        translation()


@pytest.mark.parametrize(
    "factory",
    [
        lambda: KestrelRepresentation(motifs="", position=1, reference_allele="A", alternate_allele="G"),
        lambda: KestrelRepresentation(motifs="M1-M2", position=0, reference_allele="A", alternate_allele="G"),
        lambda: KestrelRepresentation(motifs="M1-M2", position=1, reference_allele="a", alternate_allele="G"),
        lambda: AdvntrRepresentation(state="", repeat_unit=None, position=None),
        lambda: AdvntrRepresentation(state="3-1", repeat_unit="", position=None),
        lambda: AdvntrRepresentation(state="3-1", repeat_unit=None, position=0),
        lambda: EvidenceDisposition(value="unknown"),  # type: ignore[arg-type]
    ],
)
def test_representation_and_disposition_values_reject_invalid_fields(factory: Callable[[], object]) -> None:
    """Caller-local representations and dispositions retain closed value vocabularies."""
    with pytest.raises(ValueError):
        factory()


def test_value_types_reject_inconsistent_observation_and_decision_states() -> None:
    """The supporting identity values preserve source and decision boundaries."""
    identity = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
    kestrel = KestrelRepresentation(motifs="M1-M2", position=1, reference_allele="A", alternate_allele="G")
    advntr = AdvntrRepresentation(state="3-1", repeat_unit="X", position=1)
    translation = IdentityTranslation(identity=identity, status="resolved", failure=None, context_diverges=False)
    consequence = FrameConsequence(net_length_change=1, is_frameshift=True)

    observation = IdentityObservation(
        source="kestrel",
        representation=kestrel,
        translation=translation,
        frame_consequence=consequence,
    )
    disposition = EvidenceDisposition(value="admissible")
    decision = IdentityDecision(
        identity=identity,
        tier="A",
        molecular_agreement=True,
        abstention_reason=None,
    )

    assert observation.source == "kestrel"
    assert disposition.value == "admissible"
    assert decision.identity == identity

    with pytest.raises(ValueError, match="does not match representation"):
        IdentityObservation(
            source="advntr",
            representation=kestrel,
            translation=translation,
            frame_consequence=consequence,
        )
    with pytest.raises(ValueError, match="inconsistent"):
        FrameConsequence(net_length_change=1, is_frameshift=False)
    with pytest.raises(ValueError, match="requires an abstention reason"):
        IdentityDecision(identity=None, tier="abstained", molecular_agreement=False, abstention_reason=None)
    with pytest.raises(ValueError, match="does not accept an identity"):
        IdentityDecision(identity=None, tier="A", molecular_agreement=True, abstention_reason="unresolved")
    with pytest.raises(ValueError, match="cannot include an abstention reason"):
        IdentityDecision(identity=identity, tier="A", molecular_agreement=True, abstention_reason="unresolved")
    assert advntr.state == "3-1"
