"""Identity candidate capture before legacy caller projection."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import replace
from typing import Any, cast

import numpy as np
import pandas as pd
import pytest

from vntyper.scripts.identity_candidate_persistence import (
    IDENTITY_CAPTURE_COLUMNS,
    IDENTITY_SELECTION_COLUMNS,
    candidate_capture_cells,
    complete_candidate_projection_cells,
    parse_selected_candidate_cells,
    selected_candidate_cells,
)
from vntyper.scripts.identity_candidates import (
    IdentityCandidateSet,
    IdentityTranslationComponent,
    RawRepresentationKey,
    capture_advntr_observations,
    capture_kestrel_observations,
    overlay_legacy_projection,
    translation_component_from_config,
    with_candidate_evidence,
)
from vntyper.scripts.molecular_identity import (
    AdvntrRepresentation,
    IdentityTranslation,
    make_coding_edit,
    make_molecular_identity,
)

pytestmark = pytest.mark.unit

_X = "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC"
_FIVE = "TGGGGGGGCGGTGGAGCCCGGGGCTGGCTTGTTGTCCGGGGCTGAGGTGACATCGTGGGC"
_S = "TGCGGGGGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGGGCTGAGGTGACACCGTGGGC"
_C = "TTGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGC"
_A = "TGCGGGCGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGGGCCGAGGTGACACCGTGGGC"
_J = "TGCGGGCGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGTGCCGAGGTGACACCGTGGGC"

TEST_CONFIG = {
    "motifs": {"X": _X, "5": _FIVE, "S": _S, "C": _C, "A": _A, "J": _J},
    "advntr": {"mappable_repeat_units": {"2": "X", "5": "V", "6": "C"}, "rotation_offset": 39},
}
REAL_COMPONENT = translation_component_from_config(TEST_CONFIG)
ALL_GATES = (
    "is_frameshift",
    "is_valid_frameshift",
    "depth_confidence_pass",
    "alt_filter_pass",
    "motif_filter_pass",
    "flag_filter_pass",
)

CAPTURE_ORDINAL_COLUMN = "__Identity_Observation_Ordinal"
SELECTED_ORDINAL_COLUMN = "__Identity_Selected_Observation_Ordinal"
GROUP_CONTEXT_COLUMN = "__Identity_Group_Context_Diverges"


def _translation(inserted: str, *, context_diverges: bool = False) -> IdentityTranslation:
    identity = make_molecular_identity((make_coding_edit(60, 59, "", inserted),))
    return IdentityTranslation(identity, "resolved", None, context_diverges)


class _MotifSensitiveComponent:
    """Test translator making motif context visibly part of independent capture."""

    def __init__(self, *, equivalent: bool = False, divergent_motif: str | None = None) -> None:
        self.equivalent = equivalent
        self.divergent_motif = divergent_motif

    def translate_kestrel(self, representation: Any) -> IdentityTranslation:
        inserted = "C" if self.equivalent or representation.motifs == "X-5" else "A"
        return _translation(inserted, context_diverges=representation.motifs == self.divergent_motif)

    def translate_advntr(self, representation: Any) -> IdentityTranslation:
        inserted = "C" if representation.state.endswith("G_LEN1") else "A"
        return _translation(inserted)


class _UnresolvedComponent(_MotifSensitiveComponent):
    """Translator closing every representation without an identity."""

    def translate_kestrel(self, representation: Any) -> IdentityTranslation:
        return IdentityTranslation(None, "unresolved", "motif-anchor-mismatch", False)


def _kestrel_row(
    motifs: str,
    pair_sequence: str,
    *,
    ref: str = "G",
    alt: str = "GG",
    support: int = 4,
    depth: int = 400,
    flag: str = "Not flagged",
    **gates: bool,
) -> dict[str, object]:
    row: dict[str, object] = {
        "Motifs": motifs,
        "POS": 67,
        "REF": ref,
        "ALT": alt,
        "Motif_sequence": pair_sequence,
        "Estimated_Depth_AlternateVariant": support,
        "Estimated_Depth_Variant_ActiveRegion": depth,
        "Flag": flag,
    }
    row.update(dict.fromkeys(ALL_GATES, True))
    row.update(gates)
    return row


def _captured(
    rows: Iterable[Mapping[str, object]],
    component: Any = REAL_COMPONENT,
) -> IdentityCandidateSet:
    return capture_kestrel_observations(rows, component)


def _overlay_all(captured: IdentityCandidateSet, selected_index: int = 0) -> IdentityCandidateSet:
    ordinals = tuple(candidate.observation_ordinal for candidate in captured.candidates)
    evidenced = with_candidate_evidence(captured, [_row_for_candidate(candidate) for candidate in captured.candidates])
    return overlay_legacy_projection(evidenced, ordinals, ordinals[selected_index])


def _row_for_candidate(candidate: Any, **changes: object) -> dict[str, object]:
    values = candidate.row_key.values
    row = _kestrel_row(
        str(values[0]),
        str(candidate.observation.representation.pair_sequence),
        ref=str(values[2]),
        alt=str(values[3]),
        support=int(candidate.support),
        depth=int(candidate.depth),
    )
    row["POS"] = int(values[1])
    row[CAPTURE_ORDINAL_COLUMN] = str(candidate.observation_ordinal)
    row.update(changes)
    return row


def _persisted_cells() -> dict[str, str]:
    captured = _captured([_kestrel_row("S-C", _C + _S)])
    selected = _overlay_all(captured)
    return {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}


def test_same_raw_tuple_under_two_motifs_retains_two_hypotheses() -> None:
    """Capture before three-field projection retains both resolved identities."""
    rows = [
        _kestrel_row("X-5", _FIVE + _X),
        _kestrel_row("S-C", _C + _S),
    ]
    captured = _captured(rows, _MotifSensitiveComponent())
    passing_ordinals = tuple(candidate.observation_ordinal for candidate in captured.candidates)
    legacy_selected_ordinal = passing_ordinals[0]

    candidates = overlay_legacy_projection(captured, passing_ordinals, legacy_selected_ordinal)

    assert candidates.identity_hypothesis_count == 2
    assert candidates.selected_row_key == captured.candidates[0].row_key
    assert {candidate.row_key.values for candidate in candidates.candidates} == {
        ("X-5", 67, "G", "GG"),
        ("S-C", 67, "G", "GG"),
    }


def test_complete_candidate_projection_persists_every_passing_hypothesis_before_legacy_selection() -> None:
    """Mutation caught: only the later legacy-selected candidate receives authoritative metadata."""
    captured = _captured(
        [_kestrel_row("X-5", _FIVE + _X), _kestrel_row("S-C", _C + _S)],
        _MotifSensitiveComponent(),
    )
    rows = [_row_for_candidate(candidate) for candidate in captured.candidates]
    evidenced = with_candidate_evidence(captured, rows)

    projections = complete_candidate_projection_cells(evidenced, (0, 1))

    assert tuple(sorted(projections)) == (0, 1)
    assert projections[0][SELECTED_ORDINAL_COLUMN] == "0"
    assert projections[1][SELECTED_ORDINAL_COLUMN] == "1"
    assert projections[0]["__Identity_Hypothesis_Count"] == "2"
    assert projections[1]["__Identity_Hypothesis_Count"] == "2"


def test_equivalent_representations_do_not_sum_support() -> None:
    """The exact legacy row supplies support for an equivalent identity group."""
    rows = [
        _kestrel_row("S-C", _C + _S, support=4),
        _kestrel_row("A-J", _J + _A, ref="C", alt="CG", support=7),
    ]

    candidates = _overlay_all(_captured(rows), selected_index=0)

    assert candidates.selected_support == 4
    assert candidates.equivalent_representation_count == 2
    assert len(candidates.selected_group.candidates) == 2


def test_group_union_cannot_unblock_tier_a() -> None:
    """One divergent equivalent representation blocks the complete group."""
    rows = [_kestrel_row("X-5", _FIVE + _X), _kestrel_row("S-C", _C + _S)]
    candidates = _overlay_all(
        _captured(rows, _MotifSensitiveComponent(equivalent=True, divergent_motif="S-C")),
        selected_index=0,
    )

    assert candidates.selected_group.context_diverges is True


def test_group_unions_every_false_gate_and_flag() -> None:
    """Equivalent observations retain every blocker instead of borrowing a clean row."""
    captured = _captured(
        [_kestrel_row("X-5", _FIVE + _X), _kestrel_row("S-C", _C + _S)],
        _MotifSensitiveComponent(equivalent=True),
    )
    evidence = [
        _row_for_candidate(captured.candidates[0]),
        _row_for_candidate(
            captured.candidates[1],
            motif_filter_pass=False,
            flag_filter_pass=False,
            Flag="Context_Doubt, Technical_Artifact",
        ),
    ]
    evidenced = with_candidate_evidence(captured, evidence)
    selected = overlay_legacy_projection(
        evidenced,
        (captured.candidates[0].observation_ordinal,),
        captured.candidates[0].observation_ordinal,
    )

    assert selected.selected_group.blocking_gates == frozenset({"motif_filter_pass", "flag_filter_pass"})
    assert selected.selected_group.flags == frozenset({"Context_Doubt", "Technical_Artifact"})


def test_advntr_positive_states_are_translated_before_reconciliation() -> None:
    """Two complete positive State rows expose two caller-local hypotheses."""
    rows = [
        {
            "Variant": "I22_2_G_LEN1",
            "RU": "2",
            "POS": "22",
            "Net_indel_length": 1,
            "NumberOfSupportingReads": 6,
            "MeanCoverage": 25,
            "Flag": "Not flagged",
        },
        {
            "Variant": "I22_2_A_LEN1",
            "RU": "2",
            "POS": "22",
            "Net_indel_length": 1,
            "NumberOfSupportingReads": 9,
            "MeanCoverage": 30,
            "Flag": "Not flagged",
        },
    ]

    candidates = capture_advntr_observations(rows, REAL_COMPONENT)
    representations = tuple(candidate.observation.representation for candidate in candidates.candidates)
    assert all(isinstance(representation, AdvntrRepresentation) for representation in representations)
    advntr_representations = cast(tuple[AdvntrRepresentation, ...], representations)

    assert candidates.identity_hypothesis_count == 2
    assert tuple(representation.repeat_units for representation in advntr_representations) == (
        ("2",),
        ("2",),
    )
    assert tuple(representation.positions for representation in advntr_representations) == (
        (22,),
        (22,),
    )


def test_below_floor_resolved_row_does_not_enter_hypothesis_count() -> None:
    """Only keys surviving all six unchanged gates count as considered identities."""
    captured = _captured(
        [_kestrel_row("X-5", _FIVE + _X), _kestrel_row("S-C", _C + _S)],
        _MotifSensitiveComponent(),
    )
    evidenced = with_candidate_evidence(
        captured,
        [
            _row_for_candidate(captured.candidates[0]),
            _row_for_candidate(captured.candidates[1], depth_confidence_pass=False),
        ],
    )
    selected = overlay_legacy_projection(
        evidenced,
        (captured.candidates[0].observation_ordinal,),
        captured.candidates[0].observation_ordinal,
    )

    assert selected.identity_hypothesis_count == 1
    assert selected.candidates[1].eligible is False
    assert selected.candidates[1].blocking_gates == frozenset({"depth_confidence_pass"})


def test_raw_key_excludes_support_flags_and_display_name() -> None:
    """Only Motifs/POS/REF/ALT define a retained Kestrel representation."""
    left = _captured([_kestrel_row("X-5", _FIVE + _X, support=4, flag="first")]).candidates[0]
    right_row = _kestrel_row("X-5", _FIVE + _X, support=99, flag="second")
    right_row["Nomenclature"] = "different display"
    right = _captured([right_row]).candidates[0]

    assert left.row_key == right.row_key
    assert left.row_key.values == ("X-5", 67, "G", "GG")


@pytest.mark.parametrize("position", [67, np.int64(67), "67"])
def test_kestrel_capture_normalizes_dataframe_scalar_forms(position: object) -> None:
    """DataFrame-backed Kestrel rows retain strict built-in numeric value types."""
    row = _kestrel_row("X-5", _FIVE + _X)
    row["POS"] = position
    row["Estimated_Depth_AlternateVariant"] = np.int64(4)
    row["Estimated_Depth_Variant_ActiveRegion"] = np.float64(400.0)

    candidate = _captured([row]).candidates[0]

    assert candidate.row_key.values == ("X-5", 67, "G", "GG")
    assert type(candidate.row_key.values[1]) is int
    assert candidate.support == 4
    assert type(candidate.support) is int
    assert candidate.depth == 400.0
    assert type(candidate.depth) is float


@pytest.mark.parametrize("position", ["067", "+67", "67.0", " 67", "67 ", "", "A"])
def test_kestrel_capture_rejects_noncanonical_position_strings(position: str) -> None:
    """Only the unsigned decimal spelling emitted by the VCF dataframe is accepted."""
    row = _kestrel_row("X-5", _FIVE + _X)
    row["POS"] = position

    with pytest.raises(ValueError, match="POS must be an integer"):
        _captured([row])


def test_selected_row_alone_supplies_depth_as_well_as_support() -> None:
    """Selecting the second equivalent row does not retain the first row's depth."""
    rows = [
        _kestrel_row("S-C", _C + _S, support=4, depth=400),
        _kestrel_row("A-J", _J + _A, ref="C", alt="CG", support=7, depth=700),
    ]

    candidates = _overlay_all(_captured(rows), selected_index=1)

    assert candidates.selected_support == 7
    assert candidates.selected_depth == 700


def test_selected_candidate_metadata_survives_tsv_round_trip(tmp_path: Any) -> None:
    """A4 can replay selected identity evidence from TSV scalar columns alone."""
    rows = [
        _kestrel_row("S-C", _C + _S, support=4),
        _kestrel_row("A-J", _J + _A, ref="C", alt="CG", support=7, motif_filter_pass=False),
    ]
    captured = _captured(rows)
    for candidate, row in zip(captured.candidates, rows, strict=True):
        row[CAPTURE_ORDINAL_COLUMN] = str(candidate.observation_ordinal)
    evidenced = with_candidate_evidence(captured, rows)
    selected = overlay_legacy_projection(
        evidenced,
        (captured.candidates[0].observation_ordinal,),
        captured.candidates[0].observation_ordinal,
    )
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}
    path = tmp_path / "candidate.tsv"
    pd.DataFrame([cells]).to_csv(path, sep="\t", index=False)

    row = pd.read_csv(path, sep="\t", dtype=str).iloc[0].to_dict()
    replayed = parse_selected_candidate_cells(row)

    assert replayed.translation == selected.selected_candidate.observation.translation
    assert replayed.selected_row_key == selected.selected_row_key
    assert replayed.equivalent_representation_count == 2
    assert replayed.identity_hypothesis_count == 1
    assert replayed.blocking_gates == frozenset({"motif_filter_pass"})
    assert replayed.flags == frozenset()
    assert set(cells) == set(IDENTITY_CAPTURE_COLUMNS) | set(IDENTITY_SELECTION_COLUMNS)


def test_unresolved_selected_row_persists_zero_equivalent_representations() -> None:
    """An unresolved positive row remains replayable without inventing a group."""
    captured = _captured([_kestrel_row("X-5", _FIVE + _X)], _UnresolvedComponent())
    evidenced = with_candidate_evidence(captured, [_row_for_candidate(captured.candidates[0])])
    selected = overlay_legacy_projection(
        evidenced,
        (captured.candidates[0].observation_ordinal,),
        captured.candidates[0].observation_ordinal,
    )

    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}
    replayed = parse_selected_candidate_cells(cells)

    assert replayed.translation.status == "unresolved"
    assert replayed.translation.identity is None
    assert replayed.equivalent_representation_count == 0
    assert replayed.identity_hypothesis_count == 0


def test_component_is_built_only_from_explicit_config() -> None:
    """Checked-in motif, RU-map, and rotation values are copied into one component."""
    config: dict[str, Any] = {
        "motifs": {"X": _X},
        "advntr": {"mappable_repeat_units": {"2": "X"}, "rotation_offset": 39},
    }

    component = translation_component_from_config(config)
    config["motifs"]["X"] = "A" * 60
    config["advntr"]["rotation_offset"] = 1

    assert isinstance(component, IdentityTranslationComponent)
    assert component.kestrel_motifs["X"] == _X
    assert component.advntr_rotation_offset == 39


@pytest.mark.parametrize(
    "config",
    [
        {},
        {"motifs": {"X": _X}},
        {"motifs": {"X": _X}, "advntr": {"mappable_repeat_units": {"2": "X"}}},
    ],
)
def test_incomplete_component_config_fails_closed(config: Mapping[str, object]) -> None:
    """Missing injected authorities cannot fall back to module constants."""
    with pytest.raises((KeyError, ValueError, TypeError)):
        translation_component_from_config(config)


def test_overlay_rejects_a_selected_key_that_did_not_pass() -> None:
    """Metadata cannot claim selection of a row outside the six-gate projection."""
    captured = _captured([_kestrel_row("X-5", _FIVE + _X), _kestrel_row("S-C", _C + _S)])

    with pytest.raises(ValueError, match="selected Kestrel row did not pass"):
        overlay_legacy_projection(
            captured,
            (captured.candidates[0].observation_ordinal,),
            captured.candidates[1].observation_ordinal,
        )


def test_evidence_overlay_requires_all_six_gate_columns() -> None:
    """Candidate eligibility cannot be computed while any legacy gate is absent."""
    captured = _captured([_kestrel_row("X-5", _FIVE + _X)])
    incomplete = _row_for_candidate(captured.candidates[0])
    incomplete.pop("flag_filter_pass")

    with pytest.raises(ValueError, match="flag_filter_pass"):
        with_candidate_evidence(captured, [incomplete])


def test_raw_representation_key_rejects_wrong_source_shape() -> None:
    """A malformed serialized key cannot masquerade as a candidate lookup key."""
    with pytest.raises(ValueError, match="Kestrel raw representation"):
        RawRepresentationKey("kestrel", ("X-5", "67", "G", "GG"))  # type: ignore[arg-type]


@pytest.mark.parametrize("rotation", [True, 0, 61])
def test_component_rejects_invalid_rotation(rotation: object) -> None:
    """The injected rotation authority is a closed one-through-sixty integer."""
    with pytest.raises(ValueError, match="rotation offset"):
        IdentityTranslationComponent({"X": _X}, {"2": "X"}, rotation)  # type: ignore[arg-type]


@pytest.mark.parametrize("motifs", [{}, {"": _X}, {"X": ""}])
def test_component_rejects_invalid_motif_maps(motifs: dict[str, str]) -> None:
    """Empty keys, values, and maps cannot become translation authorities."""
    with pytest.raises(ValueError, match="motif map"):
        IdentityTranslationComponent(motifs, {"2": "X"}, 39)


@pytest.mark.parametrize(
    "config",
    [
        {"motifs": [], "advntr": {}},
        {"motifs": {"X": _X}, "advntr": {"mappable_repeat_units": [], "rotation_offset": 39}},
    ],
)
def test_component_config_rejects_non_mapping_sections(config: Mapping[str, object]) -> None:
    """The injection boundary does not coerce sequence-like configuration."""
    with pytest.raises(TypeError, match="mapping"):
        translation_component_from_config(config)


def test_candidate_value_objects_reject_inconsistent_evidence() -> None:
    """Candidate records fail closed on source, evidence, gate, flag, and eligibility drift."""
    candidate = _captured([_kestrel_row("X-5", _FIVE + _X)]).candidates[0]
    advntr_key = RawRepresentationKey("advntr", ("I22_2_G_LEN1", ("2",), (22,)))
    mutations = (
        {"row_key": advntr_key},
        {"support": -1},
        {"blocking_gates": frozenset({"not_a_gate"})},
        {"flags": frozenset({""})},
        {"eligible": 1},
    )

    for changes in mutations:
        with pytest.raises(ValueError):
            replace(candidate, **changes)


def test_candidate_sets_reject_inconsistent_selection_state() -> None:
    """Candidate sets reject mixed callers, absent selections, and ambiguous selected keys."""
    candidate = _captured([_kestrel_row("X-5", _FIVE + _X)]).candidates[0]
    with pytest.raises(ValueError, match="Unsupported candidate-set source"):
        IdentityCandidateSet("other", ())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="must be a tuple"):
        IdentityCandidateSet("kestrel", [])  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="cannot mix"):
        IdentityCandidateSet("advntr", (candidate,))
    with pytest.raises(ValueError, match="absent"):
        IdentityCandidateSet("kestrel", (), candidate.observation_ordinal)
    with pytest.raises(ValueError, match="no selected row"):
        _ = IdentityCandidateSet("kestrel", (candidate,)).selected_candidate


def test_evidence_and_projection_reject_incomplete_or_unknown_raw_keys() -> None:
    """The six-gate overlay cannot silently lose, coerce, or invent a raw key."""
    captured = _captured([_kestrel_row("X-5", _FIVE + _X)])
    row = _row_for_candidate(captured.candidates[0])
    row["is_frameshift"] = 1
    with pytest.raises(ValueError, match="must be a boolean"):
        with_candidate_evidence(captured, [row])
    with pytest.raises(ValueError, match="missing a captured"):
        with_candidate_evidence(captured, [])
    with pytest.raises(ValueError, match="only for Kestrel"):
        with_candidate_evidence(IdentityCandidateSet("advntr", ()), [])
    with pytest.raises(ValueError, match="unknown observation ordinal"):
        overlay_legacy_projection(captured, (99,), 99)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        (IDENTITY_SELECTION_COLUMNS[0], "not-json", "valid JSON"),
        (IDENTITY_SELECTION_COLUMNS[0], "{}", "invalid shape"),
        (IDENTITY_SELECTION_COLUMNS[3], '["z","a"]', "sorted unique"),
        (IDENTITY_SELECTION_COLUMNS[4], "[1]", "non-empty strings"),
        (IDENTITY_CAPTURE_COLUMNS[4], 1, "canonical boolean"),
        (IDENTITY_SELECTION_COLUMNS[1], -1, "non-negative integer"),
        (IDENTITY_SELECTION_COLUMNS[2], "01", "non-negative integer"),
        (IDENTITY_SELECTION_COLUMNS[1], "0", "positive equivalent"),
    ],
)
def test_persisted_candidate_codec_rejects_malformed_scalars(column: str, value: object, message: str) -> None:
    """Malformed TSV scalar cells fail closed instead of being reinterpreted."""
    cells: dict[str, object] = dict(_persisted_cells())
    cells[column] = value

    with pytest.raises(ValueError, match=message):
        parse_selected_candidate_cells(cells)


def test_persisted_candidate_codec_parses_complete_advntr_raw_key() -> None:
    """The persistence codec reconstructs nested adVNTR RU/POS tuples exactly."""
    cells = _persisted_cells()
    raw_key = '{"source":"advntr","values":["I22_2_G_LEN1",["2"],[22]]}'
    cells[IDENTITY_CAPTURE_COLUMNS[0]] = raw_key
    cells[IDENTITY_SELECTION_COLUMNS[0]] = raw_key

    replayed = parse_selected_candidate_cells(cells)

    assert replayed.selected_row_key == RawRepresentationKey("advntr", ("I22_2_G_LEN1", ("2",), (22,)))


def test_duplicate_raw_keys_retain_distinct_ordinal_bound_observations() -> None:
    """Exact duplicate raw keys remain separate and selection binds the chosen observation."""
    captured = _captured(
        [
            _kestrel_row("S-C", _C + _S, support=4, depth=400),
            _kestrel_row("S-C", _C + _S, support=7, depth=700),
        ]
    )
    assert [candidate.observation_ordinal for candidate in captured.candidates] == [0, 1]
    assert captured.candidates[0].row_key == captured.candidates[1].row_key
    evidence = []
    for candidate in captured.candidates:
        row = _row_for_candidate(candidate)
        row[CAPTURE_ORDINAL_COLUMN] = str(candidate.observation_ordinal)
        evidence.append(row)

    evidenced = with_candidate_evidence(captured, evidence)
    selected = overlay_legacy_projection(evidenced, (1,), 1)

    assert len(selected.candidates) == 2
    assert selected.selected_row_key == captured.candidates[1].row_key
    assert selected.selected_candidate.observation_ordinal == 1
    assert selected.selected_support == 7
    assert selected.selected_depth == 700
    assert selected.equivalent_representation_count == 2

    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}
    replayed = parse_selected_candidate_cells(cells)
    assert replayed.selected_observation_ordinal == 1


def test_candidate_set_requires_unique_nonnegative_observation_ordinals() -> None:
    """An ordinal is stable row identity and therefore cannot be negative or duplicated."""
    captured = _captured([_kestrel_row("S-C", _C + _S), _kestrel_row("A-J", _J + _A, ref="C", alt="CG")])
    first, second = captured.candidates

    with pytest.raises(ValueError, match="non-negative"):
        replace(first, observation_ordinal=-1)
    with pytest.raises(ValueError, match="unique"):
        IdentityCandidateSet("kestrel", (first, replace(second, observation_ordinal=first.observation_ordinal)))
    with pytest.raises(ValueError, match="non-negative"):
        IdentityCandidateSet("kestrel", (first, second), False)
    with pytest.raises(ValueError, match="non-negative"):
        overlay_legacy_projection(captured, (False,), False)


def test_evidence_ordinal_must_bind_its_original_raw_key() -> None:
    """A projected row cannot borrow another captured observation's ordinal."""
    captured = _captured([_kestrel_row("S-C", _C + _S)])
    row = _row_for_candidate(captured.candidates[0], Motifs="A-J", REF="C", ALT="CG")
    row[CAPTURE_ORDINAL_COLUMN] = "0"

    with pytest.raises(ValueError, match="raw key does not match"):
        with_candidate_evidence(captured, [row])


def test_production_comma_joined_flags_are_trimmed_and_unioned() -> None:
    """Candidate flags consume the exact comma-and-space format emitted by add_flags."""
    captured = _captured([_kestrel_row("S-C", _C + _S)])
    row = _row_for_candidate(captured.candidates[0], Flag="Context_Doubt, Technical_Artifact")
    row[CAPTURE_ORDINAL_COLUMN] = "0"

    evidenced = with_candidate_evidence(captured, [row])

    assert evidenced.candidates[0].flags == frozenset({"Context_Doubt", "Technical_Artifact"})


@pytest.mark.parametrize("flag", [",Context_Doubt", "Context_Doubt,", "Context_Doubt,,Technical_Artifact"])
def test_malformed_comma_joined_flags_fail_closed(flag: str) -> None:
    """Empty flag elements cannot silently disappear during candidate capture."""
    with pytest.raises(ValueError, match="empty element"):
        _captured([_kestrel_row("S-C", _C + _S, flag=flag)])


@pytest.mark.parametrize(
    "flag",
    [
        "Context_Doubt, Not flagged",
        "Not applicable, Context_Doubt",
        "Context_Doubt, Context_Doubt",
    ],
)
def test_reserved_or_duplicate_comma_joined_flag_tokens_fail_closed(flag: str) -> None:
    """Joined flags cannot mix absence placeholders or collapse duplicate evidence."""
    with pytest.raises(ValueError, match="reserved|duplicate"):
        _captured([_kestrel_row("S-C", _C + _S, flag=flag)])


def test_group_context_divergence_survives_tsv_round_trip() -> None:
    """One divergent equivalent representation persists a conservative group scalar."""
    captured = _captured(
        [_kestrel_row("X-5", _FIVE + _X), _kestrel_row("S-C", _C + _S)],
        _MotifSensitiveComponent(equivalent=True, divergent_motif="S-C"),
    )
    evidence = []
    for candidate in captured.candidates:
        row = _row_for_candidate(candidate)
        row[CAPTURE_ORDINAL_COLUMN] = str(candidate.observation_ordinal)
        evidence.append(row)
    selected = overlay_legacy_projection(with_candidate_evidence(captured, evidence), (0, 1), 0)
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}

    assert cells[GROUP_CONTEXT_COLUMN] == "true"
    assert parse_selected_candidate_cells(cells).group_context_diverges is True


@pytest.mark.parametrize("unresolved", [False, True])
def test_absent_translation_cells_survive_default_string_tsv_loading(tmp_path: Any, unresolved: bool) -> None:
    """Canonical nonempty absence tokens survive pandas' default NA recognition."""
    component = _UnresolvedComponent() if unresolved else REAL_COMPONENT
    captured = _captured([_kestrel_row("S-C", _C + _S)], component)
    row = _row_for_candidate(captured.candidates[0])
    row[CAPTURE_ORDINAL_COLUMN] = "0"
    selected = overlay_legacy_projection(with_candidate_evidence(captured, [row]), (0,), 0)
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}
    path = tmp_path / "identity.tsv"
    pd.DataFrame([cells]).to_csv(path, sep="\t", index=False)

    loaded = pd.read_csv(path, sep="\t", dtype=str).iloc[0].to_dict()
    replayed = parse_selected_candidate_cells(loaded)

    assert loaded[IDENTITY_CAPTURE_COLUMNS[1]] != "" and loaded[IDENTITY_CAPTURE_COLUMNS[3]] != ""
    assert replayed.translation == selected.selected_candidate.observation.translation


def test_replay_rejects_capture_and_selection_raw_key_mismatch() -> None:
    """Selected metadata cannot be combined with another observation's translation."""
    first = _overlay_all(_captured([_kestrel_row("S-C", _C + _S)]))
    second = _overlay_all(_captured([_kestrel_row("A-J", _J + _A, ref="C", alt="CG")]))
    cells = {**candidate_capture_cells(first.selected_candidate), **selected_candidate_cells(second)}

    with pytest.raises(ValueError, match="capture and selected raw keys"):
        parse_selected_candidate_cells(cells)


def test_replay_rejects_capture_and_selection_ordinal_mismatch_for_duplicate_key() -> None:
    """Equal raw keys cannot conceal selection metadata from another observation."""
    captured = _captured([_kestrel_row("S-C", _C + _S), _kestrel_row("S-C", _C + _S, support=7)])
    evidence = []
    for candidate in captured.candidates:
        row = _row_for_candidate(candidate)
        row[CAPTURE_ORDINAL_COLUMN] = str(candidate.observation_ordinal)
        evidence.append(row)
    evidenced = with_candidate_evidence(captured, evidence)
    first = overlay_legacy_projection(evidenced, (0, 1), 0)
    second = overlay_legacy_projection(evidenced, (0, 1), 1)
    cells = {**candidate_capture_cells(first.selected_candidate), **selected_candidate_cells(second)}

    with pytest.raises(ValueError, match="capture and selected observation ordinals"):
        parse_selected_candidate_cells(cells)


@pytest.mark.parametrize(
    "raw_json",
    [
        '{"source":"kestrel","values":["S-C",67,"G","GG"]} ',
        '{"values":["S-C",67,"G","GG"],"source":"kestrel"}',
        '{"source":"kestrel","values":["\\u0053-C",67,"G","GG"]}',
        '{"source":"kestrel","source":"kestrel","values":["S-C",67,"G","GG"]}',
    ],
)
def test_replay_rejects_noncanonical_raw_key_json(raw_json: str) -> None:
    """Equivalent JSON spellings cannot create multiple raw-key encodings."""
    cells = _persisted_cells()
    cells[IDENTITY_CAPTURE_COLUMNS[0]] = raw_json
    cells[IDENTITY_SELECTION_COLUMNS[0]] = raw_json

    with pytest.raises(ValueError, match="canonical"):
        parse_selected_candidate_cells(cells)


def test_replay_rejects_noncanonical_boolean_tokens() -> None:
    """Case variants are not accepted when the codec emits lowercase booleans."""
    cells = _persisted_cells()
    cells[IDENTITY_CAPTURE_COLUMNS[4]] = "TRUE"

    with pytest.raises(ValueError, match="canonical boolean"):
        parse_selected_candidate_cells(cells)


@pytest.mark.parametrize(
    ("column", "value"),
    [
        (IDENTITY_SELECTION_COLUMNS[2], "0"),
        (IDENTITY_SELECTION_COLUMNS[3], '["not_a_legacy_gate"]'),
        (CAPTURE_ORDINAL_COLUMN, "01"),
        (SELECTED_ORDINAL_COLUMN, 0),
    ],
)
def test_replay_rejects_closed_state_and_noncanonical_ordinal_values(column: str, value: object) -> None:
    """Counts, blockers, and ordinals must form a possible canonical selected state."""
    cells: dict[str, object] = dict(_persisted_cells())
    cells[column] = value

    with pytest.raises(ValueError):
        parse_selected_candidate_cells(cells)


def test_unresolved_replay_may_report_other_resolved_hypotheses() -> None:
    """An unresolved selection has zero equivalents but may coexist with resolved hypotheses."""
    captured = _captured([_kestrel_row("S-C", _C + _S)], _UnresolvedComponent())
    row = _row_for_candidate(captured.candidates[0])
    row[CAPTURE_ORDINAL_COLUMN] = "0"
    selected = overlay_legacy_projection(with_candidate_evidence(captured, [row]), (0,), 0)
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}
    cells[IDENTITY_SELECTION_COLUMNS[2]] = "1"

    assert parse_selected_candidate_cells(cells).identity_hypothesis_count == 1


def test_replay_rejects_impossible_context_divergence_states() -> None:
    """Group context must conservatively contain selected context and cannot exist unresolved."""
    unresolved = _captured([_kestrel_row("S-C", _C + _S)], _UnresolvedComponent())
    unresolved_selected = overlay_legacy_projection(
        with_candidate_evidence(unresolved, [_row_for_candidate(unresolved.candidates[0])]),
        (0,),
        0,
    )
    unresolved_cells = {
        **candidate_capture_cells(unresolved_selected.selected_candidate),
        **selected_candidate_cells(unresolved_selected),
    }
    unresolved_cells[GROUP_CONTEXT_COLUMN] = "true"
    with pytest.raises(ValueError, match="unresolved selected identity"):
        parse_selected_candidate_cells(unresolved_cells)

    divergent = _captured(
        [_kestrel_row("S-C", _C + _S)],
        _MotifSensitiveComponent(equivalent=True, divergent_motif="S-C"),
    )
    divergent_selected = overlay_legacy_projection(
        with_candidate_evidence(divergent, [_row_for_candidate(divergent.candidates[0])]),
        (0,),
        0,
    )
    divergent_cells = {
        **candidate_capture_cells(divergent_selected.selected_candidate),
        **selected_candidate_cells(divergent_selected),
    }
    divergent_cells[GROUP_CONTEXT_COLUMN] = "false"
    with pytest.raises(ValueError, match="Selected context divergence"):
        parse_selected_candidate_cells(divergent_cells)


def test_replay_rejects_unresolved_selected_evidence_with_blocking_gates() -> None:
    """An unresolved selected observation cannot carry a group-level gate union."""
    captured = _captured([_kestrel_row("S-C", _C + _S)], _UnresolvedComponent())
    selected = _overlay_all(captured)
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}
    cells[IDENTITY_SELECTION_COLUMNS[3]] = '["motif_filter_pass"]'

    with pytest.raises(ValueError, match="unresolved selected identity.*blocking"):
        parse_selected_candidate_cells(cells)


def test_replay_rejects_resolved_singleton_evidence_with_blocking_gates() -> None:
    """A selected singleton passed all gates, so its group cannot report a blocker."""
    cells = _persisted_cells()
    cells[IDENTITY_SELECTION_COLUMNS[3]] = '["motif_filter_pass"]'

    with pytest.raises(ValueError, match="singleton.*blocking"):
        parse_selected_candidate_cells(cells)


def test_replay_rejects_resolved_singleton_group_context_mismatch() -> None:
    """A singleton group has exactly the selected observation's context state."""
    captured = _captured(
        [_kestrel_row("S-C", _C + _S)],
        _MotifSensitiveComponent(equivalent=True),
    )
    selected = _overlay_all(captured)
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}
    assert cells[IDENTITY_CAPTURE_COLUMNS[4]] == "false"
    cells[GROUP_CONTEXT_COLUMN] = "true"

    with pytest.raises(ValueError, match="singleton.*context"):
        parse_selected_candidate_cells(cells)


@pytest.mark.parametrize("unresolved", [False, True])
def test_selected_candidate_cells_reject_impossible_selected_blockers(unresolved: bool) -> None:
    """The producer cannot serialize blockers for unresolved or singleton selections."""
    component = _UnresolvedComponent() if unresolved else REAL_COMPONENT
    captured = _captured([_kestrel_row("S-C", _C + _S)], component)
    candidate = replace(
        captured.candidates[0],
        blocking_gates=frozenset({"motif_filter_pass"}),
        eligible=True,
    )
    selected = IdentityCandidateSet("kestrel", (candidate,), candidate.observation_ordinal)

    with pytest.raises(ValueError, match="unresolved selected identity|singleton"):
        selected_candidate_cells(selected)


def test_multi_representation_group_may_union_blockers_and_context_divergence() -> None:
    """Conservative blocker and context unions remain valid for equivalent representations."""
    captured = _captured(
        [_kestrel_row("X-5", _FIVE + _X), _kestrel_row("S-C", _C + _S)],
        _MotifSensitiveComponent(equivalent=True, divergent_motif="S-C"),
    )
    evidence = [
        _row_for_candidate(captured.candidates[0]),
        _row_for_candidate(captured.candidates[1], motif_filter_pass=False),
    ]
    selected = overlay_legacy_projection(with_candidate_evidence(captured, evidence), (0,), 0)
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}

    replayed = parse_selected_candidate_cells(cells)

    assert replayed.equivalent_representation_count == 2
    assert replayed.blocking_gates == frozenset({"motif_filter_pass"})
    assert replayed.translation.context_diverges is False
    assert replayed.group_context_diverges is True


def test_resolved_singleton_may_retain_nonblocking_flags() -> None:
    """Flags remain dynamic annotations and do not imply an impossible gate blocker."""
    captured = _captured([_kestrel_row("S-C", _C + _S, flag="Context_Doubt")])
    selected = _overlay_all(captured)
    cells = {**candidate_capture_cells(selected.selected_candidate), **selected_candidate_cells(selected)}

    replayed = parse_selected_candidate_cells(cells)

    assert replayed.equivalent_representation_count == 1
    assert replayed.blocking_gates == frozenset()
    assert replayed.flags == frozenset({"Context_Doubt"})
