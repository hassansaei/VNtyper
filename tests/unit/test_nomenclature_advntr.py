"""Translating adVNTR HMM states into MUC1 nomenclature.

adVNTR reports states like ``I22_2_G_LEN1`` (insertion after HMM match 22 of repeat
unit 2, first inserted base G, one base long) and ``D17_2`` (deletion of match 17 in
unit 2), joined with ``&`` into compound states.

Two properties of the caller constrain what may be named, both read from the pinned
adVNTR source rather than assumed:

* ``vntr_finder.py`` records only the *first* inserted base -- "If there are run of
  insertions, the sequence might differ, but we just take the first base". So
  ``LEN1`` carries a complete allele and ``LEN>1`` carries a length only.
* Components are joined only when consecutive within one repeat unit, or when a
  deletion and insertion share a position within one repeat unit.

Every state in these tests was emitted by a real run over the benchmark.

Research use only.
"""

import pytest

from vntyper.scripts.nomenclature import MAPPABLE_RUS, from_advntr

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Which repeat units carry a canonical coordinate at all
# ---------------------------------------------------------------------------


def test_only_the_rotation_units_are_mappable() -> None:
    """Of adVNTR's nine RUs only three are rotations of a known motif.

    RU1, RU4 and RU8 are 60 bp but match no motif; RU3, RU7 and RU9 are 78, 48 and
    54 bp. Projecting any of them through the RU2 map would fabricate a coordinate.
    """
    assert set(MAPPABLE_RUS) == {2, 5, 6}


def test_a_state_in_an_unmappable_unit_never_gets_a_number() -> None:
    """``ins25bp`` really does emit RU3 and RU7 states."""
    for state in ("I49_3_C_LEN5", "I69_3_G_LEN1", "D49_7"):
        calls = from_advntr(state)
        assert all(call.name is None for call in calls), state
        assert all(call.tier == "C" for call in calls), state


# ---------------------------------------------------------------------------
# Single states
# ---------------------------------------------------------------------------


def test_the_canonical_dupc_state() -> None:
    """``I22_2_G_LEN1`` is adVNTR's commonest call on the dupC cohort (79/100)."""
    (call,) = from_advntr("I22_2_G_LEN1")
    assert call.name == "59dupC"
    assert call.source == "advntr"
    assert call.net_length == 1


def test_the_insg_state_separates_from_dupc() -> None:
    """adVNTR distinguishes what Kestrel collapses: I22, I23 and I24 are distinct."""
    assert from_advntr("I23_2_C_LEN1")[0].name == "58_59insG"
    assert from_advntr("I24_2_C_LEN1")[0].name == "57_58insG"
    assert from_advntr("I28_2_C_LEN1")[0].name == "53_54insG"


def test_a_multi_base_insertion_states_its_length_but_never_its_sequence() -> None:
    """``LEN4`` means "four bases, the first was G" -- not GGGG.

    Naming it 56_59dupCCCC would assert three bases the caller never reported.
    """
    (call,) = from_advntr("I22_2_G_LEN4")
    assert call.name is None
    assert call.net_length == 4
    assert "sequence-undetermined" in call.flags
    assert call.tier != "A"


def test_a_state_in_a_non_x_unit_is_flagged() -> None:
    """RU6 is a rotation of motif C, not X, so its projection is not free."""
    (call,) = from_advntr("I23_6_G_LEN1")
    assert "motif-context-diverges" in call.flags


# ---------------------------------------------------------------------------
# Compound states
# ---------------------------------------------------------------------------


def test_consecutive_deletions_in_one_unit_merge_into_one_event() -> None:
    """``delGCCCA`` arrives as five single-base deletions; it is one 5 bp deletion."""
    calls = from_advntr("D17_2&D18_2&D19_2&D20_2&D21_2")
    assert len(calls) == 1
    assert calls[0].event == "deletion"
    assert calls[0].net_length == -5
    assert calls[0].name == "1_5delGCCCA"


def test_a_deletion_and_insertion_at_one_position_merge_into_a_delins() -> None:
    """The event Kestrel structurally cannot represent.

    The *shape* is a delins; the inserted sequence is not knowable from LEN2, so
    no bare number is emitted.
    """
    calls = from_advntr("D27_2&I27_2_A_LEN2")
    assert len(calls) == 1
    assert calls[0].event == "delins"
    assert "sequence-undetermined" in calls[0].flags
    assert calls[0].name is None


def test_components_in_different_units_are_never_merged() -> None:
    """adVNTR only joins components sharing a repeat unit; merging across units
    would fuse two independent events into one fabricated allele."""
    calls = from_advntr("D17_2&D18_6")
    assert len(calls) == 2


def test_non_consecutive_deletions_are_never_merged() -> None:
    calls = from_advntr("D17_2&D20_2")
    assert len(calls) == 2


# ---------------------------------------------------------------------------
# Non-states
# ---------------------------------------------------------------------------


def test_a_negative_run_yields_no_call() -> None:
    """All 200 normal controls returned this."""
    assert from_advntr("Not applicable") == ()


def test_unparseable_input_yields_no_call_rather_than_a_guess() -> None:
    for state in ("", "   ", "banana", "I_2_G_LEN1", "Q22_2_G_LEN1"):
        assert from_advntr(state) == (), state
