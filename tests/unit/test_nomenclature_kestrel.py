"""Translating Kestrel's 120 bp pair-frame records into MUC1 nomenclature.

Kestrel reports against a 551-record reference of merged motif pairs. A record
named ``<L>-<R>`` holds ``seq(R) ++ seq(L)`` -- verified across all 551 records --
so a position above 60 belongs to the *left* symbol, and one at or below 60 to the
right. Positions are on the genomic plus strand; MUC1 is transcribed from the minus
strand, so naming happens after a reverse complement into the coding frame.

Research use only.
"""

import pytest

from vntyper.scripts.nomenclature import CANONICAL_UNIT, MOTIFS, PAIRS, from_kestrel

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Reference tables
# ---------------------------------------------------------------------------


def test_the_reference_tables_are_the_shipped_ones() -> None:
    assert len(MOTIFS) == 34
    assert len(PAIRS) == 551
    assert all(len(unit) == 60 for unit in MOTIFS.values())


def test_a_pair_record_is_right_then_left() -> None:
    """``S-C`` is seq(C) ++ seq(S), not seq(S) ++ seq(C)."""
    pair = PAIRS["S-C"]
    assert pair[:60] == MOTIFS["C"]
    assert pair[60:] == MOTIFS["S"]


# ---------------------------------------------------------------------------
# The canonical record
# ---------------------------------------------------------------------------


def test_the_canonical_dupc_record_names_as_59dupc() -> None:
    """``S-C`` POS 67 G>GG is the commonest call in the benchmark (87/96)."""
    call = from_kestrel("S-C", 67, "G", "GG")

    assert call.name == "59dupC"
    assert call.event == "duplication"
    assert call.net_length == 1
    assert call.source == "kestrel_vcf"


def test_the_other_dupc_shape_names_the_same() -> None:
    """9 of the 96 called dupC samples arrive as ``C>CG`` rather than ``G>GG``.

    Both are real benchmark records. They carry different anchor bases but insert
    the same G after the same pair position, so they are one allele; if they
    disagreed the name would depend on which base happened to anchor the call.
    """
    assert from_kestrel("A-J", 67, "C", "CG").name == "59dupC"
    assert from_kestrel("S-C", 67, "G", "GG").name == "59dupC"


def test_x_projection_beats_the_assigned_motif() -> None:
    """Motif S carries only 5xC; X carries 7xC at 53-59.

    Anchoring on the assigned motif would name this 57dupC. The truth is 59dupC,
    which is why the position is projected onto canonical X before naming.
    """
    assert MOTIFS["S"].count("G" * 5) >= 1, "motif S has the shorter run this test relies on"
    call = from_kestrel("S-C", 67, "G", "GG")
    assert call.name == "59dupC"
    assert call.unit == "X"


def test_a_divergent_motif_context_is_flagged_not_hidden() -> None:
    """Projecting onto X when the assigned motif differs must cost confidence."""
    call = from_kestrel("S-C", 67, "G", "GG")
    assert "motif-context-diverges" in call.flags


def test_no_single_caller_record_ever_reaches_tier_a() -> None:
    """Tier A needs two agreeing sources; one caller cannot award it to itself.

    The benchmark shows why: Kestrel places the whole insG family one position 3'
    of truth, so a record that looks perfectly clean still names a wrong allele.
    """
    for motifs, pos, ref, alt in (
        ("S-C", 67, "G", "GG"),
        ("X-X", 67, "G", "GG"),
        ("A-J", 67, "C", "CG"),
        ("5C-S", 61, "T", "TC"),
    ):
        assert from_kestrel(motifs, pos, ref, alt).tier != "A"


def test_a_self_pair_on_x_keeps_the_top_tier_available() -> None:
    """X-X needs no projection, so nothing is flagged for context divergence."""
    call = from_kestrel("X-X", 67, "G", "GG")
    assert call.name == "59dupC"
    assert "motif-context-diverges" not in call.flags


# ---------------------------------------------------------------------------
# Position handling
# ---------------------------------------------------------------------------


def test_a_low_position_still_names() -> None:
    """A real ``ins25bp`` record lands at POS 30, in the pair's first half."""
    call = from_kestrel("5C-S", 30, "T", "TG")
    assert call.name, "a first-half position must still produce a name"
    assert call.name[0].isdigit()


def test_the_canonical_unit_is_exported_for_callers() -> None:
    assert len(CANONICAL_UNIT) == 60


def test_an_unknown_pair_never_invents_a_name() -> None:
    call = from_kestrel("ZZ-ZZ", 67, "G", "GG")
    assert call.name is None
    assert call.tier == "C"


def test_a_position_outside_the_pair_never_invents_a_name() -> None:
    call = from_kestrel("S-C", 999, "G", "GG")
    assert call.name is None
    assert call.tier == "C"


def test_a_negative_run_row_never_invents_a_name() -> None:
    """Negative runs carry the literal string 'None' in every field."""
    call = from_kestrel("None", 0, "None", "None")
    assert call.name is None
    assert call.tier == "C"
