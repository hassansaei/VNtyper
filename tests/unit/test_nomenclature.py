"""Round trip: simulator edit -> translator name.

The MUC1 simulator works in the same coding frame as the MUC1 literature, so one of
its mutation definitions *is* a reference edit on a 60 bp repeat unit. Applying a
definition and naming the result must reproduce the name that definition describes,
for every repeat unit the definition is allowed on.

The benchmark tree is supplied out of band via ``VNTYPER_SIM_ROOT`` so that no path
to it is committed; without it these tests skip.

Research use only.
"""

import json
import os
from pathlib import Path

import pytest

from vntyper.scripts.nomenclature import (
    CANONICAL_UNIT,
    ambiguity_interval,
    name_edit,
    normalise,
    repeat_form,
)

pytestmark = pytest.mark.unit

SIM_ROOT = os.environ.get("VNTYPER_SIM_ROOT")

#: Expected name per (mutation, repeat symbol), derived from the simulator's own
#: definitions. ``insG_pos54`` is the reason this is keyed by motif and not by
#: mutation alone: one edit takes two different correct names depending on the unit
#: it landed on.
EXPECTED: dict[tuple[str, str], str] = {
    ("dupC", "X"): "59dupC",
    ("insCCCC", "X"): "56_59dupCCCC",
    ("insG", "X"): "58_59insG",
    ("dupA", "X"): "60dupA",
    ("delinsAT", "X"): "54_56delinsAT",
    ("ins25bp", "X"): "30_31insCAGGCCGGCCCCGGGCTCCGGACAC",
    ("insC_pos23", "A"): "23dupC",
    ("insC_pos23", "E"): "23dupC",
    ("insG_pos58", "B"): "57_58insG",
    ("insG_pos58", "X"): "57_58insG",
    ("insG_pos54", "B"): "53_54insG",
    ("insG_pos54", "J"): "54dupG",
    ("insA_pos54", "A"): "53_54insA",
    ("insA_pos54", "H"): "53_54insA",
}


def _config() -> dict:
    if not SIM_ROOT:
        pytest.skip("VNTYPER_SIM_ROOT unset; simulator config unavailable")
    path = Path(SIM_ROOT) / "muconeup_config_experiment.snapshot.json"
    if not path.exists():
        pytest.skip(f"simulator config not found at {path}")
    return json.loads(path.read_text())


def _cases() -> list[tuple[str, str, dict]]:
    """Every (mutation, repeat symbol) pair the simulator can produce."""
    if not SIM_ROOT:
        return []
    path = Path(SIM_ROOT) / "muconeup_config_experiment.snapshot.json"
    if not path.exists():
        return []
    data = json.loads(path.read_text())
    return [
        (mutation, symbol, definition)
        for mutation, definition in data["mutations"].items()
        for symbol in definition["allowed_repeats"]
    ]


def _edit(change: dict) -> tuple[int, int, str]:
    """Turn one simulator change into a 1-based inclusive edit.

    An insertion is expressed as ``end == start - 1`` -- the empty span between two
    bases -- so that insertions and deletions share one coordinate convention.
    """
    kind = change["type"]
    if kind == "insert":
        return change["start"], change["start"] - 1, change["sequence"]
    if kind == "delete":
        return change["start"], change["end"], ""
    if kind in ("delete_insert", "replace"):
        return change["start"], change["end"], change["sequence"]
    raise AssertionError(f"unhandled simulator change type {kind!r}")


def test_canonical_unit_is_the_simulators_unit() -> None:
    """The translator's canonical X must be the simulator's X, byte for byte."""
    data = _config()
    assert CANONICAL_UNIT == data["repeats"]["X"]


def test_canonical_unit_carries_the_published_c_tract() -> None:
    """Wenzel's tract: 7xC at 53-59, A at 60."""
    assert CANONICAL_UNIT[52:59] == "C" * 7
    assert CANONICAL_UNIT[59] == "A"


def test_the_round_trip_covers_exactly_fifty_pairs() -> None:
    """19 mutation definitions x their allowed repeats = 50. This is the 50/50."""
    _config()
    assert len(_cases()) == 50


@pytest.mark.parametrize(
    "mutation,symbol,definition",
    _cases(),
    ids=lambda value: value if isinstance(value, str) else "",
)
def test_round_trip_names_every_simulator_edit(mutation: str, symbol: str, definition: dict) -> None:
    """Every simulator edit must yield a name, on every unit it is allowed on."""
    data = _config()
    unit = data["repeats"][symbol]

    (change,) = definition["changes"]
    start, end, inserted = _edit(change)
    name = name_edit(unit, start, end, inserted)

    assert name, f"{mutation}/{symbol}: produced no name"
    assert not name.startswith("c."), "spec §2.2 forbids a c. prefix"
    assert name[0].isdigit(), f"{mutation}/{symbol}: name must open with a position, got {name!r}"

    expected = EXPECTED.get((mutation, symbol))
    if expected is not None:
        assert name == expected, f"{mutation}/{symbol}"


# ---------------------------------------------------------------------------
# Normalisation invariants
# ---------------------------------------------------------------------------


def test_every_anchor_inside_the_window_yields_one_name() -> None:
    """Anchoring anywhere in the C-tract must give the same name.

    This is the whole point of normalising: 59dupC, 53dupC and the older 27dupC
    are one event, not three.
    """
    names = {name_edit(CANONICAL_UNIT, anchor, anchor - 1, "C") for anchor in range(53, 61)}
    assert names == {"59dupC"}


def test_normalisation_is_idempotent() -> None:
    """Normalising an already-normalised edit must not move it again."""
    once = normalise(CANONICAL_UNIT, 55, 54, "C")
    twice = normalise(CANONICAL_UNIT, *once)
    assert once == twice


def test_a_deletion_anywhere_in_the_tract_yields_one_name() -> None:
    """The same invariant for deletions, which roll by a different rule."""
    names = {name_edit(CANONICAL_UNIT, anchor, anchor, "") for anchor in range(53, 60)}
    assert names == {"59delC"}


def test_junction_spanning_never_yields_a_coordinate_below_one() -> None:
    """A span at the 5' edge stays on the unit; no negative or zero coordinate."""
    name = name_edit(CANONICAL_UNIT, 1, 5, "")
    assert name == "1_5delGCCCA"
    assert "-" not in name
    assert not name.startswith("0")


# ---------------------------------------------------------------------------
# Ambiguity interval -- nullable, emitted only when wider than 1 bp
# ---------------------------------------------------------------------------


def test_ambiguity_interval_spans_the_whole_c_tract() -> None:
    """dupC is ambiguous across the 7xC tract: every anchor in 53-59 is the allele."""
    assert ambiguity_interval(CANONICAL_UNIT, 59, 58, "C") == (53, 59)


def test_ambiguity_interval_is_none_when_the_edit_cannot_shift() -> None:
    """insG into an all-C tract is unambiguous: G never equals C, so it cannot roll."""
    assert ambiguity_interval(CANONICAL_UNIT, 59, 58, "G") is None


def test_ambiguity_interval_is_none_for_a_delins() -> None:
    """A delins is anchored by definition, so it has no ambiguity window."""
    assert ambiguity_interval(CANONICAL_UNIT, 54, 56, "AT") is None


# ---------------------------------------------------------------------------
# Repeat form -- nullable, emitted only inside a detectable tract
# ---------------------------------------------------------------------------


def test_repeat_form_counts_the_tract_before_and_after() -> None:
    """States what was measured -- the tract went 7 -> 8 -- not which base was added."""
    assert repeat_form(CANONICAL_UNIT, 59, 58, "C") == "53C[7]>53C[8]"


def test_repeat_form_scales_to_a_longer_insertion() -> None:
    """insCCCC is far clearer as 53C[11] than as 56_59dupCCCC."""
    assert repeat_form(CANONICAL_UNIT, 59, 58, "CCCC") == "53C[7]>53C[11]"


def test_repeat_form_counts_down_for_a_deletion() -> None:
    assert repeat_form(CANONICAL_UNIT, 59, 59, "") == "53C[7]>53C[6]"


def test_repeat_form_is_none_outside_any_tract() -> None:
    """Position 2 sits in GCCCACGG..., no homopolymer run to report."""
    assert repeat_form(CANONICAL_UNIT, 2, 2, "A") is None
