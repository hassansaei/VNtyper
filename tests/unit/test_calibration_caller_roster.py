"""Frozen caller populations cannot lose unavailable observations."""

from dataclasses import replace

import pytest

from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_roster import (
    bind_caller_observations,
    caller_eligible_roster_document,
    decode_caller_eligible_roster,
)

pytestmark = pytest.mark.unit


def _raw():
    return [
        {"key": "a", "group_key": "g1", "strata": ["assay-a", "nominal"]},
        {"key": "b", "group_key": "g2", "strata": ["nominal"]},
    ]


def _rows():
    return (
        CallerObservation("a", "g1", True, ("v",), True, ("v",), ()),
        CallerObservation("b", "g2", True, ("v",), None, (), ()),
    )


def test_roster_roundtrip_is_immutable_and_binding_keeps_the_no_call():
    raw = _raw()
    roster = decode_caller_eligible_roster(raw)
    assert caller_eligible_roster_document(roster) == raw
    raw[0]["strata"].append("other")
    assert roster.members[0].strata == ("assay-a", "nominal")
    assert bind_caller_observations(tuple(reversed(_rows())), roster) == _rows()


@pytest.mark.parametrize("mode", ["omit", "extra", "swap", "duplicate"])
def test_outcome_roster_must_match_exactly(mode):
    rows = _rows()
    if mode == "omit":
        rows = rows[:1]
    elif mode == "extra":
        rows += (replace(rows[0], key="c", group_key="g3"),)
    elif mode == "swap":
        rows = (replace(rows[0], group_key="g2"), replace(rows[1], group_key="g1"))
    else:
        rows += (rows[0],)
    with pytest.raises(ValueError):
        bind_caller_observations(rows, decode_caller_eligible_roster(_raw()))


@pytest.mark.parametrize(
    "mode", ["empty", "fields", "key", "group", "strata-empty", "strata-duplicate", "strata-order", "order"]
)
def test_closed_canonical_roster_rejects_malformed_members(mode):
    raw = _raw()
    if mode == "empty":
        raw = []
    elif mode == "fields":
        raw[0]["unexpected"] = True
    elif mode in {"key", "group"}:
        raw[0]["key" if mode == "key" else "group_key"] = " "
    elif mode == "strata-empty":
        raw[0]["strata"] = []
    elif mode == "strata-duplicate":
        raw[0]["strata"] = ["nominal", "nominal"]
    elif mode == "strata-order":
        raw[0]["strata"] = ["z", "a"]
    else:
        raw.reverse()
    with pytest.raises(ValueError):
        decode_caller_eligible_roster(raw)


@pytest.mark.parametrize("field", ["key", "group_key"])
def test_repeated_representatives_fail(field):
    raw = _raw()
    raw[1][field] = raw[0][field]
    with pytest.raises(ValueError):
        decode_caller_eligible_roster(raw)


def test_typed_roster_integrity_is_rechecked():
    roster = decode_caller_eligible_roster(_raw())
    for bad in (
        None,
        replace(roster, sha256="0" * 64),
        replace(roster, members=list(roster.members)),
        replace(roster, members=(object(),)),
        replace(roster, members=(replace(roster.members[0], strata=["nominal"]),)),
    ):
        with pytest.raises(ValueError):
            caller_eligible_roster_document(bad)
