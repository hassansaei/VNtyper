"""Outcome-independent caller populations and exact observation binding."""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import NoReturn

from vntyper.scripts.calibration_caller_metrics import CallerObservation, validate_caller_observations
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class EligibleCallerMember:
    """One predeclared independent representative with overlapping metadata strata."""

    key: str
    group_key: str
    strata: tuple[str, ...]


@dataclass(frozen=True)
class CallerEligibleRoster:
    """Immutable population frozen before observing candidate performance."""

    members: tuple[EligibleCallerMember, ...]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _text(value: object) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        _fail("caller roster identifiers must be non-empty trimmed strings")
    return value


def decode_caller_eligible_roster(value: object) -> CallerEligibleRoster:
    """Decode exact metadata-only members before predictions are collected.

    Args:
        value: Closed rows with key, group_key and sorted unique strata.

    Returns:
        One representative per group and the canonical roster digest.

    Raises:
        ValueError: If content is empty, malformed, unordered or duplicated.
    """
    if not isinstance(value, list) or not value:
        _fail("caller eligible roster must be a non-empty list")
    members: list[EligibleCallerMember] = []
    for row in value:
        if not isinstance(row, Mapping) or set(row) != {"key", "group_key", "strata"}:
            _fail("caller eligible roster row fields differ from the closed contract")
        strata_raw = row["strata"]
        if not isinstance(strata_raw, list) or not strata_raw:
            _fail("caller roster strata must be a non-empty list")
        strata = tuple(_text(item) for item in strata_raw)
        if strata != tuple(sorted(set(strata))):
            _fail("caller roster strata must be sorted and unique")
        members.append(EligibleCallerMember(_text(row["key"]), _text(row["group_key"]), strata))
    if (
        len({member.key for member in members}) != len(members)
        or len({member.group_key for member in members}) != len(members)
        or members != sorted(members, key=lambda member: member.group_key)
    ):
        _fail("caller roster requires unique specimen/group keys sorted by group")
    return CallerEligibleRoster(tuple(members), canonical_sha256(value))


def caller_eligible_roster_document(roster: CallerEligibleRoster) -> list[dict[str, object]]:
    """Project a roster after validating immutable content and its digest.

    Args:
        roster: Previously decoded metadata-only population.

    Returns:
        Fresh JSON-compatible canonical roster rows.

    Raises:
        ValueError: If typed content is mutable, forged or inconsistent.
    """
    if (
        not isinstance(roster, CallerEligibleRoster)
        or not isinstance(roster.members, tuple)
        or any(not isinstance(member, EligibleCallerMember) for member in roster.members)
        or any(not isinstance(member.strata, tuple) for member in roster.members)
    ):
        _fail("caller eligible roster requires decoded immutable content")
    raw: list[dict[str, object]] = [
        {"key": member.key, "group_key": member.group_key, "strata": list(member.strata)} for member in roster.members
    ]
    if decode_caller_eligible_roster(raw) != roster:
        _fail("caller eligible roster differs from its canonical content or digest")
    return raw


def bind_caller_observations(
    rows: Sequence[CallerObservation], roster: CallerEligibleRoster
) -> tuple[CallerObservation, ...]:
    """Reject omitted, extra or swapped outcomes before caller metric evaluation.

    Args:
        rows: Complete candidate outcomes, including unavailable calls.
        roster: Frozen independent primary representatives and strata.

    Returns:
        Validated complete outcomes in deterministic group order.

    Raises:
        ValueError: If any member is missing, additional, duplicated or swapped.
    """
    caller_eligible_roster_document(roster)
    observations = validate_caller_observations(rows)
    if {(row.key, row.group_key) for row in observations} != {
        (member.key, member.group_key) for member in roster.members
    }:
        _fail("caller outcome set must match the frozen eligible roster exactly")
    return observations
