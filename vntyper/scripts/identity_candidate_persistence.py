"""Pure TSV codecs for retained identity-candidate metadata."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral
from typing import cast

from vntyper.scripts.identity_candidates import (
    CandidateSource,
    IdentityCandidate,
    IdentityCandidateSet,
    RawKeyValues,
    RawRepresentationKey,
)
from vntyper.scripts.molecular_identity import (
    IdentityTranslation,
    TranslationFailure,
    TranslationStatus,
    parse_molecular_identity,
    serialize_molecular_identity,
)

IDENTITY_CAPTURE_COLUMNS: tuple[str, ...] = (
    "__Identity_Raw_Representation_Key",
    "__Identity_Molecular_Identity",
    "__Identity_Translation_Status",
    "__Identity_Translation_Failure",
    "__Identity_Context_Diverges",
)

IDENTITY_SELECTION_COLUMNS: tuple[str, ...] = (
    "__Identity_Selected_Raw_Representation_Key",
    "__Identity_Equivalent_Representation_Count",
    "__Identity_Hypothesis_Count",
    "__Identity_Group_Blocking_Gates",
    "__Identity_Group_Flags",
)


@dataclass(frozen=True)
class PersistedIdentityCandidate:
    """Selected candidate evidence decoded from TSV-persistent scalar columns."""

    translation: IdentityTranslation
    selected_row_key: RawRepresentationKey
    equivalent_representation_count: int
    identity_hypothesis_count: int
    blocking_gates: frozenset[str]
    flags: frozenset[str]


def candidate_capture_cells(candidate: IdentityCandidate) -> dict[str, str]:
    """Serialize one captured translation into internal scalar columns.

    Args:
        candidate: Independently captured caller representation.

    Returns:
        Internal columns carrying raw key and closed translation state.
    """
    translation = candidate.observation.translation
    return {
        IDENTITY_CAPTURE_COLUMNS[0]: _serialize_raw_key(candidate.row_key),
        IDENTITY_CAPTURE_COLUMNS[1]: serialize_molecular_identity(translation.identity) if translation.identity else "",
        IDENTITY_CAPTURE_COLUMNS[2]: translation.status,
        IDENTITY_CAPTURE_COLUMNS[3]: translation.failure or "",
        IDENTITY_CAPTURE_COLUMNS[4]: "true" if translation.context_diverges else "false",
    }


def selected_candidate_cells(candidates: IdentityCandidateSet) -> dict[str, str]:
    """Serialize selected-row counts and conservative group evidence.

    Args:
        candidates: Candidate set carrying an exact legacy selection.

    Returns:
        Internal scalar columns sufficient for later reconciliation replay.
    """
    selected = candidates.selected_candidate
    if selected.identity is None:
        equivalent_count = 0
        blocking_gates = selected.blocking_gates
        flags = selected.flags
    else:
        group = candidates.selected_group
        equivalent_count = len(group.candidates)
        blocking_gates = group.blocking_gates
        flags = group.flags
    return {
        IDENTITY_SELECTION_COLUMNS[0]: _serialize_raw_key(selected.row_key),
        IDENTITY_SELECTION_COLUMNS[1]: str(equivalent_count),
        IDENTITY_SELECTION_COLUMNS[2]: str(candidates.identity_hypothesis_count),
        IDENTITY_SELECTION_COLUMNS[3]: _serialize_strings(blocking_gates),
        IDENTITY_SELECTION_COLUMNS[4]: _serialize_strings(flags),
    }


def parse_selected_candidate_cells(row: Mapping[str, object]) -> PersistedIdentityCandidate:
    """Decode selected candidate metadata after a TSV round trip.

    Args:
        row: Mapping containing every internal capture and selection column.

    Returns:
        Validated selected identity evidence for reconciliation.

    Raises:
        KeyError: If a required internal column is absent.
        ValueError: If serialized values are inconsistent or malformed.
    """
    identity_text = _required_string(row, IDENTITY_CAPTURE_COLUMNS[1], allow_empty=True)
    status = _required_string(row, IDENTITY_CAPTURE_COLUMNS[2])
    failure_text = _required_string(row, IDENTITY_CAPTURE_COLUMNS[3], allow_empty=True)
    identity = parse_molecular_identity(identity_text) if identity_text else None
    translation = IdentityTranslation(
        identity,
        cast(TranslationStatus, status),
        cast(TranslationFailure | None, failure_text or None),
        _parse_bool(row[IDENTITY_CAPTURE_COLUMNS[4]]),
    )
    equivalent_count = _nonnegative_int(row, IDENTITY_SELECTION_COLUMNS[1])
    hypothesis_count = _nonnegative_int(row, IDENTITY_SELECTION_COLUMNS[2])
    if (translation.identity is None) != (equivalent_count == 0):
        raise ValueError("Equivalent representation count must be zero exactly when identity is unresolved")
    return PersistedIdentityCandidate(
        translation=translation,
        selected_row_key=_parse_raw_key(_required_string(row, IDENTITY_SELECTION_COLUMNS[0])),
        equivalent_representation_count=equivalent_count,
        identity_hypothesis_count=hypothesis_count,
        blocking_gates=_parse_strings(_required_string(row, IDENTITY_SELECTION_COLUMNS[3])),
        flags=_parse_strings(_required_string(row, IDENTITY_SELECTION_COLUMNS[4])),
    )


def _serialize_raw_key(key: RawRepresentationKey) -> str:
    return json.dumps({"source": key.source, "values": key.values}, separators=(",", ":"), sort_keys=True)


def _parse_raw_key(value: str) -> RawRepresentationKey:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise ValueError("Raw representation key must be valid JSON") from error
    if not isinstance(parsed, dict) or set(parsed) != {"source", "values"} or not isinstance(parsed["values"], list):
        raise ValueError("Raw representation key serialization has an invalid shape")
    source = parsed["source"]
    values = parsed["values"]
    if source == "advntr" and len(values) == 3 and isinstance(values[1], list) and isinstance(values[2], list):
        typed_values: RawKeyValues = (values[0], tuple(values[1]), tuple(values[2]))
    else:
        typed_values = cast(RawKeyValues, tuple(values))
    return RawRepresentationKey(cast(CandidateSource, source), typed_values)


def _serialize_strings(values: frozenset[str]) -> str:
    return json.dumps(sorted(values), separators=(",", ":"))


def _parse_strings(value: str) -> frozenset[str]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise ValueError("Candidate string set must be valid JSON") from error
    if not isinstance(parsed, list) or any(not isinstance(item, str) or not item for item in parsed):
        raise ValueError("Candidate string set must contain non-empty strings")
    if parsed != sorted(set(parsed)):
        raise ValueError("Candidate string set must use sorted unique canonical form")
    return frozenset(parsed)


def _parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.casefold() in {"true", "false"}:
        return value.casefold() == "true"
    raise ValueError("Identity context divergence must be a boolean scalar")


def _required_string(row: Mapping[str, object], column: str, *, allow_empty: bool = False) -> str:
    value = row[column]
    if not isinstance(value, str) or (not allow_empty and not value):
        raise ValueError(f"{column} must be a {'possibly empty ' if allow_empty else 'non-empty '}string")
    return value


def _nonnegative_int(row: Mapping[str, object], column: str) -> int:
    scalar = row[column]
    if isinstance(scalar, Integral) and not isinstance(scalar, bool):
        if scalar < 0:
            raise ValueError(f"{column} must be a canonical non-negative integer")
        return int(scalar)
    text = _required_string(row, column)
    if not text.isascii() or not text.isdecimal() or (len(text) > 1 and text.startswith("0")):
        raise ValueError(f"{column} must be a canonical non-negative integer")
    return int(text)
