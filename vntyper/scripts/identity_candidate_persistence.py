"""Pure TSV codecs for retained identity-candidate metadata."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

from vntyper.scripts.identity_candidates import (
    LEGACY_GATE_COLUMNS,
    OBSERVATION_ORDINAL_COLUMN,
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
    OBSERVATION_ORDINAL_COLUMN,
)

IDENTITY_SELECTION_COLUMNS: tuple[str, ...] = (
    "__Identity_Selected_Raw_Representation_Key",
    "__Identity_Equivalent_Representation_Count",
    "__Identity_Hypothesis_Count",
    "__Identity_Group_Blocking_Gates",
    "__Identity_Group_Flags",
    "__Identity_Selected_Observation_Ordinal",
    "__Identity_Group_Context_Diverges",
)

ABSENT_TOKEN = "absent"


@dataclass(frozen=True)
class PersistedIdentityCandidate:
    """Selected candidate evidence decoded from TSV-persistent scalar columns."""

    translation: IdentityTranslation
    selected_row_key: RawRepresentationKey
    selected_observation_ordinal: int
    equivalent_representation_count: int
    identity_hypothesis_count: int
    blocking_gates: frozenset[str]
    flags: frozenset[str]
    group_context_diverges: bool


@dataclass(frozen=True)
class PersistedIdentityCapture:
    """One complete identity capture decoded from a pre-result TSV row."""

    translation: IdentityTranslation
    row_key: RawRepresentationKey
    observation_ordinal: int


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
        IDENTITY_CAPTURE_COLUMNS[1]: (
            serialize_molecular_identity(translation.identity) if translation.identity else ABSENT_TOKEN
        ),
        IDENTITY_CAPTURE_COLUMNS[2]: translation.status,
        IDENTITY_CAPTURE_COLUMNS[3]: translation.failure or ABSENT_TOKEN,
        IDENTITY_CAPTURE_COLUMNS[4]: "true" if translation.context_diverges else "false",
        IDENTITY_CAPTURE_COLUMNS[5]: str(candidate.observation_ordinal),
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
        group_context_diverges = False
    else:
        group = candidates.selected_group
        equivalent_count = len(group.candidates)
        blocking_gates = group.blocking_gates
        flags = group.flags
        group_context_diverges = group.context_diverges
    _validate_group_invariants(
        selected.observation.translation,
        equivalent_count,
        blocking_gates,
        group_context_diverges,
    )
    return {
        IDENTITY_SELECTION_COLUMNS[0]: _serialize_raw_key(selected.row_key),
        IDENTITY_SELECTION_COLUMNS[1]: str(equivalent_count),
        IDENTITY_SELECTION_COLUMNS[2]: str(candidates.identity_hypothesis_count),
        IDENTITY_SELECTION_COLUMNS[3]: _serialize_strings(blocking_gates),
        IDENTITY_SELECTION_COLUMNS[4]: _serialize_strings(flags),
        IDENTITY_SELECTION_COLUMNS[5]: str(selected.observation_ordinal),
        IDENTITY_SELECTION_COLUMNS[6]: "true" if group_context_diverges else "false",
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
    capture_raw_key = _parse_raw_key(_required_string(row, IDENTITY_CAPTURE_COLUMNS[0]))
    selected_raw_key = _parse_raw_key(_required_string(row, IDENTITY_SELECTION_COLUMNS[0]))
    if capture_raw_key != selected_raw_key:
        raise ValueError("Identity capture and selected raw keys must match exactly")
    capture_ordinal = _nonnegative_int(row, IDENTITY_CAPTURE_COLUMNS[5])
    selected_ordinal = _nonnegative_int(row, IDENTITY_SELECTION_COLUMNS[5])
    if capture_ordinal != selected_ordinal:
        raise ValueError("Identity capture and selected observation ordinals must match exactly")
    identity_text = _required_string(row, IDENTITY_CAPTURE_COLUMNS[1])
    status = _required_string(row, IDENTITY_CAPTURE_COLUMNS[2])
    failure_text = _required_string(row, IDENTITY_CAPTURE_COLUMNS[3])
    identity = None if identity_text == ABSENT_TOKEN else parse_molecular_identity(identity_text)
    failure = None if failure_text == ABSENT_TOKEN else cast(TranslationFailure, failure_text)
    translation = IdentityTranslation(
        identity,
        cast(TranslationStatus, status),
        failure,
        _parse_bool(row[IDENTITY_CAPTURE_COLUMNS[4]]),
    )
    equivalent_count = _nonnegative_int(row, IDENTITY_SELECTION_COLUMNS[1])
    hypothesis_count = _nonnegative_int(row, IDENTITY_SELECTION_COLUMNS[2])
    blocking_gates = _parse_strings(_required_string(row, IDENTITY_SELECTION_COLUMNS[3]))
    if not blocking_gates <= frozenset(LEGACY_GATE_COLUMNS):
        raise ValueError("Persisted candidate blockers must use the six legacy gate names")
    group_context_diverges = _parse_bool(row[IDENTITY_SELECTION_COLUMNS[6]])
    if translation.identity is None:
        if equivalent_count != 0:
            raise ValueError("An unresolved selected identity requires zero equivalent representations")
        if group_context_diverges:
            raise ValueError("An unresolved selected identity cannot have divergent group context")
    elif equivalent_count < 1 or hypothesis_count < 1:
        raise ValueError("A resolved selected identity requires positive equivalent and hypothesis counts")
    _validate_group_invariants(translation, equivalent_count, blocking_gates, group_context_diverges)
    return PersistedIdentityCandidate(
        translation=translation,
        selected_row_key=selected_raw_key,
        selected_observation_ordinal=selected_ordinal,
        equivalent_representation_count=equivalent_count,
        identity_hypothesis_count=hypothesis_count,
        blocking_gates=blocking_gates,
        flags=_parse_strings(_required_string(row, IDENTITY_SELECTION_COLUMNS[4])),
        group_context_diverges=group_context_diverges,
    )


def parse_candidate_capture_cells(row: Mapping[str, object]) -> PersistedIdentityCapture:
    """Decode one pre-selection identity capture after a TSV round trip.

    Args:
        row: Mapping containing every internal identity-capture column.

    Returns:
        Validated raw-key, translation, and stable observation ordinal.

    Raises:
        KeyError: If a required internal column is absent.
        ValueError: If serialized values are inconsistent or malformed.
    """
    raw_key = _parse_raw_key(_required_string(row, IDENTITY_CAPTURE_COLUMNS[0]))
    identity_text = _required_string(row, IDENTITY_CAPTURE_COLUMNS[1])
    failure_text = _required_string(row, IDENTITY_CAPTURE_COLUMNS[3])
    identity = None if identity_text == ABSENT_TOKEN else parse_molecular_identity(identity_text)
    failure = None if failure_text == ABSENT_TOKEN else cast(TranslationFailure, failure_text)
    translation = IdentityTranslation(
        identity,
        cast(TranslationStatus, _required_string(row, IDENTITY_CAPTURE_COLUMNS[2])),
        failure,
        _parse_bool(row[IDENTITY_CAPTURE_COLUMNS[4]]),
    )
    return PersistedIdentityCapture(
        translation,
        raw_key,
        _nonnegative_int(row, IDENTITY_CAPTURE_COLUMNS[5]),
    )


def _serialize_raw_key(key: RawRepresentationKey) -> str:
    return json.dumps({"source": key.source, "values": key.values}, separators=(",", ":"), sort_keys=True)


def _validate_group_invariants(
    translation: IdentityTranslation,
    equivalent_count: int,
    blocking_gates: frozenset[str],
    group_context_diverges: bool,
) -> None:
    """Reject group evidence that cannot arise from a selected passing observation."""
    if translation.identity is None:
        if blocking_gates:
            raise ValueError("An unresolved selected identity cannot carry blocking gates")
        return
    if equivalent_count == 1:
        if blocking_gates:
            raise ValueError("A resolved singleton identity group cannot carry blocking gates")
        if group_context_diverges != translation.context_diverges:
            raise ValueError("Resolved singleton context must equal Selected context divergence")
    if translation.context_diverges and not group_context_diverges:
        raise ValueError("Selected context divergence must be present in its identity group")


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
    raw_key = RawRepresentationKey(cast(CandidateSource, source), typed_values)
    if _serialize_raw_key(raw_key) != value:
        raise ValueError("Raw representation key must use its canonical JSON serialization")
    return raw_key


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
    result = frozenset(parsed)
    if _serialize_strings(result) != value:
        raise ValueError("Candidate string set must use canonical JSON serialization")
    return result


def _parse_bool(value: object) -> bool:
    if not isinstance(value, str) or value not in {"true", "false"}:
        raise ValueError("Identity context divergence must use a canonical boolean token")
    return value == "true"


def _required_string(row: Mapping[str, object], column: str) -> str:
    value = row[column]
    if not isinstance(value, str) or not value:
        raise ValueError(f"{column} must be a non-empty string")
    return value


def _nonnegative_int(row: Mapping[str, object], column: str) -> int:
    value = row[column]
    if (
        not isinstance(value, str)
        or not value
        or not value.isascii()
        or not value.isdecimal()
        or (len(value) > 1 and value.startswith("0"))
    ):
        raise ValueError(f"{column} must be a canonical non-negative integer")
    return int(value)
