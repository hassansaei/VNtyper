"""Pure public and diagnostic projections of typed molecular identity."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Any

from vntyper.scripts.molecular_identity import (
    AdvntrRepresentation,
    IdentityTranslation,
    serialize_molecular_identity,
)

if TYPE_CHECKING:
    from vntyper.scripts.identity_candidate_persistence import PersistedIdentityCandidate
    from vntyper.scripts.identity_candidates import IdentityTranslator

IDENTITY_COLUMNS: tuple[str, ...] = (
    "Molecular_Identity",
    "Molecular_Identity_Status",
    "Equivalent_Representation_Count",
    "Identity_Hypothesis_Count",
)

LEGACY_IDENTITY_NOT_RECORDED = "legacy identity not recorded"

IDENTITY_COLUMN_HELP: dict[str, str] = {
    "Molecular_Identity": "Stable molecular identity recorded by the emitting run.",
    "Molecular_Identity_Status": "Whether the recorded identity was unique, selected among alternatives, or unresolved.",
    "Equivalent_Representation_Count": "Caller representations recorded as equivalent to this molecular identity.",
    "Identity_Hypothesis_Count": "Distinct resolved molecular identities considered for this caller result.",
}

IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS: tuple[str, ...] = (
    "Molecular_Identity",
    "Molecular_Identity_Translation_Status",
    "Molecular_Identity_Translation_Failure",
    "Molecular_Identity_Context_Diverges",
)


def identity_compatibility_cells(
    row: Mapping[str, object],
    *,
    schema_version: object,
) -> dict[str, object]:
    """Project recorded or explicit legacy identity cells from a summary row.

    Schema 1 and schema 2 introduced no required identity container. Field presence,
    rather than allele-shaped legacy cells, is therefore the only compatibility
    discriminator. A complete quartet is copied exactly; any absent member makes all
    four cells explicit legacy text. Schema 3 instead requires the complete quartet;
    rendering a malformed current row as legacy would hide a broken run artifact.

    Args:
        row: One caller-positive result row read from a pipeline summary.
        schema_version: Summary schema recorded with the row, or ``None`` for a summary
            written before schema provenance existed.

    Returns:
        The exact four compatibility cells in public column order.
    """
    if all(column in row for column in IDENTITY_COLUMNS):
        return {column: row[column] for column in IDENTITY_COLUMNS}
    if schema_version == 3:
        raise ValueError("summary schema 3 positive row requires the complete molecular identity quartet")
    return dict.fromkeys(IDENTITY_COLUMNS, LEGACY_IDENTITY_NOT_RECORDED)


def identity_compatible_result_row(
    row: Mapping[str, object],
    *,
    schema_version: object,
    positive: bool,
) -> dict[str, object]:
    """Build one downstream result row with identity compatibility presentation.

    Args:
        row: One row from a caller result recorded in a pipeline summary.
        schema_version: Summary schema recorded with the row, or ``None``.
        positive: Whether the caller row is a result rather than its frozen negative
            placeholder.

    Returns:
        A new row. Positive rows end in the exact identity quartet; negative rows are
        copied without widening their caller schema.
    """
    if not positive:
        return dict(row)
    projected = {key: value for key, value in row.items() if key not in IDENTITY_COLUMNS}
    projected.update(identity_compatibility_cells(row, schema_version=schema_version))
    return projected


def identity_result_cells(
    selected: IdentityTranslation,
    *other_positive_translations: IdentityTranslation,
) -> dict[str, str | int]:
    """Project one emitted row against its complete caller-local positive set.

    Args:
        selected: Translation of the emitted row receiving the cells.
        *other_positive_translations: Translations of every other emitted positive
            row from the same caller. Detection-ineligible and below-floor rows must
            be excluded by the caller before this presentation boundary.

    Returns:
        The exact four public molecular-identity cells.

    Raises:
        ValueError: If any input is not a validated identity translation.
    """
    translations = (selected, *other_positive_translations)
    if any(not isinstance(translation, IdentityTranslation) for translation in translations):
        raise ValueError("Identity presentation requires IdentityTranslation values")
    identities = tuple(translation.identity for translation in translations if translation.identity is not None)
    distinct_identities = frozenset(identities)
    if selected.identity is None:
        return _result_cells(selected, 0, len(distinct_identities))
    equivalent_count = sum(identity == selected.identity for identity in identities)
    return _result_cells(selected, equivalent_count, len(distinct_identities))


def persisted_identity_result_cells(candidate: PersistedIdentityCandidate) -> dict[str, str | int]:
    """Project strict selected-row metadata persisted by the Kestrel stage.

    Args:
        candidate: Decoded selected candidate with caller-local aggregate counts.

    Returns:
        The exact four public molecular-identity cells.

    Raises:
        ValueError: If the supplied object or its counts describe an impossible state.
    """
    from vntyper.scripts.identity_candidate_persistence import PersistedIdentityCandidate

    if not isinstance(candidate, PersistedIdentityCandidate):
        raise ValueError("Persisted identity presentation requires PersistedIdentityCandidate")
    equivalent_count = _nonnegative_count(candidate.equivalent_representation_count, "equivalent representation")
    hypothesis_count = _nonnegative_count(candidate.identity_hypothesis_count, "identity hypothesis")
    if candidate.translation.identity is None:
        if equivalent_count != 0:
            raise ValueError("An unresolved identity requires zero equivalent representations")
    else:
        if equivalent_count == 0:
            raise ValueError("A resolved identity requires a positive equivalent representation count")
        if hypothesis_count == 0:
            raise ValueError("A resolved identity requires a positive hypothesis count")
    return _result_cells(candidate.translation, equivalent_count, hypothesis_count)


def identity_translation_diagnostic_cells(translation: IdentityTranslation) -> dict[str, str | bool]:
    """Project one typed translation onto the Kestrel pre-result diagnostics.

    Args:
        translation: Complete-context translation captured before caller filtering.

    Returns:
        The exact four non-decision diagnostic cells.

    Raises:
        ValueError: If ``translation`` is not a validated identity translation.
    """
    if not isinstance(translation, IdentityTranslation):
        raise ValueError("Identity diagnostics require an IdentityTranslation")
    return {
        "Molecular_Identity": serialize_molecular_identity(translation.identity) if translation.identity else "",
        "Molecular_Identity_Translation_Status": translation.status,
        "Molecular_Identity_Translation_Failure": translation.failure or "",
        "Molecular_Identity_Context_Diverges": translation.context_diverges,
    }


def persisted_identity_result_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, str | int]]:
    """Decode and project selected Kestrel rows without pandas coupling.

    Args:
        rows: Positive Kestrel rows containing the complete internal identity codec.

    Returns:
        One public identity-cell mapping per input row, in input order.
    """
    from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells

    return [persisted_identity_result_cells(parse_selected_candidate_cells(row)) for row in rows]


def advntr_identity_result_rows(
    rows: Iterable[Mapping[str, object]],
    identity_component: IdentityTranslator,
) -> list[dict[str, str | int]]:
    """Translate and project the complete emitted adVNTR positive set.

    Args:
        rows: Every emitted positive adVNTR row.
        identity_component: Current-run complete-context translation component.

    Returns:
        One public identity-cell mapping per input row, in input order.
    """
    from vntyper.scripts.identity_candidates import capture_advntr_observations

    translations: list[IdentityTranslation] = []
    for row in rows:
        repeat_units = row["RU"]
        positions = row["POS"]
        repeat_units_absent = repeat_units == "."
        positions_absent = positions == "."
        if repeat_units_absent != positions_absent:
            raise ValueError("adVNTR RU and POS must both use the absent-context sentinel or neither use it")
        if repeat_units_absent:
            state_values = [row[column] for column in ("State", "Variant") if column in row]
            if len(state_values) != 1 or not isinstance(state_values[0], str) or not state_values[0]:
                raise ValueError("adVNTR row must contain exactly one non-empty State or Variant field")
            representation = AdvntrRepresentation(state_values[0], None, None)
            translations.append(identity_component.translate_advntr(representation))
            continue
        candidate = capture_advntr_observations((row,), identity_component).candidates[0]
        translations.append(candidate.observation.translation)
    translation_tuple = tuple(translations)
    return [
        identity_result_cells(selected, *(translation_tuple[:index] + translation_tuple[index + 1 :]))
        for index, selected in enumerate(translation_tuple)
    ]


def _result_cells(
    selected: IdentityTranslation,
    equivalent_count: int,
    hypothesis_count: int,
) -> dict[str, str | int]:
    if selected.identity is None:
        status = "unresolved"
        identity = ""
    else:
        status = "unique" if hypothesis_count == 1 else "legacy-selected-among-multiple"
        identity = serialize_molecular_identity(selected.identity)
    return {
        "Molecular_Identity": identity,
        "Molecular_Identity_Status": status,
        "Equivalent_Representation_Count": equivalent_count,
        "Identity_Hypothesis_Count": hypothesis_count,
    }


def _nonnegative_count(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"Persisted {name} count must be a non-negative integer")
    return value
