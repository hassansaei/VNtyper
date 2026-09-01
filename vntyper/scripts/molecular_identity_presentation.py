"""Pure public and diagnostic projections of typed molecular identity."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Any

from vntyper.scripts.molecular_identity import IdentityTranslation, serialize_molecular_identity

if TYPE_CHECKING:
    from vntyper.scripts.identity_candidate_persistence import PersistedIdentityCandidate

IDENTITY_COLUMNS: tuple[str, ...] = (
    "Molecular_Identity",
    "Molecular_Identity_Status",
    "Equivalent_Representation_Count",
    "Identity_Hypothesis_Count",
)

IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS: tuple[str, ...] = (
    "Molecular_Identity",
    "Molecular_Identity_Translation_Status",
    "Molecular_Identity_Translation_Failure",
    "Molecular_Identity_Context_Diverges",
)


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
    identity_component: Any,
) -> list[dict[str, str | int]]:
    """Translate and project the complete emitted adVNTR positive set.

    Args:
        rows: Every emitted positive adVNTR row.
        identity_component: Current-run complete-context translation component.

    Returns:
        One public identity-cell mapping per input row, in input order.
    """
    from vntyper.scripts.identity_candidates import capture_advntr_observations

    candidates = capture_advntr_observations(rows, identity_component).candidates
    translations = tuple(candidate.observation.translation for candidate in candidates)
    return [
        identity_result_cells(selected, *(translations[:index] + translations[index + 1 :]))
        for index, selected in enumerate(translations)
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
