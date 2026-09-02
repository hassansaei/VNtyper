"""Pure frame projections for generated whole-locus dominance candidates."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import pandas as pd

from vntyper.scripts.identity_candidate_persistence import IDENTITY_SELECTION_COLUMNS
from vntyper.scripts.identity_candidates import OBSERVATION_ORDINAL_COLUMN
from vntyper.scripts.molecular_identity_presentation import IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS


def passing_candidate_frame(frame: pd.DataFrame, filter_columns: Sequence[str]) -> pd.DataFrame:
    """Return every row surviving the fixed legacy gates in original order.

    Args:
        frame: Complete pre-selection Kestrel frame.
        filter_columns: Frozen legacy gate column names.

    Returns:
        A copy containing every positive candidate hypothesis.

    Raises:
        ValueError: If gate or authoritative projection columns are absent or malformed.
    """
    required = {*filter_columns, OBSERVATION_ORDINAL_COLUMN, *IDENTITY_SELECTION_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"complete Kestrel dominance candidates lack required columns: {missing}")
    mask = pd.Series(True, index=frame.index)
    for column in filter_columns:
        values = frame[column]
        valid = values.map(lambda value: value is True or value == "True")
        malformed = ~(valid | values.map(lambda value: value is False or value == "False"))
        if malformed.any():
            raise ValueError(f"complete Kestrel dominance gate {column} must contain Boolean values")
        mask &= valid
    passing = frame.loc[mask].copy()
    ordinals = tuple(int(value) for value in passing[OBSERVATION_ORDINAL_COLUMN])
    if ordinals != tuple(sorted(set(ordinals))):
        raise ValueError("complete Kestrel dominance candidate ordinals must be unique and increasing")
    return passing


def legacy_result_candidates(frame: pd.DataFrame) -> pd.DataFrame:
    """Drop the pre-result-only translation diagnostics before final annotation.

    The legacy selection path removes these four columns from the selected row before
    the public identity quartet is appended. The enabled-dominance path annotates rows
    projected straight out of the retained pre-result, so it must remove them too, or
    its ``kestrel_result.tsv`` gains three columns and moves ``Molecular_Identity`` out
    of the quartet position the legacy result publishes.

    Args:
        frame: Passing candidate rows projected from the complete pre-result.

    Returns:
        A copy without the four diagnostic columns, in the legacy column order.

    Raises:
        ValueError: If any diagnostic column is absent.
    """
    missing = sorted(set(IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"complete Kestrel dominance candidates lack translation diagnostics: {missing}")
    return frame.drop(columns=list(IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS))


def merge_candidate_annotations(pre_result: pd.DataFrame, annotated: pd.DataFrame) -> pd.DataFrame:
    """Merge annotated passing rows back into the complete pre-result frame.

    Args:
        pre_result: Complete pre-selection frame, including failed evidence rows.
        annotated: Passing subset with fixed nomenclature projections.

    Returns:
        Complete frame retaining every original row and adding annotations by index.

    Raises:
        ValueError: If annotated rows are not a subset of the pre-result index.
    """
    if not set(annotated.index) <= set(pre_result.index):
        raise ValueError("annotated dominance candidates must belong to the complete pre-result")
    merged = pre_result.copy()
    for column in annotated.columns:
        if column not in merged.columns:
            merged[column] = ""
        merged.loc[annotated.index, column] = annotated[column]
    return merged


def publish_dominance_selection(
    pre_result: pd.DataFrame,
    filter_columns: Sequence[str],
    selected_ordinal: int,
    annotate: Callable[[pd.DataFrame], pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the fixed enabled-dominance publication sequence on a complete pre-result.

    This is the one sequence production uses when a profile enables dominance: every
    passing hypothesis is projected out of the retained pre-result, stripped of the
    pre-result-only translation diagnostics exactly as the legacy path strips them,
    annotated, merged back so the pre-result keeps every row and its diagnostics, and
    the legacy-selected ordinal is published unchanged.

    Args:
        pre_result: Complete pre-selection Kestrel frame as read from disk.
        filter_columns: Frozen legacy gate column names.
        selected_ordinal: Exact ordinal the legacy ranking selected.
        annotate: The fixed nomenclature annotation for the passing candidates.

    Returns:
        The merged pre-result to persist and the one annotated legacy-selected row.

    Raises:
        ValueError: If any projection step rejects the frame.
    """
    candidates = legacy_result_candidates(passing_candidate_frame(pre_result, filter_columns))
    annotated = annotate(candidates)
    merged = merge_candidate_annotations(pre_result, annotated)
    return merged, selected_candidate_frame(annotated, selected_ordinal)


def selected_candidate_frame(annotated: pd.DataFrame, selected_ordinal: int) -> pd.DataFrame:
    """Project the unchanged legacy-selected row from complete annotations.

    Args:
        annotated: Complete passing candidate annotations.
        selected_ordinal: Exact ordinal selected by legacy ranking.

    Returns:
        The one fixed annotated legacy selection.

    Raises:
        ValueError: If the selected ordinal is missing or duplicated.
    """
    if isinstance(selected_ordinal, bool) or not isinstance(selected_ordinal, int) or selected_ordinal < 0:
        raise ValueError("legacy-selected dominance ordinal must be a non-negative integer")
    selected = annotated.loc[annotated[OBSERVATION_ORDINAL_COLUMN].map(int) == selected_ordinal].copy()
    if len(selected) != 1:
        raise ValueError("legacy-selected dominance ordinal must identify exactly one complete candidate")
    return selected
