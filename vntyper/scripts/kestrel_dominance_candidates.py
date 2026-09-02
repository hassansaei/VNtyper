"""Pure frame projections for generated whole-locus dominance candidates."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from vntyper.scripts.identity_candidate_persistence import IDENTITY_SELECTION_COLUMNS
from vntyper.scripts.identity_candidates import OBSERVATION_ORDINAL_COLUMN


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
