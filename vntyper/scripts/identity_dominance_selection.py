"""Pure candidate grouping, metric aggregation, and identity-dominance selection decisions."""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Final

import pandas as pd

from vntyper.scripts.kestrel_decision_config import KestrelSelection
from vntyper.scripts.molecular_identity import MolecularIdentity, parse_molecular_identity

logger = logging.getLogger(__name__)

DEFAULT_CONFIDENCE_PRIORITY: Final[dict[str, int]] = {
    "High_Precision*": 3,
    "High_Precision": 2,
    "Low_Precision": 1,
    "Negative": 0,
}

UNFLAGGED_DEFAULT: Final[str] = "Not flagged"


@dataclass(frozen=True)
class IdentityCandidateGroup:
    """Group of passing Kestrel candidates sharing one canonical molecular identity."""

    identity_key: str
    identity: MolecularIdentity | None
    max_confidence_priority: int
    context_fidelity: int
    assembly_cardinality: int
    peak_alternate_depth: float
    peak_depth_score: float
    candidates: pd.DataFrame

    def __post_init__(self) -> None:
        """Validate group fields and invariants."""
        if not self.identity_key:
            raise ValueError("IdentityCandidateGroup requires a non-empty identity_key")
        if self.assembly_cardinality < 1:
            raise ValueError("IdentityCandidateGroup requires at least one candidate row")
        if self.context_fidelity not in (0, 1):
            raise ValueError("context_fidelity must be binary (0 or 1)")


def group_candidates_by_identity(
    df: pd.DataFrame,
    confidence_priority: Mapping[str, int] | None = None,
) -> list[IdentityCandidateGroup]:
    """Partition passing Kestrel candidate rows into distinct identity groups.

    Args:
        df: DataFrame of Kestrel candidate rows that have passed initial quality gates.
        confidence_priority: Mapping from confidence string to integer priority.

    Returns:
        List of IdentityCandidateGroup objects with computed group-level metrics.
    """
    if df.empty:
        return []

    priority_map = confidence_priority if confidence_priority is not None else DEFAULT_CONFIDENCE_PRIORITY

    grouped_rows: dict[str, list[pd.Series]] = defaultdict(list)
    for _, row in df.iterrows():
        raw_identity = row.get("Molecular_Identity")
        molid_str = str(raw_identity).strip() if pd.notna(raw_identity) else ""
        if molid_str:
            key = molid_str
        else:
            pos = row.get("POS", "")
            ref = row.get("REF", "")
            alt = row.get("ALT", "")
            key = f"unresolved:{pos}:{ref}:{alt}"
        grouped_rows[key].append(row)

    groups: list[IdentityCandidateGroup] = []
    for key, row_list in grouped_rows.items():
        group_df = pd.DataFrame(row_list)
        parsed_identity = None
        if not key.startswith("unresolved:"):
            try:
                parsed_identity = parse_molecular_identity(key)
            except ValueError:
                parsed_identity = None

        conf_series = group_df["Confidence"].map(lambda c: priority_map.get(str(c), 0))
        max_priority = int(conf_series.max()) if not conf_series.empty else 0

        # Context fidelity: 1 if at least one candidate has canonical X context
        if "Molecular_Identity_Context_Diverges" in group_df.columns:
            diverges_vals = group_df["Molecular_Identity_Context_Diverges"]
        else:
            diverges_vals = pd.Series([True] * len(group_df), index=group_df.index)
        has_canonical = any(v is False or str(v).strip().lower() == "false" for v in diverges_vals)
        context_fidelity = 1 if has_canonical else 0

        cardinality = len(group_df)
        if "Estimated_Depth_AlternateVariant" in group_df.columns:
            alt_depths = pd.to_numeric(group_df["Estimated_Depth_AlternateVariant"], errors="coerce").fillna(0.0)
        else:
            alt_depths = pd.Series(0.0, index=group_df.index)
        peak_alt_depth = float(alt_depths.max()) if not alt_depths.empty else 0.0

        if "Depth_Score" in group_df.columns:
            depth_scores = pd.to_numeric(group_df["Depth_Score"], errors="coerce").fillna(0.0)
        else:
            depth_scores = pd.Series(0.0, index=group_df.index)
        peak_depth_score = float(depth_scores.max()) if not depth_scores.empty else 0.0

        groups.append(
            IdentityCandidateGroup(
                identity_key=key,
                identity=parsed_identity,
                max_confidence_priority=max_priority,
                context_fidelity=context_fidelity,
                assembly_cardinality=cardinality,
                peak_alternate_depth=peak_alt_depth,
                peak_depth_score=peak_depth_score,
                candidates=group_df,
            )
        )

    return groups


def rank_candidate_groups(groups: Sequence[IdentityCandidateGroup]) -> list[IdentityCandidateGroup]:
    """Rank candidate groups lexicographically by clinical priority and support.

    Sort hierarchy:
        1. max_confidence_priority DESC (High_Precision* > High_Precision > Low_Precision)
        2. context_fidelity DESC (canonical X repeat context preferred over divergent motifs)
        3. peak_alternate_depth DESC (highest alternate k-mer depth)
        4. assembly_cardinality DESC (most supporting motif pair assemblies)
        5. peak_depth_score DESC (highest allele fraction)

    Args:
        groups: Sequence of candidate groups to rank.

    Returns:
        Sorted list of candidate groups, highest ranked first.
    """
    return sorted(
        groups,
        key=lambda g: (
            g.max_confidence_priority,
            g.context_fidelity,
            g.peak_alternate_depth,
            g.assembly_cardinality,
            g.peak_depth_score,
        ),
        reverse=True,
    )


def select_group_representative_row(
    group: IdentityCandidateGroup,
    confidence_priority: Mapping[str, int] | None = None,
    unflagged_value: str = UNFLAGGED_DEFAULT,
) -> pd.DataFrame:
    """Select the single best representative row from an identity candidate group.

    Tie-breaking hierarchy within the winning identity group:
        1. Confidence DESC
        2. is_unflagged DESC
        3. haplo_count DESC
        4. Depth_Score DESC
        5. POS ASC

    Args:
        group: The winning IdentityCandidateGroup.
        confidence_priority: Mapping from confidence string to integer priority.
        unflagged_value: Value in Flag column considered unflagged.

    Returns:
        DataFrame containing exactly one row representing the winning group.
    """
    df = group.candidates.copy()
    if df.empty:
        return df
    if len(df) == 1:
        return df

    priority_map = confidence_priority if confidence_priority is not None else DEFAULT_CONFIDENCE_PRIORITY

    df["_priority"] = df["Confidence"].map(lambda c: priority_map.get(str(c), 0))
    if "Flag" in df.columns:
        df["_is_unflagged"] = (df["Flag"] == unflagged_value).astype(int)
    else:
        df["_is_unflagged"] = 1

    if "haplo_count" in df.columns:
        df["_haplo_count"] = pd.to_numeric(df["haplo_count"], errors="coerce").fillna(0)
    else:
        df["_haplo_count"] = 0

    if "Depth_Score" in df.columns:
        df["_depth_score"] = pd.to_numeric(df["Depth_Score"], errors="coerce").fillna(0.0)
    else:
        df["_depth_score"] = 0.0

    if "POS" in df.columns:
        df["_pos"] = pd.to_numeric(df["POS"], errors="coerce").fillna(0)
    else:
        df["_pos"] = 0

    sorted_df = df.sort_values(
        by=["_priority", "_is_unflagged", "_haplo_count", "_depth_score", "_pos"],
        ascending=[False, False, False, False, True],
    )

    added_cols = ["_priority", "_is_unflagged", "_haplo_count", "_depth_score", "_pos"]
    return sorted_df.head(1).drop(columns=added_cols)


def select_identity_dominance_variant(
    df: pd.DataFrame,
    selection: KestrelSelection | Mapping[str, object] | None = None,
) -> pd.DataFrame:
    """Select the single dominant variant by grouping passing rows on molecular identity.

    Args:
        df: Pre-filtered DataFrame of Kestrel candidate variants.
        selection: Optional selection configuration or typed KestrelSelection.

    Returns:
        DataFrame with exactly 1 row (the dominant variant), or empty if input empty.
    """
    if df.empty:
        return df
    if len(df) == 1:
        return df

    confidence_priority = None
    unflagged_value = UNFLAGGED_DEFAULT
    if isinstance(selection, KestrelSelection):
        confidence_priority = selection.confidence_priority
        unflagged_value = selection.unflagged_value
    elif isinstance(selection, Mapping):
        conf_prio = selection.get("confidence_priority")
        if isinstance(conf_prio, Mapping):
            confidence_priority = {str(k): int(v) for k, v in conf_prio.items()}
        unflagged_val = selection.get("unflagged_value")
        if isinstance(unflagged_val, str):
            unflagged_value = unflagged_val

    groups = group_candidates_by_identity(df, confidence_priority=confidence_priority)
    if not groups:
        return df.head(1)

    ranked_groups = rank_candidate_groups(groups)
    winning_group = ranked_groups[0]

    result = select_group_representative_row(
        winning_group,
        confidence_priority=confidence_priority,
        unflagged_value=unflagged_value,
    )

    logger.info(
        "Selected dominant identity '%s' with %d supporting assemblies, peak AltDepth=%.1f",
        winning_group.identity_key,
        winning_group.assembly_cardinality,
        winning_group.peak_alternate_depth,
    )
    return result
