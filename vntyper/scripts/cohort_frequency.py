"""Pure module for computing cohort call frequency tables.

This module groups Kestrel variant calls across a cohort, calculates call frequencies
relative to the full cohort roster (including samples without positive calls), marks
calls at or below a configured maximum frequency threshold, and produces a summary
DataFrame for presentation and export.

In accordance with VNtyper design, this feature is named 'cohort call frequency'
(never 'rare allele detection') because cohort mode aggregates one final call per
sample and does not evaluate a pre-filter allele spectrum.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd

from vntyper.scripts.report_formatting import is_empty_result_row

logger = logging.getLogger(__name__)

CALL_FREQUENCY_COLUMNS: list[str] = [
    "Grouping_Key",
    "Grouping_Key_Kind",
    "Molecular_Identity",
    "Motifs",
    "POS",
    "REF",
    "ALT",
    "Variant",
    "Sample_Count",
    "Frequency",
    "Below_Cutoff",
    "Samples",
    "Min_Depth_Score",
    "Max_Depth_Score",
]

_EXCLUDED_CALLER_METADATA = frozenset(
    {
        "Sample",
        "decision_profile_name",
        "decision_profile_revision",
        "decision_profile_sha256",
        "advntr_evidence_digest",
    }
)


def _is_placeholder(row: pd.Series | Mapping[str, Any]) -> bool:
    """Check if a row is the negative empty-result placeholder."""
    items = row.items() if hasattr(row, "items") else dict(row).items()
    cells = {k: "" if pd.isna(v) or v is None else str(v) for k, v in items if k not in _EXCLUDED_CALLER_METADATA}
    return is_empty_result_row(cells)


def _extract_grouping_key(row: pd.Series | Mapping[str, Any]) -> tuple[str, str]:
    """Determine the grouping key and key kind for one call row.

    Uses `Molecular_Identity` if `Molecular_Identity_Status` is `unique` or
    `legacy-selected-among-multiple` and the value is non-empty and not the
    legacy sentinel; otherwise falls back to the caller representation `(Motifs, POS, REF, ALT)`.

    Returns:
        tuple[str, str]: (Grouping_Key, Grouping_Key_Kind)
    """
    status = str(row.get("Molecular_Identity_Status") or "").strip()
    mol_id = str(row.get("Molecular_Identity") or "").strip()
    if status in ("unique", "legacy-selected-among-multiple") and mol_id and mol_id != "legacy identity not recorded":
        return mol_id, "molecular-identity"

    motifs = str(row.get("Motifs") or row.get("Motif") or "").strip()
    pos = str(row.get("POS") or "").strip()
    ref = str(row.get("REF") or "").strip()
    alt = str(row.get("ALT") or "").strip()
    return f"{motifs}:{pos}:{ref}:{alt}", "caller-representation"


def _first_non_empty(group: pd.DataFrame, *columns: str) -> str:
    """Extract first non-empty string value from given columns in group."""
    for col in columns:
        if col in group.columns:
            for val in group[col]:
                if pd.notna(val):
                    val_str = str(val).strip()
                    if val_str and val_str != "None":
                        return val_str
    return ""


def call_frequency_frame(
    kestrel_df: pd.DataFrame,
    cohort_size: int,
    max_frequency: float,
) -> pd.DataFrame:
    """Compute the cohort call frequency summary table.

    Args:
        kestrel_df: Aggregated Kestrel calls across the cohort.
        cohort_size: Total number of samples in the cohort roster (denominator).
        max_frequency: Frequency threshold at or below which calls are marked as Below_Cutoff.

    Returns:
        pd.DataFrame: Grouped call frequency table sorted ascending by Frequency, then Grouping_Key.
    """
    if kestrel_df.empty:
        return pd.DataFrame(columns=CALL_FREQUENCY_COLUMNS)

    # 1. Exclude placeholders using is_empty_result_row over row minus Sample
    placeholder_mask = kestrel_df.apply(_is_placeholder, axis=1)
    filtered_df = kestrel_df.loc[~placeholder_mask].copy()

    if filtered_df.empty:
        return pd.DataFrame(columns=CALL_FREQUENCY_COLUMNS)

    # 2. Coerce numeric columns after placeholder exclusion
    if "Depth_Score" in filtered_df.columns:
        filtered_df["_numeric_depth_score"] = pd.to_numeric(filtered_df["Depth_Score"], errors="coerce")
    else:
        filtered_df["_numeric_depth_score"] = np.nan

    # 3. Compute grouping keys
    key_tuples = filtered_df.apply(_extract_grouping_key, axis=1)
    filtered_df["_grouping_key"] = [kt[0] for kt in key_tuples]
    filtered_df["_grouping_kind"] = [kt[1] for kt in key_tuples]

    # 4. Group calls
    records: list[dict[str, Any]] = []
    for (key, kind), group in filtered_df.groupby(["_grouping_key", "_grouping_kind"]):
        samples = (
            sorted({str(s) for s in group["Sample"].dropna() if str(s).strip()}) if "Sample" in group.columns else []
        )
        sample_count = len(samples)
        raw_frequency = float(sample_count) / float(cohort_size) if cohort_size > 0 else 0.0
        below_cutoff = "yes" if raw_frequency <= max_frequency else "no"
        frequency = round(raw_frequency, 4)

        # Determine presentation representations
        mol_id = ""
        if kind == "molecular-identity":
            mol_id = str(key)
        else:
            if "Molecular_Identity" in group.columns:
                for val in group["Molecular_Identity"]:
                    val_str = str(val or "").strip()
                    if val_str and val_str not in ("None", "legacy identity not recorded"):
                        mol_id = val_str
                        break

        motifs = _first_non_empty(group, "Motifs", "Motif")
        pos = _first_non_empty(group, "POS")
        ref = _first_non_empty(group, "REF")
        alt = _first_non_empty(group, "ALT")
        variant = _first_non_empty(group, "Variant")

        numeric_scores = group["_numeric_depth_score"].dropna()
        min_depth = float(numeric_scores.min()) if not numeric_scores.empty else None
        max_depth = float(numeric_scores.max()) if not numeric_scores.empty else None

        records.append(
            {
                "Grouping_Key": str(key),
                "Grouping_Key_Kind": str(kind),
                "Molecular_Identity": mol_id,
                "Motifs": motifs,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "Variant": variant,
                "Sample_Count": sample_count,
                "Frequency": frequency,
                "Below_Cutoff": below_cutoff,
                "Samples": "; ".join(samples),
                "Min_Depth_Score": min_depth,
                "Max_Depth_Score": max_depth,
            }
        )

    result_df = pd.DataFrame(records, columns=CALL_FREQUENCY_COLUMNS)
    result_df = result_df.sort_values(by=["Frequency", "Grouping_Key"], ascending=[True, True]).reset_index(drop=True)
    return result_df
