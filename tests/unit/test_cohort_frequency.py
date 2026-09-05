"""Unit tests for pure cohort call frequency aggregation."""

from __future__ import annotations

import pandas as pd
import pytest

from vntyper.scripts.cohort_frequency import (
    CALL_FREQUENCY_COLUMNS,
    call_frequency_frame,
)

pytestmark = pytest.mark.unit


def test_call_frequency_frame_columns_and_empty_input() -> None:
    """An empty input DataFrame yields an empty frame with exact declared columns."""
    empty_df = pd.DataFrame()
    result = call_frequency_frame(empty_df, cohort_size=5, max_frequency=0.05)

    assert list(result.columns) == CALL_FREQUENCY_COLUMNS
    assert result.empty
    assert "Flag" not in result.columns


def test_placeholder_rows_are_completely_excluded() -> None:
    """Empty-result placeholder rows are excluded and do not form a 'None' group."""
    rows = [
        {
            "Sample": "sample_1",
            "Motif": "None",
            "Variant": "None",
            "POS": "None",
            "REF": "None",
            "ALT": "None",
            "Motif_sequence": "None",
            "Estimated_Depth_AlternateVariant": "None",
            "Estimated_Depth_Variant_ActiveRegion": "None",
            "Depth_Score": "None",
            "Confidence": "Negative",
        },
        {
            "Sample": "sample_2",
            "Molecular_Identity": "chr1|155160963|155160963|C|CC",
            "Molecular_Identity_Status": "unique",
            "Motifs": "5-5",
            "POS": "155160963",
            "REF": "C",
            "ALT": "CC",
            "Variant": "insC",
            "Depth_Score": "25.5",
            "Confidence": "High_Precision",
        },
    ]
    df = pd.DataFrame(rows)
    result = call_frequency_frame(df, cohort_size=2, max_frequency=0.05)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["Grouping_Key"] == "chr1|155160963|155160963|C|CC"
    assert row["Grouping_Key_Kind"] == "molecular-identity"
    assert row["Sample_Count"] == 1
    assert row["Frequency"] == 0.5
    assert row["Below_Cutoff"] == "no"
    assert row["Samples"] == "sample_2"
    assert row["Min_Depth_Score"] == 25.5
    assert row["Max_Depth_Score"] == 25.5


def test_molecular_identity_grouping_when_unique_or_legacy_selected() -> None:
    """Unique and legacy-selected-among-multiple molecular identities are grouped by identity."""
    rows = [
        {
            "Sample": "sample_1",
            "Molecular_Identity": "chr1|100|100|A|G",
            "Molecular_Identity_Status": "unique",
            "Motifs": "1-1",
            "POS": "100",
            "REF": "A",
            "ALT": "G",
            "Variant": "snv1",
            "Depth_Score": "10.0",
        },
        {
            "Sample": "sample_2",
            "Molecular_Identity": "chr1|100|100|A|G",
            "Molecular_Identity_Status": "legacy-selected-among-multiple",
            "Motifs": "1-1",
            "POS": "100",
            "REF": "A",
            "ALT": "G",
            "Variant": "snv1",
            "Depth_Score": "12.0",
        },
    ]
    df = pd.DataFrame(rows)
    result = call_frequency_frame(df, cohort_size=10, max_frequency=0.25)

    assert len(result) == 1
    row = result.iloc[0]
    assert row["Grouping_Key"] == "chr1|100|100|A|G"
    assert row["Grouping_Key_Kind"] == "molecular-identity"
    assert row["Sample_Count"] == 2
    assert row["Frequency"] == 0.2
    assert row["Below_Cutoff"] == "yes"
    assert row["Samples"] == "sample_1; sample_2"
    assert row["Min_Depth_Score"] == 10.0
    assert row["Max_Depth_Score"] == 12.0


def test_legacy_and_unresolved_rows_fall_back_to_caller_representation() -> None:
    """Rows with unresolved or 'legacy identity not recorded' fall back to caller representation without collapsing."""
    rows = [
        {
            "Sample": "sample_1",
            "Molecular_Identity": "legacy identity not recorded",
            "Molecular_Identity_Status": "legacy identity not recorded",
            "Motifs": "1-1",
            "POS": "100",
            "REF": "A",
            "ALT": "T",
            "Variant": "snv1",
            "Depth_Score": "8.0",
        },
        {
            "Sample": "sample_2",
            "Molecular_Identity": "legacy identity not recorded",
            "Molecular_Identity_Status": "legacy identity not recorded",
            "Motifs": "2-2",
            "POS": "200",
            "REF": "C",
            "ALT": "G",
            "Variant": "snv2",
            "Depth_Score": "9.0",
        },
        {
            "Sample": "sample_3",
            "Molecular_Identity": "",
            "Molecular_Identity_Status": "unresolved",
            "Motifs": "1-1",
            "POS": "100",
            "REF": "A",
            "ALT": "T",
            "Variant": "snv1",
            "Depth_Score": "11.0",
        },
    ]
    df = pd.DataFrame(rows)
    result = call_frequency_frame(df, cohort_size=5, max_frequency=0.5)

    assert len(result) == 2
    keys = list(result["Grouping_Key"])
    assert "1-1:100:A:T" in keys
    assert "2-2:200:C:G" in keys

    group_1 = result[result["Grouping_Key"] == "1-1:100:A:T"].iloc[0]
    assert group_1["Grouping_Key_Kind"] == "caller-representation"
    assert group_1["Sample_Count"] == 2
    assert group_1["Frequency"] == 0.4
    assert group_1["Samples"] == "sample_1; sample_3"
    assert group_1["Min_Depth_Score"] == 8.0
    assert group_1["Max_Depth_Score"] == 11.0

    group_2 = result[result["Grouping_Key"] == "2-2:200:C:G"].iloc[0]
    assert group_2["Grouping_Key_Kind"] == "caller-representation"
    assert group_2["Sample_Count"] == 1
    assert group_2["Frequency"] == 0.2
    assert group_2["Samples"] == "sample_2"


def test_numeric_coercion_and_ordering() -> None:
    """Depth scores are compared numerically, not lexicographically; results sorted by Frequency then key."""
    rows = [
        {
            "Sample": "sample_1",
            "Molecular_Identity": "chr1|10|10|A|T",
            "Molecular_Identity_Status": "unique",
            "Motifs": "1-1",
            "POS": "10",
            "REF": "A",
            "ALT": "T",
            "Variant": "varB",
            "Depth_Score": "9.5",
        },
        {
            "Sample": "sample_2",
            "Molecular_Identity": "chr1|10|10|A|T",
            "Molecular_Identity_Status": "unique",
            "Motifs": "1-1",
            "POS": "10",
            "REF": "A",
            "ALT": "T",
            "Variant": "varB",
            "Depth_Score": "100.0",
        },
        {
            "Sample": "sample_3",
            "Molecular_Identity": "chr1|5|5|C|G",
            "Molecular_Identity_Status": "unique",
            "Motifs": "1-1",
            "POS": "5",
            "REF": "C",
            "ALT": "G",
            "Variant": "varA",
            "Depth_Score": "15.0",
        },
    ]
    df = pd.DataFrame(rows)
    result = call_frequency_frame(df, cohort_size=4, max_frequency=0.3)

    assert len(result) == 2
    assert result.iloc[0]["Grouping_Key"] == "chr1|5|5|C|G"
    assert result.iloc[0]["Frequency"] == 0.25
    assert result.iloc[0]["Below_Cutoff"] == "yes"

    assert result.iloc[1]["Grouping_Key"] == "chr1|10|10|A|T"
    assert result.iloc[1]["Frequency"] == 0.5
    assert result.iloc[1]["Below_Cutoff"] == "no"
    assert result.iloc[1]["Min_Depth_Score"] == 9.5
    assert result.iloc[1]["Max_Depth_Score"] == 100.0


def test_denominator_comes_from_cohort_size_parameter() -> None:
    """Cohort size parameter controls denominator even when few samples have rows."""
    rows = [
        {
            "Sample": "s1",
            "Molecular_Identity": "id1",
            "Molecular_Identity_Status": "unique",
            "Motifs": "1-1",
            "POS": "10",
            "REF": "A",
            "ALT": "C",
            "Variant": "snv",
            "Depth_Score": "20.0",
        }
    ]
    df = pd.DataFrame(rows)
    result = call_frequency_frame(df, cohort_size=100, max_frequency=0.05)
    assert result.iloc[0]["Frequency"] == 0.01
    assert result.iloc[0]["Below_Cutoff"] == "yes"


def test_zero_or_negative_cohort_size_handled_safely() -> None:
    """Zero cohort size sets frequency to 0.0 and marks Below_Cutoff appropriately without ZeroDivisionError."""
    rows = [
        {
            "Sample": "s1",
            "Molecular_Identity": "id1",
            "Molecular_Identity_Status": "unique",
            "Motifs": "1-1",
            "POS": "10",
            "REF": "A",
            "ALT": "C",
            "Variant": "snv",
            "Depth_Score": "20.0",
        }
    ]
    df = pd.DataFrame(rows)
    result = call_frequency_frame(df, cohort_size=0, max_frequency=0.05)
    assert result.iloc[0]["Frequency"] == 0.0
    assert result.iloc[0]["Below_Cutoff"] == "yes"


def test_all_placeholder_rows_returns_empty_frame() -> None:
    """DataFrame with only placeholder rows returns empty frame with declared columns."""
    rows = [
        {
            "Sample": "s1",
            "Motif": "None",
            "Variant": "None",
            "POS": "None",
            "REF": "None",
            "ALT": "None",
            "Motif_sequence": "None",
            "Estimated_Depth_AlternateVariant": "None",
            "Estimated_Depth_Variant_ActiveRegion": "None",
            "Depth_Score": "None",
            "Confidence": "Negative",
        }
    ]
    df = pd.DataFrame(rows)
    result = call_frequency_frame(df, cohort_size=5, max_frequency=0.05)
    assert result.empty
    assert list(result.columns) == CALL_FREQUENCY_COLUMNS
