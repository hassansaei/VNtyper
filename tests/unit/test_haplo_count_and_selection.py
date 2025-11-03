#!/usr/bin/env python3
"""
tests/unit/test_haplo_count_and_selection.py

Unit tests for the new functions:
- add_haplo_count: Calculate supporting evidence count for variants
- select_single_best_variant: Select single best variant with deterministic tie-breaking

These functions are part of the Issue #136 fix implementation.

References:
    - Issue #136: https://github.com/hassansaei/VNtyper/issues/136
    - PR #140 code review
    - Streamlined implementation plan (Day 5)
"""

import pandas as pd
import pytest

from vntyper.scripts.kestrel_genotyping import add_haplo_count, select_single_best_variant


class TestAddHaploCount:
    """Test the add_haplo_count function."""

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty with haplo_count=0."""
        df = pd.DataFrame()
        result = add_haplo_count(df)
        assert result.empty
        assert "haplo_count" in result.columns

    def test_basic_haplo_count(self):
        """Basic haplo_count calculation."""
        df = pd.DataFrame(
            {
                "POS": [67, 67, 67, 54],
                "REF": ["G", "G", "G", "C"],
                "ALT": ["GG", "GG", "GG", "GC"],
            }
        )

        result = add_haplo_count(df)

        assert result.iloc[0]["haplo_count"] == 3  # GG at 67: 3 times
        assert result.iloc[1]["haplo_count"] == 3
        assert result.iloc[2]["haplo_count"] == 3
        assert result.iloc[3]["haplo_count"] == 1  # GC at 54: 1 time

    def test_haplo_count_position_specific(self):
        """haplo_count must be position-specific."""
        df = pd.DataFrame(
            {
                "POS": [67, 67, 54],
                "REF": ["G", "C", "G"],
                "ALT": ["GG", "GG", "GG"],  # Same ALT, different POS/REF
            }
        )
        result = add_haplo_count(df)
        # Each is unique by (POS, REF, ALT) combination
        assert all(result["haplo_count"] == 1)

    def test_haplo_count_different_positions_same_alt(self):
        """Different positions with same ALT should be counted separately."""
        df = pd.DataFrame(
            {
                "POS": [54, 54, 67, 67],
                "REF": ["C", "C", "G", "G"],
                "ALT": ["GC", "GC", "GC", "GC"],
            }
        )
        result = add_haplo_count(df)
        # POS=54: 2 GC variants
        assert result.iloc[0]["haplo_count"] == 2
        assert result.iloc[1]["haplo_count"] == 2
        # POS=67: 2 GC variants
        assert result.iloc[2]["haplo_count"] == 2
        assert result.iloc[3]["haplo_count"] == 2

    def test_haplo_count_missing_columns(self, caplog):
        """Missing POS/REF/ALT columns should log warning and set haplo_count=0."""
        df = pd.DataFrame({"SomeColumn": [1, 2, 3]})
        result = add_haplo_count(df)
        assert "haplo_count" in result.columns
        assert all(result["haplo_count"] == 0)
        assert "Missing POS/REF/ALT columns" in caplog.text

    def test_haplo_count_preserves_other_columns(self):
        """add_haplo_count should preserve all other columns."""
        df = pd.DataFrame(
            {
                "POS": [67, 67],
                "REF": ["G", "G"],
                "ALT": ["GG", "GG"],
                "Depth_Score": [0.01, 0.009],
                "Confidence": ["High_Precision", "High_Precision"],
            }
        )
        result = add_haplo_count(df)
        assert "Depth_Score" in result.columns
        assert "Confidence" in result.columns
        assert len(result) == 2

    def test_haplo_count_realistic_example(self):
        """Realistic example from insG_pos54 sample."""
        # Simulated data: 389 GG variants at POS 67, 176 GC variants at POS 54
        df = pd.DataFrame(
            {
                "POS": [67] * 389 + [54] * 176,
                "REF": ["G"] * 389 + ["C"] * 176,
                "ALT": ["GG"] * 389 + ["GC"] * 176,
            }
        )
        result = add_haplo_count(df)
        # Check first GG variant
        assert result.iloc[0]["haplo_count"] == 389
        # Check first GC variant
        assert result.iloc[389]["haplo_count"] == 176


class TestSelectSingleBestVariant:
    """Test the select_single_best_variant function."""

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty."""
        df = pd.DataFrame()
        result = select_single_best_variant(df)
        assert result.empty

    def test_single_variant_returns_as_is(self):
        """Single variant should be returned unchanged."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision"],
                "haplo_count": [100],
                "Depth_Score": [0.01],
                "POS": [67],
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        assert result.iloc[0]["Confidence"] == "High_Precision"

    def test_confidence_priority(self):
        """Higher confidence level wins regardless of other metrics."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision*", "High_Precision", "Low_Precision"],
                "haplo_count": [50, 389, 15],  # Middle has highest haplo_count
                "Depth_Score": [0.008, 0.010, 0.006],
                "POS": [60, 67, 54],
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        assert result.iloc[0]["Confidence"] == "High_Precision*"  # Wins on priority

    def test_haplo_count_tie_breaker(self):
        """When Confidence tied, highest haplo_count wins."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision", "High_Precision", "High_Precision"],
                "haplo_count": [100, 389, 50],  # Middle has highest
                "Depth_Score": [0.01, 0.009, 0.011],
                "POS": [60, 67, 54],
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        assert result.iloc[0]["haplo_count"] == 389
        assert result.iloc[0]["POS"] == 67

    def test_depth_score_tie_breaker(self):
        """When Confidence and haplo_count tied, highest Depth_Score wins."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision"] * 3,
                "haplo_count": [100, 100, 100],  # All tied
                "Depth_Score": [0.009, 0.011, 0.010],  # Middle highest
                "POS": [60, 67, 54],
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        assert result.iloc[0]["Depth_Score"] == 0.011
        assert result.iloc[0]["POS"] == 67

    def test_pos_tie_breaker(self):
        """When all metrics tied, lowest POS wins (reproducibility)."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision"] * 3,
                "haplo_count": [100, 100, 100],  # All tied
                "Depth_Score": [0.010, 0.010, 0.010],  # All tied
                "POS": [67, 54, 60],  # Lowest is 54
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        assert result.iloc[0]["POS"] == 54  # Lowest POS wins

    def test_deterministic_tie_breaking(self):
        """Complete deterministic tie-breaking scenario."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision"] * 3,
                "haplo_count": [100, 100, 100],  # Tie
                "Depth_Score": [0.010, 0.009, 0.010],  # Two tied at 0.010
                "POS": [67, 54, 60],  # Among 0.010 group, lowest is 60
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        # First row wins (POS=67, Depth_Score=0.010, haplo_count=100)
        # because multi-key sort gives first match
        assert result.iloc[0]["POS"] == 60

    def test_missing_haplo_count_column(self):
        """Missing haplo_count should default to 0."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision", "High_Precision"],
                "Depth_Score": [0.01, 0.009],
                "POS": [67, 54],
                # No haplo_count column
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        # Should select based on Depth_Score (0.01 > 0.009)
        assert result.iloc[0]["POS"] == 67

    def test_negative_confidence_lowest_priority(self):
        """Negative confidence should have lowest priority."""
        df = pd.DataFrame(
            {
                "Confidence": ["Low_Precision", "Negative", "High_Precision"],
                "haplo_count": [10, 1000, 50],  # Negative has highest haplo_count
                "Depth_Score": [0.005, 0.020, 0.008],
                "POS": [54, 60, 67],
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        assert result.iloc[0]["Confidence"] == "High_Precision"  # Wins on priority

    def test_unknown_confidence_defaults_to_zero(self):
        """Unknown confidence levels should get priority 0."""
        df = pd.DataFrame(
            {
                "Confidence": ["UnknownLevel", "High_Precision"],
                "haplo_count": [1000, 10],
                "Depth_Score": [0.02, 0.008],
                "POS": [54, 67],
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        # High_Precision should win (priority=2 vs 0)
        assert result.iloc[0]["Confidence"] == "High_Precision"

    def test_no_priority_column_in_result(self):
        """Result should not contain temporary _priority column."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision", "Low_Precision"],
                "haplo_count": [100, 50],
                "Depth_Score": [0.01, 0.005],
                "POS": [67, 54],
            }
        )
        result = select_single_best_variant(df)
        assert "_priority" not in result.columns

    def test_preserves_all_original_columns(self):
        """Result should preserve all original columns."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision", "Low_Precision"],
                "haplo_count": [100, 50],
                "Depth_Score": [0.01, 0.005],
                "POS": [67, 54],
                "REF": ["G", "C"],
                "ALT": ["GG", "GC"],
                "Motif": ["X", "Y"],
                "ExtraColumn": ["A", "B"],
            }
        )
        result = select_single_best_variant(df)
        assert "REF" in result.columns
        assert "ALT" in result.columns
        assert "Motif" in result.columns
        assert "ExtraColumn" in result.columns

    def test_numeric_coercion(self):
        """Non-numeric values should be coerced to numeric."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision", "High_Precision"],
                "haplo_count": ["100", "50"],  # String numbers
                "Depth_Score": ["0.01", "0.005"],  # String numbers
                "POS": ["67", "54"],  # String numbers
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        # Should convert and select correctly
        assert int(result.iloc[0]["haplo_count"]) == 100

    def test_realistic_scenario_gg_vs_gc(self):
        """Realistic scenario: GG at POS 67 vs GC at POS 54."""
        df = pd.DataFrame(
            {
                "Confidence": ["High_Precision", "High_Precision"],
                "haplo_count": [389, 176],  # GG has more support
                "Depth_Score": [0.010233, 0.008842],  # GG has higher depth
                "POS": [67, 54],
                "REF": ["G", "C"],
                "ALT": ["GG", "GC"],
            }
        )
        result = select_single_best_variant(df)
        assert len(result) == 1
        # GG should win (same Confidence, higher haplo_count)
        assert result.iloc[0]["ALT"] == "GG"
        assert result.iloc[0]["POS"] == 67
        assert result.iloc[0]["haplo_count"] == 389


class TestIntegration:
    """Integration tests combining both functions."""

    def test_add_haplo_count_then_select(self):
        """Typical workflow: add haplo_count, then select single best."""
        # Simulate variants from Kestrel output
        df = pd.DataFrame(
            {
                "POS": [67, 67, 67, 54, 54],
                "REF": ["G", "G", "G", "C", "C"],
                "ALT": ["GG", "GG", "GG", "GC", "GC"],
                "Depth_Score": [0.010, 0.009, 0.011, 0.008, 0.007],
                "Confidence": ["High_Precision"] * 5,
            }
        )

        # Step 1: Add haplo_count
        df = add_haplo_count(df)
        assert df.iloc[0]["haplo_count"] == 3  # GG: 3 times
        assert df.iloc[3]["haplo_count"] == 2  # GC: 2 times

        # Step 2: Select single best
        result = select_single_best_variant(df)
        assert len(result) == 1
        # GG should win (higher haplo_count)
        assert result.iloc[0]["ALT"] == "GG"
        assert result.iloc[0]["haplo_count"] == 3

    def test_condition_6_integration(self):
        """Test that Condition 6 fix integrates correctly."""
        # Variants with intermediate depth scores
        df = pd.DataFrame(
            {
                "POS": [67, 54],
                "REF": ["G", "C"],
                "ALT": ["GG", "GC"],
                "Depth_Score": [0.00490, 0.00480],  # Both in intermediate range
                "Confidence": ["Low_Precision", "Low_Precision"],  # Both assigned via cond6
                "haplo_count": [389, 176],
            }
        )

        result = select_single_best_variant(df)
        assert len(result) == 1
        # GG should win (same Confidence, higher haplo_count)
        assert result.iloc[0]["ALT"] == "GG"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
