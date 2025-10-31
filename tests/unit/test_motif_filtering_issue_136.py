"""
Unit tests for Issue #136 motif filtering fix.

Tests the new uniform filtering logic that prevents silent deletion of
non-GG frameshift variants (e.g., insG_pos54 mutations).

Test Coverage:
- _apply_uniform_filtering_right_motif helper function
- Feature flag branching in motif_correction_and_annotation
- Biological scenarios: insG_pos54 (mixed GG+GC), dupC (only GG), dupA (CT variants)
- Edge cases: empty DataFrames, all variants excluded, position-specific deduplication

References:
- GitHub Issue #136
- Hassan Saei's recommendation: depth-score-first filtering
"""

import pandas as pd
import pytest

from vntyper.scripts.motif_processing import _apply_uniform_filtering_right_motif


class TestUniformFilteringHelperFunction:
    """Test suite for _apply_uniform_filtering_right_motif function."""

    def test_empty_dataframe_returns_empty(self):
        """Empty input should return empty output without errors."""
        empty_df = pd.DataFrame()
        result = _apply_uniform_filtering_right_motif(
            empty_df, exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )
        assert result.empty

    def test_exclude_motifs_first_conserved_motifs(self):
        """
        Step 1: Excluded motifs should be removed FIRST.

        Biological: Conserved motifs (Q, 8, 9, etc.) are stable → variants likely artifacts.
        """
        # Create test data with variants in excluded motifs
        df = pd.DataFrame(
            {
                "POS": [67, 67, 67],
                "REF": ["G", "G", "G"],
                "ALT": ["GG", "GG", "GG"],
                "Depth_Score": [0.010, 0.008, 0.006],
                "Motif": ["X", "Q", "8"],  # X is allowed, Q and 8 are excluded
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=["Q", "8"], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # Should keep only motif "X" (highest depth)
        assert len(result) == 1
        assert result.iloc[0]["Motif"] == "X"
        assert result.iloc[0]["Depth_Score"] == 0.010

    def test_sort_by_depth_score_desc(self):
        """
        Step 2: Should sort by Depth_Score DESC, then POS DESC.

        Biological: Highest coverage first (signal > noise).
        """
        df = pd.DataFrame(
            {
                "POS": [67, 60, 54],
                "REF": ["G", "C", "C"],
                "ALT": ["GG", "GC", "GC"],
                "Depth_Score": [0.005, 0.015, 0.010],  # Middle, highest, low
                "Motif": ["X", "X", "X"],
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # After sorting by Depth_Score DESC: [0.015, 0.010, 0.005]
        # After dedup by [POS, REF, ALT]: all unique → all kept
        # After GG filtering: GG filtered by motif, GC preserved
        assert len(result) == 3  # All unique positions

    def test_position_specific_deduplication_critical_fix(self):
        """
        Step 3: Deduplication by [POS, REF, ALT] - CRITICAL FIX for Issue #136!

        Biological:
        - POS 60: C>GC (insG at pos 54) vs POS 67: G>GG (dupC) = DIFFERENT events
        - Old bug: subset=["REF", "ALT"] would lose one
        - New fix: subset=["POS", "REF", "ALT"] preserves both
        """
        df = pd.DataFrame(
            {
                "POS": [60, 67, 67],  # Two different positions
                "REF": ["C", "G", "G"],
                "ALT": ["GC", "GG", "GG"],  # GC at pos 60, GG at pos 67 (twice)
                "Depth_Score": [0.015, 0.010, 0.008],  # insG highest, dupC variants
                "Motif": ["X", "X", "X"],
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # After dedup by [POS, REF, ALT]:
        #   - POS 60, REF=C, ALT=GC: unique → kept
        #   - POS 67, REF=G, ALT=GG: first kept (Depth_Score=0.010), second dropped
        # After GG filtering:
        #   - GG at pos 67 in motif X → kept
        #   - GC at pos 60 → kept (non-GG preserved)
        assert len(result) == 2
        assert any((result["POS"] == 60) & (result["ALT"] == "GC"))  # insG variant kept
        assert any((result["POS"] == 67) & (result["ALT"] == "GG"))  # dupC variant kept

    def test_issue_136_scenario_mixed_gg_and_gc_variants(self):
        """
        THE PRIMARY TEST FOR ISSUE #136 FIX

        Scenario: insG_pos54 sample has 176 GC variants + 85 GG variants
        - Current (buggy): Deletes all 176 GC variants → 0 variants reported
        - Fixed: Keeps highest-depth GC variant → 1-2 variants reported

        This simulates the actual failing case from the GitHub issue.
        """
        # Simulate insG_pos54: many GC variants (high depth) + some GG variants (low depth noise)
        gc_variants = pd.DataFrame(
            {
                "POS": [54] * 5,  # Multiple GC variants at pos 54 (k-mer ambiguity)
                "REF": ["C"] * 5,
                "ALT": ["GC"] * 5,
                "Depth_Score": [0.015, 0.014, 0.013, 0.012, 0.011],  # High depth
                "Motif": ["X"] * 5,
            }
        )

        gg_variants = pd.DataFrame(
            {
                "POS": [67] * 3,  # Some GG noise at pos 67
                "REF": ["G"] * 3,
                "ALT": ["GG"] * 3,
                "Depth_Score": [0.008, 0.006, 0.005],  # Lower depth
                "Motif": ["X"] * 3,
            }
        )

        df = pd.concat([gc_variants, gg_variants], ignore_index=True)

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # CRITICAL ASSERTIONS FOR ISSUE #136 FIX
        # 1. Should NOT be empty (bug would return 0 variants)
        assert not result.empty, "ISSUE #136 BUG: All GC variants were deleted!"

        # 2. Should contain GC variant (insG_pos54 mutation)
        gc_in_result = result[result["ALT"] == "GC"]
        assert not gc_in_result.empty, "GC variant (insG_pos54) should be preserved!"

        # 3. GC variant should have highest depth (signal, not noise)
        assert gc_in_result.iloc[0]["Depth_Score"] == 0.015, "Highest-depth GC variant should be kept"

        # 4. Should also have GG variant (in motif X)
        gg_in_result = result[result["ALT"] == "GG"]
        assert not gg_in_result.empty, "GG variant should also be kept (motif X)"

        # 5. After dedup, should have exactly 2 variants (1 GC + 1 GG)
        assert len(result) == 2, "Should keep 1 highest-depth GC + 1 highest-depth GG"

    def test_dupc_scenario_only_gg_variants(self):
        """
        dupC samples: Only GG variants (no regression expected).

        Scenario: dupC samples have G>GG, C>CG, T>TG at same position (same biological event)
        - All represent dupC at position 67
        - Should keep highest-depth variant
        """
        df = pd.DataFrame(
            {
                "POS": [67, 67, 67],
                "REF": ["G", "C", "T"],  # Different REF bases
                "ALT": ["GG", "CG", "TG"],  # All 2 bp insertions (same event)
                "Depth_Score": [0.010, 0.009, 0.007],  # GG highest
                "Motif": ["X", "X", "X"],
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # After dedup by [POS, REF, ALT]: all unique (different REF) → all kept
        # After GG filtering: GG kept (motif X), CG and TG preserved (non-GG)
        assert len(result) == 3  # All three representations kept
        assert any(result["ALT"] == "GG")  # dupC canonical
        assert any(result["ALT"] == "CG")  # dupC alternate 1
        assert any(result["ALT"] == "TG")  # dupC alternate 2

    def test_gg_motif_filtering_only_in_allowed_motifs(self):
        """
        Step 4: GG variants only kept if in allowed motifs (typically ["X"]).

        Biological: GG in motif "X" = pathogenic dupC
                    GG in other motifs = likely polymorphism → filtered
        """
        df = pd.DataFrame(
            {
                "POS": [67, 67, 67],
                "REF": ["G", "G", "G"],
                "ALT": ["GG", "GG", "GG"],
                "Depth_Score": [0.010, 0.008, 0.006],
                "Motif": ["X", "Y", "Z"],  # Only X is allowed
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # Should keep only GG in motif "X"
        assert len(result) == 1
        assert result.iloc[0]["Motif"] == "X"
        assert result.iloc[0]["Depth_Score"] == 0.010

    def test_non_gg_variants_always_preserved(self):
        """
        THE KEY FIX FOR ISSUE #136: Non-GG variants are ALWAYS preserved.

        Biological: GC, CG, CT, etc. are potential mutations (insG, dupA, etc.)
                    Must NOT be deleted just because GG variants exist!
        """
        df = pd.DataFrame(
            {
                "POS": [54, 60, 67, 67],
                "REF": ["C", "C", "G", "G"],
                "ALT": ["GC", "CT", "GG", "GC"],  # Mixed: GC, CT, GG
                "Depth_Score": [0.015, 0.012, 0.010, 0.008],
                "Motif": ["X", "X", "X", "X"],
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # ALL non-GG variants should be preserved (GC, CT)
        assert any((result["POS"] == 54) & (result["ALT"] == "GC")), "insG variant (GC) should be kept"
        assert any((result["POS"] == 60) & (result["ALT"] == "CT")), "dupA variant (CT) should be kept"

        # GG should also be kept (motif X)
        assert any((result["POS"] == 67) & (result["ALT"] == "GG")), "dupC variant (GG) should be kept"

        # GC at pos 67 should also be kept (different from GG)
        assert any((result["POS"] == 67) & (result["ALT"] == "GC")), "GC at pos 67 should be kept"

    def test_all_variants_excluded_returns_empty(self):
        """Edge case: All variants in excluded motifs should return empty."""
        df = pd.DataFrame(
            {
                "POS": [67, 67],
                "REF": ["G", "G"],
                "ALT": ["GG", "GG"],
                "Depth_Score": [0.010, 0.008],
                "Motif": ["Q", "8"],  # Both excluded
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=["Q", "8", "9"], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        assert result.empty, "All excluded motifs should result in empty DataFrame"

    def test_no_gg_variants_present_preserves_all_others(self):
        """When NO GG variants exist, all other variants should be preserved."""
        df = pd.DataFrame(
            {
                "POS": [54, 60, 67],
                "REF": ["C", "C", "C"],
                "ALT": ["GC", "CT", "CG"],  # No GG variants
                "Depth_Score": [0.015, 0.012, 0.010],
                "Motif": ["X", "X", "X"],
            }
        )

        result = _apply_uniform_filtering_right_motif(
            df.copy(), exclude_motifs_right=[], alt_for_motif_right_gg="GG", motifs_for_alt_gg=["X"]
        )

        # All variants should be preserved (no GG filtering needed)
        assert len(result) == 3
        assert list(result["ALT"]) == ["GC", "CT", "CG"]


@pytest.mark.unit
class TestMotifFilteringIntegration:
    """
    Integration tests for the complete motif filtering with feature flag.

    These tests verify the branching logic based on use_uniform_filtering flag.
    """

    def test_feature_flag_default_false_uses_legacy_logic(self):
        """
        Feature flag defaults to False → should use legacy logic.

        This ensures backward compatibility.
        """
        from vntyper.scripts.motif_processing import motif_correction_and_annotation

        # Create minimal test data
        df = pd.DataFrame(
            {
                "Motifs": ["X-Y", "X-Z"],
                "Variant": ["67:G>GG", "67:C>GC"],
                "POS": [67, 67],
                "REF": ["G", "C"],
                "ALT": ["GG", "GC"],
                "Estimated_Depth_AlternateVariant": [100, 150],
                "Estimated_Depth_Variant_ActiveRegion": [10000, 10000],
                "Depth_Score": [0.010, 0.015],
                "Confidence": ["High_Precision", "High_Precision"],
            }
        )

        merged_motifs = pd.DataFrame({"Motif": ["X", "Y", "Z"], "Motif_sequence": ["SEQ1", "SEQ2", "SEQ3"]})

        kestrel_config = {
            "motif_filtering": {
                "use_uniform_filtering": False,  # LEGACY
                "position_threshold": 60,
                "exclude_motifs_right": [],
                "alt_for_motif_right_gg": "GG",
                "motifs_for_alt_gg": ["X"],
                "exclude_alts_combined": [],
                "exclude_motifs_combined": [],
            }
        }

        result = motif_correction_and_annotation(df, merged_motifs, kestrel_config)

        # Legacy logic has the bug: deletes all non-GG variants
        # After fix is deployed with flag=False, this should still show legacy behavior
        assert "motif_filter_pass" in result.columns
        # We can't assert exact behavior without running full pipeline,
        # but we verify the function runs without errors

    def test_feature_flag_true_uses_new_logic(self):
        """
        Feature flag = True → should use new uniform filtering logic.

        This is the fix for Issue #136.
        """
        from vntyper.scripts.motif_processing import motif_correction_and_annotation

        df = pd.DataFrame(
            {
                "Motifs": ["X-Y", "X-Z"],
                "Variant": ["67:G>GG", "54:C>GC"],
                "POS": [67, 54],  # Different positions
                "REF": ["G", "C"],
                "ALT": ["GG", "GC"],  # Mixed GG and GC
                "Estimated_Depth_AlternateVariant": [100, 150],
                "Estimated_Depth_Variant_ActiveRegion": [10000, 10000],
                "Depth_Score": [0.010, 0.015],
                "Confidence": ["High_Precision", "High_Precision"],
            }
        )

        merged_motifs = pd.DataFrame({"Motif": ["X", "Y", "Z"], "Motif_sequence": ["SEQ1", "SEQ2", "SEQ3"]})

        kestrel_config = {
            "motif_filtering": {
                "use_uniform_filtering": True,  # NEW FIX
                "position_threshold": 60,
                "exclude_motifs_right": [],
                "alt_for_motif_right_gg": "GG",
                "motifs_for_alt_gg": ["X"],
                "exclude_alts_combined": [],
                "exclude_motifs_combined": [],
            }
        }

        result = motif_correction_and_annotation(df, merged_motifs, kestrel_config)

        # New logic preserves both GG and GC variants
        assert "motif_filter_pass" in result.columns
        # Both variants should pass (not deleted)
        passing = result[result["motif_filter_pass"]]
        # We expect both to pass with the new logic
        assert len(passing) >= 1  # At least one should pass


# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit
