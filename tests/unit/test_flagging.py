#!/usr/bin/env python3
# tests/unit/test_flagging.py

"""
Unit tests for the flagging module.
Validates flagging rules, condition evaluation, and duplicate detection
using vntyper/scripts/kestrel_config.json.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts.flagging import add_flags, evaluate_condition, regex_match

pytestmark = pytest.mark.unit

# Conserved motifs defined in kestrel_config.json flagging_rules
CONSERVED_MOTIFS = ["1", "2", "3", "4", "6", "7", "8", "9"]
NON_CONSERVED_MOTIFS = ["5", "6p", "D", "E", "Q", "X"]


@pytest.fixture(scope="session")
def kestrel_config():
    """Load the Kestrel configuration from vntyper/scripts/kestrel_config.json."""
    this_file = Path(__file__).resolve()
    config_path = this_file.parents[2] / "vntyper" / "scripts" / "kestrel_config.json"
    if not config_path.exists():
        pytest.exit(f"kestrel_config.json not found at {config_path}", returncode=1)
    with config_path.open("r") as f:
        return json.load(f)


@pytest.fixture
def flagging_rules(kestrel_config):
    """Extract flagging rules from the Kestrel config."""
    return kestrel_config["flagging_rules"]


# --- regex_match tests ---


class TestRegexMatch:
    def test_simple_match(self):
        assert regex_match(r"^D", "D5") is True

    def test_no_match(self):
        assert regex_match(r"^D", "E5") is False

    def test_non_string_value(self):
        assert regex_match(r"^\d+", 42) is True

    def test_invalid_pattern_returns_false(self):
        assert regex_match(r"[invalid", "test") is False


# --- evaluate_condition tests ---


class TestEvaluateCondition:
    def test_simple_comparison(self):
        row = pd.Series({"Depth_Score": 0.3, "Motif": "2"})
        assert evaluate_condition(row, "Depth_Score < 0.4") is True

    def test_in_operator(self):
        row = pd.Series({"Motif": "2"})
        assert evaluate_condition(row, "Motif in ['1', '2', '3']") is True

    def test_not_in_operator(self):
        row = pd.Series({"Motif": "5"})
        assert evaluate_condition(row, "Motif in ['1', '2', '3']") is False

    def test_missing_column_returns_false(self):
        row = pd.Series({"Motif": "2"})
        assert evaluate_condition(row, "NonExistent > 0") is False

    def test_combined_and_condition(self):
        row = pd.Series({"Depth_Score": 0.3, "Motif": "2"})
        condition = "(Depth_Score < 0.4) and (Motif in ['1', '2', '3'])"
        assert evaluate_condition(row, condition) is True

    def test_combined_and_condition_fails_one(self):
        row = pd.Series({"Depth_Score": 0.5, "Motif": "2"})
        condition = "(Depth_Score < 0.4) and (Motif in ['1', '2', '3'])"
        assert evaluate_condition(row, condition) is False

    def test_regex_match_in_condition(self):
        row = pd.Series({"REF": "C", "ALT": "CGGCA"})
        condition = "regex_match('^C', REF) and ALT == 'CGGCA'"
        assert evaluate_condition(row, condition) is True

    def test_pd_na_in_list_returns_false(self):
        """pd.NA in an 'in' check should return False, not raise TypeError (#154)."""
        row = pd.Series({"Depth_Score": 0.015, "Motif": pd.NA})
        condition = "(Depth_Score < 0.4) and (Motif in ['1', '2', '3'])"
        assert evaluate_condition(row, condition) is False

    def test_pd_na_equality_returns_false(self):
        """pd.NA in an '==' check should return False, not raise TypeError."""
        row = pd.Series({"REF": pd.NA, "ALT": "CG"})
        assert evaluate_condition(row, "REF == 'C'") is False

    def test_none_in_list_returns_false(self):
        """Explicit None in an 'in' check should also return False."""
        row = pd.Series({"Motif": None, "Depth_Score": 0.1})
        condition = "(Depth_Score < 0.4) and (Motif in ['1', '2'])"
        assert evaluate_condition(row, condition) is False


# --- Low_Depth_Conserved_Motifs flagging rule tests ---


class TestLowDepthConservedMotifsFlagging:
    """Tests for the Low_Depth_Conserved_Motifs rule with threshold 0.4."""

    def test_conserved_motif_below_threshold_flagged(self, flagging_rules):
        """Depth_Score=0.3 in conserved motif '2' should be flagged."""
        df = pd.DataFrame({"Depth_Score": [0.3], "Motif": ["2"], "REF": ["C"], "ALT": ["CG"]})
        result = add_flags(df, flagging_rules)
        assert "Low_Depth_Conserved_Motifs" in result.loc[0, "Flag"]

    def test_conserved_motif_above_threshold_not_flagged(self, flagging_rules):
        """Depth_Score=0.5 in conserved motif '2' should NOT be flagged."""
        df = pd.DataFrame({"Depth_Score": [0.5], "Motif": ["2"], "REF": ["C"], "ALT": ["CG"]})
        result = add_flags(df, flagging_rules)
        assert "Low_Depth_Conserved_Motifs" not in result.loc[0, "Flag"]

    def test_boundary_just_below_threshold(self, flagging_rules):
        """Depth_Score=0.39 (< 0.4) in conserved motif should be flagged."""
        df = pd.DataFrame({"Depth_Score": [0.39], "Motif": ["1"], "REF": ["C"], "ALT": ["CG"]})
        result = add_flags(df, flagging_rules)
        assert "Low_Depth_Conserved_Motifs" in result.loc[0, "Flag"]

    def test_boundary_at_threshold_not_flagged(self, flagging_rules):
        """Depth_Score=0.4 (not < 0.4) in conserved motif should NOT be flagged."""
        df = pd.DataFrame({"Depth_Score": [0.4], "Motif": ["1"], "REF": ["C"], "ALT": ["CG"]})
        result = add_flags(df, flagging_rules)
        assert "Low_Depth_Conserved_Motifs" not in result.loc[0, "Flag"]

    @pytest.mark.parametrize("motif", CONSERVED_MOTIFS)
    def test_each_conserved_motif_flagged(self, flagging_rules, motif):
        """Each conserved motif should be flagged when Depth_Score is below threshold."""
        df = pd.DataFrame({"Depth_Score": [0.1], "Motif": [motif], "REF": ["C"], "ALT": ["CG"]})
        result = add_flags(df, flagging_rules)
        assert "Low_Depth_Conserved_Motifs" in result.loc[0, "Flag"]

    @pytest.mark.parametrize("motif", NON_CONSERVED_MOTIFS)
    def test_non_conserved_motif_not_flagged(self, flagging_rules, motif):
        """Non-conserved motifs should NOT be flagged regardless of Depth_Score."""
        df = pd.DataFrame({"Depth_Score": [0.01], "Motif": [motif], "REF": ["C"], "ALT": ["CG"]})
        result = add_flags(df, flagging_rules)
        assert "Low_Depth_Conserved_Motifs" not in result.loc[0, "Flag"]


# --- False_Positive_4bp_Insertion flagging rule tests ---


class TestFalsePositive4bpInsertionFlagging:
    def test_matching_variant_flagged(self, flagging_rules):
        """REF='C', ALT='CGGCA' should be flagged as False_Positive_4bp_Insertion."""
        df = pd.DataFrame({"REF": ["C"], "ALT": ["CGGCA"], "Depth_Score": [0.5], "Motif": ["5"]})
        result = add_flags(df, flagging_rules)
        assert "False_Positive_4bp_Insertion" in result.loc[0, "Flag"]

    def test_non_matching_variant_not_flagged(self, flagging_rules):
        """Different ALT should NOT be flagged as False_Positive_4bp_Insertion."""
        df = pd.DataFrame({"REF": ["C"], "ALT": ["CG"], "Depth_Score": [0.5], "Motif": ["5"]})
        result = add_flags(df, flagging_rules)
        assert "False_Positive_4bp_Insertion" not in result.loc[0, "Flag"]


# --- Multiple flags and "Not flagged" tests ---


class TestMultipleFlags:
    def test_variant_matching_both_rules(self, flagging_rules):
        """A variant matching both rules should have both flags comma-separated."""
        df = pd.DataFrame({"REF": ["C"], "ALT": ["CGGCA"], "Depth_Score": [0.1], "Motif": ["2"]})
        result = add_flags(df, flagging_rules)
        flag_value = result.loc[0, "Flag"]
        assert "False_Positive_4bp_Insertion" in flag_value
        assert "Low_Depth_Conserved_Motifs" in flag_value

    def test_no_flags_applied(self, flagging_rules):
        """A variant matching no rules should get 'Not flagged'."""
        df = pd.DataFrame({"REF": ["A"], "ALT": ["AT"], "Depth_Score": [0.5], "Motif": ["5"]})
        result = add_flags(df, flagging_rules)
        assert result.loc[0, "Flag"] == "Not flagged"

    def test_empty_dataframe(self, flagging_rules):
        """Empty input should return empty output with Flag column."""
        df = pd.DataFrame(columns=["REF", "ALT", "Depth_Score", "Motif"])
        result = add_flags(df, flagging_rules)
        assert result.empty
        assert "Flag" in result.columns


# --- Duplicate flagging tests ---


class TestDuplicateFlagging:
    def test_duplicate_rows_flagged(self):
        """Duplicate rows (same REF, ALT) should have the lower-priority one flagged."""
        df = pd.DataFrame(
            {
                "REF": ["C", "C", "A"],
                "ALT": ["CG", "CG", "AT"],
                "Depth_Score": [0.8, 0.5, 0.6],
                "Motif": ["5", "5", "5"],
            }
        )
        duplicates_config = {
            "enabled": True,
            "flag_name": "Potential_Duplicate",
            "group_by": ["REF", "ALT"],
            "sort_by": [{"column": "Depth_Score", "ascending": False}],
        }
        result = add_flags(df, {}, duplicates_config=duplicates_config)
        assert "Potential_Duplicate" not in result.loc[0, "Flag"]
        assert "Potential_Duplicate" in result.loc[1, "Flag"]
        assert "Potential_Duplicate" not in result.loc[2, "Flag"]

    def test_duplicates_disabled(self):
        """When duplicates config is disabled, no duplicate flags should be added."""
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.5],
                "Motif": ["5", "5"],
            }
        )
        duplicates_config = {"enabled": False}
        result = add_flags(df, {}, duplicates_config=duplicates_config)
        assert all(result["Flag"] == "Not flagged")
