# tests/unit/test_flagging.py

"""
Unit tests for the flagging module.
Validates flagging rules, condition evaluation, and duplicate detection
using vntyper/scripts/kestrel_config.json.
"""

import json
import logging
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


# --- Mutation-killing tests (Refs #179) ---
#
# Every test below was written against a surviving mutant recorded in
# docs/development/mutation-testing.md, and each was confirmed to fail with the
# mutation applied by hand and to pass with it reverted. The mutant it kills is named
# in the docstring so the link survives a later edit.


def _dup_config(sort_by=None, flag_name="Potential_Duplicate"):
    """Build a duplicate_flagging config, omitting 'sort_by' entirely when None."""
    config = {
        "enabled": True,
        "flag_name": flag_name,
        "group_by": ["REF", "ALT"],
    }
    if sort_by is not None:
        config["sort_by"] = sort_by
    return config


class TestEvaluateConditionErrorPath:
    """flagging.py:89 - the non-NameError error path must fail closed, not open."""

    def test_type_error_in_condition_returns_false(self):
        """A condition raising TypeError evaluates to False, never True (kills 89).

        The whole point of the except-Exception handler is that a malformed row
        cannot invent a flag. A non-numeric Depth_Score makes `<` raise TypeError,
        which is not a NameError and so reaches the second handler.
        """
        row = pd.Series({"Depth_Score": "n/a", "Motif": "2"})
        assert evaluate_condition(row, "Depth_Score < 0.4") is False

    def test_zero_division_in_condition_returns_false(self):
        """Any other evaluation error also evaluates to False (kills 89)."""
        row = pd.Series({"Depth_Score": 0.0, "Motif": "2"})
        assert evaluate_condition(row, "1 / Depth_Score > 1") is False

    def test_error_path_is_logged_at_error_level(self, caplog):
        """The failing condition is reported at ERROR, not merely DEBUG (kills 89)."""
        row = pd.Series({"Depth_Score": "n/a"})
        with caplog.at_level(logging.ERROR, logger="vntyper.scripts.flagging"):
            assert evaluate_condition(row, "Depth_Score < 0.4") is False
        errors = [r for r in caplog.records if r.levelno >= logging.ERROR]
        assert errors, "an unevaluatable condition must be logged at ERROR"

    def test_a_broken_rule_cannot_invent_a_flag(self):
        """End to end: a rule that cannot be evaluated leaves the row unflagged (kills 89)."""
        df = pd.DataFrame({"Depth_Score": ["n/a"], "Motif": ["2"], "REF": ["C"], "ALT": ["CG"]})
        result = add_flags(df, {"Low_Depth": "Depth_Score < 0.4"})
        assert result.loc[0, "Flag"] == "Not flagged"


class TestDuplicateSortSelection:
    """flagging.py:154/157 - which columns and directions the duplicate sort uses."""

    def test_configured_sort_by_is_used_instead_of_the_fallback(self):
        """A supplied sort_by wins over the Depth_Score fallback (kills 154).

        Position ascending and Depth_Score descending pick opposite primaries here,
        so a mutant that takes the fallback branch when sort_by *is* present flags
        the wrong row.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Position": [200, 100],
                "Depth_Score": [0.9, 0.1],
                "Motif": ["5", "5"],
            }
        )
        result = add_flags(df, {}, duplicates_config=_dup_config(sort_by=[{"column": "Position", "ascending": True}]))
        assert result.loc[1, "Flag"] == "Not flagged"
        assert result.loc[0, "Flag"] == "Potential_Duplicate"

    def test_missing_sort_by_falls_back_to_depth_score_descending(self):
        """With no sort_by, the highest Depth_Score is the primary record (kills 157).

        The fallback is `["Depth_Score"], [False]`; flipping that False to True would
        keep the *weakest* call and flag the strongest one as the duplicate.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.2],
                "Motif": ["5", "5"],
            }
        )
        result = add_flags(df, {}, duplicates_config=_dup_config())
        assert result.loc[0, "Flag"] == "Not flagged"
        assert result.loc[1, "Flag"] == "Potential_Duplicate"

    def test_empty_sort_by_list_also_falls_back(self):
        """An explicit empty sort_by takes the same fallback as a missing one (kills 154/157)."""
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.2],
                "Motif": ["5", "5"],
            }
        )
        result = add_flags(df, {}, duplicates_config=_dup_config(sort_by=[]))
        assert result.loc[0, "Flag"] == "Not flagged"
        assert result.loc[1, "Flag"] == "Potential_Duplicate"

    def test_priority_sort_actually_reorders_before_grouping(self):
        """The primary is chosen after sorting, not by input order (kills 217).

        The strongest call is the *second* row here, so a mutant that discards the
        sort (inplace=False) keeps the first row and flags the strongest call.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.2, 0.8],
                "Motif": ["5", "5"],
            }
        )
        result = add_flags(
            df, {}, duplicates_config=_dup_config(sort_by=[{"column": "Depth_Score", "ascending": False}])
        )
        assert result.loc[1, "Flag"] == "Not flagged"
        assert result.loc[0, "Flag"] == "Potential_Duplicate"


# --- #197: duplicate sort key falls back to Depth_Score alone ---


def test_duplicate_ordering_uses_depth_score_only_and_no_motif_column(kestrel_config):
    """#197: Motif is not a duplicate sort key; Depth_Score descending is.

    Specification. @hassansaei: "Fall back to the 1.3 Depth_Score-only rule
    [...] Do not use Motifs or Motif as a sort key." The previous config named
    the plural `Motifs`, which does not exist at step 6.5 and raised KeyError.
    """
    sort_by = kestrel_config["duplicate_flagging"]["sort_by"]

    assert sort_by == [{"column": "Depth_Score", "ascending": False}]
    assert all(entry["column"] not in ("Motif", "Motifs") for entry in sort_by)


def test_duplicate_flagging_stays_disabled_in_the_shipped_config(kestrel_config):
    """@hassansaei on #197: "Leave duplicate_flagging.enabled as false (as now)."."""
    assert kestrel_config["duplicate_flagging"]["enabled"] is False


def test_the_flagged_row_is_deterministic_when_depth_scores_tie():
    """A single-key sort_values is not stable; ties must not reorder.

    Implementation detail, not part of the #197 decision: with Depth_Score as the only
    sort key, pandas' default quicksort does not guarantee tied rows keep their input
    order, so which tied row is left unflagged (and therefore preferred by
    select_single_best_variant, which prefers unflagged rows) would be arbitrary.

    A handful of tied rows is not enough to observe this on the numpy build used here -
    argsort(kind="quicksort") happens to leave short equal-value runs untouched, and
    only starts reordering ties once the array is large enough (~260+ rows on this
    build) for numpy's introspective quicksort to switch strategy mid-sort. 300 rows
    is used to sit safely past that threshold and actually exercise the guarantee
    `kind="stable"` gives.

    Goes through add_flags with duplicate_flagging enabled, the path this would
    actually break in production if the toggle were ever turned on (it stays off in
    the shipped config; see test_duplicate_flagging_stays_disabled_in_the_shipped_config).
    """
    n = 300
    df = pd.DataFrame(
        {
            "REF": ["C"] * n,
            "ALT": ["CC"] * n,
            "Depth_Score": [0.5] * n,
            "POS": list(range(n)),
        }
    )
    duplicates_config = _dup_config(sort_by=[{"column": "Depth_Score", "ascending": False}])

    first = add_flags(df.copy(), {}, duplicates_config=duplicates_config)
    second = add_flags(df.copy(), {}, duplicates_config=duplicates_config)

    assert first["Flag"].tolist() == second["Flag"].tolist()
    assert first.loc[first["POS"] == 0, "Flag"].iloc[0] == "Not flagged"


class TestDuplicateSortColumnMustExist:
    """Characterisation of a real trap: a sort column absent from the frame raises.

    Before #197, `kestrel_config.json` shipped `duplicate_flagging.sort_by` naming the
    columns `Motifs` and `POS`. By the time flagging runs - step 6.5 of
    `process_kmer_results`, immediately after `motif_correction_and_annotation` - the
    `Motifs` column has been dropped: step 6 projects onto an explicit `keep_cols` list
    that carries `Motif` and `POS` but not `Motifs`.

    So flipping that toggle on with that stale `sort_by` would not silently mis-sort,
    it would raise `KeyError` from `sort_values` on the first frame it saw. That is the
    *good* failure mode and the opposite of AGENTS.md trap 3, where a config string
    naming a column that does not exist merely logs a warning and turns the rule off.

    #197 (@hassansaei) replaced that three-key `sort_by` with a single `Depth_Score`
    descending key - see `test_duplicate_ordering_uses_depth_score_only_and_no_motif_column`
    below - so the shipped config no longer has this shape. These tests reconstruct the
    retired three-key form explicitly via `_dup_config` instead of reading it off
    `kestrel_config.json`, so the trap stays documented for any future `sort_by` that
    reintroduces a since-dropped column name.
    """

    def test_unknown_sort_column_raises_rather_than_mis_sorting(self):
        """A sort_by naming a column the frame lacks raises KeyError, loudly."""
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.5],
                "Motif": ["5", "5"],
            }
        )
        config = _dup_config(sort_by=[{"column": "Motifs", "ascending": True}])
        with pytest.raises(KeyError):
            add_flags(df, {}, duplicates_config=config)

    def test_retired_three_key_sort_by_would_still_raise(self):
        """The pre-#197 sort_by shape (Depth_Score, Motifs, POS) still raises KeyError.

        #197 removed this exact three-key shape from the shipped config (it is no
        longer readable from `kestrel_config.json`), so it is reconstructed here
        explicitly to show the trap it documents is a property of the code, not of
        whatever the config currently happens to ship.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.5],
                "POS": [100, 200],
                "Motif": ["5", "5"],
            }
        )
        stale_sort_by = [
            {"column": "Depth_Score", "ascending": False},
            {"column": "Motifs", "ascending": True},
            {"column": "POS", "ascending": True},
        ]
        config = _dup_config(sort_by=stale_sort_by)
        with pytest.raises(KeyError):
            add_flags(df, {}, duplicates_config=config)


class TestDuplicateFlagCombination:
    """flagging.py:244/250 - how an existing flag and a duplicate flag combine."""

    def test_existing_flag_is_kept_and_extended_on_a_duplicate_row(self):
        """A duplicate that already carries a flag keeps both (kills 244 and both 250s).

        The '+' concatenation on line 250 is the only place the existing flag and the
        duplicate flag are joined; replacing either '+' with '-' raises TypeError, and
        inverting the '==' on line 244 drops the existing flag entirely.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.5],
                "Motif": ["2", "2"],
            }
        )
        result = add_flags(
            df,
            {"Low_Depth": "Depth_Score < 0.9"},
            duplicates_config=_dup_config(sort_by=[{"column": "Depth_Score", "ascending": False}]),
        )
        assert result.loc[0, "Flag"] == "Low_Depth"
        assert result.loc[1, "Flag"] == "Low_Depth, Potential_Duplicate"

    def test_unflagged_duplicate_gets_only_the_duplicate_flag(self):
        """'Not flagged' is replaced, never appended to (kills 244).

        Inverting the '==' produces the literal string 'Not flagged, Potential_Duplicate',
        which reads as both flagged and unflagged at once.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.5],
                "Motif": ["5", "5"],
            }
        )
        result = add_flags(
            df, {}, duplicates_config=_dup_config(sort_by=[{"column": "Depth_Score", "ascending": False}])
        )
        assert result.loc[0, "Flag"] == "Not flagged"
        assert result.loc[1, "Flag"] == "Potential_Duplicate"


class TestDuplicateMarkingLeavesNoTrace:
    """flagging.py:258/259 - the frame handed back is the caller's frame plus 'Flag'."""

    def test_original_row_order_is_restored(self):
        """Duplicate marking must not reorder the caller's rows (kills 258).

        Sorting by Depth_Score descending permutes these three rows, so a mutant that
        discards the restoring sort returns them in scoring order. That order reaches
        kestrel_pre_result.tsv, where row order is the debugging aid.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C", "C"],
                "ALT": ["CG", "CG", "CG"],
                "Depth_Score": [0.2, 0.8, 0.5],
                "Motif": ["5", "5", "5"],
            }
        )
        result = add_flags(
            df, {}, duplicates_config=_dup_config(sort_by=[{"column": "Depth_Score", "ascending": False}])
        )
        assert list(result.index) == [0, 1, 2]
        assert list(result["Depth_Score"]) == [0.2, 0.8, 0.5]

    def test_internal_bookkeeping_columns_are_dropped(self):
        """The three scratch columns never reach the caller (kills 259).

        They would otherwise be written straight into the output TSV.
        """
        df = pd.DataFrame(
            {
                "REF": ["C", "C"],
                "ALT": ["CG", "CG"],
                "Depth_Score": [0.8, 0.5],
                "Motif": ["5", "5"],
            }
        )
        result = add_flags(
            df, {}, duplicates_config=_dup_config(sort_by=[{"column": "Depth_Score", "ascending": False}])
        )
        assert "__original_index" not in result.columns
        assert "__dup_indicator" not in result.columns
        assert "__is_duplicate" not in result.columns
        assert list(result.columns) == ["REF", "ALT", "Depth_Score", "Motif", "Flag"]
