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

from vntyper.scripts.flagging import add_artifact_gate, add_flags, evaluate_condition, regex_match

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


def test_invalid_regex_is_false_and_observable(caplog: pytest.LogCaptureFixture) -> None:
    """An invalid configured regex fails closed and leaves an ERROR breadcrumb."""
    with caplog.at_level(logging.ERROR, logger="vntyper.scripts.flagging"):
        assert regex_match("[", "value") is False
    assert "Error in regex_match" in caplog.text


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

    def test_missing_column_returns_false(self, caplog: pytest.LogCaptureFixture) -> None:
        """A missing condition name fails closed and records its name at WARNING."""
        row = pd.Series({"Motif": "2"})
        with caplog.at_level(logging.WARNING, logger="vntyper.scripts.flagging"):
            assert evaluate_condition(row, "NonExistent > 0") is False
        matching_records = [record for record in caplog.records if "NonExistent" in record.getMessage()]
        assert len(matching_records) == 1
        assert matching_records[0].levelno == logging.WARNING

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

    def test_error_path_is_logged_at_error_level(self, caplog: pytest.LogCaptureFixture) -> None:
        """The failing condition is reported at ERROR, not merely DEBUG (kills 89)."""
        row = pd.Series({"Depth_Score": "n/a"})
        with caplog.at_level(logging.ERROR, logger="vntyper.scripts.flagging"):
            assert evaluate_condition(row, "Depth_Score < 0.4") is False
        assert "Error evaluating condition" in caplog.text

    def test_evaluation_error_is_false_and_observable(self, caplog: pytest.LogCaptureFixture) -> None:
        """A non-NameError condition failure stays false and logs at ERROR."""
        row = pd.Series({"Depth_Score": "n/a"})
        with caplog.at_level(logging.ERROR, logger="vntyper.scripts.flagging"):
            assert evaluate_condition(row, "Depth_Score < 0.4") is False
        matching_records = [record for record in caplog.records if "Error evaluating condition" in record.getMessage()]
        assert len(matching_records) == 1
        assert matching_records[0].levelno == logging.ERROR

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
    [...] Do not use Motifs or Motif as a sort key."

    This is a **policy** change -- which record wins when two calls share a REF/ALT --
    not a repair for a broken column reference. An earlier revision of this docstring
    claimed the retired three-key `sort_by` named a column that "does not exist at step
    6.5 and raised KeyError". It did exist:
    `test_the_frame_flagging_sees_still_carries_the_motifs_column` below measures the
    frame `add_flags` is actually handed and finds `Motifs`, `POS` and `Depth_Score` all
    present, so the retired shape would have sorted, not raised. The reason to drop
    `Motifs` from the key is that motif identity is not a duplicate-priority criterion,
    which is what @hassansaei decided; it is not that the column had disappeared.
    """
    sort_by = kestrel_config["duplicate_flagging"]["sort_by"]

    assert sort_by == [{"column": "Depth_Score", "ascending": False}]
    assert all(entry["column"] not in ("Motif", "Motifs") for entry in sort_by)


def test_the_frame_flagging_sees_still_carries_the_motifs_column(kestrel_config):
    """`Motifs` survives step 6 and is present when `add_flags` runs at step 6.5.

    The claim two docstrings in this module used to make -- that
    `motif_correction_and_annotation` drops `Motifs` -- is false, and this measures it
    rather than reasoning about it. The narrowing to `keep_cols` is applied to the
    *working* frames (`motif_left` / `motif_right`), not to what the function returns:
    it deep-copies its input into `original_df`, does the motif work on the narrow
    frames, copies the annotations back by index, and returns `original_df` at the
    original width. `process_kmer_results` then hands exactly that frame to `add_flags`.

    Asserted for `POS` and `Depth_Score` too, because the retired three-key `sort_by`
    named all three; every one of them is present.
    """
    from vntyper.scripts.motif_processing import motif_correction_and_annotation

    df = pd.DataFrame(
        {
            "Motifs": ["5-9", "5-9"],
            "Variant": ["v1", "v2"],
            "POS": [10, 70],
            "REF": ["C", "C"],
            "ALT": ["CG", "CG"],
            "Motif_sequence": ["ACGT", "ACGT"],
            "Depth_Score": [0.8, 0.5],
            "Confidence": ["High_Precision", "High_Precision"],
            "Estimated_Depth_AlternateVariant": [10, 10],
            "Estimated_Depth_Variant_ActiveRegion": [100, 100],
            "is_frameshift": [True, True],
            "is_valid_frameshift": [True, True],
        }
    )
    merged_motifs = pd.DataFrame({"Motif": ["9", "5"], "Motif_sequence": ["ACGT", "ACGT"]})

    annotated = motif_correction_and_annotation(df, merged_motifs, kestrel_config)

    assert {"Motifs", "POS", "Depth_Score"} <= set(annotated.columns)
    assert len(annotated) == len(df), "step 6 marks rows, it does not drop them"

    # And the retired three-key sort_by therefore sorts this frame rather than raising.
    stale_sort_by = [
        {"column": "Depth_Score", "ascending": False},
        {"column": "Motifs", "ascending": True},
        {"column": "POS", "ascending": True},
    ]
    flagged = add_flags(annotated, {}, duplicates_config=_dup_config(sort_by=stale_sort_by))

    assert "Flag" in flagged.columns


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
    """A sort column absent from the frame raises, loudly. A property of `add_flags`.

    `sort_values` raises `KeyError` for a column the frame does not carry, so a
    `duplicate_flagging.sort_by` naming one fails fast rather than silently mis-sorting.
    That is the *good* failure mode and the opposite of AGENTS.md trap 3, where a config
    string naming a missing column merely logs a warning and turns the rule off.

    **The frames below omit `Motifs` deliberately, and that is a synthetic condition, not
    a reproduction of production.** An earlier revision of this docstring justified them
    by claiming step 6 had dropped `Motifs` before flagging runs, so that the pre-#197
    three-key `sort_by` (`Depth_Score`, `Motifs`, `POS`) would have raised on the first
    real frame it saw. That is false, and
    `test_the_frame_flagging_sees_still_carries_the_motifs_column` above measures it:
    `motif_correction_and_annotation` returns its input frame at the original width, so
    `Motifs`, `POS` and `Depth_Score` are all present at step 6.5 and the retired shape
    would have sorted normally.

    #197 (@hassansaei) replaced that three-key `sort_by` with a single `Depth_Score`
    descending key because motif identity is not a duplicate-priority criterion - an
    authorised policy change, not a repair. The shipped config no longer has the
    three-key shape, so these tests reconstruct it via `_dup_config` rather than reading
    it off `kestrel_config.json`; what they document is the `add_flags` behaviour that
    would catch any *future* `sort_by` naming a column its frame really does lack.
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

    def test_retired_three_key_sort_by_raises_on_a_frame_that_lacks_motifs(self):
        """The pre-#197 sort_by shape (Depth_Score, Motifs, POS) raises on THIS frame.

        The frame is built without `Motifs` on purpose. A production frame at step 6.5
        carries it (see `test_the_frame_flagging_sees_still_carries_the_motifs_column`),
        so this is not a reproduction of a bug the shipped config had - it pins that
        `add_flags` fails fast rather than mis-sorting, using the retired shape only
        because it is the multi-key example already at hand.
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


# --- #174: the artifact gate ---
#
# `add_artifact_gate` splits the flags into two classes. Most flags are advisory: they
# say "this call warrants scrutiny", and their only effect is to deprioritise the row
# during selection. A few name a known technical artifact - a row that is not a candidate
# variant at all - and those must not reach `kestrel_result.tsv`, where a consumer would
# read them as a call.
#
# Which flags are artifacts is configuration, never a name written into Python. Every test
# below passes the list in explicitly for that reason; only
# `test_the_shipped_config_declares_the_4bp_insertion_as_an_artifact` reads the shipped
# value, and it is the single place the name appears.

ARTIFACT_FLAG = "False_Positive_4bp_Insertion"


class TestArtifactGate:
    """`flag_filter_pass`: False iff the row carries a declared artifact flag."""

    def test_a_declared_artifact_flag_fails_the_gate(self):
        """The whole point: a declared artifact is marked for exclusion."""
        df = pd.DataFrame({"Flag": [ARTIFACT_FLAG]})

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].tolist() == [False]

    def test_an_advisory_flag_does_not_fail_the_gate(self):
        """#174 excludes artifacts, not flags.

        `Low_Depth_Conserved_Motifs` marks a call that warrants scrutiny, not one that is
        technically impossible. Excluding it would delete real low-depth calls, which is
        the failure mode this milestone exists to prevent.
        """
        df = pd.DataFrame({"Flag": ["Low_Depth_Conserved_Motifs", "Not flagged"]})

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].tolist() == [True, True]

    def test_a_row_carrying_both_an_artifact_and_an_advisory_flag_is_excluded(self):
        """`Flag` is comma-joined, so membership is per element."""
        df = pd.DataFrame({"Flag": [f"Low_Depth_Conserved_Motifs, {ARTIFACT_FLAG}"]})

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].tolist() == [False]

    def test_the_artifact_may_appear_in_any_position_of_the_joined_value(self):
        """Order within `Flag` is incidental; `add_flags` emits rules in config order."""
        df = pd.DataFrame({"Flag": [f"{ARTIFACT_FLAG}, Potential_Duplicate"]})

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].tolist() == [False]

    def test_a_flag_name_is_matched_whole_and_not_as_a_substring(self):
        """A substring test would exclude `..._Insertion_v2` given `..._Insertion`."""
        df = pd.DataFrame({"Flag": [f"{ARTIFACT_FLAG}_v2"]})

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].tolist() == [True]

    def test_a_declared_artifact_that_is_a_substring_of_a_present_flag_does_not_match(self):
        """The mirror of the case above: the artifact is the shorter of the two names."""
        df = pd.DataFrame({"Flag": ["Low_Depth_Conserved_Motifs"]})

        assert add_artifact_gate(df, ["Low_Depth"])["flag_filter_pass"].tolist() == [True]

    def test_a_frame_without_a_flag_column_still_gains_the_gate(self):
        """A negative run's frame legitimately carries no `Flag` (report_formatting.py:68).

        The gate must still be present, or `filter_final_dataframe` raises on a missing
        required column - turning a safety gate into an abort.
        """
        out = add_artifact_gate(pd.DataFrame({"POS": [1, 2]}), [ARTIFACT_FLAG])

        assert out["flag_filter_pass"].tolist() == [True, True]

    @pytest.mark.parametrize("value", [None, float("nan"), pd.NA])
    def test_an_unknown_flag_value_is_not_a_declared_artifact(self, value):
        """`Flag=None` is an accepted input - test_haplo_count_and_selection.py:398 relies on it.

        Selection still deprioritises it. But absence of evidence is not evidence of an
        artifact, so it passes this gate rather than being excluded, and a bare `.split()`
        would have raised.
        """
        df = pd.DataFrame({"Flag": [value]})

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].tolist() == [True]

    def test_an_all_nan_float_flag_column_passes_rather_than_raising(self):
        """The column-wide version of the case above.

        A frame whose `Flag` column is entirely missing values is typed float64 by pandas,
        not object, so the gate meets `numpy.float64` NaN rather than `None`. Elementwise
        coercion is what keeps that from raising.
        """
        df = pd.DataFrame({"Flag": [float("nan"), float("nan")]})
        assert df["Flag"].dtype == "float64", "the premise of this test is the float64 dtype"

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].tolist() == [True, True]

    def test_an_empty_frame_gains_the_gate_column(self):
        """`filter_final_dataframe` raises on a missing gate; an empty frame is no exception."""
        assert "flag_filter_pass" in add_artifact_gate(pd.DataFrame(), [ARTIFACT_FLAG]).columns

    def test_an_empty_frame_that_does_carry_a_flag_column_also_gains_the_gate(self):
        """The other empty shape: zero rows but the full column set."""
        out = add_artifact_gate(pd.DataFrame({"Flag": pd.Series(dtype=object)}), [ARTIFACT_FLAG])

        assert out["flag_filter_pass"].tolist() == []

    def test_the_gate_is_a_boolean_column(self):
        """`filter_final_dataframe` ANDs it into a boolean mask, so the dtype is contractual."""
        df = pd.DataFrame({"Flag": [ARTIFACT_FLAG, "Not flagged"]})

        assert add_artifact_gate(df, [ARTIFACT_FLAG])["flag_filter_pass"].dtype == bool

    def test_an_empty_artifact_list_excludes_nothing(self):
        """Emptying `artifact_flags` in config restores the pre-#174 behaviour, with no code change.

        This is the reversibility property #174 rests on: the decision lives in
        `kestrel_config.json`, so narrowing or withdrawing it is a config edit.
        """
        df = pd.DataFrame({"Flag": [ARTIFACT_FLAG, "Low_Depth_Conserved_Motifs"]})

        assert add_artifact_gate(df, [])["flag_filter_pass"].tolist() == [True, True]

    def test_several_declared_artifacts_are_all_honoured(self):
        """The config key is a list, so a second entry must gate as well as the first."""
        df = pd.DataFrame({"Flag": [ARTIFACT_FLAG, "Another_Artifact", "Not flagged"]})

        gate = add_artifact_gate(df, [ARTIFACT_FLAG, "Another_Artifact"])["flag_filter_pass"]

        assert gate.tolist() == [False, False, True]

    def test_the_callers_frame_is_not_mutated(self):
        """Stages return copies; a leaked column would reach `kestrel_pre_result.tsv` twice."""
        df = pd.DataFrame({"Flag": [ARTIFACT_FLAG]})

        add_artifact_gate(df, [ARTIFACT_FLAG])

        assert "flag_filter_pass" not in df.columns

    def test_the_other_columns_and_the_row_order_survive(self):
        """The gate marks; it does not drop, reorder or project."""
        df = pd.DataFrame({"POS": [67, 54], "Flag": [ARTIFACT_FLAG, "Not flagged"]})

        out = add_artifact_gate(df, [ARTIFACT_FLAG])

        assert list(out.columns) == ["POS", "Flag", "flag_filter_pass"]
        assert out["POS"].tolist() == [67, 54]

    def test_the_flag_column_itself_is_left_untouched(self):
        """The evidence stays readable in `kestrel_pre_result.tsv` beside the False gate."""
        df = pd.DataFrame({"Flag": [f"Low_Depth_Conserved_Motifs, {ARTIFACT_FLAG}"]})

        out = add_artifact_gate(df, [ARTIFACT_FLAG])

        assert out["Flag"].iloc[0] == f"Low_Depth_Conserved_Motifs, {ARTIFACT_FLAG}"


def test_the_shipped_config_declares_the_4bp_insertion_as_an_artifact(kestrel_config):
    """Ties the code to configuration: the flag name is never written inline in Python.

    The declared artifact must also be a flag some rule can actually raise, or the gate
    would be inert - a typo in either half is invisible without this second assertion.
    """
    assert kestrel_config["artifact_flags"] == ["False_Positive_4bp_Insertion"]
    assert "False_Positive_4bp_Insertion" in kestrel_config["flagging_rules"]


def test_no_artifact_flag_name_is_written_into_the_flagging_module(kestrel_config):
    """#174's non-negotiable constraint, asserted rather than trusted.

    `add_artifact_gate` reads the list from `kestrel_config.json`. If a flag name ever
    appears in the module's source, the reversibility property above is a fiction:
    emptying `artifact_flags` would no longer restore the previous behaviour.
    """
    source = (Path(__file__).resolve().parents[2] / "vntyper" / "scripts" / "flagging.py").read_text(encoding="utf-8")

    for flag in kestrel_config["artifact_flags"]:
        assert flag not in source, f"{flag!r} is hardcoded in flagging.py; it must come from configuration"


def test_a_string_artifact_flags_value_is_refused_rather_than_iterated():
    """`"artifact_flags": "Foo"` is valid JSON and an easy thing to write instead of a list.

    A bare string satisfies `Sequence[str]`, and `set("Foo")` is a set of *characters* - so
    the gate would silently degrade to matching single letters and let every artifact
    through. The whole point of #174 is that a known artifact is not reported as a call, so
    a config typo must not quietly undo it. Found by adversarial review of the PR.
    """
    df = pd.DataFrame({"Flag": ["False_Positive_4bp_Insertion"]})

    with pytest.raises(ValueError, match="must be a list of flag names"):
        add_artifact_gate(df, "False_Positive_4bp_Insertion")


def test_an_artifact_flag_containing_a_comma_is_refused():
    """`Flag` is comma-joined, so such a name could never match and would be inert."""
    df = pd.DataFrame({"Flag": ["Artifact,Comma"]})

    with pytest.raises(ValueError, match="must not contain a comma"):
        add_artifact_gate(df, ["Artifact,Comma"])


def test_the_shipped_artifact_flags_survive_both_guards():
    """The guards must not reject the real configuration."""
    config = json.loads(Path("vntyper/scripts/kestrel_config.json").read_text(encoding="utf-8"))
    df = pd.DataFrame({"Flag": ["False_Positive_4bp_Insertion", "Not flagged"]})

    out = add_artifact_gate(df, config["artifact_flags"])

    assert out["flag_filter_pass"].tolist() == [False, True]
