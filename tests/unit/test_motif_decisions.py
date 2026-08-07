"""The motif filtering decision layer, and the end-to-end oracle that guards it.

``motif_correction_and_annotation`` was 182 lines holding 28 of
``motif_processing.py``'s 46 surviving mutants, because it fused the config
read, the column plumbing, the left/right positional split and the motif
filtering decisions into one function that could not be called without a full
frame. This file pins the whole of it end to end (``TestTheEndToEndOracle``)
and then exercises the extracted pure functions one rule at a time.

The oracle came first and is the invariant the extraction had to preserve: it
was written and recorded against the pre-extraction code, and it passed
unmodified afterwards. If it ever needs editing alongside a change to the
decision layer, that change is not a refactor.

Everything here is CHARACTERISATION unless a test's own docstring says
otherwise.
"""

import json
import logging
from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts.motif_decisions import (
    apply_combined_exclusions,
    apply_gg_alt_rule,
    apply_right_motif_exclusions,
    has_gg_alternate,
    split_left_right,
)
from vntyper.scripts.motif_processing import motif_correction_and_annotation

pytestmark = pytest.mark.unit

CONFIG_PATH = Path("vntyper/scripts/kestrel_config.json")


def _shipped_config():
    """The real ``kestrel_config.json``, not a hand-written stand-in.

    The decisions under test are driven entirely by ``motif_filtering``, so a
    fixture config would pin a policy nobody ships.
    """
    return json.loads(CONFIG_PATH.read_text())


def _merged_motifs():
    """The motif annotation table ``motif_correction_and_annotation`` merges on."""
    motifs = ["1", "2", "3", "4", "5", "6", "7", "8", "X", "Y", "Z"]
    return pd.DataFrame({"Motif": motifs, "Motif_sequence": [f"SEQ_{m}" for m in motifs]})


def _kestrel_frame(**overrides):
    """A frame shaped like the one ``process_kmer_results`` hands over.

    Covers, in row order: a left motif (POS < 60), a plain right motif, a right
    motif on the ``exclude_motifs_right`` list, a GG alternate, an ALT on
    ``exclude_alts_combined``, and a motif on ``exclude_motifs_combined``.
    """
    frame = pd.DataFrame(
        {
            "Motifs": ["1-2", "X-3", "8-4", "Y-5", "Z-6", "6-7"],
            "Variant": ["Insertion", "Insertion", "Insertion", "Insertion", "Insertion", "Deletion"],
            "POS": [54, 67, 67, 67, 68, 55],
            "REF": ["C", "G", "C", "G", "C", "C"],
            "ALT": ["CC", "GC", "CG", "GG", "CCGCC", "CT"],
            "Estimated_Depth_AlternateVariant": [50, 51, 52, 53, 54, 55],
            "Estimated_Depth_Variant_ActiveRegion": [5000] * 6,
            "Depth_Score": [0.010, 0.011, 0.012, 0.013, 0.014, 0.015],
            "Confidence": ["High_Precision"] * 6,
        }
    )
    for column, values in overrides.items():
        frame[column] = values
    return frame


# ---------------------------------------------------------------------------
# The end-to-end oracle
# ---------------------------------------------------------------------------


class TestTheEndToEndOracle:
    """What the shipped config does to a frame covering every filtering rule.

    Recorded before the decision layer was extracted and unchanged since. These
    are not claims that the policy is right - they are the byte-for-byte record
    of what it currently is.
    """

    def test_the_exact_pass_fail_verdict_for_every_rule(self):
        out = motif_correction_and_annotation(_kestrel_frame(), _merged_motifs(), _shipped_config())

        assert out["motif_filter_pass"].tolist() == [
            True,  # 1-2  POS 54  left motif, nothing excludes it
            True,  # X-3  POS 67  right motif, plain ALT
            False,  # 8-4  POS 67  Motif '8' is on exclude_motifs_right
            True,  # Y-5  POS 67  GG alternate, kept: motifs_for_alt_gg is empty so
            #                     the allowlist narrowing never fires
            False,  # Z-6  POS 68  ALT 'CCGCC' is on exclude_alts_combined
            False,  # 6-7  POS 55  left motif '7' is on exclude_motifs_combined
        ]

    def test_the_input_rows_are_all_returned_the_stage_marks_it_does_not_filter(self):
        """``filter_final_dataframe`` gates on the column; this stage must not drop rows."""
        frame = _kestrel_frame()
        out = motif_correction_and_annotation(frame.copy(deep=True), _merged_motifs(), _shipped_config())

        assert len(out) == len(frame)
        assert out["Motifs"].tolist() == frame["Motifs"].tolist()
        assert "original_index" not in out.columns

    def test_only_passing_rows_carry_motif_motif_fasta_and_pos_fasta(self):
        """A failing row keeps its input columns and gets NA annotations, not garbage.

        ``POS_fasta`` equals ``POS`` for every passing row, including the ones at 67 that
        sit above ``position_threshold``: it is a verbatim copy, not a rebase (#203).
        """
        out = motif_correction_and_annotation(_kestrel_frame(), _merged_motifs(), _shipped_config())

        assert out["Motif"].tolist()[:2] == ["2", "X"]
        assert out["Motif"].iloc[3] == "Y"
        assert out["Motif"].isna().tolist() == [False, False, True, False, True, True]
        assert out["Motif_fasta"].tolist()[:2] == ["1-2", "X-3"]
        assert out["Motif_fasta"].isna().tolist() == [False, False, True, False, True, True]
        assert out["POS_fasta"].tolist()[:2] == [54, 67]
        assert out["POS_fasta"].iloc[3] == 67
        assert out["POS_fasta"].isna().tolist() == [False, False, True, False, True, True]

    def test_the_right_motif_is_the_left_half_and_the_left_motif_is_the_right_half(self):
        """The naming is inverted on purpose and is the easiest thing to 'fix' wrongly.

        A row at POS >= 60 sits in the *right* half of the repeat unit, so the
        motif it belongs to is the one named before the dash ('X' of 'X-3'). A
        row below the threshold takes the half after the dash ('2' of '1-2').
        """
        out = motif_correction_and_annotation(_kestrel_frame(), _merged_motifs(), _shipped_config())

        assert out.loc[0, "Motifs"] == "1-2" and out.loc[0, "Motif"] == "2", "POS 54 takes the right half"
        assert out.loc[1, "Motifs"] == "X-3" and out.loc[1, "Motif"] == "X", "POS 67 takes the left half"

    def test_pos_and_pos_fasta_are_both_the_vcf_position_in_the_pair_record(self):
        """#203: ``Motif_fasta`` names a 120 bp *pair* record, so ``POS`` is already right.

        Every record of ``All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`` is two
        60 bp motifs concatenated and named ``<left>-<right>``; ``1-2`` is a record ID,
        ``2`` is not. ``Motif_fasta`` is a verbatim copy of the VCF ``#CHROM``, i.e. the
        pair name, so the coordinate beside it in ``output.bed`` lives in the 120 bp
        record's space and the raw ``POS`` is the correct value. Subtracting
        ``position_threshold`` would make every row at or above 60 wrong by exactly 60 bp.

        ``Motif`` - the half-motif name - is what is repeat-unit-relative, and nothing
        writes a coordinate beside it.
        """
        out = motif_correction_and_annotation(_kestrel_frame(), _merged_motifs(), _shipped_config())

        assert out["POS"].tolist() == [54, 67, 67, 67, 68, 55]
        assert out.loc[1, "POS_fasta"] == 67, "not 7: no rebase, by design"
        assert out.loc[0, "POS_fasta"] == 54

    def test_the_position_threshold_boundary_decides_which_half_of_the_id_is_the_motif(self):
        """POS == 60 exactly is a RIGHT motif, so it takes the half before the dash.

        Nothing else in this file puts a row on the boundary, and ``<`` -> ``<=``
        on the left slice is invisible without one: the row simply changes which
        motif it is reported against.
        """
        frame = _kestrel_frame().iloc[:2].copy()
        frame["Motifs"] = ["1-2", "3-4"]
        frame["POS"] = [59, 60]
        frame["ALT"] = ["CC", "CC"]

        out = motif_correction_and_annotation(frame, _merged_motifs(), _shipped_config())

        assert out["motif_filter_pass"].tolist() == [True, True]
        assert out["Motif"].tolist() == ["2", "3"], "POS 59 takes the right half, POS 60 the left"

    def test_an_empty_frame_still_gains_the_gate_column(self):
        out = motif_correction_and_annotation(pd.DataFrame(), _merged_motifs(), _shipped_config())

        for column in ("motif_filter_pass", "Motif_fasta", "POS_fasta", "Motif"):
            assert column in out.columns

    def test_a_row_that_is_not_a_valid_frameshift_cannot_pass(self):
        """``motif_filter_pass`` is ANDed with ``is_valid_frameshift`` when present."""
        frame = _kestrel_frame(is_valid_frameshift=[True, False, True, pd.NA, True, True])

        out = motif_correction_and_annotation(frame, _merged_motifs(), _shipped_config())

        assert out["motif_filter_pass"].tolist() == [True, False, False, False, False, False]


# ---------------------------------------------------------------------------
# split_left_right
# ---------------------------------------------------------------------------


class TestSplitLeftRight:
    """The positional split: everything below the threshold is a left motif."""

    def _frame(self, positions):
        return pd.DataFrame({"Motifs": [f"L{i}-R{i}" for i in range(len(positions))], "POS": positions})

    def test_the_threshold_is_exclusive_below_and_inclusive_at_and_above(self):
        """POS == threshold is a RIGHT motif. ``<`` -> ``<=`` moves the boundary silently."""
        left, right = split_left_right(self._frame([59, 60, 61]), 60)

        assert left["POS"].tolist() == [59]
        assert right["POS"].tolist() == [60, 61]

    def test_every_row_lands_in_exactly_one_half(self):
        """``<`` -> ``>=`` on the left slice would put rows in both halves."""
        frame = self._frame([10, 59, 60, 200])
        left, right = split_left_right(frame, 60)

        assert len(left) + len(right) == len(frame)
        assert set(left["POS"]).isdisjoint(set(right["POS"]))

    def test_the_motif_id_is_split_on_the_dash_into_left_and_right_halves(self):
        left, right = split_left_right(pd.DataFrame({"Motifs": ["A-B"], "POS": [10]}), 60)

        assert left["Motif_left"].tolist() == ["A"]
        assert left["Motif_right"].tolist() == ["B"]
        assert right.empty

    def test_a_non_numeric_position_becomes_minus_one_and_lands_left(self):
        """``to_numeric(errors='coerce').fillna(-1)`` is a real branch, not defensive noise."""
        left, right = split_left_right(pd.DataFrame({"Motifs": ["A-B"], "POS": ["not-a-number"]}), 60)

        assert left["POS"].tolist() == [-1]
        assert right.empty

    def test_a_string_position_is_coerced_to_int_before_comparison(self):
        """Kestrel VCF positions arrive as strings; '9' must not compare as text."""
        left, right = split_left_right(pd.DataFrame({"Motifs": ["A-B", "C-D"], "POS": ["9", "600"]}), 60)

        assert left["POS"].tolist() == [9]
        assert right["POS"].tolist() == [600]

    def test_the_caller_s_frame_is_not_mutated(self):
        frame = pd.DataFrame({"Motifs": ["A-B"], "POS": ["67"]})

        split_left_right(frame, 60)

        assert frame["POS"].tolist() == ["67"]
        assert "Motif_left" not in frame.columns


# ---------------------------------------------------------------------------
# has_gg_alternate / apply_right_motif_exclusions / apply_gg_alt_rule
# ---------------------------------------------------------------------------


def _right_frame(motifs, alts):
    return pd.DataFrame({"Motif": motifs, "ALT": alts, "POS": [67] * len(motifs)})


class TestHasGgAlternate:
    """The gate on the whole legacy right-motif branch."""

    def test_a_gg_alternate_anywhere_in_the_frame_opens_the_branch(self):
        assert has_gg_alternate(_right_frame(["X", "Y"], ["GC", "GG"]), "GG") is True

    def test_no_gg_alternate_leaves_the_branch_shut(self):
        assert has_gg_alternate(_right_frame(["X", "Y"], ["GC", "CT"]), "GG") is False

    def test_the_match_is_word_bounded_so_gg_inside_a_longer_alt_does_not_count(self):
        """``\\bGG\\b`` - 'CGGCG' contains 'GG' but is a different allele."""
        assert has_gg_alternate(_right_frame(["X"], ["CGGCG"]), "GG") is False

    def test_the_alt_label_comes_from_config_and_is_not_hardcoded(self):
        assert has_gg_alternate(_right_frame(["X"], ["CT"]), "CT") is True
        assert has_gg_alternate(_right_frame(["X"], ["CT"]), "GG") is False


class TestApplyRightMotifExclusions:
    """Conserved motifs are dropped: a variant called there is likely an artifact."""

    def test_listed_motifs_are_removed_and_unlisted_ones_survive(self):
        out = apply_right_motif_exclusions(_right_frame(["X", "8", "9", "Y"], ["GG"] * 4), ["8", "9", "7", "6p", "6"])

        assert out["Motif"].tolist() == ["X", "Y"]

    def test_an_empty_exclusion_list_removes_nothing(self):
        """``~isin([])`` is all-True. ``not`` deleted here would remove everything."""
        out = apply_right_motif_exclusions(_right_frame(["X", "8"], ["GG", "GG"]), [])

        assert out["Motif"].tolist() == ["X", "8"]

    def test_excluding_every_motif_yields_an_empty_frame_not_an_error(self):
        out = apply_right_motif_exclusions(_right_frame(["8", "9"], ["GG", "GG"]), ["8", "9"])

        assert out.empty
        assert out.columns.tolist() == ["Motif", "ALT", "POS"]

    def test_the_caller_s_frame_is_not_mutated(self):
        frame = _right_frame(["X", "8"], ["GG", "GG"])

        apply_right_motif_exclusions(frame, ["8"])

        assert frame["Motif"].tolist() == ["X", "8"]


class TestApplyGgAltRule:
    """The GG allowlist narrowing - inert on the shipped config, armed by #186."""

    def test_an_empty_allowlist_leaves_the_frame_untouched(self):
        """This is the shipped config. The `.any()` guard is what makes it inert.

        Without the guard an empty allowlist would delete every row, which is
        the #186 failure mode the config guard in ``motif_processing`` refuses.
        """
        frame = _right_frame(["X", "Y"], ["GG", "GC"])

        out = apply_gg_alt_rule(frame, [])

        assert out["Motif"].tolist() == ["X", "Y"]

    def test_a_populated_allowlist_narrows_to_the_allowed_motifs(self):
        out = apply_gg_alt_rule(_right_frame(["X", "Y", "Z"], ["GG", "GC", "CT"]), ["X"])

        assert out["Motif"].tolist() == ["X"]

    def test_an_allowlist_that_matches_nothing_leaves_the_frame_untouched(self):
        """``.any()`` false -> no narrowing. ``in`` -> ``not in`` would invert this."""
        out = apply_gg_alt_rule(_right_frame(["X", "Y"], ["GG", "GC"]), ["NOTHERE"])

        assert out["Motif"].tolist() == ["X", "Y"]

    def test_the_caller_s_frame_is_not_mutated(self):
        frame = _right_frame(["X", "Y"], ["GG", "GC"])

        apply_gg_alt_rule(frame, ["X"])

        assert frame["Motif"].tolist() == ["X", "Y"]


# ---------------------------------------------------------------------------
# apply_combined_exclusions
# ---------------------------------------------------------------------------


class TestApplyCombinedExclusions:
    """The last gate, applied to left and right motifs together."""

    def _combined(self):
        return pd.DataFrame(
            {
                "Motif": ["X", "Y", "6", "7"],
                "ALT": ["GG", "CCGCC", "GC", "CT"],
            }
        )

    def test_listed_alts_are_removed(self):
        out = apply_combined_exclusions(self._combined(), ["CCGCC", "CGGCG", "CGGCC"], [])

        assert out["ALT"].tolist() == ["GG", "GC", "CT"]

    def test_listed_motifs_are_removed(self):
        out = apply_combined_exclusions(self._combined(), [], ["6", "6p", "7"])

        assert out["Motif"].tolist() == ["X", "Y"]

    def test_both_lists_apply_and_a_row_can_be_caught_by_either(self):
        out = apply_combined_exclusions(self._combined(), ["CCGCC", "CGGCG", "CGGCC"], ["6", "6p", "7"])

        assert out["Motif"].tolist() == ["X"]

    def test_two_empty_lists_remove_nothing(self):
        """``not`` deleted from either ``~`` would invert the whole gate to a keeplist."""
        out = apply_combined_exclusions(self._combined(), [], [])

        assert len(out) == 4

    def test_the_alt_match_is_exact_not_a_substring(self):
        frame = pd.DataFrame({"Motif": ["X"], "ALT": ["CCGCCA"]})

        out = apply_combined_exclusions(frame, ["CCGCC"], [])

        assert len(out) == 1, "'CCGCCA' is a different allele from 'CCGCC'"

    def test_the_caller_s_frame_is_not_mutated(self):
        frame = self._combined()

        apply_combined_exclusions(frame, ["CCGCC"], ["6"])

        assert len(frame) == 4


# ---------------------------------------------------------------------------
# Defect W8 (#179 audit): the whole-column .max() on the dash count
# ---------------------------------------------------------------------------


def _dash_frame(motifs):
    """One row per motif ID, all otherwise identical and all otherwise passing."""
    n = len(motifs)
    return pd.DataFrame(
        {
            "Motifs": motifs,
            "Variant": ["Insertion"] * n,
            "POS": [67] * n,
            "REF": ["C"] * n,
            "ALT": ["CC"] * n,
            "Estimated_Depth_AlternateVariant": [50] * n,
            "Estimated_Depth_Variant_ActiveRegion": [5000] * n,
            "Depth_Score": [0.01] * n,
            "Confidence": ["High_Precision"] * n,
        }
    )


class TestMalformedMotifIdsAreContainedToTheirOwnRow:
    """Defect W8 from the #179 audit, fixed here.

    ``motif_processing.py`` used to gate the split on an aggregate::

        if "Motifs" not in working_df.columns or working_df["Motifs"].str.count("-").max() != 1:

    ``.max()`` reduces the whole column, so the verdict for every row was decided
    by whichever row had the most dashes. That failed in both directions:

    * one row with two dashes drove ``max()`` to 2, the branch returned an empty
      frame, and **every** row in the sample got ``motif_filter_pass=False`` -
      the report says ``Negative`` and the pipeline exits 0;
    * a row with *no* dash rode through unnoticed whenever some other row had
      exactly one, because ``max()`` was still 1. It was split into
      ``Motif='BROKEN', Motif_right=None``, merged to nothing, and reported as
      a passing call against a motif that is not in the reference.

    Both are now decided per row. The disposition of the malformed row follows
    what this stage already does with every row it rejects (see the oracle
    above): the row stays in the returned frame, ``motif_filter_pass`` is False
    and the annotation columns are NA. This stage marks, it does not filter -
    ``filter_final_dataframe`` in ``kestrel_genotyping.py`` is what drops rows,
    and it needs the gate column present on every input row.
    """

    def test_an_unsplittable_motif_id_fails_only_its_own_row(self):
        """SPECIFICATION. One malformed motif ID may not suppress the sample."""
        out = motif_correction_and_annotation(_dash_frame(["X-5", "X-5-6"]), _merged_motifs(), _shipped_config())

        assert out["motif_filter_pass"].tolist() == [True, False]
        assert out["Motif"].iloc[0] == "X"
        assert pd.isna(out["Motif"].iloc[1]), "the malformed row is failed, not annotated with a guess"
        assert len(out) == 2, "and it is still returned - this stage marks, it does not filter"

    def test_a_motif_id_with_no_dash_fails_rather_than_riding_through(self):
        """SPECIFICATION. The under-suppressing half of W8.

        Under ``.max()`` this row passed whenever any sibling row had exactly one
        dash, and was annotated ``Motif='BROKEN'``.
        """
        out = motif_correction_and_annotation(_dash_frame(["X-5", "BROKEN"]), _merged_motifs(), _shipped_config())

        assert out["motif_filter_pass"].tolist() == [True, False]
        assert pd.isna(out["Motif"].iloc[1])

    def test_a_frame_where_no_row_can_be_split_still_fails_every_row(self):
        """Unchanged behaviour: nothing splittable means nothing passes."""
        out = motif_correction_and_annotation(_dash_frame(["BROKEN", "A-B-2"]), _merged_motifs(), _shipped_config())

        assert out["motif_filter_pass"].tolist() == [False, False]
        assert len(out) == 2

    def test_a_frame_where_every_row_splits_is_untouched_by_the_per_row_guard(self):
        """The no-malformed-rows path - the only one the 58 real samples exercise."""
        out = motif_correction_and_annotation(_dash_frame(["X-5", "Y-6"]), _merged_motifs(), _shipped_config())

        assert out["motif_filter_pass"].tolist() == [True, True]

    def test_a_frame_with_no_motifs_column_at_all_fails_every_row(self):
        """Unchanged behaviour: without the column there is nothing to split on."""
        out = motif_correction_and_annotation(
            _dash_frame(["X-5"]).drop(columns=["Motifs"]), _merged_motifs(), _shipped_config()
        )

        assert out["motif_filter_pass"].tolist() == [False]
        assert len(out) == 1

    def test_the_dropped_rows_are_named_in_the_log(self, caplog):
        """A silently suppressed row is the failure mode; the log has to name it."""
        caplog.set_level(logging.WARNING, logger="vntyper.scripts.motif_processing")

        motif_correction_and_annotation(_dash_frame(["X-5", "X-5-6"]), _merged_motifs(), _shipped_config())

        assert "X-5-6" in caplog.text
        assert "1" in caplog.text, "the count of rows dropped"
