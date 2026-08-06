# tests/unit/test_advntr_frameshift_filter.py

"""
Characterisation of the adVNTR frameshift filter
(``advntr_processing_del`` / ``advntr_processing_ins``).

**These tests pin current behaviour, including behaviour that is wrong.** They are
deliberately not a fix. The filter decides which adVNTR calls reach
``output_adVNTR_result.tsv``, and that file is a reported genotype: it drives the adVNTR
column of the HTML report, the cohort summary, and the Kestrel/adVNTR cross-match.
Changing which rows survive changes what a clinician reads, so the repair belongs to a
human decision, not to a test-coverage pass.

What the code does
------------------
Both functions compute ``frame = abs(Insertion_len - Deletion_length)`` and then keep a
row only if ``frame`` lands in a per-function arithmetic series::

    del_frame = arange(100) * 3 + 2   ->  {2, 5, 8, ...}   (== 2 mod 3)
    ins_frame = arange(100) * 3 + 1   ->  {1, 4, 7, ...}   (== 1 mod 3)

Between them those two series cover every net length change that is *not* a multiple of
three -- that is, every frameshift. But the two halves are gated by different guards:
``advntr_processing_del`` additionally requires ``Deletion_length >= 1`` and
``advntr_processing_ins`` requires ``Insertion_len >= 1``. A *pure* deletion has
``Insertion_len == 0``, so it can never satisfy the insertion guard, and a pure insertion
has ``Deletion_length == 0``, so it can never satisfy the deletion guard.

The consequence, which the parametrised cases below pin exactly:

* a pure deletion whose length is ``1 mod 3`` (1, 4, 7, ... bp) is a frameshift and is
  **silently dropped**;
* a pure insertion whose length is ``2 mod 3`` (2, 5, 8, ... bp) is a frameshift and is
  **silently dropped**.

The 1 bp deletion is the headline case: it is the single most common frameshift class and
it never reaches the report. The pathogenic MUC1 dupC is a 1 bp *insertion*, which does
survive -- which is why this has gone unnoticed.
"""

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr

pytestmark = pytest.mark.unit


def variant_frame(state: str) -> pd.DataFrame:
    """A single-row adVNTR frame in the shape ``process_advntr_output`` hands over."""
    return pd.DataFrame(
        {
            "VID": [25561],
            "State": [state],
            "NumberOfSupportingReads": [11],
            "MeanCoverage": [153.98],
            "Pvalue": [0.0001],
        }
    )


def pure_deletion(length: int) -> str:
    """An adVNTR state naming ``length`` consecutive single-base deletions."""
    return "&".join(f"D{8 + offset}_2" for offset in range(length))


def pure_insertion(length: int) -> str:
    """An adVNTR state naming a single insertion of ``length`` bases."""
    return f"I10_2_A_LEN{length}"


def surviving_rows(state: str) -> int:
    """How many rows survive the full del+ins filter, exactly as the parser combines them."""
    frame = variant_frame(state)
    kept_del = advntr.advntr_processing_del(frame.copy())
    kept_ins = advntr.advntr_processing_ins(frame.copy())
    return len(kept_del) + len(kept_ins)


# ---------------------------------------------------------------------------
# Column derivation -- what the two functions add before they filter
# ---------------------------------------------------------------------------


class TestDerivedColumns:
    def test_state_is_renamed_to_variant(self):
        result = advntr.advntr_processing_ins(variant_frame(pure_insertion(1)))

        assert "Variant" in result.columns
        assert "State" not in result.columns
        assert list(result["Variant"]) == ["I10_2_A_LEN1"]

    def test_deletion_and_insertion_counts_come_from_the_letters_in_the_state(self):
        # net -1, so this row survives the insertion half and its columns are observable.
        result = advntr.advntr_processing_ins(variant_frame("D8_2&D9_2&I9_2_A_LEN1"))

        assert list(result["Deletion_length"]) == [2]
        assert list(result["Insertion_length"]) == [1]

    def test_a_state_without_a_len_token_yields_insertion_len_zero(self):
        result = advntr.advntr_processing_del(variant_frame(pure_deletion(2)))

        assert list(result["Insertion_len"]) == [0]

    def test_frame_is_the_absolute_difference_and_is_carried_as_a_string(self):
        result = advntr.advntr_processing_del(variant_frame("D8_2&D9_2&D10_2&I10_2_A_LEN1"))

        # 1 inserted base, 3 deleted bases -> net -2.
        assert list(result["frame"]) == ["2"]
        assert result["frame"].dtype == object


# ---------------------------------------------------------------------------
# E4 -- the dropped frameshifts
# ---------------------------------------------------------------------------


class TestPureDeletionsAreFilteredByLength:
    """``length % 3`` fully determines the outcome; ``1 mod 3`` is dropped in error."""

    @pytest.mark.parametrize(
        ("length", "shifts_the_frame", "kept_today"),
        [
            (1, True, False),  # BUG: a 1 bp deletion is a frameshift and is dropped
            (2, True, True),
            (3, False, False),  # correct: in frame
            (4, True, False),  # BUG
            (5, True, True),
            (6, False, False),  # correct: in frame
            (7, True, False),  # BUG
            (8, True, True),
            (9, False, False),  # correct: in frame
        ],
    )
    def test_a_pure_deletion_survives_only_when_its_length_is_two_mod_three(self, length, shifts_the_frame, kept_today):
        """
        Specification (#182, decided 2026-08-06). @hassansaei: "keep the same 3n+1 / 3n+2
        rule for adVNTR as for Kestrel (#181). This is intentional shared convention, not
        something to relax independently." For a pure deletion, ``Deletion_length`` is
        always >= 1, so ``advntr_processing_del``'s guard never decides the outcome here --
        only ``del_frame`` membership (``length % 3 == 2``) does, which this table pins as
        decided behaviour. It does not rule on the separate guard-interaction gap visible
        in the 1/4/7 rows -- see
        ``test_a_frameshifting_deletion_of_one_mod_three_bases_is_lost``, which stays a
        characterisation.
        """
        assert (length % 3 != 0) is shifts_the_frame, "test table disagrees with arithmetic"

        assert surviving_rows(pure_deletion(length)) == (1 if kept_today else 0)

    @pytest.mark.parametrize("length", [1, 4, 7])
    def test_a_frameshifting_deletion_of_one_mod_three_bases_is_lost(self, length):
        """
        The defect, stated on its own. These deletions shift the reading frame and are
        therefore candidate pathogenic calls, yet neither half of the filter accepts them:
        ``advntr_processing_del`` rejects the frame, and ``advntr_processing_ins`` rejects
        the row for having no inserted bases.
        """
        state = pure_deletion(length)

        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0, "rejected on frame"
        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 0, "rejected on Insertion_len >= 1"


class TestPureInsertionsAreFilteredByLength:
    """The mirror image: ``2 mod 3`` insertions are dropped in error."""

    @pytest.mark.parametrize(
        ("length", "shifts_the_frame", "kept_today"),
        [
            (1, True, True),  # the pathogenic MUC1 dupC lives here, which is why this went unseen
            (2, True, False),  # BUG
            (3, False, False),  # correct: in frame
            (4, True, True),
            (5, True, False),  # BUG
            (6, False, False),  # correct: in frame
            (7, True, True),
            (8, True, False),  # BUG
            (9, False, False),  # correct: in frame
        ],
    )
    def test_a_pure_insertion_survives_only_when_its_length_is_one_mod_three(
        self, length, shifts_the_frame, kept_today
    ):
        """
        Specification (#182, decided 2026-08-06). @hassansaei: "keep the same 3n+1 / 3n+2
        rule for adVNTR as for Kestrel (#181). This is intentional shared convention, not
        something to relax independently." For a pure insertion, ``Insertion_len`` is
        always >= 1, so ``advntr_processing_ins``'s guard never decides the outcome here --
        only ``ins_frame`` membership (``length % 3 == 1``) does, which this table pins as
        decided behaviour. It does not rule on the separate guard-interaction gap visible
        in the 2/5/8 rows -- see
        ``test_a_frameshifting_insertion_of_two_mod_three_bases_is_lost``, which stays a
        characterisation.
        """
        assert (length % 3 != 0) is shifts_the_frame, "test table disagrees with arithmetic"

        assert surviving_rows(pure_insertion(length)) == (1 if kept_today else 0)

    @pytest.mark.parametrize("length", [2, 5, 8])
    def test_a_frameshifting_insertion_of_two_mod_three_bases_is_lost(self, length):
        state = pure_insertion(length)

        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 0, "rejected on frame"
        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0, "rejected on Deletion_length >= 1"

    def test_the_canonical_integration_fixture_variant_still_survives(self):
        """``I22_2_G_LEN1`` is what the shipped adVNTR fixture reports; it must keep passing."""
        assert surviving_rows("I22_2_G_LEN1") == 1


class TestMixedIndels:
    """A state naming both deletions and insertions can satisfy either guard."""

    def test_a_deletion_with_a_compensating_insertion_is_judged_on_the_net_change(self):
        # 2 deleted bases, 1 inserted base -> net -1, frame 1, which is the insertion series.
        state = "D8_2&D9_2&I9_2_A_LEN1"

        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 1
        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0

    def test_a_net_in_frame_indel_is_dropped_by_both_halves(self):
        # 3 deleted bases, 6 inserted -> net +3, in frame.
        state = "D8_2&D9_2&D10_2&I10_2_A_LEN6"

        assert surviving_rows(state) == 0


# ---------------------------------------------------------------------------
# H1 -- which module global the filter actually reads
# ---------------------------------------------------------------------------


class TestSettingsComeFromTheDerivedGlobal:
    """
    ``advntr_genotyping`` builds two globals at import: ``advntr_config`` (the whole JSON)
    and ``advntr_settings`` (its ``advntr_settings`` sub-dict). The frameshift filter reads
    the **derived** one. Patching ``advntr_config`` after import is a silent no-op, so a
    test that patches it proves nothing.
    """

    def test_patching_the_derived_settings_changes_the_accepted_frames(self, monkeypatch):
        # multiplier 1 makes del_frame {2, 3, 4, ...}, so a 3 bp deletion starts passing.
        monkeypatch.setattr(advntr, "advntr_settings", {"max_frameshift": 100, "frameshift_multiplier": 1})

        assert len(advntr.advntr_processing_del(variant_frame(pure_deletion(3)))) == 1

    def test_patching_the_raw_config_global_has_no_effect(self, monkeypatch):
        monkeypatch.setattr(
            advntr,
            "advntr_config",
            {"advntr_settings": {"max_frameshift": 100, "frameshift_multiplier": 1}},
        )

        assert len(advntr.advntr_processing_del(variant_frame(pure_deletion(3)))) == 0

    def test_max_frameshift_bounds_the_accepted_series(self, monkeypatch):
        monkeypatch.setattr(advntr, "advntr_settings", {"max_frameshift": 1, "frameshift_multiplier": 3})

        # With max_frameshift 1 the only accepted deletion frame is 2.
        assert len(advntr.advntr_processing_del(variant_frame(pure_deletion(2)))) == 1
        assert len(advntr.advntr_processing_del(variant_frame(pure_deletion(5)))) == 0
