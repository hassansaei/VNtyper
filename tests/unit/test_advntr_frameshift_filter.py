# tests/unit/test_advntr_frameshift_filter.py

"""
The adVNTR pathogenic-frame filter
(``advntr_processing_del`` / ``advntr_processing_ins``).

The filter decides which adVNTR calls reach ``output_adVNTR_result.tsv``, and that file is
a reported genotype: it drives the adVNTR column of the HTML report, the cohort summary,
and the Kestrel/adVNTR cross-match. Which rows survive is therefore what a clinician
reads, and the rule below is a domain decision (#182, @hassansaei), not an implementation
detail to be relaxed for tidiness.

Two different questions, and VNtyper answers the narrower one
-------------------------------------------------------------
Let ``Delta = Insertion_len - Deletion_length`` be the **signed** net change in bases
(``Net_indel_length`` in the code). Two distinct predicates can be built on it:

1. *"Does this change the reading frame at all?"* -- ``Delta % 3 != 0``. True for a net
   change of 1, 2, 4, 5, 7, 8, ... bases in **either** direction.
2. *"Does this enter the pathogenic ADTKD-MUC1 frame?"* -- ``Delta % 3 == 1``, i.e.
   ``Delta = +1 (mod 3)``. That is the frame yielding the toxic MUC1-fs neo-protein,
   classically exemplified by dupC (a 1 bp insertion, ``Delta = +1``).

**VNtyper filters on (2), the narrower one.** Written out over the sign:

===========================  ============  ==================================
net change                   ``Delta % 3``  verdict
===========================  ============  ==================================
insertion of 3n+1 bases      1              **pathogenic frame** -- reported
deletion of 3n+2 bases       1              **pathogenic frame** -- reported
insertion of 3n+2 bases      2              frameshift, other frame -- dropped
deletion of 3n+1 bases       2              frameshift, other frame -- dropped
net change of 3n bases       0              in frame -- dropped
===========================  ============  ==================================

A ``Delta = 2 (mod 3)`` state really is a frameshift. It is dropped on purpose: that frame
has not been established as pathogenic in ADTKD-MUC1 patients, and reporting it would
broaden a clinical filter on no evidence. So the 1/4/7 bp pure deletions and the 2/5/8 bp
pure insertions below are **not** missing pathogenic calls -- they are correctly-rejected
non-pathogenic frameshifts. (An earlier revision of this file labelled them "BUG" and
called the 1 bp deletion "the headline case"; that reading conflated question (1) with
question (2). It is corrected here.)

How the two functions implement it
----------------------------------
``advntr_processing_del`` keeps ``Delta < 0`` with ``|Delta|`` in the ``3n+2`` series;
``advntr_processing_ins`` keeps ``Delta > 0`` with ``Delta`` in the ``3n+1`` series. Both
series come from ``advntr_config.json``'s ``frameshift_multiplier`` and ``max_frameshift``
(see :class:`TestSettingsComeFromTheDerivedGlobal`). Their union is exactly
``Delta % 3 == 1``, because Python's ``%`` returns a non-negative residue: ``-2 % 3 == 1``,
``-5 % 3 == 1``, while ``-1 % 3 == 2``.

**The sign test is the repair this file also guards.** Both arms used to test
``abs(Delta)`` against their series while guarding only on "names at least one deletion" /
"names at least one insertion". A *mixed* state satisfies both guards, so either arm could
admit it on the magnitude alone -- and the arm that did so was accepting the opposite,
non-pathogenic frame. See :class:`TestMixedIndelsAreJudgedOnTheSignedNetChange`.
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


def mixed_state(inserted: int, deleted: int) -> str:
    """A state naming ``deleted`` single-base deletions and one ``inserted``-bp insertion.

    ``Deletion_length`` is ``Variant.str.count("D")``, so ``deleted`` deletion parts give
    ``Deletion_length == deleted``; ``Insertion_len`` is the sum over ``LEN`` tokens, so
    the single insertion part gives ``Insertion_len == inserted``. Hence
    ``Net_indel_length == inserted - deleted`` exactly.

    Args:
        inserted (int): Inserted bases, >= 1.
        deleted (int): Deleted bases, >= 1.

    Returns:
        str: e.g. ``mixed_state(3, 1) == "I9_2_A_LEN3&D50_2"``.
    """
    parts = [f"I9_2_A_LEN{inserted}"] + [f"D{50 + offset}_2" for offset in range(deleted)]
    return "&".join(parts)


def surviving_rows(state: str) -> int:
    """How many rows survive the full del+ins filter, exactly as the parser combines them."""
    frame = variant_frame(state)
    kept_del = advntr.advntr_processing_del(frame.copy())
    kept_ins = advntr.advntr_processing_ins(frame.copy())
    return len(kept_del) + len(kept_ins)


def net_indel_length(state: str) -> int:
    """``Insertion_len - Deletion_length`` for ``state``, from the production derivation."""
    return int(advntr.derive_indel_columns(variant_frame(state))["Net_indel_length"].iloc[0])


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
        # 4 inserted, 3 deleted -> Delta = +1, the pathogenic frame, so this row survives
        # the insertion arm and its derived columns are observable on the result.
        result = advntr.advntr_processing_ins(variant_frame("D8_2&D9_2&D10_2&I10_2_A_LEN4"))

        assert list(result["Deletion_length"]) == [3]
        assert list(result["Insertion_length"]) == [1]

    def test_a_state_without_a_len_token_yields_insertion_len_zero(self):
        result = advntr.advntr_processing_del(variant_frame(pure_deletion(2)))

        assert list(result["Insertion_len"]) == [0]

    def test_net_indel_length_is_signed_and_frame_is_its_magnitude(self):
        """``Net_indel_length`` keeps the direction; ``frame`` is the magnitude only.

        This split is the whole repair: ``frame`` is what the configurable ``3n+1`` /
        ``3n+2`` series are matched against, and it cannot carry a sign, so the sign has
        to live in its own column for each arm to test it.
        """
        # 1 inserted base, 3 deleted bases -> Delta = -2.
        result = advntr.advntr_processing_del(variant_frame("D8_2&D9_2&D10_2&I10_2_A_LEN1"))

        assert list(result["Net_indel_length"]) == [-2]
        assert list(result["frame"]) == ["2"]
        assert result["frame"].dtype == object

    def test_derive_indel_columns_drops_no_rows(self):
        """The derivation is total; only the two arms filter."""
        frame = pd.DataFrame(
            {
                "VID": [25561, 25561, 25561],
                "State": [pure_insertion(1), pure_deletion(3), mixed_state(3, 1)],
                "NumberOfSupportingReads": [11, 11, 11],
                "MeanCoverage": [153.98, 153.98, 153.98],
                "Pvalue": [0.0001, 0.0001, 0.0001],
            }
        )

        derived = advntr.derive_indel_columns(frame)

        assert len(derived) == 3
        assert list(derived["Net_indel_length"]) == [1, -3, 2]
        assert list(derived["frame"]) == ["1", "3", "2"]


# ---------------------------------------------------------------------------
# Pure states -- one direction only, so the sign is never in doubt
# ---------------------------------------------------------------------------


class TestPureDeletionsAreFilteredByLength:
    """``length % 3`` fully determines the outcome for a pure deletion."""

    @pytest.mark.parametrize(
        ("length", "shifts_the_frame", "pathogenic_frame"),
        [
            (1, True, False),  # Delta = -1, -1 % 3 == 2: a frameshift, the other frame
            (2, True, True),  # Delta = -2, -2 % 3 == 1: the pathogenic frame
            (3, False, False),  # Delta = -3: in frame
            (4, True, False),  # Delta = -4 == 2 mod 3
            (5, True, True),  # Delta = -5 == 1 mod 3
            (6, False, False),  # Delta = -6: in frame
            (7, True, False),  # Delta = -7 == 2 mod 3
            (8, True, True),  # Delta = -8 == 1 mod 3
            (9, False, False),  # Delta = -9: in frame
        ],
    )
    def test_a_pure_deletion_survives_only_when_its_length_is_two_mod_three(
        self, length, shifts_the_frame, pathogenic_frame
    ):
        """
        Specification (#182, decided 2026-08-06). @hassansaei: "keep the same 3n+1 / 3n+2
        rule for adVNTR as for Kestrel (#181). This is intentional shared convention, not
        something to relax independently."

        ``shifts_the_frame`` and ``pathogenic_frame`` are the two predicates the module
        docstring separates, and the third column -- not the second -- is what the filter
        implements. The 1/4/7 bp rows are frameshifts that VNtyper deliberately does not
        report; they are not lost pathogenic calls.
        """
        assert (length % 3 != 0) is shifts_the_frame, "test table disagrees with arithmetic"
        assert (-length % 3 == 1) is pathogenic_frame, "test table disagrees with the signed rule"

        assert surviving_rows(pure_deletion(length)) == (1 if pathogenic_frame else 0)

    @pytest.mark.parametrize("length", [1, 4, 7])
    def test_a_one_mod_three_deletion_is_a_frameshift_but_not_the_pathogenic_frame(self, length):
        """
        Stated on its own because it is the case most often misread as a defect.

        A 1, 4 or 7 bp deletion gives ``Delta = -1, -4, -7``, every one of which is
        ``2 (mod 3)`` -- the reading frame moves, but into the frame VNtyper does not
        report. Neither arm accepts it, and each rejects it for the right reason: the
        deletion arm because the magnitude is not in its ``3n+2`` series, the insertion
        arm because the net change is not an insertion at all.
        """
        state = pure_deletion(length)

        assert net_indel_length(state) == -length
        assert -length % 3 == 2, "these lengths are the non-pathogenic frame"
        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0, "magnitude not in the 3n+2 series"
        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 0, "net change is not an insertion"


class TestPureInsertionsAreFilteredByLength:
    """The mirror image: for a pure insertion, ``1 mod 3`` is the pathogenic frame."""

    @pytest.mark.parametrize(
        ("length", "shifts_the_frame", "pathogenic_frame"),
        [
            (1, True, True),  # Delta = +1: the pathogenic frame -- this is where dupC lives
            (2, True, False),  # Delta = +2 == 2 mod 3: a frameshift, the other frame
            (3, False, False),  # in frame
            (4, True, True),
            (5, True, False),
            (6, False, False),  # in frame
            (7, True, True),
            (8, True, False),
            (9, False, False),  # in frame
        ],
    )
    def test_a_pure_insertion_survives_only_when_its_length_is_one_mod_three(
        self, length, shifts_the_frame, pathogenic_frame
    ):
        """
        Specification (#182, decided 2026-08-06). @hassansaei: "keep the same 3n+1 / 3n+2
        rule for adVNTR as for Kestrel (#181). This is intentional shared convention, not
        something to relax independently."

        The pathogenic MUC1 dupC is the ``length == 1`` row: a 1 bp insertion,
        ``Delta = +1``. The 2/5/8 bp rows are frameshifts in the other frame and are
        dropped on purpose, exactly as the 1/4/7 bp deletions are.
        """
        assert (length % 3 != 0) is shifts_the_frame, "test table disagrees with arithmetic"
        assert (length % 3 == 1) is pathogenic_frame, "test table disagrees with the signed rule"

        assert surviving_rows(pure_insertion(length)) == (1 if pathogenic_frame else 0)

    @pytest.mark.parametrize("length", [2, 5, 8])
    def test_a_two_mod_three_insertion_is_a_frameshift_but_not_the_pathogenic_frame(self, length):
        state = pure_insertion(length)

        assert net_indel_length(state) == length
        assert length % 3 == 2, "these lengths are the non-pathogenic frame"
        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 0, "magnitude not in the 3n+1 series"
        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0, "net change is not a deletion"

    def test_the_canonical_integration_fixture_variant_still_survives(self):
        """``I22_2_G_LEN1`` is what the shipped adVNTR fixture reports; it must keep passing."""
        assert surviving_rows("I22_2_G_LEN1") == 1


# ---------------------------------------------------------------------------
# Mixed states -- where the sign of the net change is decided, and used to be lost
# ---------------------------------------------------------------------------


class TestMixedIndelsAreJudgedOnTheSignedNetChange:
    """A state naming both deletions and insertions is judged on ``Delta``, sign included.

    This is the class the repair exists for. Every previous formulation guarded only on
    presence -- ``Deletion_length >= 1`` for the deletion arm, ``Insertion_len >= 1`` for
    the insertion arm -- and a mixed state satisfies **both**, so a magnitude in either
    series was enough to be reported by the arm that happened to match it, whatever the
    real direction of the change was.
    """

    #: ``(inserted, deleted)``, the residue class of ``Delta``, and the verdict. Every
    #: residue appears with a net-insertion and a net-deletion sign, which is the whole
    #: point: only the sign distinguishes the two rows of each residue pair.
    #:
    #: ============  =====  ===========  =========================================
    #: state         Delta  ``Delta%3``  verdict
    #: ============  =====  ===========  =========================================
    #: (4, 1)         +3     0           in frame -- dropped
    #: (1, 4)         -3     0           in frame -- dropped
    #: (5, 1)         +4     1           net insertion of 3n+1 -- REPORTED (ins)
    #: (1, 3)         -2     1           net deletion of 3n+2 -- REPORTED (del)
    #: (3, 1)         +2     2           net insertion of 3n+2 -- dropped
    #: (1, 2)         -1     2           net deletion of 3n+1 -- dropped
    #: ============  =====  ===========  =========================================
    SIGNED_CASES = [
        pytest.param(4, 1, 3, 0, False, False, id="delta+3-in-frame"),
        pytest.param(1, 4, -3, 0, False, False, id="delta-3-in-frame"),
        pytest.param(5, 1, 4, 1, True, False, id="delta+4-pathogenic-via-ins"),
        pytest.param(1, 3, -2, 1, False, True, id="delta-2-pathogenic-via-del"),
        pytest.param(3, 1, 2, 2, False, False, id="delta+2-other-frame"),
        pytest.param(1, 2, -1, 2, False, False, id="delta-1-other-frame"),
    ]

    @pytest.mark.parametrize(("inserted", "deleted", "delta", "residue", "kept_by_ins", "kept_by_del"), SIGNED_CASES)
    def test_each_signed_residue_class_reaches_the_right_arm(
        self, inserted, deleted, delta, residue, kept_by_ins, kept_by_del
    ):
        """Both arms are asserted for every case, so an arm cannot claim a row it must not.

        The arithmetic is restated in the parameters and cross-checked against Python here,
        so a wrong expectation cannot be smuggled in as a plausible-looking id.
        """
        assert inserted - deleted == delta, "case table disagrees with its own arithmetic"
        assert delta % 3 == residue, "case table disagrees with Python's signed modulo"
        assert (residue == 1 and delta > 0) is kept_by_ins, "case table disagrees with the signed rule"
        assert (residue == 1 and delta < 0) is kept_by_del, "case table disagrees with the signed rule"

        state = mixed_state(inserted, deleted)
        assert net_indel_length(state) == delta

        assert len(advntr.advntr_processing_ins(variant_frame(state))) == (1 if kept_by_ins else 0)
        assert len(advntr.advntr_processing_del(variant_frame(state))) == (1 if kept_by_del else 0)
        assert surviving_rows(state) == (1 if (kept_by_ins or kept_by_del) else 0)

    def test_a_net_insertion_of_two_mod_three_is_not_admitted_by_the_deletion_arm(self):
        """Regression for the state this branch made newly reachable.

        ``I9_2_A_LEN3&D50_2``: 3 inserted bases against 1 deleted base, so ``Delta = +2``,
        which is ``2 (mod 3)`` -- a frameshift into the frame VNtyper does not report.

        On ``main`` the state never got this far: ``Insertion_len`` collapsed to 0 for any
        state with material after its first ``LEN``, giving ``Delta = -1`` and ``|Delta| =
        1``, which is in neither series. Commit ``2a267aa`` (#192) made ``Insertion_len``
        the sum of every ``LEN`` token, raising it to 3 -- and because the filter then
        tested ``|Delta| = 2`` against the deletion arm's ``3n+2`` series while guarding
        only on ``Deletion_length >= 1``, the state was **reported via the deletion arm**
        despite being a net insertion. The sign error is older; #192 is what made it reach
        a real state shape.
        """
        state = "I9_2_A_LEN3&D50_2"

        assert net_indel_length(state) == 2
        assert 2 % 3 == 2, "not the pathogenic frame"
        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0, "net insertion, so not a deletion call"
        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 0, "magnitude 2 is not in the 3n+1 series"
        assert surviving_rows(state) == 0

    def test_a_net_deletion_of_one_mod_three_is_not_admitted_by_the_insertion_arm(self):
        """The mirror case, and the one a previous test asserted backwards.

        ``D8_2&D9_2&I9_2_A_LEN1``: 1 inserted base against 2 deleted bases, so
        ``Delta = -1``, which is ``2 (mod 3)`` -- a frameshift, but not the pathogenic
        one. ``|Delta| = 1`` is in the insertion arm's ``3n+1`` series, and the arm's only
        other guard used to be ``Insertion_len >= 1``, which this state satisfies. It was
        therefore reported as an insertion call although its net change is a deletion.
        """
        state = "D8_2&D9_2&I9_2_A_LEN1"

        assert net_indel_length(state) == -1
        assert -1 % 3 == 2, "not the pathogenic frame"
        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 0, "net deletion, so not an insertion call"
        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0, "magnitude 1 is not in the 3n+2 series"
        assert surviving_rows(state) == 0

    def test_a_compensated_deletion_in_the_pathogenic_frame_is_still_reported(self):
        """The repair narrows the filter; it must not close it.

        ``I9_2_A_LEN3&D50_2&D51_2``: 3 inserted against 2 deleted, ``Delta = +1`` -- the
        pathogenic frame, reached by a compound state. It is reported, via the insertion
        arm, exactly as it was before.
        """
        state = "I9_2_A_LEN3&D50_2&D51_2"

        assert net_indel_length(state) == 1
        assert len(advntr.advntr_processing_ins(variant_frame(state))) == 1
        assert len(advntr.advntr_processing_del(variant_frame(state))) == 0

    def test_a_net_in_frame_indel_is_dropped_by_both_halves(self):
        # 3 deleted bases, 6 inserted -> Delta = +3, in frame.
        state = "D8_2&D9_2&D10_2&I10_2_A_LEN6"

        assert net_indel_length(state) == 3
        assert surviving_rows(state) == 0

    def test_no_state_is_ever_reported_by_both_arms_at_once(self):
        """The two arms are disjoint, so ``pd.concat`` cannot double-count a row.

        ``process_advntr_output`` concatenates the two halves without de-duplicating on
        the state, so an overlap would emit the same call twice. Disjointness follows from
        the sign tests being strict and mutually exclusive; this asserts it over the whole
        mixed grid rather than trusting the argument.
        """
        for inserted in range(1, 10):
            for deleted in range(1, 10):
                state = mixed_state(inserted, deleted)
                by_ins = len(advntr.advntr_processing_ins(variant_frame(state)))
                by_del = len(advntr.advntr_processing_del(variant_frame(state)))

                assert by_ins + by_del <= 1, f"{state} reported by both arms"
                assert (by_ins + by_del == 1) == ((inserted - deleted) % 3 == 1), state


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
        # multiplier 1 makes the deletion series {2, 3, 4, ...}, so a 3 bp deletion starts
        # passing. The sign test is unaffected -- Delta = -3 is still a net deletion.
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

        # With max_frameshift 1 the only accepted deletion magnitude is 2.
        assert len(advntr.advntr_processing_del(variant_frame(pure_deletion(2)))) == 1
        assert len(advntr.advntr_processing_del(variant_frame(pure_deletion(5)))) == 0

    def test_the_settings_still_govern_the_series_after_the_sign_test_was_added(self, monkeypatch):
        """The sign test must not have hardcoded what stayed configurable.

        A multiplier of 1 makes the insertion series ``{1, 2, 3, ...}`` and the deletion
        series ``{2, 3, 4, ...}``, so magnitudes the shipped ``3n+1``/``3n+2`` series
        reject start passing -- for both signs, and each only through the arm matching its
        sign. That is the property to hold: the sign test is unconditional, the series is
        configuration.
        """
        net_insertion = mixed_state(3, 1)  # Delta = +2, magnitude 2, not in the shipped 3n+1 series
        net_deletion = mixed_state(1, 4)  # Delta = -3, magnitude 3, not in the shipped 3n+2 series

        # Under the shipped settings both are dropped, which is what makes the contrast real.
        assert surviving_rows(net_insertion) == 0
        assert surviving_rows(net_deletion) == 0

        monkeypatch.setattr(advntr, "advntr_settings", {"max_frameshift": 100, "frameshift_multiplier": 1})

        assert len(advntr.advntr_processing_ins(variant_frame(net_insertion))) == 1
        assert len(advntr.advntr_processing_del(variant_frame(net_insertion))) == 0, "a net insertion, either way"

        assert len(advntr.advntr_processing_del(variant_frame(net_deletion))) == 1
        assert len(advntr.advntr_processing_ins(variant_frame(net_deletion))) == 0, "a net deletion, either way"

    def test_the_two_offsets_are_the_ones_the_convention_names(self):
        """The insertion arm is the ``3n+1`` series and the deletion arm the ``3n+2`` one."""
        assert advntr.INSERTION_FRAME_OFFSET == 1
        assert advntr.DELETION_FRAME_OFFSET == 2

        ins_frame = advntr.accepted_frame_magnitudes(advntr.INSERTION_FRAME_OFFSET)
        del_frame = advntr.accepted_frame_magnitudes(advntr.DELETION_FRAME_OFFSET)

        assert list(ins_frame[:4]) == ["1", "4", "7", "10"]
        assert list(del_frame[:4]) == ["2", "5", "8", "11"]
        assert set(ins_frame).isdisjoint(del_frame)
