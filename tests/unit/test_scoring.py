# tests/unit/test_scoring.py

"""
Unit tests for the scoring functionality in vntyper/scripts/scoring.py.

Validates frame-score calculations, depth splitting, and frameshift extraction.
"""

import pandas as pd
import pytest

from vntyper.scripts.confidence_assignment import calculate_depth_score_and_assign_confidence
from vntyper.scripts.scoring import (
    extract_frameshifts,
    split_depth_and_calculate_frame_score,
    split_frame_score,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


def _score_chain(ref: str, alt: str, sample: str = "Del:10:100") -> pd.Series:
    """Run one variant through the three scoring stages in pipeline order.

    Args:
        ref: Reference allele.
        alt: Alternate allele.
        sample: The colon-delimited ``Sample`` field Kestrel emits.

    Returns:
        pd.Series: The single output row after
        ``split_depth_and_calculate_frame_score`` ->
        ``split_frame_score`` -> ``extract_frameshifts``.
    """
    frame = pd.DataFrame(
        {
            "Sample": [sample],
            "REF": [ref],
            "ALT": [alt],
            "POS": [67],
            "Motifs": ["X-5"],
            "Variant": ["v"],
            "Motif_sequence": ["SEQ1"],
        }
    )
    return extract_frameshifts(split_frame_score(split_depth_and_calculate_frame_score(frame))).loc[0]


@pytest.mark.parametrize(
    "df_input,expected_len",
    [
        (pd.DataFrame(), 0),
    ],
)
def test_split_depth_and_calculate_frame_score_empty_df(df_input, expected_len):
    """Verify that an empty input DataFrame remains empty."""
    out = split_depth_and_calculate_frame_score(df_input)
    assert len(out) == expected_len, (
        "Empty input should yield empty output after split_depth_and_calculate_frame_score."
    )


def test_split_depth_and_calculate_frame_score_no_frameshift():
    """
    Test that non-frameshift variants (multiple of 3 difference) are correctly flagged.

    When the difference (ALT length - REF length) is a multiple of 3,
    the variant is marked with is_frameshift=False but retained in the DataFrame.
    This follows the flag-and-defer pattern where filtering is done downstream.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "Sample": ["Del:10:100"],  # Only the first 'Del' part is not used, but we keep format for test
            "REF": ["ATG"],  # length 3
            "ALT": ["ATGATG"],  # length 6  -> difference = 3 -> multiple of 3
            "Motifs": ["mock_motif"],
            "Variant": ["mock_variant"],
            "POS": [123],
            "Motif_sequence": ["mock_sequence"],
        }
    )

    # Act
    out = split_depth_and_calculate_frame_score(df)

    # Assert - row is retained with is_frameshift flag
    assert len(out) == 1, "Row should be retained in the DataFrame"
    assert "is_frameshift" in out.columns, "Should have 'is_frameshift' column"
    assert not out.loc[0, "is_frameshift"], (
        "Variant with multiple-of-3 difference should be marked as is_frameshift=False"
    )

    # Assert - Frame_Score calculation is correct
    assert "Frame_Score" in out.columns, "Should have 'Frame_Score' column"
    expected_frame_score = (6 - 3) / 3  # (alt_len - ref_len) / 3 = 1.0
    assert out.loc[0, "Frame_Score"] == expected_frame_score, f"Frame_Score should be {expected_frame_score}"


def test_split_depth_and_calculate_frame_score_frameshift():
    """
    Test frameshift variant detection and Frame_Score calculation.

    If the difference (ALT length - REF length) is not a multiple of 3,
    the variant should be retained and a 'Frame_Score' should be added.
    """
    df = pd.DataFrame(
        {
            "Sample": ["Del:50:500"],
            "REF": ["ATG"],  # length 3
            "ALT": ["ATGA"],  # length 4 -> difference = 1 -> frameshift
            "Motifs": ["mock_motif"],
            "Variant": ["mock_variant"],
            "POS": [456],
            "Motif_sequence": ["mock_sequence"],
        }
    )
    out = split_depth_and_calculate_frame_score(df)
    assert not out.empty, "Expected to retain a frameshift variant (difference not multiple of 3)."
    assert "Frame_Score" in out.columns, "Output should have a 'Frame_Score' column."
    # Check that is_frameshift was True
    assert "is_frameshift" in out.columns, "Output should have 'is_frameshift' marking frameshift or not."
    assert all(out["is_frameshift"]), "All retained rows should be frameshift variants."


def test_split_frame_score_empty_df():
    """Verify that an empty input DataFrame remains empty when split_frame_score is called."""
    df = pd.DataFrame()
    out = split_frame_score(df)
    assert out.empty, "Empty input should yield empty output after split_frame_score."


def test_split_frame_score_basic():
    """Test basic splitting of frame score into 'direction' and 'frameshift_amount'."""
    df = pd.DataFrame(
        {
            "Frame_Score": [1.0, -2.0],  # not directly used, but indicates frameshift
            "ref_len": [3, 6],
            "alt_len": [4, 4],  # alt_len - ref_len => [1, -2]
            "is_frameshift": [
                True,
                True,
            ],  # frameshift is assumed True from previous step
        }
    )
    out = split_frame_score(df)

    # We drop 'is_frameshift', 'ref_len', 'alt_len'
    # We keep 'direction', 'frameshift_amount', 'Frame_Score', etc.
    expected_columns = {"Frame_Score", "direction", "frameshift_amount"}
    assert expected_columns.issubset(set(out.columns)), f"Output must contain at least: {expected_columns}"

    # direction = sign(alt_len - ref_len)
    # frameshift_amount = abs(alt_len - ref_len) % 3
    # For row0: alt_len - ref_len = 1 => direction=1, frameshift_amount=1
    # For row1: alt_len - ref_len = -2 => direction < 0 => -1, frameshift_amount=2
    assert out.loc[0, "direction"] == 1, "Expected direction=1 for alt_len-ref_len=1."
    assert out.loc[0, "frameshift_amount"] == 1, "Expected frameshift_amount=1 for difference=1."

    assert out.loc[1, "direction"] == -1, "Expected direction=-1 for alt_len-ref_len=-2."
    assert out.loc[1, "frameshift_amount"] == 2, "Expected frameshift_amount=2 for difference=-2."


def test_extract_frameshifts_empty_df():
    """Verify that an empty input DataFrame remains empty in extract_frameshifts."""
    df = pd.DataFrame()
    out = extract_frameshifts(df)
    assert out.empty, "Empty input should yield empty output after extract_frameshifts."


def test_extract_frameshifts_mixed():
    """
    Test that frameshift variants are correctly flagged based on insertion/deletion patterns.

    Valid frameshift patterns:
    - Insertion frameshift: direction > 0 AND frameshift_amount == 1 (3n+1)
    - Deletion frameshift: direction < 0 AND frameshift_amount == 2 (3n+2)

    All rows are retained but marked with is_valid_frameshift flag.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "direction": [1, 1, -1, -1, 1],
            "frameshift_amount": [1, 2, 2, 1, 1],
            "Frame_Score": [0.33, 0.66, -0.66, -0.33, 0.33],
            "Variant": ["ins_ok", "ins_wrong", "del_ok", "del_wrong", "ins_ok2"],
        }
    )

    # Act
    out = extract_frameshifts(df)

    # Assert - all rows are retained
    assert len(out) == 5, "All 5 rows should be retained in the DataFrame"

    # Assert - is_valid_frameshift column exists
    assert "is_valid_frameshift" in out.columns, "Should have 'is_valid_frameshift' column"

    # Assert - correct rows are marked as valid frameshifts
    # Row 0: direction=1, frameshift_amount=1 => valid insertion (3n+1)
    assert out.loc[0, "is_valid_frameshift"], "Row 0 (ins_ok) should be valid: direction>0 and frameshift_amount=1"

    # Row 1: direction=1, frameshift_amount=2 => invalid (wrong frameshift amount for insertion)
    assert not out.loc[1, "is_valid_frameshift"], (
        "Row 1 (ins_wrong) should be invalid: direction>0 but frameshift_amount=2"
    )

    # Row 2: direction=-1, frameshift_amount=2 => valid deletion (3n+2)
    assert out.loc[2, "is_valid_frameshift"], "Row 2 (del_ok) should be valid: direction<0 and frameshift_amount=2"

    # Row 3: direction=-1, frameshift_amount=1 => invalid (wrong frameshift amount for deletion)
    assert not out.loc[3, "is_valid_frameshift"], (
        "Row 3 (del_wrong) should be invalid: direction<0 but frameshift_amount=1"
    )

    # Row 4: direction=1, frameshift_amount=1 => valid insertion (3n+1)
    assert out.loc[4, "is_valid_frameshift"], "Row 4 (ins_ok2) should be valid: direction>0 and frameshift_amount=1"

    # Assert - verify the count of valid frameshifts
    valid_count = out["is_valid_frameshift"].sum()
    assert valid_count == 3, "Should have exactly 3 valid frameshifts (rows 0, 2, 4)"


# ---------------------------------------------------------------------------
# Characterisation tests (issue #179)
#
# Everything below records what the code does today so that a change becomes
# visible. None of it is a specification, and none of it should be read as a
# claim that the current behaviour is correct.
# ---------------------------------------------------------------------------


def test_the_sample_field_splits_into_alt_depth_then_region_depth():
    """``Sample`` is ``Del:<alt depth>:<active-region depth>``, in that order.

    The two depth column names are assigned positionally from the split, and
    nothing downstream re-checks the assignment. Swapping the two names in
    ``split_depth_and_calculate_frame_score`` leaves every line covered, every
    existing test green, and inverts ``Depth_Score`` for every variant in the
    run -- which inverts every confidence label the pipeline reports. This test
    is the only thing standing between that swap and a clinical report.

    It also pins that both columns arrive as **strings**: the split produces
    object dtype, and ``calculate_depth_score_and_assign_confidence`` is what
    coerces them to numbers.
    """
    out = split_depth_and_calculate_frame_score(
        pd.DataFrame(
            {
                "Sample": ["Del:10:100"],
                "REF": ["G"],
                "ALT": ["GG"],
                "POS": [67],
                "Motifs": ["X-5"],
                "Variant": ["Insertion"],
                "Motif_sequence": ["SEQ1"],
            }
        )
    )

    alt_depth = out.loc[0, "Estimated_Depth_AlternateVariant"]
    region_depth = out.loc[0, "Estimated_Depth_Variant_ActiveRegion"]

    assert out.loc[0, "Del"] == "Del", "the first colon-delimited field is the literal 'Del' tag"
    assert alt_depth == "10", "the second field is the alternate-allele depth"
    assert region_depth == "100", "the third field is the variant-active-region depth"
    assert int(alt_depth) == 10
    assert int(region_depth) == 100
    assert int(alt_depth) < int(region_depth), (
        "the alternate depth is a subset of the active-region depth, so it can never be the larger of the two"
    )


def test_the_sample_split_order_decides_the_direction_of_depth_score():
    """The consequence of the split order, asserted end to end.

    ``Depth_Score`` is ``alt / region``. With ``Sample='Del:10:100'`` that is
    0.1; with the two column names swapped it would be 10. The two land on
    opposite sides of every threshold in ``kestrel_config.json``, so this is the
    same defect as the test above seen from the other end.
    """
    scored = split_depth_and_calculate_frame_score(
        pd.DataFrame(
            {
                "Sample": ["Del:10:100"],
                "REF": ["G"],
                "ALT": ["GG"],
                "POS": [67],
                "Motifs": ["X-5"],
                "Variant": ["Insertion"],
                "Motif_sequence": ["SEQ1"],
            }
        )
    )
    from tests.builders import kestrel_config

    out = calculate_depth_score_and_assign_confidence(scored, kestrel_config())
    assert out.loc[0, "Depth_Score"] == pytest.approx(0.1)
    assert out.loc[0, "Depth_Score"] < 1.0, "a depth score above 1 means the split order was inverted"


@pytest.mark.parametrize("frameshift_amount", [0, 1, 2])
def test_direction_zero_is_never_a_valid_frameshift(frameshift_amount):
    """``direction == 0`` fails both patterns, whatever the frameshift amount.

    ``extract_frameshifts`` requires ``direction > 0`` for an insertion and
    ``direction < 0`` for a deletion. Zero satisfies neither, so a length-preserving
    change can never be a valid frameshift. Widening either comparison to ``>=`` or
    ``<=`` would admit substitutions; this pins that it does not happen.

    Args:
        frameshift_amount: The modulo-3 remainder to pair with ``direction == 0``.
    """
    out = extract_frameshifts(
        pd.DataFrame({"direction": [0], "frameshift_amount": [frameshift_amount]}),
    )
    assert not out.loc[0, "is_valid_frameshift"]


def test_a_substitution_yields_direction_zero_and_is_discarded():
    """The end-to-end form of the case above: ``REF == ALT`` reaches direction 0."""
    row = _score_chain(ref="AT", alt="AT")

    assert row["direction"] == 0
    assert row["frameshift_amount"] == 0
    assert not row["is_frameshift"]
    assert not row["is_valid_frameshift"]


@pytest.mark.parametrize("deleted_bases", [1, 4, 7])
def test_a_3n_plus_1_bp_deletion_is_a_frameshift_that_the_rule_discards(deleted_bases):
    """Characterisation: a (3n+1)-bp deletion is a frameshift but not a *valid* one.

    ``REF='AT'``, ``ALT='A'`` deletes one base. ``is_frameshift`` is True, because
    ``(alt_len - ref_len) % 3 != 0``. ``is_valid_frameshift`` is **False**, because
    ``extract_frameshifts`` requires ``frameshift_amount == 2`` in the deletion
    direction and a 1-bp deletion has ``frameshift_amount == 1``.

    **This contradicts the general definition of a frameshift.** Any indel whose
    length is not a multiple of three shifts the reading frame; the deletion arm of
    this rule accepts only half of those. ``filter_final_dataframe`` requires
    ``is_valid_frameshift``, so every (3n+1)-bp deletion is dropped from the final
    result with no warning and no row in the report. That is a candidate
    silent-false-negative source, and no fixture in this repository exercised a
    1-bp or (3n+1)-bp deletion before this test.

    The rule is deliberately left alone: the thresholds and patterns are calibrated
    against Saei et al. (2023), and changing which variants survive needs a
    golden-cohort diff, not a unit test. This records the behaviour so the decision
    is visible.

    Args:
        deleted_bases: How many bases the deletion removes; each is ``3n + 1``.
    """
    assert deleted_bases % 3 == 1

    ref = "A" * (deleted_bases + 1)
    row = _score_chain(ref=ref, alt="A")

    assert row["Frame_Score"] == pytest.approx(-deleted_bases / 3)
    assert row["direction"] == -1
    assert row["frameshift_amount"] == 1
    assert row["is_frameshift"], "a non-multiple-of-three deletion does shift the reading frame"
    assert not row["is_valid_frameshift"], (
        "and yet the deletion arm of the rule requires frameshift_amount == 2, so this row is discarded"
    )


@pytest.mark.parametrize("deleted_bases", [2, 5, 8])
def test_a_3n_plus_2_bp_deletion_is_the_only_deletion_the_rule_accepts(deleted_bases):
    """The other half of the pair above: ``frameshift_amount == 2`` deletions survive.

    Asserting both halves is what makes the asymmetry a finding rather than an
    observation -- the rule is not "deletions are dropped", it is "half of the
    frameshifting deletions are dropped".

    Args:
        deleted_bases: How many bases the deletion removes; each is ``3n + 2``.
    """
    assert deleted_bases % 3 == 2

    ref = "A" * (deleted_bases + 1)
    row = _score_chain(ref=ref, alt="A")

    assert row["direction"] == -1
    assert row["frameshift_amount"] == 2
    assert row["is_frameshift"]
    assert row["is_valid_frameshift"]


@pytest.mark.parametrize(
    ("inserted_bases", "expected_valid"),
    [(1, True), (2, False), (4, True), (5, False)],
)
def test_the_insertion_arm_accepts_3n_plus_1_only(inserted_bases, expected_valid):
    """The insertion arm mirrors the deletion arm, on the opposite remainder.

    Insertions are accepted at ``frameshift_amount == 1`` and deletions at
    ``frameshift_amount == 2``, so the two arms are not symmetric in the remainder
    they accept. Both are pinned so a change to either is visible.

    Args:
        inserted_bases: How many bases the insertion adds.
        expected_valid: Today's ``is_valid_frameshift`` for that insertion.
    """
    alt = "A" * (inserted_bases + 1)
    row = _score_chain(ref="A", alt=alt)

    assert row["direction"] == 1
    assert bool(row["is_frameshift"]) is (inserted_bases % 3 != 0)
    assert bool(row["is_valid_frameshift"]) is expected_valid
