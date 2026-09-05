"""Boundary characterisation for ``confidence_assignment``.

Why this file exists
--------------------
``calculate_depth_score_and_assign_confidence`` produces the ``Confidence``
string a clinician reads. Before this file it had 100% line coverage and a 21%
mutation score: inverting the central condition (``cond1``'s ``==`` to ``!=``)
left the whole suite green, and five of its six threshold comparisons could be
changed without a single test failing. Line coverage says the code ran; it says
nothing about whether the *boundary* was pinned.

The 54-cell boundary matrix below is **characterisation**: it records the
threshold arithmetic so a changed comparison operator becomes visible, and
makes no claim that any individual cutoff is right.

The *ordering* is different. @hassansaei decided on #183 (2026-08-06) that
2.x's last-wins sequential assignment is the intended behaviour and that
1.3's absolute region-depth cap must not be restored. Tests that pin the
order of the ``df.loc`` assignments are therefore **specification**, and
changing that order requires a new decision on #183.

What is pinned
--------------
* A 54-cell matrix over the cross product of the three alt-depth thresholds
  (each probed at threshold-1, threshold and threshold+1) and the two
  depth-score thresholds (each probed just below, exactly at and just above).
  Every threshold is read from the shipped config through
  ``tests.builders.kestrel_config``; none is hardcoded. This part is
  characterisation. #184 rewrote five of its cells -- see the comment above
  ``_EXPECTED_MATRIX`` -- but the matrix as a whole still makes no claim that
  any individual cutoff is right.
* The closed mid-band ``[0.00469, 0.00515]``, which #184 decided must be
  Low_Precision at every alternate depth. That part **is** specification: the
  tests under the "#184" banner at the foot of this file each quote
  @hassansaei and cite the issue. The rest of the matrix is not covered by
  that decision.
* The region-depth tier interaction -- see
  ``test_region_depth_demotion_is_overwritten_by_a_later_high_precision_tier``
  and ``test_a_low_region_depth_row_is_deliberately_not_capped_at_low_precision``.
  Its *ordering* (last-wins) is specification, per #183.

Notes for anyone editing this file
----------------------------------
* ``High_Precision*`` is a **literal label**. The trailing star is part of the
  value, not a wildcard. ``test_the_star_in_high_precision_star_is_literal``
  pins that.
* ``calculate_depth_score_and_assign_confidence`` reads its ``kestrel_config``
  argument, unlike ``run_kestrel`` which reads the module-level global. Passing
  a config in here therefore does something; passing one into ``run_kestrel``
  does not.
"""

import math
from fractions import Fraction

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.scripts.confidence_assignment import (
    NEGATIVE_LABEL,
    calculate_depth_score_and_assign_confidence,
)
from vntyper.scripts.kestrel_genotyping import select_single_best_variant
from vntyper.scripts.scoring import split_depth_and_calculate_frame_score

pytestmark = pytest.mark.unit

_CONF = kestrel_config()["confidence_assignment"]
_LOW = _CONF["depth_score_thresholds"]["low"]
_HIGH = _CONF["depth_score_thresholds"]["high"]
_ALT_LOW = _CONF["alt_depth_thresholds"]["low"]
_ALT_MID_LOW = _CONF["alt_depth_thresholds"]["mid_low"]
_ALT_MID_HIGH = _CONF["alt_depth_thresholds"]["mid_high"]
_REGION_THRESHOLD = _CONF["var_active_region_threshold"]

_LOW_LABEL = _CONF["confidence_levels"]["low_precision"]
_HIGH_LABEL = _CONF["confidence_levels"]["high_precision"]
_HIGH_STAR_LABEL = _CONF["confidence_levels"]["high_precision_star"]

# Step used for the "just below" / "just above" probes. It has to be small
# enough to stay inside the band between the two depth-score thresholds and
# large enough to survive the float round trip through `alt / region` -- the gap
# is ~4.6e-4 and the round-trip error is ~1e-18, so 1e-6 clears both by orders
# of magnitude.
_EPS = 1e-6

# Column order of the matrix below.
_ALT_DEPTHS = {
    "alt_low-1": _ALT_LOW - 1,
    "alt_low": _ALT_LOW,
    "alt_low+1": _ALT_LOW + 1,
    "mid_low-1": _ALT_MID_LOW - 1,
    "mid_low": _ALT_MID_LOW,
    "mid_low+1": _ALT_MID_LOW + 1,
    "mid_high-1": _ALT_MID_HIGH - 1,
    "mid_high": _ALT_MID_HIGH,
    "mid_high+1": _ALT_MID_HIGH + 1,
}

# Row order of the matrix below.
_DEPTH_SCORES = {
    "low-eps": _LOW - _EPS,
    "low": _LOW,
    "low+eps": _LOW + _EPS,
    "high-eps": _HIGH - _EPS,
    "high": _HIGH,
    "high+eps": _HIGH + _EPS,
}

_N = NEGATIVE_LABEL
_L = _LOW_LABEL
_H = _HIGH_LABEL
_S = _HIGH_STAR_LABEL

# The snapshot. Rows are Depth_Score probes, columns are alt-depth probes, both
# in the declaration order above. This is what the code emits today.
#
# The `high` row was rewritten by #184: the mid-band demotion now covers the CLOSED
# interval [low, high] and is applied last, so `Depth_Score == high` exactly is
# Low_Precision at every alt depth. Five cells moved -- alt_low+1, mid_low, mid_low+1,
# mid_high and mid_high+1. The sixth, mid_high-1, is the one cell in the row whose
# constructed region depth cannot land exactly on the threshold (see below): it sits
# one ULP ABOVE `high`, outside the closed band, so it is still promoted. That is not
# an inconsistency -- it is the tightest evidence in this file that the band's upper
# edge is exactly at `high` and not one step past it.
#
#                  alt_low-1  alt_low  alt_low+1  mid_low-1  mid_low  mid_low+1  mid_high-1  mid_high  mid_high+1
_EXPECTED_MATRIX = {
    "low-eps": (_N, _N, _N, _N, _N, _N, _N, _N, _N),
    "low": (_L, _L, _L, _L, _L, _L, _L, _L, _L),
    "low+eps": (_L, _L, _L, _L, _L, _L, _L, _L, _L),
    "high-eps": (_L, _L, _L, _L, _L, _L, _L, _L, _L),
    "high": (_L, _L, _L, _L, _L, _L, _H, _L, _L),
    "high+eps": (_L, _L, _H, _L, _H, _H, _H, _S, _S),
}

# `Depth_Score` is computed by production as `alt / region`, so each cell supplies
# a region depth of `alt / target`. That round trip lands on the target exactly for
# 53 of the 54 cells. The exception is alt=99 at Depth_Score == high: no double `r`
# satisfies `99 / r == 0.00515`, because one ULP of `r` near 19223 is wider than the
# pre-image interval of the target. That cell lands one ULP above `high`, which is
# on the same side of every comparison the function makes there. The row assertion
# below is written as `>= high` for that reason and only for that reason.
_CELLS_THAT_CANNOT_LAND_EXACTLY_ON_THE_THRESHOLD = {("mid_high-1", "high")}


def _assert_depth_score_landed_where_intended(achieved: float, ds_id: str, alt_id: str) -> None:
    """Fail if the constructed region depth did not put Depth_Score where the cell claims.

    Without this the matrix could pass for the wrong reason: a region depth that
    silently lands one side of a threshold instead of on it would still produce
    *a* label, and the snapshot would then record a boundary that was never
    probed.

    Args:
        achieved: The ``Depth_Score`` the production function computed.
        ds_id: The row key from :data:`_DEPTH_SCORES`.
        alt_id: The column key from :data:`_ALT_DEPTHS`.
    """
    if ds_id == "low-eps":
        assert achieved < _LOW, f"{alt_id}/{ds_id}: {achieved!r} should be below the low threshold {_LOW!r}"
    elif ds_id == "low":
        assert achieved == _LOW, f"{alt_id}/{ds_id}: {achieved!r} should be exactly the low threshold {_LOW!r}"
    elif ds_id in ("low+eps", "high-eps"):
        assert _LOW < achieved < _HIGH, f"{alt_id}/{ds_id}: {achieved!r} should sit strictly between the thresholds"
    elif ds_id == "high":
        assert achieved >= _HIGH, f"{alt_id}/{ds_id}: {achieved!r} should reach the high threshold {_HIGH!r}"
        if (alt_id, ds_id) not in _CELLS_THAT_CANNOT_LAND_EXACTLY_ON_THE_THRESHOLD:
            assert achieved == _HIGH, f"{alt_id}/{ds_id}: {achieved!r} should be exactly {_HIGH!r}"
        else:
            assert achieved != _HIGH, (
                f"{alt_id}/{ds_id} is recorded as unreachable, but {achieved!r} landed exactly on "
                f"{_HIGH!r}. Drop it from _CELLS_THAT_CANNOT_LAND_EXACTLY_ON_THE_THRESHOLD."
            )
            assert abs(achieved - _HIGH) <= math.ulp(_HIGH), (
                f"{alt_id}/{ds_id}: {achieved!r} is more than one ULP from {_HIGH!r}"
            )
    elif ds_id == "high+eps":
        assert achieved > _HIGH, f"{alt_id}/{ds_id}: {achieved!r} should be above the high threshold {_HIGH!r}"
    else:  # pragma: no cover - guards against an unlabelled row being added
        raise AssertionError(f"no landing assertion defined for depth-score probe {ds_id!r}")


def _assign(alt_depth: float, region_depth: float) -> pd.Series:
    """Run one row through the production function and return it.

    Args:
        alt_depth: ``Estimated_Depth_AlternateVariant``.
        region_depth: ``Estimated_Depth_Variant_ActiveRegion``.

    Returns:
        pd.Series: The single output row.
    """
    frame = pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [alt_depth],
            "Estimated_Depth_Variant_ActiveRegion": [region_depth],
        }
    )
    return calculate_depth_score_and_assign_confidence(frame, kestrel_config()).loc[0]


def test_the_star_in_high_precision_star_is_literal():
    """The trailing ``*`` is part of the label, not a wildcard.

    A helper in this repository used to treat it as a prefix wildcard, which made
    every assertion of ``High_Precision*`` silently accept plain
    ``High_Precision``. The three labels are written out here so a config edit
    that renames one is caught here rather than in a report a clinician reads.
    """
    assert _LOW_LABEL == "Low_Precision"
    assert _HIGH_LABEL == "High_Precision"
    assert _HIGH_STAR_LABEL == "High_Precision*"
    assert _HIGH_STAR_LABEL != _HIGH_LABEL
    assert NEGATIVE_LABEL == "Negative"


def test_the_matrix_covers_every_cell():
    """The snapshot must describe all 54 cells, or the parametrisation lies.

    A row with the wrong width would make ``_EXPECTED_MATRIX[ds][index]`` read a
    neighbouring cell instead of failing.
    """
    assert len(_ALT_DEPTHS) == 9
    assert len(_DEPTH_SCORES) == 6
    assert set(_EXPECTED_MATRIX) == set(_DEPTH_SCORES)
    for ds_id, row in _EXPECTED_MATRIX.items():
        assert len(row) == len(_ALT_DEPTHS), f"row {ds_id!r} has {len(row)} cells, expected {len(_ALT_DEPTHS)}"


@pytest.mark.parametrize("ds_id", list(_DEPTH_SCORES))
@pytest.mark.parametrize("alt_id", list(_ALT_DEPTHS))
def test_confidence_boundary_matrix(alt_id, ds_id):
    """Pin the ``Confidence`` label at every threshold crossing.

    Each cell fixes an alt depth one step below, exactly at, or one step above one
    of the three alt-depth thresholds, and a region depth chosen so ``Depth_Score``
    lands just below, exactly on, or just above one of the two depth-score
    thresholds. The region depth is always far above
    ``var_active_region_threshold``, so ``cond1``'s region clause is provably out
    of play and each cell isolates the depth-score and alt-depth comparisons.

    Args:
        alt_id: Column key from :data:`_ALT_DEPTHS`.
        ds_id: Row key from :data:`_DEPTH_SCORES`.
    """
    alt_depth = _ALT_DEPTHS[alt_id]
    target_depth_score = _DEPTH_SCORES[ds_id]
    region_depth = alt_depth / target_depth_score

    assert region_depth > _REGION_THRESHOLD, (
        f"{alt_id}/{ds_id}: region depth {region_depth!r} is not above {_REGION_THRESHOLD!r}, so cond1's "
        "region clause would fire and this cell would no longer isolate the depth-score comparison"
    )

    row = _assign(alt_depth, region_depth)
    _assert_depth_score_landed_where_intended(row["Depth_Score"], ds_id, alt_id)

    expected = _EXPECTED_MATRIX[ds_id][list(_ALT_DEPTHS).index(alt_id)]
    assert row["Confidence"] == expected, (
        f"alt_depth={alt_depth} ({alt_id}), Depth_Score={row['Depth_Score']!r} ({ds_id}): "
        f"expected {expected!r}, got {row['Confidence']!r}"
    )
    assert row["depth_confidence_pass"] == (expected != NEGATIVE_LABEL)


def test_depth_score_below_the_low_threshold_is_always_negative():
    """The whole ``low-eps`` row is Negative, whatever the alt depth.

    Stated once as a whole-row claim so that a change turning the guard into a
    per-alt-depth rule fails here with one clear message rather than as nine
    separate cell failures.
    """
    for alt_id, alt_depth in _ALT_DEPTHS.items():
        row = _assign(alt_depth, alt_depth / (_LOW - _EPS))
        assert row["Confidence"] == NEGATIVE_LABEL, f"{alt_id}: expected Negative below the low threshold"
        assert not row["depth_confidence_pass"]


def test_region_depth_demotion_is_overwritten_by_a_later_high_precision_tier():
    """Specification: a low region depth does **not** cap the confidence label.

    ``cond1`` demotes a row whose ``Estimated_Depth_Variant_ActiveRegion`` is at or
    below ``var_active_region_threshold`` to Low_Precision, but it is applied
    *first*; ``cond2`` (High_Precision*) and ``cond5`` (High_Precision) are applied
    after it and overwrite the demotion outright.

    **Tier precedence is specified.** @hassansaei decided on #183 (2026-08-06)
    that 2.x's last-wins sequential assignment -- the six conditions applied in
    source order, later conditions overwriting earlier ones -- is the intended
    behaviour going forward, and that 1.3's absolute region-depth cap must not
    be restored. This test pins that order: do read it as a claim that a
    High_Precision call on a 150-read active region is correct. Changing the
    order requires a new decision on #183.

    ``docs/pipeline/scoring-and-confidence.md`` used to state that High_Precision
    requires a region depth above 200. The code is correct and the doc table was
    the artefact that was wrong; the doc has been corrected to describe the code.
    """
    below = _REGION_THRESHOLD - 50
    assert below <= _REGION_THRESHOLD

    # alt depth inside [mid_low, mid_high) with a depth score far above `high`.
    demoted_then_promoted = _assign(50, below)
    assert demoted_then_promoted["Depth_Score"] >= _HIGH
    assert demoted_then_promoted["Confidence"] == _HIGH_LABEL

    # alt depth at or above mid_high: the demotion is overwritten by cond2 instead.
    demoted_then_starred = _assign(150, below)
    assert demoted_then_starred["Depth_Score"] >= _HIGH
    assert demoted_then_starred["Confidence"] == _HIGH_STAR_LABEL

    # The same alt depth on an ample region gets the same label, which is the
    # point: the region-depth tier is invisible in the emitted label.
    ample = _assign(50, 5000)
    assert ample["Estimated_Depth_Variant_ActiveRegion"] > _REGION_THRESHOLD
    assert ample["Confidence"] == _HIGH_LABEL
    assert ample["Confidence"] == demoted_then_promoted["Confidence"]


def test_the_region_depth_clause_is_only_reachable_through_the_alt_threshold_gap():
    """Characterisation: ``cond1``'s region clause decides a label in one narrow case.

    ``cond1`` is applied before every promotion, so its region-depth branch only
    survives to the output when no later condition fires at all. That needs an alt
    depth strictly between ``alt_depth_thresholds.low`` (20) and
    ``alt_depth_thresholds.mid_low`` (21) together with a depth score at or above
    ``high``: ``cond4`` wants ``<= 20``, ``cond3``/``cond5`` want ``>= 21``, and
    ``cond6`` wants a score below ``high``, so all four miss.

    The config leaves that gap open. Alt depths arrive as integers from the
    ``Sample`` split, so nothing in the shipped pipeline lands inside it -- but the
    gap is real in the config, and this is the only place the region threshold is
    observable in the final label at all.
    """
    assert _ALT_MID_LOW - _ALT_LOW == 1, "the gap this test probes only exists while mid_low is alt_low + 1"

    in_the_gap = (_ALT_LOW + _ALT_MID_LOW) / 2
    at_threshold = _assign(in_the_gap, _REGION_THRESHOLD)
    assert at_threshold["Depth_Score"] >= _HIGH
    assert at_threshold["Confidence"] == _LOW_LABEL

    just_above_threshold = _assign(in_the_gap, _REGION_THRESHOLD + 1)
    assert just_above_threshold["Depth_Score"] >= _HIGH
    assert just_above_threshold["Confidence"] == NEGATIVE_LABEL, (
        "with the region clause off and every alt-depth condition missing, no condition assigns a label "
        "and the row keeps the default"
    )


def test_a_low_region_depth_row_is_deliberately_not_capped_at_low_precision():
    """alt=50 on a 150-read region reports High_Precision, and that is intended.

    Specification, not characterisation. v1.3 made region depth an absolute
    first-wins cap; 2.x applies the tiers as sequential ``df.loc`` assignments,
    so ``cond1``'s demotion is overwritten by a later match. @hassansaei on
    #183 decided this is the intended behaviour going forward and that the
    1.3 cap must not be restored: the pattern is rare, and the cases where it
    appears are already caught by the Depth_Score flagging rule.
    """
    df = pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [50],
            "Estimated_Depth_Variant_ActiveRegion": [150],
        }
    )

    out = calculate_depth_score_and_assign_confidence(df, kestrel_config())

    assert out["Confidence"].iloc[0] == "High_Precision"
    assert out["Estimated_Depth_Variant_ActiveRegion"].iloc[0] <= 200


# ---------------------------------------------------------------------------
# #184: the mid-band Depth_Score demotion takes precedence over the High tiers.
# ---------------------------------------------------------------------------

#: Every integer ``(alt, region)`` pair, up to alt=1030, whose float quotient is
#: *exactly* the shipped ``high`` threshold. ``0.00515 == 103 / 20000`` in IEEE 754
#: (the division is exact for this pair), so ``alt`` is always a multiple of 103 and
#: therefore always >= 103 > ``alt_mid_high`` (100). ``cond5``'s ``[21, 100)`` window
#: is arithmetically unreachable here; only ``cond2`` can fire, so the only label the
#: #184 change moves on integer depths is ``High_Precision*``.
#:
#: Enumerated exhaustively, not guessed: for each alt in 1..12000 the only region
#: depths that can quotient to within one ULP of the threshold are floor(alt/high)-1,
#: floor(alt/high) and floor(alt/high)+1, and the search over all three found 116
#: pairs, every one of them with alt a multiple of 103. The first ten are listed.
EXACT_HIGH_THRESHOLD_PAIRS = [
    (103, 20000),
    (206, 40000),
    (309, 60000),
    (412, 80000),
    (515, 100000),
    (618, 120000),
    (721, 140000),
    (824, 160000),
    (927, 180000),
    (1030, 200000),
]


def test_the_enumerated_boundary_pairs_all_land_exactly_on_the_high_threshold():
    """The table above is only evidence if every pair in it really hits the boundary.

    A pair whose quotient lands one ULP off would still produce *a* label, and the
    parametrised tests below would then be pinning a neighbouring cell while
    claiming to pin the boundary.
    """
    for alt, region in EXACT_HIGH_THRESHOLD_PAIRS:
        assert alt / region == _HIGH, f"({alt}, {region}) gives {alt / region!r}, not {_HIGH!r}"
        assert region > _REGION_THRESHOLD, f"({alt}, {region}) would also trip cond1's region clause"


@pytest.mark.parametrize(("alt", "region"), EXACT_HIGH_THRESHOLD_PAIRS)
def test_the_top_of_the_mid_band_is_low_precision(alt, region):
    """``Depth_Score == 0.00515`` exactly must be Low_Precision (#184).

    **Specification.** @hassansaei on #184: "Any variant with Depth_Score between
    0.00469 and 0.00515 (inclusive) must be Low_Precision, even when alternate
    depth is >= 21 (or higher). That mid-range Depth_Score demotion from 1.3 is
    still the desired behaviour."

    ``cond6`` already demoted the OPEN interval at every alt depth because it was
    applied last, so this boundary was the only divergence from his intent. Before
    this change ``cond2`` promoted these rows to ``High_Precision*``.

    Args:
        alt: ``Estimated_Depth_AlternateVariant`` from :data:`EXACT_HIGH_THRESHOLD_PAIRS`.
        region: The matching ``Estimated_Depth_Variant_ActiveRegion``.
    """
    assert alt / region == _HIGH, "this pair does not reach the boundary in IEEE 754"

    row = _assign(alt, region)

    assert row["Depth_Score"] == _HIGH
    assert row["Confidence"] == _LOW_LABEL
    assert row["depth_confidence_pass"]


def test_the_bottom_of_the_mid_band_is_low_precision():
    """``Depth_Score == 0.00469`` exactly must be Low_Precision (#184).

    **Specification**, same quote as above -- the band is inclusive at *both*
    ends. ``cond1``'s ``Depth_Score == low_threshold`` clause already handled this
    end and nothing overwrote it, so this test pins existing behaviour rather than
    changing it. It is here because the new mid-band rule now also covers this
    point: if the rule were ever written with ``inclusive="right"`` the label would
    not move today (``cond1`` still fires) but would move the moment ``cond1``
    changed, and this is what would notice.
    """
    assert _LOW == 469 / 100000, "this pair does not reach the low threshold in IEEE 754"

    row = _assign(469, 100000)

    assert row["Depth_Score"] == _LOW
    assert row["Confidence"] == _LOW_LABEL


def test_one_ulp_above_the_mid_band_is_still_promoted():
    """The change must not swallow the High tiers -- only the top boundary moves.

    ``5151 / 1000000`` is the smallest ratio above the threshold that a plausible
    integer depth pair produces, and it is still above ``math.nextafter(high, 1)``
    -- the smallest representable step up from the threshold. That makes this the
    tightest available proof that the mid-band's ``inclusive="both"`` upper edge
    does not extend past the threshold itself.
    """
    alt, region = 5151, 1000000
    assert alt / region > math.nextafter(_HIGH, 1)

    row = _assign(alt, region)

    assert row["Confidence"] == _HIGH_STAR_LABEL


@pytest.mark.parametrize(("alt", "region"), EXACT_HIGH_THRESHOLD_PAIRS)
def test_no_row_at_the_boundary_falls_through_to_negative(alt, region):
    """Guards the trap in the second implementation suggested on #184.

    Changing ``cond2``/``cond5`` to ``> high`` *alone* would leave these rows
    matching no condition at all, so they would keep ``NEGATIVE_LABEL`` -- turning
    a reported call into a non-call, which is strictly worse than the bug being
    fixed. The chosen implementation (mid-band applied last, inclusive at both
    ends) cannot do that, and this test is what holds it to that.

    Args:
        alt: ``Estimated_Depth_AlternateVariant`` from :data:`EXACT_HIGH_THRESHOLD_PAIRS`.
        region: The matching ``Estimated_Depth_Variant_ActiveRegion``.
    """
    row = _assign(alt, region)

    assert row["Confidence"] != NEGATIVE_LABEL
    assert row["depth_confidence_pass"]


def test_a_boundary_demotion_changes_which_variant_is_reported():
    """The label is not the only thing that moves -- selection ranks on it.

    ``select_single_best_variant`` sorts by ``Confidence`` first
    (``kestrel_genotyping.py:758``, priority map at :723), so demoting a boundary
    row from ``High_Precision*`` to ``Low_Precision`` hands the reported genotype
    to a different candidate whenever a sample carries several. This is the part of
    #184's blast radius that is not visible in the label alone, and it is why the
    change is gated on the golden cohort rather than treated as cosmetic.

    Before the change the boundary row won on confidence priority 3 against the
    other candidate's 2; after it, it loses with priority 1 and POS 68 is reported
    instead of POS 67.
    """
    frame = kestrel_stage_frame("final", rows=2)
    # Row 0 (POS 67): exactly at the boundary -- High_Precision* before, Low_Precision after.
    frame.loc[0, "Estimated_Depth_AlternateVariant"] = 103
    frame.loc[0, "Estimated_Depth_Variant_ActiveRegion"] = 20000
    # Row 1 (POS 68): alt inside cond5's window with a score well above the boundary,
    # so it is High_Precision both before and after and is unaffected by the change.
    frame.loc[1, "Estimated_Depth_AlternateVariant"] = 50
    frame.loc[1, "Estimated_Depth_Variant_ActiveRegion"] = 1000

    scored = calculate_depth_score_and_assign_confidence(frame, kestrel_config())
    assert scored.loc[0, "Depth_Score"] == _HIGH
    assert scored.loc[0, "Confidence"] == _LOW_LABEL
    assert scored.loc[1, "Confidence"] == _HIGH_LABEL

    best = select_single_best_variant(scored)

    assert len(best) == 1
    assert best["Confidence"].iloc[0] == _HIGH_LABEL
    assert best["POS"].iloc[0] == 68, "the boundary row must no longer win selection"
    assert best.index[0] == 1


def test_cond5_is_unreachable_at_the_boundary_on_integer_depths():
    """Records why only ``High_Precision*`` can move, and proves it arithmetically.

    ``0.00515 == 103 / 20000`` in lowest terms, so any integer pair whose *exact*
    ratio is the threshold has alt as a multiple of 103. The smallest is 103, which
    is already above ``alt_mid_high`` (100), so ``cond5``'s ``[21, 100)`` window can
    never contain a boundary row.

    The exhaustive half of the claim is checked directly against IEEE 754 rather
    than against exact rationals, because the production comparison is a float
    division: for every integer alt at or below ``alt_mid_high`` the only region
    depths whose quotient could round onto the threshold are within one of
    ``alt / high``, and none of them does.

    If this ever fails, the config thresholds changed and the #184 blast-radius
    analysis must be redone from scratch.
    """
    assert Fraction(_HIGH).limit_denominator(10**6) == Fraction(103, 20000)
    assert Fraction(_HIGH).limit_denominator(10**6).numerator > _ALT_MID_HIGH

    for alt in range(1, _ALT_MID_HIGH + 1):
        approx = math.floor(alt / _HIGH)
        for region in (approx - 1, approx, approx + 1):
            if region > 0:
                assert alt / region != _HIGH, (
                    f"alt={alt}, region={region} reaches the boundary inside cond5's window; "
                    "the #184 blast radius is wider than High_Precision* and must be re-derived"
                )


def test_a_fractional_region_depth_reaches_the_boundary_inside_cond5s_window():
    """Characterisation of the scope limit on the #184 blast-radius argument.

    ``calculate_depth_score_and_assign_confidence`` coerces its depth columns with
    ``pd.to_numeric(...).fillna(0)`` and never casts to ``int``, so a fractional
    depth is accepted. With one, ``Depth_Score == high`` *is* reachable at alt
    depths 21-99, where ``cond5`` used to promote to plain ``High_Precision`` --
    a case the "only High_Precision* moves" argument does not cover.

    Nothing in this module enforces integrality; it is a property of the input,
    pinned one layer up by
    :func:`test_production_depths_arrive_as_whole_numbers`. This test exists so
    the limitation is recorded in the suite rather than only in a commit message.
    """
    alt = _ALT_MID_LOW  # 21, the bottom of cond5's window
    region = alt / _HIGH
    assert region != int(region), "this test is only meaningful with a fractional region depth"

    row = _assign(alt, region)

    assert row["Depth_Score"] == _HIGH
    assert row["Confidence"] == _LOW_LABEL


def test_production_depths_arrive_as_whole_numbers():
    """The 103/20000 blast-radius argument depends on this and nothing enforces it.

    Both depth columns are produced by splitting Kestrel's ``Sample`` field
    (``Del:<alt>:<region>``), which carries read counts, so they are integral in
    practice. ``confidence_assignment.py`` itself never casts to ``int``, so if a
    fractional depth ever reached it, ``Depth_Score`` could hit 0.00515 at alt
    depths 21-99 and plain ``High_Precision`` would move as well -- see
    :func:`test_a_fractional_region_depth_reaches_the_boundary_inside_cond5s_window`.

    This pins the property at the boundary where it is actually established: the
    production split, run over the repository's shared notion of a raw Kestrel
    frame.
    """
    frame = split_depth_and_calculate_frame_score(kestrel_stage_frame("raw", rows=3))

    for column in ("Estimated_Depth_AlternateVariant", "Estimated_Depth_Variant_ActiveRegion"):
        values = pd.to_numeric(frame[column], errors="coerce").dropna()
        assert not values.empty, f"{column} produced no usable depths"
        assert (values == values.astype(int)).all(), f"{column} carries a fractional depth"


def _make_df(alt_depth: float, region_depth: float) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [alt_depth],
            "Estimated_Depth_Variant_ActiveRegion": [region_depth],
        }
    )


def test_moving_reporting_floor_alone_leaves_low_precision_band_edge_intact() -> None:
    """Lowering reporting_floor alone admits variants without widening Low_Precision band (#311)."""
    conf = kestrel_config()
    conf["confidence_assignment"]["reporting_floor"] = 0.00200
    # Depth_Score between 0.00200 and 0.00469:
    # Above reporting floor (0.003 >= 0.002), but below Low_Precision band (0.00469).
    # With alt=10, region=3333 -> Depth_Score = 0.00300
    df = _make_df(10, 3333)
    scored = calculate_depth_score_and_assign_confidence(df, conf)
    assert scored.loc[0, "Depth_Score"] < conf["confidence_assignment"]["depth_score_thresholds"]["low"]
    assert scored.loc[0, "Confidence"] == NEGATIVE_LABEL


def test_moving_low_precision_band_edge_alone_leaves_reporting_floor_intact() -> None:
    """Raising depth_score_thresholds.low alone preserves reporting floor guard (#311)."""
    conf = kestrel_config()
    conf["confidence_assignment"]["depth_score_thresholds"]["low"] = 0.00490
    # Variant with Depth_Score 0.00460: below reporting floor (0.00469) -> remains Negative
    df = _make_df(46, 10000)
    scored = calculate_depth_score_and_assign_confidence(df, conf)
    assert scored.loc[0, "Confidence"] == NEGATIVE_LABEL
