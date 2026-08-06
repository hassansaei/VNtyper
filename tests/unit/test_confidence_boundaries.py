"""Boundary characterisation for ``confidence_assignment``.

Why this file exists
--------------------
``calculate_depth_score_and_assign_confidence`` produces the ``Confidence``
string a clinician reads. Before this file it had 100% line coverage and a 21%
mutation score: inverting the central condition (``cond1``'s ``==`` to ``!=``)
left the whole suite green, and five of its six threshold comparisons could be
changed without a single test failing. Line coverage says the code ran; it says
nothing about whether the *boundary* was pinned.

Everything here is **characterisation**. It records what the code does today so
that a future change to a threshold, a comparison operator or the order the
conditions are applied in becomes visible. It is not a specification, and none
of it should be read as a claim that the current behaviour is correct.

What is pinned
--------------
* A 54-cell matrix over the cross product of the three alt-depth thresholds
  (each probed at threshold-1, threshold and threshold+1) and the two
  depth-score thresholds (each probed just below, exactly at and just above).
  Every threshold is read from the shipped config through
  ``tests.builders.kestrel_config``; none is hardcoded.
* The region-depth tier interaction, whose precedence is **unspecified** -- see
  ``test_region_depth_demotion_is_overwritten_by_a_later_high_precision_tier``.

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

import pandas as pd
import pytest

from tests.builders import kestrel_config
from vntyper.scripts.confidence_assignment import (
    NEGATIVE_LABEL,
    calculate_depth_score_and_assign_confidence,
)

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
#                  alt_low-1  alt_low  alt_low+1  mid_low-1  mid_low  mid_low+1  mid_high-1  mid_high  mid_high+1
_EXPECTED_MATRIX = {
    "low-eps": (_N, _N, _N, _N, _N, _N, _N, _N, _N),
    "low": (_L, _L, _L, _L, _L, _L, _L, _L, _L),
    "low+eps": (_L, _L, _L, _L, _L, _L, _L, _L, _L),
    "high-eps": (_L, _L, _L, _L, _L, _L, _L, _L, _L),
    "high": (_L, _L, _H, _L, _H, _H, _H, _S, _S),
    "high+eps": (_L, _L, _H, _L, _H, _H, _H, _S, _S),
}

# `Depth_Score` is computed by production as `alt / region`, so each cell supplies
# a region depth of `alt / target`. That round trip lands on the target exactly for
# 53 of the 54 cells. The exception is alt=99 at Depth_Score == high: no double `r`
# satisfies `99 / r == 0.00515`, because one ULP of `r` near 19223 is wider than the
# pre-image interval of the target. That cell lands one ULP above `high`, which is
# on the same side of every comparison the function makes there, so its expected
# label is unchanged. The row assertion below is written as `>= high` for that
# reason and only for that reason.
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
    """Characterisation: a low region depth does **not** cap the confidence label.

    ``cond1`` demotes a row whose ``Estimated_Depth_Variant_ActiveRegion`` is at or
    below ``var_active_region_threshold`` to Low_Precision, but it is applied
    *first*; ``cond2`` (High_Precision*) and ``cond5`` (High_Precision) are applied
    after it and overwrite the demotion outright.

    **Tier precedence is unspecified.** Nothing in the code states whether an
    earlier demotion or a later promotion should win, and the six conditions are
    simply applied in source order. This test records the order that ships today.
    It is characterisation, not specification: do not read it as a claim that a
    High_Precision call on a 150-read active region is correct.

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
