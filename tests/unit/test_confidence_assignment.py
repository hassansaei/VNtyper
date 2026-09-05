# tests/unit/test_confidence_assignment.py

"""
Unit tests for confidence assignment functionality.
Validates confidence levels based on depth scores and thresholds
from vntyper/scripts/kestrel_config.json.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from tests.builders import kestrel_config as shipped_kestrel_config
from vntyper.scripts.confidence_assignment import (
    NEGATIVE_LABEL,
    calculate_depth_score_and_assign_confidence,
)

pytestmark = pytest.mark.unit

# Fixed region depth used across tests to isolate Depth_Score effects.
# Must be above var_active_region_threshold (200) to avoid triggering cond1's region check.
_REGION_DEPTH = 10000


@pytest.fixture(scope="session")
def kestrel_config():
    """Load the Kestrel configuration from vntyper/scripts/kestrel_config.json."""
    this_file = Path(__file__).resolve()
    config_path = this_file.parents[2] / "vntyper" / "scripts" / "kestrel_config.json"
    if not config_path.exists():
        pytest.exit(f"kestrel_config.json not found at {config_path}", returncode=1)
    with config_path.open("r") as f:
        return json.load(f)


def _make_df(alt_depth: float, region_depth: float = _REGION_DEPTH) -> pd.DataFrame:
    """Helper to create a single-row DataFrame with the required depth columns."""
    return pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [alt_depth],
            "Estimated_Depth_Variant_ActiveRegion": [region_depth],
        }
    )


def test_calculate_depth_score_empty_df(kestrel_config):
    """Empty input should yield empty output."""
    df = pd.DataFrame()
    out = calculate_depth_score_and_assign_confidence(df, kestrel_config)
    assert out.empty


def test_depth_score_below_threshold_is_negative(kestrel_config):
    """Variants with Depth_Score < low_threshold get Negative (filtered out)."""
    low = kestrel_config["confidence_assignment"]["depth_score_thresholds"]["low"]
    # Derive alt_depth that produces Depth_Score just below low_threshold
    alt_depth = low * _REGION_DEPTH * 0.5  # half the threshold
    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == NEGATIVE_LABEL
    assert not out.loc[0, "depth_confidence_pass"]


def test_depth_score_at_threshold_boundary_not_negative(kestrel_config):
    """Variants with Depth_Score == low_threshold should NOT be Negative."""
    low = kestrel_config["confidence_assignment"]["depth_score_thresholds"]["low"]
    # Derive alt_depth that produces Depth_Score exactly at low_threshold
    alt_depth = low * _REGION_DEPTH
    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth), kestrel_config)
    assert out.loc[0, "Confidence"] != NEGATIVE_LABEL, (
        "Depth_Score at exactly low_threshold must not be Negative; "
        "the Negative filter applies only to strictly lower scores."
    )
    assert out.loc[0, "depth_confidence_pass"]


def test_depth_score_in_low_precision_band(kestrel_config):
    """Variants with Depth_Score in (low, high) band get Low_Precision."""
    conf = kestrel_config["confidence_assignment"]
    low = conf["depth_score_thresholds"]["low"]
    high = conf["depth_score_thresholds"]["high"]
    mid_low = conf["alt_depth_thresholds"]["mid_low"]
    mid_high = conf["alt_depth_thresholds"]["mid_high"]
    low_label = conf["confidence_levels"]["low_precision"]

    # Target Depth_Score in (low, high) and alt in [mid_low, mid_high] => cond3 Low_Precision
    depth_score = (low + high) / 2
    alt_depth = (mid_low + mid_high) // 2
    region_depth = int(alt_depth / depth_score) + 1

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == low_label


def test_depth_score_high_precision_star(kestrel_config):
    """High alt depth + high Depth_Score => High_Precision*."""
    conf = kestrel_config["confidence_assignment"]
    high = conf["depth_score_thresholds"]["high"]
    alt_mid_high = conf["alt_depth_thresholds"]["mid_high"]
    high_star_label = conf["confidence_levels"]["high_precision_star"]

    # alt >= mid_high AND Depth_Score >= high_threshold => cond2
    alt_depth = alt_mid_high
    region_depth = int(alt_depth / high)  # produces Depth_Score >= high

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == high_star_label


def test_low_region_depth_with_sufficient_depth_score(kestrel_config):
    """Low region depth (below var_active_region_threshold) should get Low_Precision, not Negative."""
    conf = kestrel_config["confidence_assignment"]
    low = conf["depth_score_thresholds"]["low"]
    var_region_threshold = conf["var_active_region_threshold"]
    low_label = conf["confidence_levels"]["low_precision"]

    # Region depth at threshold, but Depth_Score above low_threshold
    region_depth = var_region_threshold
    alt_depth = low * region_depth * 2  # Depth_Score = 2 * low_threshold (well above)

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == low_label


def test_low_region_depth_with_insufficient_depth_score(kestrel_config):
    """Low region depth AND Depth_Score below threshold should still be Negative."""
    conf = kestrel_config["confidence_assignment"]
    low = conf["depth_score_thresholds"]["low"]
    var_region_threshold = conf["var_active_region_threshold"]

    # Both below thresholds
    region_depth = var_region_threshold
    alt_depth = low * region_depth * 0.5  # Depth_Score = 0.5 * low_threshold

    out = calculate_depth_score_and_assign_confidence(_make_df(alt_depth, region_depth), kestrel_config)
    assert out.loc[0, "Confidence"] == NEGATIVE_LABEL
    assert not out.loc[0, "depth_confidence_pass"]


def test_both_depths_unparseable_stays_negative_with_a_nan_score(kestrel_config):
    """Unparseable depths coerce to 0, so 0/0 is NaN and the row is Negative.

    Changing the coercion default from ``fillna(0)`` to ``fillna(1)`` instead
    manufactures a 1.0 score and a false scored call.
    """
    df = pd.DataFrame(
        [
            {
                "Estimated_Depth_AlternateVariant": "not_a_number",
                "Estimated_Depth_Variant_ActiveRegion": "also_unparseable",
            }
        ]
    )

    row = calculate_depth_score_and_assign_confidence(df, kestrel_config).iloc[0]

    assert row["Estimated_Depth_AlternateVariant"] == 0.0
    assert row["Estimated_Depth_Variant_ActiveRegion"] == 0.0
    assert pd.isna(row["Depth_Score"]), "0/0 must stay NaN, never a manufactured score"
    assert row["Confidence"] == NEGATIVE_LABEL
    assert not row["depth_confidence_pass"]


def test_an_unparseable_alt_depth_scores_zero_and_stays_negative(kestrel_config):
    """An unparseable alternate depth yields a zero score, not a scored call."""
    df = pd.DataFrame(
        [
            {
                "Estimated_Depth_AlternateVariant": "corrupt",
                "Estimated_Depth_Variant_ActiveRegion": _REGION_DEPTH,
            }
        ]
    )

    row = calculate_depth_score_and_assign_confidence(df, kestrel_config).iloc[0]

    assert row["Estimated_Depth_AlternateVariant"] == 0.0
    assert row["Estimated_Depth_Variant_ActiveRegion"] == _REGION_DEPTH
    assert row["Depth_Score"] == 0.0
    assert row["Confidence"] == NEGATIVE_LABEL
    assert not row["depth_confidence_pass"]


def test_an_unparseable_region_depth_keeps_a_parseable_alt_negative(kestrel_config):
    """A malformed region depth must not turn a parseable 120 alt depth into a call.

    At HEAD, 120/0 becomes a NaN score and remains Negative. The ``fillna(1)``
    mutant turns it into 120/1, a score of 120 and a false High_Precision*.
    """
    df = pd.DataFrame(
        [
            {
                "Estimated_Depth_AlternateVariant": 120,
                "Estimated_Depth_Variant_ActiveRegion": "corrupt",
            }
        ]
    )

    row = calculate_depth_score_and_assign_confidence(df, kestrel_config).iloc[0]

    assert row["Estimated_Depth_AlternateVariant"] == 120.0
    assert row["Estimated_Depth_Variant_ActiveRegion"] == 0.0
    assert pd.isna(row["Depth_Score"]), "120/0 must stay NaN, never become 120"
    assert row["Confidence"] == NEGATIVE_LABEL
    assert not row["depth_confidence_pass"]


# ---------------------------------------------------------------------------
# The calibration constants are required, not defaulted.
# ---------------------------------------------------------------------------

#: The six calibration constants, as dotted paths inside ``confidence_assignment``,
#: paired with the ``.get()`` default that used to stand in for a missing key.
#:
#: Those defaults were not neutral placeholders -- they encoded a SECOND, WRONG
#: calibration. ``depth_score_thresholds.high`` defaulted to 0.4 where the shipped
#: config supplies 0.00515, a factor of 78. Dropping a single key was enough to
#: genotype a whole cohort against thresholds nobody chose, silently and with no
#: log line. ``--config-path`` replaces the whole config rather than merging it
#: (AGENTS.md trap 2) and a partial config already fails with ``KeyError``
#: elsewhere in the pipeline, so a partial config is not supported input and the
#: right behaviour is to fail loudly here too.
CALIBRATION_CONSTANTS = {
    "reporting_floor": (0.00469, 0.0),
    "depth_score_thresholds.low": (0.00469, 0.2),
    "depth_score_thresholds.high": (0.00515, 0.4),
    "alt_depth_thresholds.low": (20, 5),
    "alt_depth_thresholds.mid_low": (21, 10),
    "alt_depth_thresholds.mid_high": (100, 20),
    "var_active_region_threshold": (200, 0),
}


def test_the_shipped_config_supplies_every_calibration_constant():
    """Requiring the keys only leaves behaviour unchanged because all six are shipped.

    This is the whole safety argument for removing the defaults, so it is asserted
    rather than assumed. It also pins the shipped values against the defaults they
    replaced: a config edit moving one of them back towards its old fallback would
    be a recalibration, and would surface here.
    """
    conf = shipped_kestrel_config()["confidence_assignment"]

    for dotted, (shipped, old_default) in CALIBRATION_CONSTANTS.items():
        node = conf
        parts = dotted.split(".")
        for part in parts[:-1]:
            assert part in node, f"{dotted}: the shipped config is missing the {part!r} section"
            node = node[part]
        assert parts[-1] in node, f"{dotted}: the shipped config is missing this key"
        assert node[parts[-1]] == shipped, f"{dotted}: shipped value moved"
        assert node[parts[-1]] != old_default, f"{dotted}: the shipped value now equals the old fallback"


@pytest.mark.parametrize("dotted", sorted(CALIBRATION_CONSTANTS))
def test_dropping_a_calibration_constant_raises_instead_of_recalibrating(dotted):
    """A missing calibration key must raise, not fall back to a different calibration.

    Before this change every one of these was read as ``.get(key, <default>)``, so a
    config missing ``depth_score_thresholds.high`` scored against 0.4 instead of
    0.00515 and turned an ample-depth ``High_Precision*`` call into
    ``Low_Precision`` -- a changed genotype, produced silently, from a config
    nothing validates. Labels under ``confidence_levels`` deliberately keep their
    fallbacks: those are cosmetic strings, not calibration.

    Args:
        dotted: Dotted path of the key to remove, from
            :data:`CALIBRATION_CONSTANTS`.
    """
    config = shipped_kestrel_config()
    node = config["confidence_assignment"]
    parts = dotted.split(".")
    for part in parts[:-1]:
        node = node[part]
    del node[parts[-1]]

    with pytest.raises(KeyError, match=parts[-1]):
        calculate_depth_score_and_assign_confidence(_make_df(120, 12000), config)


def test_an_empty_frame_still_short_circuits_before_the_config_is_read():
    """The empty-input fast path must not need a config at all.

    ``calculate_depth_score_and_assign_confidence`` returns before touching
    ``kestrel_config``, and callers rely on that: a run with no Kestrel variants
    reaches this function with an empty frame. Requiring the calibration keys must
    not turn that into a ``KeyError``.
    """
    out = calculate_depth_score_and_assign_confidence(pd.DataFrame(), {})

    assert out.empty
