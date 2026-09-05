"""Cross-file invariant: calibration.json must agree with kestrel_config.json.

This is provenance, not validation. It catches a threshold changed in one file
and not the other, and it forces every calibrated constant to carry a stated
source and rationale. It does NOT prove a value is scientifically correct --
the behavioural tests do that.
"""

import json
from pathlib import Path
from typing import Any

import pytest

pytestmark = pytest.mark.unit

SCRIPTS = Path(__file__).resolve().parents[2] / "vntyper" / "scripts"


@pytest.fixture(scope="module")
def calibration() -> dict:
    return json.loads((SCRIPTS / "calibration.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def kestrel_config() -> dict:
    return json.loads((SCRIPTS / "kestrel_config.json").read_text(encoding="utf-8"))


def _resolve(config: dict, dotted: str) -> Any:
    node: Any = config
    for part in dotted.split("."):
        node = node[part]
    return node


def test_every_calibrated_value_matches_the_config(calibration, kestrel_config) -> None:
    """The manifest and the live config must not drift apart."""
    mismatches = []
    for dotted, entry in calibration["constants"].items():
        actual = _resolve(kestrel_config, dotted)
        if actual != entry["value"]:
            mismatches.append(f"{dotted}: config={actual!r} manifest={entry['value']!r}")
    assert not mismatches, "calibration.json disagrees with kestrel_config.json: " + "; ".join(mismatches)


def test_every_entry_has_provenance(calibration) -> None:
    """A calibrated number without a stated source and rationale is unreviewable."""
    incomplete = [
        dotted
        for dotted, entry in calibration["constants"].items()
        if not entry.get("source") or not entry.get("rationale")
    ]
    assert not incomplete, f"These constants lack source or rationale: {incomplete}"


def test_every_value_is_within_its_declared_bounds(calibration) -> None:
    """Bounds encode a reviewed safe range; a value outside it is a red flag."""
    violations = []
    for dotted, entry in calibration["constants"].items():
        bounds = entry.get("bounds")
        if bounds is None:
            continue
        low, high = bounds
        if not low <= entry["value"] <= high:
            violations.append(f"{dotted}={entry['value']} outside {bounds}")
    assert not violations, "; ".join(violations)


def test_linked_constants_hold_the_same_value(calibration, kestrel_config) -> None:
    """Constants declared as aliases of one another must not drift.

    ``alt_filtering.gg_depth_score_threshold`` and
    ``confidence_assignment.reporting_floor`` are the same calibrated
    number stored twice; nothing in the config prevents them diverging.
    """
    for group in calibration.get("linked", []):
        values = {dotted: _resolve(kestrel_config, dotted) for dotted in group}
        assert len(set(values.values())) == 1, f"Linked constants have drifted apart: {values}"


def test_state_space_limits_are_equal(kestrel_config) -> None:
    """max_align_states and max_hap_states are always set together."""
    settings = kestrel_config["kestrel_settings"]
    assert settings["max_align_states"] == settings["max_hap_states"]


def test_threshold_ordering(kestrel_config) -> None:
    """Ordering invariants that make the confidence bands coherent."""
    conf = kestrel_config["confidence_assignment"]
    assert conf["depth_score_thresholds"]["low"] < conf["depth_score_thresholds"]["high"]
    alt = conf["alt_depth_thresholds"]
    assert alt["low"] < alt["mid_low"] < alt["mid_high"]
