"""Round-trip tests for the packaged decision-profile renderer."""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.render_decision_profile import build_decision_profile, build_decision_projection, render
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile_schema import validate_complete_inventory

pytestmark = pytest.mark.unit


def test_renderer_round_trip_matches_the_committed_artifacts(tmp_path: Path) -> None:
    render(tmp_path)
    committed = Path(__file__).parents[2] / "vntyper" / "profiles"

    for name in ("decision_projection.json", "decision_profile.json", "decision_profile.sha256"):
        assert (tmp_path / name).read_bytes() == (committed / name).read_bytes()


def test_builder_produces_a_valid_canonical_profile() -> None:
    projection = build_decision_projection()
    profile = build_decision_profile(projection)

    validate_complete_inventory(profile)
    assert load_strict_json_object(canonical_json_bytes(profile)) == profile
    assert len(canonical_sha256(profile)) == 64
