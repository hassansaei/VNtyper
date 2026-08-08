"""Pure config-path coverage for the in-container image probe."""

from __future__ import annotations

import runpy
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

PROBE = Path(__file__).parents[1] / "docker" / "image_probe.py"


def _collect_references(config: dict) -> dict:
    """Load the probe definitions without running its container-only ``main``."""
    namespace = runpy.run_path(str(PROBE), run_name="image_probe_test")
    return namespace["_collect_references"]("/absent/package", config)


def test_unset_reference_candidates_are_not_treated_as_paths(tmp_path: Path) -> None:
    """A shipped null CRAM candidate means not configured, not a filesystem path."""
    reference = tmp_path / "configured.fa"
    reference.write_text(">chr1\nA\n", encoding="utf-8")

    references = _collect_references(
        {"reference_data": {"configured": str(reference), "cram_reference_hg38": None}},
    )

    assert set(references) == {"configured"}
    assert references["configured"]["exists"] is True


def test_non_null_non_string_reference_candidate_fails_closed() -> None:
    """Malformed config cannot silently disappear from image structure validation."""
    with pytest.raises(ValueError, match="reference_data.bad must be a path string or null"):
        _collect_references({"reference_data": {"bad": 42}})
