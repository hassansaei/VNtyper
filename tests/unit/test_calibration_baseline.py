"""Closed row alignment for independently derived shipped baselines."""

from __future__ import annotations

import pytest

from vntyper.scripts.calibration_baseline import project_baseline
from vntyper.scripts.calibration_features import decode_label_artifact

pytestmark = pytest.mark.unit


def _labels():
    artifact = decode_label_artifact(
        {
            "schema_version": "calibration-labels-v1",
            "rows": [
                {
                    "label_key": f"label-{key}",
                    "manifest_key": key,
                    "truth_status": "control",
                    "expected_identity": None,
                    "expected_display_name": None,
                    "mutation_class": "duplication",
                }
                for key in ("a", "b")
            ],
        }
    )
    return {row.manifest_key: row for row in artifact.rows}


def _row(key: str) -> dict[str, object]:
    return {
        "manifest_key": key,
        "order": 0,
        "canonical_identity": None,
        "name": None,
        "confidence": None,
        "flag": None,
        "tier": None,
        "support": None,
        "tie": False,
        "abstention": None,
        "identity_projection": {},
    }


@pytest.mark.parametrize(
    "expected,observed",
    [
        ([_row("a"), _row("a")], [_row("a"), _row("b")]),
        ([_row("a"), _row("b")], [_row("b"), _row("a")]),
        ([_row("a"), _row("b")], [_row("a")]),
    ],
)
def test_projection_rejects_duplicate_reordered_or_missing_row_key_alignment(
    expected: list[dict[str, object]], observed: list[dict[str, object]]
) -> None:
    with pytest.raises(ValueError, match="keys.*align|duplicate|labels"):
        project_baseline(expected, observed, _labels())
