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


@pytest.mark.parametrize("missing", ["name", "canonical_identity", "tier", "manifest_key"])
def test_projection_reports_a_missing_row_field_as_a_typed_error(missing: str) -> None:
    """Mutation caught: a truncated baseline row escapes as a bare ``KeyError``."""
    expected = [_row("a"), _row("b")]
    observed = [_row("a"), _row("b")]
    del expected[1][missing]

    with pytest.raises(ValueError, match=f"lacks required field: {missing}|row keys must be non-empty strings"):
        project_baseline(expected, observed, _labels())


def test_projection_counts_control_findings_for_representation_limited_call() -> None:
    """Control findings are counted whenever a control row has a canonical identity, even if unrepresentable."""
    row_a = _row("a")
    row_a["canonical_identity"] = "c.54_56delinsAT"
    row_a["name"] = "frameshift +1, representation-limited"
    row_a["tier"] = "A"

    row_b = _row("b")
    row_b["canonical_identity"] = None
    row_b["name"] = None

    expected = [row_a, row_b]
    observed = [dict(row_a), dict(row_b)]

    projected = project_baseline(expected, observed, _labels())
    expected_data = projected["expected"]
    assert isinstance(expected_data, dict)
    aggregate = expected_data["aggregate"]
    assert isinstance(aggregate, dict)

    assert aggregate["control_findings"] == 1
    assert aggregate["displayed"] == 0
    assert aggregate["wrong"] == 0
