"""Strict replay validation for custodian-supplied locked baselines."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence

from vntyper.scripts.calibration_baseline import project_baseline
from vntyper.scripts.calibration_features import FeatureArtifact, LabelArtifact

logger = logging.getLogger(__name__)

_ROW_FIELDS = {
    "manifest_key",
    "order",
    "canonical_identity",
    "name",
    "confidence",
    "flag",
    "tier",
    "support",
    "tie",
    "abstention",
    "identity_projection",
}


def validate_locked_baseline(baseline: Mapping[str, object], features: FeatureArtifact, labels: LabelArtifact) -> None:
    """Require a closed, exactly replayed locked baseline in canonical row order."""
    if set(baseline) != {"schema_version", "expected", "observed"} or baseline["schema_version"] != (
        "calibration-baseline-replay-v1"
    ):
        raise ValueError("locked calibration baseline fields or schema version differ")
    expected = _projection(baseline["expected"], "expected")
    observed = _projection(baseline["observed"], "observed")
    if expected != observed:
        raise ValueError("locked calibration baseline expected and observed replay differ")
    rows = _rows(expected["rows"])
    feature_keys = tuple(row.manifest_key for row in features.rows)
    label_keys = tuple(row.manifest_key for row in labels.rows)
    row_keys = tuple(row["manifest_key"] for row in rows)
    if row_keys != feature_keys or row_keys != label_keys:
        raise ValueError("locked calibration baseline, feature, and label row keys or order differ")
    if tuple(row["order"] for row in rows) != tuple(range(len(rows))):
        raise ValueError("locked calibration baseline row order is not canonical")
    recomputed = project_baseline(rows, rows, {row.manifest_key: row for row in labels.rows}, require_equal=True)
    if recomputed != dict(baseline):
        raise ValueError("locked calibration baseline aggregate, tier, row, or order replay differs")


def _projection(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != {"aggregate", "per_tier", "rows"}:
        raise ValueError(f"locked calibration {label} baseline projection fields differ")
    if not isinstance(value["aggregate"], Mapping) or set(value["aggregate"]) != {
        "displayed",
        "exact",
        "wrong",
        "control_findings",
    }:
        raise ValueError(f"locked calibration {label} baseline aggregate fields differ")
    _counts(value["aggregate"], f"{label} aggregate")
    if not isinstance(value["per_tier"], Mapping):
        raise ValueError(f"locked calibration {label} baseline per-tier value must be an object")
    for tier, counts in value["per_tier"].items():
        if (
            not isinstance(tier, str)
            or not tier
            or not isinstance(counts, Mapping)
            or set(counts)
            != {
                "displayed",
                "exact",
                "wrong",
            }
        ):
            raise ValueError(f"locked calibration {label} baseline per-tier fields differ")
        _counts(counts, f"{label} per-tier")
    return {"aggregate": dict(value["aggregate"]), "per_tier": dict(value["per_tier"]), "rows": value["rows"]}


def _rows(value: object) -> list[dict[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError("locked calibration baseline rows must be a non-empty array")
    rows: list[dict[str, object]] = []
    for raw in value:
        if not isinstance(raw, Mapping) or set(raw) != _ROW_FIELDS:
            raise ValueError("locked calibration baseline row fields differ")
        projection = raw["identity_projection"]
        if not isinstance(projection, Mapping):
            raise ValueError("locked calibration baseline identity projection must be an object")
        if tuple(projection) != tuple(sorted(projection)):
            raise ValueError("locked calibration baseline identity projection order is not canonical")
        for identity, rendered in projection.items():
            if (
                not isinstance(identity, str)
                or not identity
                or not isinstance(rendered, Mapping)
                or set(rendered)
                != {
                    "name",
                    "tier",
                }
            ):
                raise ValueError("locked calibration baseline identity projection fields differ")
            if not isinstance(rendered["name"], str) or not rendered["name"]:
                raise ValueError("locked calibration baseline identity projection name must be non-empty")
            if rendered["tier"] is not None and (not isinstance(rendered["tier"], str) or not rendered["tier"]):
                raise ValueError("locked calibration baseline identity projection tier must be non-empty or null")
        _row_types(raw)
        rows.append(dict(raw))
    return rows


def _counts(value: Mapping[str, object], label: str) -> None:
    if any(isinstance(count, bool) or not isinstance(count, int) or count < 0 for count in value.values()):
        raise ValueError(f"locked calibration baseline {label} counts must be non-negative integers")


def _row_types(row: Mapping[str, object]) -> None:
    if not isinstance(row["manifest_key"], str) or not row["manifest_key"]:
        raise ValueError("locked calibration baseline manifest key must be non-empty")
    if isinstance(row["order"], bool) or not isinstance(row["order"], int) or row["order"] < 0:
        raise ValueError("locked calibration baseline order must be a non-negative integer")
    for field in ("canonical_identity", "name", "confidence", "flag", "tier", "abstention"):
        if row[field] is not None and (not isinstance(row[field], str) or not row[field]):
            raise ValueError(f"locked calibration baseline {field} must be non-empty or null")
    support = row["support"]
    if support is not None and (
        isinstance(support, bool) or not isinstance(support, (int, float)) or not math.isfinite(support) or support < 0
    ):
        raise ValueError("locked calibration baseline support must be finite, non-negative, or null")
    if not isinstance(row["tie"], bool):
        raise ValueError("locked calibration baseline tie must be boolean")
