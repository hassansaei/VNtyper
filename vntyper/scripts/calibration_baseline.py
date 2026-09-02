"""Pure independent projection of shipped calibration baselines."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Protocol

from vntyper.scripts.calibration_features import LabelRow


class BaselineExtraction(Protocol):
    """Structural view of independently parsed run projection rows."""

    @property
    def expected_row(self) -> Mapping[str, object]:
        """Return the summary-derived projection row."""
        ...

    @property
    def observed_row(self) -> Mapping[str, object]:
        """Return the fixed-result-file-derived projection row."""
        ...


def build_baseline(
    extractions: Sequence[BaselineExtraction], labels_by_key: Mapping[str, LabelRow]
) -> dict[str, object]:
    """Build independently aggregated expected and observed shipped projections.

    Args:
        extractions: Completed run extractions in canonical member order.
        labels_by_key: Independently decoded truth labels in the same order.

    Returns:
        Closed expected/observed baseline replay document.

    Raises:
        ValueError: If rows, keys, labels, or independent projections differ.
    """
    expected_rows = [dict(_mapping(extraction.expected_row, "expected extraction row")) for extraction in extractions]
    observed_rows = [dict(_mapping(extraction.observed_row, "observed extraction row")) for extraction in extractions]
    return project_baseline(expected_rows, observed_rows, labels_by_key, require_equal=True)


def project_baseline(
    expected_rows: Sequence[object],
    observed_rows: Sequence[object],
    labels_by_key: Mapping[str, LabelRow],
    *,
    require_equal: bool = False,
) -> dict[str, object]:
    """Aggregate independently decoded expected and observed decision rows.

    Args:
        expected_rows: Summary-derived baseline rows.
        observed_rows: Fixed-result-file-derived baseline rows.
        labels_by_key: Independently decoded labels in canonical row order.
        require_equal: Whether all derived aggregate and row values must match.

    Returns:
        Closed expected/observed baseline replay document.

    Raises:
        ValueError: If row keys do not align uniquely and exactly, or required
            independent projections differ.
    """
    expected_values = [dict(_mapping(row, "expected baseline row")) for row in expected_rows]
    observed_values = [dict(_mapping(row, "observed baseline row")) for row in observed_rows]
    expected_keys = _row_keys(expected_values, "expected")
    observed_keys = _row_keys(observed_values, "observed")
    label_keys = tuple(labels_by_key)
    if expected_keys != observed_keys or expected_keys != label_keys:
        raise ValueError("calibration baseline expected, observed, and label row keys must align exactly")
    expected = _projection(expected_values, labels_by_key)
    observed = _projection(observed_values, labels_by_key)
    if expected != observed and require_equal:
        raise ValueError("calibration extraction did not reproduce the shipped decision projection baseline")
    return {"schema_version": "calibration-baseline-replay-v1", "expected": expected, "observed": observed}


def _row_keys(rows: Sequence[Mapping[str, object]], source: str) -> tuple[str, ...]:
    keys = tuple(row.get("manifest_key") for row in rows)
    if any(not isinstance(key, str) or not key for key in keys):
        raise ValueError(f"calibration {source} baseline row keys must be non-empty strings")
    if len(keys) != len(set(keys)):
        raise ValueError(f"calibration {source} baseline row keys must not contain duplicates")
    return keys  # type: ignore[return-value]


def _projection(rows: list[dict[str, object]], labels_by_key: Mapping[str, LabelRow]) -> dict[str, object]:
    aggregate = {"displayed": 0, "exact": 0, "wrong": 0, "control_findings": 0}
    per_tier: dict[str, dict[str, int]] = {}
    for order, row in enumerate(rows):
        row["order"] = order
        key = _row_field(row, "manifest_key")
        if not isinstance(key, str):
            raise ValueError("calibration baseline row manifest key must be a string")
        label = labels_by_key.get(key)
        if label is None:
            raise ValueError(f"calibration baseline row has no independent label: {key}")
        name = _row_field(row, "name")
        identity = _row_field(row, "canonical_identity")
        tier = _row_field(row, "tier")
        if name is not None:
            aggregate["displayed"] += 1
            if label.truth_status == "control":
                aggregate["control_findings"] += 1
            if name != label.expected_display_name:
                aggregate["wrong"] += 1
        if identity is not None and identity == label.expected_identity:
            aggregate["exact"] += 1
        if isinstance(tier, str):
            counts = per_tier.setdefault(tier, {"displayed": 0, "exact": 0, "wrong": 0})
            if name is not None:
                counts["displayed"] += 1
                if name != label.expected_display_name:
                    counts["wrong"] += 1
            if identity is not None and identity == label.expected_identity:
                counts["exact"] += 1
    return {"aggregate": aggregate, "per_tier": per_tier, "rows": rows}


def _row_field(row: Mapping[str, object], key: str) -> object:
    if key not in row:
        raise ValueError(f"calibration baseline row lacks required field: {key}")
    return row[key]


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"calibration {label} must be an object")
    return value
