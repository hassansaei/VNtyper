"""Pure role-bundle projection for calibration evidence artifacts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from vntyper.scripts.calibration_baseline import project_baseline
from vntyper.scripts.calibration_features import decode_label_artifact


@dataclass(frozen=True)
class RoleInputs:
    """One role's closed feature, label, baseline, and run-hash payloads."""

    features: dict[str, object]
    labels: dict[str, object]
    baseline: dict[str, object]
    run_hashes: dict[str, dict[str, str]]


def build_role_inputs(
    keys: tuple[str, ...],
    feature_rows: Sequence[object],
    label_rows: Sequence[object],
    baseline: Mapping[str, object],
    run_hashes: Mapping[str, Mapping[str, str]],
) -> RoleInputs:
    """Project the exact inputs belonging to a predeclared calibration role."""
    key_set = set(keys)
    features: dict[str, object] = {
        "schema_version": "calibration-features-v1",
        "rows": [row for row in feature_rows if _mapping(row, "feature row").get("manifest_key") in key_set],
    }
    labels: dict[str, object] = {
        "schema_version": "calibration-labels-v1",
        "rows": [row for row in label_rows if _mapping(row, "label row").get("manifest_key") in key_set],
    }
    decoded_labels = decode_label_artifact(labels)
    expected_rows = _sequence(_mapping(baseline.get("expected"), "baseline expected").get("rows"), "expected rows")
    observed_rows = _sequence(_mapping(baseline.get("observed"), "baseline observed").get("rows"), "observed rows")
    selected_expected = [
        row for row in expected_rows if _mapping(row, "baseline expected row").get("manifest_key") in key_set
    ]
    selected_observed = [
        row for row in observed_rows if _mapping(row, "baseline observed row").get("manifest_key") in key_set
    ]
    return RoleInputs(
        features,
        labels,
        project_baseline(
            selected_expected,
            selected_observed,
            {row.manifest_key: row for row in decoded_labels.rows},
        ),
        {key: dict(run_hashes[key]) for key in keys},
    )


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"calibration {label} must be an object")
    return value


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"calibration {label} must be a sequence")
    return value
