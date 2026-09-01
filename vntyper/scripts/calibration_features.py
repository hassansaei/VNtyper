"""Allowlisted calibration features and separately keyed truth labels."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.calibration_manifest import PartitionManifest
from vntyper.scripts.canonical_json import canonical_sha256

TruthStatus = Literal["mutated", "control"]

APPROVED_RUNTIME_FEATURES = frozenset(
    {
        "canonical_identity",
        "motif_context",
        "pair_context",
        "context_local_representation_count",
        "alternate_kmer_path_depth",
        "active_region_depth",
        "depth_score",
        "gdp",
        "dp",
        "structural_gate",
        "cooccurring_identity_count",
        "haplotype_record_count",
        "haplotype_record_share",
        "haplotype_record_count_margin",
        "haplotype_record_share_margin",
        "haplotype_record_tie",
        "xd_availability_count",
        "xd_availability_fraction",
        "xd_median",
        "xd_interquartile_range",
        "advntr_evidence_disposition",
        "advntr_sequencing_read_support",
        "advntr_p_value",
        "advntr_coverage",
        "tool_version",
        "reference_version",
        "assay_class",
    }
)
_RECIPE_FIELDS = {
    "transform_fit_role",
    "candidate_comparison_role",
    "missing_value_fit_role",
    "cap_fit_role",
    "feature_selection_fit_role",
}


@dataclass(frozen=True)
class FeatureRow:
    """One opaque-keyed runtime feature row without truth or identifiers."""

    feature_key: str
    manifest_key: str
    features: Mapping[str, object]


@dataclass(frozen=True)
class FeatureArtifact:
    """Strict feature rows stored independently from labels."""

    rows: tuple[FeatureRow, ...]
    sha256: str


@dataclass(frozen=True)
class LabelRow:
    """One independently established truth row in a separate artifact."""

    label_key: str
    manifest_key: str
    truth_status: TruthStatus
    expected_identity: str | None
    expected_display_name: str | None
    mutation_class: str


@dataclass(frozen=True)
class LabelArtifact:
    """Strict separately keyed truth labels."""

    rows: tuple[LabelRow, ...]
    sha256: str


@dataclass(frozen=True)
class TrainingRecipe:
    """Roles authorized to fit transforms and compare frozen candidates."""

    transform_fit_role: str
    candidate_comparison_role: str
    missing_value_fit_role: str
    cap_fit_role: str
    feature_selection_fit_role: str


def decode_feature_artifact(value: object) -> FeatureArtifact:
    """Decode feature rows and reject every field outside the runtime allowlist."""
    root = _artifact_root(value, "calibration-features-v1", "calibration feature artifact")
    raw_rows = root["rows"]
    assert isinstance(raw_rows, list)
    rows = tuple(_decode_feature_row(raw) for raw in raw_rows)
    _unique_increasing(tuple(row.feature_key for row in rows), "calibration feature keys")
    _unique_increasing(tuple(row.manifest_key for row in rows), "calibration feature manifest keys")
    return FeatureArtifact(rows, canonical_sha256(root))


def decode_label_artifact(value: object) -> LabelArtifact:
    """Decode separately keyed truth labels under a closed schema."""
    root = _artifact_root(value, "calibration-labels-v1", "calibration label artifact")
    raw_rows = root["rows"]
    assert isinstance(raw_rows, list)
    rows = tuple(_decode_label_row(raw) for raw in raw_rows)
    _unique_increasing(tuple(row.label_key for row in rows), "calibration label keys")
    _unique_increasing(tuple(row.manifest_key for row in rows), "calibration label manifest keys")
    return LabelArtifact(rows, canonical_sha256(root))


def validate_training_recipe(value: object) -> TrainingRecipe:
    """Require transforms from training and comparisons from policy selection."""
    raw = _exact_object(value, _RECIPE_FIELDS, "calibration training recipe")
    training_fields = (
        "transform_fit_role",
        "missing_value_fit_role",
        "cap_fit_role",
        "feature_selection_fit_role",
    )
    for field in training_fields:
        if raw[field] != "training":
            raise ValueError(f"calibration {field} must use training only")
    if raw["candidate_comparison_role"] != "policy-selection":
        raise ValueError("calibration candidates must be compared on policy-selection only")
    return TrainingRecipe(
        "training",
        "policy-selection",
        "training",
        "training",
        "training",
    )


def validate_artifact_alignment(
    features: FeatureArtifact,
    labels: LabelArtifact,
    manifest: PartitionManifest,
) -> None:
    """Require one separate feature and label row for every manifest member."""
    if not isinstance(features, FeatureArtifact) or not isinstance(labels, LabelArtifact):
        raise ValueError("calibration feature and label artifacts must be decoded values")
    if not isinstance(manifest, PartitionManifest):
        raise ValueError("calibration partition manifest must be a PartitionManifest")
    manifest_keys = tuple(member.key for member in manifest.members)
    feature_manifest_keys = tuple(row.manifest_key for row in features.rows)
    label_manifest_keys = tuple(row.manifest_key for row in labels.rows)
    if feature_manifest_keys != manifest_keys or label_manifest_keys != manifest_keys:
        raise ValueError("calibration feature and label rows must align exactly with the partition manifest")
    if {row.feature_key for row in features.rows} & {row.label_key for row in labels.rows}:
        raise ValueError("calibration feature and label keys must use separate namespaces")


def _artifact_root(value: object, version: str, label: str) -> Mapping[str, object]:
    root = _exact_object(value, {"schema_version", "rows"}, label)
    if root["schema_version"] != version:
        raise ValueError(f"{label} schema version must be {version}")
    if not isinstance(root["rows"], list) or not root["rows"]:
        raise ValueError(f"{label} rows must be a non-empty list")
    return root


def _decode_feature_row(value: object) -> FeatureRow:
    raw = _exact_object(value, {"feature_key", "manifest_key", "features"}, "calibration feature row")
    feature_key = _text(raw["feature_key"], "calibration feature key")
    manifest_key = _text(raw["manifest_key"], "calibration feature manifest key")
    features = raw["features"]
    if not isinstance(features, Mapping) or not features:
        raise ValueError("calibration features must be a non-empty object")
    unknown = set(features) - APPROVED_RUNTIME_FEATURES
    if unknown:
        raise ValueError(f"calibration feature fields are not allowlisted: {sorted(unknown)}")
    parsed = {str(name): _feature_value(str(name), item) for name, item in features.items()}
    return FeatureRow(feature_key, manifest_key, MappingProxyType(parsed))


def _decode_label_row(value: object) -> LabelRow:
    raw = _exact_object(
        value,
        {
            "label_key",
            "manifest_key",
            "truth_status",
            "expected_identity",
            "expected_display_name",
            "mutation_class",
        },
        "calibration label row",
    )
    status = raw["truth_status"]
    if status not in {"mutated", "control"}:
        raise ValueError(f"unsupported calibration truth status: {status!r}")
    identity = _optional_text(raw["expected_identity"], "expected identity")
    display_name = _optional_text(raw["expected_display_name"], "expected display name")
    mutation_class = _text(raw["mutation_class"], "calibration mutation class")
    if status == "mutated" and (identity is None or display_name is None):
        raise ValueError("mutated calibration labels require expected identity and displayed name")
    if status == "control" and (identity is not None or display_name is not None):
        raise ValueError("control calibration labels require null expected identity and displayed name")
    return LabelRow(
        _text(raw["label_key"], "calibration label key"),
        _text(raw["manifest_key"], "calibration label manifest key"),
        cast(TruthStatus, status),
        identity,
        display_name,
        mutation_class,
    )


def _feature_value(name: str, value: object) -> object:
    if isinstance(value, (str, bool)) or value is None:
        return value
    if isinstance(value, int):
        if value < 0:
            raise ValueError(f"calibration feature {name} cannot be negative")
        return value
    if isinstance(value, float):
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"calibration feature {name} must be finite and non-negative")
        if name == "xd_availability_fraction" and value > 1:
            raise ValueError("calibration feature XD availability fraction cannot exceed one")
        return value
    raise ValueError(f"calibration feature {name} has an unsupported value type")


def _unique_increasing(values: tuple[str, ...], label: str) -> None:
    if values != tuple(sorted(values)) or len(values) != len(set(values)):
        raise ValueError(f"{label} must be unique and increasing")


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a non-empty string")
    return value


def _optional_text(value: object, label: str) -> str | None:
    if value is None:
        return None
    return _text(value, label)


def _exact_object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    return value
