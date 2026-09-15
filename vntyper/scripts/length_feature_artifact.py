"""Closed feature bundles carrying the full measurement provenance needed for reuse."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from typing import NoReturn, cast

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_feature_provenance import (
    decode_length_feature_provenance,
    encode_length_feature_provenance,
)
from vntyper.scripts.length_features import (
    LENGTH_FEATURES_SCHEMA_VERSION,
    REGION_NAMES,
    LengthFeatures,
    LengthFeatureStatus,
    RegionFeatures,
    encode_length_features,
)

logger = logging.getLogger(__name__)
_SCHEMA = "length-feature-artifact-v1"
_ROW_FIELDS = {
    "manifest_key",
    "assembly",
    "assay_class",
    "input_scope",
    "annotation_sha256",
    "counting_policy_sha256",
    "provenance_sha256",
    "regions",
    "A",
    "F",
    "status",
    "reasons",
}
_REGION_FIELDS = {"length_bp", "depth_sum", "mean_depth", "covered_fraction", "supporting_fragment_count"}


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"length feature artifact {label} fields differ from the closed contract")
    return value


def _region(value: object) -> RegionFeatures | None:
    if value is None:
        return None
    raw = _object(value, _REGION_FIELDS, "region")
    # The shared encoder checks numeric types, ranges, derived means and support.
    return RegionFeatures(
        cast(int, raw["length_bp"]),
        cast(float, raw["depth_sum"]),
        cast(float, raw["mean_depth"]),
        cast(float, raw["covered_fraction"]),
        cast(int | None, raw["supporting_fragment_count"]),
    )


def _ratio(value: object) -> float | None:
    if value is not None and (
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
    ):
        _fail("length feature artifact ratios must be finite numbers or null")
    return cast(float | None, value)


def decode_length_feature_artifact(value: object) -> LengthFeatures:
    """Restore one exact feature row and revalidate its full provenance binding.

    Args:
        value: Parsed closed feature bundle, with both feature and provenance bodies.

    Returns:
        Immutable features with their original feature and provenance digests.

    Raises:
        ValueError: If schema, fields, types, derived values or provenance differ.
    """
    raw = _object(value, {"schema_version", "features", "provenance"}, "bundle")
    if raw["schema_version"] != _SCHEMA:
        _fail("length feature artifact schema version is unsupported")
    document = _object(raw["features"], {"schema_version", "rows"}, "features")
    if document["schema_version"] != LENGTH_FEATURES_SCHEMA_VERSION:
        _fail("length feature artifact feature schema version is unsupported")
    rows = document["rows"]
    if not isinstance(rows, list) or len(rows) != 1:
        _fail("length feature artifact requires exactly one feature row")
    row = _object(rows[0], _ROW_FIELDS, "row")
    regions = _object(row["regions"], set(REGION_NAMES), "regions")
    reasons, status = row["reasons"], row["status"]
    if not isinstance(reasons, list) or any(not isinstance(reason, str) for reason in reasons):
        _fail("length feature artifact reasons must be an array of strings")
    if not isinstance(status, str) or status not in {"measured", "partial", "unavailable"}:
        _fail("length feature artifact status is unsupported")
    provenance = decode_length_feature_provenance(raw["provenance"])
    measured = LengthFeatures(
        manifest_key=cast(str, row["manifest_key"]),
        assembly=cast(str, row["assembly"]),
        assay_class=cast(str, row["assay_class"]),
        input_scope=cast(str, row["input_scope"]),
        annotation_sha256=cast(str, row["annotation_sha256"]),
        counting_policy_sha256=cast(str, row["counting_policy_sha256"]),
        provenance_sha256=cast(str, row["provenance_sha256"]),
        regions={name: _region(regions[name]) for name in REGION_NAMES},
        a=_ratio(row["A"]),
        f=_ratio(row["F"]),
        status=cast(LengthFeatureStatus, status),
        reasons=tuple(reasons),
        provenance=provenance,
        sha256=canonical_sha256(document),
    )
    encode_length_features(measured)
    return measured


def encode_length_feature_artifact(features: LengthFeatures) -> dict[str, object]:
    """Bundle a validated feature row with its complete measurement provenance.

    Args:
        features: Existing immutable measurement with bound provenance.

    Returns:
        New closed JSON-compatible bundle, without changing either component digest.

    Raises:
        ValueError: If typed features or their provenance binding are invalid.
    """
    return {
        "schema_version": _SCHEMA,
        "features": encode_length_features(features),
        "provenance": encode_length_feature_provenance(features.provenance),
    }
