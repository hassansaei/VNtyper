"""Pure depth summaries and ratios for VNTR total-length research."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_annotation import (
    Interval,
    LengthAnnotation,
    decode_length_annotation,
    encode_length_annotation,
)
from vntyper.scripts.length_feature_provenance import (
    DENOMINATOR_NAMES,
    DenominatorQc,
    LengthFeatureContext,
    LengthFeatureProvenance,
    QueryInterval,
    bind_length_feature_provenance,
    encode_length_feature_context,
    encode_length_feature_provenance,
)

LENGTH_FEATURES_SCHEMA_VERSION = "length-features-v1"
REGION_NAMES = ("CORE", "INVARIANT", "ARRAY", "LEFT_FLANK", "RIGHT_FLANK")
LengthFeatureStatus = Literal["measured", "partial", "unavailable"]
_SHA256_LENGTH = 64


@dataclass(frozen=True)
class DepthPosition:
    """One explicitly zero-based depth observation with optional transient support IDs."""

    contig: str
    position_zero_based: int
    depth: int
    supporting_fragment_ids: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        _text(self.contig, "depth contig")
        if not _is_int(self.position_zero_based) or self.position_zero_based < 0:
            raise ValueError("depth position must be a zero-based non-negative integer")
        if not _is_int(self.depth) or self.depth < 0:
            raise ValueError("depth must be a non-negative integer and cannot be boolean")
        fragment_ids = self.supporting_fragment_ids
        if fragment_ids is not None:
            if not isinstance(fragment_ids, tuple) or any(
                not isinstance(item, str) or not item or item.strip() != item for item in fragment_ids
            ):
                raise ValueError("supporting fragment IDs must be null or a tuple of non-empty strings")
            if len(fragment_ids) != len(set(fragment_ids)):
                raise ValueError("supporting fragment IDs must be unique at each depth position")
            if self.depth == 0 and fragment_ids:
                raise ValueError("supporting fragment IDs must be empty at zero depth")
            if self.depth > 0 and not fragment_ids:
                raise ValueError("positive depth requires non-empty supporting fragment IDs when evidence is declared")
            if len(fragment_ids) > self.depth:
                raise ValueError("distinct supporting fragment count cannot exceed depth")


@dataclass(frozen=True)
class RegionFeatures:
    """Depth statistics for every base of one annotated region."""

    length_bp: int
    depth_sum: float
    mean_depth: float
    covered_fraction: float
    supporting_fragment_count: int | None


@dataclass(frozen=True)
class LengthFeatures:
    """Immutable feature row and digest of its exact serialized artifact."""

    manifest_key: str
    assembly: str
    assay_class: str
    input_scope: str
    annotation_sha256: str
    counting_policy_sha256: str
    provenance_sha256: str
    regions: Mapping[str, RegionFeatures | None]
    a: float | None
    f: float | None
    status: LengthFeatureStatus
    reasons: tuple[str, ...]
    provenance: LengthFeatureProvenance
    sha256: str

    def __post_init__(self) -> None:
        if isinstance(self.regions, Mapping):
            object.__setattr__(self, "regions", MappingProxyType(dict(self.regions)))


def extract_length_features(
    depth: Sequence[DepthPosition],
    annotation: LengthAnnotation,
    context: LengthFeatureContext,
) -> LengthFeatures:
    """Compute strict region summaries and the A/F depth ratios.

    Args:
        depth: Exactly one observation for each position in the union of declared regions.
        annotation: Validated immutable region annotation.
        context: Validated input identity and counting policy from the depth adapter.

    Returns:
        An immutable feature row. Missing annotation or zero denominators produce null ratios.

    Raises:
        ValueError: If inputs, positions, contigs, depth values, or context bindings are incompatible.
    """
    if not isinstance(annotation, LengthAnnotation):
        raise ValueError("length feature annotation must be a LengthAnnotation")
    if not isinstance(context, LengthFeatureContext):
        raise ValueError("length feature context must be a LengthFeatureContext")
    encode_length_feature_context(context)
    _validate_annotation(annotation)
    _validate_context(annotation, context)
    if not isinstance(depth, Sequence) or isinstance(depth, (str, bytes)):
        raise ValueError("length feature depth must be a sequence of DepthPosition values")
    if any(not isinstance(item, DepthPosition) for item in depth):
        raise ValueError("length feature depth must contain only DepthPosition values")

    contigs = {item.contig for item in depth}
    if len(contigs) > 1:
        raise ValueError("length feature depth contains mixed contigs")
    if contigs and contigs != {context.original_contig}:
        raise ValueError("length feature depth contig does not match context original contig")
    positions = [item.position_zero_based for item in depth]
    if len(positions) != len(set(positions)):
        raise ValueError("length feature depth contains a duplicate position")

    intervals = _annotation_intervals(annotation)
    expected_positions = {position for interval in intervals for position in range(interval.start, interval.end)}
    if context.counting_policy.queried_intervals != _positions_to_intervals(expected_positions):
        raise ValueError("length feature queried interval union does not match annotation positions")
    observed_positions = set(positions)
    outside = observed_positions - expected_positions
    if outside:
        raise ValueError("length feature depth contains an out-of-range position")
    missing = expected_positions - observed_positions
    if missing:
        raise ValueError("length feature depth has a missing declared position")
    fragment_availability = {item.supporting_fragment_ids is not None for item in depth}
    if len(fragment_availability) > 1:
        raise ValueError("length feature depth has incomplete fragment metadata")
    fragments_available = fragment_availability == {True}

    by_position = {item.position_zero_based: item for item in depth}
    region_intervals = _region_intervals(annotation)
    summaries = {
        name: _summarize_region(region_intervals[name], by_position, fragments_available) for name in REGION_NAMES
    }
    denominator_qc = _denominator_qc(region_intervals, summaries, by_position, fragments_available)
    provenance = bind_length_feature_provenance(context, denominator_qc)
    a, f, status, reasons = _ratios_status_and_reasons(summaries)

    features = LengthFeatures(
        manifest_key=provenance.manifest_key,
        assembly=provenance.assembly,
        assay_class=provenance.assay_class,
        input_scope=provenance.input_scope,
        annotation_sha256=provenance.annotation_sha256,
        counting_policy_sha256=provenance.counting_policy_sha256,
        provenance_sha256=provenance.sha256,
        regions=MappingProxyType(summaries),
        a=a,
        f=f,
        status=status,
        reasons=reasons,
        provenance=provenance,
        sha256="",
    )
    return replace(features, sha256=canonical_sha256(_encode_length_features_unchecked(features)))


def _ratios_status_and_reasons(
    summaries: Mapping[str, RegionFeatures | None],
) -> tuple[float | None, float | None, LengthFeatureStatus, tuple[str, ...]]:
    reasons: list[str] = []
    for name, reason in (
        ("CORE", "missing_core_annotation"),
        ("INVARIANT", "missing_invariant_annotation"),
        ("ARRAY", "missing_array_annotation"),
        ("LEFT_FLANK", "missing_left_flank_annotation"),
        ("RIGHT_FLANK", "missing_right_flank_annotation"),
    ):
        if summaries[name] is None:
            reasons.append(reason)

    core = summaries["CORE"]
    invariant = summaries["INVARIANT"]
    array = summaries["ARRAY"]
    left_flank = summaries["LEFT_FLANK"]
    right_flank = summaries["RIGHT_FLANK"]
    a: float | None = None
    if core is not None and invariant is not None:
        if invariant.mean_depth == 0:
            reasons.append("zero_invariant_mean_depth")
        else:
            a = core.mean_depth / invariant.mean_depth
    f: float | None = None
    if array is not None and left_flank is not None and right_flank is not None:
        flank_length = left_flank.length_bp + right_flank.length_bp
        flank_mean = (left_flank.depth_sum + right_flank.depth_sum) / flank_length
        if flank_mean == 0:
            reasons.append("zero_combined_flank_mean_depth")
        else:
            f = array.mean_depth / flank_mean

    present_count = sum(summary is not None for summary in summaries.values())
    if present_count == 0:
        status: LengthFeatureStatus = "unavailable"
    elif a is not None and f is not None and present_count == len(REGION_NAMES):
        status = "measured"
    else:
        status = "partial"
    return a, f, status, tuple(reasons)


def encode_length_features(features: LengthFeatures) -> dict[str, object]:
    """Encode one feature row in the exact ``length-features-v1`` artifact schema.

    Args:
        features: Feature row returned by :func:`extract_length_features`.

    Returns:
        A canonical JSON-compatible artifact containing one row.

    Raises:
        ValueError: If the input is not a ``LengthFeatures`` value.
    """
    if not isinstance(features, LengthFeatures):
        raise ValueError("length features must be a LengthFeatures value")
    _validate_length_features(features)
    return _encode_length_features_unchecked(features)


def _encode_length_features_unchecked(features: LengthFeatures) -> dict[str, object]:
    row: dict[str, object] = {
        "manifest_key": features.manifest_key,
        "assembly": features.assembly,
        "assay_class": features.assay_class,
        "input_scope": features.input_scope,
        "annotation_sha256": features.annotation_sha256,
        "counting_policy_sha256": features.counting_policy_sha256,
        "provenance_sha256": features.provenance_sha256,
        "regions": {name: _encode_region(features.regions[name]) for name in REGION_NAMES},
        "A": features.a,
        "F": features.f,
        "status": features.status,
        "reasons": list(features.reasons),
    }
    return {"schema_version": LENGTH_FEATURES_SCHEMA_VERSION, "rows": [row]}


def _validate_annotation(annotation: LengthAnnotation) -> None:
    try:
        encoded = encode_length_annotation(annotation)
        decoded = decode_length_annotation(encoded)
    except (AttributeError, KeyError, TypeError) as error:
        raise ValueError("length feature annotation contains invalid typed content") from error
    if not isinstance(annotation.sha256, str) or annotation.sha256 != decoded.sha256:
        raise ValueError("length feature annotation digest does not match its canonical content")


def _validate_context(annotation: LengthAnnotation, context: LengthFeatureContext) -> None:
    if context.assembly != annotation.assembly:
        raise ValueError("length feature assembly mismatch")
    if context.original_contig not in annotation.accepted_contigs:
        raise ValueError("length feature contig mismatch")
    if context.reference_fasta_sha256 != annotation.reference_fasta_sha256:
        raise ValueError("length feature reference digest mismatch")
    if context.annotation_sha256 != annotation.sha256:
        raise ValueError("length feature annotation digest mismatch")


def _validate_length_features(features: LengthFeatures) -> None:
    provenance_document = encode_length_feature_provenance(features.provenance)
    provenance = features.provenance
    for value, label in (
        (features.manifest_key, "length feature manifest key"),
        (features.assembly, "length feature assembly"),
        (features.assay_class, "length feature assay class"),
        (features.input_scope, "length feature input scope"),
    ):
        _text(value, label)
    if features.input_scope not in {"full", "regional"}:
        raise ValueError("length feature input scope must be full or regional")
    for value, label in (
        (features.annotation_sha256, "length feature annotation digest"),
        (features.counting_policy_sha256, "length feature counting policy digest"),
        (features.provenance_sha256, "length feature provenance digest"),
    ):
        _digest(value, label)
    if (
        features.manifest_key,
        features.assembly,
        features.assay_class,
        features.input_scope,
        features.annotation_sha256,
        features.counting_policy_sha256,
        features.provenance_sha256,
    ) != (
        provenance.manifest_key,
        provenance.assembly,
        provenance.assay_class,
        provenance.input_scope,
        provenance.annotation_sha256,
        provenance.counting_policy_sha256,
        canonical_sha256(provenance_document),
    ):
        raise ValueError("length feature row identity does not match its bound provenance")
    if not isinstance(features.regions, Mapping) or set(features.regions) != set(REGION_NAMES):
        raise ValueError("length feature regions must contain exactly the five declared region names")
    for name in REGION_NAMES:
        _validate_region(features.regions[name], name)
    _validate_feature_qc(features.regions, provenance.denominator_qc)
    expected_a, expected_f, expected_status, expected_reasons = _ratios_status_and_reasons(features.regions)
    if (features.a, features.f, features.status, features.reasons) != (
        expected_a,
        expected_f,
        expected_status,
        expected_reasons,
    ):
        raise ValueError("length feature ratios, status, or reasons do not match region values")
    if not isinstance(features.sha256, str) or features.sha256 != canonical_sha256(
        _encode_length_features_unchecked(features)
    ):
        raise ValueError("length features digest does not match its canonical content")


def _validate_region(region: RegionFeatures | None, name: str) -> None:
    if region is None:
        return
    if not isinstance(region, RegionFeatures):
        raise ValueError(f"length feature {name} region must be RegionFeatures or null")
    if not _is_int(region.length_bp) or region.length_bp <= 0:
        raise ValueError(f"length feature {name} length must be a positive integer")
    for value, label in (
        (region.depth_sum, "depth sum"),
        (region.mean_depth, "mean depth"),
        (region.covered_fraction, "covered fraction"),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"length feature {name} {label} must be finite numeric")
    if region.depth_sum < 0 or region.mean_depth < 0 or not 0 <= region.covered_fraction <= 1:
        raise ValueError(f"length feature {name} numeric values are out of range")
    if region.mean_depth != region.depth_sum / region.length_bp:
        raise ValueError(f"length feature {name} mean does not match its sum and length")
    support = region.supporting_fragment_count
    if support is not None and (not _is_int(support) or support < 0 or support > region.depth_sum):
        raise ValueError(f"length feature {name} support count is invalid")


def _validate_feature_qc(regions: Mapping[str, RegionFeatures | None], qc: DenominatorQc) -> None:
    counts = qc.supporting_fragment_counts
    support_values = [region.supporting_fragment_count for region in regions.values() if region is not None]
    if qc.evidence_kind == "unavailable" and any(value is not None for value in support_values):
        raise ValueError("unavailable length feature fragment evidence requires null region support")
    if qc.evidence_kind == "read-pair-identity-qc-proxy" and any(value is None for value in support_values):
        raise ValueError("declared length feature fragment evidence requires region support counts")
    for name in ("CORE", "INVARIANT", "ARRAY"):
        region = regions[name]
        expected = None if region is None else region.supporting_fragment_count
        if counts[name] != expected:
            raise ValueError(f"length feature {name} support does not match bound denominator QC")
    left = regions["LEFT_FLANK"]
    right = regions["RIGHT_FLANK"]
    both = counts["BOTH_FLANKS"]
    if left is None or right is None:
        if both is not None:
            raise ValueError("length feature BOTH_FLANKS support requires both flank regions")
        return
    left_count = left.supporting_fragment_count
    right_count = right.supporting_fragment_count
    if left_count is None or right_count is None:
        if both is not None:
            raise ValueError("length feature BOTH_FLANKS support must be null when fragment evidence is unavailable")
    elif both is None or not max(left_count, right_count) <= both <= left_count + right_count:
        raise ValueError("length feature BOTH_FLANKS support is inconsistent with individual flank support")


def _annotation_intervals(annotation: LengthAnnotation) -> tuple[Interval, ...]:
    result: list[Interval] = []
    for intervals in _region_intervals(annotation).values():
        if intervals is not None:
            result.extend(intervals)
    return tuple(result)


def _positions_to_intervals(positions: set[int]) -> tuple[QueryInterval, ...]:
    if not positions:
        return ()
    sorted_positions = sorted(positions)
    starts_and_ends: list[tuple[int, int]] = []
    start = sorted_positions[0]
    previous = start
    for position in sorted_positions[1:]:
        if position != previous + 1:
            starts_and_ends.append((start, previous + 1))
            start = position
        previous = position
    starts_and_ends.append((start, previous + 1))
    return tuple(QueryInterval(start, end) for start, end in starts_and_ends)


def _region_intervals(annotation: LengthAnnotation) -> dict[str, tuple[Interval, ...] | None]:
    return {
        "CORE": annotation.core,
        "INVARIANT": annotation.invariant,
        "ARRAY": _as_interval_tuple(annotation.array),
        "LEFT_FLANK": _as_interval_tuple(annotation.left_flank),
        "RIGHT_FLANK": _as_interval_tuple(annotation.right_flank),
    }


def _as_interval_tuple(interval: Interval | None) -> tuple[Interval, ...] | None:
    if interval is None:
        return None
    return (interval,)


def _summarize_region(
    intervals: tuple[Interval, ...] | None,
    by_position: Mapping[int, DepthPosition],
    fragments_available: bool,
) -> RegionFeatures | None:
    if intervals is None:
        return None
    observations = [by_position[position] for interval in intervals for position in range(interval.start, interval.end)]
    length_bp = len(observations)
    depth_sum = float(sum(item.depth for item in observations))
    covered_fraction = sum(item.depth > 0 for item in observations) / length_bp
    supporting_fragment_count: int | None = None
    if fragments_available:
        supporting_fragment_count = len(
            {
                fragment_id
                for item in observations
                for fragment_id in cast(tuple[str, ...], item.supporting_fragment_ids)
            }
        )
    return RegionFeatures(
        length_bp=length_bp,
        depth_sum=depth_sum,
        mean_depth=depth_sum / length_bp,
        covered_fraction=covered_fraction,
        supporting_fragment_count=supporting_fragment_count,
    )


def _denominator_qc(
    region_intervals: Mapping[str, tuple[Interval, ...] | None],
    summaries: Mapping[str, RegionFeatures | None],
    by_position: Mapping[int, DepthPosition],
    fragments_available: bool,
) -> DenominatorQc:
    if not fragments_available:
        return DenominatorQc(
            evidence_kind="unavailable",
            supporting_fragment_counts=MappingProxyType(dict.fromkeys(DENOMINATOR_NAMES)),
        )
    counts: dict[str, int | None] = {}
    for name in ("CORE", "INVARIANT", "ARRAY"):
        summary = summaries[name]
        counts[name] = None if summary is None else summary.supporting_fragment_count
    left = region_intervals["LEFT_FLANK"]
    right = region_intervals["RIGHT_FLANK"]
    counts["BOTH_FLANKS"] = None if left is None or right is None else _fragment_count((*left, *right), by_position)
    return DenominatorQc(
        evidence_kind="read-pair-identity-qc-proxy",
        supporting_fragment_counts=MappingProxyType(counts),
    )


def _fragment_count(intervals: tuple[Interval, ...], by_position: Mapping[int, DepthPosition]) -> int:
    return len(
        {
            fragment_id
            for interval in intervals
            for position in range(interval.start, interval.end)
            for fragment_id in cast(tuple[str, ...], by_position[position].supporting_fragment_ids)
        }
    )


def _encode_region(region: RegionFeatures | None) -> dict[str, object] | None:
    if region is None:
        return None
    return {
        "length_bp": region.length_bp,
        "depth_sum": region.depth_sum,
        "mean_depth": region.mean_depth,
        "covered_fraction": region.covered_fraction,
        "supporting_fragment_count": region.supporting_fragment_count,
    }


def _digest(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != _SHA256_LENGTH
        or any(char not in "0123456789abcdef" for char in value)
    ):
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a non-empty trimmed string")
    return value


def _is_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)
