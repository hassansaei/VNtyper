"""Pure depth summaries and ratios for VNTR total-length research."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_annotation import Interval, LengthAnnotation

LENGTH_FEATURES_SCHEMA_VERSION = "length-features-v1"
REGION_NAMES = ("CORE", "INVARIANT", "ARRAY", "LEFT_FLANK", "RIGHT_FLANK")
LengthFeatureStatus = Literal["measured", "partial", "unavailable"]
_SHA256_LENGTH = 64


@dataclass(frozen=True)
class DepthPosition:
    """One explicitly zero-based depth observation with optional transient support IDs."""

    contig: str
    position_zero_based: int
    depth: float
    supporting_fragment_ids: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        _text(self.contig, "depth contig")
        if not _is_int(self.position_zero_based) or self.position_zero_based < 0:
            raise ValueError("depth position must be a zero-based non-negative integer")
        if isinstance(self.depth, bool) or not isinstance(self.depth, (int, float)):
            raise ValueError("depth must be numeric and cannot be boolean")
        if not math.isfinite(self.depth) or self.depth < 0:
            raise ValueError("depth must be finite and non-negative")
        object.__setattr__(self, "depth", float(self.depth))
        fragment_ids = self.supporting_fragment_ids
        if fragment_ids is not None:
            if not isinstance(fragment_ids, tuple) or any(
                not isinstance(item, str) or not item for item in fragment_ids
            ):
                raise ValueError("supporting fragment IDs must be null or a tuple of non-empty strings")
            if len(fragment_ids) != len(set(fragment_ids)):
                raise ValueError("supporting fragment IDs must be unique at each depth position")


@dataclass(frozen=True)
class LengthFeatureProvenance:
    """Typed identity and digest bindings supplied by a future depth adapter."""

    manifest_key: str
    assembly: str
    assay_class: str
    input_scope: str
    contig: str
    reference_fasta_sha256: str
    annotation_sha256: str
    counting_policy_sha256: str
    provenance_sha256: str
    aligner_name: str
    aligner_version: str
    aligner_arguments_sha256: str
    primary_secondary_marking: str
    preprocessing_id: str

    def __post_init__(self) -> None:
        for value, label in (
            (self.manifest_key, "length feature manifest key"),
            (self.assembly, "length feature assembly"),
            (self.assay_class, "length feature assay class"),
            (self.input_scope, "length feature input scope"),
            (self.contig, "length feature contig"),
            (self.aligner_name, "length feature aligner name"),
            (self.aligner_version, "length feature aligner version"),
            (self.primary_secondary_marking, "length feature primary/secondary marking"),
            (self.preprocessing_id, "length feature preprocessing ID"),
        ):
            _text(value, label)
        if self.input_scope not in {"full", "regional"}:
            raise ValueError("length feature input scope must be full or regional")
        for value, label in (
            (self.reference_fasta_sha256, "length feature reference FASTA digest"),
            (self.annotation_sha256, "length feature annotation digest"),
            (self.counting_policy_sha256, "length feature counting policy digest"),
            (self.provenance_sha256, "length feature provenance digest"),
            (self.aligner_arguments_sha256, "length feature aligner arguments digest"),
        ):
            _digest(value, label)


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
    sha256: str


def extract_length_features(
    depth: Sequence[DepthPosition],
    annotation: LengthAnnotation,
    provenance: LengthFeatureProvenance,
) -> LengthFeatures:
    """Compute strict region summaries and the A/F depth ratios.

    Args:
        depth: Exactly one observation for each position in the union of declared regions.
        annotation: Validated immutable region annotation.
        provenance: Input identity and digest bindings from the depth-producing boundary.

    Returns:
        An immutable feature row. Missing annotation or zero denominators produce null ratios.

    Raises:
        ValueError: If inputs, positions, contigs, depth values, or provenance bindings are incompatible.
    """
    if not isinstance(annotation, LengthAnnotation):
        raise ValueError("length feature annotation must be a LengthAnnotation")
    if not isinstance(provenance, LengthFeatureProvenance):
        raise ValueError("length feature provenance must be a LengthFeatureProvenance")
    _validate_provenance(annotation, provenance)
    if not isinstance(depth, Sequence) or isinstance(depth, (str, bytes)):
        raise ValueError("length feature depth must be a sequence of DepthPosition values")
    if any(not isinstance(item, DepthPosition) for item in depth):
        raise ValueError("length feature depth must contain only DepthPosition values")

    contigs = {item.contig for item in depth}
    if len(contigs) > 1:
        raise ValueError("length feature depth contains mixed contigs")
    if contigs and contigs != {provenance.contig}:
        raise ValueError("length feature depth contig does not match provenance contig")
    positions = [item.position_zero_based for item in depth]
    if len(positions) != len(set(positions)):
        raise ValueError("length feature depth contains a duplicate position")

    intervals = _annotation_intervals(annotation)
    expected_positions = {position for interval in intervals for position in range(interval.start, interval.end)}
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
    features = LengthFeatures(
        manifest_key=provenance.manifest_key,
        assembly=provenance.assembly,
        assay_class=provenance.assay_class,
        input_scope=provenance.input_scope,
        annotation_sha256=provenance.annotation_sha256,
        counting_policy_sha256=provenance.counting_policy_sha256,
        provenance_sha256=provenance.provenance_sha256,
        regions=MappingProxyType(summaries),
        a=a,
        f=f,
        status=status,
        reasons=tuple(reasons),
        sha256="",
    )
    return replace(features, sha256=canonical_sha256(encode_length_features(features)))


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


def _validate_provenance(annotation: LengthAnnotation, provenance: LengthFeatureProvenance) -> None:
    if provenance.assembly != annotation.assembly:
        raise ValueError("length feature assembly mismatch")
    if provenance.contig not in annotation.accepted_contigs:
        raise ValueError("length feature contig mismatch")
    if provenance.reference_fasta_sha256 != annotation.reference_fasta_sha256:
        raise ValueError("length feature reference digest mismatch")
    if provenance.annotation_sha256 != annotation.sha256:
        raise ValueError("length feature annotation digest mismatch")


def _annotation_intervals(annotation: LengthAnnotation) -> tuple[Interval, ...]:
    result: list[Interval] = []
    for intervals in _region_intervals(annotation).values():
        if intervals is not None:
            result.extend(intervals)
    return tuple(result)


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
    depth_sum = sum(item.depth for item in observations)
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
