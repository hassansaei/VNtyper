"""Closed standard 13-feature measurement contract for MUC1 length research."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from statistics import fmean, pstdev
from types import MappingProxyType

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_features import DepthPosition

STANDARD_FEATURE_ORDER = (
    "A",
    "F",
    "core_zero_fraction",
    "core_depth_cv",
    "core_bin_cv",
    "core_half_log_ratio",
    "invariant_end_log_ratio",
    "flank_end_log_ratio",
    "log_invariant_depth",
    "mapq_zero_fraction",
    "mean_mapq",
    "soft_clipped_read_fraction",
    "query_sequence_gc_fraction",
)
STANDARD_ASSEMBLY = "GRCh38"
STANDARD_ACCEPTED_CONTIGS = ("chr1", "1")
STANDARD_LOCUS_START = 155188296
STANDARD_LOCUS_END = 155192429
STANDARD_REFERENCE_LOCUS_SHA256 = "e7aed8bef66426a3839471ec5f5e96c4a870d387fa286ebe9a245ee0c6f86e2d"
STANDARD_ANNOTATION_SHA256 = "32aa3a8e59bfacadd8d604b27ab6167d5d6a71f82aeae24a1e67c32f8cab5dd7"

_REGIONS = MappingProxyType(
    {
        "CORE": ((155188726, 155191939),),
        "INVARIANT_LEFT": ((155188486, 155188726),),
        "INVARIANT_RIGHT": ((155191939, 155192239),),
        "ARRAY": ((155188529, 155192010),),
        "LEFT_FLANK": ((155188296, 155188486),),
        "RIGHT_FLANK": ((155192239, 155192429),),
    }
)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_MEASUREMENT_FIELDS = {
    "schema_version",
    "assembly",
    "contig",
    "feature_order",
    "values",
    "qc",
    "reasons",
    "annotation_sha256",
    "reference_locus_sha256",
    "feature_definition_sha256",
    "sha256",
}
_QC_FIELDS = {
    "eligible_read_count",
    "invariant_mean_depth",
    "invariant_covered_fraction",
    "invariant_supporting_fragments",
    "left_flank_mean_depth",
    "left_flank_covered_fraction",
    "left_flank_supporting_fragments",
    "right_flank_mean_depth",
    "right_flank_covered_fraction",
    "right_flank_supporting_fragments",
    "combined_flank_mean_depth",
    "combined_flank_covered_fraction",
    "combined_flank_supporting_fragments",
}


def standard_feature_definition_document() -> dict[str, object]:
    """Return the frozen public geometry and mathematical feature definitions."""
    return {
        "schema_version": "standard-length-feature-definition-v1",
        "assembly": STANDARD_ASSEMBLY,
        "accepted_contigs": list(STANDARD_ACCEPTED_CONTIGS),
        "locus": {"start": STANDARD_LOCUS_START, "end": STANDARD_LOCUS_END},
        "regions": {name: [list(interval) for interval in intervals] for name, intervals in _REGIONS.items()},
        "record_policy": {
            "excluded_sam_flags": 0x704,
            "supplementary_included": True,
            "mates_counted_separately": True,
            "read_statistics_require_nonempty_query_sequence": True,
            "base_depth": "aligned-query-and-reference-pairs-base-quality-at-least-zero",
        },
        "feature_order": list(STANDARD_FEATURE_ORDER),
        "formula_version": "standard-length-formulas-v1",
        "core_bin_count": 8,
        "log_pseudocount": 1,
    }


STANDARD_FEATURE_DEFINITION_SHA256 = canonical_sha256(standard_feature_definition_document())


@dataclass(frozen=True)
class StandardReadSummary:
    """Aggregate read-level evidence from one indexed locus traversal."""

    eligible_read_count: int
    mapq_zero_count: int
    mapq_sum: int
    soft_clipped_read_count: int
    query_gc_bases: int
    query_bases: int

    def __post_init__(self) -> None:
        values = (
            self.eligible_read_count,
            self.mapq_zero_count,
            self.mapq_sum,
            self.soft_clipped_read_count,
            self.query_gc_bases,
            self.query_bases,
        )
        if any(not _is_int(value) or value < 0 for value in values):
            raise ValueError("standard length read summary values must be non-negative integers")
        if self.mapq_zero_count > self.eligible_read_count:
            raise ValueError("standard length MAPQ-zero count exceeds eligible reads")
        if self.soft_clipped_read_count > self.eligible_read_count:
            raise ValueError("standard length soft-clipped count exceeds eligible reads")
        if self.query_gc_bases > self.query_bases:
            raise ValueError("standard length GC bases exceed query bases")
        if self.eligible_read_count == 0 and any(values[1:]):
            raise ValueError("empty standard length read evidence requires zero aggregates")
        if self.eligible_read_count > 0 and self.query_bases == 0:
            raise ValueError("eligible standard length reads require query bases")


@dataclass(frozen=True)
class StandardLengthQc:
    """Observable denominator and read evidence used by the standard model QC."""

    eligible_read_count: int
    invariant_mean_depth: float
    invariant_covered_fraction: float
    invariant_supporting_fragments: int | None
    left_flank_mean_depth: float
    left_flank_covered_fraction: float
    left_flank_supporting_fragments: int | None
    right_flank_mean_depth: float
    right_flank_covered_fraction: float
    right_flank_supporting_fragments: int | None
    combined_flank_mean_depth: float
    combined_flank_covered_fraction: float
    combined_flank_supporting_fragments: int | None


@dataclass(frozen=True)
class StandardLengthMeasurement:
    """Immutable standard predictor vector, QC evidence, and canonical identity."""

    assembly: str
    contig: str
    feature_order: tuple[str, ...]
    values: Mapping[str, float | None]
    qc: StandardLengthQc
    reasons: tuple[str, ...]
    annotation_sha256: str
    reference_locus_sha256: str
    feature_definition_sha256: str
    sha256: str

    def __post_init__(self) -> None:
        if isinstance(self.values, Mapping):
            object.__setattr__(self, "values", MappingProxyType(dict(self.values)))


def _positions(intervals: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    return tuple(position for start, end in intervals for position in range(start, end))


def _values(name: str, by_position: Mapping[int, DepthPosition]) -> tuple[DepthPosition, ...]:
    return tuple(by_position[position] for position in _positions(_REGIONS[name]))


def _mean(rows: Sequence[DepthPosition]) -> float:
    return fmean(row.depth for row in rows)


def _covered(rows: Sequence[DepthPosition]) -> float:
    return sum(row.depth > 0 for row in rows) / len(rows)


def _support(rows: Sequence[DepthPosition]) -> int | None:
    if any(row.supporting_fragment_ids is None for row in rows):
        return None
    return len({fragment for row in rows for fragment in (row.supporting_fragment_ids or ())})


def _bins(values: Sequence[int], count: int) -> tuple[tuple[int, ...], ...]:
    quotient, remainder = divmod(len(values), count)
    result = []
    offset = 0
    for index in range(count):
        length = quotient + (1 if index < remainder else 0)
        result.append(tuple(values[offset : offset + length]))
        offset += length
    return tuple(result)


def extract_standard_length_features(
    depth: Sequence[DepthPosition],
    read_summary: StandardReadSummary,
    *,
    assembly: str,
    contig: str,
    reference_locus_sha256: str,
) -> StandardLengthMeasurement:
    """Derive the frozen 13-feature vector from exact locus depth and read aggregates."""
    if assembly != STANDARD_ASSEMBLY:
        raise ValueError("standard length features require GRCh38")
    if contig not in STANDARD_ACCEPTED_CONTIGS:
        raise ValueError("standard length feature contig is unsupported")
    if reference_locus_sha256 != STANDARD_REFERENCE_LOCUS_SHA256:
        raise ValueError("standard length reference locus digest differs")
    if not isinstance(read_summary, StandardReadSummary):
        raise ValueError("standard length read summary must be typed")
    if not isinstance(depth, Sequence) or isinstance(depth, (str, bytes)):
        raise ValueError("standard length depth must be a sequence")
    if any(not isinstance(row, DepthPosition) for row in depth):
        raise ValueError("standard length depth contains an invalid row")
    positions = [row.position_zero_based for row in depth]
    expected = set(range(STANDARD_LOCUS_START, STANDARD_LOCUS_END))
    if len(positions) != len(set(positions)) or set(positions) != expected:
        raise ValueError("standard length depth must cover the exact locus once")
    if {row.contig for row in depth} != {contig}:
        raise ValueError("standard length depth contig differs")
    availability = {row.supporting_fragment_ids is None for row in depth}
    if len(availability) > 1:
        raise ValueError("standard length fragment evidence is incomplete")

    by_position = {row.position_zero_based: row for row in depth}
    core = _values("CORE", by_position)
    invariant_left = _values("INVARIANT_LEFT", by_position)
    invariant_right = _values("INVARIANT_RIGHT", by_position)
    invariant = (*invariant_left, *invariant_right)
    array = _values("ARRAY", by_position)
    left_flank = _values("LEFT_FLANK", by_position)
    right_flank = _values("RIGHT_FLANK", by_position)
    flanks = (*left_flank, *right_flank)
    core_depth = tuple(row.depth for row in core)
    core_mean = fmean(core_depth)
    invariant_mean = _mean(invariant)
    flank_mean = _mean(flanks)
    bin_means = tuple(fmean(part) for part in _bins(core_depth, 8))
    half = len(core_depth) // 2

    reasons: list[str] = []
    a = None if invariant_mean == 0 else core_mean / invariant_mean
    if a is None:
        reasons.append("zero_invariant_mean_depth")
    f = None if flank_mean == 0 else _mean(array) / flank_mean
    if f is None:
        reasons.append("zero_combined_flank_mean_depth")
    if read_summary.eligible_read_count == 0:
        reasons.append("no_eligible_reads")
    values: dict[str, float | None] = {
        "A": a,
        "F": f,
        "core_zero_fraction": sum(value == 0 for value in core_depth) / len(core_depth),
        "core_depth_cv": None if core_mean == 0 else pstdev(core_depth) / core_mean,
        "core_bin_cv": None if fmean(bin_means) == 0 else pstdev(bin_means) / fmean(bin_means),
        "core_half_log_ratio": math.log((fmean(core_depth[:half]) + 1) / (fmean(core_depth[half:]) + 1)),
        "invariant_end_log_ratio": math.log((_mean(invariant_left) + 1) / (_mean(invariant_right) + 1)),
        "flank_end_log_ratio": math.log((_mean(left_flank) + 1) / (_mean(right_flank) + 1)),
        "log_invariant_depth": math.log1p(invariant_mean),
        "mapq_zero_fraction": (
            None
            if read_summary.eligible_read_count == 0
            else read_summary.mapq_zero_count / read_summary.eligible_read_count
        ),
        "mean_mapq": (
            None if read_summary.eligible_read_count == 0 else read_summary.mapq_sum / read_summary.eligible_read_count
        ),
        "soft_clipped_read_fraction": (
            None
            if read_summary.eligible_read_count == 0
            else read_summary.soft_clipped_read_count / read_summary.eligible_read_count
        ),
        "query_sequence_gc_fraction": (
            None if read_summary.query_bases == 0 else read_summary.query_gc_bases / read_summary.query_bases
        ),
    }
    reasons.extend(
        f"zero_{name.removesuffix('_cv')}_mean_depth"
        for name in ("core_depth_cv", "core_bin_cv")
        if values[name] is None
    )

    qc = StandardLengthQc(
        read_summary.eligible_read_count,
        invariant_mean,
        _covered(invariant),
        _support(invariant),
        _mean(left_flank),
        _covered(left_flank),
        _support(left_flank),
        _mean(right_flank),
        _covered(right_flank),
        _support(right_flank),
        flank_mean,
        _covered(flanks),
        _support(flanks),
    )
    measurement = StandardLengthMeasurement(
        assembly,
        contig,
        STANDARD_FEATURE_ORDER,
        MappingProxyType(values),
        qc,
        tuple(reasons),
        STANDARD_ANNOTATION_SHA256,
        reference_locus_sha256,
        STANDARD_FEATURE_DEFINITION_SHA256,
        "",
    )
    return replace(measurement, sha256=canonical_sha256(_measurement_payload(measurement)))


def _qc_document(qc: StandardLengthQc) -> dict[str, object]:
    return {name: getattr(qc, name) for name in sorted(_QC_FIELDS)}


def _measurement_payload(measurement: StandardLengthMeasurement) -> dict[str, object]:
    return {
        "schema_version": "standard-length-measurement-v1",
        "assembly": measurement.assembly,
        "contig": measurement.contig,
        "feature_order": list(measurement.feature_order),
        "values": {name: measurement.values[name] for name in STANDARD_FEATURE_ORDER},
        "qc": _qc_document(measurement.qc),
        "reasons": list(measurement.reasons),
        "annotation_sha256": measurement.annotation_sha256,
        "reference_locus_sha256": measurement.reference_locus_sha256,
        "feature_definition_sha256": measurement.feature_definition_sha256,
    }


def _finite(value: object, label: str, *, nonnegative: bool = False, fraction: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be finite numeric")
    try:
        result = float(value)
    except OverflowError:
        raise ValueError(f"{label} must be finite numeric") from None
    if not math.isfinite(result) or (nonnegative and result < 0) or (fraction and not 0 <= result <= 1):
        raise ValueError(f"{label} is out of range")
    return result


def _optional_count(value: object, label: str) -> int | None:
    if value is None:
        return None
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{label} must be a non-negative integer or null")
    return value


def _required_count(value: object, label: str) -> int:
    result = _optional_count(value, label)
    if result is None:
        raise ValueError(f"{label} must be a non-negative integer")
    return result


def encode_standard_length_measurement(measurement: StandardLengthMeasurement) -> dict[str, object]:
    """Project and revalidate a standard measurement as closed JSON content."""
    if not isinstance(measurement, StandardLengthMeasurement):
        raise ValueError("standard length measurement must be typed")
    decoded = decode_standard_length_measurement({**_measurement_payload(measurement), "sha256": measurement.sha256})
    if decoded != measurement:
        raise ValueError("standard length measurement typed content differs")
    return {**_measurement_payload(measurement), "sha256": measurement.sha256}


def decode_standard_length_measurement(value: object) -> StandardLengthMeasurement:
    """Decode a closed standard measurement and verify its canonical digest."""
    if not isinstance(value, Mapping) or set(value) != _MEASUREMENT_FIELDS:
        raise ValueError("standard length measurement fields differ")
    if value["schema_version"] != "standard-length-measurement-v1":
        raise ValueError("standard length measurement schema is unsupported")
    if value["assembly"] != STANDARD_ASSEMBLY or value["contig"] not in STANDARD_ACCEPTED_CONTIGS:
        raise ValueError("standard length measurement reference identity differs")
    if value["feature_order"] != list(STANDARD_FEATURE_ORDER):
        raise ValueError("standard length measurement feature order differs")
    raw_values = value["values"]
    if not isinstance(raw_values, Mapping) or set(raw_values) != set(STANDARD_FEATURE_ORDER):
        raise ValueError("standard length measurement values differ from feature order")
    values = {
        name: None if raw_values[name] is None else _finite(raw_values[name], f"standard feature {name}")
        for name in STANDARD_FEATURE_ORDER
    }
    raw_qc = value["qc"]
    if not isinstance(raw_qc, Mapping) or set(raw_qc) != _QC_FIELDS:
        raise ValueError("standard length measurement QC fields differ")
    qc = StandardLengthQc(
        eligible_read_count=_required_count(raw_qc["eligible_read_count"], "eligible reads"),
        invariant_mean_depth=_finite(raw_qc["invariant_mean_depth"], "invariant mean depth", nonnegative=True),
        invariant_covered_fraction=_finite(
            raw_qc["invariant_covered_fraction"], "invariant covered fraction", fraction=True
        ),
        invariant_supporting_fragments=_optional_count(raw_qc["invariant_supporting_fragments"], "invariant support"),
        left_flank_mean_depth=_finite(raw_qc["left_flank_mean_depth"], "left flank mean depth", nonnegative=True),
        left_flank_covered_fraction=_finite(
            raw_qc["left_flank_covered_fraction"], "left flank covered fraction", fraction=True
        ),
        left_flank_supporting_fragments=_optional_count(
            raw_qc["left_flank_supporting_fragments"], "left flank support"
        ),
        right_flank_mean_depth=_finite(raw_qc["right_flank_mean_depth"], "right flank mean depth", nonnegative=True),
        right_flank_covered_fraction=_finite(
            raw_qc["right_flank_covered_fraction"], "right flank covered fraction", fraction=True
        ),
        right_flank_supporting_fragments=_optional_count(
            raw_qc["right_flank_supporting_fragments"], "right flank support"
        ),
        combined_flank_mean_depth=_finite(
            raw_qc["combined_flank_mean_depth"], "combined flank mean depth", nonnegative=True
        ),
        combined_flank_covered_fraction=_finite(
            raw_qc["combined_flank_covered_fraction"], "combined flank covered fraction", fraction=True
        ),
        combined_flank_supporting_fragments=_optional_count(
            raw_qc["combined_flank_supporting_fragments"], "combined flank support"
        ),
    )
    reasons = value["reasons"]
    if not isinstance(reasons, list) or any(not isinstance(reason, str) or not reason for reason in reasons):
        raise ValueError("standard length measurement reasons must be strings")
    if len(reasons) != len(set(reasons)):
        raise ValueError("standard length measurement reasons must be unique")
    for field, expected in (
        ("annotation_sha256", STANDARD_ANNOTATION_SHA256),
        ("reference_locus_sha256", STANDARD_REFERENCE_LOCUS_SHA256),
        ("feature_definition_sha256", STANDARD_FEATURE_DEFINITION_SHA256),
    ):
        if value[field] != expected:
            raise ValueError(f"standard length measurement {field} differs")
    digest = value["sha256"]
    if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
        raise ValueError("standard length measurement digest is invalid")
    measurement = StandardLengthMeasurement(
        STANDARD_ASSEMBLY,
        str(value["contig"]),
        STANDARD_FEATURE_ORDER,
        MappingProxyType(values),
        qc,
        tuple(reasons),
        STANDARD_ANNOTATION_SHA256,
        STANDARD_REFERENCE_LOCUS_SHA256,
        STANDARD_FEATURE_DEFINITION_SHA256,
        digest,
    )
    if canonical_sha256(_measurement_payload(measurement)) != digest:
        raise ValueError("standard length measurement digest differs")
    return measurement


def _is_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)
