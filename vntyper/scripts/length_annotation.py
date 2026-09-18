"""Strict immutable geometry for VNTR length measurements."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TypeGuard

from vntyper.scripts.canonical_json import canonical_sha256

LENGTH_ANNOTATION_SCHEMA_VERSION = "length-annotation-v1"
_REGION_NAMES = {"CORE", "INVARIANT", "ARRAY", "LEFT_FLANK", "RIGHT_FLANK"}
_SHA256_LENGTH = 64


@dataclass(frozen=True)
class Interval:
    """One zero-based, half-open interval."""

    start: int
    end: int

    def __post_init__(self) -> None:
        if not _is_int(self.start) or not _is_int(self.end) or self.start < 0 or self.end <= self.start:
            raise ValueError("length annotation interval must be a non-empty zero-based half-open interval")

    @property
    def length_bp(self) -> int:
        """Return the number of bases in the interval."""
        return self.end - self.start


@dataclass(frozen=True)
class LengthAnnotation:
    """Validated immutable annotation and its canonical content digest."""

    schema_version: str
    assembly: str
    contig: str
    accepted_contig_aliases: tuple[str, ...]
    reference_fasta_sha256: str
    coordinate_system: str
    boundary_definition: str
    repeat_unit_bp: int
    core: tuple[Interval, ...] | None
    invariant: tuple[Interval, ...] | None
    array: Interval | None
    left_flank: Interval | None
    right_flank: Interval | None
    array_only_bp: int | None
    target_only_bp: int | None
    target_boundary_conversion_sha256: str | None
    physical_a_compatible: bool
    physical_f_compatible: bool
    annotation_provenance: str
    annotation_version: str
    sha256: str

    @property
    def accepted_contigs(self) -> tuple[str, ...]:
        """Return the primary contig followed by its accepted aliases."""
        return (self.contig, *self.accepted_contig_aliases)

    @property
    def reference_core_repeat_count(self) -> int | None:
        """Return the verified complete CORE unit count when annotated."""
        if self.core is None:
            return None
        length_bp = sum(interval.length_bp for interval in self.core)
        if length_bp % self.repeat_unit_bp:
            return None
        return length_bp // self.repeat_unit_bp

    @property
    def reference_invariant_repeat_count(self) -> int | None:
        """Return the verified complete INVARIANT unit count when annotated."""
        if self.invariant is None:
            return None
        return sum(interval.length_bp for interval in self.invariant) // self.repeat_unit_bp


def one_based_closed_to_zero_based_half_open(start: object, end: object) -> Interval:
    """Convert a one-based closed interval to the annotation coordinate system.

    Args:
        start: Inclusive one-based start.
        end: Inclusive one-based end.

    Returns:
        The equivalent zero-based half-open interval.

    Raises:
        ValueError: If either boundary is not an integer or the interval is invalid.
    """
    if not _is_int(start) or not _is_int(end) or start < 1 or end < start:
        raise ValueError("one-based closed interval boundaries must be positive integers with start <= end")
    return Interval(start - 1, end)


def decode_length_annotation(value: object) -> LengthAnnotation:
    """Decode and validate a closed ``length-annotation-v1`` object.

    Args:
        value: Parsed JSON-compatible annotation object.

    Returns:
        An immutable annotation carrying its canonical digest.

    Raises:
        ValueError: If the schema, geometry, provenance, or compatibility claims are invalid.
    """
    root = _exact_object(
        value,
        {
            "schema_version",
            "assembly",
            "contig",
            "accepted_contig_aliases",
            "reference_fasta_sha256",
            "coordinate_system",
            "boundary_definition",
            "repeat_unit_bp",
            "regions",
            "array_boundary_geometry",
            "target_boundary_conversion_sha256",
            "physical_hypothesis_compatibility",
            "annotation_provenance",
            "annotation_version",
        },
        "length annotation",
    )
    if root["schema_version"] != LENGTH_ANNOTATION_SCHEMA_VERSION:
        raise ValueError(f"length annotation schema version must be {LENGTH_ANNOTATION_SCHEMA_VERSION}")
    if root["coordinate_system"] != "zero-based-half-open":
        raise ValueError("length annotation coordinate system must be zero-based-half-open")
    assembly = _text(root["assembly"], "length annotation assembly")
    contig = _text(root["contig"], "length annotation contig")
    aliases = _aliases(root["accepted_contig_aliases"], contig)
    reference_digest = _digest(root["reference_fasta_sha256"], "reference FASTA digest")
    boundary_definition = _text(root["boundary_definition"], "length annotation boundary definition")
    repeat_unit_bp = root["repeat_unit_bp"]
    if not _is_int(repeat_unit_bp) or repeat_unit_bp <= 0:
        raise ValueError("length annotation repeat unit bp must be a positive integer")

    regions = _exact_object(root["regions"], _REGION_NAMES, "length annotation regions")
    core = _optional_interval_list(regions["CORE"], "CORE")
    invariant = _optional_interval_list(regions["INVARIANT"], "INVARIANT")
    array = _optional_interval(regions["ARRAY"], "ARRAY")
    left_flank = _optional_interval(regions["LEFT_FLANK"], "LEFT_FLANK")
    right_flank = _optional_interval(regions["RIGHT_FLANK"], "RIGHT_FLANK")
    # A complete biological repeat can contain a sequence insertion or deletion.
    # Such a CORE remains measurable, but its bp span cannot establish a physical
    # reference repeat count by integer division.
    _validate_repeat_geometry(invariant, repeat_unit_bp, "INVARIANT")
    if core is not None and invariant is not None and _overlap_length(core, invariant):
        raise ValueError("length annotation CORE and INVARIANT intervals must be disjoint")
    _validate_flanks(left_flank, right_flank, core, invariant, array)

    geometry = _exact_object(
        root["array_boundary_geometry"], {"array_only_bp", "target_only_bp"}, "length annotation boundary geometry"
    )
    array_only_bp = _optional_nonnegative_int(geometry["array_only_bp"], "array-only bp")
    target_only_bp = _optional_nonnegative_int(geometry["target_only_bp"], "target-only bp")
    _validate_boundary_geometry(core, invariant, array, array_only_bp, target_only_bp)

    conversion_digest = _optional_digest(root["target_boundary_conversion_sha256"], "target boundary conversion digest")
    compatibility = _exact_object(
        root["physical_hypothesis_compatibility"],
        {"physical_A", "physical_F"},
        "length annotation physical hypothesis compatibility",
    )
    physical_a = _strict_bool(compatibility["physical_A"], "physical A compatibility")
    physical_f = _strict_bool(compatibility["physical_F"], "physical F compatibility")
    if (physical_a or physical_f) and conversion_digest is None:
        raise ValueError("physical length hypotheses require a target boundary conversion digest")
    if physical_a and (core is None or invariant is None):
        raise ValueError("physical A compatibility requires CORE and INVARIANT annotation")
    if physical_a and sum(interval.length_bp for interval in core or ()) % repeat_unit_bp:
        raise ValueError("physical A compatibility requires a fixed-width CORE reference count")
    if physical_f and (array is None or left_flank is None or right_flank is None):
        raise ValueError("physical F compatibility requires ARRAY and both flank annotations")

    return LengthAnnotation(
        schema_version=LENGTH_ANNOTATION_SCHEMA_VERSION,
        assembly=assembly,
        contig=contig,
        accepted_contig_aliases=aliases,
        reference_fasta_sha256=reference_digest,
        coordinate_system="zero-based-half-open",
        boundary_definition=boundary_definition,
        repeat_unit_bp=repeat_unit_bp,
        core=core,
        invariant=invariant,
        array=array,
        left_flank=left_flank,
        right_flank=right_flank,
        array_only_bp=array_only_bp,
        target_only_bp=target_only_bp,
        target_boundary_conversion_sha256=conversion_digest,
        physical_a_compatible=physical_a,
        physical_f_compatible=physical_f,
        annotation_provenance=_text(root["annotation_provenance"], "length annotation provenance"),
        annotation_version=_text(root["annotation_version"], "length annotation version"),
        sha256=canonical_sha256(root),
    )


def encode_length_annotation(annotation: LengthAnnotation) -> dict[str, object]:
    """Encode a validated annotation as exact JSON-compatible primitives.

    Args:
        annotation: Validated length annotation.

    Returns:
        The closed ``length-annotation-v1`` object.

    Raises:
        ValueError: If the value is not a validated annotation.
    """
    if not isinstance(annotation, LengthAnnotation):
        raise ValueError("length annotation must be a LengthAnnotation")
    return {
        "schema_version": annotation.schema_version,
        "assembly": annotation.assembly,
        "contig": annotation.contig,
        "accepted_contig_aliases": list(annotation.accepted_contig_aliases),
        "reference_fasta_sha256": annotation.reference_fasta_sha256,
        "coordinate_system": annotation.coordinate_system,
        "boundary_definition": annotation.boundary_definition,
        "repeat_unit_bp": annotation.repeat_unit_bp,
        "regions": {
            "CORE": _encode_interval_list(annotation.core),
            "INVARIANT": _encode_interval_list(annotation.invariant),
            "ARRAY": _encode_interval(annotation.array),
            "LEFT_FLANK": _encode_interval(annotation.left_flank),
            "RIGHT_FLANK": _encode_interval(annotation.right_flank),
        },
        "array_boundary_geometry": {
            "array_only_bp": annotation.array_only_bp,
            "target_only_bp": annotation.target_only_bp,
        },
        "target_boundary_conversion_sha256": annotation.target_boundary_conversion_sha256,
        "physical_hypothesis_compatibility": {
            "physical_A": annotation.physical_a_compatible,
            "physical_F": annotation.physical_f_compatible,
        },
        "annotation_provenance": annotation.annotation_provenance,
        "annotation_version": annotation.annotation_version,
    }


def _aliases(value: object, contig: str) -> tuple[str, ...]:
    if not isinstance(value, list) or any(not isinstance(alias, str) or not alias for alias in value):
        raise ValueError("length annotation accepted contig aliases must be a string list")
    aliases = tuple(value)
    if aliases != tuple(sorted(aliases)) or len(aliases) != len(set(aliases)) or contig in aliases:
        raise ValueError(
            "length annotation accepted contig aliases must be unique, increasing, and exclude the primary"
        )
    return aliases


def _optional_interval_list(value: object, name: str) -> tuple[Interval, ...] | None:
    if value is None:
        return None
    if not isinstance(value, list) or not value:
        raise ValueError(f"length annotation {name} must be null or a non-empty interval list")
    intervals = tuple(_decode_interval(item, name) for item in value)
    if intervals != tuple(sorted(intervals, key=lambda item: (item.start, item.end))):
        raise ValueError(f"length annotation {name} intervals must be increasing")
    if any(left.end > right.start for left, right in zip(intervals, intervals[1:], strict=False)):
        raise ValueError(f"length annotation {name} intervals must be disjoint")
    return intervals


def _optional_interval(value: object, name: str) -> Interval | None:
    if value is None:
        return None
    return _decode_interval(value, name)


def _decode_interval(value: object, name: str) -> Interval:
    raw = _exact_object(value, {"start", "end"}, f"length annotation {name} interval")
    try:
        return Interval(raw["start"], raw["end"])  # type: ignore[arg-type]
    except ValueError as exc:
        raise ValueError(f"length annotation {name} interval is invalid: {exc}") from exc


def _validate_repeat_geometry(intervals: tuple[Interval, ...] | None, unit_bp: int, name: str) -> None:
    if intervals is not None and any(interval.length_bp % unit_bp for interval in intervals):
        raise ValueError(f"length annotation {name} intervals must contain complete repeat units")


def _validate_flanks(
    left: Interval | None,
    right: Interval | None,
    core: tuple[Interval, ...] | None,
    invariant: tuple[Interval, ...] | None,
    array: Interval | None,
) -> None:
    measured = [interval for intervals in (core, invariant) if intervals is not None for interval in intervals]
    if array is not None:
        measured.append(array)
    if left is not None and right is not None and left.end > right.start:
        raise ValueError("length annotation flank intervals must be disjoint and ordered")
    if measured and left is not None and left.end > min(interval.start for interval in measured):
        raise ValueError("length annotation left flank must precede and not overlap measured regions")
    if measured and right is not None and right.start < max(interval.end for interval in measured):
        raise ValueError("length annotation right flank must follow and not overlap measured regions")


def _validate_boundary_geometry(
    core: tuple[Interval, ...] | None,
    invariant: tuple[Interval, ...] | None,
    array: Interval | None,
    array_only_bp: int | None,
    target_only_bp: int | None,
) -> None:
    if core is None or invariant is None or array is None:
        if array_only_bp is not None or target_only_bp is not None:
            raise ValueError("length annotation boundary geometry must be null when target or ARRAY is unavailable")
        return
    target = (*core, *invariant)
    intersection_bp = sum(_intersection_length(interval, array) for interval in target)
    expected_array_only = array.length_bp - intersection_bp
    expected_target_only = sum(interval.length_bp for interval in target) - intersection_bp
    if (array_only_bp, target_only_bp) != (expected_array_only, expected_target_only):
        raise ValueError("length annotation boundary geometry does not match ARRAY and target intervals")


def _overlap_length(left: tuple[Interval, ...], right: tuple[Interval, ...]) -> int:
    return sum(
        _intersection_length(left_interval, right_interval) for left_interval in left for right_interval in right
    )


def _intersection_length(left: Interval, right: Interval) -> int:
    return max(0, min(left.end, right.end) - max(left.start, right.start))


def _optional_nonnegative_int(value: object, label: str) -> int | None:
    if value is None:
        return None
    if not _is_int(value) or value < 0:
        raise ValueError(f"length annotation {label} must be null or a non-negative integer")
    return value


def _strict_bool(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"length annotation {label} must be boolean")
    return value


def _optional_digest(value: object, label: str) -> str | None:
    if value is None:
        return None
    return _digest(value, label)


def _digest(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != _SHA256_LENGTH
        or any(char not in "0123456789abcdef" for char in value)
    ):
        raise ValueError(f"length annotation {label} must be a lowercase SHA-256 digest")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a non-empty trimmed string")
    return value


def _exact_object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    return value


def _encode_interval(interval: Interval | None) -> dict[str, int] | None:
    if interval is None:
        return None
    return {"start": interval.start, "end": interval.end}


def _encode_interval_list(intervals: tuple[Interval, ...] | None) -> list[dict[str, int]] | None:
    if intervals is None:
        return None
    return [{"start": interval.start, "end": interval.end} for interval in intervals]


def _is_int(value: object) -> TypeGuard[int]:
    return isinstance(value, int) and not isinstance(value, bool)
