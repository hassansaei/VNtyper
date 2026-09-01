"""Closed schema and leaf inventory for complete VNtyper decision profiles."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import TypeGuard

from vntyper.scripts.canonical_json import canonical_sha256

_ROOT_FIELDS = {
    "schema_version",
    "profile_id",
    "profile_revision",
    "profile_kind",
    "generated_metadata",
    "inventory",
}
_PROFILE_KINDS = frozenset({"packaged", "explicit-custom", "generated"})
_COMPONENTS = frozenset({"advntr", "cross_match", "dominance", "kestrel", "nomenclature", "shark"})
_GENERATED_METADATA_FIELDS = {
    "packaged_base_hash",
    "generator_name",
    "generator_version",
    "objective",
    "dataset_manifest_hash",
    "partition_manifest_hash",
    "seed",
}
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_NUMERIC_FIELD_KEYS = {"class", "value", "unit", "comparator", "inclusive"}
_NONNUMERIC_FIELD_KEYS = {"class", "value"}
_POINTER_TOKEN_RE = re.compile(r"(?:[^~/]|~[01])+\Z")

_CRITICAL_NUMERIC_METADATA: dict[str, tuple[object, str, str, bool]] = {
    "/components/kestrel/confidence_assignment/depth_score_thresholds/low": (
        0.00469,
        "depth-score-ratio",
        "gte",
        True,
    ),
    "/components/kestrel/alt_filtering/gg_depth_score_threshold": (
        0.00469,
        "depth-score-ratio",
        "gte",
        True,
    ),
    "/components/kestrel/confidence_assignment/depth_score_thresholds/high": (
        0.00515,
        "depth-score-ratio",
        "gte",
        True,
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": (
        20,
        "alternate-kmer-path-depth",
        "lte",
        True,
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": (
        21,
        "alternate-kmer-path-depth",
        "gte",
        True,
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": (
        100,
        "alternate-kmer-path-depth",
        "gte-and-upper-exclusive",
        False,
    ),
    "/components/kestrel/confidence_assignment/var_active_region_threshold": (
        200,
        "active-region-kmer-depth",
        "lte",
        True,
    ),
    "/components/nomenclature/thresholds/bam_flank": (8, "base-pairs-per-side", "eq", True),
    "/components/nomenclature/thresholds/bam_thin_haplotype_record_support": (
        3,
        "resolved-haplotype-records",
        "lt",
        False,
    ),
    "/components/nomenclature/identity_reconciliation/kestrel_min_alternate_kmer_path_depth": (
        5,
        "alternate-kmer-path-depth",
        "gte",
        True,
    ),
    "/components/nomenclature/identity_reconciliation/advntr_min_sequencing_read_support": (
        5,
        "sequencing-reads",
        "gte",
        True,
    ),
}

_GENERATED_BOUNDS: dict[str, tuple[type[object], object]] = {
    "/components/dominance/enabled": (bool, None),
    "/components/dominance/minimum_record_count_margin": (int, (0, None)),
    "/components/dominance/minimum_record_share": (float, (0.0, 1.0)),
    "/components/dominance/minimum_record_share_margin": (float, (0.0, 1.0)),
    "/components/dominance/xd_veto": (
        str,
        frozenset({"disabled", "missingness", "concentration", "discordance"}),
    ),
    "/components/dominance/abstain_on_inadmissible_advntr": (bool, None),
}


class ValidationClass(str, Enum):
    """How one complete-profile decision leaf may vary."""

    FIXED_SAFETY = "fixed-safety"
    EXPLICIT_CUSTOM = "explicit-custom"
    GENERATED_MUTABLE = "generated-mutable"


@dataclass(frozen=True)
class DecisionField:
    """One validated decision leaf and its comparison semantics."""

    pointer: str
    validation_class: ValidationClass
    value: object
    unit: str | None = None
    comparator: str | None = None
    inclusive: bool | None = None


def _require_exact_fields(value: Mapping[str, object], expected: set[str], *, label: str) -> None:
    actual = set(value)
    if actual != expected:
        raise ValueError(f"{label} fields differ: expected {sorted(expected)}, got {sorted(actual)}")


def _decode_pointer(pointer: object) -> tuple[str, ...]:
    if not isinstance(pointer, str) or not pointer.startswith("/"):
        raise ValueError(f"decision inventory key is not a JSON Pointer: {pointer!r}")
    raw_tokens = pointer[1:].split("/")
    if not raw_tokens or any(_POINTER_TOKEN_RE.fullmatch(token) is None for token in raw_tokens):
        raise ValueError(f"decision inventory key is not a canonical JSON Pointer: {pointer!r}")
    return tuple(token.replace("~1", "/").replace("~0", "~") for token in raw_tokens)


def _encode_token(token: str) -> str:
    return token.replace("~", "~0").replace("/", "~1")


def flatten_decision_projection(projection: Mapping[str, object]) -> dict[str, object]:
    """Flatten a component projection into JSON-Pointer-addressed leaf values.

    Empty arrays and objects are retained as leaves. Array indices are retained in
    pointers, so array ordering remains semantic while mapping key order does not.

    Args:
        projection: Complete component mapping without profile metadata.

    Returns:
        A pointer-to-value mapping for every scalar or empty-container leaf.
    """
    flattened: dict[str, object] = {}

    def visit(value: object, tokens: tuple[str, ...]) -> None:
        if isinstance(value, Mapping) and value:
            for key, child in value.items():
                if not isinstance(key, str):
                    raise ValueError("decision projection object keys must be strings")
                visit(child, (*tokens, key))
            return
        if isinstance(value, list) and value:
            for index, child in enumerate(value):
                visit(child, (*tokens, str(index)))
            return
        pointer = "/" + "/".join(_encode_token(token) for token in tokens)
        flattened[pointer] = value

    visit({"components": dict(projection)}, ())
    return flattened


def _is_number(value: object) -> TypeGuard[int | float]:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _parse_field(pointer: str, raw: object) -> DecisionField:
    if not isinstance(raw, Mapping):
        raise ValueError(f"decision field {pointer} must be an object")
    value = raw.get("value")
    numeric = _is_number(value)
    _require_exact_fields(raw, _NUMERIC_FIELD_KEYS if numeric else _NONNUMERIC_FIELD_KEYS, label=pointer)
    raw_class = raw["class"]
    try:
        validation_class = ValidationClass(raw_class)
    except (TypeError, ValueError) as error:
        raise ValueError(f"decision field {pointer} has unsupported class: {raw_class!r}") from error
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"decision field {pointer} must be finite")
    if not numeric:
        return DecisionField(pointer, validation_class, value)
    unit = raw["unit"]
    comparator = raw["comparator"]
    inclusive = raw["inclusive"]
    if not isinstance(unit, str) or not unit:
        raise ValueError(f"numeric decision field {pointer} unit must be a non-empty string")
    if not isinstance(comparator, str) or not comparator:
        raise ValueError(f"numeric decision field {pointer} comparator must be a non-empty string")
    if not isinstance(inclusive, bool):
        raise ValueError(f"numeric decision field {pointer} inclusive must be Boolean")
    return DecisionField(pointer, validation_class, value, unit, comparator, inclusive)


def _validate_profile_metadata(profile: Mapping[str, object]) -> str:
    _require_exact_fields(profile, _ROOT_FIELDS, label="decision profile root")
    if profile["schema_version"] != 1 or isinstance(profile["schema_version"], bool):
        raise ValueError("decision profile schema_version must be 1")
    for key in ("profile_id", "profile_revision"):
        if not isinstance(profile[key], str) or not profile[key]:
            raise ValueError(f"decision profile {key} must be a non-empty string")
    kind = profile["profile_kind"]
    if not isinstance(kind, str) or kind not in _PROFILE_KINDS:
        raise ValueError(f"decision profile kind is unsupported: {kind!r}")
    generated_metadata = profile["generated_metadata"]
    if kind != "generated":
        if generated_metadata is not None:
            raise ValueError("only a generated profile may carry generated_metadata")
        return kind
    if not isinstance(generated_metadata, Mapping):
        raise ValueError("generated profile requires generated_metadata")
    _require_exact_fields(generated_metadata, _GENERATED_METADATA_FIELDS, label="generated_metadata")
    for hash_field in ("packaged_base_hash", "dataset_manifest_hash", "partition_manifest_hash"):
        value = generated_metadata[hash_field]
        if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
            raise ValueError(f"generated_metadata.{hash_field} must be a lowercase SHA-256 digest")
    for text_field in ("generator_name", "generator_version", "objective"):
        if not isinstance(generated_metadata[text_field], str) or not generated_metadata[text_field]:
            raise ValueError(f"generated_metadata.{text_field} must be a non-empty string")
    seed = generated_metadata["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("generated_metadata.seed must be a non-negative integer")
    return kind


def _validate_critical_fields(fields: Mapping[str, DecisionField]) -> None:
    for pointer, expected in _CRITICAL_NUMERIC_METADATA.items():
        field = fields.get(pointer)
        if field is None:
            raise ValueError(f"decision profile is missing critical fixed-safety field {pointer}")
        actual = (field.value, field.unit, field.comparator, field.inclusive)
        if field.validation_class is not ValidationClass.FIXED_SAFETY or actual != expected:
            raise ValueError(f"decision profile critical fixed-safety field differs: {pointer}")
    low = fields["/components/kestrel/confidence_assignment/depth_score_thresholds/low"]
    gg = fields["/components/kestrel/alt_filtering/gg_depth_score_threshold"]
    if low.value != gg.value:
        raise ValueError("independent GG depth-score minimum must equal the reporting floor")


def _same_json_type(left: object, right: object) -> bool:
    if _is_number(left) and _is_number(right):
        return True
    return type(left) is type(right)


def _validate_generated_value(field: DecisionField) -> None:
    expected = _GENERATED_BOUNDS.get(field.pointer)
    if expected is None:
        raise ValueError(f"unsupported generated-mutable field: {field.pointer}")
    expected_type, bounds = expected
    value = field.value
    if expected_type is float:
        if not _is_number(value):
            raise ValueError(f"generated-mutable field {field.pointer} must be numeric")
        numeric_value = float(value)
        assert isinstance(bounds, tuple)
        if not bounds[0] <= numeric_value <= bounds[1]:
            raise ValueError(f"generated-mutable field {field.pointer} is outside its frozen range")
        return
    if expected_type is int:
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"generated-mutable field {field.pointer} must be an integer")
        assert isinstance(bounds, tuple)
        if value < bounds[0]:
            raise ValueError(f"generated-mutable field {field.pointer} is outside its frozen range")
        return
    if not isinstance(value, expected_type):
        raise ValueError(f"generated-mutable field {field.pointer} has the wrong type")
    if isinstance(bounds, frozenset) and value not in bounds:
        raise ValueError(f"generated-mutable field {field.pointer} is outside its frozen enum")


def validate_complete_inventory(
    profile: Mapping[str, object], *, packaged_profile: Mapping[str, object] | None = None
) -> tuple[DecisionField, ...]:
    """Validate a complete decision inventory and its immutable-field contract.

    Args:
        profile: Decoded complete profile.
        packaged_profile: Verified packaged baseline required for explicit or
            generated profiles.

    Returns:
        Validated fields in sorted JSON-Pointer order.

    Raises:
        ValueError: If metadata, field coverage, values, or mutability violate
            the closed profile contract.
    """
    if not isinstance(profile, Mapping):
        raise ValueError("decision profile must be an object")
    kind = _validate_profile_metadata(profile)
    raw_inventory = profile["inventory"]
    if not isinstance(raw_inventory, Mapping) or not raw_inventory:
        raise ValueError("decision profile inventory must be a non-empty object")
    fields: dict[str, DecisionField] = {}
    for raw_pointer, raw_field in raw_inventory.items():
        tokens = _decode_pointer(raw_pointer)
        if len(tokens) < 2 or tokens[0] != "components" or tokens[1] not in _COMPONENTS:
            raise ValueError(f"decision inventory field is outside the closed components: {raw_pointer}")
        assert isinstance(raw_pointer, str)
        fields[raw_pointer] = _parse_field(raw_pointer, raw_field)
    _validate_critical_fields(fields)
    generated_pointers = {
        pointer for pointer, field in fields.items() if field.validation_class is ValidationClass.GENERATED_MUTABLE
    }
    if generated_pointers != set(_GENERATED_BOUNDS):
        raise ValueError(
            f"generated-mutable fields differ: expected {sorted(_GENERATED_BOUNDS)}, got {sorted(generated_pointers)}"
        )
    for pointer in generated_pointers:
        _validate_generated_value(fields[pointer])

    if kind == "packaged":
        if packaged_profile is not None:
            raise ValueError("packaged profile validation does not accept a packaged_profile baseline")
    else:
        if packaged_profile is None:
            raise ValueError("custom decision profile validation requires the packaged baseline")
        if packaged_profile.get("profile_kind") != "packaged":
            raise ValueError("custom decision profile baseline must be packaged")
        packaged_fields = {field.pointer: field for field in validate_complete_inventory(packaged_profile)}
        if set(fields) != set(packaged_fields):
            raise ValueError(f"inventory fields differ: expected {sorted(packaged_fields)}, got {sorted(fields)}")
        if profile["profile_id"] == packaged_profile["profile_id"]:
            raise ValueError("custom decision profile must use a distinct profile_id")
        if kind == "generated":
            metadata = profile["generated_metadata"]
            assert isinstance(metadata, Mapping)
            if metadata["packaged_base_hash"] != canonical_sha256(packaged_profile):
                raise ValueError("generated profile packaged_base_hash does not match the packaged profile")
        for pointer, field in fields.items():
            baseline = packaged_fields[pointer]
            if field.validation_class is not baseline.validation_class:
                raise ValueError(f"decision field class differs from packaged profile: {pointer}")
            if (field.unit, field.comparator, field.inclusive) != (
                baseline.unit,
                baseline.comparator,
                baseline.inclusive,
            ):
                raise ValueError(f"decision field semantics differ from packaged profile: {pointer}")
            if not _same_json_type(field.value, baseline.value):
                raise ValueError(f"decision field type differs from packaged profile: {pointer}")
            if field.validation_class is ValidationClass.FIXED_SAFETY and field.value != baseline.value:
                raise ValueError(f"immutable fixed-safety field differs: {pointer}")
            if (
                kind == "generated"
                and field.validation_class is ValidationClass.EXPLICIT_CUSTOM
                and field.value != baseline.value
            ):
                raise ValueError(f"generated profile must copy explicit-custom field: {pointer}")
    components = _projection_from_fields(fields)
    if set(components) != _COMPONENTS:
        raise ValueError(
            f"decision profile components differ: expected {sorted(_COMPONENTS)}, got {sorted(components)}"
        )
    from vntyper.scripts.decision_profile_semantics import validate_component_semantics

    validate_component_semantics(components)
    return tuple(fields[pointer] for pointer in sorted(fields))


def _projection_from_fields(fields: Mapping[str, DecisionField]) -> dict[str, object]:
    children_by_prefix: dict[tuple[str, ...], set[str]] = {}
    decoded: dict[str, tuple[str, ...]] = {}
    for pointer in fields:
        tokens = _decode_pointer(pointer)
        decoded[pointer] = tokens
        for index, token in enumerate(tokens):
            children_by_prefix.setdefault(tokens[:index], set()).add(token)

    def new_container(prefix: tuple[str, ...]) -> object:
        children = children_by_prefix.get(prefix, set())
        if children and all(token.isdigit() for token in children):
            indices = sorted(int(token) for token in children)
            if indices == list(range(len(indices))):
                return []
        return {}

    root: dict[str, object] = {}
    for pointer in sorted(fields):
        tokens = decoded[pointer]
        current: object = root
        for index, token in enumerate(tokens):
            final = index == len(tokens) - 1
            if isinstance(current, dict):
                if final:
                    current[token] = fields[pointer].value
                else:
                    current = current.setdefault(token, new_container(tokens[: index + 1]))
            elif isinstance(current, list):
                if not token.isdigit():
                    raise ValueError(f"decision pointer uses a non-index token inside an array: {pointer}")
                array_index = int(token)
                while len(current) <= array_index:
                    current.append(None)
                if final:
                    current[array_index] = fields[pointer].value
                else:
                    if current[array_index] is None:
                        current[array_index] = new_container(tokens[: index + 1])
                    current = current[array_index]
            else:
                raise ValueError(f"decision pointers conflict at {pointer}")
    components = root.get("components")
    if not isinstance(components, dict):
        raise ValueError("decision inventory does not reconstruct a components object")
    return components


def component_projection(
    profile: Mapping[str, object],
    component: str,
    *,
    packaged_profile: Mapping[str, object] | None = None,
) -> object:
    """Return one reconstructed immutable decision component.

    Args:
        profile: Validated or packaged complete decision profile.
        component: Closed component name.
        packaged_profile: Verified packaged baseline for a custom profile.

    Returns:
        The component's JSON-compatible decision mapping.

    Raises:
        ValueError: If the profile or component is invalid.
    """
    if component not in _COMPONENTS:
        raise ValueError(f"unsupported decision component: {component}")
    fields = {field.pointer: field for field in validate_complete_inventory(profile, packaged_profile=packaged_profile)}
    return _projection_from_fields(fields)[component]
