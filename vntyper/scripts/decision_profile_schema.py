"""Closed schema and leaf inventory for complete VNtyper decision profiles."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import TypeGuard

from vntyper.scripts.calibration_caller_policy import (
    ADVNTR_CALLER_POLICY_POINTERS,
    KESTREL_CALLER_POLICY_POINTERS,
    caller_policy_values_document,
    decode_caller_policy_values,
)
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
_CALLER_GENERATED_METADATA_FIELDS = _GENERATED_METADATA_FIELDS | {
    "caller_policy_sha256",
    "generated_pointers",
    "generation_target",
    "required_callers",
}
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_NUMERIC_FIELD_KEYS = {"class", "value", "unit", "comparator", "inclusive"}
_NONNUMERIC_FIELD_KEYS = {"class", "value"}
_POINTER_TOKEN_RE = re.compile(r"(?:[^~/]|~[01])+\Z")

_CRITICAL_NUMERIC_METADATA: dict[str, tuple[object, str, str, bool]] = {
    "/components/kestrel/confidence_assignment/reporting_floor": (
        0.00469,
        "depth-score-ratio",
        "gte",
        True,
    ),
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

CALLER_ADVNTR_FIELD_METADATA: dict[str, tuple[str | None, str | None, bool | None]] = {
    "/components/advntr/calibrated_calling/adapter_filter": (None, None, None),
    "/components/advntr/calibrated_calling/cutoff": ("probability", "lt", False),
    "/components/advntr/calibrated_calling/minimum_read_match_ratio": ("matched-read-fraction", "gte", True),
    "/components/advntr/calibrated_calling/minimum_read_support": ("sequencing-reads", "gte", True),
    "/components/advntr/calibrated_calling/mode": (None, None, None),
    "/components/advntr/calibrated_calling/prune_reverse": (None, None, None),
    "/components/advntr/calibrated_calling/rare_unit_fraction": (
        "eligible-repeat-unit-fraction",
        "gte-when-enabled",
        True,
    ),
}
_NULLABLE_NUMERIC_POINTERS = frozenset({"/components/advntr/calibrated_calling/rare_unit_fraction"})


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
    numeric = _is_number(value) or (value is None and pointer in _NULLABLE_NUMERIC_POINTERS)
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


def _validate_profile_metadata(profile: Mapping[str, object]) -> tuple[str, int]:
    _require_exact_fields(profile, _ROOT_FIELDS, label="decision profile root")
    schema_version = profile["schema_version"]
    if isinstance(schema_version, bool) or not isinstance(schema_version, int) or schema_version not in {1, 2}:
        raise ValueError("decision profile schema_version must be 1 or 2")
    for key in ("profile_id", "profile_revision"):
        if not isinstance(profile[key], str) or not profile[key]:
            raise ValueError(f"decision profile {key} must be a non-empty string")
    kind = profile["profile_kind"]
    if not isinstance(kind, str) or kind not in _PROFILE_KINDS:
        raise ValueError(f"decision profile kind is unsupported: {kind!r}")
    generated_metadata = profile["generated_metadata"]
    if kind != "generated":
        if schema_version != 1:
            raise ValueError("decision profile schema_version 2 is reserved for generated caller profiles")
        if generated_metadata is not None:
            raise ValueError("only a generated profile may carry generated_metadata")
        return kind, schema_version
    if not isinstance(generated_metadata, Mapping):
        raise ValueError("generated profile requires generated_metadata")
    expected_metadata = _GENERATED_METADATA_FIELDS if schema_version == 1 else _CALLER_GENERATED_METADATA_FIELDS
    _require_exact_fields(generated_metadata, expected_metadata, label="generated_metadata")
    hash_fields = ["packaged_base_hash", "dataset_manifest_hash", "partition_manifest_hash"]
    if schema_version == 2:
        hash_fields.append("caller_policy_sha256")
    for hash_field in hash_fields:
        value = generated_metadata[hash_field]
        if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
            raise ValueError(f"generated_metadata.{hash_field} must be a lowercase SHA-256 digest")
    for text_field in ("generator_name", "generator_version", "objective"):
        if not isinstance(generated_metadata[text_field], str) or not generated_metadata[text_field]:
            raise ValueError(f"generated_metadata.{text_field} must be a non-empty string")
    seed = generated_metadata["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("generated_metadata.seed must be a non-negative integer")
    if schema_version == 2:
        if generated_metadata["generation_target"] != "callers":
            raise ValueError("generated_metadata.generation_target must be callers")
        if generated_metadata["objective"] != "caller-safety-v1":
            raise ValueError("generated_metadata.objective must be caller-safety-v1")
        required_callers = generated_metadata["required_callers"]
        if required_callers not in (["kestrel"], ["advntr", "kestrel"]):
            raise ValueError("generated_metadata.required_callers has unsupported caller membership")
        generated_pointers = generated_metadata["generated_pointers"]
        if not isinstance(generated_pointers, list) or any(not isinstance(item, str) for item in generated_pointers):
            raise ValueError("generated_metadata.generated_pointers must be an array of strings")
    return kind, schema_version


def _validate_critical_fields(
    fields: Mapping[str, DecisionField], *, revision: str = "2", caller_generated: bool = False
) -> None:
    expected_critical = dict(_CRITICAL_NUMERIC_METADATA)
    if revision == "1":
        expected_critical.pop("/components/kestrel/confidence_assignment/reporting_floor", None)
    for pointer, expected in expected_critical.items():
        field = fields.get(pointer)
        if field is None:
            raise ValueError(f"decision profile is missing critical fixed-safety field {pointer}")
        actual = (field.value, field.unit, field.comparator, field.inclusive)
        if caller_generated and pointer in KESTREL_CALLER_POLICY_POINTERS:
            if field.validation_class is not ValidationClass.GENERATED_MUTABLE or actual[1:] != expected[1:]:
                raise ValueError(f"decision profile caller field semantics differ: {pointer}")
            continue
        if field.validation_class is not ValidationClass.FIXED_SAFETY or actual != expected:
            raise ValueError(f"decision profile critical fixed-safety field differs: {pointer}")
    if revision == "1":
        low = fields["/components/kestrel/confidence_assignment/depth_score_thresholds/low"]
        gg = fields["/components/kestrel/alt_filtering/gg_depth_score_threshold"]
        if low.value != gg.value:
            raise ValueError("independent GG depth-score minimum must equal the reporting floor")
    elif not caller_generated:
        floor = fields["/components/kestrel/confidence_assignment/reporting_floor"]
        gg = fields["/components/kestrel/alt_filtering/gg_depth_score_threshold"]
        if floor.value != gg.value:
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
    kind, schema_version = _validate_profile_metadata(profile)
    caller_generated = kind == "generated" and schema_version == 2
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

    if kind == "packaged":
        if packaged_profile is not None:
            raise ValueError("packaged profile validation does not accept a packaged_profile baseline")
    else:
        if packaged_profile is None:
            raise ValueError("custom decision profile validation requires the packaged baseline")
        if packaged_profile.get("profile_kind") != "packaged":
            raise ValueError("custom decision profile baseline must be packaged")
        packaged_fields = {field.pointer: field for field in validate_complete_inventory(packaged_profile)}
        metadata = profile["generated_metadata"] if caller_generated else None
        required_callers = ()
        caller_pointers: set[str] = set()
        if caller_generated:
            assert isinstance(metadata, Mapping)
            raw_required_callers = metadata["required_callers"]
            assert isinstance(raw_required_callers, list)
            required_callers = tuple(raw_required_callers)
            caller_pointers = set(KESTREL_CALLER_POLICY_POINTERS)
            if "advntr" in required_callers:
                caller_pointers.update(ADVNTR_CALLER_POLICY_POINTERS)
        expected_fields = set(packaged_fields) | (
            set(ADVNTR_CALLER_POLICY_POINTERS) if "advntr" in required_callers else set()
        )
        if set(fields) != expected_fields:
            raise ValueError(f"inventory fields differ: expected {sorted(expected_fields)}, got {sorted(fields)}")
        if profile["profile_id"] == packaged_profile["profile_id"]:
            raise ValueError("custom decision profile must use a distinct profile_id")
        if kind == "generated":
            metadata = profile["generated_metadata"]
            assert isinstance(metadata, Mapping)
            if metadata["packaged_base_hash"] != canonical_sha256(packaged_profile):
                raise ValueError("generated profile packaged_base_hash does not match the packaged profile")
        for pointer, field in fields.items():
            if pointer in ADVNTR_CALLER_POLICY_POINTERS:
                if not caller_generated or pointer not in caller_pointers:
                    raise ValueError(f"decision field is outside the caller-generated capability: {pointer}")
                expected_semantics = CALLER_ADVNTR_FIELD_METADATA[pointer]
                if field.validation_class is not ValidationClass.GENERATED_MUTABLE:
                    raise ValueError(f"decision field class differs from caller-generated policy: {pointer}")
                if (field.unit, field.comparator, field.inclusive) != expected_semantics:
                    raise ValueError(f"decision field semantics differ from caller-generated policy: {pointer}")
                continue
            baseline = packaged_fields[pointer]
            caller_mutable = caller_generated and pointer in KESTREL_CALLER_POLICY_POINTERS
            expected_class = ValidationClass.GENERATED_MUTABLE if caller_mutable else baseline.validation_class
            if field.validation_class is not expected_class:
                raise ValueError(f"decision field class differs from packaged profile: {pointer}")
            if (field.unit, field.comparator, field.inclusive) != (
                baseline.unit,
                baseline.comparator,
                baseline.inclusive,
            ):
                raise ValueError(f"decision field semantics differ from packaged profile: {pointer}")
            if not _same_json_type(field.value, baseline.value):
                raise ValueError(f"decision field type differs from packaged profile: {pointer}")
            if (
                not caller_mutable
                and field.validation_class is ValidationClass.FIXED_SAFETY
                and field.value != baseline.value
            ):
                raise ValueError(f"immutable fixed-safety field differs: {pointer}")
            if (
                kind == "generated"
                and field.validation_class is ValidationClass.EXPLICIT_CUSTOM
                and field.value != baseline.value
            ):
                raise ValueError(f"generated profile must copy explicit-custom field: {pointer}")

    _validate_critical_fields(
        fields,
        revision=str(profile.get("profile_revision", "2")),
        caller_generated=caller_generated,
    )
    generated_pointers = {
        pointer for pointer, field in fields.items() if field.validation_class is ValidationClass.GENERATED_MUTABLE
    }
    expected_generated_pointers = set(_GENERATED_BOUNDS)
    if caller_generated:
        metadata = profile["generated_metadata"]
        assert isinstance(metadata, Mapping)
        raw_required_callers = metadata["required_callers"]
        assert isinstance(raw_required_callers, list)
        caller_pointers = set(KESTREL_CALLER_POLICY_POINTERS)
        if "advntr" in raw_required_callers:
            caller_pointers.update(ADVNTR_CALLER_POLICY_POINTERS)
        if metadata["generated_pointers"] != sorted(caller_pointers):
            raise ValueError("generated_metadata.generated_pointers differs from the exact caller capability")
        policy_document = {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": list(raw_required_callers),
            "values": {pointer: fields[pointer].value for pointer in sorted(caller_pointers)},
        }
        policy = decode_caller_policy_values(policy_document)
        if metadata["caller_policy_sha256"] != policy.sha256:
            raise ValueError("generated profile caller policy digest does not match its inventory")
        if caller_policy_values_document(policy) != policy_document:
            raise ValueError("generated profile caller policy does not round-trip canonically")
        expected_generated_pointers.update(caller_pointers)
    if generated_pointers != expected_generated_pointers:
        raise ValueError(
            "generated-mutable fields differ: "
            f"expected {sorted(expected_generated_pointers)}, got {sorted(generated_pointers)}"
        )
    for pointer in set(_GENERATED_BOUNDS):
        _validate_generated_value(fields[pointer])
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
