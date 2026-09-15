"""Pure exact paired truth targets and explicit boundary conversions."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.calibration_intake_contract import LengthTruth, TruthRecord
from vntyper.scripts.canonical_json import canonical_sha256

ConversionSourceUnit = Literal["bp", "repeat-count"]

_ROOT_FIELDS = {"schema_version", "conversions"}
_COMMON_FIELDS = {
    "conversion_id",
    "source_unit",
    "source_boundary_definition",
    "target_boundary_definition",
    "repeat_unit_bp",
}
_BP_FIELDS = _COMMON_FIELDS | {"allele_1_nonrepeat_bp", "allele_2_nonrepeat_bp"}
_COUNT_FIELDS = _COMMON_FIELDS | {"allele_1_offset", "allele_2_offset"}


@dataclass(frozen=True)
class LengthConversion:
    """One registered exact conversion into a target repeat-count boundary."""

    conversion_id: str
    source_unit: ConversionSourceUnit
    source_boundary_definition: str
    target_boundary_definition: str
    repeat_unit_bp: int
    allele_1_subtract: Fraction
    allele_2_subtract: Fraction


@dataclass(frozen=True)
class ConversionRegistry:
    """Immutable explicit conversion registry and its canonical digest."""

    conversions: Mapping[str, LengthConversion]
    sha256: str


def decode_conversion_registry(value: object) -> ConversionRegistry:
    """Decode a strict length-boundary conversion registry.

    Args:
        value: Parsed conversion registry JSON value.

    Returns:
        An immutable registry keyed by explicit conversion ID.

    Raises:
        ValueError: If fields, values, or conversion IDs are invalid.
    """
    root = _exact_object(value, _ROOT_FIELDS, "calibration conversion registry")
    if root["schema_version"] != "calibration-length-conversions-v1":
        raise ValueError("calibration conversion registry schema version must be calibration-length-conversions-v1")
    raw_conversions = root["conversions"]
    if not isinstance(raw_conversions, list) or not raw_conversions:
        raise ValueError("calibration conversion registry conversions must be a non-empty list")
    decoded = tuple(_decode_conversion(row) for row in raw_conversions)
    ids = tuple(row.conversion_id for row in decoded)
    if len(ids) != len(set(ids)):
        raise ValueError("calibration conversion ids must be unique")
    conversions = MappingProxyType(
        {row.conversion_id: row for row in sorted(decoded, key=lambda row: row.conversion_id)}
    )
    normalized = {
        "schema_version": "calibration-length-conversions-v1",
        "conversions": [_encode_conversion(row) for row in conversions.values()],
    }
    return ConversionRegistry(conversions, canonical_sha256(normalized))


def length_target(truth: TruthRecord) -> float | None:
    """Return a direct exact paired repeat-count target when one is available.

    Args:
        truth: Validated intake truth record.

    Returns:
        Total diploid repeat count, or ``None`` for unavailable fitting truth.

    Raises:
        ValueError: If the record requires an explicit unit or boundary conversion.
    """
    length = _eligible_exact_pair(truth)
    if length is None:
        return None
    if length.unit != "repeat-count" or length.conversion_id is not None:
        raise ValueError("calibration length truth requires an explicit registered conversion")
    allele_1, allele_2 = _required_pair(length)
    return float(allele_1 + allele_2)


def converted_length_target(
    truth: TruthRecord,
    registry: ConversionRegistry,
    *,
    target_boundary_definition: str,
) -> float | None:
    """Return an exact paired target through one explicitly registered conversion.

    Args:
        truth: Validated intake truth record.
        registry: Separately decoded and hash-bound conversion registry.
        target_boundary_definition: Required target boundary identity.

    Returns:
        Total diploid repeat count, or ``None`` for unavailable fitting truth.

    Raises:
        ValueError: If registry bindings mismatch or conversion is non-integral.
    """
    if not isinstance(registry, ConversionRegistry):
        raise ValueError("calibration conversions must be a ConversionRegistry")
    if not isinstance(target_boundary_definition, str) or not target_boundary_definition:
        raise ValueError("calibration target boundary definition must be a non-empty string")
    length = _eligible_exact_pair(truth)
    if length is None:
        return None
    if length.conversion_id is None:
        if length.unit == "repeat-count" and length.boundary_definition == target_boundary_definition:
            return length_target(truth)
        raise ValueError("calibration length boundary mismatch requires a registered conversion")
    conversion = registry.conversions.get(length.conversion_id)
    if conversion is None:
        raise ValueError(f"unknown calibration length conversion id: {length.conversion_id}")
    if conversion.source_unit != length.unit:
        raise ValueError("calibration length conversion source unit does not match truth")
    if conversion.source_boundary_definition != length.boundary_definition:
        raise ValueError("calibration length conversion source boundary does not match truth")
    if conversion.target_boundary_definition != target_boundary_definition:
        raise ValueError("calibration length conversion target boundary does not match requested target boundary")
    if conversion.repeat_unit_bp != length.repeat_unit_bp:
        raise ValueError("calibration length conversion repeat unit does not match truth")
    allele_1, allele_2 = _required_pair(length)
    converted_1 = allele_1 - conversion.allele_1_subtract
    converted_2 = allele_2 - conversion.allele_2_subtract
    if conversion.source_unit == "bp":
        converted_1 /= conversion.repeat_unit_bp
        converted_2 /= conversion.repeat_unit_bp
    if converted_1 <= 0 or converted_2 <= 0:
        raise ValueError("calibration length conversion must produce positive allele counts")
    if converted_1.denominator != 1 or converted_2.denominator != 1:
        raise ValueError("calibration length conversion must produce integral exact repeat counts")
    return float(converted_1 + converted_2)


def _eligible_exact_pair(truth: TruthRecord) -> LengthTruth | None:
    if not isinstance(truth, TruthRecord):
        raise ValueError("calibration truth must be a TruthRecord")
    length = truth.length
    if truth.status != "confirmed" or length is None or length.measurement != "exact":
        return None
    if length.allele_1 is None or length.allele_2 is None:
        return None
    return length


def _required_pair(length: LengthTruth) -> tuple[Fraction, Fraction]:
    if length.allele_1 is None or length.allele_2 is None:
        raise ValueError("calibration length truth does not contain an exact pair")
    return length.allele_1, length.allele_2


def _decode_conversion(value: object) -> LengthConversion:
    if not isinstance(value, Mapping):
        raise ValueError("calibration length conversion must be an object")
    source_unit = value.get("source_unit")
    if source_unit == "bp":
        raw = _exact_object(value, _BP_FIELDS, "calibration bp conversion")
        repeat_unit_bp = _positive_integer(raw["repeat_unit_bp"], "calibration conversion repeat unit bp")
        subtract_1 = _nonnegative_number(raw["allele_1_nonrepeat_bp"], "calibration allele 1 nonrepeat bp")
        subtract_2 = _nonnegative_number(raw["allele_2_nonrepeat_bp"], "calibration allele 2 nonrepeat bp")
    elif source_unit == "repeat-count":
        raw = _exact_object(value, _COUNT_FIELDS, "calibration repeat-count conversion")
        repeat_unit_bp = _positive_integer(raw["repeat_unit_bp"], "calibration conversion repeat unit bp")
        subtract_1 = _nonnegative_number(raw["allele_1_offset"], "calibration allele 1 offset")
        subtract_2 = _nonnegative_number(raw["allele_2_offset"], "calibration allele 2 offset")
        if subtract_1.denominator != 1 or subtract_2.denominator != 1:
            raise ValueError("calibration repeat-count conversion offsets must be integral")
    else:
        raise ValueError(f"unsupported calibration conversion source unit: {source_unit!r}")
    return LengthConversion(
        _nonempty_string(raw["conversion_id"], "calibration conversion id"),
        cast(ConversionSourceUnit, source_unit),
        _nonempty_string(raw["source_boundary_definition"], "calibration conversion source boundary"),
        _nonempty_string(raw["target_boundary_definition"], "calibration conversion target boundary"),
        repeat_unit_bp,
        subtract_1,
        subtract_2,
    )


def _encode_conversion(conversion: LengthConversion) -> dict[str, object]:
    common: dict[str, object] = {
        "conversion_id": conversion.conversion_id,
        "source_unit": conversion.source_unit,
        "source_boundary_definition": conversion.source_boundary_definition,
        "target_boundary_definition": conversion.target_boundary_definition,
    }
    if conversion.source_unit == "bp":
        return {
            **common,
            "repeat_unit_bp": conversion.repeat_unit_bp,
            "allele_1_nonrepeat_bp": _json_number(conversion.allele_1_subtract),
            "allele_2_nonrepeat_bp": _json_number(conversion.allele_2_subtract),
        }
    return {
        **common,
        "repeat_unit_bp": conversion.repeat_unit_bp,
        "allele_1_offset": _json_number(conversion.allele_1_subtract),
        "allele_2_offset": _json_number(conversion.allele_2_subtract),
    }


def _exact_object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    return value


def _nonempty_string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a non-empty string")
    return value


def _positive_integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _nonnegative_number(value: object, label: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f"{label} must be a finite nonnegative number")
    return Fraction(str(value))


def _json_number(value: Fraction) -> int | float:
    if value.denominator == 1:
        return value.numerator
    return float(value)
