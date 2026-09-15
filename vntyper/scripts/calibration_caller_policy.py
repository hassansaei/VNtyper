"""Immutable complete caller-policy values keyed by binding JSON pointers."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, NoReturn, TypeAlias, cast

from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

CallerName = Literal["advntr", "kestrel"]
CallerPolicyScalar: TypeAlias = str | bool | int | float | None

KESTREL_CALLER_POLICY_POINTERS = tuple(
    sorted(
        (
            "/components/kestrel/alt_filtering/gg_depth_score_threshold",
            "/components/kestrel/confidence_assignment/reporting_floor",
            "/components/kestrel/confidence_assignment/var_active_region_threshold",
            "/components/kestrel/confidence_assignment/depth_score_thresholds/low",
            "/components/kestrel/confidence_assignment/depth_score_thresholds/high",
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/low",
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low",
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high",
        )
    )
)
ADVNTR_CALLER_POLICY_POINTERS = tuple(
    sorted(
        (
            "/components/advntr/calibrated_calling/mode",
            "/components/advntr/calibrated_calling/cutoff",
            "/components/advntr/calibrated_calling/minimum_read_support",
            "/components/advntr/calibrated_calling/rare_unit_fraction",
            "/components/advntr/calibrated_calling/adapter_filter",
            "/components/advntr/calibrated_calling/minimum_read_match_ratio",
            "/components/advntr/calibrated_calling/prune_reverse",
        )
    )
)

_ROOT_FIELDS = {"schema_version", "required_callers", "values"}
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))

_GG_DEPTH = "/components/kestrel/alt_filtering/gg_depth_score_threshold"
_REPORTING_FLOOR = "/components/kestrel/confidence_assignment/reporting_floor"
_ACTIVE_DEPTH = "/components/kestrel/confidence_assignment/var_active_region_threshold"
_DEPTH_LOW = "/components/kestrel/confidence_assignment/depth_score_thresholds/low"
_DEPTH_HIGH = "/components/kestrel/confidence_assignment/depth_score_thresholds/high"
_ALT_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/low"
_ALT_MID_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low"
_ALT_MID_HIGH = "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high"

_ADVNTR_MODE = "/components/advntr/calibrated_calling/mode"
_ADVNTR_CUTOFF = "/components/advntr/calibrated_calling/cutoff"
_ADVNTR_SUPPORT = "/components/advntr/calibrated_calling/minimum_read_support"
_ADVNTR_RARE = "/components/advntr/calibrated_calling/rare_unit_fraction"
_ADVNTR_ADAPTER = "/components/advntr/calibrated_calling/adapter_filter"
_ADVNTR_MATCH = "/components/advntr/calibrated_calling/minimum_read_match_ratio"
_ADVNTR_PRUNE = "/components/advntr/calibrated_calling/prune_reverse"


@dataclass(frozen=True)
class CallerPolicyValues:
    """One complete policy-value inventory without transported model/build assets."""

    required_callers: tuple[CallerName, ...]
    values: Mapping[str, CallerPolicyScalar]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"{label} fields differ from the closed contract")
    return value


def _callers(value: object) -> tuple[CallerName, ...]:
    if value == ["kestrel"]:
        return ("kestrel",)
    if value == ["advntr", "kestrel"]:
        return ("advntr", "kestrel")
    _fail("required_callers must be ['kestrel'] or sorted ['advntr', 'kestrel']")


def _unit_number(value: object, pointer: str) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        _fail(f"{pointer} must be a finite number between zero and one")
    if not 0 <= value <= 1:
        _fail(f"{pointer} must be between zero and one")
    return value


def _nonnegative_integer(value: object, pointer: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        _fail(f"{pointer} must be a nonnegative integer")
    return value


def _positive_integer(value: object, pointer: str) -> int:
    result = _nonnegative_integer(value, pointer)
    if result < 1:
        _fail(f"{pointer} must be a positive integer")
    return result


def _bounded_float(value: object, pointer: str, *, include_one: bool) -> float:
    if not isinstance(value, float) or not math.isfinite(value):
        _fail(f"{pointer} must be a finite float")
    if value <= 0 or (value > 1 if include_one else value >= 1):
        boundary = "greater than zero through one" if include_one else "strictly between zero and one"
        _fail(f"{pointer} must be a finite float {boundary}")
    return value


def _decode_values(value: object, callers: tuple[CallerName, ...]) -> Mapping[str, CallerPolicyScalar]:
    expected = set(KESTREL_CALLER_POLICY_POINTERS)
    if "advntr" in callers:
        expected.update(ADVNTR_CALLER_POLICY_POINTERS)
    raw = _object(value, expected, "caller policy pointers")
    parsed: dict[str, CallerPolicyScalar] = {
        _GG_DEPTH: _unit_number(raw[_GG_DEPTH], _GG_DEPTH),
        _REPORTING_FLOOR: _unit_number(raw[_REPORTING_FLOOR], _REPORTING_FLOOR),
        _ACTIVE_DEPTH: _nonnegative_integer(raw[_ACTIVE_DEPTH], _ACTIVE_DEPTH),
        _DEPTH_LOW: _unit_number(raw[_DEPTH_LOW], _DEPTH_LOW),
        _DEPTH_HIGH: _unit_number(raw[_DEPTH_HIGH], _DEPTH_HIGH),
        _ALT_LOW: _nonnegative_integer(raw[_ALT_LOW], _ALT_LOW),
        _ALT_MID_LOW: _nonnegative_integer(raw[_ALT_MID_LOW], _ALT_MID_LOW),
        _ALT_MID_HIGH: _nonnegative_integer(raw[_ALT_MID_HIGH], _ALT_MID_HIGH),
    }
    depth_low = cast(int | float, parsed[_DEPTH_LOW])
    depth_high = cast(int | float, parsed[_DEPTH_HIGH])
    if depth_low > depth_high:
        _fail("Kestrel depth-score low must not exceed high")
    alt_low = cast(int, parsed[_ALT_LOW])
    alt_mid_low = cast(int, parsed[_ALT_MID_LOW])
    alt_mid_high = cast(int, parsed[_ALT_MID_HIGH])
    if alt_mid_low != alt_low + 1 or alt_mid_high < alt_mid_low + 1:
        _fail("Kestrel alternate-depth partition must satisfy mid_low=low+1 and mid_high>=mid_low+1")
    if "advntr" in callers:
        mode = raw[_ADVNTR_MODE]
        if mode not in {"legacy", "exact"} or not isinstance(mode, str):
            _fail(f"{_ADVNTR_MODE} must be mode legacy or exact")
        rare = raw[_ADVNTR_RARE]
        parsed.update(
            {
                _ADVNTR_MODE: mode,
                _ADVNTR_CUTOFF: _bounded_float(raw[_ADVNTR_CUTOFF], _ADVNTR_CUTOFF, include_one=False),
                _ADVNTR_SUPPORT: _positive_integer(raw[_ADVNTR_SUPPORT], _ADVNTR_SUPPORT),
                _ADVNTR_RARE: (None if rare is None else _bounded_float(rare, _ADVNTR_RARE, include_one=True)),
                _ADVNTR_ADAPTER: _boolean(raw[_ADVNTR_ADAPTER], _ADVNTR_ADAPTER),
                _ADVNTR_MATCH: _bounded_float(raw[_ADVNTR_MATCH], _ADVNTR_MATCH, include_one=True),
                _ADVNTR_PRUNE: _boolean(raw[_ADVNTR_PRUNE], _ADVNTR_PRUNE),
            }
        )
    return MappingProxyType(dict(sorted(parsed.items())))


def _boolean(value: object, pointer: str) -> bool:
    if not isinstance(value, bool):
        _fail(f"{pointer} must be boolean")
    return value


def _document(policy: CallerPolicyValues) -> dict[str, object]:
    return {
        "schema_version": "calibration-caller-policy-values-v1",
        "required_callers": list(policy.required_callers),
        "values": dict(policy.values),
    }


def decode_caller_policy_values(value: object) -> CallerPolicyValues:
    """Decode one complete finite caller policy without bundle transport assets.

    Args:
        value: Closed caller policy value document keyed by exact JSON pointers.

    Returns:
        Immutable validated values and their canonical digest.

    Raises:
        ValueError: If fields, pointer membership, types or numeric relationships are invalid.
    """
    raw = _object(value, _ROOT_FIELDS, "caller policy")
    if raw["schema_version"] != "calibration-caller-policy-values-v1":
        _fail("caller policy schema_version must be calibration-caller-policy-values-v1")
    callers = _callers(raw["required_callers"])
    values = _decode_values(raw["values"], callers)
    policy = CallerPolicyValues(callers, values, "")
    return CallerPolicyValues(callers, values, canonical_sha256(_document(policy)))


def _require_policy(policy: CallerPolicyValues) -> CallerPolicyValues:
    if not isinstance(policy, CallerPolicyValues):
        _fail("caller policy projection requires CallerPolicyValues")
    if not isinstance(policy.required_callers, tuple) or type(policy.values) is not _MAPPING_PROXY_TYPE:
        _fail("caller policy projection requires immutable decoded content")
    decoded = decode_caller_policy_values(_document(policy))
    if decoded != policy:
        _fail("caller policy differs from its canonical content or digest")
    return policy


def caller_policy_values_document(policy: CallerPolicyValues) -> dict[str, object]:
    """Project validated policy values as fresh canonical JSON-compatible content.

    Args:
        policy: A previously decoded immutable caller policy.

    Returns:
        Fresh closed policy values without transport asset hashes.

    Raises:
        ValueError: If typed content or its canonical digest was forged.
    """
    return _document(_require_policy(policy))
