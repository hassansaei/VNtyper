"""Closed, outcome-independent declarations for finite total-length studies."""

from __future__ import annotations

import logging
import re
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_FIELDS = {
    "schema_version",
    "objective",
    "seed",
    "fold_count",
    "grouping_rule",
    "representative_preprocessing_priority",
    "maximum_free_parameters",
    "maximum_candidate_count",
    "candidate_grid",
    "baseline",
    "tie_margin_repeat_units",
    "required_strata",
    "declared_exclusions",
    "uncertainty",
    "multiplicity",
    "acceptance",
    "qc",
    "feature_bounds_rule",
    "truth_distribution_sha256",
    "eligible_roster_sha256",
}
_ACCEPTANCE = {
    "minimum_independent_count",
    "maximum_mae",
    "minimum_relative_mae_improvement",
    "paired_error_difference_upper_limit",
    "tolerance_absolute",
    "tolerance_relative",
    "minimum_tolerance_lower_bound",
    "minimum_availability_lower_bound",
}
_QC = {
    "minimum_denominator_mean_depth",
    "minimum_denominator_covered_fraction",
    "minimum_denominator_supporting_fragments",
}
_UNCERTAINTY = {
    "bootstrap_iterations": 10000,
    "bootstrap_interval": "percentile",
    "confidence": 0.95,
    "binomial_bound": "one-sided-exact",
}
_MULTIPLICITY = {"mandatory": "intersection-union", "exploratory": "holm"}
_KINDS = {"affine-A": 2, "affine-F": 2, "physical-A": 0, "physical-F": 0}
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))


@dataclass(frozen=True)
class LengthHypothesis:
    """One predeclared promotable model family, without outcome-fitted coefficients."""

    candidate_id: str
    model_kind: str
    free_parameters: int


@dataclass(frozen=True)
class LengthProtocol:
    """Immutable study rules; decoding never grants scientific approval."""

    seed: int
    fold_count: int
    representative_preprocessing_priority: tuple[str, ...]
    maximum_free_parameters: int
    maximum_candidate_count: int
    candidates: tuple[LengthHypothesis, ...]
    required_strata: tuple[str, ...]
    declared_exclusions: tuple[str, ...]
    acceptance: Mapping[str, int | float]
    qc: Mapping[str, int | float]
    qc_sha256: str
    truth_distribution_sha256: str
    eligible_roster_sha256: str
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"length protocol {name} must contain exactly its declared fields")
    return value


def _integer(value: object, name: str, *, minimum: int = 0, maximum: int = 2**53 - 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        _fail(f"length protocol {name} must be an integer in [{minimum}, {maximum}]")
    return value


def _number(value: object, name: str, *, positive: bool = False, maximum: float = sys.float_info.max) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= maximum:
        _fail(f"length protocol {name} must be finite, nonnegative and within its allowed range")
    if positive and value == 0:
        _fail(f"length protocol {name} must be positive")
    return value


def _nonpositive_number(value: object, name: str) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not -sys.float_info.max <= value <= 0:
        _fail(f"length protocol {name} must be finite and nonpositive")
    return value


def _strings(value: object, name: str, *, ordered: bool = True, allow_empty: bool = False) -> tuple[str, ...]:
    if not isinstance(value, list) or (not value and not allow_empty):
        _fail(f"length protocol {name} must be a list with the required members")
    if any(not isinstance(item, str) or not item or item != item.strip() for item in value):
        _fail(f"length protocol {name} must contain non-empty text")
    if len(value) != len(set(value)) or (ordered and value != sorted(value)):
        _fail(f"length protocol {name} must have unique{' increasing' if ordered else ''} members")
    return tuple(value)


def _fixed(value: object, expected: object, name: str) -> None:
    # Exact type checks also keep True from matching a frozen numeric one.
    if type(value) is not type(expected) or value != expected:
        _fail(f"length protocol {name} must equal the frozen rule {expected}")


def _candidates(value: object, maximum_count: int, maximum_parameters: int) -> tuple[LengthHypothesis, ...]:
    if not isinstance(value, list) or not 1 <= len(value) <= maximum_count:
        _fail("length protocol candidate grid is empty or exceeds its declared maximum")
    result: list[LengthHypothesis] = []
    for item in value:
        row = _object(item, {"candidate_id", "model_kind"}, "candidate")
        candidate_id = _strings([row["candidate_id"]], "candidate ID")[0]
        kind = row["model_kind"]
        if not isinstance(kind, str) or kind not in _KINDS:
            _fail("length protocol candidate kind is not a supported promotable hypothesis")
        parameters = _KINDS[kind]
        if parameters > maximum_parameters:
            _fail("length protocol candidate exceeds declared free parameters")
        result.append(LengthHypothesis(candidate_id, kind, parameters))
    keys = [item.candidate_id for item in result]
    if keys != sorted(set(keys)) or len({item.model_kind for item in result}) != len(result):
        _fail("length protocol candidates must have unique increasing IDs and unique model kinds")
    return tuple(result)


def _acceptance(value: object) -> Mapping[str, int | float]:
    raw = _object(value, _ACCEPTANCE, "acceptance")
    result: dict[str, int | float] = {}
    for key in sorted(raw):
        if key == "minimum_independent_count":
            result[key] = _integer(raw[key], key, minimum=1)
        elif key == "paired_error_difference_upper_limit":
            result[key] = _nonpositive_number(raw[key], key)
        elif key in {"maximum_mae", "tolerance_absolute"}:
            result[key] = _number(raw[key], key, positive=True)
        else:
            result[key] = _number(raw[key], key, positive=True, maximum=1)
    return MappingProxyType(result)


def _qc(value: object) -> Mapping[str, int | float]:
    raw = _object(value, _QC, "QC")
    return MappingProxyType(
        {
            "minimum_denominator_mean_depth": _number(
                raw["minimum_denominator_mean_depth"], "mean depth", positive=True
            ),
            "minimum_denominator_covered_fraction": _number(
                raw["minimum_denominator_covered_fraction"],
                "covered fraction",
                positive=True,
                maximum=1,
            ),
            "minimum_denominator_supporting_fragments": _integer(
                raw["minimum_denominator_supporting_fragments"],
                "supporting fragments",
                minimum=1,
            ),
        }
    )


def decode_length_protocol(value: object) -> LengthProtocol:
    """Validate a finite declaration before training or selection outcomes are read.

    Args:
        value: Decoded strict JSON protocol document.

    Returns:
        Immutable candidates, declared limits and the complete canonical digest.
        Custom limits remain study declarations, not claims of meeting the
        separately frozen reference simulation gate.

    Raises:
        ValueError: If a rule, hypothesis, numeric limit or field is invalid.
    """
    raw = _object(value, _FIELDS, "root")
    for key, expected in {
        "schema_version": "length-protocol-v1",
        "objective": "length-total-v1",
        "grouping_rule": "connected-leakage-groups-v1",
        "baseline": "training-mean-v1",
        "tie_margin_repeat_units": 1,
        "feature_bounds_rule": "training-min-max-expand-10-percent-v1",
    }.items():
        _fixed(raw[key], expected, key)
    for name, expected in (("uncertainty", _UNCERTAINTY), ("multiplicity", _MULTIPLICITY)):
        section = _object(raw[name], set(expected), name)
        for key, rule in expected.items():
            _fixed(section[key], rule, key)
    seed = _integer(raw["seed"], "seed")
    folds = _integer(raw["fold_count"], "fold count", minimum=2)
    maximum_parameters = _integer(raw["maximum_free_parameters"], "maximum free parameters", maximum=2)
    maximum_count = _integer(raw["maximum_candidate_count"], "maximum candidate count", minimum=1, maximum=4)
    candidates = _candidates(raw["candidate_grid"], maximum_count, maximum_parameters)
    priority = _strings(raw["representative_preprocessing_priority"], "preprocessing priority", ordered=False)
    strata = _strings(raw["required_strata"], "required strata")
    exclusions = _strings(raw["declared_exclusions"], "declared exclusions", allow_empty=True)
    digests = {}
    for field in ("truth_distribution_sha256", "eligible_roster_sha256"):
        value_digest = raw[field]
        if not isinstance(value_digest, str) or not re.fullmatch(r"[0-9a-f]{64}", value_digest):
            _fail(f"length protocol {field} must be lowercase SHA256")
        digests[field] = value_digest
    qc = _qc(raw["qc"])
    return LengthProtocol(
        seed,
        folds,
        priority,
        maximum_parameters,
        maximum_count,
        candidates,
        strata,
        exclusions,
        _acceptance(raw["acceptance"]),
        qc,
        canonical_sha256(dict(qc)),
        digests["truth_distribution_sha256"],
        digests["eligible_roster_sha256"],
        canonical_sha256(raw),
    )


def length_protocol_document(protocol: LengthProtocol) -> dict[str, object]:
    """Project a typed protocol while verifying its content and stored digest.

    Args:
        protocol: Previously decoded immutable protocol.

    Returns:
        An independent JSON-compatible document.

    Raises:
        ValueError: If direct construction or replacement broke any invariant.
    """
    if not isinstance(protocol, LengthProtocol):
        _fail("length protocol must be a LengthProtocol")
    if (
        not isinstance(protocol.representative_preprocessing_priority, tuple)
        or not isinstance(protocol.candidates, tuple)
        or any(not isinstance(item, LengthHypothesis) for item in protocol.candidates)
        or not isinstance(protocol.required_strata, tuple)
        or not isinstance(protocol.declared_exclusions, tuple)
        or not isinstance(protocol.acceptance, _MAPPING_PROXY_TYPE)
        or not isinstance(protocol.qc, _MAPPING_PROXY_TYPE)
    ):
        _fail("length protocol must use decoded immutable collections")
    raw: dict[str, object] = {
        "schema_version": "length-protocol-v1",
        "objective": "length-total-v1",
        "seed": protocol.seed,
        "fold_count": protocol.fold_count,
        "grouping_rule": "connected-leakage-groups-v1",
        "representative_preprocessing_priority": list(protocol.representative_preprocessing_priority),
        "maximum_free_parameters": protocol.maximum_free_parameters,
        "maximum_candidate_count": protocol.maximum_candidate_count,
        "candidate_grid": [
            {"candidate_id": item.candidate_id, "model_kind": item.model_kind} for item in protocol.candidates
        ],
        "baseline": "training-mean-v1",
        "tie_margin_repeat_units": 1,
        "required_strata": list(protocol.required_strata),
        "declared_exclusions": list(protocol.declared_exclusions),
        "uncertainty": dict(_UNCERTAINTY),
        "multiplicity": dict(_MULTIPLICITY),
        "acceptance": dict(protocol.acceptance),
        "qc": dict(protocol.qc),
        "feature_bounds_rule": "training-min-max-expand-10-percent-v1",
        "truth_distribution_sha256": protocol.truth_distribution_sha256,
        "eligible_roster_sha256": protocol.eligible_roster_sha256,
    }
    if decode_length_protocol(raw) != protocol:
        _fail("length protocol content or digest does not match its decoded contract")
    return raw
