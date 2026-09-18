"""Closed outcome-independent protocols for finite caller policy studies."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from typing import NoReturn

from vntyper.scripts.calibration_caller_acceptance import CallerGateRules
from vntyper.scripts.calibration_caller_policy import (
    CallerName,
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_FIELDS = {
    "schema_version",
    "objective",
    "baseline_policy_sha256",
    "seed",
    "fold_count",
    "grouping_rule",
    "maximum_candidate_count",
    "maximum_free_parameters",
    "search_dimension_count",
    "candidate_grid",
    "required_strata",
    "declared_exclusions",
    "uncertainty",
    "multiplicity",
    "acceptance",
}
_CANDIDATE_FIELDS = {"candidate_id", "policy", "free_parameters"}
_ACCEPTANCE_FIELDS = {
    "minimum_positive_groups",
    "minimum_negative_groups",
    "maximum_fpr_upper",
    "minimum_sensitivity_delta_lower",
    "maximum_no_call_increase",
    "minimum_exact_benefit",
    "minimum_fpr_benefit",
}
_UNCERTAINTY = {
    "bootstrap_iterations": 10_000,
    "bootstrap_interval": "percentile",
    "confidence": 0.95,
    "binomial_bound": "one-sided-exact",
    "paired_difference": "group-percentile",
}
_MULTIPLICITY = {"mandatory": "intersection-union", "exploratory": "holm"}
_ALT_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/low"
_ALT_MID_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low"
_ALT_BOUNDARY_COORDINATE = "kestrel-alt-low-boundary"


@dataclass(frozen=True)
class CallerPolicyCandidate:
    """One full protocol-enumerated policy and its baseline-relative complexity."""

    candidate_id: str
    policy: CallerPolicyValues
    free_parameters: int


@dataclass(frozen=True)
class CallerProtocol:
    """Finite caller study rules frozen before any candidate outcome is opened.

    ``free_parameters`` on each candidate is its policy-complexity tie-breaker
    relative to the baseline. ``search_dimension_count`` is instead the union
    of coordinates explored anywhere in the grid. Neither value adjusts for
    multiple comparisons; the frozen multiplicity rules govern that separately.
    """

    baseline_policy: CallerPolicyValues
    baseline_policy_sha256: str
    required_callers: tuple[CallerName, ...]
    seed: int
    fold_count: int
    maximum_candidate_count: int
    maximum_free_parameters: int
    search_dimension_count: int
    candidates: tuple[CallerPolicyCandidate, ...]
    required_strata: tuple[str, ...]
    declared_exclusions: tuple[str, ...]
    gate_rules: CallerGateRules
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"caller protocol {label} fields differ from the closed contract")
    return value


def _fixed(value: object, expected: object, field: str) -> None:
    if type(value) is not type(expected) or value != expected:
        _fail(f"caller protocol {field} must equal the frozen rule {expected}")


def _integer(value: object, field: str, *, minimum: int = 0, maximum: int = 2**53 - 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        _fail(f"caller protocol {field} must be an integer in [{minimum}, {maximum}]")
    return value


def _strings(value: object, field: str, *, allow_empty: bool = False) -> tuple[str, ...]:
    if not isinstance(value, list) or (not value and not allow_empty):
        _fail(f"caller protocol {field} must be a list with the required members")
    if any(not isinstance(item, str) or not item or item.strip() != item for item in value):
        _fail(f"caller protocol {field} must contain non-empty trimmed text")
    if value != sorted(set(value)):
        _fail(f"caller protocol {field} must contain sorted unique text")
    return tuple(value)


def _rate(
    value: object,
    field: str,
    *,
    minimum: float,
    maximum: float,
    positive: bool = False,
) -> Fraction:
    # Canonical JSON writes integral rates such as 1.0 as 1; both numeric forms
    # must decode identically while booleans remain invalid protocol values.
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not minimum <= value <= maximum
        or not math.isfinite(value)
    ):
        _fail(f"caller protocol acceptance {field} must be a finite number in [{minimum}, {maximum}]")
    if positive and value == 0:
        _fail(f"caller protocol acceptance {field} must be positive")
    return Fraction(str(value))


def _require_policy(policy: CallerPolicyValues, label: str) -> CallerPolicyValues:
    try:
        caller_policy_values_document(policy)
    except ValueError as error:
        _fail(f"caller protocol {label} requires decoded immutable CallerPolicyValues: {error}")
    return policy


def _changed_coordinates(candidate: CallerPolicyValues, baseline: CallerPolicyValues) -> frozenset[str]:
    current = _require_policy(candidate, "candidate policy")
    old = _require_policy(baseline, "baseline policy")
    if current.required_callers != old.required_callers:
        _fail("caller protocol candidate required_callers must equal the frozen baseline")
    changed: set[str] = set()
    for pointer, value in current.values.items():
        if value == old.values[pointer]:
            continue
        changed.add(_ALT_BOUNDARY_COORDINATE if pointer in {_ALT_LOW, _ALT_MID_LOW} else pointer)
    return frozenset(changed)


def caller_policy_free_parameter_count(candidate: CallerPolicyValues, baseline: CallerPolicyValues) -> int:
    """Count changed policy coordinates for the selection simplicity tie-breaker.

    The linked alternate-depth ``low`` and ``mid_low`` values define one
    boundary coordinate. This candidate-relative count is not the number of
    dimensions explored by the complete grid and is not a multiplicity adjustment.

    Args:
        candidate: One complete validated candidate policy.
        baseline: The separately frozen complete baseline policy.

    Returns:
        Number of independent policy coordinates changed from baseline.

    Raises:
        ValueError: If either policy is forged or their caller sets differ.
    """
    return len(_changed_coordinates(candidate, baseline))


def _acceptance(value: object, strata: tuple[str, ...], seed: int) -> CallerGateRules:
    raw = _object(value, _ACCEPTANCE_FIELDS, "acceptance")
    return CallerGateRules(
        required_strata=strata,
        minimum_positive_groups=_integer(raw["minimum_positive_groups"], "minimum_positive_groups", minimum=1),
        minimum_negative_groups=_integer(raw["minimum_negative_groups"], "minimum_negative_groups", minimum=1),
        seed=seed,
        maximum_fpr_upper=_rate(raw["maximum_fpr_upper"], "maximum_fpr_upper", minimum=0, maximum=1),
        minimum_sensitivity_delta_lower=_rate(
            raw["minimum_sensitivity_delta_lower"],
            "minimum_sensitivity_delta_lower",
            minimum=-1,
            maximum=1,
        ),
        maximum_no_call_increase=_rate(
            raw["maximum_no_call_increase"], "maximum_no_call_increase", minimum=0, maximum=1
        ),
        minimum_exact_benefit=_rate(
            raw["minimum_exact_benefit"], "minimum_exact_benefit", minimum=0, maximum=1, positive=True
        ),
        minimum_fpr_benefit=_rate(
            raw["minimum_fpr_benefit"], "minimum_fpr_benefit", minimum=0, maximum=1, positive=True
        ),
    )


def _candidates(
    value: object,
    baseline: CallerPolicyValues,
    maximum_count: int,
    maximum_parameters: int,
) -> tuple[tuple[CallerPolicyCandidate, ...], int]:
    if not isinstance(value, list) or not 1 <= len(value) <= maximum_count:
        _fail("caller protocol candidate grid is empty or exceeds maximum_candidate_count")
    result: list[CallerPolicyCandidate] = []
    searched: set[str] = set()
    for item in value:
        row = _object(item, _CANDIDATE_FIELDS, "candidate")
        policy = decode_caller_policy_values(row["policy"])
        coordinates = _changed_coordinates(policy, baseline)
        actual = len(coordinates)
        declared = _integer(row["free_parameters"], "candidate free_parameters")
        if declared != actual:
            _fail("caller protocol candidate free_parameters differs from changed policy coordinates")
        if actual > maximum_parameters:
            _fail("caller protocol candidate exceeds maximum_free_parameters")
        candidate_id = row["candidate_id"]
        if not isinstance(candidate_id, str) or candidate_id != policy.sha256:
            _fail("caller protocol candidate_id must equal the full policy digest")
        result.append(CallerPolicyCandidate(candidate_id, policy, actual))
        searched.update(coordinates)
    identities = [candidate.candidate_id for candidate in result]
    if identities != sorted(set(identities)):
        _fail("caller protocol candidates must have sorted unique candidate IDs and policies")
    if len(searched) > maximum_parameters:
        _fail("caller protocol search dimension exceeds maximum_free_parameters")
    return tuple(result), len(searched)


def _validate_fixed_sections(raw: Mapping[str, object]) -> None:
    for field, expected in {
        "schema_version": "calibration-caller-protocol-v1",
        "objective": "caller-safety-v1",
        "grouping_rule": "connected-leakage-groups-v1",
    }.items():
        _fixed(raw[field], expected, field)
    for name, rules in (("uncertainty", _UNCERTAINTY), ("multiplicity", _MULTIPLICITY)):
        section = _object(raw[name], set(rules), name)
        for field, rule in rules.items():
            _fixed(section[field], rule, field)


def decode_caller_protocol(value: object, *, baseline_policy: CallerPolicyValues) -> CallerProtocol:
    """Decode a finite caller study declaration against its baseline asset.

    Args:
        value: Closed protocol document containing sorted full-policy candidates.
        baseline_policy: Separately loaded complete baseline policy bound by digest.

    Returns:
        Immutable protocol with verified candidate and grid-wide complexity.

    Raises:
        ValueError: If any field, policy, limit, fixed method or digest is invalid.
    """
    raw = _object(value, _FIELDS, "root")
    _validate_fixed_sections(raw)
    baseline = _require_policy(baseline_policy, "baseline_policy")
    if raw["baseline_policy_sha256"] != baseline.sha256:
        _fail("caller protocol baseline_policy_sha256 differs from the opened baseline policy")
    seed = _integer(raw["seed"], "seed")
    folds = _integer(raw["fold_count"], "fold_count", minimum=2)
    maximum_count = _integer(raw["maximum_candidate_count"], "maximum_candidate_count", minimum=1)
    coordinate_capacity = len(baseline.values) - 1
    maximum_parameters = _integer(
        raw["maximum_free_parameters"],
        "maximum_free_parameters",
        maximum=coordinate_capacity,
    )
    declared_search_dimension = _integer(raw["search_dimension_count"], "search_dimension_count")
    candidates, search_dimension = _candidates(raw["candidate_grid"], baseline, maximum_count, maximum_parameters)
    if declared_search_dimension != search_dimension:
        _fail("caller protocol search_dimension_count differs from the derived grid-wide coordinate union")
    strata = _strings(raw["required_strata"], "required_strata")
    exclusions = _strings(raw["declared_exclusions"], "declared_exclusions", allow_empty=True)
    rules = _acceptance(raw["acceptance"], strata, seed)
    return CallerProtocol(
        baseline,
        baseline.sha256,
        baseline.required_callers,
        seed,
        folds,
        maximum_count,
        maximum_parameters,
        search_dimension,
        candidates,
        strata,
        exclusions,
        rules,
        canonical_sha256(raw),
    )


def _acceptance_document(rules: CallerGateRules) -> dict[str, object]:
    return {
        "minimum_positive_groups": rules.minimum_positive_groups,
        "minimum_negative_groups": rules.minimum_negative_groups,
        "maximum_fpr_upper": float(rules.maximum_fpr_upper),
        "minimum_sensitivity_delta_lower": float(rules.minimum_sensitivity_delta_lower),
        "maximum_no_call_increase": float(rules.maximum_no_call_increase),
        "minimum_exact_benefit": float(rules.minimum_exact_benefit),
        "minimum_fpr_benefit": float(rules.minimum_fpr_benefit),
    }


def _document(protocol: CallerProtocol) -> dict[str, object]:
    return {
        "schema_version": "calibration-caller-protocol-v1",
        "objective": "caller-safety-v1",
        "baseline_policy_sha256": protocol.baseline_policy_sha256,
        "seed": protocol.seed,
        "fold_count": protocol.fold_count,
        "grouping_rule": "connected-leakage-groups-v1",
        "maximum_candidate_count": protocol.maximum_candidate_count,
        "maximum_free_parameters": protocol.maximum_free_parameters,
        "search_dimension_count": protocol.search_dimension_count,
        "candidate_grid": [
            {
                "candidate_id": candidate.candidate_id,
                "policy": caller_policy_values_document(candidate.policy),
                "free_parameters": candidate.free_parameters,
            }
            for candidate in protocol.candidates
        ],
        "required_strata": list(protocol.required_strata),
        "declared_exclusions": list(protocol.declared_exclusions),
        "uncertainty": dict(_UNCERTAINTY),
        "multiplicity": dict(_MULTIPLICITY),
        "acceptance": _acceptance_document(protocol.gate_rules),
    }


def _require_protocol(protocol: CallerProtocol) -> CallerProtocol:
    if not isinstance(protocol, CallerProtocol):
        _fail("caller protocol projection requires CallerProtocol")
    if (
        not isinstance(protocol.required_callers, tuple)
        or not isinstance(protocol.candidates, tuple)
        or any(not isinstance(candidate, CallerPolicyCandidate) for candidate in protocol.candidates)
        or not isinstance(protocol.required_strata, tuple)
        or not isinstance(protocol.declared_exclusions, tuple)
        or not isinstance(protocol.gate_rules, CallerGateRules)
    ):
        _fail("caller protocol projection requires decoded immutable collections")
    decoded = decode_caller_protocol(_document(protocol), baseline_policy=protocol.baseline_policy)
    if decoded != protocol:
        _fail("caller protocol differs from its canonical content or digest")
    return protocol


def caller_protocol_document(protocol: CallerProtocol) -> dict[str, object]:
    """Project a typed caller protocol after full content and digest revalidation.

    Args:
        protocol: Previously decoded immutable caller protocol.

    Returns:
        Fresh JSON-compatible protocol with verified derived complexity counts.

    Raises:
        ValueError: If direct construction or replacement broke any invariant.
    """
    return _document(_require_protocol(protocol))
