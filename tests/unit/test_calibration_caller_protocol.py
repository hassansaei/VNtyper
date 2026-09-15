"""Finite caller protocols bind full policies and two distinct complexity counts."""

from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
from importlib import import_module
from types import MappingProxyType

import pytest

from tests.unit.test_calibration_caller_policy import policy_document, policy_values
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, decode_caller_policy_values
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def changed_policy(**changes: object) -> CallerPolicyValues:
    raw = policy_document()
    policy_values(raw).update(changes)
    return decode_caller_policy_values(raw)


def candidate_row(policy: CallerPolicyValues, free_parameters: int) -> dict[str, object]:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    return {
        "candidate_id": policy.sha256,
        "policy": policies.caller_policy_values_document(policy),
        "free_parameters": free_parameters,
    }


def protocol_document(
    baseline: CallerPolicyValues,
    candidates: list[tuple[CallerPolicyValues, int]] | None = None,
) -> dict[str, object]:
    if candidates is None:
        candidates = [
            (
                changed_policy(**{"/components/kestrel/confidence_assignment/reporting_floor": 0.004}),
                1,
            ),
            (
                changed_policy(
                    **{
                        "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": 25,
                        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": 26,
                    }
                ),
                1,
            ),
        ]
    rows = sorted(
        (candidate_row(policy, count) for policy, count in candidates), key=lambda row: str(row["candidate_id"])
    )
    changed_coordinates: set[str] = set()
    for policy, _count in candidates:
        for pointer, value in policy.values.items():
            if value != baseline.values[pointer]:
                changed_coordinates.add(
                    "kestrel-alt-low-boundary"
                    if pointer.endswith(("alt_depth_thresholds/low", "alt_depth_thresholds/mid_low"))
                    else pointer
                )
    return {
        "schema_version": "calibration-caller-protocol-v1",
        "objective": "caller-safety-v1",
        "baseline_policy_sha256": baseline.sha256,
        "seed": 731,
        "fold_count": 5,
        "grouping_rule": "connected-leakage-groups-v1",
        "maximum_candidate_count": 8,
        "maximum_free_parameters": 4,
        "search_dimension_count": len(changed_coordinates),
        "candidate_grid": rows,
        "required_strata": ["fallback", "nominal"],
        "declared_exclusions": ["predeclared-ineligible-input"],
        "uncertainty": {
            "bootstrap_iterations": 10000,
            "bootstrap_interval": "percentile",
            "confidence": 0.95,
            "binomial_bound": "one-sided-exact",
            "paired_difference": "group-percentile",
        },
        "multiplicity": {"mandatory": "intersection-union", "exploratory": "holm"},
        "acceptance": {
            "minimum_positive_groups": 20,
            "minimum_negative_groups": 300,
            "maximum_fpr_upper": 0.01,
            "minimum_sensitivity_delta_lower": -0.02,
            "maximum_no_call_increase": 0.02,
            "minimum_exact_benefit": 0.01,
            "minimum_fpr_benefit": 0.005,
        },
    }


def test_protocol_roundtrip_binds_full_policies_and_gate_defaults() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    expected = deepcopy(raw)

    protocol = protocols.decode_caller_protocol(raw, baseline_policy=baseline)

    assert protocol.baseline_policy_sha256 == baseline.sha256
    assert protocol.required_callers == ("advntr", "kestrel")
    assert protocol.search_dimension_count == 2
    assert tuple(candidate.free_parameters for candidate in protocol.candidates) == (1, 1)
    assert protocol.gate_rules == protocols.CallerGateRules(
        required_strata=("fallback", "nominal"),
        minimum_positive_groups=20,
        minimum_negative_groups=300,
        seed=731,
    )
    assert protocol.gate_rules.maximum_fpr_upper == Fraction(1, 100)
    assert protocol.gate_rules.minimum_sensitivity_delta_lower == Fraction(-1, 50)
    assert protocol.gate_rules.maximum_no_call_increase == Fraction(1, 50)
    assert protocol.gate_rules.minimum_exact_benefit == Fraction(1, 100)
    assert protocol.gate_rules.minimum_fpr_benefit == Fraction(1, 200)
    assert protocol.sha256 == canonical_sha256(expected)
    raw["candidate_grid"] = []
    assert protocols.caller_protocol_document(protocol) == expected
    with pytest.raises(FrozenInstanceError):
        protocol.seed = 1


def test_changed_coordinate_count_treats_coupled_alt_boundary_as_one() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    alt_pair = changed_policy(
        **{
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": 25,
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": 26,
        }
    )
    alt_pair_and_high = changed_policy(
        **{
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": 25,
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": 26,
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": 120,
        }
    )

    assert protocols.caller_policy_free_parameter_count(baseline, baseline) == 0
    assert protocols.caller_policy_free_parameter_count(alt_pair, baseline) == 1
    assert protocols.caller_policy_free_parameter_count(alt_pair_and_high, baseline) == 2


def test_per_candidate_complexity_is_not_the_grid_search_dimension() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    raw["maximum_free_parameters"] = 1

    with pytest.raises(ValueError, match="search dimension"):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


def test_candidate_complexity_must_match_and_fit_the_declared_limit() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    candidate = changed_policy(
        **{
            "/components/kestrel/confidence_assignment/reporting_floor": 0.004,
            "/components/kestrel/confidence_assignment/var_active_region_threshold": 250,
        }
    )
    wrong = protocol_document(baseline, [(candidate, 1)])
    too_small = protocol_document(baseline, [(candidate, 2)])
    too_small["maximum_free_parameters"] = 1

    with pytest.raises(ValueError, match="free_parameters"):
        protocols.decode_caller_protocol(wrong, baseline_policy=baseline)
    with pytest.raises(ValueError, match="maximum_free_parameters"):
        protocols.decode_caller_protocol(too_small, baseline_policy=baseline)


def test_search_dimension_count_is_recomputed_instead_of_trusted() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    raw["search_dimension_count"] = 1

    with pytest.raises(ValueError, match="search_dimension_count"):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


def test_baseline_policy_is_a_separate_hash_bound_asset() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    raw["baseline_policy_sha256"] = "f" * 64

    with pytest.raises(ValueError, match="baseline_policy_sha256"):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)
    with pytest.raises(ValueError, match="CallerPolicyValues"):
        protocols.decode_caller_protocol(protocol_document(baseline), baseline_policy=policy_document())


def test_every_candidate_uses_the_baseline_caller_composition() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    kestrel_only = decode_caller_policy_values(policy_document(include_advntr=False))
    raw = protocol_document(baseline, [(kestrel_only, 0)])

    with pytest.raises(ValueError, match="required_callers"):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


@pytest.mark.parametrize("change", ["identifier", "order", "duplicate", "overflow", "empty"])
def test_candidate_grid_is_explicit_unique_sorted_and_capped(change: str) -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    rows = raw["candidate_grid"]
    assert isinstance(rows, list)
    if change == "identifier":
        rows[0]["candidate_id"] = "f" * 64
    elif change == "order":
        rows.reverse()
    elif change == "duplicate":
        rows[1] = deepcopy(rows[0])
    elif change == "overflow":
        raw["maximum_candidate_count"] = 1
    else:
        raw["candidate_grid"] = []

    with pytest.raises(ValueError, match="candidate"):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


@pytest.mark.parametrize("section", ["root", "candidate", "uncertainty", "multiplicity", "acceptance"])
def test_protocol_is_closed_at_every_object(section: str) -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    target: object
    if section == "root":
        target = raw
    elif section == "candidate":
        rows = raw["candidate_grid"]
        assert isinstance(rows, list)
        target = rows[0]
    else:
        target = raw[section]
    assert isinstance(target, dict)
    target["extra"] = True

    with pytest.raises(ValueError, match="fields"):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", "calibration-caller-protocol-v2"),
        ("objective", "lexicographic-safety-v1"),
        ("grouping_rule", "per-library-v1"),
    ],
)
def test_protocol_rejects_unversioned_rule_changes(field: str, value: object) -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    raw[field] = value

    with pytest.raises(ValueError, match=field):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


@pytest.mark.parametrize(
    "field,value",
    [
        ("seed", True),
        ("seed", -1),
        ("fold_count", 1),
        ("maximum_candidate_count", 0),
        ("maximum_free_parameters", True),
        ("search_dimension_count", -1),
    ],
)
def test_protocol_integer_limits_exclude_bool_and_invalid_bounds(field: str, value: object) -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    raw[field] = value

    with pytest.raises(ValueError, match=field):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


@pytest.mark.parametrize(
    "section,field,value",
    [
        ("uncertainty", "bootstrap_iterations", 9999),
        ("uncertainty", "confidence", True),
        ("uncertainty", "paired_difference", "library-percentile"),
        ("multiplicity", "mandatory", "holm"),
        ("multiplicity", "exploratory", "none"),
    ],
)
def test_uncertainty_and_multiplicity_rules_are_frozen(section: str, field: str, value: object) -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    values = raw[section]
    assert isinstance(values, dict)
    values[field] = value

    with pytest.raises(ValueError, match=field):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


@pytest.mark.parametrize(
    "field,value",
    [
        ("minimum_positive_groups", True),
        ("minimum_negative_groups", 0),
        ("maximum_fpr_upper", float("nan")),
        ("maximum_fpr_upper", 1),
        ("minimum_sensitivity_delta_lower", -1.01),
        ("maximum_no_call_increase", 1.01),
        ("minimum_exact_benefit", 0.0),
        ("minimum_fpr_benefit", -0.01),
    ],
)
def test_acceptance_limits_have_strict_types_and_domains(field: str, value: object) -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    acceptance = raw["acceptance"]
    assert isinstance(acceptance, dict)
    acceptance[field] = value

    with pytest.raises(ValueError, match=field):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


@pytest.mark.parametrize(
    "field,value",
    [
        ("required_strata", []),
        ("required_strata", ["nominal", "fallback"]),
        ("required_strata", ["nominal", "nominal"]),
        ("declared_exclusions", ["b", "a"]),
        ("declared_exclusions", [True]),
    ],
)
def test_strata_and_exclusions_are_explicit_sorted_unique_text(field: str, value: object) -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    raw[field] = value

    with pytest.raises(ValueError, match=field):
        protocols.decode_caller_protocol(raw, baseline_policy=baseline)


def test_public_projection_revalidates_nested_typed_content_and_digest() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    protocol = protocols.decode_caller_protocol(protocol_document(baseline), baseline_policy=baseline)

    with pytest.raises(ValueError, match="CallerProtocol"):
        protocols.caller_protocol_document(protocol_document(baseline))
    with pytest.raises(ValueError, match="immutable"):
        protocols.caller_protocol_document(replace(protocol, candidates=list(protocol.candidates)))
    with pytest.raises(ValueError, match="canonical"):
        protocols.caller_protocol_document(replace(protocol, sha256="f" * 64))
    with pytest.raises(ValueError, match="search_dimension_count"):
        protocols.caller_protocol_document(replace(protocol, search_dimension_count=1))
    forged = replace(protocol.candidates[0], free_parameters=4)
    with pytest.raises(ValueError, match="free_parameters"):
        protocols.caller_protocol_document(replace(protocol, candidates=(forged, *protocol.candidates[1:])))
    mutable_policy = replace(protocol.candidates[0].policy, values=dict(protocol.candidates[0].policy.values))
    mutable_candidate = replace(protocol.candidates[0], policy=mutable_policy)
    with pytest.raises(ValueError, match="immutable"):
        protocols.caller_protocol_document(replace(protocol, candidates=(mutable_candidate, *protocol.candidates[1:])))


def test_protocol_collections_do_not_accept_mutation() -> None:
    protocols = import_module("vntyper.scripts.calibration_caller_protocol")
    baseline = decode_caller_policy_values(policy_document())
    protocol = protocols.decode_caller_protocol(protocol_document(baseline), baseline_policy=baseline)

    assert isinstance(protocol.gate_rules, protocols.CallerGateRules)
    assert isinstance(protocol.baseline_policy.values, MappingProxyType)
    with pytest.raises(FrozenInstanceError):
        protocol.candidates[0].free_parameters = 10
