"""Closed caller policy values preserve exact tunable-pointer semantics."""

from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from importlib import import_module
from types import MappingProxyType

import pytest

from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


KESTREL_VALUES: dict[str, object] = {
    "/components/kestrel/alt_filtering/gg_depth_score_threshold": 0.00469,
    "/components/kestrel/confidence_assignment/reporting_floor": 0.003,
    "/components/kestrel/confidence_assignment/var_active_region_threshold": 200,
    "/components/kestrel/confidence_assignment/depth_score_thresholds/low": 0.00469,
    "/components/kestrel/confidence_assignment/depth_score_thresholds/high": 0.00515,
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": 20,
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": 21,
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": 100,
}
ADVNTR_VALUES: dict[str, object] = {
    "/components/advntr/calibrated_calling/mode": "exact",
    "/components/advntr/calibrated_calling/cutoff": 0.001,
    "/components/advntr/calibrated_calling/minimum_read_support": 3,
    "/components/advntr/calibrated_calling/rare_unit_fraction": None,
    "/components/advntr/calibrated_calling/adapter_filter": False,
    "/components/advntr/calibrated_calling/minimum_read_match_ratio": 0.9,
    "/components/advntr/calibrated_calling/prune_reverse": False,
}


def policy_document(*, include_advntr: bool = True) -> dict[str, object]:
    values: dict[str, object] = dict(KESTREL_VALUES)
    callers = ["kestrel"]
    if include_advntr:
        values.update(ADVNTR_VALUES)
        callers = ["advntr", "kestrel"]
    return {
        "schema_version": "calibration-caller-policy-values-v1",
        "required_callers": callers,
        "values": values,
    }


def policy_values(document: dict[str, object]) -> dict[str, object]:
    values = document["values"]
    assert isinstance(values, dict)
    return values


@pytest.mark.parametrize("include_advntr", [False, True])
def test_policy_roundtrip_is_hash_bound_and_immutable(include_advntr: bool) -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document(include_advntr=include_advntr)
    expected = deepcopy(raw)

    policy = policies.decode_caller_policy_values(raw)

    assert policy.required_callers == (("advntr", "kestrel") if include_advntr else ("kestrel",))
    assert tuple(policy.values) == tuple(sorted(policy_values(expected)))
    assert policy.sha256 == canonical_sha256(expected)
    policy_values(raw)[next(iter(KESTREL_VALUES))] = 0.9
    assert policies.caller_policy_values_document(policy) == expected
    with pytest.raises(TypeError):
        policy.values[next(iter(KESTREL_VALUES))] = 0.9
    with pytest.raises(FrozenInstanceError):
        policy.sha256 = "f" * 64


def test_pointer_constants_are_the_exact_binding_inventory() -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")

    assert tuple(sorted(KESTREL_VALUES)) == policies.KESTREL_CALLER_POLICY_POINTERS
    assert tuple(sorted(ADVNTR_VALUES)) == policies.ADVNTR_CALLER_POLICY_POINTERS


@pytest.mark.parametrize("where", ["root", "values"])
@pytest.mark.parametrize("change", ["missing", "extra"])
def test_policy_rejects_missing_and_unknown_fields(where: str, change: str) -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    target = raw if where == "root" else policy_values(raw)
    assert isinstance(target, dict)
    if change == "missing":
        del target[next(iter(target))]
    else:
        target["extra"] = True

    with pytest.raises(ValueError, match="fields|pointers"):
        policies.decode_caller_policy_values(raw)


@pytest.mark.parametrize(
    "callers",
    [[], ["advntr"], ["other", "kestrel"], ["kestrel", "advntr"], ["kestrel", "kestrel"]],
)
def test_required_callers_are_sorted_unique_and_always_include_kestrel(callers: list[str]) -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    raw["required_callers"] = callers

    with pytest.raises(ValueError, match="required_callers"):
        policies.decode_caller_policy_values(raw)


def test_advntr_pointers_are_present_exactly_when_advntr_is_required() -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    missing = policy_document()
    del policy_values(missing)[next(iter(ADVNTR_VALUES))]
    unexpected = policy_document(include_advntr=False)
    policy_values(unexpected).update(ADVNTR_VALUES)

    with pytest.raises(ValueError, match="pointers"):
        policies.decode_caller_policy_values(missing)
    with pytest.raises(ValueError, match="pointers"):
        policies.decode_caller_policy_values(unexpected)


@pytest.mark.parametrize(
    "pointer,value",
    [
        ("/components/kestrel/alt_filtering/gg_depth_score_threshold", -0.01),
        ("/components/kestrel/confidence_assignment/reporting_floor", 1.01),
        ("/components/kestrel/confidence_assignment/depth_score_thresholds/low", float("nan")),
        ("/components/kestrel/confidence_assignment/depth_score_thresholds/high", float("inf")),
        ("/components/kestrel/confidence_assignment/depth_score_thresholds/high", True),
        ("/components/kestrel/confidence_assignment/depth_score_thresholds/high", 10**1_000),
    ],
)
def test_kestrel_fractional_thresholds_are_finite_unit_numbers(pointer: str, value: object) -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    policy_values(raw)[pointer] = value

    with pytest.raises(ValueError, match="finite|between"):
        policies.decode_caller_policy_values(raw)


def test_kestrel_depth_score_thresholds_are_ordered() -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    policy_values(raw)["/components/kestrel/confidence_assignment/depth_score_thresholds/low"] = 0.8
    policy_values(raw)["/components/kestrel/confidence_assignment/depth_score_thresholds/high"] = 0.7

    with pytest.raises(ValueError, match="low.*high"):
        policies.decode_caller_policy_values(raw)


@pytest.mark.parametrize("value", [True, -1, 1.5, "20"])
@pytest.mark.parametrize(
    "pointer",
    [
        "/components/kestrel/confidence_assignment/var_active_region_threshold",
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/low",
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low",
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high",
    ],
)
def test_kestrel_depth_thresholds_are_nonnegative_integers(pointer: str, value: object) -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    policy_values(raw)[pointer] = value

    with pytest.raises(ValueError, match="integer"):
        policies.decode_caller_policy_values(raw)


@pytest.mark.parametrize("low,mid_low,mid_high", [(20, 20, 100), (20, 22, 100), (20, 21, 21)])
def test_alternate_depths_preserve_the_v1_partition(low: int, mid_low: int, mid_high: int) -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    values = policy_values(raw)
    values["/components/kestrel/confidence_assignment/alt_depth_thresholds/low"] = low
    values["/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low"] = mid_low
    values["/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high"] = mid_high

    with pytest.raises(ValueError, match="alternate-depth partition"):
        policies.decode_caller_policy_values(raw)


@pytest.mark.parametrize(
    "pointer,value",
    [
        ("/components/advntr/calibrated_calling/mode", "future"),
        ("/components/advntr/calibrated_calling/mode", []),
        ("/components/advntr/calibrated_calling/cutoff", 0.0),
        ("/components/advntr/calibrated_calling/cutoff", 1.0),
        ("/components/advntr/calibrated_calling/cutoff", True),
        ("/components/advntr/calibrated_calling/minimum_read_support", 0),
        ("/components/advntr/calibrated_calling/minimum_read_support", True),
        ("/components/advntr/calibrated_calling/rare_unit_fraction", 0.0),
        ("/components/advntr/calibrated_calling/rare_unit_fraction", float("inf")),
        ("/components/advntr/calibrated_calling/adapter_filter", 0),
        ("/components/advntr/calibrated_calling/minimum_read_match_ratio", 0.0),
        ("/components/advntr/calibrated_calling/minimum_read_match_ratio", 1),
        ("/components/advntr/calibrated_calling/prune_reverse", "false"),
    ],
)
def test_advntr_policy_values_have_exact_r3_types_and_domains(pointer: str, value: object) -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    policy_values(raw)[pointer] = value

    with pytest.raises(ValueError, match="mode|float|integer|boolean"):
        policies.decode_caller_policy_values(raw)


def test_valid_domain_boundaries_remain_available() -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    raw = policy_document()
    values = policy_values(raw)
    values["/components/kestrel/alt_filtering/gg_depth_score_threshold"] = 0
    values["/components/kestrel/confidence_assignment/reporting_floor"] = 1
    values["/components/kestrel/confidence_assignment/depth_score_thresholds/low"] = 0
    values["/components/kestrel/confidence_assignment/depth_score_thresholds/high"] = 1
    values["/components/kestrel/confidence_assignment/var_active_region_threshold"] = 0
    values["/components/advntr/calibrated_calling/rare_unit_fraction"] = 1.0
    values["/components/advntr/calibrated_calling/minimum_read_match_ratio"] = 1.0

    assert policies.caller_policy_values_document(policies.decode_caller_policy_values(raw)) == raw


def test_public_projection_revalidates_typed_content_and_digest() -> None:
    policies = import_module("vntyper.scripts.calibration_caller_policy")
    policy = policies.decode_caller_policy_values(policy_document())

    with pytest.raises(ValueError, match="CallerPolicyValues"):
        policies.caller_policy_values_document(policy_document())
    with pytest.raises(ValueError, match="immutable"):
        policies.caller_policy_values_document(replace(policy, values=dict(policy.values)))
    with pytest.raises(ValueError, match="canonical"):
        policies.caller_policy_values_document(replace(policy, sha256="f" * 64))
    forged_values = dict(policy.values)
    forged_values["/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low"] = 22
    with pytest.raises(ValueError, match="alternate-depth partition"):
        policies.caller_policy_values_document(replace(policy, values=MappingProxyType(forged_values)))


def test_expected_refusals_silence_only_the_error_log_and_only_inside_the_block(caplog):
    import logging

    from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values, expected_refusals

    caplog.set_level(logging.INFO)
    document = policy_document(include_advntr=True)
    document["values"]["/components/advntr/calibrated_calling/cutoff"] = 2.0

    with expected_refusals(), pytest.raises(ValueError, match="strictly between zero and one"):
        decode_caller_policy_values(document)
    assert [record for record in caplog.records if record.levelno >= logging.ERROR] == []

    with pytest.raises(ValueError, match="strictly between zero and one"):
        decode_caller_policy_values(document)
    assert [record.levelno for record in caplog.records if record.levelno >= logging.ERROR] == [logging.ERROR]
