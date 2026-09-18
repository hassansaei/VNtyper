"""Finite simulation declarations bind resource limits and leakage families."""

from copy import deepcopy
from dataclasses import replace
from importlib import import_module

import pytest

pytestmark = pytest.mark.unit


def simulation_case(key="case-a", **changes):
    haplotype = {
        "repeat_units": ["ACGT", "TGCA"],
        "left_flank": "AACCGG",
        "right_flank": "TTGGCC",
        "edits": [],
    }
    result = {
        "case_id": key,
        "group_key": "group-a",
        "backbone_family": "backbone-a",
        "pair_family": "pair-a",
        "seed_family": "seed-a",
        "role": "development-assessment",
        "primary": True,
        "strata": ["nominal"],
        "scenario": "nominal",
        "caller_positive": False,
        "truth_variant_ids": [],
        "haplotypes": [deepcopy(haplotype), deepcopy(haplotype)],
        "reads": {
            "pair_count": 3,
            "read_length": 8,
            "fragment_length": 12,
            "seed": 42,
            "substitution_rate": 0.01,
            "quality_score": 30,
            "allele_copy_weights": [1.0, 1.0],
        },
    }
    result.update(changes)
    return result


def simulation_protocol(*cases):
    return {
        "schema_version": "calibration-simulation-protocol-v1",
        "purpose": "development-smoke",
        "repeat_unit_bp": 4,
        "maximum_haplotype_bp": 100,
        "maximum_case_count": 10,
        "maximum_generated_bases": 1000,
        "acceptance_protocol_sha256": {"callers": "a" * 64, "length": "b" * 64},
        "cases": list(cases or (simulation_case(),)),
    }


def decode(value):
    return import_module("calibration_sim.protocol").decode_simulation_protocol(value)


def test_protocol_derives_truth_and_exact_resource_counts_without_generating_reads():
    result = decode(simulation_protocol())
    assert result.generated_bases == 48
    assert result.pair_count == 3
    assert result.independent_group_count == 1
    assert result.cases[0].truth.total_repeat_count == 4
    assert result.cases[0].haplotypes[0].sequence == "AACCGGACGTTGCATTGGCC"
    module = import_module("calibration_sim.protocol")
    assert module.simulation_protocol_document(result) == simulation_protocol()


@pytest.mark.parametrize(
    "field,value",
    [
        ("maximum_case_count", 0),
        ("maximum_generated_bases", 47),
        ("maximum_haplotype_bp", 19),
        ("repeat_unit_bp", True),
        ("purpose", "validated"),
        ("cases", []),
        ("extra", "ignored"),
    ],
)
def test_invalid_protocol_limits_and_unknown_fields_are_rejected(field, value):
    raw = simulation_protocol()
    raw[field] = value
    with pytest.raises(ValueError):
        decode(raw)


@pytest.mark.parametrize(
    "changes",
    [
        {"case_id": "../escape"},
        {"case_id": "manifest.json"},
        {"case_id": "protocol.json"},
        {"role": "heldout"},
        {"primary": 1},
        {"strata": ["z", "a"]},
        {"scenario": "unexplained"},
        {"caller_positive": False, "truth_variant_ids": ["insertion"]},
        {"caller_positive": True, "truth_variant_ids": []},
        {"caller_positive": 1},
        {"unknown": None},
    ],
)
def test_case_contract_is_closed_and_truth_is_explicit(changes):
    with pytest.raises(ValueError):
        decode(simulation_protocol(simulation_case(**changes)))


@pytest.mark.parametrize(
    "reads",
    [
        {"pair_count": -1},
        {"read_length": 13},
        {"fragment_length": 30},
        {"substitution_rate": float("nan")},
        {"quality_score": 94},
        {"allele_copy_weights": [0, 0]},
        {"seed": True},
        {"extra": 0},
    ],
)
def test_read_arguments_are_validated_before_generation(reads):
    case = simulation_case()
    case["reads"].update(reads)
    with pytest.raises(ValueError):
        decode(simulation_protocol(case))


@pytest.mark.parametrize("family", ["group_key", "backbone_family", "pair_family", "seed_family"])
def test_shared_families_cannot_cross_roles_or_count_as_independent_groups(family):
    first = simulation_case()
    second = simulation_case(
        "case-b",
        group_key="group-b",
        backbone_family="backbone-b",
        pair_family="pair-b",
        seed_family="seed-b",
        role="validation",
    )
    second["reads"]["seed"] = 99
    second["haplotypes"][0]["repeat_units"] = ["AAAA", "CCCC"]
    second["haplotypes"][1]["repeat_units"] = ["CCCC", "AAAA"]
    second[family] = first[family]
    with pytest.raises(ValueError, match="famil|group|role"):
        decode(simulation_protocol(first, second))


def test_technical_arms_do_not_inflate_independent_groups_or_primary_counts():
    first = simulation_case()
    second = simulation_case("case-b", primary=False, scenario="stress")
    second["reads"]["substitution_rate"] = 0.1
    result = decode(simulation_protocol(first, second))
    assert result.independent_group_count == 1
    assert result.pair_count == 6
    assert result.generated_bases == 96
    assert len([case for case in result.cases if case.primary]) == 1


@pytest.mark.parametrize("mode", ["two-primary", "no-primary", "duplicate-id", "reversed"])
def test_every_group_has_one_frozen_primary_and_case_order_is_canonical(mode):
    first = simulation_case()
    second = simulation_case("case-b", primary=False)
    if mode == "two-primary":
        second["primary"] = True
    elif mode == "no-primary":
        first["primary"] = False
    elif mode == "duplicate-id":
        second["case_id"] = first["case_id"]
    cases = (second, first) if mode == "reversed" else (first, second)
    with pytest.raises(ValueError):
        decode(simulation_protocol(*cases))


def test_smoke_cannot_be_relabelled_as_confirmation_by_changing_roles():
    with pytest.raises(ValueError, match="smoke"):
        decode(simulation_protocol(simulation_case(role="validation")))


def test_protocol_identity_changes_with_parameters_and_projection_rejects_forgery():
    module = import_module("calibration_sim.protocol")
    raw = simulation_protocol()
    first = decode(raw)
    raw["cases"][0]["reads"]["seed"] = 43
    second = decode(raw)
    assert first.sha256 != second.sha256
    with pytest.raises(ValueError):
        module.simulation_protocol_document(replace(first, generated_bases=999))
    with pytest.raises(ValueError):
        module.simulation_protocol_document(replace(first, sha256=second.sha256))


def test_study_design_declares_disjoint_roles_without_attesting_to_power_or_validity():
    first = simulation_case(role="training")
    second = simulation_case(
        "case-b",
        role="validation",
        group_key="group-b",
        backbone_family="backbone-b",
        seed_family="seed-b",
        pair_family="pair-b",
    )
    second["reads"]["seed"] = 99
    second["haplotypes"][0]["repeat_units"] = ["AAAA", "CCCC"]
    second["haplotypes"][1]["repeat_units"] = ["CCCC", "AAAA"]
    raw = simulation_protocol(first, second)
    raw["purpose"] = "study-design"
    result = decode(raw)
    assert result.independent_group_count == 2
    assert tuple(case.role for case in result.cases) == ("training", "validation")


@pytest.mark.parametrize("collision", ["seed", "backbone"])
def test_relabelled_identical_seed_or_unedited_backbone_does_not_create_independence(collision):
    first = simulation_case()
    second = simulation_case(
        "case-b", group_key="group-b", backbone_family="backbone-b", seed_family="seed-b", pair_family="pair-b"
    )
    if collision == "seed":
        second["haplotypes"][0]["repeat_units"] = ["AAAA", "CCCC"]
        second["haplotypes"][1]["repeat_units"] = ["CCCC", "AAAA"]
    else:
        second["reads"]["seed"] = 99
        second["haplotypes"][0]["edits"] = [
            {"repeat_index": 0, "offset": 1, "deleted_bases": "", "inserted_bases": "C"}
        ]
    with pytest.raises(ValueError, match="families"):
        decode(simulation_protocol(first, second))


def test_one_shared_ancestral_haplotype_cannot_be_relabelled_as_independent():
    first = simulation_case()
    second = simulation_case(
        "case-b", group_key="group-b", backbone_family="backbone-b", seed_family="seed-b", pair_family="pair-b"
    )
    second["reads"]["seed"] = 99
    second["haplotypes"][1]["repeat_units"] = ["AAAA", "CCCC"]
    raw = simulation_protocol(first, second)
    raw["purpose"] = "study-design"

    with pytest.raises(ValueError, match="families"):
        decode(raw)


def test_aggregate_generation_budget_counts_all_technical_arms():
    raw = simulation_protocol(simulation_case(), simulation_case("case-b", primary=False))
    raw["maximum_generated_bases"] = 80
    with pytest.raises(ValueError, match="total read bases"):
        decode(raw)


def test_zero_coverage_arm_and_explicit_positive_identity_are_retained():
    case = simulation_case(caller_positive=True, truth_variant_ids=["insertion-1"])
    case["reads"]["pair_count"] = 0
    result = decode(simulation_protocol(case))
    assert result.generated_bases == 0
    assert result.cases[0].caller_positive is True
    assert result.cases[0].truth_variant_ids == ("insertion-1",)


@pytest.mark.parametrize(
    "field,value",
    [
        ("repeat_units", "ACGT"),
        ("left_flank", []),
        ("edits", {}),
        ("edits", [{"repeat_index": 0, "offset": 0, "deleted_bases": "X", "inserted_bases": "C"}]),
    ],
)
def test_haplotype_construction_rejects_malformed_units_and_edits(field, value):
    case = simulation_case()
    case["haplotypes"][0][field] = value
    with pytest.raises(ValueError):
        decode(simulation_protocol(case))


@pytest.mark.parametrize("value", [{"callers": "a" * 64}, {"callers": "z" * 64, "length": "b" * 64}])
def test_both_exact_acceptance_protocol_identities_are_required(value):
    raw = simulation_protocol()
    raw["acceptance_protocol_sha256"] = value
    with pytest.raises(ValueError):
        decode(raw)


def test_projection_rejects_noncanonical_or_invalid_json_bytes():
    module = import_module("calibration_sim.protocol")
    protocol = decode(simulation_protocol())
    for forged in (None, replace(protocol, canonical_json=b"{"), replace(protocol, canonical_json="{}")):
        with pytest.raises(ValueError):
            module.simulation_protocol_document(forged)
