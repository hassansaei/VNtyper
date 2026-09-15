"""Target study declarations freeze search, role, applicability and exposure bindings."""

from copy import deepcopy
from dataclasses import replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_caller_policy import policy_document
from tests.unit.test_calibration_caller_protocol import protocol_document as caller_protocol
from tests.unit.test_calibration_candidate import candidate_document
from tests.unit.test_calibration_length_protocol import protocol_document as length_protocol
from tests.unit.test_calibration_manifest import _manifest, _member
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def study_document(target="length"):
    candidate = candidate_document(target)
    baseline = {
        "schema_version": "calibration-length-baseline-plan-v1",
        "baseline_kind": "training-mean-v1",
        "annotation_sha256": "3" * 64,
        "counting_policy_sha256": "2" * 64,
        "maximum_condition_number": 1000000,
        "producer": candidate["producer"],
    }
    protocol = length_protocol()
    if target == "callers":
        policy = policy_document()
        baseline = {
            "schema_version": "calibration-caller-baseline-plan-v1",
            "policy": policy,
            "assets_sha256": "4" * 64,
            "producer": candidate["producer"],
        }
        protocol = caller_protocol(decode_caller_policy_values(policy))
    members = []
    for index, role in enumerate(("training", "policy-selection", "validation", "locked-heldout")):
        member = _member(f"artifact-{index}", role, assay_class="capture")
        if role == "locked-heldout":
            member["provenance"] = "external-custodian"
        members.append(member)
    return {
        "schema_version": "calibration-study-v2",
        "target": target,
        "protocol": protocol,
        "partitions": _manifest(*members),
        "baseline": baseline,
        "applicability": candidate["applicability"],
        "exposure_ledger_id": "5" * 64,
    }


@pytest.mark.parametrize("target", ["length", "callers"])
def test_target_study_roundtrip_freezes_all_bindings(target):
    module = import_module("vntyper.scripts.calibration_target_contract")
    document = study_document(target)
    expected = deepcopy(document)
    study = module.decode_target_study(document)
    assert study.target == target
    assert study.sha256 == canonical_sha256(document)
    assert study.baseline.sha256 == canonical_sha256(document["baseline"])
    assert study.exposure_ledger_id == "5" * 64
    document["baseline"]["producer"]["name"] = "changed"
    assert module.target_study_document(study) == expected
    for forged in (replace(study, sha256="0" * 64), replace(study, exposure_ledger_id="0" * 64)):
        with pytest.raises(ValueError):
            module.target_study_document(forged)


@pytest.mark.parametrize("target", ["length", "callers"])
@pytest.mark.parametrize(
    "change", ["schema", "retarget", "ledger", "extra", "missing-role", "locked-provenance", "assay"]
)
def test_study_refuses_wrong_target_roles_and_bindings(target, change):
    document = study_document(target)
    if change == "schema":
        document["schema_version"] = "calibration-study-v1"
    elif change == "retarget":
        document["target"] = "callers" if target == "length" else "length"
    elif change == "ledger":
        document["exposure_ledger_id"] = "operator-name"
    elif change == "extra":
        document["unknown"] = None
    elif change == "missing-role":
        document["partitions"]["members"].pop()
    elif change == "locked-provenance":
        document["partitions"]["members"][-1]["provenance"] = "development"
    else:
        document["partitions"]["members"][0]["assay_class"] = "undeclared"
    with pytest.raises(ValueError):
        import_module("vntyper.scripts.calibration_target_contract").decode_target_study(document)


@pytest.mark.parametrize(
    "field,value",
    [
        ("baseline_kind", "global-mean"),
        ("maximum_condition_number", True),
        ("maximum_condition_number", float("inf")),
        ("maximum_condition_number", 1),
        ("counting_policy_sha256", "0" * 64),
        ("annotation_sha256", "bad"),
    ],
)
def test_length_baseline_is_a_bound_training_only_recipe(field, value):
    document = study_document()
    document["baseline"][field] = value
    with pytest.raises(ValueError):
        import_module("vntyper.scripts.calibration_target_contract").decode_target_study(document)


def test_caller_composition_must_equal_study_applicability():
    document = study_document("callers")
    document["applicability"]["required_callers"] = ["kestrel"]
    with pytest.raises(ValueError, match="caller"):
        import_module("vntyper.scripts.calibration_target_contract").decode_target_study(document)


def test_study_requires_typed_projection_and_known_target():
    module = import_module("vntyper.scripts.calibration_target_contract")
    for value in (None, {}, {**study_document(), "target": []}):
        with pytest.raises(ValueError):
            module.decode_target_study(value)
    with pytest.raises(ValueError):
        module.target_study_document(study_document())


@pytest.mark.parametrize("target", ["length", "callers"])
def test_projection_rejects_forged_nested_baseline_and_partition(target):
    module = import_module("vntyper.scripts.calibration_target_contract")
    study = module.decode_target_study(study_document(target))
    for forged in (
        replace(study, baseline=replace(study.baseline, sha256="0" * 64)),
        replace(study, baseline=None),
        replace(study, partitions=replace(study.partitions, sha256="0" * 64)),
    ):
        with pytest.raises(ValueError):
            module.target_study_document(forged)


def test_caller_baseline_schema_cannot_be_silently_reinterpreted():
    document = study_document("callers")
    document["baseline"]["schema_version"] = "calibration-length-baseline-plan-v1"
    with pytest.raises(ValueError, match="baseline schema"):
        import_module("vntyper.scripts.calibration_target_contract").decode_target_study(document)
