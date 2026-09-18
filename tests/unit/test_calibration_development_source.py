"""Closed previously examined evidence source for target assessment."""

from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_candidate import candidate_document
from tests.unit.test_calibration_target_runs import run_document
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def source_document(target="length"):
    candidate = import_module("vntyper.scripts.calibration_candidate").decode_candidate(candidate_document(target))
    runs = import_module("vntyper.scripts.calibration_target_runs").decode_target_runs(run_document(target=target))
    roster = [{"key": "artifact-001", "group_key": "group-001", "strata": ["all"]}]
    raw = {
        "schema_version": "calibration-development-source-v1",
        "target": target,
        "evidence_role": "development-assessment",
        "promotion_eligible": False,
        "candidate_sha256": candidate.sha256,
        "study_sha256": candidate.study_sha256,
        "run_manifest_sha256": runs.sha256,
        "partition_sha256": canonical_sha256({"members": roster}),
        "roster": roster,
        "identities_by_key": {
            "artifact-001": [
                {"namespace": "physical-readset", "sha256": "b" * 64},
                {"namespace": "specimen", "sha256": "c" * 64},
            ]
        },
        "evidence_domains": {"artifact-001": "external"},
        "previously_examined": {"artifact-001": True},
        "truth_asset": {"path": "/invented/truth.json", "sha256": "e" * 64, "size_bytes": 17},
    }
    return raw, candidate, runs


@pytest.mark.parametrize("target", ["length", "callers"])
def test_development_source_roundtrip_binds_candidate_runs_members_and_all_identities(target):
    module = import_module("vntyper.scripts.calibration_development_source")
    raw, candidate, runs = source_document(target)
    source = module.decode_development_source(raw, candidate=candidate, runs=runs)
    assert source.keys == ("artifact-001",)
    assert source.identities == (
        ("physical-readset", "b" * 64),
        ("specimen", "c" * 64),
    )
    assert source.promotion_eligible is False
    assert module.development_source_document(source, candidate=candidate, runs=runs) == raw
    with pytest.raises(FrozenInstanceError):
        source.partition_sha256 = "0" * 64


@pytest.mark.parametrize(
    "change,match",
    [
        (("evidence_role", "validation"), "development-assessment"),
        (("promotion_eligible", True), "promotion"),
        (("candidate_sha256", "0" * 64), "candidate"),
        (("study_sha256", "0" * 64), "study"),
        (("run_manifest_sha256", "0" * 64), "run"),
        (("target", "callers"), "target"),
    ],
)
def test_development_source_rejects_promotion_role_and_context_drift(change, match):
    module = import_module("vntyper.scripts.calibration_development_source")
    raw, candidate, runs = source_document()
    raw[change[0]] = change[1]
    with pytest.raises(ValueError, match=match):
        module.decode_development_source(raw, candidate=candidate, runs=runs)


@pytest.mark.parametrize("mode", ["missing-flag", "false-flag", "identity-key", "domain", "physical"])
def test_development_source_requires_every_member_to_be_previously_examined_and_audited(mode):
    module = import_module("vntyper.scripts.calibration_development_source")
    raw, candidate, runs = source_document()
    if mode == "missing-flag":
        raw["previously_examined"] = {}
    elif mode == "false-flag":
        raw["previously_examined"]["artifact-001"] = False
    elif mode == "identity-key":
        raw["identities_by_key"] = {}
    elif mode == "domain":
        raw["evidence_domains"]["artifact-001"] = "guessed"
    else:
        raw["identities_by_key"]["artifact-001"][0]["sha256"] = "9" * 64
    with pytest.raises(ValueError):
        module.decode_development_source(raw, candidate=candidate, runs=runs)


def test_development_source_rejects_unknown_fields_mutability_and_stale_digest():
    module = import_module("vntyper.scripts.calibration_development_source")
    raw, candidate, runs = source_document()
    source = module.decode_development_source(raw, candidate=candidate, runs=runs)
    changed = deepcopy(raw)
    changed["unknown"] = True
    with pytest.raises(ValueError, match="fields"):
        module.decode_development_source(changed, candidate=candidate, runs=runs)
    with pytest.raises(ValueError, match="immutable|identity"):
        module.development_source_document(
            replace(source, evidence_domains=dict(source.evidence_domains)), candidate=candidate, runs=runs
        )
    with pytest.raises(ValueError, match="identity"):
        module.development_source_document(replace(source, sha256="0" * 64), candidate=candidate, runs=runs)
