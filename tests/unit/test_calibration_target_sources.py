"""Trusted curation splits normalized truth into sealed target role sources."""

from copy import deepcopy
from dataclasses import replace
from importlib import import_module
from pathlib import Path
from typing import cast

import pytest

from tests.unit.test_calibration_identity import audit, intake_for
from tests.unit.test_calibration_target_contract import study_document
from tests.unit.test_calibration_target_runs import run_document
from vntyper.scripts.calibration_caller_observations import decode_caller_truth
from vntyper.scripts.calibration_identity import build_partition_members
from vntyper.scripts.calibration_intake_contract import decode_intake
from vntyper.scripts.calibration_length_controller import _truth
from vntyper.scripts.calibration_manifest import GROUP_NAMESPACES
from vntyper.scripts.calibration_role_source import decode_role_source
from vntyper.scripts.calibration_target_contract import decode_target_study, target_study_document
from vntyper.scripts.calibration_target_runs import decode_target_runs, target_runs_document
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION

pytestmark = pytest.mark.unit


def curator_fixture(target, tmp_path):
    raw = intake_for(3)
    roles = ("training", "policy-selection", "validation")
    for index, role in enumerate(roles):
        raw["assignments"][index]["role"] = role
        raw["artifacts"][index].update(assay_class="capture", expected_sha256=None)
        raw["truth"].append(
            {
                "specimen_key": f"sample-{index}",
                "genotype": "positive" if index != 1 else "negative",
                "variants": [f"variant-{index}"] if index != 1 else [],
                "length": {
                    "allele_1": 40 + index,
                    "allele_2": 70 + index,
                    "unit": "repeat-count",
                    "repeat_unit_bp": 60,
                    "boundary_definition": TARGET_BOUNDARY_DEFINITION,
                    "measurement": "exact",
                    "lower_bound": None,
                    "upper_bound": None,
                    "conversion_id": None,
                },
                "method": "synthetic-truth-v1",
                "source_digest": canonical_sha256(f"truth-{index}"),
                "source_row": index + 1,
                "status": "confirmed",
            }
        )
    declaration = decode_intake(raw)
    identity_audit = audit(raw)
    partition_members = build_partition_members(declaration, identity_audit)
    study_raw = study_document(target)
    study_raw["partitions"]["members"][:3] = [
        {
            "key": member.key,
            "role": member.role,
            "provenance": member.provenance,
            "assay_class": member.assay_class,
            "groups": {name: list(member.groups[name]) for name in GROUP_NAMESPACES},
        }
        for member in partition_members
    ]
    study = decode_target_study(study_raw)
    runs_raw = {"schema_version": "calibration-runs-v2", "target": target, "runs": []}
    for index, artifact in enumerate(declaration.artifacts):
        row = deepcopy(run_document(target=target)["runs"][0])
        row["manifest_key"] = artifact.key
        row["input_sha256"] = canonical_sha256(f"run-input-{index}")
        runs_raw["runs"].append(row)
    runs = decode_target_runs(runs_raw)
    return declaration, identity_audit, study, runs


@pytest.mark.parametrize("target", ["length", "callers"])
def test_curator_prepares_one_sealed_role_source_from_revalidated_intake(target, tmp_path):
    module = import_module("vntyper.scripts.calibration_target_sources")
    declaration, identity_audit, study, runs = curator_fixture(target, tmp_path)
    staging = tmp_path / "staging"
    staging.mkdir()
    installed = (tmp_path / "installed").resolve()
    source = module.prepare_target_role_source(
        declaration,
        identity_audit,
        study,
        runs,
        staging,
        installed_root=installed,
        role="training",
        evidence_domains={"artifact-0": "synthetic"},
        strata_by_key={"artifact-0": ("all",)},
    )
    source_raw = load_strict_json_object((staging / "roles/training/source.json").read_bytes())
    assert decode_role_source(source_raw, study=study, runs=runs, expected_role="training") == source
    truth_raw = load_strict_json_object((staging / "roles/training/truth.json").read_bytes())
    if target == "length":
        assert _truth(truth_raw, source.keys) == {"artifact-0": 110.0}
    else:
        assert decode_caller_truth(truth_raw, source.keys).by_key["artifact-0"].variants == ("variant-0",)
    curation = load_strict_json_object((staging / "curation/training.json").read_bytes())
    assert curation["previously_examined"] == {"artifact-0": False}
    assert curation["source_sha256"] == source.sha256
    assert source.truth_asset.path == installed / "roles/training/truth.json"
    second = module.prepare_target_role_source(
        declaration,
        identity_audit,
        study,
        runs,
        staging,
        installed_root=installed,
        role="policy-selection",
        evidence_domains={"artifact-1": "synthetic"},
        strata_by_key={"artifact-1": ("all",)},
    )
    assert second.role == "policy-selection"
    assert {path.name for path in (staging / "roles").iterdir()} == {"training", "policy-selection"}
    for path in staging.rglob("*.json"):
        assert path.stat().st_mode & 0o777 == 0o600


def test_curator_applies_explicit_registered_length_boundary_conversion(tmp_path):
    module = import_module("vntyper.scripts.calibration_target_sources")
    truth_module = import_module("vntyper.scripts.calibration_truth")
    declaration, identity_audit, study, runs = curator_fixture("length", tmp_path)
    encoded = import_module("vntyper.scripts.calibration_intake_contract").encode_intake(declaration)
    encoded["truth"][0]["length"].update(
        allele_1=49,
        allele_2=79,
        boundary_definition="reported-including-nine-terminal-v1",
        conversion_id="subtract-nine-terminal-per-allele-v1",
    )
    declaration = decode_intake(encoded)
    identity_audit = audit(encoded)
    registry = truth_module.decode_conversion_registry(
        {
            "schema_version": "calibration-length-conversions-v2",
            "conversions": [
                {
                    "conversion_id": "subtract-nine-terminal-per-allele-v1",
                    "source_unit": "repeat-count",
                    "source_boundary_definition": "reported-including-nine-terminal-v1",
                    "target_boundary_definition": TARGET_BOUNDARY_DEFINITION,
                    "repeat_unit_bp": 60,
                    "allele_1_adjustment": -9,
                    "allele_2_adjustment": -9,
                }
            ],
        }
    )
    staging = tmp_path / "staging"
    staging.mkdir()
    source = module.prepare_target_role_source(
        declaration,
        identity_audit,
        study,
        runs,
        staging,
        installed_root=(tmp_path / "installed").resolve(),
        role="training",
        evidence_domains={"artifact-0": "external"},
        strata_by_key={"artifact-0": ("all",)},
        conversion_registry=registry,
    )
    assert _truth(load_strict_json_object((staging / "roles/training/truth.json").read_bytes()), source.keys) == {
        "artifact-0": 110.0
    }


def test_curator_preserves_confirmed_unknown_caller_truth_as_unassessable_identity(tmp_path):
    module = import_module("vntyper.scripts.calibration_target_sources")
    declaration, _, study, runs = curator_fixture("callers", tmp_path)
    encoded = import_module("vntyper.scripts.calibration_intake_contract").encode_intake(declaration)
    encoded["truth"][0].update(genotype="unknown", variants=[])
    declaration = decode_intake(encoded)
    identity_audit = audit(encoded)
    staging = tmp_path / "staging"
    staging.mkdir()
    source = module.prepare_target_role_source(
        declaration,
        identity_audit,
        study,
        runs,
        staging,
        installed_root=(tmp_path / "installed").resolve(),
        role="training",
        evidence_domains={"artifact-0": "external"},
        strata_by_key={"artifact-0": ("all",)},
    )
    truth = decode_caller_truth(
        load_strict_json_object((staging / "roles/training/truth.json").read_bytes()), source.keys
    )
    assert truth.by_key["artifact-0"].genotype == "unknown"
    assert truth.by_key["artifact-0"].variants is None


@pytest.mark.parametrize("mode", ["locked", "stale-audit", "domain", "strata", "strata-type", "output"])
def test_curator_rejects_unauthorized_or_incomplete_preparation_without_writing(tmp_path, mode):
    module = import_module("vntyper.scripts.calibration_target_sources")
    declaration, identity_audit, study, runs = curator_fixture("length", tmp_path)
    staging = tmp_path / "staging"
    staging.mkdir()
    role = "locked-heldout" if mode == "locked" else "training"
    domains = {} if mode == "domain" else {"artifact-0": "synthetic"}
    strata = {} if mode == "strata" else {"artifact-0": ("all",)}
    if mode == "strata-type":
        strata = {"artifact-0": cast(tuple[str, ...], ["all"])}
    installed = Path("relative") if mode == "output" else (tmp_path / "installed").resolve()
    if mode == "stale-audit":
        identity_audit = replace(identity_audit, sha256="0" * 64)
    with pytest.raises(ValueError):
        module.prepare_target_role_source(
            declaration,
            identity_audit,
            study,
            runs,
            staging,
            installed_root=installed,
            role=role,
            evidence_domains=domains,
            strata_by_key=strata,
        )
    assert not tuple(staging.iterdir())


def test_curator_refuses_to_mix_existing_sources_from_another_context(tmp_path):
    module = import_module("vntyper.scripts.calibration_target_sources")
    declaration, identity_audit, study, runs = curator_fixture("length", tmp_path)
    staging = tmp_path / "staging"
    staging.mkdir()
    module.prepare_target_role_source(
        declaration,
        identity_audit,
        study,
        runs,
        staging,
        installed_root=(tmp_path / "installed").resolve(),
        role="training",
        evidence_domains={"artifact-0": "synthetic"},
        strata_by_key={"artifact-0": ("all",)},
    )
    audit_path = staging / "curation/training.json"
    raw = load_strict_json_object(audit_path.read_bytes())
    raw["study_sha256"] = "0" * 64
    audit_path.write_bytes(import_module("vntyper.scripts.canonical_json").canonical_json_bytes(raw))
    with pytest.raises(ValueError, match="current trust context"):
        module.prepare_target_role_source(
            declaration,
            identity_audit,
            study,
            runs,
            staging,
            installed_root=(tmp_path / "installed").resolve(),
            role="policy-selection",
            evidence_domains={"artifact-1": "synthetic"},
            strata_by_key={"artifact-1": ("all",)},
        )
    assert not (staging / "roles/policy-selection").exists()


@pytest.mark.parametrize(
    "mode", ["unconverted", "no-truth", "target", "unrelated-output", "missing-run", "membership", "quarantine"]
)
def test_curator_fails_closed_for_missing_truth_or_context(mode, tmp_path):
    module = import_module("vntyper.scripts.calibration_target_sources")
    declaration, identity_audit, study, runs = curator_fixture("length", tmp_path)
    if mode == "unconverted":
        encoded = import_module("vntyper.scripts.calibration_intake_contract").encode_intake(declaration)
        encoded["truth"][0]["length"].update(boundary_definition="other-boundary", conversion_id="required-conversion")
        declaration = decode_intake(encoded)
        identity_audit = audit(encoded)
    elif mode == "no-truth":
        encoded = import_module("vntyper.scripts.calibration_intake_contract").encode_intake(declaration)
        encoded["truth"] = [row for row in encoded["truth"] if row["specimen_key"] != "sample-0"]
        declaration = decode_intake(encoded)
        identity_audit = audit(encoded)
    elif mode == "target":
        runs = replace(runs, target="callers")
    elif mode == "missing-run":
        raw_runs = target_runs_document(runs)
        raw_runs["runs"] = [row for row in raw_runs["runs"] if row["manifest_key"] != "artifact-0"]
        runs = decode_target_runs(raw_runs)
    elif mode == "membership":
        raw_study = target_study_document(study)
        raw_study["partitions"]["members"][0]["groups"]["batch"] = ["changed-membership"]
        study = decode_target_study(raw_study)
    elif mode == "quarantine":
        encoded = import_module("vntyper.scripts.calibration_intake_contract").encode_intake(declaration)
        encoded["specimens"][0].update(individual_key=None, identity_status="unresolved")
        declaration = decode_intake(encoded)
        identity_audit = audit(encoded)
        raw_study = target_study_document(study)
        replacement = build_partition_members(declaration, identity_audit)[0]
        raw_study["partitions"]["members"][0] = {
            "key": replacement.key,
            "role": replacement.role,
            "provenance": replacement.provenance,
            "assay_class": replacement.assay_class,
            "groups": {name: list(replacement.groups[name]) for name in GROUP_NAMESPACES},
        }
        study = decode_target_study(raw_study)
    staging = tmp_path / "staging"
    staging.mkdir()
    if mode == "unrelated-output":
        (staging / "unrelated").write_text("forbidden")
    with pytest.raises(ValueError):
        module.prepare_target_role_source(
            declaration,
            identity_audit,
            study,
            runs,
            staging,
            installed_root=(tmp_path / "installed").resolve(),
            role="training",
            evidence_domains={"artifact-0": "external"},
            strata_by_key={"artifact-0": ("all",)},
        )
    assert not (staging / "roles").exists()
