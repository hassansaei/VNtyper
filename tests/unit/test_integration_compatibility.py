"""Tests for the append-only real-integration compatibility contract."""

import copy
import importlib
import json
import subprocess
from pathlib import Path
from typing import Any

import pytest

pytestmark = pytest.mark.unit

BOOTSTRAP_REVISION = "52c4146596fef2d1e2402991fbab062ba8021889^"
BOOTSTRAP_IDENTITIES = {
    ("bam_tests", "example_b178_hg19_subset_fast"),
    ("bam_tests", "example_a5c1_hg19_subset_fast"),
    ("bam_tests", "example_66bf_hg19_subset_fast"),
    ("bam_tests", "example_7a61_hg19_subset_fast"),
    ("bam_tests", "example_dfc3_hg19_subset_fast"),
    ("bam_tests", "example_dfc3_GRCh37_fast"),
    ("bam_tests", "example_dfc3_GRCh38_fast"),
    ("bam_tests", "example_dfc3_hg38_ensembl_fast"),
    ("bam_tests", "example_40cf_hg38_subset_fast_gdp_guard"),
    ("advntr_tests", "example_a5c1_hg19_subset_advntr"),
}
BOOTSTRAP_ROUTING_COUNTS = {
    ("bam_tests", "example_b178_hg19_subset_fast"): (14690, 14690, 0, 1),
    ("bam_tests", "example_a5c1_hg19_subset_fast"): (20888, 20888, 0, 3),
    ("bam_tests", "example_66bf_hg19_subset_fast"): (19841, 19841, 0, 3),
    ("bam_tests", "example_7a61_hg19_subset_fast"): (2359, 2359, 0, 8),
    ("bam_tests", "example_dfc3_hg19_subset_fast"): (31596, 31596, 0, 53),
    ("bam_tests", "example_dfc3_GRCh37_fast"): (31593, 31593, 0, 2),
    ("bam_tests", "example_dfc3_GRCh38_fast"): (31603, 31603, 0, 2),
    ("bam_tests", "example_dfc3_hg38_ensembl_fast"): (31603, 31603, 0, 2),
    ("bam_tests", "example_40cf_hg38_subset_fast_gdp_guard"): (3474, 3474, 0, 93),
    ("advntr_tests", "example_a5c1_hg19_subset_advntr"): (20888, 20888, 0, 3),
}
CURRENT_IDENTITIES = BOOTSTRAP_IDENTITIES | {
    ("single_end_bam_tests", "example_b178_hg19_single_end"),
    ("single_end_bam_tests", "example_b178_hg19_single_end_keep"),
    ("single_end_bam_tests", "example_b178_hg19_single_end_delete"),
    ("single_end_bam_tests", "example_b178_hg19_single_end_delete_overrides_keep"),
    ("single_end_bam_tests", "example_b178_hg19_single_end_archive"),
    ("advntr_tests", "example_b178_hg19_bwa_advntr"),
    ("fastq_tests", "example_6449_hg19_subset_fastq_shark"),
    ("fastq_tests", "example_6449_hg19_subset_single_fastq"),
    ("fastq_tests", "example_b178_hg19_subset_paired_fastq_no_shark"),
}


def _module() -> Any:
    """Import the production module while keeping the initial RED collectable."""
    try:
        return importlib.import_module("scripts.integration_compatibility")
    except ModuleNotFoundError:
        return None


def _assertion(value: object, tolerance: dict[str, object] | None = None) -> dict[str, object]:
    return {"value": value, "tolerance": tolerance}


def _contract(*, suite: str = "bam_tests", name: str = "case-a", path: str = "tests/data/input.bam") -> dict:
    return {
        "suite": suite,
        "test_name": name,
        "fixture": None,
        "inputs": [{"path": path, "md5": "a" * 32, "record_digest": None}],
        "reference": {"assembly": "hg19", "path": None, "sha256": None},
        "execution": {
            "threads": 2,
            "log_level": "DEBUG",
            "cli_options": ["--fast-mode", "--keep-intermediates"],
            "modules": [],
        },
        "routing": {
            "records": {"r1": 4, "r2": 4, "other": 0, "single": 1},
            "selected_basenames": ["output_R1.fastq.gz", "output_R2.fastq.gz", "output_single.fastq.gz"],
        },
        "artifacts": {
            "required_present": [
                "summary_report.html",
                "kestrel/kestrel_result.tsv",
                "coverage/coverage_summary.tsv",
                "pipeline_summary.json",
            ],
            "required_absent": [],
            "archive": False,
        },
        "outcomes": {
            "kestrel": {
                "Confidence": _assertion("High_Precision*"),
                "Estimated_Depth_AlternateVariant": _assertion(5),
                "Estimated_Depth_Variant_ActiveRegion": _assertion(100),
                "Depth_Score": _assertion(0.05),
            },
            "advntr": {},
            "cross_match": {},
        },
        "compatibility_since": "2.0.9",
        "provenance_commit": "b" * 40,
    }


def _manifest(*contracts: dict) -> dict:
    return {"schema_version": 1, "contracts": list(contracts)}


def _resources(*paths: str) -> dict:
    return {
        "file_resources": [
            {"local_path": str(Path(path).parent), "filename": Path(path).name, "md5sum": "a" * 32} for path in paths
        ]
    }


def _live_case(contract: dict) -> dict:
    case = {
        "test_name": contract["test_name"],
        "bam": contract["inputs"][0]["path"],
        "reference_assembly": contract["reference"]["assembly"],
        "threads": contract["execution"]["threads"],
        "log_level": contract["execution"]["log_level"],
        "cli_options": [*contract["execution"]["cli_options"]],
        "expected_files": [*contract["artifacts"]["required_present"]],
        "expected_absent": [*contract["artifacts"]["required_absent"]],
        "expected_archive": contract["artifacts"]["archive"],
        "expected_fastq_records": copy.deepcopy(contract["routing"]["records"]),
        "expected_selected_fastqs": [*contract["routing"]["selected_basenames"]],
        "kestrel_assertions": {
            field: value["value"]
            if value["tolerance"] is None
            else {
                "value": value["value"],
                ("tolerance_percentage" if value["tolerance"]["kind"] == "percentage" else "log10_tolerance"): value[
                    "tolerance"
                ]["value"],
            }
            for field, value in contract["outcomes"]["kestrel"].items()
        },
        "advntr_assertions": {field: value["value"] for field, value in contract["outcomes"]["advntr"].items()},
        "cross_match_assertions": {
            "comments": [],
            "data": [{field: str(value["value"]) for field, value in contract["outcomes"]["cross_match"].items()}],
        }
        if contract["outcomes"]["cross_match"]
        else {},
    }
    modules = contract["execution"]["modules"]
    if modules:
        case["cli_options"].extend(["--extra-modules", ",".join(modules)])
    return case


def _advntr_contract() -> dict:
    contract = _contract(suite="advntr_tests", name="case-advntr")
    contract["execution"]["modules"] = ["advntr"]
    contract["outcomes"]["advntr"] = {
        "VID": _assertion("25561"),
        "State": _assertion("I22_2_G_LEN1"),
        "NumberOfSupportingReads": _assertion(11),
        "MeanCoverage": _assertion(153.986111111),
        "Pvalue": _assertion(6.78296229901e-7),
    }
    contract["outcomes"]["cross_match"] = {
        "Kestrel_POS": _assertion("67"),
        "Kestrel_REF": _assertion("G"),
        "Kestrel_ALT": _assertion("GG"),
        "Kestrel_Allele_Change": _assertion("G"),
        "Kestrel_Variant_Type": _assertion("Insertion"),
        "Advntr_POS": _assertion("22"),
        "Advntr_REF": _assertion("T"),
        "Advntr_ALT": _assertion("TG"),
        "Advntr_Allele_Change": _assertion("G"),
        "Advntr_Variant_Type": _assertion("Insertion"),
        "Match": _assertion("Yes"),
    }
    return contract


def test_validate_manifest_accepts_complete_contract() -> None:
    """Catch rejecting the exact, fully resolved contract schema."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()

    indexed = module.validate_manifest(_manifest(contract), _resources("tests/data/input.bam"))

    assert indexed == {("bam_tests", "case-a"): contract}


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda value: value.update(routing=None), "requires routing"),
        (lambda value: value["routing"].update(selected_basenames=[]), "selected routing evidence"),
        (lambda value: value["outcomes"].update(kestrel={}), "Kestrel fields"),
        (lambda value: value["artifacts"].update(required_present=[]), "required artifacts"),
    ],
)
def test_validate_manifest_requires_complete_kestrel_success_content(mutation: Any, message: str) -> None:
    """Catch accepting a success row without routing, value, or artifact evidence."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    mutation(contract)

    with pytest.raises(ValueError, match=message):
        module.validate_manifest(_manifest(contract), _resources("tests/data/input.bam"))


@pytest.mark.parametrize(
    ("outcome", "mutation"),
    [
        ("advntr", lambda value: value.clear()),
        ("cross_match", lambda value: value.clear()),
        ("cross_match", lambda value: value.pop("Match")),
    ],
)
def test_validate_manifest_requires_advntr_specific_outcomes(outcome: str, mutation: Any) -> None:
    """Catch accepting an adVNTR module success without its module-specific value map."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _advntr_contract()
    mutation(contract["outcomes"][outcome])

    with pytest.raises(ValueError, match=f"{outcome} fields"):
        module.validate_manifest(_manifest(contract), _resources("tests/data/input.bam"))


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda value: value.update(extra=True), "extra keys"),
        (lambda value: value.pop("reference"), "missing keys"),
        (lambda value: value["execution"].update(threads=True), "threads"),
        (lambda value: value["execution"].update(threads=0), "threads"),
        (lambda value: value["inputs"][0].update(path="/tmp/input.bam"), "relative"),
        (lambda value: value["inputs"][0].update(path="tests/../input.bam"), "relative"),
        (lambda value: value["execution"].update(cli_options=["--fast-mode", "--fast-mode"]), "duplicate"),
        (
            lambda value: value["artifacts"].update(required_present=["summary_report.html", "summary_report.html"]),
            "duplicate",
        ),
        (
            lambda value: value["outcomes"]["kestrel"].update(
                Depth_Score=_assertion(0.1, {"kind": "percentage", "value": float("inf")})
            ),
            "finite",
        ),
    ],
)
def test_validate_manifest_rejects_malformed_contract(mutation: Any, message: str) -> None:
    """Catch accepting schema drift, unsafe paths, Booleans, duplicates, or unbounded tolerances."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    mutation(contract)

    with pytest.raises(ValueError, match=message):
        module.validate_manifest(_manifest(contract), _resources("tests/data/input.bam"))


def test_validate_manifest_rejects_duplicate_suite_name_identity() -> None:
    """Catch collapsing or silently replacing duplicate identities."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()

    with pytest.raises(ValueError, match="duplicate identity"):
        module.validate_manifest(_manifest(contract, copy.deepcopy(contract)), _resources("tests/data/input.bam"))


def test_validate_manifest_allows_same_input_in_distinct_suites() -> None:
    """Catch treating an input path, rather than suite/name, as identity."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    bam = _contract()
    advntr = _advntr_contract()

    indexed = module.validate_manifest(_manifest(bam, advntr), _resources("tests/data/input.bam"))

    assert set(indexed) == {("bam_tests", "case-a"), ("advntr_tests", "case-advntr")}


def test_validate_manifest_rejects_unresolved_or_wrong_resource_digest() -> None:
    """Catch accepting an input whose archive identity is absent or differs."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"

    with pytest.raises(ValueError, match="unresolved resource digest"):
        module.validate_manifest(_manifest(_contract()), {"file_resources": []})
    resources = _resources("tests/data/input.bam")
    resources["file_resources"][0]["md5sum"] = "c" * 32
    with pytest.raises(ValueError, match="resource digest mismatch"):
        module.validate_manifest(_manifest(_contract()), resources)


def test_validate_manifest_binds_generated_fixture_to_source_and_recipe() -> None:
    """Catch treating a generated BAM as an unregistered archive member or unpinned recipe."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract(suite="single_end_bam_tests", path="tests/data/source.bam")
    contract["fixture"] = {
        "name": "derived-single",
        "kind": "single_end_bam",
        "source_path": "tests/data/source.bam",
        "output_path": "tests/data/derived/single.bam",
        "recipe_digest": "2527b10c934fd45d51c184b82898687d67972b6dc089291d2e171a33ee12487f",
    }
    resources = _resources("tests/data/source.bam")
    resources["derived_fixtures"] = [
        {
            "name": "derived-single",
            "kind": "single_end_bam",
            "source_bam": "tests/data/source.bam",
            "output_bam": "tests/data/derived/single.bam",
        }
    ]

    module.validate_manifest(_manifest(contract), resources)
    contract["fixture"]["recipe_digest"] = "0" * 64
    with pytest.raises(ValueError, match="fixture recipe digest mismatch"):
        module.validate_manifest(_manifest(contract), resources)


def test_validate_manifest_resolves_derived_record_and_reference_digests() -> None:
    """Catch accepting an unregistered derived-record or reference fingerprint."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract(path="tests/data/derived/input.cram")
    contract["inputs"][0] = {
        "path": "tests/data/derived/input.cram",
        "md5": None,
        "record_digest": "d" * 64,
    }
    contract["reference"] = {"assembly": "hg19", "path": "reference/ref.fa", "sha256": "e" * 64}
    resources = {
        "file_resources": [],
        "derived_fixtures": [
            {
                "output_cram": "tests/data/derived/input.cram",
                "fixture_record_digest": "d" * 64,
                "reference_fasta": "reference/ref.fa",
                "fixture_reference_sha256": "e" * 64,
            }
        ],
    }

    module.validate_manifest(_manifest(contract), resources)
    resources["derived_fixtures"][0]["fixture_reference_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="reference digest mismatch"):
        module.validate_manifest(_manifest(contract), resources)


def test_check_compatibility_rejects_base_deletion_and_mutation() -> None:
    """Catch permitting removal or any in-place change to an established success."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    resources = _resources("tests/data/input.bam")
    live = {"integration_tests": {"bam_tests": [_live_case(contract)]}}

    with pytest.raises(ValueError, match="removed"):
        module.check_compatibility(_manifest(contract), _manifest(), live, resources)
    changed = copy.deepcopy(contract)
    changed["outcomes"]["kestrel"]["Confidence"]["value"] = "Negative"
    with pytest.raises(ValueError, match="mutated"):
        module.check_compatibility(_manifest(contract), _manifest(changed), live, resources)


@pytest.mark.parametrize("field", ["threads", "log_level"])
def test_check_compatibility_requires_explicit_live_invocation_metadata(field: str) -> None:
    """Catch restoring an implicit invocation default that can hide Task 5 declaration drift."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    case = _live_case(contract)
    case.pop(field)
    live = {"integration_tests": {"bam_tests": [case]}}

    with pytest.raises(ValueError, match=f"live case .* {field}"):
        module.check_compatibility(_manifest(contract), _manifest(contract), live, _resources("tests/data/input.bam"))


def test_check_compatibility_accepts_explicit_task5_invocation_metadata() -> None:
    """Catch accidentally rejecting the planned explicit two-thread DEBUG success declaration."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    live = {"integration_tests": {"bam_tests": [_live_case(contract)]}}

    module.check_compatibility(_manifest(contract), _manifest(contract), live, _resources("tests/data/input.bam"))


def test_check_compatibility_accepts_task5_cross_match_result_shape() -> None:
    """Catch requiring a fabricated overall-match key instead of the exact summary row."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _advntr_contract()
    live = {"integration_tests": {"advntr_tests": [_live_case(contract)]}}

    module.check_compatibility(_manifest(contract), _manifest(contract), live, _resources("tests/data/input.bam"))


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda value: value.update(expected_exit_code=False), "expected_exit_code"),
        (lambda value: value.update(expected_exit_code="0"), "expected_exit_code"),
        (lambda value: value.update(threads=True), "threads"),
        (lambda value: value.pop("test_name"), "test_name"),
        (lambda value: value.update(test_name=None), "test_name"),
        (lambda value: value.update(bam="/tmp/input.bam"), "normalized relative path"),
        (lambda value: value.update(cli_options=[1]), "non-empty string"),
        (
            lambda value: value["cli_options"].extend(["--extra-modules", "shark,shark"]),
            "duplicate module",
        ),
    ],
)
def test_check_compatibility_rejects_malformed_live_declarations(mutation: Any, message: str) -> None:
    """Catch Boolean coercion, unnamed identities, and normalized-away module duplication."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    case = _live_case(contract)
    mutation(case)
    live = {"integration_tests": {"bam_tests": [case]}}

    with pytest.raises(ValueError, match=message):
        module.check_compatibility(_manifest(contract), _manifest(contract), live, _resources("tests/data/input.bam"))


@pytest.mark.parametrize(
    "path",
    [
        ("suite",),
        ("test_name",),
        ("inputs", 0, "path"),
        ("reference", "assembly"),
        ("execution", "cli_options"),
        ("execution", "modules"),
        ("artifacts", "required_present"),
        ("outcomes", "kestrel", "Confidence", "value"),
        ("outcomes", "kestrel", "Confidence", "tolerance"),
    ],
)
def test_check_compatibility_rejects_live_contract_mutation(path: tuple[object, ...]) -> None:
    """Catch normalizing a regression by changing any protected live declaration field."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    live_contract = copy.deepcopy(contract)
    target: Any = live_contract
    for part in path[:-1]:
        target = target[part]
    leaf = path[-1]
    current = target[leaf]
    if isinstance(current, list):
        target[leaf] = [*current, "changed"]
    elif current is None:
        target[leaf] = {"kind": "percentage", "value": 5}
    else:
        target[leaf] = f"changed-{current}"
    live = {"integration_tests": {live_contract["suite"]: [_live_case(live_contract)]}}

    message = "unsupported" if path == ("suite",) else r"does not (?:resolve|match live declaration)"
    with pytest.raises(ValueError, match=message):
        module.check_compatibility(_manifest(contract), _manifest(contract), live, _resources("tests/data/input.bam"))


def test_check_compatibility_rejects_success_without_contract_and_accepts_append() -> None:
    """Catch one-way validation that overlooks newly added real successes."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    first = _contract()
    second = _contract(name="case-b", path="tests/data/second.bam")
    live = {"integration_tests": {"bam_tests": [_live_case(first), _live_case(second)]}}
    resources = _resources("tests/data/input.bam", "tests/data/second.bam")

    with pytest.raises(ValueError, match="missing from manifest"):
        module.check_compatibility(_manifest(first), _manifest(first), live, resources)
    module.check_compatibility(_manifest(first), _manifest(first, second), live, resources)


def test_check_compatibility_rejects_missing_bootstrap_evidence_and_duplicate_live_identity() -> None:
    """Catch unchecked empty bootstrap and ambiguous identity resolution."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    contract = _contract()
    resources = _resources("tests/data/input.bam")
    live_case = _live_case(contract)
    live = {"integration_tests": {"bam_tests": [live_case]}}

    with pytest.raises(ValueError, match="requires authoritative bootstrap history"):
        module.check_compatibility(None, _manifest(contract), live, resources)
    duplicate_live = {"integration_tests": {"bam_tests": [live_case, copy.deepcopy(live_case)]}}
    with pytest.raises(ValueError, match="resolve exactly once"):
        module.check_compatibility(_manifest(contract), _manifest(contract), duplicate_live, resources)


def test_manifest_contains_every_required_inactive_success_identity() -> None:
    """Catch omitting any historical, single-end, FASTQ, or pure-paired adVNTR success from the seed."""
    manifest = json.loads(Path("tests/compatibility/real_success_baseline.json").read_text())

    assert {(row["suite"], row["test_name"]) for row in manifest["contracts"]} == CURRENT_IDENTITIES


def test_bootstrap_seed_is_reconstructed_from_authoritative_git_history() -> None:
    """Catch omitting, renaming, or mutating any one of the ten pre-regression successes."""
    module = _module()
    assert module is not None, "integration compatibility module is not implemented"
    historical = json.loads(
        subprocess.run(
            ["git", "show", f"{BOOTSTRAP_REVISION}:tests/test_data_config.json"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )
    manifest = json.loads(Path("tests/compatibility/real_success_baseline.json").read_text())
    indexed = module.validate_bootstrap_manifest(
        manifest, historical, json.loads(Path("tests/test_data_config.json").read_text())
    )

    assert set(indexed) >= BOOTSTRAP_IDENTITIES
    assert {
        (suite, case["test_name"])
        for suite in ("bam_tests", "advntr_tests")
        for case in historical["integration_tests"][suite]
        if case.get("expected_exit_code", 0) == 0
    } == BOOTSTRAP_IDENTITIES

    missing = copy.deepcopy(manifest)
    missing["contracts"] = [
        row
        for row in missing["contracts"]
        if (row["suite"], row["test_name"]) != ("bam_tests", "example_b178_hg19_subset_fast")
    ]
    with pytest.raises(ValueError, match="missing authoritative identities"):
        module.validate_bootstrap_manifest(
            missing, historical, json.loads(Path("tests/test_data_config.json").read_text())
        )

    mutated_history = copy.deepcopy(historical)
    mutated_history["integration_tests"]["bam_tests"][0]["test_name"] = "renamed"
    with pytest.raises(ValueError, match="exact ten-ID seed"):
        module.validate_bootstrap_manifest(
            manifest, mutated_history, json.loads(Path("tests/test_data_config.json").read_text())
        )

    for identity, expected in BOOTSTRAP_ROUTING_COUNTS.items():
        records = indexed[identity]["routing"]["records"]
        assert tuple(records[key] for key in ("r1", "r2", "other", "single")) == expected
    mutated_routing = copy.deepcopy(manifest)
    mutated_routing["contracts"][0]["routing"]["records"]["r1"] += 1
    with pytest.raises(ValueError, match="authoritative routing"):
        module.validate_bootstrap_manifest(
            mutated_routing, historical, json.loads(Path("tests/test_data_config.json").read_text())
        )
