"""Validate immutable real-integration success contracts.

This module is deliberately independent of pytest and CI.  Its public core consumes
already-loaded JSON objects so callers and unit tests can provide authoritative inputs
without hidden filesystem or Git access.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import PurePosixPath
from typing import Any

Manifest = dict[str, Any]
Identity = tuple[str, str]

SCHEMA_VERSION = 1
BOOTSTRAP_REVISION = "52c4146596fef2d1e2402991fbab062ba8021889^"
BOOTSTRAP_COMMIT = "47ce0b1168210e55841d07f280ba8134a073e96e"
BOOTSTRAP_IDENTITIES: frozenset[Identity] = frozenset(
    {
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
)
_BOOTSTRAP_ROUTING_COUNTS: dict[Identity, tuple[int, int, int, int]] = {
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
_BOOTSTRAP_SELECTED_BASENAMES = ["output_R1.fastq.gz", "output_R2.fastq.gz", "output_single.fastq.gz"]

_TOP_KEYS = {"schema_version", "contracts"}
_CONTRACT_KEYS = {
    "suite",
    "test_name",
    "fixture",
    "inputs",
    "reference",
    "execution",
    "routing",
    "artifacts",
    "outcomes",
    "compatibility_since",
    "provenance_commit",
}
_INPUT_KEYS = {"path", "md5", "record_digest"}
_FIXTURE_KEYS = {"name", "kind", "source_path", "output_path", "recipe_digest"}
_REFERENCE_KEYS = {"assembly", "path", "sha256"}
_EXECUTION_KEYS = {"threads", "log_level", "cli_options", "modules"}
_ROUTING_KEYS = {"records", "selected_basenames"}
_COUNT_KEYS = {"r1", "r2", "other", "single"}
_ARTIFACT_KEYS = {"required_present", "required_absent", "archive"}
_OUTCOME_KEYS = {"kestrel", "advntr", "cross_match"}
_ASSERTION_KEYS = {"value", "tolerance"}
_TOLERANCE_KEYS = {"kind", "value"}
_SHA_RE = re.compile(r"[0-9a-f]{40}")
_MD5_RE = re.compile(r"[0-9a-f]{32}")
_SHA256_RE = re.compile(r"[0-9a-f]{64}")


def _mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _exact_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing:
        raise ValueError(f"{label} has missing keys: {missing}")
    if extra:
        raise ValueError(f"{label} has extra keys: {extra}")


def _string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a non-empty string")
    return value


def _integer(value: object, label: str, *, positive: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or (positive and value <= 0) or value < 0:
        qualifier = "positive " if positive else "non-negative "
        raise ValueError(f"{label} must be a {qualifier}integer")
    return value


def _relative_path(value: object, label: str) -> str:
    text = _string(value, label)
    path = PurePosixPath(text)
    if path.is_absolute() or ".." in path.parts or "." in path.parts or text != path.as_posix():
        raise ValueError(f"{label} must be a normalized relative path")
    return text


def _string_list(value: object, label: str, *, paths: bool = False, basenames: bool = False) -> list[str]:
    if not isinstance(value, list):
        raise ValueError(f"{label} must be a list")
    result: list[str] = []
    for index, item in enumerate(value):
        text = _relative_path(item, f"{label}[{index}]") if paths else _string(item, f"{label}[{index}]")
        if basenames and PurePosixPath(text).name != text:
            raise ValueError(f"{label}[{index}] must be a basename")
        result.append(text)
    if len(set(result)) != len(result):
        raise ValueError(f"{label} contains a duplicate")
    return result


def _resource_digests(resource_config: dict[str, Any]) -> dict[str, str]:
    resources = resource_config.get("file_resources", [])
    if not isinstance(resources, list):
        raise ValueError("resource config file_resources must be a list")
    result: dict[str, str] = {}
    for index, resource_value in enumerate(resources):
        resource = _mapping(resource_value, f"file_resources[{index}]")
        try:
            path = (PurePosixPath(str(resource["local_path"])) / str(resource["filename"])).as_posix()
            digest = str(resource["md5sum"])
        except KeyError as error:
            raise ValueError(f"file_resources[{index}] is missing {error.args[0]}") from error
        if path in result and result[path] != digest:
            raise ValueError(f"resource config has inconsistent digests for {path}")
        result[path] = digest
    return result


def _derived_record_digests(resource_config: dict[str, Any]) -> dict[str, str]:
    result: dict[str, str] = {}
    fixtures = resource_config.get("derived_fixtures", [])
    if not isinstance(fixtures, list):
        raise ValueError("resource config derived_fixtures must be a list")
    for fixture_value in fixtures:
        fixture = _mapping(fixture_value, "derived fixture")
        path = fixture.get("output_cram") or fixture.get("output_bam")
        digest = fixture.get("fixture_record_digest") or fixture.get("record_digest")
        if isinstance(path, str) and isinstance(digest, str):
            result[path] = digest
    return result


def _reference_digests(resource_config: dict[str, Any]) -> dict[str, str]:
    """Return registered reference-file SHA-256 digests."""
    result: dict[str, str] = {}
    fixtures = resource_config.get("derived_fixtures", [])
    if not isinstance(fixtures, list):
        raise ValueError("resource config derived_fixtures must be a list")
    for fixture_value in fixtures:
        fixture = _mapping(fixture_value, "derived fixture")
        path = fixture.get("reference_fasta")
        digest = fixture.get("fixture_reference_sha256")
        if not isinstance(path, str) or not isinstance(digest, str):
            continue
        normalized_path = _relative_path(path, "derived fixture reference_fasta")
        if _SHA256_RE.fullmatch(digest) is None:
            raise ValueError(f"derived fixture reference digest for {normalized_path} must be a lowercase SHA-256")
        if normalized_path in result and result[normalized_path] != digest:
            raise ValueError(f"resource config has inconsistent reference digests for {normalized_path}")
        result[normalized_path] = digest
    return result


def _validate_input(value: object, label: str, resources: dict[str, str], derived: dict[str, str]) -> None:
    item = _mapping(value, label)
    _exact_keys(item, _INPUT_KEYS, label)
    path = _relative_path(item["path"], f"{label}.path")
    md5 = item["md5"]
    record_digest = item["record_digest"]
    if (md5 is None) == (record_digest is None):
        raise ValueError(f"{label} must declare exactly one resource digest")
    if md5 is not None:
        if not isinstance(md5, str) or _MD5_RE.fullmatch(md5) is None:
            raise ValueError(f"{label}.md5 must be a lowercase MD5")
        if path not in resources:
            raise ValueError(f"{label} has unresolved resource digest: {path}")
        if resources[path] != md5:
            raise ValueError(f"{label} resource digest mismatch for {path}")
    else:
        if not isinstance(record_digest, str) or _SHA256_RE.fullmatch(record_digest) is None:
            raise ValueError(f"{label}.record_digest must be a lowercase SHA-256")
        if path not in derived:
            raise ValueError(f"{label} has unresolved resource digest: {path}")
        if derived[path] != record_digest:
            raise ValueError(f"{label} resource digest mismatch for {path}")


def _validate_reference(value: object, label: str, reference_digests: dict[str, str]) -> None:
    reference = _mapping(value, label)
    _exact_keys(reference, _REFERENCE_KEYS, label)
    _string(reference["assembly"], f"{label}.assembly")
    path, digest = reference["path"], reference["sha256"]
    if (path is None) != (digest is None):
        raise ValueError(f"{label}.path and sha256 must be declared together")
    if path is not None:
        normalized_path = _relative_path(path, f"{label}.path")
        if not isinstance(digest, str) or _SHA256_RE.fullmatch(digest) is None:
            raise ValueError(f"{label}.sha256 must be a lowercase SHA-256")
        registered_digest = reference_digests.get(normalized_path)
        if registered_digest is None:
            raise ValueError(f"unresolved reference digest for {normalized_path}")
        if registered_digest != digest:
            raise ValueError(f"reference digest mismatch for {normalized_path}")


def _fixture_recipe_digest(fixture: dict[str, Any]) -> str:
    encoded = json.dumps(fixture, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _validate_fixture(value: object, label: str, resource_config: dict[str, Any], input_paths: list[str]) -> None:
    if value is None:
        return
    fixture = _mapping(value, label)
    _exact_keys(fixture, _FIXTURE_KEYS, label)
    for key in ("name", "kind"):
        _string(fixture[key], f"{label}.{key}")
    source = _relative_path(fixture["source_path"], f"{label}.source_path")
    _relative_path(fixture["output_path"], f"{label}.output_path")
    if source not in input_paths:
        raise ValueError(f"{label}.source_path must identify a contract input")
    digest = fixture["recipe_digest"]
    if not isinstance(digest, str) or _SHA256_RE.fullmatch(digest) is None:
        raise ValueError(f"{label}.recipe_digest must be a lowercase SHA-256")
    matches = [
        item
        for item in resource_config.get("derived_fixtures", [])
        if isinstance(item, dict) and item.get("name") == fixture["name"]
    ]
    if len(matches) != 1:
        raise ValueError(f"{label} does not resolve exactly once")
    registered = matches[0]
    if registered.get("kind") != fixture["kind"] or registered.get("source_bam") != source:
        raise ValueError(f"{label} does not match registered fixture identity")
    output = registered.get("output_bam") or registered.get("output_cram")
    if output != fixture["output_path"]:
        raise ValueError(f"{label} does not match registered fixture output")
    if _fixture_recipe_digest(registered) != digest:
        raise ValueError(f"{label} fixture recipe digest mismatch")


def _validate_execution(value: object, label: str) -> None:
    execution = _mapping(value, label)
    _exact_keys(execution, _EXECUTION_KEYS, label)
    _integer(execution["threads"], f"{label}.threads", positive=True)
    _string(execution["log_level"], f"{label}.log_level")
    _string_list(execution["cli_options"], f"{label}.cli_options")
    modules = _string_list(execution["modules"], f"{label}.modules")
    if modules != sorted(modules):
        raise ValueError(f"{label}.modules must be a normalized module set")


def _validate_routing(value: object, label: str) -> None:
    if value is None:
        return
    routing = _mapping(value, label)
    _exact_keys(routing, _ROUTING_KEYS, label)
    records = _mapping(routing["records"], f"{label}.records")
    _exact_keys(records, _COUNT_KEYS, f"{label}.records")
    for key in sorted(_COUNT_KEYS):
        _integer(records[key], f"{label}.records.{key}")
    _string_list(routing["selected_basenames"], f"{label}.selected_basenames", basenames=True)


def _validate_artifacts(value: object, label: str) -> None:
    artifacts = _mapping(value, label)
    _exact_keys(artifacts, _ARTIFACT_KEYS, label)
    present = _string_list(artifacts["required_present"], f"{label}.required_present", paths=True)
    absent = _string_list(artifacts["required_absent"], f"{label}.required_absent", paths=True)
    if set(present) & set(absent):
        raise ValueError(f"{label} present and absent artifacts overlap")
    if not isinstance(artifacts["archive"], bool):
        raise ValueError(f"{label}.archive must be a Boolean")


def _validate_assertions(value: object, label: str) -> None:
    assertions = _mapping(value, label)
    for field, assertion_value in assertions.items():
        _string(field, f"{label} field")
        assertion = _mapping(assertion_value, f"{label}.{field}")
        _exact_keys(assertion, _ASSERTION_KEYS, f"{label}.{field}")
        expected = assertion["value"]
        if expected is None or isinstance(expected, (dict, list)):
            raise ValueError(f"{label}.{field}.value must be a scalar")
        tolerance_value = assertion["tolerance"]
        if tolerance_value is None:
            continue
        tolerance = _mapping(tolerance_value, f"{label}.{field}.tolerance")
        _exact_keys(tolerance, _TOLERANCE_KEYS, f"{label}.{field}.tolerance")
        if tolerance["kind"] not in {"percentage", "log10"}:
            raise ValueError(f"{label}.{field}.tolerance kind is unsupported")
        number = tolerance["value"]
        if isinstance(number, bool) or not isinstance(number, (int, float)) or not math.isfinite(number):
            raise ValueError(f"{label}.{field}.tolerance value must be finite")
        if number < 0 or (tolerance["kind"] == "percentage" and number > 100):
            raise ValueError(f"{label}.{field}.tolerance is unbounded")


def _validate_contract(
    value: object,
    index: int,
    resources: dict[str, str],
    derived: dict[str, str],
    reference_digests: dict[str, str],
    resource_config: dict[str, Any],
) -> tuple[Identity, dict[str, Any]]:
    label = f"contracts[{index}]"
    contract = _mapping(value, label)
    _exact_keys(contract, _CONTRACT_KEYS, label)
    suite = _string(contract["suite"], f"{label}.suite")
    test_name = _string(contract["test_name"], f"{label}.test_name")
    inputs = contract["inputs"]
    if not isinstance(inputs, list) or not inputs:
        raise ValueError(f"{label}.inputs must be a non-empty list")
    paths: list[str] = []
    for input_index, item in enumerate(inputs):
        _validate_input(item, f"{label}.inputs[{input_index}]", resources, derived)
        paths.append(item["path"])
    if len(set(paths)) != len(paths):
        raise ValueError(f"{label}.inputs contains a duplicate")
    _validate_fixture(contract["fixture"], f"{label}.fixture", resource_config, paths)
    _validate_reference(contract["reference"], f"{label}.reference", reference_digests)
    _validate_execution(contract["execution"], f"{label}.execution")
    _validate_routing(contract["routing"], f"{label}.routing")
    _validate_artifacts(contract["artifacts"], f"{label}.artifacts")
    outcomes = _mapping(contract["outcomes"], f"{label}.outcomes")
    _exact_keys(outcomes, _OUTCOME_KEYS, f"{label}.outcomes")
    for outcome in sorted(_OUTCOME_KEYS):
        _validate_assertions(outcomes[outcome], f"{label}.outcomes.{outcome}")
    _string(contract["compatibility_since"], f"{label}.compatibility_since")
    provenance = contract["provenance_commit"]
    if not isinstance(provenance, str) or _SHA_RE.fullmatch(provenance) is None:
        raise ValueError(f"{label}.provenance_commit must be a full lowercase Git SHA")
    return (suite, test_name), contract


def validate_manifest(manifest: Manifest, resource_config: dict[str, Any]) -> dict[Identity, dict[str, Any]]:
    """Validate and index one strict compatibility manifest.

    Args:
        manifest: Parsed manifest JSON.
        resource_config: Parsed test resource configuration.

    Returns:
        Contracts indexed by ``(suite, test_name)``.

    Raises:
        ValueError: If the manifest or any resource binding is invalid.
    """
    root = _mapping(manifest, "manifest")
    _exact_keys(root, _TOP_KEYS, "manifest")
    if isinstance(root["schema_version"], bool) or root["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"manifest schema_version must be {SCHEMA_VERSION}")
    contracts = root["contracts"]
    if not isinstance(contracts, list):
        raise ValueError("manifest contracts must be a list")
    resources = _resource_digests(resource_config)
    derived = _derived_record_digests(resource_config)
    reference_digests = _reference_digests(resource_config)
    result: dict[Identity, dict[str, Any]] = {}
    for index, value in enumerate(contracts):
        identity, contract = _validate_contract(value, index, resources, derived, reference_digests, resource_config)
        if identity in result:
            raise ValueError(f"manifest contains duplicate identity {identity}")
        result[identity] = contract
    return result


def _assertions_from_live(value: object) -> dict[str, dict[str, object]]:
    if not isinstance(value, dict):
        return {}
    result: dict[str, dict[str, object]] = {}
    for field, raw in value.items():
        if isinstance(raw, dict) and "value" in raw:
            tolerance = None
            if "tolerance_percentage" in raw:
                tolerance = {"kind": "percentage", "value": raw["tolerance_percentage"]}
            elif "log10_tolerance" in raw:
                tolerance = {"kind": "log10", "value": raw["log10_tolerance"]}
            result[field] = {"value": raw["value"], "tolerance": tolerance}
        else:
            result[field] = {"value": raw, "tolerance": None}
    return result


def _split_options(case: dict[str, Any], suite: str) -> tuple[list[str], list[str]]:
    options = case.get("cli_options", [])
    if not isinstance(options, list):
        raise ValueError(f"live case {case.get('test_name')} cli_options must be a list")
    modules = ["advntr"] if suite == "advntr_tests" else []
    normalized: list[str] = []
    index = 0
    while index < len(options):
        option = options[index]
        if option == "--extra-modules":
            if index + 1 >= len(options) or not isinstance(options[index + 1], str):
                raise ValueError(f"live case {case.get('test_name')} has malformed --extra-modules")
            modules.extend(module for module in options[index + 1].split(",") if module not in modules)
            index += 2
        else:
            normalized.append(option)
            index += 1
    return normalized, sorted(modules)


def _live_fixture(case: dict[str, Any], resource_config: dict[str, Any]) -> dict[str, str] | None:
    fixture_name = case.get("fixture_name")
    if fixture_name is None:
        return None
    matches = [
        item
        for item in resource_config.get("derived_fixtures", [])
        if isinstance(item, dict) and item.get("name") == fixture_name
    ]
    if len(matches) != 1:
        return None
    registered = matches[0]
    output = registered.get("output_bam") or registered.get("output_cram")
    source = registered.get("source_bam")
    if not isinstance(output, str) or not isinstance(source, str):
        return None
    return {
        "name": str(registered.get("name")),
        "kind": str(registered.get("kind")),
        "source_path": source,
        "output_path": output,
        "recipe_digest": _fixture_recipe_digest(registered),
    }


def _live_inputs(case: dict[str, Any], suite: str, resource_config: dict[str, Any]) -> list[str]:
    if suite == "fastq_tests":
        return [case[key] for key in ("fastq1", "fastq2") if isinstance(case.get(key), str) and case[key]]
    for key in ("bam", "cram"):
        if isinstance(case.get(key), str):
            return [case[key]]
    fixture = _live_fixture(case, resource_config)
    return [fixture["source_path"]] if fixture else []


def _default_artifacts(case: dict[str, Any], suite: str) -> list[str]:
    if suite == "advntr_tests":
        defaults = ["summary_report.html", "kestrel/kestrel_result.tsv", "advntr/output_adVNTR_result.tsv"]
    elif suite == "fastq_tests":
        defaults = list(case.get("expected_files", []))
    else:
        defaults = [
            "summary_report.html",
            "kestrel/kestrel_result.tsv",
            "coverage/coverage_summary.tsv",
            "pipeline_summary.json",
        ]
        if case.get("check_igv_report"):
            defaults.append("igv_report.html")
    for path in case.get("expected_present", []):
        if path not in defaults:
            defaults.append(path)
    return defaults


def _live_projection(case: dict[str, Any], suite: str, resource_config: dict[str, Any]) -> dict[str, Any]:
    options, modules = _split_options(case, suite)
    records = case.get("expected_fastq_records") or case.get("expected_mixed_fastq_records")
    selected = case.get("expected_selected_fastqs")
    routing = None if records is None else {"records": records, "selected_basenames": selected or []}
    return {
        "suite": suite,
        "test_name": case.get("test_name"),
        "fixture": _live_fixture(case, resource_config),
        "input_paths": _live_inputs(case, suite, resource_config),
        "reference": {
            "assembly": case.get("reference_assembly", "hg19"),
            "path": case.get("reference_fasta"),
            "sha256": case.get("fixture_reference_sha256"),
        },
        "execution": {
            "threads": case.get("threads", 4),
            "log_level": case.get("log_level", "DEBUG"),
            "cli_options": options,
            "modules": modules,
        },
        "routing": routing,
        "artifacts": {
            "required_present": _default_artifacts(case, suite),
            "required_absent": case.get("expected_absent", []),
            "archive": case.get("expected_archive", "--archive-results" in options),
        },
        "outcomes": {
            "kestrel": _assertions_from_live(case.get("kestrel_assertions", {})),
            "advntr": _assertions_from_live(case.get("advntr_assertions", {})),
            "cross_match": _assertions_from_live(case.get("cross_match_assertions", {})),
        },
    }


def _contract_projection(contract: dict[str, Any]) -> dict[str, Any]:
    return {
        "suite": contract["suite"],
        "test_name": contract["test_name"],
        "fixture": contract["fixture"],
        "input_paths": [item["path"] for item in contract["inputs"]],
        "reference": contract["reference"],
        "execution": contract["execution"],
        "routing": contract["routing"],
        "artifacts": contract["artifacts"],
        "outcomes": contract["outcomes"],
    }


def _live_successes(test_config: dict[str, Any]) -> dict[Identity, list[dict[str, Any]]]:
    suites = test_config.get("integration_tests", {})
    if not isinstance(suites, dict):
        raise ValueError("live test config integration_tests must be an object")
    result: dict[Identity, list[dict[str, Any]]] = {}
    for suite, cases in suites.items():
        if not isinstance(cases, list):
            raise ValueError(f"live suite {suite} must be a list")
        for case_value in cases:
            case = _mapping(case_value, f"live suite {suite} case")
            if case.get("expected_exit_code", 0) == 0:
                identity = (suite, str(case.get("test_name")))
                result.setdefault(identity, []).append(case)
    return result


def validate_bootstrap_manifest(
    manifest: Manifest,
    historical_test_config: dict[str, Any],
    resource_config: dict[str, Any],
) -> dict[Identity, dict[str, Any]]:
    """Validate the missing-base bootstrap against authoritative pre-regression history.

    Args:
        manifest: Candidate current manifest.
        historical_test_config: ``tests/test_data_config.json`` loaded from the pinned revision.
        resource_config: Current resource registry used to resolve immutable inputs.

    Returns:
        The validated manifest index.

    Raises:
        ValueError: If the exact ten-case seed or any historical declaration differs.
    """
    indexed = validate_manifest(manifest, resource_config)
    historical = _live_successes(historical_test_config)
    historical_ids = {identity for identity in historical if identity[0] in {"bam_tests", "advntr_tests"}}
    if historical_ids != BOOTSTRAP_IDENTITIES:
        raise ValueError("pinned bootstrap history does not contain the exact ten-ID seed")
    if not BOOTSTRAP_IDENTITIES.issubset(indexed):
        missing = sorted(BOOTSTRAP_IDENTITIES - indexed.keys())
        raise ValueError(f"bootstrap manifest is missing authoritative identities: {missing}")
    for identity in BOOTSTRAP_IDENTITIES:
        contract = indexed[identity]
        if contract["provenance_commit"] != BOOTSTRAP_COMMIT:
            raise ValueError(f"bootstrap contract {identity} has wrong provenance commit")
        r1, r2, other, single = _BOOTSTRAP_ROUTING_COUNTS[identity]
        authoritative_routing = {
            "records": {"r1": r1, "r2": r2, "other": other, "single": single},
            "selected_basenames": _BOOTSTRAP_SELECTED_BASENAMES,
        }
        if contract["routing"] != authoritative_routing:
            raise ValueError(f"bootstrap contract {identity} mutates authoritative routing")
        case = historical[identity][0]
        projection = _contract_projection(contract)
        live = _live_projection(case, identity[0], resource_config)
        for field in ("suite", "test_name", "input_paths", "reference", "execution", "artifacts", "outcomes"):
            if projection[field] != live[field]:
                raise ValueError(f"bootstrap contract {identity} mutates historical field {field}")
    return indexed


def check_compatibility(
    base_manifest: Manifest | None,
    current_manifest: Manifest,
    live_test_config: dict[str, Any],
    resource_config: dict[str, Any],
    *,
    historical_test_config: dict[str, Any] | None = None,
) -> None:
    """Check schema, append-only history, and bidirectional live declarations.

    Args:
        base_manifest: Manifest loaded from the event base, or ``None`` for the one-time bootstrap.
        current_manifest: Manifest from the current checkout.
        live_test_config: Current integration declarations.
        resource_config: Current immutable resource registry.
        historical_test_config: Pinned historical declarations, required when ``base_manifest`` is absent.

    Raises:
        ValueError: If any compatibility invariant fails.
    """
    current = validate_manifest(current_manifest, resource_config)
    if base_manifest is None:
        if historical_test_config is None:
            raise ValueError("missing base manifest requires authoritative bootstrap history")
        validate_bootstrap_manifest(current_manifest, historical_test_config, resource_config)
    else:
        base = validate_manifest(base_manifest, resource_config)
        removed = sorted(base.keys() - current.keys())
        if removed:
            raise ValueError(f"compatibility contracts were removed: {removed}")
        mutated = [identity for identity in base if base[identity] != current[identity]]
        if mutated:
            raise ValueError(f"compatibility contracts were mutated: {mutated}")

    successes = _live_successes(live_test_config)
    duplicate = sorted(identity for identity, cases in successes.items() if len(cases) != 1)
    if duplicate:
        raise ValueError(f"manifest identities do not resolve exactly once: {duplicate}")
    for identity, contract in current.items():
        cases = successes.get(identity, [])
        if len(cases) != 1:
            raise ValueError(f"contract {identity} does not resolve exactly once to an exit-0 live case")
        if _contract_projection(contract) != _live_projection(cases[0], identity[0], resource_config):
            raise ValueError(f"contract {identity} does not match live declaration")
    missing = sorted(successes.keys() - current.keys())
    if missing:
        raise ValueError(f"qualifying live successes are missing from manifest: {missing}")
