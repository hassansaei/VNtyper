"""Schema, identity, and live-policy tests for the BLE001 adapter."""

import ast
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from scripts.ble001_policy import (
    Diagnostic,
    FailOpenPolicy,
    HandlerPolicy,
    Measurement,
    Policy,
    enclosing_symbol,
    load_policy,
    measure_ble001,
    read_ruff_paths,
    validate_policy,
)
from scripts.ble001_policy_validation import behavior_node_error

pytestmark = pytest.mark.unit
REPO_ROOT = Path(__file__).resolve().parents[2]
EXPECTED_RUFF_PATHS = ("vntyper/", "docker/app/", "tests/", "scripts/", "docs/")
EXPECTED_UNRESOLVED_BEHAVIOR_NODES = {
    "tests/unit/web/test_result_expiry.py::test_cleanup_continues_after_one_delete_error",
    "tests/unit/web/test_tasks.py::test_input_cleanup_logs_one_error_and_attempts_every_owned_path",
    "tests/unit/test_coverage_gate.py::test_coverage_read_failure_returns_none_and_fails_gate",
    "tests/unit/test_golden_cohort_launcher.py::test_resolve_returns_structured_import_failure",
    "tests/unit/test_advntr_output_parsing.py::test_unreadable_advntr_output_logs_and_returns_none",
    "tests/unit/test_cohort_exports.py::test_pseudonym_table_write_failure_is_logged",
    "tests/unit/test_cohort_inputs.py::test_identity_read_failure_uses_directory_fallback",
    "tests/unit/test_cohort_inputs.py::test_cleanup_attempts_all_directories_after_failure",
    "tests/unit/test_cohort_inputs.py::test_bad_archive_is_skipped_and_other_inputs_continue",
    "tests/unit/test_cohort_inputs.py::test_summary_read_failure_returns_three_empty_results",
    "tests/unit/test_cohort_summary_oracle.py::test_image_encoding_failure_returns_empty_string",
    "tests/unit/test_cohort_summary_oracle.py::test_donut_failure_returns_empty_string",
    "tests/unit/test_cohort_summary_oracle.py::test_report_config_failure_returns_empty_mapping",
    "tests/unit/test_flagging.py::test_invalid_regex_is_false_and_observable",
    "tests/unit/test_generate_report.py::test_kestrel_conversion_failure_preserves_both_frames",
    "tests/unit/test_generate_report.py::test_igv_extraction_failure_returns_empty_fragment",
    "tests/unit/test_generate_report.py::test_fastp_failure_returns_empty_mapping",
    "tests/unit/test_generate_report.py::test_pipeline_log_failure_returns_failure_message",
    "tests/unit/test_generate_report.py::test_pipeline_summary_failure_returns_empty_mapping",
    "tests/unit/test_install_references.py::test_executable_probe_error_returns_false",
    "tests/unit/test_screening_summary.py::test_report_config_failure_returns_empty_mapping",
    "tests/unit/test_summary_parsers.py::test_md5_failure_returns_none",
    "tests/unit/test_summary_parsers.py::test_csv_failure_returns_error_comment",
    "tests/unit/test_summary_parsers.py::test_json_failure_returns_error_mapping",
    "tests/unit/test_summary_parsers.py::test_tsv_failure_returns_error_comment",
    "tests/unit/test_summary_record_step.py::test_parser_failure_is_recorded_and_step_is_appended",
    "tests/unit/test_utils.py::test_tool_version_unexpected_failure_returns_unknown",
}


def _ruff_lint_values() -> tuple[list[str], list[str]]:
    """Read only Ruff's selected rules and per-file ignore values."""
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    lint_match = re.search(r"(?ms)^\[tool\.ruff\.lint\]\n(.*?)(?=^\[|\Z)", text)
    ignores_match = re.search(r"(?ms)^\[tool\.ruff\.lint\.per-file-ignores\]\n(.*?)(?=^\[|\Z)", text)
    if lint_match is None or ignores_match is None:
        raise AssertionError("Ruff lint configuration sections are missing")
    select_match = re.search(r"(?ms)^select\s*=\s*(\[.*?\])", lint_match.group(1))
    if select_match is None:
        raise AssertionError("Ruff select list is missing")
    selected = ast.literal_eval(select_match.group(1))
    ignored = [
        rule
        for value in re.findall(r"(?m)^\s*[^#\n=]+\s*=\s*(\[[^\n]*\])", ignores_match.group(1))
        for rule in ast.literal_eval(value)
    ]
    return selected, ignored


def _track(root: Path, *relative_paths: str) -> None:
    """Create a temporary Git index containing the requested ordinary paths."""
    subprocess.run(["git", "init", "--quiet"], cwd=root, check=True)
    for relative_path in relative_paths:
        path = root / relative_path.rstrip("/")
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            path.write_text("tracked\n", encoding="utf-8")
    subprocess.run(["git", "add", "--", *[path.rstrip("/") for path in relative_paths]], cwd=root, check=True)


def _policy_payload(
    *, handlers: list[dict[str, Any]] | None = None, fail_open: list[dict[str, Any]] | None = None
) -> dict[str, Any]:
    """Build a literal valid policy payload for strict-schema mutation tests."""
    handler_rows = handlers
    if handler_rows is None:
        handler_rows = [
            {
                "path": "module.py",
                "symbol": "run",
                "normal_count": 1,
                "all_count": 1,
                "category": "A",
                "rationale": "Terminal boundary.",
            }
        ]
    categories = Counter(row["category"] for row in handler_rows)
    return {
        "reviewed_ruff_version": "ruff 0.16.1",
        "expected_counts": {
            "normal": sum(row["normal_count"] for row in handler_rows),
            "ignore_noqa": sum(row["all_count"] for row in handler_rows),
            "identities": len(handler_rows),
            "categories": {category: categories[category] for category in ("A", "B", "C")},
        },
        "handlers": handler_rows,
        "fail_open": fail_open if fail_open is not None else [],
    }


def _write_policy(tmp_path: Path, payload: dict[str, Any]) -> Path:
    """Write one temporary policy fixture and return its path."""
    path = tmp_path / "policy.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _fail_open_record(behavior_test: str = "tests/unit/test_module.py::test_exact_behavior") -> dict[str, Any]:
    """Return a complete literal fail-open record."""
    return {
        "path": "module.py",
        "symbol": "run",
        "disposition": "preserved-contract",
        "outcome": "Returns False.",
        "rationale": "The caller treats the fallback as unavailable evidence.",
        "behavior_test": behavior_test,
        "observable_via": "error log",
    }


def _category_c_handler() -> dict[str, Any]:
    """Return a complete literal category-C handler record."""
    return {
        "path": "module.py",
        "symbol": "run",
        "normal_count": 1,
        "all_count": 1,
        "category": "C",
        "rationale": "Fallback boundary.",
    }


def test_load_policy_returns_exact_frozen_dataclasses(tmp_path: Path) -> None:
    """Schema loading preserves every reviewed field without coercion or defaults."""
    payload = _policy_payload(handlers=[_category_c_handler()], fail_open=[_fail_open_record()])
    assert load_policy(_write_policy(tmp_path, payload)) == Policy(
        reviewed_ruff_version="ruff 0.16.1",
        expected_normal=1,
        expected_all=1,
        expected_identities=1,
        expected_categories=(0, 0, 1),
        handlers=(HandlerPolicy("module.py", "run", 1, 1, "C", "Fallback boundary."),),
        fail_open=(
            FailOpenPolicy(
                "module.py",
                "run",
                "preserved-contract",
                "Returns False.",
                "The caller treats the fallback as unavailable evidence.",
                "tests/unit/test_module.py::test_exact_behavior",
                "error log",
            ),
        ),
    )


def test_load_policy_rejects_unreadable_or_invalid_json(tmp_path: Path) -> None:
    """An absent or syntactically invalid artifact cannot become an empty policy."""
    with pytest.raises(ValueError, match="read policy"):
        load_policy(tmp_path / "missing.json")
    invalid = tmp_path / "invalid.json"
    invalid.write_text("not-json", encoding="utf-8")
    with pytest.raises(ValueError, match="parse policy JSON"):
        load_policy(invalid)


@pytest.mark.parametrize(
    ("case", "expected"),
    [
        ("missing reviewed version", "reviewed_ruff_version"),
        ("empty policy", "at least one handler"),
        ("zero aggregate count", "strictly positive"),
        ("zero handler count", "all_count.*strictly positive"),
        ("negative count", "non-negative"),
        ("normal greater than all", "normal_count.*all_count"),
        ("invalid category", "category"),
        ("empty rationale", "rationale"),
        ("duplicate handler", "duplicate handler"),
        ("duplicate fail-open", "duplicate fail_open"),
        ("absolute path", "relative"),
        ("parent path", "parent traversal"),
        ("category C without fail-open", "category C.*fail_open"),
        ("fail-open for non-C", "non-category-C"),
        ("invalid disposition", "disposition"),
        ("empty test node", "behavior_test"),
        ("category total mismatch", "handler totals"),
        ("category identity mismatch", "category identity counts"),
        ("unknown top-level key", "unknown keys"),
        ("unknown handler key", "unknown keys"),
    ],
)
def test_load_policy_rejects_malformed_schema(tmp_path: Path, case: str, expected: str) -> None:
    """Malformed or ambiguous reviewed metadata always fails closed."""
    payload = _policy_payload()
    if case == "missing reviewed version":
        del payload["reviewed_ruff_version"]
    elif case == "empty policy":
        payload["expected_counts"] = {
            "normal": 0,
            "ignore_noqa": 0,
            "identities": 0,
            "categories": {"A": 0, "B": 0, "C": 0},
        }
        payload["handlers"] = []
    elif case == "zero aggregate count":
        payload["expected_counts"]["normal"] = 0
        payload["handlers"][0]["normal_count"] = 0
    elif case == "zero handler count":
        payload["handlers"][0]["normal_count"] = 0
        payload["handlers"][0]["all_count"] = 0
    elif case == "negative count":
        payload["expected_counts"]["normal"] = -1
    elif case == "normal greater than all":
        payload["handlers"][0]["normal_count"] = 2
    elif case == "invalid category":
        payload["handlers"][0]["category"] = "D"
    elif case == "empty rationale":
        payload["handlers"][0]["rationale"] = " "
    elif case == "duplicate handler":
        payload["expected_counts"] = {
            "normal": 2,
            "ignore_noqa": 2,
            "identities": 2,
            "categories": {"A": 2, "B": 0, "C": 0},
        }
        payload["handlers"].append(dict(payload["handlers"][0]))
    elif case == "duplicate fail-open":
        payload = _policy_payload(
            handlers=[_category_c_handler()],
            fail_open=[_fail_open_record(), _fail_open_record()],
        )
    elif case == "absolute path":
        payload["handlers"][0]["path"] = "/module.py"
    elif case == "parent path":
        payload["handlers"][0]["path"] = "../module.py"
    elif case == "category C without fail-open":
        payload = _policy_payload(handlers=[_category_c_handler()], fail_open=[])
    elif case == "fail-open for non-C":
        payload["fail_open"] = [_fail_open_record()]
    elif case == "invalid disposition":
        row = _fail_open_record()
        row["disposition"] = "later"
        payload = _policy_payload(handlers=[_category_c_handler()], fail_open=[row])
    elif case == "empty test node":
        payload = _policy_payload(handlers=[_category_c_handler()], fail_open=[_fail_open_record("")])
    elif case == "category total mismatch":
        payload["expected_counts"]["ignore_noqa"] = 2
    elif case == "category identity mismatch":
        payload["handlers"][0]["category"] = "B"
    elif case == "unknown top-level key":
        payload["misspelled"] = []
    elif case == "unknown handler key":
        payload["handlers"][0]["misspelled"] = True
    with pytest.raises(ValueError, match=expected):
        load_policy(_write_policy(tmp_path, payload))


@pytest.mark.parametrize(
    "node_id",
    [
        "tests/unit/test_module.py::test_exact_behavior",
        "tests/unit/web/test_tasks.py::test_input_cleanup",
        "tests/unit/test_module.py::TestPolicy::test_exact_behavior",
        "tests/unit/test_module.py::test_exact_behavior[param-value]",
    ],
)
def test_load_policy_accepts_supported_behavior_node_shapes(tmp_path: Path, node_id: str) -> None:
    """Top-level, nested, class-qualified, and parametrized unit nodes are valid."""
    payload = _policy_payload(handlers=[_category_c_handler()], fail_open=[_fail_open_record(node_id)])
    assert load_policy(_write_policy(tmp_path, payload)).fail_open[0].behavior_test == node_id


@pytest.mark.parametrize(
    "node_id",
    [
        "tests/integration/test_module.py::test_exact_behavior",
        "tests/unit/module.py::test_exact_behavior",
        "tests/unit/../unit/test_module.py::test_exact_behavior",
        "tests/unit/test_module.py::test_exact_behavior\nsecond",
        "tests/unit/test_module.py",
    ],
)
def test_load_policy_rejects_unsupported_behavior_node_shapes(tmp_path: Path, node_id: str) -> None:
    """Behavior evidence cannot escape or ambiguously name the unit tier."""
    payload = _policy_payload(handlers=[_category_c_handler()], fail_open=[_fail_open_record(node_id)])
    with pytest.raises(ValueError, match="behavior_test"):
        load_policy(_write_policy(tmp_path, payload))


def test_behavior_node_resolution_accepts_top_level_class_and_parametrized_nodes(tmp_path: Path) -> None:
    """Static resolution recognizes every supported pytest node shape without importing tests."""
    test_path = tmp_path / "tests/unit/test_module.py"
    test_path.parent.mkdir(parents=True)
    test_path.write_text(
        "def test_top():\n    pass\n\nclass TestPolicy:\n    def test_method(self):\n        pass\n",
        encoding="utf-8",
    )
    _track(tmp_path, "tests/unit/test_module.py")
    assert behavior_node_error(tmp_path, "tests/unit/test_module.py::test_top[param]") is None
    assert behavior_node_error(tmp_path, "tests/unit/test_module.py::TestPolicy::test_method") is None


def test_behavior_node_resolution_reports_unreadable_unparseable_and_missing_paths(tmp_path: Path) -> None:
    """Each unresolved-node failure names whether reading, parsing, or symbol lookup failed."""
    missing = behavior_node_error(tmp_path, "tests/unit/test_missing.py::test_absent")
    assert missing is not None and "could not read tests/unit/test_missing.py" in missing

    broken_path = tmp_path / "tests/unit/test_broken.py"
    broken_path.parent.mkdir(parents=True)
    broken_path.write_text("def broken(:\n", encoding="utf-8")
    _track(tmp_path, "tests/unit/test_broken.py")
    broken = behavior_node_error(tmp_path, "tests/unit/test_broken.py::test_absent")
    assert broken is not None and "could not parse tests/unit/test_broken.py" in broken

    valid_path = tmp_path / "tests/unit/test_valid.py"
    valid_path.write_text("def test_present():\n    pass\n", encoding="utf-8")
    subprocess.run(["git", "add", "--", "tests/unit/test_valid.py"], cwd=tmp_path, check=True)
    unresolved = behavior_node_error(tmp_path, "tests/unit/test_valid.py::TestMissing::test_absent")
    assert unresolved is not None and "could not resolve TestMissing" in unresolved


def test_behavior_node_resolution_rejects_symlink_and_untracked_source_evidence(tmp_path: Path) -> None:
    """A node resolves only from an ordinary Git-tracked source file inside the repository."""
    outside = tmp_path.parent / f"{tmp_path.name}-outside.py"
    outside.write_text("def test_external():\n    pass\n", encoding="utf-8")
    linked = tmp_path / "tests/unit/test_linked.py"
    linked.parent.mkdir(parents=True)
    linked.symlink_to(outside)
    _track(tmp_path, "tests/unit/test_linked.py")
    symlink_error = behavior_node_error(tmp_path, "tests/unit/test_linked.py::test_external")
    assert symlink_error is not None and "symlink" in symlink_error

    untracked = tmp_path / "tests/unit/test_untracked.py"
    untracked.write_text("def test_present():\n    pass\n", encoding="utf-8")
    untracked_error = behavior_node_error(tmp_path, "tests/unit/test_untracked.py::test_present")
    assert untracked_error is not None and "tracked" in untracked_error


def test_policy_rejects_a_category_c_handler_without_behavior_evidence(tmp_path: Path) -> None:
    """Every category-C symbol must have one complete fail-open record."""
    payload = _policy_payload(handlers=[_category_c_handler()], fail_open=[])
    with pytest.raises(ValueError, match="category C.*fail_open"):
        load_policy(_write_policy(tmp_path, payload))


def _synthetic_policy(*, normal_count: int = 1, all_count: int = 1, symbol: str = "run") -> Policy:
    """Return literal category-A synthetic policy data for identity validation tests."""
    return Policy(
        "ruff 0.16.1",
        normal_count,
        all_count,
        1,
        (1, 0, 0),
        (HandlerPolicy("module.py", symbol, normal_count, all_count, "A", "Terminal boundary."),),
        (),
    )


def _write_synthetic_handler(tmp_path: Path, symbol: str = "run") -> Diagnostic:
    """Write one broad handler and return its exact synthetic Ruff diagnostic."""
    (tmp_path / "module.py").write_text(
        f"def {symbol}():\n    try:\n        pass\n    except Exception:\n        raise\n",
        encoding="utf-8",
    )
    return Diagnostic("module.py", 4, 5, "BLE001", "blind-except")


def test_validate_policy_accepts_exact_identity_counts_and_reviewed_version(tmp_path: Path) -> None:
    """Exact reviewed Ruff 0.16.1 diagnostics produce no validation errors."""
    diagnostic = _write_synthetic_handler(tmp_path)
    measurement = Measurement("ruff 0.16.1", (diagnostic,))
    assert validate_policy(tmp_path, measurement, measurement, _synthetic_policy()) == []


def test_validate_policy_accepts_version_only_drift_with_identical_diagnostics(tmp_path: Path) -> None:
    """Reviewed Ruff metadata is attribution, not a dependency pin."""
    diagnostic = _write_synthetic_handler(tmp_path)
    measurement = Measurement("ruff 0.16.0", (diagnostic,))
    assert validate_policy(tmp_path, measurement, measurement, _synthetic_policy()) == []


def test_validate_policy_attributes_identity_drift_to_reviewed_and_actual_versions(tmp_path: Path) -> None:
    """A tool-version change plus identity drift names both measuring contexts."""
    diagnostic = _write_synthetic_handler(tmp_path)
    errors = validate_policy(
        tmp_path,
        Measurement("ruff 0.16.0", (diagnostic,)),
        Measurement("ruff 0.16.0", (diagnostic,)),
        _synthetic_policy(symbol="stale"),
    )
    assert any(
        "reviewed Ruff ruff 0.16.1" in error and "actual Ruff ruff 0.16.0" in error and "reclassify" in error
        for error in errors
    )


def test_validate_policy_rejects_different_versions_between_measurement_modes(tmp_path: Path) -> None:
    """Suppression comparisons cannot combine inventories emitted by different Ruff versions."""
    diagnostic = _write_synthetic_handler(tmp_path)
    errors = validate_policy(
        tmp_path,
        Measurement("ruff 0.16.1", (diagnostic,)),
        Measurement("ruff 0.16.0", (diagnostic,)),
        _synthetic_policy(),
    )
    assert any("normal Ruff ruff 0.16.1" in error and "ignore-noqa Ruff ruff 0.16.0" in error for error in errors)


@pytest.mark.parametrize(
    ("policy", "normal_count", "all_count", "expected"),
    [
        (Policy("ruff 0.16.1", 0, 0, 0, (0, 0, 0), (), ()), 1, 1, "module.py::run.*expected 0.*actual 1"),
        (_synthetic_policy(normal_count=0, all_count=1), 1, 1, "normal count.*expected 0.*actual 1"),
        (_synthetic_policy(symbol="stale"), 1, 1, "removed.*module.py::stale"),
        (_synthetic_policy(), 2, 2, "module.py::run.*expected 1.*actual 2"),
    ],
)
def test_validate_policy_rejects_unclassified_changed_stale_or_multiply_counted_symbols(
    tmp_path: Path,
    policy: Policy,
    normal_count: int,
    all_count: int,
    expected: str,
) -> None:
    """Stable identity counts must match exactly rather than merely preserve aggregate totals."""
    diagnostic = _write_synthetic_handler(tmp_path)
    normal = Measurement("ruff 0.16.1", (diagnostic,) * normal_count)
    all_handlers = Measurement("ruff 0.16.1", (diagnostic,) * all_count)
    assert any(re.search(expected, error) for error in validate_policy(tmp_path, normal, all_handlers, policy))


def test_validate_policy_rejects_a_normal_diagnostic_absent_from_all_mode(tmp_path: Path) -> None:
    """Normal diagnostics must be an exact multiset subset of ignore-noqa diagnostics."""
    diagnostic = _write_synthetic_handler(tmp_path)
    errors = validate_policy(
        tmp_path,
        Measurement("ruff 0.16.1", (diagnostic,)),
        Measurement("ruff 0.16.1", ()),
        _synthetic_policy(normal_count=1, all_count=0),
    )
    assert any("normal diagnostics are not a subset" in error for error in errors)


def test_validate_policy_reports_unreadable_or_unparseable_source(tmp_path: Path) -> None:
    """An unreadable source names its path and cannot produce a guessed symbol."""
    diagnostic = Diagnostic("missing.py", 1, 1, "BLE001", "blind-except")
    errors = validate_policy(
        tmp_path,
        Measurement("ruff 0.16.1", (diagnostic,)),
        Measurement("ruff 0.16.1", (diagnostic,)),
        Policy("ruff 0.16.1", 1, 1, 0, (0, 0, 0), (), ()),
    )
    assert any("missing.py" in error and "read" in error for error in errors)


def test_validate_policy_reports_an_unparseable_source_path(tmp_path: Path) -> None:
    """A syntactically invalid diagnostic source reaches and reports AST parsing."""
    (tmp_path / "broken.py").write_text("def broken(:\n", encoding="utf-8")
    diagnostic = Diagnostic("broken.py", 1, 1, "BLE001", "blind-except")
    errors = validate_policy(
        tmp_path,
        Measurement("ruff 0.16.1", (diagnostic,)),
        Measurement("ruff 0.16.1", (diagnostic,)),
        Policy("ruff 0.16.1", 1, 1, 0, (0, 0, 0), (), ()),
    )
    assert any("broken.py" in error and "parse" in error for error in errors)


def test_live_ble001_diagnostics_match_reviewed_handler_counts() -> None:
    """Every live BLE001 diagnostic maps exactly once to the reviewed inventory."""
    scope = read_ruff_paths(REPO_ROOT / "Makefile")
    assert scope == EXPECTED_RUFF_PATHS
    normal = measure_ble001(REPO_ROOT, scope, ignore_noqa=False)
    all_handlers = measure_ble001(REPO_ROOT, scope, ignore_noqa=True)
    policy = load_policy(REPO_ROOT / "scripts/ble001_policy.json")

    assert policy.expected_normal == 103
    assert policy.expected_all == 108
    assert policy.expected_all - policy.expected_normal == 5
    assert policy.expected_identities == 79
    assert policy.expected_categories == (30, 16, 33)
    assert len(policy.handlers) == 79
    assert Counter(row.category for row in policy.handlers) == Counter({"A": 30, "B": 16, "C": 33})
    assert sum(row.normal_count for row in policy.handlers) == 103
    assert sum(row.all_count for row in policy.handlers) == 108
    assert len(normal.diagnostics) == policy.expected_normal
    assert len(all_handlers.diagnostics) == policy.expected_all
    expected_normal = {(row.path, row.symbol): row.normal_count for row in policy.handlers}
    expected_all = {(row.path, row.symbol): row.all_count for row in policy.handlers}

    def counts(measurement: Measurement) -> dict[tuple[str, str], int]:
        found: Counter[tuple[str, str]] = Counter()
        for diagnostic in measurement.diagnostics:
            source = (REPO_ROOT / diagnostic.path).read_text(encoding="utf-8")
            found[(diagnostic.path, enclosing_symbol(source, diagnostic.row))] += 1
        return dict(found)

    assert counts(normal) == {key: value for key, value in expected_normal.items() if value}
    assert counts(all_handlers) == {key: value for key, value in expected_all.items() if value}
    assert {(row.path, row.symbol) for row in policy.handlers if row.category == "C"} == {
        (row.path, row.symbol) for row in policy.fail_open
    }


def test_live_policy_cli_lists_unresolved_behavior_nodes_without_success() -> None:
    """Tasks 0-3 remain visibly incomplete until every frozen behavior node exists."""
    result = subprocess.run(
        [sys.executable, "scripts/ble001_policy.py", "--repo-root", ".", "--policy", "scripts/ble001_policy.json"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    actual = set(re.findall(r"(?m)^unresolved behavior-test node ([^:]+::[^:]+):", result.stdout))
    assert actual == EXPECTED_UNRESOLVED_BEHAVIOR_NODES
    assert result.stdout.count("unresolved behavior-test node") == len(EXPECTED_UNRESOLVED_BEHAVIOR_NODES)
    assert "categories A/B/C:" not in result.stdout


def test_ble001_and_g004_remain_deliberately_unselected() -> None:
    """The frozen inventory must not silently become an enabled or ignored lint rule."""
    selected, ignored = _ruff_lint_values()
    assert "BLE001" not in selected
    assert "BLE001" not in ignored
    assert "G004" not in selected
    assert "G004" not in ignored
