"""Unit tests for the reviewed BLE001 exception-policy adapter."""

import ast
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import ble001_policy as policy_module  # noqa: E402
from ble001_policy import (  # noqa: E402
    Diagnostic,
    FailOpenPolicy,
    HandlerPolicy,
    Measurement,
    Policy,
    enclosing_symbol,
    load_policy,
    main,
    measure_ble001,
    normalize_diagnostics,
    read_ruff_paths,
    validate_policy,
)
from ble001_policy_validation import behavior_node_error  # noqa: E402

pytestmark = pytest.mark.unit
REPO_ROOT = Path(__file__).resolve().parents[2]

SOURCE = """\
try:
    pass
except Exception:
    pass

def outer():
    def inner():
        try:
            pass
        except Exception:
            return None
    return inner()

class Worker:
    def run(self):
        try:
            pass
        except Exception:
            return 1
"""


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


def test_read_ruff_paths_reads_the_single_make_authority(tmp_path: Path) -> None:
    """A valid single assignment returns the exact ordered path scope."""
    for relative_path in ("vntyper", "docker/app", "tests", "scripts", "docs"):
        (tmp_path / relative_path).mkdir(parents=True)
    makefile = tmp_path / "Makefile"
    makefile.write_text("RUFF_PATHS := vntyper/ docker/app/ tests/ scripts/ docs/\n", encoding="utf-8")
    assert read_ruff_paths(makefile) == ("vntyper/", "docker/app/", "tests/", "scripts/", "docs/")


def test_read_ruff_paths_rejects_missing_or_duplicate_assignments(tmp_path: Path) -> None:
    """A missing or ambiguous path authority fails closed."""
    makefile = tmp_path / "Makefile"
    makefile.write_text("all:\n\ttrue\n", encoding="utf-8")
    with pytest.raises(ValueError, match="exactly one RUFF_PATHS"):
        read_ruff_paths(makefile)
    makefile.write_text("RUFF_PATHS := a/\nRUFF_PATHS := b/\n", encoding="utf-8")
    with pytest.raises(ValueError, match="exactly one RUFF_PATHS"):
        read_ruff_paths(makefile)


@pytest.mark.parametrize(
    "scope",
    [
        "--output-file captured.json module.py",
        "$(shell-touch-marker)",
        "missing.py",
        "../outside.py",
        "/absolute.py",
        "folder\\file.py",
    ],
)
def test_read_ruff_paths_rejects_options_nonpaths_and_repository_escapes(tmp_path: Path, scope: str) -> None:
    """Only existing ordinary paths contained by the repository can enter Ruff argv."""
    makefile = tmp_path / "Makefile"
    makefile.write_text(f"RUFF_PATHS := {scope}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="RUFF_PATHS token"):
        read_ruff_paths(makefile)


def test_read_ruff_paths_rejects_a_symlink_escape(tmp_path: Path) -> None:
    """A pathname spelled below the repository cannot resolve to an outside scope."""
    outside = tmp_path.parent / f"{tmp_path.name}-outside.py"
    outside.write_text("", encoding="utf-8")
    (tmp_path / "escaped.py").symlink_to(outside)
    makefile = tmp_path / "Makefile"
    makefile.write_text("RUFF_PATHS := escaped.py\n", encoding="utf-8")
    with pytest.raises(ValueError, match="escapes repository"):
        read_ruff_paths(makefile)


def test_normalize_diagnostics_returns_sorted_relative_records(tmp_path: Path) -> None:
    """Ruff output order and absolute in-root paths do not affect normalized identity."""
    second = tmp_path / "z.py"
    first = tmp_path / "a.py"
    payload = json.dumps(
        [
            {
                "filename": str(second),
                "code": "BLE001",
                "message": "second",
                "location": {"row": 9, "column": 3},
            },
            {
                "filename": str(first),
                "code": "BLE001",
                "message": "first",
                "location": {"row": 2, "column": 1},
            },
        ]
    )
    assert normalize_diagnostics(payload, tmp_path) == (
        Diagnostic("a.py", 2, 1, "BLE001", "first"),
        Diagnostic("z.py", 9, 3, "BLE001", "second"),
    )


def test_normalize_diagnostics_accepts_the_cell_start_location_spelling(tmp_path: Path) -> None:
    """The Ruff compatibility location spelling produces the same normalized record."""
    payload = json.dumps(
        [
            {
                "filename": "module.py",
                "code": "BLE001",
                "message": "blind-except",
                "location": {"cell": {"start": {"row": 4, "column": 7}}},
            }
        ]
    )
    assert normalize_diagnostics(payload, tmp_path) == (Diagnostic("module.py", 4, 7, "BLE001", "blind-except"),)


@pytest.mark.parametrize("payload", ["not-json", "{}"])
def test_normalize_diagnostics_rejects_malformed_payloads(tmp_path: Path, payload: str) -> None:
    """Malformed top-level Ruff output fails closed."""
    with pytest.raises(ValueError):
        normalize_diagnostics(payload, tmp_path)


def test_normalize_diagnostics_reaches_repository_containment_validation(tmp_path: Path) -> None:
    """A complete diagnostic with an escaping filename fails at the containment boundary."""
    payload = json.dumps(
        [
            {
                "filename": "../outside.py",
                "code": "BLE001",
                "message": "blind-except",
                "location": {"row": 1, "column": 1},
            }
        ]
    )
    with pytest.raises(ValueError, match="escapes repository root"):
        normalize_diagnostics(payload, tmp_path)


def test_normalize_diagnostics_reaches_malformed_location_validation(tmp_path: Path) -> None:
    """A complete diagnostic with a malformed location names that exact defect."""
    payload = json.dumps(
        [
            {
                "filename": "module.py",
                "code": "BLE001",
                "message": "blind-except",
                "location": {"row": 0, "column": 1},
            }
        ]
    )
    with pytest.raises(ValueError, match="invalid row or column"):
        normalize_diagnostics(payload, tmp_path)


def test_enclosing_symbol_uses_qualified_names_not_lines() -> None:
    """Blank-line movement preserves a handler's qualified symbol identity."""
    assert enclosing_symbol(SOURCE, 3) == "<module>"
    assert enclosing_symbol(SOURCE, 10) == "outer.inner"
    assert enclosing_symbol(SOURCE, 18) == "Worker.run"
    assert enclosing_symbol("\n\n" + SOURCE, 12) == "outer.inner"


def test_enclosing_symbol_rejects_invalid_source_or_row() -> None:
    """Invalid AST input and impossible row numbers cannot produce an identity."""
    with pytest.raises(ValueError, match="parse"):
        enclosing_symbol("def broken(:\n", 1)
    with pytest.raises(ValueError, match="positive"):
        enclosing_symbol(SOURCE, 0)


def test_measure_ble001_runs_independent_versioned_normal_and_all_measurements(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Dropping no-cache, version attribution, scope, or all-mode isolation breaks the command contract."""
    scope = ("vntyper/", "docker/app/", "tests/", "scripts/", "docs/")
    normal_payload = json.dumps(
        [
            {
                "filename": "vntyper/module.py",
                "code": "BLE001",
                "message": "blind-except",
                "location": {"row": 3, "column": 1},
            }
        ]
    )
    responses = iter(
        [
            subprocess.CompletedProcess(["ruff", "--version"], 0, "ruff 0.16.1\n", ""),
            subprocess.CompletedProcess(["ruff", "check"], 1, normal_payload, ""),
            subprocess.CompletedProcess(["ruff", "--version"], 0, "ruff 0.16.1\n", ""),
            subprocess.CompletedProcess(["ruff", "check"], 0, "[]", ""),
        ]
    )
    calls: list[tuple[list[str], dict[str, Any]]] = []

    def run(argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append((argv, kwargs))
        return next(responses)

    monkeypatch.setattr(policy_module.subprocess, "run", run)

    normal = measure_ble001(tmp_path, scope, ignore_noqa=False)
    all_handlers = measure_ble001(tmp_path, scope, ignore_noqa=True)

    assert normal == Measurement(
        "ruff 0.16.1",
        (Diagnostic("vntyper/module.py", 3, 1, "BLE001", "blind-except"),),
    )
    assert all_handlers == Measurement("ruff 0.16.1", ())
    base = ["ruff", "check", "--no-cache", "--select", "BLE001", "--output-format", "json"]
    assert [call[0] for call in calls] == [
        ["ruff", "--version"],
        [*base, "--", *scope],
        ["ruff", "--version"],
        [*base, "--ignore-noqa", "--", *scope],
    ]
    assert all(
        kwargs
        == {
            "cwd": tmp_path,
            "capture_output": True,
            "text": True,
            "check": False,
        }
        for _, kwargs in calls
    )


def test_measure_ble001_reports_ruff_check_failure_with_actual_version(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A Ruff operational failure cannot be mistaken for an empty inventory."""
    responses = iter(
        [
            subprocess.CompletedProcess(["custom-ruff", "--version"], 0, "ruff 0.16.1\n", ""),
            subprocess.CompletedProcess(["custom-ruff", "check"], 2, "", "configuration broke"),
        ]
    )
    monkeypatch.setattr(policy_module.subprocess, "run", lambda *_args, **_kwargs: next(responses))
    with pytest.raises(RuntimeError, match=r"ruff 0\.16\.1.*configuration broke"):
        measure_ble001(tmp_path, ("vntyper/",), ignore_noqa=False, ruff_executable="custom-ruff")


def test_measure_ble001_option_shaped_scope_cannot_write_output(tmp_path: Path) -> None:
    """The argv separator makes a malicious output option an inert path operand."""
    (tmp_path / "module.py").write_text("try:\n    pass\nexcept Exception:\n    pass\n", encoding="utf-8")
    captured = tmp_path / "captured.json"
    measurement = measure_ble001(
        tmp_path,
        ("--output-file", captured.name, "module.py"),
        ignore_noqa=False,
    )
    assert not captured.exists()
    assert len(measurement.diagnostics) == 1


@pytest.mark.parametrize(
    ("version_result", "expected"),
    [
        (subprocess.CompletedProcess(["ruff", "--version"], 1, "", "cannot execute"), "cannot execute"),
        (subprocess.CompletedProcess(["ruff", "--version"], 0, "", ""), "empty version"),
    ],
)
def test_measure_ble001_rejects_unusable_version_attribution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    version_result: subprocess.CompletedProcess[str],
    expected: str,
) -> None:
    """Missing Ruff version evidence fails before any diagnostic command runs."""
    calls = 0

    def run(*_args: Any, **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        return version_result

    monkeypatch.setattr(policy_module.subprocess, "run", run)
    with pytest.raises(RuntimeError, match=expected):
        measure_ble001(tmp_path, ("vntyper/",), ignore_noqa=False)
    assert calls == 1


def test_measure_ble001_reports_an_unavailable_ruff_executable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An unavailable executable fails closed with its configured name."""

    def run(*_args: Any, **_kwargs: Any) -> subprocess.CompletedProcess[str]:
        raise FileNotFoundError("not found")

    monkeypatch.setattr(policy_module.subprocess, "run", run)
    with pytest.raises(RuntimeError, match="custom-ruff.*not found"):
        measure_ble001(tmp_path, ("vntyper/",), ignore_noqa=False, ruff_executable="custom-ruff")


def _cli_policy() -> Policy:
    """Return a literal two-handler policy for CLI output tests."""
    return Policy(
        reviewed_ruff_version="ruff 0.16.1",
        expected_normal=1,
        expected_all=2,
        expected_identities=2,
        expected_categories=(1, 1, 0),
        handlers=(
            HandlerPolicy("module.py", "run", 1, 1, "A", "Terminal boundary."),
            HandlerPolicy("module.py", "candidate", 0, 1, "B", "Audit candidate."),
        ),
        fail_open=(),
    )


def test_main_returns_one_and_prints_every_validation_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """A validation failure cannot be hidden behind a successful CLI status."""
    policy_path = tmp_path / "policy.json"
    monkeypatch.setattr(policy_module, "read_ruff_paths", lambda _path: ("module.py",))
    monkeypatch.setattr(policy_module, "load_policy", lambda _path: _cli_policy())
    measurements = iter((Measurement("ruff 0.16.1", ()), Measurement("ruff 0.16.1", ())))
    monkeypatch.setattr(policy_module, "measure_ble001", lambda *_args, **_kwargs: next(measurements))
    monkeypatch.setattr(policy_module, "validate_policy", lambda *_args: ["first mismatch", "second mismatch"])

    assert main(["--repo-root", str(tmp_path), "--policy", str(policy_path), "--ruff", "custom-ruff"]) == 1
    output = capsys.readouterr().out
    assert "first mismatch" in output
    assert "second mismatch" in output


def test_main_prints_version_counts_and_categories_for_a_clean_policy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """A clean CLI run reports enough exact evidence to reproduce its result."""
    policy_path = tmp_path / "policy.json"
    monkeypatch.setattr(policy_module, "read_ruff_paths", lambda _path: ("module.py",))
    monkeypatch.setattr(policy_module, "load_policy", lambda _path: _cli_policy())
    measurements = iter((Measurement("ruff 0.16.1", ()), Measurement("ruff 0.16.1", ())))
    monkeypatch.setattr(policy_module, "measure_ble001", lambda *_args, **_kwargs: next(measurements))
    monkeypatch.setattr(policy_module, "validate_policy", lambda *_args: [])

    assert main(["--repo-root", str(tmp_path), "--policy", str(policy_path)]) == 0
    output = capsys.readouterr().out
    assert "reviewed Ruff: ruff 0.16.1" in output
    assert "actual Ruff: ruff 0.16.1" in output
    assert "normal/all: 1/2" in output
    assert "categories A/B/C: 1/1/0" in output


def test_main_fails_closed_on_an_operational_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """An unreadable authority yields status 1 and an actionable error."""
    monkeypatch.setattr(policy_module, "read_ruff_paths", lambda _path: (_ for _ in ()).throw(ValueError("bad scope")))
    assert main(["--repo-root", str(tmp_path), "--policy", str(tmp_path / "policy.json")]) == 1
    assert "bad scope" in capsys.readouterr().out


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
    assert behavior_node_error(tmp_path, "tests/unit/test_module.py::test_top[param]") is None
    assert behavior_node_error(tmp_path, "tests/unit/test_module.py::TestPolicy::test_method") is None


def test_behavior_node_resolution_reports_unreadable_unparseable_and_missing_paths(tmp_path: Path) -> None:
    """Each unresolved-node failure names whether reading, parsing, or symbol lookup failed."""
    missing = behavior_node_error(tmp_path, "tests/unit/test_missing.py::test_absent")
    assert missing is not None and "could not read tests/unit/test_missing.py" in missing

    broken_path = tmp_path / "tests/unit/test_broken.py"
    broken_path.parent.mkdir(parents=True)
    broken_path.write_text("def broken(:\n", encoding="utf-8")
    broken = behavior_node_error(tmp_path, "tests/unit/test_broken.py::test_absent")
    assert broken is not None and "could not parse tests/unit/test_broken.py" in broken

    valid_path = tmp_path / "tests/unit/test_valid.py"
    valid_path.write_text("def test_present():\n    pass\n", encoding="utf-8")
    unresolved = behavior_node_error(tmp_path, "tests/unit/test_valid.py::TestMissing::test_absent")
    assert unresolved is not None and "could not resolve TestMissing" in unresolved


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
    policy = load_policy(REPO_ROOT / "scripts/ble001_policy.json")
    unresolved = [
        record.behavior_test
        for record in policy.fail_open
        if behavior_node_error(REPO_ROOT, record.behavior_test) is not None
    ]
    assert unresolved
    result = subprocess.run(
        [sys.executable, "scripts/ble001_policy.py", "--repo-root", ".", "--policy", "scripts/ble001_policy.json"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    assert result.stdout.count("unresolved behavior-test node") == len(unresolved)
    assert all(f"unresolved behavior-test node {node}:" in result.stdout for node in unresolved)
    assert "categories A/B/C:" not in result.stdout


def test_ble001_and_g004_remain_deliberately_unselected() -> None:
    """The frozen inventory must not silently become an enabled or ignored lint rule."""
    selected, ignored = _ruff_lint_values()
    assert "BLE001" not in selected
    assert "BLE001" not in ignored
    assert "G004" not in selected
    assert "G004" not in ignored
