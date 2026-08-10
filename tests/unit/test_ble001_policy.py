"""Unit tests for the reviewed BLE001 exception-policy adapter."""

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import ble001_policy as policy_module  # noqa: E402
from ble001_policy import (  # noqa: E402
    Diagnostic,
    HandlerPolicy,
    Measurement,
    Policy,
    enclosing_symbol,
    main,
    measure_ble001,
    normalize_diagnostics,
    read_ruff_paths,
)

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


def test_read_ruff_paths_reads_the_single_make_authority(tmp_path: Path) -> None:
    """A valid single assignment returns the exact ordered path scope."""
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


@pytest.mark.parametrize("payload", ["not-json", "{}", '[{"filename": "../outside.py"}]'])
def test_normalize_diagnostics_rejects_malformed_or_out_of_root_payloads(tmp_path: Path, payload: str) -> None:
    """Malformed structures and paths escaping the repository fail closed."""
    with pytest.raises(ValueError):
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
        [*base, *scope],
        ["ruff", "--version"],
        [*base, "--ignore-noqa", *scope],
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
