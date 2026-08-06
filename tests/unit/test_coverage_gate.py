"""Unit tests for the coverage ratchet gate."""

import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import coverage_gate  # noqa: E402


def test_read_floor_returns_the_configured_value() -> None:
    """The floor is parsed from the real pyproject.toml."""
    assert coverage_gate.read_floor() >= 25.0


def test_read_floor_raises_when_pyproject_is_unreadable(monkeypatch, tmp_path) -> None:
    """A gate that cannot find its threshold must fail, not silently pass."""
    monkeypatch.setattr(coverage_gate, "PYPROJECT", tmp_path / "nope.toml")
    with pytest.raises(SystemExit):
        coverage_gate.read_floor()


def test_read_floor_raises_when_fail_under_is_absent(monkeypatch, tmp_path) -> None:
    """A pyproject without fail_under must fail loudly."""
    stub = tmp_path / "pyproject.toml"
    stub.write_text("[tool.coverage.report]\nshow_missing = true\n", encoding="utf-8")
    monkeypatch.setattr(coverage_gate, "PYPROJECT", stub)
    with pytest.raises(SystemExit):
        coverage_gate.read_floor()


def _run_gate(monkeypatch, floor: float, total: float) -> str:
    """Run ``main()`` against a stubbed floor and coverage total.

    Args:
        monkeypatch: Pytest monkeypatch fixture.
        floor: The floor ``read_floor`` should report.
        total: The coverage ``read_total`` should report.

    Returns:
        str: Everything the gate printed.
    """
    monkeypatch.setattr(coverage_gate, "read_floor", lambda: floor)
    monkeypatch.setattr(coverage_gate, "read_total", lambda: total)
    monkeypatch.setattr(sys, "argv", ["coverage_gate.py"])
    printed: list[str] = []
    monkeypatch.setattr(coverage_gate, "emit", lambda line, **_: printed.append(line))
    monkeypatch.setattr("builtins.print", lambda *args, **_: printed.append(" ".join(str(a) for a in args)))
    coverage_gate.main()
    return "\n".join(printed)


def test_the_ratchet_stays_quiet_when_the_floor_is_already_the_rounded_coverage(monkeypatch) -> None:
    """36.81% against a floor of 36 has no integer left to gain.

    Suggesting `fail_under = 36` when it is already 36 is noise that trains
    readers to ignore the notice.
    """
    output = _run_gate(monkeypatch, floor=36.0, total=36.81)
    assert "Raise it" not in output


def test_the_ratchet_names_the_new_floor_once_coverage_crosses_an_integer(monkeypatch) -> None:
    """37.20% against a floor of 36 is a real, lockable gain."""
    output = _run_gate(monkeypatch, floor=36.0, total=37.20)
    assert "set `fail_under = 37`" in output


def test_the_gate_fails_below_the_floor(monkeypatch) -> None:
    """Coverage under the floor is an error exit, not a warning."""
    monkeypatch.setattr(coverage_gate, "read_floor", lambda: 36.0)
    monkeypatch.setattr(coverage_gate, "read_total", lambda: 35.0)
    monkeypatch.setattr(sys, "argv", ["coverage_gate.py"])
    assert coverage_gate.main() == 1
