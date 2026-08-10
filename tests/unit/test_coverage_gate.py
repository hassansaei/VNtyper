"""Unit tests for the coverage ratchet gate."""

import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import coverage_gate  # noqa: E402

# The floor is a ratchet: it only ever moves up. Pinning the literal here is the point -
# a `>=` bound would let someone lower `fail_under` in pyproject.toml all the way to the
# bound and still leave this test green, which is exactly the regression the gate exists
# to prevent. When you legitimately raise the floor (use the number `make test-unit-cov`
# prints, never the rounded TOTAL column), raise this literal in the same commit.
#
# 86 is a BRANCH-INCLUSIVE figure: `branch = true` was enabled in #196, taking the floor
# from 70 to the 74.22% that run measured; splitting cohort_summary.py into five focused
# modules took the same suite to 79.02%; and the fixes and tests of the gated phase took
# it to 80.24%. Statement-only coverage is higher still, so the two measurements are not
# interchangeable. Milestone 4's focused modules and failure-path tests first measured
# 84.68% and moved the integer floor to 84; the completed input-safety closure measured
# 86.23% and moved it to 86. See `test_branch_coverage_is_enabled` below for why the gap
# has to be guarded.
CURRENT_COVERAGE_FLOOR = 86.0


def test_read_floor_returns_the_configured_value() -> None:
    """The floor parsed from the real pyproject.toml is exactly the pinned value.

    Raises:
        AssertionError: If `fail_under` drifts from ``CURRENT_COVERAGE_FLOOR`` in either
            direction - lowered (a silent weakening of CI) or raised without updating
            this pin.
    """
    floor = coverage_gate.read_floor()
    assert floor == CURRENT_COVERAGE_FLOOR, (
        f"pyproject.toml [tool.coverage.report] fail_under is {floor:.0f}, but this test "
        f"pins it at {CURRENT_COVERAGE_FLOOR:.0f}. The floor is a ratchet: never lower it "
        "to make a build pass - add tests instead. If you raised it, update "
        "CURRENT_COVERAGE_FLOOR here in the same commit."
    )
    assert floor >= CURRENT_COVERAGE_FLOOR, "The coverage floor must never be lowered."


def test_branch_coverage_is_enabled() -> None:
    """`[tool.coverage.run] branch = true` must stay on, because nothing else notices.

    A config assertion earns a test here for one specific reason: **the ratchet cannot
    catch this regression, because deleting the setting moves the number the WRONG WAY.**
    Both figures below are measured, on the same suite and the same commit:

        branch-inclusive (branch = true)   74.22%   <- what fail_under = 74 was set from
        statement-only   (setting removed) 76.60%   <- clears the floor by 2.60 points

    So removing `branch = true` *raises* the reported total, the coverage gate goes
    green, CI goes green, and the measurement has silently been weakened back to the one
    #196 replaced - an `if` that is entered but never taken counted as fully covered.
    A floor can only catch coverage falling; this makes it rise. That is why the
    assertion has to exist, and it is the only thing standing between that edit and a
    passing build.

    The config is read through coverage.py's own parser rather than by grepping the
    TOML, so a commented-out or overridden setting fails here exactly as it would in
    a real run.

    Raises:
        AssertionError: If branch coverage is off in pyproject.toml.
    """
    import coverage

    config = coverage.Coverage(config_file=str(coverage_gate.PYPROJECT)).config

    assert config.branch is True, (
        "pyproject.toml [tool.coverage.run] no longer sets `branch = true`. The floor "
        f"of {CURRENT_COVERAGE_FLOOR:.0f} is a branch-inclusive figure; measuring "
        "statements only would clear it while covering strictly less. Re-enable branch "
        "coverage rather than relying on the floor to catch this - it cannot."
    )


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


@pytest.mark.parametrize("corrupt", [False, True])
def test_read_total_returns_none_for_missing_or_corrupt_coverage_data(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, corrupt: bool
) -> None:
    monkeypatch.chdir(tmp_path)
    if corrupt:
        (tmp_path / ".coverage").write_bytes(b"not a coverage database")
    assert coverage_gate.read_total() is None


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


def test_unreadable_coverage_fails_the_gate(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(sys, "argv", ["coverage_gate.py"])
    monkeypatch.setattr(coverage_gate, "read_floor", lambda: 86.0)
    monkeypatch.setattr(coverage_gate, "read_total", lambda: None)
    assert coverage_gate.main() == 1
    assert "could not read" in capsys.readouterr().out.lower()


def test_total_below_the_soft_target_is_warning_only(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(sys, "argv", ["coverage_gate.py", "--target", "90", "--ratchet-margin", "100"])
    monkeypatch.setattr(coverage_gate, "read_floor", lambda: 80.0)
    monkeypatch.setattr(coverage_gate, "read_total", lambda: 85.0)
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)

    assert coverage_gate.main() == 0
    output = capsys.readouterr().out
    assert "Coverage 85.00% is below the 90% target" in output
    assert "Not a failure" in output


def test_write_summary_is_a_noop_without_a_github_summary_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    summary = tmp_path / "summary.md"
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    coverage_gate.write_summary(87.5, 86.0, 90.0)
    assert not summary.exists()


def test_write_summary_appends_exact_github_markdown(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    summary = tmp_path / "summary.md"
    summary.write_text("existing\n", encoding="utf-8")
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))

    coverage_gate.write_summary(87.5, 86.0, 90.0)

    assert summary.read_text(encoding="utf-8") == (
        "existing\n"
        "### Unit test coverage\n\n"
        "`███████████████████████████████████░░░░░` **87.50%**\n\n"
        "| | % |\n| --- | --- |\n"
        "| Measured | **87.50** |\n"
        "| Hard floor (CI fails below) | 86 |\n"
        "| Target | 90 — 🎯 below target |\n\n"
    )


def test_scripts_only_coverage_scope_matches_the_repository_tree() -> None:
    scripts = sorted(
        path.relative_to(coverage_gate.REPO_ROOT).as_posix() for path in coverage_gate.REPO_ROOT.glob("scripts/**/*.py")
    )
    makefile = (coverage_gate.REPO_ROOT / "Makefile").read_text(encoding="utf-8")
    assert scripts
    assert "--cov=scripts" in makefile
    assert "--cov-omit" not in makefile
    assert all(path.startswith("scripts/") for path in scripts)
