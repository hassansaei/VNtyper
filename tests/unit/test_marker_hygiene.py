"""
tests/unit/test_marker_hygiene.py

Guards the test suite against silently losing coverage in CI.

CI runs `pytest -m unit`. A test file under tests/unit/ that forgets
`pytestmark = pytest.mark.unit` is collected, deselected, and never run - with no
warning and a green build. That is not hypothetical: test_haplo_count_and_selection.py
sat unmarked with 30 passing tests covering the Issue #136/#145 genotype tie-breaking
logic, invisible to CI.

`--strict-markers` (pytest.ini) catches a *misspelled* marker but cannot catch a
*missing* one. This test closes that gap.

It works off the live session's collected items, so it costs no subprocess and no
measurable time: any file that contributes zero items to a `-m unit` run is either
unmarked or has no tests, and we distinguish those by source inspection.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

UNIT_DIR = Path(__file__).resolve().parent
REPO_ROOT = UNIT_DIR.parents[1]

# A module defines tests if it has a top-level `def test_*` or a `class Test*`.
_DEFINES_TESTS = re.compile(r"^\s*(?:async\s+def\s+test_|def\s+test_|class\s+Test)", re.MULTILINE)

# The `browser` marker, declared either module-wide or on an individual test.
_DECLARES_BROWSER_MARKER = re.compile(
    r"^\s*(?:pytestmark\s*=.*pytest\.mark\.browser|@pytest\.mark\.browser)",
    re.MULTILINE,
)

# The `golden` marker, declared either module-wide or on an individual test.
_DECLARES_GOLDEN_MARKER = re.compile(
    r"^\s*(?:pytestmark\s*=.*pytest\.mark\.golden|@pytest\.mark\.golden)",
    re.MULTILINE,
)


def test_every_unit_file_is_selected_by_the_unit_marker(request: pytest.FixtureRequest) -> None:
    """Every test file under tests/unit/ must be selected by `pytest -m unit`.

    Args:
        request: Pytest request, used to reach the current session's collected items.

    Raises:
        AssertionError: If a file defining tests contributed no items to this run.
    """
    selected_files = {
        Path(str(item.fspath)).resolve() for item in request.session.items if str(item.fspath).endswith(".py")
    }

    unmarked = []
    for path in sorted(UNIT_DIR.rglob("test_*.py")):
        if path.resolve() in selected_files:
            continue
        if _DEFINES_TESTS.search(path.read_text(encoding="utf-8")):
            unmarked.append(str(path.relative_to(UNIT_DIR)))

    assert not unmarked, (
        "These files under tests/unit/ define tests but contributed nothing to a "
        f"`pytest -m unit` run, so CI never executes them: {unmarked}. "
        "Add `pytestmark = pytest.mark.unit` immediately after the imports."
    )


def test_no_test_modules_outside_the_known_tiers() -> None:
    """Test modules must live under a tier some command actually runs.

    The tiers are ``tests/unit`` (pure, what CI gates on), ``tests/integration``
    (needs the Zenodo archive), ``tests/docker`` (needs a Docker daemon),
    ``tests/browser`` (needs a real browser engine; ``make test-browser``, #242),
    and ``tests/golden`` (needs the simulated known-truth cohort).
    ``tests/browser`` is the newest and is deliberately listed here rather than
    tolerated: the report is a file people open in a browser, and its DataTables
    filtering and client-side rounding do not exist until JavaScript has run, so
    they are invisible to every other tier.

    A ``test_*.py`` outside all five is either a mis-placed test that no tier
    runs, or a helper module whose name makes a real collection failure invisible.

    ``golden`` is a tier like the others, run by
    ``pytest -m golden tests/golden``. It is deliberately outside ``make
    check-all``: it compares against a known-truth simulated cohort supplied via
    ``VNTYPER_SIM_ROOT`` and ``VNTYPER_ADVNTR_ROOT``. Selected golden execution
    fails closed without complete roots; deselected collection remains safe.

    Raises:
        AssertionError: If a stray test module is found.
    """
    tests_root = UNIT_DIR.parent
    allowed = {"unit", "integration", "docker", "browser", "golden"}
    stray = [
        str(path.relative_to(tests_root))
        for path in tests_root.rglob("test_*.py")
        if path.relative_to(tests_root).parts[0] not in allowed
    ]
    assert not stray, (
        f"These test modules live outside tests/{{{','.join(sorted(allowed))}}}: {stray}. "
        "Move real tests into a tier; rename helpers so they do not start with 'test_'."
    )


@pytest.mark.parametrize("module_name", ["test_molecular_identity_golden.py", "test_nomenclature_golden.py"])
def test_golden_modules_collect_when_deselected_without_benchmark_roots(module_name: str) -> None:
    """Non-golden collection must not resolve out-of-band roots during import."""
    result = _run_golden_probe(module_name, "not golden", collect_only=True)

    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("module_name", ["test_molecular_identity_golden.py", "test_nomenclature_golden.py"])
def test_selected_golden_modules_fail_without_benchmark_roots(module_name: str) -> None:
    """Selected golden execution must fail closed rather than skip absent evidence."""
    result = _run_golden_probe(module_name, "golden", collect_only=False)
    output = result.stdout + result.stderr

    assert result.returncode == 1, output
    assert "VNTYPER_SIM_ROOT must name an explicit golden corpus root" in output
    assert " skipped" not in output


def _run_golden_probe(module_name: str, marker: str, *, collect_only: bool) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment.pop("VNTYPER_SIM_ROOT", None)
    environment.pop("VNTYPER_ADVNTR_ROOT", None)
    argv = [
        sys.executable,
        "-m",
        "pytest",
        "-m",
        marker,
        "-q",
        "-o",
        "log_cli=false",
        "--maxfail=1",
    ]
    if collect_only:
        argv.append("--collect-only")
    argv.append(str(UNIT_DIR.parent / "golden" / module_name))
    if collect_only:
        argv.append(f"{__file__}::test_root_pytest_ini_is_the_single_live_marker_authority")
    return subprocess.run(argv, cwd=REPO_ROOT, env=environment, capture_output=True, text=True, check=False)


def test_every_browser_file_declares_the_browser_marker() -> None:
    """Every test module under tests/browser/ must carry the ``browser`` marker.

    ``make test-browser`` runs ``pytest tests/browser -m browser``, so a module
    that forgets the marker is collected, deselected and never run - the same
    silent hole ``test_every_unit_file_is_selected_by_the_unit_marker`` closes for
    the unit tier, and ``--strict-markers`` cannot catch it because an omission is
    not a misspelling. Widening ``allowed`` above to admit the tier is what makes
    this necessary: the directory check alone now accepts an unmarked file.

    This one works off the source text rather than the session's collected items,
    because a ``pytest -m unit`` run never collects the browser tier and so has
    nothing to inspect.

    Raises:
        AssertionError: If a browser module defining tests declares no marker.
    """
    browser_dir = UNIT_DIR.parent / "browser"
    if not browser_dir.is_dir():
        pytest.skip("tests/browser is absent from this checkout")

    unmarked = []
    for path in sorted(browser_dir.rglob("test_*.py")):
        source = path.read_text(encoding="utf-8")
        if not _DEFINES_TESTS.search(source):
            continue
        if not _DECLARES_BROWSER_MARKER.search(source):
            unmarked.append(str(path.relative_to(browser_dir)))

    assert not unmarked, (
        f"These files under tests/browser/ define tests but never declare the `browser` marker, so "
        f"`pytest tests/browser -m browser` silently deselects them: {unmarked}. "
        "Add `pytestmark = pytest.mark.browser` immediately after the imports."
    )


def test_every_golden_file_declares_the_golden_marker() -> None:
    """Every test module under tests/golden/ must carry the ``golden`` marker.

    ``pytest -m golden tests/golden`` silently deselects an unmarked module, so
    admitting the directory as a known tier must be paired with the same
    source-level omission guard the browser tier carries.

    Raises:
        AssertionError: If a golden module defining tests declares no marker.
    """
    golden_dir = UNIT_DIR.parent / "golden"
    if not golden_dir.is_dir():
        pytest.skip("tests/golden is absent from this checkout")

    unmarked = []
    for path in sorted(golden_dir.rglob("test_*.py")):
        source = path.read_text(encoding="utf-8")
        if not _DEFINES_TESTS.search(source):
            continue
        if not _DECLARES_GOLDEN_MARKER.search(source):
            unmarked.append(str(path.relative_to(golden_dir)))

    assert not unmarked, (
        f"These files under tests/golden/ define tests but never declare the `golden` marker, so "
        f"`pytest tests/golden -m golden` silently deselects them: {unmarked}. "
        "Add `pytestmark = pytest.mark.golden` immediately after the imports."
    )


def test_collected_test_module_basenames_are_unique() -> None:
    """Every collected test module must have a repository-unique basename.

    Pytest's default import mode imports non-package test modules by basename. Two
    tier directories containing the same ``test_*.py`` basename therefore collide
    during collection before marker selection can deselect either module.

    Raises:
        AssertionError: If two collected test modules share a basename.
    """
    tests_root = UNIT_DIR.parent
    paths_by_basename: dict[str, list[str]] = {}
    for path in sorted(tests_root.rglob("test_*.py")):
        paths_by_basename.setdefault(path.name, []).append(str(path.relative_to(tests_root)))

    duplicates = {name: paths for name, paths in paths_by_basename.items() if len(paths) > 1}
    assert not duplicates, (
        f"Collected test modules share basenames under pytest's default import mode: {duplicates}. "
        "Rename each test module so its basename is unique across all tiers."
    )


def test_root_pytest_ini_is_the_single_live_marker_authority() -> None:
    pytest_ini = (REPO_ROOT / "pytest.ini").read_text(encoding="utf-8")
    pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "--strict-markers" in pytest_ini
    for marker in ("unit", "integration", "docker", "browser", "golden", "smoke", "slow"):
        assert f"{marker}:" in pytest_ini
    assert "[tool.pytest.ini_options]" not in pyproject
