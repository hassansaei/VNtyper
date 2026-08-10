"""
tests/unit/test_makefile_recipes.py

Guards the Makefile recipes that branch on an auto-detected optional tool.

Several targets locate a tool once at parse time::

    ACTIONLINT ?= $(shell command -v actionlint 2>/dev/null)

and then branch on whether the lookup found anything. The trap is that make
substitutes the *empty string* when it did not, so a recipe line consisting of
nothing but ``$(TOOL);`` becomes a bare ``;`` -- and ``;`` with no command in
front of it is a **syntax error**, not a no-op. Because the whole ``if ... fi``
reaches the shell as one compound command, that error is raised while parsing,
before the guard ``[ -n "$(TOOL)" ]`` ever runs. The fallback branch the guard
exists to select is therefore unreachable, and the target fails on exactly the
machines the fallback was written for. ``sh`` and ``bash`` both reject it.

This bit: ``make ci-local`` -- the gate AGENTS.md requires for any change under
``.github/workflows/`` -- died at its first step on every machine without
``actionlint`` installed, reporting only ``/bin/sh: 2: Syntax error: ";"
unexpected``. CI never caught it because no workflow runs actionlint; it is a
local-only target.

The fix is to quote the expansion (``"$(TOOL)"``), which makes the branch a
syntactically valid word. The quotes are load-bearing -- hence this test.

Pure unit test: it reads the Makefile, asks make to expand recipes with
``--dry-run`` (no target is executed), and syntax-checks the result with
``sh -n`` (which parses without executing). No network, no Docker, no data.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
MAKEFILE = REPO_ROOT / "Makefile"

# `TOOL ?= $(shell command -v tool ...)` -- a lookup that yields "" when absent.
_AUTO_DETECTED = re.compile(
    r"^([A-Z_][A-Z0-9_]*)\s*\?=\s*\$\(shell\s+command\s+-v\b",
    re.MULTILINE,
)

# A target line: `name:` or `name: prereqs`, but not `name := value`.
_TARGET = re.compile(r"^([a-zA-Z0-9_][a-zA-Z0-9_.-]*)\s*:(?!=)")


def _auto_detected_vars() -> list[str]:
    """Names of make variables assigned from a `command -v` lookup.

    Returns:
        list: Variable names, in the order they appear in the Makefile.
    """
    return _AUTO_DETECTED.findall(MAKEFILE.read_text())


def _recipes() -> dict[str, str]:
    """Map each target name to the raw text of its own recipe.

    Prerequisites are not followed; only the tab-indented lines that belong to
    the target itself are collected.

    Returns:
        dict: Target name -> recipe text (empty string for targets with none).
    """
    recipes: dict[str, str] = {}
    current: str | None = None
    for line in MAKEFILE.read_text().splitlines():
        if line.startswith("\t"):
            if current is not None:
                recipes[current] += line + "\n"
            continue
        match = _TARGET.match(line)
        if match:
            current = match.group(1)
            recipes.setdefault(current, "")
        elif line.strip() and not line.startswith("#"):
            # A non-recipe, non-comment line ends the current recipe.
            current = None
    return recipes


def _targets_using(var: str) -> list[str]:
    """Targets whose own recipe expands `var`.

    Args:
        var: Make variable name, without the `$(...)` wrapper.

    Returns:
        list: Matching target names.
    """
    needle = f"$({var})"
    return [name for name, body in _recipes().items() if needle in body]


def test_auto_detected_tools_are_discoverable() -> None:
    """The parser finds the auto-detect variables this test exists to cover.

    A refactor that renames or restyles the assignments would otherwise leave
    the checks below iterating over nothing and passing vacuously.
    """
    found = _auto_detected_vars()
    assert found, (
        "No `TOOL ?= $(shell command -v ...)` assignments found in the Makefile. "
        "Either they were removed (delete this test) or the pattern changed "
        "(update _AUTO_DETECTED) -- the recipe checks below are now vacuous."
    )
    assert "ACTIONLINT" in found, f"expected ACTIONLINT among auto-detected tools, got {found}"


@pytest.mark.parametrize("var", _auto_detected_vars())
def test_recipe_is_valid_shell_when_tool_is_absent(var: str) -> None:
    """Recipes must still parse when the auto-detected tool was not found.

    Renders each recipe with the tool variable forced empty -- the state on a
    machine that lacks it -- and syntax-checks the result. `make --dry-run`
    expands without running, and `sh -n` parses without executing, so this
    never invokes the tool or its fallback.

    Args:
        var: Name of the auto-detected make variable to blank out.
    """
    targets = _targets_using(var)
    assert targets, f"{var} is auto-detected but no recipe uses it; drop the variable"

    for target in targets:
        expanded = subprocess.run(
            ["make", "--dry-run", target, f"{var}="],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert expanded.returncode == 0, f"`make --dry-run {target} {var}=` failed:\n{expanded.stderr}"

        checked = subprocess.run(
            ["sh", "-n"],
            input=expanded.stdout,
            capture_output=True,
            text=True,
            check=False,
        )
        assert checked.returncode == 0, (
            f"Target `{target}` is not valid shell when {var} is empty "
            f"(i.e. when the tool is not installed):\n"
            f"  {checked.stderr.strip()}\n\n"
            f"An empty expansion in command position leaves a bare `;`. Quote it "
            f'-- `"$({var})"` -- so the branch stays a syntactically valid word.\n\n'
            f"Rendered recipe:\n{expanded.stdout}"
        )


def test_scripts_coverage_target_is_isolated_and_fixed_at_88() -> None:
    text = MAKEFILE.read_text(encoding="utf-8")
    recipe = _recipes()["test-scripts-cov"]
    assert re.search(r"^SCRIPTS_COVERAGE_TARGET\s*\?=\s*88$", text, re.MULTILINE)
    assert "--cov=scripts" in recipe
    assert "--cov-fail-under=$(SCRIPTS_COVERAGE_TARGET)" in recipe
    assert "--cov-report=term-missing" in recipe
    assert "coverage.xml" not in recipe
    assert "mktemp -d" in recipe
    assert "COVERAGE_FILE=" in recipe
    assert 'test -n "$$tmp_dir"' in recipe
    assert 'rm -rf -- "$$tmp_dir"' in recipe


def test_repository_coverage_thresholds_remain_independent() -> None:
    text = MAKEFILE.read_text(encoding="utf-8")
    assert re.search(r"^COVERAGE_TARGET\s*\?=\s*86$", text, re.MULTILINE)
    assert re.search(r"^PATCH_COVERAGE_TARGET\s*\?=\s*80$", text, re.MULTILINE)


def test_mypy_runtime_and_test_scopes_are_deliberately_separate() -> None:
    recipes = _recipes()
    assert "mypy vntyper/ docker/app/ scripts/" in recipes["type-check"]
    assert "tests/" not in recipes["type-check"]
    assert "mypy vntyper/ tests/" in recipes["type-check-all"]
    assert "scripts/" not in recipes["type-check-all"]


def test_ruff_recipes_terminate_options_before_repository_scopes() -> None:
    """Every Ruff invocation treats RUFF_PATHS entries only as path operands."""
    scoped_lines = [
        line.strip()
        for recipe in _recipes().values()
        for line in recipe.splitlines()
        if line.strip().startswith("ruff ") and "$(RUFF_PATHS)" in line
    ]
    assert scoped_lines
    assert all(re.search(r"(?:^|\s)--\s+\$\(RUFF_PATHS\)(?:\s|$)", line) for line in scoped_lines), scoped_lines


def _write_executable(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)


def _run_scripts_coverage_with_fake_pytest(
    tmp_path: Path,
    pytest_status: int,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    capture = tmp_path / "coverage-file.txt"
    _write_executable(
        fake_bin / "pytest",
        """#!/bin/sh
printf '%s\n' "$COVERAGE_FILE" > "$COVERAGE_CAPTURE"
exit "$FAKE_PYTEST_STATUS"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"
    env["COVERAGE_CAPTURE"] = str(capture)
    env["FAKE_PYTEST_STATUS"] = str(pytest_status)
    result = subprocess.run(
        ["make", "--no-print-directory", "test-scripts-cov"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    return result, Path(capture.read_text(encoding="utf-8").strip())


def test_scripts_coverage_removes_temporary_data_after_success(tmp_path: Path) -> None:
    result, coverage_file = _run_scripts_coverage_with_fake_pytest(tmp_path, pytest_status=0)

    assert result.returncode == 0
    assert coverage_file.name == ".coverage"
    assert not coverage_file.parent.exists()


def test_scripts_coverage_removes_temporary_data_and_preserves_failure(tmp_path: Path) -> None:
    result, coverage_file = _run_scripts_coverage_with_fake_pytest(tmp_path, pytest_status=17)

    assert result.returncode == 2
    assert "Error 17" in result.stderr
    assert coverage_file.name == ".coverage"
    assert not coverage_file.parent.exists()


def test_scripts_coverage_does_not_run_pytest_when_mktemp_fails(tmp_path: Path) -> None:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    capture = tmp_path / "coverage-file.txt"
    _write_executable(fake_bin / "mktemp", "#!/bin/sh\nexit 23\n")
    _write_executable(
        fake_bin / "pytest",
        """#!/bin/sh
printf '%s\n' "$COVERAGE_FILE" > "$COVERAGE_CAPTURE"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"
    env["COVERAGE_CAPTURE"] = str(capture)

    result = subprocess.run(
        ["make", "--no-print-directory", "test-scripts-cov"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 2
    assert "Error 23" in result.stderr
    assert not capture.exists(), "pytest received a COVERAGE_FILE path after mktemp failed"
