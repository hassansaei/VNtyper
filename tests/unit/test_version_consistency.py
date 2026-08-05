#!/usr/bin/env python3
"""
tests/unit/test_version_consistency.py

Keeps the PyPI package and the Docker image in agreement about versions.

The supported interpreter is asserted in six places, in four different file formats:

    conda/environment_vntyper.yml            python=X.Y.Z   (what the image RUNS)
    pyproject.toml  [project]                requires-python = ">=X.Y"
    pyproject.toml  [project] classifiers    Programming Language :: Python :: X.Y
    pyproject.toml  [tool.ruff]              target-version = "pyXY"
    pyproject.toml  [tool.mypy]              python_version = "X.Y"
    .github/workflows/ci-tests.yml           matrix python-version: [...]
    tests/docker/test_image_structure.py     expected interpreter in the image

Nothing makes them agree on its own, and disagreement is silent and expensive: shipping
an interpreter the test matrix never exercises, or letting ruff accept syntax the
declared floor cannot run. This test is the mechanism that keeps them in sync.

Deliberately dependency-free (no tomllib - the 3.10 matrix leg has none, and tomli is
not a dependency) and pure file reads, so it belongs in the fast unit tier.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = REPO_ROOT / "pyproject.toml"
CONDA_ENV = REPO_ROOT / "conda" / "environment_vntyper.yml"
CI_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ci-tests.yml"
SMOKE_TEST = REPO_ROOT / "tests" / "docker" / "test_image_structure.py"


def _read(path: Path) -> str:
    """Return a file's text, skipping the test if it is absent.

    Args:
        path: File to read.

    Returns:
        str: File contents.
    """
    if not path.is_file():
        pytest.skip(f"{path.relative_to(REPO_ROOT)} not present in this tree")
    return path.read_text(encoding="utf-8")


def _search(pattern: str, text: str, what: str) -> re.Match[str]:
    """Search for a required pattern, failing with a useful message if absent.

    Args:
        pattern: Regular expression.
        text: Text to search.
        what: Human-readable description used in the failure message.

    Returns:
        re.Match: The match object.
    """
    match = re.search(pattern, text, re.MULTILINE)
    assert match, f"could not find {what} - has the file format changed?"
    return match


def conda_python() -> tuple[int, int]:
    """Return the (major, minor) Python pinned in the vntyper conda environment.

    Returns:
        tuple: Major and minor version.
    """
    match = _search(r"^\s*-\s*python=(\d+)\.(\d+)", _read(CONDA_ENV), "the conda python pin")
    return int(match.group(1)), int(match.group(2))


def requires_python_floor() -> tuple[int, int]:
    """Return the (major, minor) floor declared by ``requires-python``.

    Returns:
        tuple: Major and minor version.
    """
    match = _search(r'^requires-python\s*=\s*">=(\d+)\.(\d+)"', _read(PYPROJECT), "requires-python")
    return int(match.group(1)), int(match.group(2))


def classifier_versions() -> set[tuple[int, int]]:
    """Return the Python versions listed in the trove classifiers.

    Returns:
        set: (major, minor) pairs.
    """
    text = _read(PYPROJECT)
    return {
        (int(major), int(minor))
        for major, minor in re.findall(r'"Programming Language :: Python :: (\d+)\.(\d+)"', text)
    }


def matrix_versions() -> set[tuple[int, int]]:
    """Return the Python versions in the CI unit-test matrix.

    Returns:
        set: (major, minor) pairs.
    """
    match = _search(r"^\s*python-version:\s*\[(.+)\]", _read(CI_WORKFLOW), "the CI python matrix")
    return {(int(major), int(minor)) for major, minor in re.findall(r"'(\d+)\.(\d+)'", match.group(1))}


def test_conda_python_satisfies_requires_python() -> None:
    """The interpreter the image ships must satisfy the package's declared floor."""
    assert conda_python() >= requires_python_floor(), (
        f"conda env ships Python {conda_python()} but requires-python declares >={requires_python_floor()}"
    )


def test_ci_matrix_covers_the_shipped_interpreter() -> None:
    """CI must test the interpreter version the Docker image actually runs.

    Otherwise the version users run is the one version nothing exercises.
    """
    assert conda_python() in matrix_versions(), (
        f"the image runs Python {conda_python()}, which is not in the CI matrix {sorted(matrix_versions())}"
    )


def test_matrix_lowest_equals_requires_python_floor() -> None:
    """The matrix must actually exercise the oldest version the package claims."""
    assert min(matrix_versions()) == requires_python_floor(), (
        f"requires-python declares >={requires_python_floor()} but the lowest tested "
        f"version is {min(matrix_versions())}"
    )


def test_classifiers_match_the_matrix() -> None:
    """Trove classifiers must advertise exactly the versions CI tests."""
    assert classifier_versions() == matrix_versions(), (
        f"classifiers {sorted(classifier_versions())} != CI matrix {sorted(matrix_versions())}"
    )


def test_ruff_target_matches_requires_python() -> None:
    """Ruff must reject syntax newer than the declared floor can run."""
    major, minor = requires_python_floor()
    match = _search(r'^target-version\s*=\s*"py(\d+)"', _read(PYPROJECT), "ruff target-version")
    assert match.group(1) == f"{major}{minor}", (
        f'ruff target-version is "py{match.group(1)}" but requires-python is >={major}.{minor}'
    )


def test_mypy_python_version_matches_requires_python() -> None:
    """mypy must type-check against the declared floor."""
    major, minor = requires_python_floor()
    match = _search(r'^python_version\s*=\s*"(\d+)\.(\d+)"', _read(PYPROJECT), "mypy python_version")
    assert (int(match.group(1)), int(match.group(2))) == (major, minor), (
        f"mypy targets {match.group(0)} but requires-python is >={major}.{minor}"
    )


def test_smoke_test_expects_the_shipped_interpreter() -> None:
    """The image smoke tier must assert the version the conda environment installs."""
    major, minor = conda_python()
    match = _search(r'\("vntyper",\s*"(\d+\.\d+)"\)', _read(SMOKE_TEST), "the smoke test's expected version")
    assert match.group(1) == f"{major}.{minor}", (
        f"the image smoke test expects Python {match.group(1)} but the conda env pins {major}.{minor}"
    )


# ---------------------------------------------------------------------------
# Package dependencies vs the conda environment the Docker image installs
# ---------------------------------------------------------------------------
# The most consequential drift this guards: pyproject pinned numpy>=1.26.0,<2.0.0 while
# the conda environment installed numpy=2.0.2, so `pip install .` inside the image
# downgraded conda's numpy and the published image silently ran 1.26.4 - a different
# numerical stack from the one the environment file declares.


def _normalise(name: str) -> str:
    """Normalise a distribution name for comparison across PyPI and conda.

    Args:
        name: Raw package name.

    Returns:
        str: Lower-cased name with underscores folded to hyphens.
    """
    return name.strip().lower().replace("_", "-")


def pyproject_dependencies() -> dict[str, str]:
    """Return the runtime dependencies declared in ``[project] dependencies``.

    Returns:
        dict: Normalised package name -> PEP 440 specifier (may be empty).
    """
    text = _read(PYPROJECT)
    # DOTALL: the array spans lines. The `^dependencies` anchor keeps this to
    # [project] dependencies and away from [project.optional-dependencies].
    match = re.search(r"^dependencies = \[(.*?)^\]", text, re.MULTILINE | re.DOTALL)
    assert match, "could not find [project] dependencies - has the file format changed?"
    block = match.group(1)
    found: dict[str, str] = {}
    for raw in re.findall(r'"([^"]+)"', block):
        match = re.match(r"^([A-Za-z0-9._-]+)\s*(.*)$", raw.strip())
        if match:
            found[_normalise(match.group(1))] = match.group(2).strip()
    return found


def conda_dependencies() -> dict[str, str]:
    """Return the pinned packages in the vntyper conda environment.

    Returns:
        dict: Normalised package name -> exact pinned version.
    """
    found: dict[str, str] = {}
    for line in _read(CONDA_ENV).splitlines():
        # e.g. "  - bioconda::fastp=0.23.4" or "  - pandas=2.2.2"
        match = re.match(r"^\s*-\s*(?:[A-Za-z0-9_-]+::)?([A-Za-z0-9._-]+)=([\w.]+)\s*$", line)
        if match:
            found[_normalise(match.group(1))] = match.group(2)
    return found


def test_conda_versions_satisfy_pyproject_specifiers() -> None:
    """Every shared dependency's conda pin must satisfy the package's own specifier.

    When it does not, `pip install .` inside the image resolves a *different* version
    than the environment declares, and the shipped image stops matching the recipe.
    """
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    declared = pyproject_dependencies()
    pinned = conda_dependencies()
    shared = sorted(set(declared) & set(pinned))
    assert shared, "no shared dependencies found - have the file formats changed?"

    conflicts = []
    for name in shared:
        specifier = declared[name]
        if not specifier:
            continue
        if Version(pinned[name]) not in SpecifierSet(specifier):
            conflicts.append(f"{name}: conda pins {pinned[name]}, pyproject requires {specifier}")

    assert not conflicts, (
        "the conda environment installs versions the package forbids, so pip would "
        "replace them inside the image:\n  " + "\n  ".join(conflicts)
    )


def test_image_required_binaries_are_in_the_conda_environment() -> None:
    """Every binary the image smoke tier requires must be declared in the environment.

    Catches a tool being dropped from the recipe before a 25-minute base build proves
    it at runtime.
    """
    required = re.findall(
        r'"([a-z0-9-]+)"', _search(r"^REQUIRED_BINARIES = \((.*?)\)", _read(SMOKE_TEST), "REQUIRED_BINARIES").group(1)
    )
    assert required, "could not parse REQUIRED_BINARIES from the smoke test"

    pinned = conda_dependencies()
    # `java` is provided by openjdk rather than a package of its own.
    provided = set(pinned) | ({"java"} if "openjdk" in pinned else set())
    missing = [tool for tool in required if _normalise(tool) not in provided]
    assert not missing, f"the image smoke tier requires {missing}, which the conda environment does not declare"
