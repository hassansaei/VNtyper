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

import ast
import re
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = REPO_ROOT / "pyproject.toml"
CONDA_ENV = REPO_ROOT / "conda" / "environment_vntyper.yml"
CI_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ci-tests.yml"
SMOKE_TEST = REPO_ROOT / "tests" / "docker" / "test_image_structure.py"
WEB_REQUIREMENTS = REPO_ROOT / "docker" / "requirements-web.txt"
DOCKER_APP = REPO_ROOT / "docker" / "app"
VERSION_MODULE = REPO_ROOT / "vntyper" / "version.py"
CITATION = REPO_ROOT / "CITATION.cff"
CHANGELOG = REPO_ROOT / "docs" / "about" / "changelog.md"


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


def test_release_metadata_matches_the_package_version() -> None:
    """Citation and changelog release versions must match ``__version__``."""
    version_tree = ast.parse(_read(VERSION_MODULE))
    package_version = next(
        node.value.value
        for node in version_tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
    )
    citation_version = _search(r'^version:\s*"([^"]+)"', _read(CITATION), "the citation version").group(1)
    changelog_version = _search(
        r"^##\s+([0-9]+\.[0-9]+\.[0-9]+)(?:\s|$)", _read(CHANGELOG), "the latest changelog version"
    ).group(1)
    assert citation_version == package_version, (
        f"CITATION.cff is {citation_version}, but vntyper/version.py is {package_version}"
    )
    assert changelog_version == package_version, (
        f"docs/about/changelog.md is {changelog_version}, but vntyper/version.py is {package_version}"
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


# ---------------------------------------------------------------------------
# The web service's dependencies vs the versions the Docker image installs
# ---------------------------------------------------------------------------
# `docker/requirements-web.txt` is what the base image installs for the FastAPI +
# Celery service in `docker/app/`, and every line of it is pinned. `pyproject.toml`
# declares the same packages twice more - in the `web` extra, so `docker/app` can be
# run outside the container, and in `dev`, which is what the web unit tier resolves
# from. When those carried lower bounds instead of pins, the tier installed whatever
# PyPI offered that day and tested a framework the deployment never runs.
#
# The direction of the fix is fixed by AGENTS.md trap 10: `requirements-web.txt` is
# an input to the base image's content hash, so editing it turns the application
# image build red until a new base is published. It is therefore the authority here,
# and `pyproject.toml` is the side that moves.

# Third-party modules `docker/app` imports that `requirements-web.txt` does not name,
# mapped to the pinned distribution that supplies them. An entry is an allowance, not
# an exemption: `test_undeclared_import_allowances_are_all_still_needed` deletes the
# excuse the moment it stops applying.
#
#   starlette / kombu - hard, non-optional dependencies of `fastapi` and `celery`.
#     Both are pinned, both bound the version of their own dependency in their
#     metadata, and neither can install without pulling these in. Naming them here
#     would pin a transitive version independently of the package that owns it.
#
#   bcrypt - NOT of the same kind, and the reason this guard was written.
#     `docker/app/utils.py` imports it directly, but `requirements-web.txt` declares
#     only `passlib[bcrypt]==1.7.4`, and nothing under `docker/app/` imports passlib
#     any more. The import therefore survives on an *extra* of a distribution that is
#     otherwise unused - one dependency edit away from an ImportError at container
#     start. The fix is to declare `bcrypt` directly and drop `passlib`, but
#     `requirements-web.txt` is an input to the base image's content hash (AGENTS.md
#     trap 10), so editing it fails every application image build until a new base is
#     published on `main`. FOLLOW-UP: make that edit as part of the next base-image
#     rebuild, then delete this entry - the test above will insist on it.
UNDECLARED_IMPORT_ALLOWANCES: dict[str, str] = {
    "starlette": "fastapi",
    "kombu": "celery",
    "bcrypt": "passlib",
}


def _requirement_strings(block: str) -> list[str]:
    """Extract the double-quoted requirement strings from a TOML array body.

    Args:
        block: The text between the array's brackets.

    Returns:
        list[str]: Requirement strings, in declaration order.
    """
    return [item.strip() for item in re.findall(r'"([^"]+)"', block)]


def web_requirements() -> list[str]:
    """Return the requirement strings pinned in ``docker/requirements-web.txt``.

    Returns:
        list[str]: One PEP 508 requirement per non-empty, non-comment line.
    """
    found = []
    for raw in _read(WEB_REQUIREMENTS).splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            found.append(line)
    assert found, "docker/requirements-web.txt declared nothing - has the file format changed?"
    return found


def optional_dependencies(extra: str) -> list[str]:
    """Return the requirement strings in one ``[project.optional-dependencies]`` group.

    Args:
        extra: The group's name, e.g. ``"web"`` or ``"dev"``.

    Returns:
        list[str]: One PEP 508 requirement per entry.
    """
    match = re.search(rf"^{re.escape(extra)} = \[(.*?)^\]", _read(PYPROJECT), re.MULTILINE | re.DOTALL)
    assert match, f"could not find the optional-dependency group {extra!r} - has the file format changed?"
    return _requirement_strings(match.group(1))


def docker_app_third_party_imports() -> dict[str, list[str]]:
    """Return the third-party modules imported anywhere under ``docker/app``.

    Relative imports are the service's own modules and are skipped; anything in
    ``sys.stdlib_module_names`` comes with the interpreter. What is left has to be
    installed by something, and ``requirements-web.txt`` is the only thing that
    installs anything into the image for this service.

    Returns:
        dict[str, list[str]]: Top-level module name -> the files importing it.
    """
    found: dict[str, set[str]] = {}
    for path in sorted(DOCKER_APP.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            else:
                continue
            for name in names:
                top = name.split(".")[0]
                if top not in sys.stdlib_module_names:
                    found.setdefault(top, set()).add(path.name)
    assert found, "no imports found under docker/app - has the layout changed?"
    return {module: sorted(files) for module, files in found.items()}


@pytest.mark.parametrize("extra", ["web", "dev"])
def test_pyproject_pins_the_web_versions_the_image_installs(extra: str) -> None:
    """Every web package must be declared identically in pyproject and in the image.

    A lower bound here rather than a pin means the unit tier resolves a newer
    framework than the deployed image ever runs, so a passing web test says nothing
    about the deployment. Extras are compared too: ``fastapi`` and
    ``fastapi[standard]`` install different sets of packages.

    Args:
        extra: The optional-dependency group being checked.
    """
    from packaging.requirements import Requirement

    declared = {_normalise(Requirement(line).name): Requirement(line) for line in optional_dependencies(extra)}

    problems = []
    for line in web_requirements():
        image = Requirement(line)
        name = _normalise(image.name)
        local = declared.get(name)
        if local is None:
            problems.append(f"{name}: the image installs '{line}', the {extra!r} extra does not declare it at all")
            continue
        if local.specifier != image.specifier:
            problems.append(
                f"{name}: the {extra!r} extra requires '{local.specifier}', the image installs '{image.specifier}'"
            )
        if local.extras != image.extras:
            problems.append(
                f"{name}: the {extra!r} extra requests extras {sorted(local.extras)}, "
                f"the image installs {sorted(image.extras)}"
            )

    assert not problems, (
        f"the {extra!r} extra does not install what docker/requirements-web.txt installs, so the "
        "web tests run against a different stack than the image:\n  " + "\n  ".join(problems)
    )


def test_web_extra_declares_nothing_the_image_does_not_install() -> None:
    """The ``web`` extra is a mirror, so it must not grow entries of its own.

    ``dev`` may legitimately add test-only packages; ``web`` exists to reproduce the
    image's service environment, and anything extra in it is drift in the other
    direction.
    """
    from packaging.requirements import Requirement

    installed = {_normalise(Requirement(line).name) for line in web_requirements()}
    surplus = sorted(
        _normalise(Requirement(line).name)
        for line in optional_dependencies("web")
        if _normalise(Requirement(line).name) not in installed
    )
    assert not surplus, (
        f"the 'web' extra declares {surplus}, which docker/requirements-web.txt does not "
        "install - either add them to the image's requirements or drop them here"
    )


def test_docker_app_imports_are_declared_in_the_web_requirements() -> None:
    """Every third-party module ``docker/app`` imports must be a declared dependency.

    An import that works only because some other package happens to pull it in is
    one dependency bump away from an ImportError at container start. Modules that
    are guaranteed by a pinned distribution rather than named by one are recorded in
    ``UNDECLARED_IMPORT_ALLOWANCES``, which states what guarantees them.
    """
    from packaging.requirements import Requirement

    declared = {_normalise(Requirement(line).name) for line in web_requirements()}

    undeclared = {}
    for module, files in docker_app_third_party_imports().items():
        if _normalise(module) in declared:
            continue
        provider = UNDECLARED_IMPORT_ALLOWANCES.get(module)
        if provider is not None and _normalise(provider) in declared:
            continue
        undeclared[module] = files

    assert not undeclared, (
        "docker/app imports modules docker/requirements-web.txt does not declare, so the "
        "image only has them by accident:\n  "
        + "\n  ".join(f"{module} (imported by {', '.join(files)})" for module, files in sorted(undeclared.items()))
    )


def test_undeclared_import_allowances_are_all_still_needed() -> None:
    """Each allowance must still describe a real, still-undeclared import.

    Keeps the allowance table from outliving the situation it documents: once a
    module is declared directly, or is no longer imported, the entry has to go.
    """
    from packaging.requirements import Requirement

    declared = {_normalise(Requirement(line).name) for line in web_requirements()}
    imported = docker_app_third_party_imports()

    stale = []
    for module, provider in UNDECLARED_IMPORT_ALLOWANCES.items():
        if module not in imported:
            stale.append(f"{module}: no longer imported under docker/app - drop the allowance")
        elif _normalise(module) in declared:
            stale.append(f"{module}: now declared in requirements-web.txt - drop the allowance")
        elif _normalise(provider) not in declared:
            stale.append(f"{module}: allowed on the strength of '{provider}', which the image no longer installs")

    assert not stale, "UNDECLARED_IMPORT_ALLOWANCES is out of date:\n  " + "\n  ".join(stale)
