#!/usr/bin/env python3
"""
tests/unit/test_docker_installer_imports.py

Guards the `refs` stage of `docker/Dockerfile.base` against a class of bug that no
other unit test can see: it copies `install_references.py` (and now `reference_bundle.py`)
into the image and runs the installer, but it copies *only* the files it names in `COPY`
instructions. If either module gains a module-scope `from vntyper....` import of a
sibling that the Dockerfile does not copy, the base image build fails at
`ModuleNotFoundError` - after conda envs, adVNTR and the reference genomes have already
been built, 25-120 minutes into the most expensive step in the whole pipeline.

Function-local imports are the house style for anything that must NOT join the base
image content-hash set (see `reference_registry` usage inside
`canonical_reference_keys`), so this only inspects imports at module (function/class)
scope - `tree.body` - not imports nested inside a function body.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCKERFILE_BASE = REPO_ROOT / "docker" / "Dockerfile.base"
SCRIPTS_DIR = REPO_ROOT / "vntyper" / "scripts"

# The modules the `refs` stage's installer invocation is known to need. Kept in sync by
# hand with the module docstrings of the two files themselves; a third module added to
# either import list must be added to both this tuple and the Dockerfile's COPY lines.
_INSTALLER_MODULES = ("install_references", "reference_bundle")


#: The `refs` stage alone. A `COPY vntyper/scripts/x.py` in `envs` or `base` provisions
#: nothing for the installer, so scanning the whole file would convince this guard that
#: an unreachable module is reachable. Matches from `FROM ... AS refs` to the next
#: top-level `FROM` or end of file.
_REFS_STAGE = re.compile(r"^FROM\s+\S+\s+AS\s+refs\b(.*?)(?=^FROM\s|\Z)", re.DOTALL | re.MULTILINE)


def _vntyper_imports_in(source: str) -> set[str]:
    """Collect the last path component of every module-scope `vntyper.` import.

    Both spellings count. `import vntyper.scripts.foo` raises the identical
    `ModuleNotFoundError` as `from vntyper.scripts.foo import bar`, so collecting only
    `ast.ImportFrom` would leave a one-character style choice able to cost a 25-120
    minute base build.

    Only `tree.body` (top-level statements) is walked, so an import nested inside a
    function or class - the house style for anything that must not join the base-image
    content hash - is deliberately invisible here.

    Args:
        source: Python source text.

    Returns:
        set[str]: Last dotted component of each module-scope `vntyper.` import target.
    """
    tree = ast.parse(source)
    imported: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("vntyper."):
            imported.add(node.module.split(".")[-1])
        elif isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[-1] for alias in node.names if alias.name.startswith("vntyper."))
    return imported


def _module_scope_vntyper_imports(module_name: str) -> set[str]:
    """Read one installer module and collect its module-scope `vntyper.` imports.

    Args:
        module_name: File stem under `vntyper/scripts/`, e.g. ``"install_references"``.

    Returns:
        set[str]: Last dotted component of each module-scope `vntyper.` import target.
    """
    return _vntyper_imports_in((SCRIPTS_DIR / f"{module_name}.py").read_text(encoding="utf-8"))


def _refs_stage_text(dockerfile: str) -> str:
    """Return the body of the `refs` stage of a Dockerfile.

    Args:
        dockerfile: Full Dockerfile text.

    Returns:
        str: Everything between `FROM ... AS refs` and the next top-level `FROM`.

    Raises:
        AssertionError: If the stage cannot be found. Returning an empty string would
            make every parametrisation of the guard below pass vacuously - the guard
            would report that nothing is copied and that nothing is therefore missing.
    """
    match = _REFS_STAGE.search(dockerfile)
    assert match, (
        "docker/Dockerfile.base has no `FROM ... AS refs` stage; this guard cannot tell "
        "which modules the installer can import and must not pass by default"
    )
    return match.group(1)


def _modules_copied_into_refs_stage() -> set[str]:
    """Return every `vntyper/scripts/<name>.py` the `refs` stage copies.

    Returns:
        set[str]: File stems copied from `vntyper/scripts/` inside the `refs` stage.
    """
    return set(re.findall(r"vntyper/scripts/(\w+)\.py", _refs_stage_text(DOCKERFILE_BASE.read_text(encoding="utf-8"))))


@pytest.mark.parametrize("module_name", _INSTALLER_MODULES)
def test_installer_module_scope_imports_are_copied_into_refs_stage(module_name: str) -> None:
    """Every module-scope `vntyper.` import the installer chain makes must be COPYed.

    This is the test a 25-120 minute Docker build would otherwise be the only thing to
    catch: `install_references.py` is run by file path inside a stage that installs no
    package, so an import the Dockerfile does not provision for raises
    `ModuleNotFoundError` at the most expensive step in the base build.

    Raises:
        AssertionError: If `module_name` imports a `vntyper.scripts.X` sibling that the
            `refs` stage does not copy.
    """
    if not DOCKERFILE_BASE.exists():
        pytest.skip("docker/Dockerfile.base not present in this checkout")

    imported = _module_scope_vntyper_imports(module_name)
    copied = _modules_copied_into_refs_stage()

    missing = imported - copied
    assert not missing, (
        f"{module_name}.py imports {sorted(missing)} at module scope, but the `refs` "
        "stage in docker/Dockerfile.base does not COPY it - the base build would fail "
        "with ModuleNotFoundError after conda envs and adVNTR are already built."
    )


class TestTheGuardItself:
    """Both defects below made the guard pass while the base build would have failed."""

    def test_a_plain_import_of_a_sibling_is_collected(self) -> None:
        """`import vntyper.scripts.foo` raises the same ModuleNotFoundError as the
        `from` spelling, and collecting only `ast.ImportFrom` could not see it."""
        assert _vntyper_imports_in("import vntyper.scripts.reference_bundle\n") == {"reference_bundle"}

    def test_an_aliased_plain_import_is_collected_by_its_real_name(self) -> None:
        assert _vntyper_imports_in("import vntyper.scripts.reference_bundle as rb\n") == {"reference_bundle"}

    def test_both_spellings_and_several_names_are_collected_together(self) -> None:
        source = "import os\nimport vntyper.scripts.a, vntyper.scripts.b\nfrom vntyper.scripts.c import thing\n"
        assert _vntyper_imports_in(source) == {"a", "b", "c"}

    def test_a_plain_import_of_the_package_itself_is_not_a_sibling(self) -> None:
        """`import vntyper` needs only `vntyper/__init__.py`, which the stage copies
        by a different COPY line; it names no `vntyper/scripts/` module."""
        assert _vntyper_imports_in("import vntyper\n") == set()

    def test_a_function_local_import_stays_invisible(self) -> None:
        """Deliberate: the house style for anything that must not join the base-image
        content-hash set (see `reference_registry` inside `canonical_reference_keys`)."""
        source = "def f():\n    from vntyper.scripts.reference_registry import keys\n    return keys\n"
        assert _vntyper_imports_in(source) == set()

    def test_only_the_refs_stage_is_scanned(self) -> None:
        """A COPY in `envs` or `base` provisions nothing for the installer, so treating
        it as evidence of reachability would let the guard bless an unreachable module."""
        dockerfile = (
            "FROM scratch AS envs\n"
            "COPY vntyper/scripts/elsewhere.py /tmp/\n"
            "FROM envs AS refs\n"
            "COPY vntyper/scripts/install_references.py /opt/ir/vntyper/scripts/\n"
            "FROM scratch AS base\n"
            "COPY vntyper/scripts/afterwards.py /tmp/\n"
        )
        body = _refs_stage_text(dockerfile)
        assert set(re.findall(r"vntyper/scripts/(\w+)\.py", body)) == {"install_references"}

    def test_the_final_stage_being_refs_is_still_captured(self) -> None:
        """`\\Z` rather than only a following `FROM`, so the last stage is not lost."""
        dockerfile = "FROM scratch AS envs\nFROM envs AS refs\nCOPY vntyper/scripts/install_references.py /opt/\n"
        assert "install_references" in _refs_stage_text(dockerfile)

    def test_a_dockerfile_with_no_refs_stage_fails_rather_than_passing_vacuously(self) -> None:
        with pytest.raises(AssertionError, match="no `FROM ... AS refs` stage"):
            _refs_stage_text("FROM scratch AS envs\nCOPY vntyper/scripts/install_references.py /opt/\n")

    def test_the_real_dockerfile_copies_the_two_installer_modules(self) -> None:
        """Anchors the scoped regex against the file it actually guards."""
        if not DOCKERFILE_BASE.exists():
            pytest.skip("docker/Dockerfile.base not present in this checkout")
        assert set(_INSTALLER_MODULES) <= _modules_copied_into_refs_stage()
