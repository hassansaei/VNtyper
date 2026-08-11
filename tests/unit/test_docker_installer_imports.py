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


def _module_scope_vntyper_imports(module_name: str) -> set[str]:
    """Collect the last path component of every module-scope `from vntyper.X import ...`.

    Only `tree.body` (top-level statements) is walked, so an import nested inside a
    function or class - the house style for anything that must not join the base-image
    content hash - is deliberately invisible here.

    Args:
        module_name: File stem under `vntyper/scripts/`, e.g. ``"install_references"``.

    Returns:
        set[str]: Last dotted component of each module-scope `vntyper.` import target.
    """
    source = (SCRIPTS_DIR / f"{module_name}.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    return {
        node.module.split(".")[-1]
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("vntyper.")
    }


def _modules_copied_into_refs_stage() -> set[str]:
    """Return every `vntyper/scripts/<name>.py` the Dockerfile copies anywhere.

    Returns:
        set[str]: File stems copied from `vntyper/scripts/`.
    """
    return set(re.findall(r"vntyper/scripts/(\w+)\.py", DOCKERFILE_BASE.read_text(encoding="utf-8")))


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
