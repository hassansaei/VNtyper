#!/usr/bin/env python3
"""
tests/unit/test_workflow_consistency.py

Guards invariants that span more than one GitHub Actions workflow file.

The base Docker image is content-addressed: docker-base.yml tags it with a hash of its
build inputs, and docker-build.yml recomputes the same hash to decide whether that base
exists and to pin the application build to it. Three separate `hashFiles(...)` calls
must therefore agree exactly. If one drifts, the tag silently changes meaning: the app
either rebuilds a base that already exists, or - worse - resolves a base built from
different inputs.

YAML has no way to share the expression, so it is duplicated by necessity. This test is
what stops the duplication from rotting. It is a plain unit test: it only reads files,
needs no Docker and no network.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

WORKFLOWS = Path(__file__).resolve().parents[2] / ".github" / "workflows"

# Matches the base-image hash expression specifically (it always starts with conda/**),
# not the unrelated test-data cache key that also uses hashFiles().
_BASE_HASH = re.compile(r"hashFiles\((\s*'conda/\*\*'.*?)\)", re.DOTALL)


def _base_hash_expressions() -> dict[str, list[str]]:
    """Collect every base-image hashFiles() argument list, keyed by file name.

    Returns:
        dict: Workflow file name -> list of normalised argument strings.
    """
    found: dict[str, list[str]] = {}
    for path in sorted(WORKFLOWS.glob("*.yml")):
        matches = _BASE_HASH.findall(path.read_text(encoding="utf-8"))
        if matches:
            found[path.name] = [" ".join(match.split()) for match in matches]
    return found


def test_base_image_hash_inputs_are_identical() -> None:
    """Every base-image hashFiles() call must list exactly the same paths.

    Raises:
        AssertionError: If the workflows disagree on the base image's inputs.
    """
    found = _base_hash_expressions()
    if not found:
        # No .github/workflows/ here: a source checkout without it, or an installed
        # sdist. Nothing to guard, so skip rather than fail on an absent directory.
        pytest.skip("no GitHub Actions workflows present in this tree")

    distinct = {expression for expressions in found.values() for expression in expressions}
    assert len(distinct) == 1, (
        "workflows disagree on the base image's hash inputs, so the content-addressed "
        f"tag means different things in different jobs:\n{found}"
    )


def test_base_hash_covers_everything_that_changes_the_base() -> None:
    """The hash must cover every input that can change the base image's contents.

    `.dockerignore` and `reference/**` are easy to forget and were genuinely missed:
    excluding all of reference/ dropped three config-declared files from the image
    without changing the tag, so no rebuild was triggered.

    Raises:
        AssertionError: If a required path is absent from the hash inputs.
    """
    required = (
        "conda/**",
        "docker/Dockerfile.base",
        "docker/requirements-web.txt",
        "vntyper/scripts/install_references.py",
        "vntyper/scripts/install_references_config.json",
        "vntyper/dependencies/advntr/**",
        "reference/**",
        ".dockerignore",
    )
    expressions = _base_hash_expressions()
    if not expressions:
        pytest.skip("no GitHub Actions workflows present in this tree")
    sample = next(iter(expressions.values()))[0]

    missing = [path for path in required if f"'{path}'" not in sample]
    assert not missing, f"base image hash does not cover: {missing}"


def test_unit_coverage_matrix_and_patch_coverage_version_are_fixed() -> None:
    workflow = (WORKFLOWS / "ci-tests.yml").read_text(encoding="utf-8")
    assert "python-version: ['3.10', '3.11', '3.12', '3.13']" in workflow
    assert "run: make test-unit-cov" in workflow
    assert "matrix.python-version == '3.12'" in workflow
    assert "PATCH_COVERAGE_BASE" in workflow
    assert "`mypy vntyper/ docker/app/ scripts/`" in workflow
