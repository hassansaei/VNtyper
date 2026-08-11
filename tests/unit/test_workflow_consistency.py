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

import os
import re
import subprocess
from pathlib import Path

import pytest
import yaml

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


def _compatibility_step() -> dict[str, object]:
    workflow = yaml.safe_load((WORKFLOWS / "ci-tests.yml").read_text(encoding="utf-8"))
    matches = [
        step
        for step in workflow["jobs"]["test-unit"]["steps"]
        if step.get("name") == "Check real integration compatibility"
    ]
    assert len(matches) == 1
    return matches[0]


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True, text=True)
    return result.stdout.strip()


def _workflow_history(tmp_path: Path) -> tuple[Path, str]:
    repo = tmp_path / "repo"
    repo.mkdir(parents=True)
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "tests@example.invalid")
    _git(repo, "config", "user.name", "Tests")
    (repo / "tracked").write_text("base\n", encoding="utf-8")
    _git(repo, "add", "tracked")
    _git(repo, "commit", "-qm", "base")
    base = _git(repo, "rev-parse", "HEAD")
    _git(repo, "update-ref", "refs/remotes/origin/main", base)
    (repo / "tracked").write_text("head\n", encoding="utf-8")
    _git(repo, "commit", "-qam", "head")
    return repo, base


def _run_compatibility_step(tmp_path: Path, *, event_name: str, base_ref: str, before: str) -> tuple[int, str]:
    repo, base = _workflow_history(tmp_path)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    capture = tmp_path / "make-argv.txt"
    fake_make = fake_bin / "make"
    fake_make.write_text(
        """#!/bin/sh
set -eu
printf '%s\n' "$*" > "$WORKFLOW_CAPTURE"
base=${2#INTEGRATION_COMPAT_BASE=}
git rev-parse --verify "$base^{commit}" >/dev/null
git merge-base --is-ancestor "$base" HEAD
""",
        encoding="utf-8",
    )
    fake_make.chmod(0o755)
    step = _compatibility_step()
    result = subprocess.run(
        ["bash", "-c", str(step["run"])],
        cwd=repo,
        env={
            **os.environ,
            "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
            "EVENT_NAME": event_name,
            "PR_BASE_REF": base_ref,
            "PUSH_BEFORE": base if before == "__BASE__" else before,
            "WORKFLOW_CAPTURE": str(capture),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    return result.returncode, capture.read_text(encoding="utf-8").strip() if capture.exists() else ""


def test_integration_compatibility_workflow_step_is_single_strict_and_full_history() -> None:
    """Catch duplicate/soft checker execution or a checkout that cannot prove ancestry."""
    step = _compatibility_step()
    assert step["if"] == "matrix.python-version == '3.12'"
    assert step["env"] == {
        "EVENT_NAME": "${{ github.event_name }}",
        "PR_BASE_REF": "${{ github.base_ref }}",
        "PUSH_BEFORE": "${{ github.event.before }}",
    }
    run = str(step["run"])
    assert 'make check-integration-compatibility INTEGRATION_COMPAT_BASE="$BASE_REVISION"' in run
    assert "continue-on-error" not in step
    assert "|| true" not in run
    workflow = yaml.safe_load((WORKFLOWS / "ci-tests.yml").read_text(encoding="utf-8"))
    checkout = workflow["jobs"]["test-unit"]["steps"][0]
    assert checkout["with"]["fetch-depth"] == 0


def test_integration_compatibility_workflow_selects_real_pr_and_push_bases(tmp_path: Path) -> None:
    """Catch swapping PR target history with the direct-push before SHA."""
    pr_status, pr_argv = _run_compatibility_step(
        tmp_path / "pr", event_name="pull_request", base_ref="main", before="unused"
    )
    assert pr_status == 0
    assert pr_argv == "check-integration-compatibility INTEGRATION_COMPAT_BASE=origin/main"

    push_status, push_argv = _run_compatibility_step(
        tmp_path / "push", event_name="push", base_ref="unused", before="__BASE__"
    )
    assert push_status == 0
    push_base = push_argv.partition("INTEGRATION_COMPAT_BASE=")[2]
    assert re.fullmatch(r"[0-9a-f]{40}", push_base)


@pytest.mark.parametrize(
    ("event_name", "base_ref", "before"),
    [("workflow_dispatch", "", ""), ("push", "", "")],
)
def test_integration_compatibility_workflow_rejects_unsupported_or_empty_base(
    tmp_path: Path, event_name: str, base_ref: str, before: str
) -> None:
    """Catch allowing an event without an authoritative comparison base."""
    status, argv = _run_compatibility_step(tmp_path, event_name=event_name, base_ref=base_ref, before=before)
    assert status != 0
    assert argv == ""


@pytest.mark.parametrize(
    ("event_name", "base_ref", "before", "expected_base"),
    [
        ("push", "", "0" * 40, "0" * 40),
        ("pull_request", "missing", "", "origin/missing"),
    ],
)
def test_integration_compatibility_workflow_rejects_unreachable_bases(
    tmp_path: Path,
    event_name: str,
    base_ref: str,
    before: str,
    expected_base: str,
) -> None:
    """Catch accepting all-zero or otherwise unreachable comparison revisions."""
    status, argv = _run_compatibility_step(tmp_path, event_name=event_name, base_ref=base_ref, before=before)
    assert status != 0
    assert argv == f"check-integration-compatibility INTEGRATION_COMPAT_BASE={expected_base}"
