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
from typing import Any

import pytest
import yaml

pytestmark = pytest.mark.unit

WORKFLOWS = Path(__file__).resolve().parents[2] / ".github" / "workflows"

# Matches the base-image hash expression specifically (it always starts with conda/**),
# not the unrelated test-data cache key that also uses hashFiles().
_BASE_HASH = re.compile(r"hashFiles\((\s*'conda/\*\*'.*?)\)", re.DOTALL)

# docker-base.yml has exactly one `paths:` key (under `on: push:`); capture everything
# between it and the next top-level `workflow_dispatch:` sibling, then pull the quoted
# path entries out of that block. Reading `on:` with a YAML parser is a known trap -
# PyYAML 1.1 treats the bare word `on` as the boolean `True`, so this stays regex-based
# like `_BASE_HASH` above rather than switching to `yaml.safe_load()` for this one field.
_PUSH_PATHS_BLOCK = re.compile(r"paths:\n(.*?)\n\s*workflow_dispatch:", re.DOTALL)


def _docker_base_push_paths() -> list[str]:
    """Collect the quoted path entries from docker-base.yml's push `paths:` filter.

    Returns:
        list[str]: Path patterns in file order, or `[]` if the file or the filter is
        absent.
    """
    path = WORKFLOWS / "docker-base.yml"
    if not path.exists():
        return []
    match = _PUSH_PATHS_BLOCK.search(path.read_text(encoding="utf-8"))
    if not match:
        return []
    return re.findall(r"'([^']+)'", match.group(1))


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

    `.dockerignore` was easy to forget and was genuinely missed once: excluding all of
    reference/ dropped three config-declared files from the image without changing the
    tag, so no rebuild was triggered. The two `__init__.py` files were the sixth defect
    found at this seam: the `refs` stage COPYs them and
    `python -m vntyper.scripts.install_references` imports both, so they are image
    content, and they were in none of the five mirrors of this list.

    `reference/**` itself is deliberately absent from `required`: milestone 5 moved
    reference data to a published, checksummed bundle in berntpopp/vntyper-data, so
    `reference/` (now holding only `README.md`, `pseudonymize.py` and
    `pseudonymize_config.json`) no longer changes what the `refs` stage installs. The
    bundle pin (`install_references_config.json`'s asset/asset_sha256 fields) is already
    required below, and is what triggers a rebuild when a reference byte changes now.

    Raises:
        AssertionError: If a required path is absent from the hash inputs.
    """
    required = (
        "conda/**",
        "docker/Dockerfile.base",
        "docker/requirements-web.txt",
        "vntyper/__init__.py",
        "vntyper/scripts/__init__.py",
        "vntyper/scripts/install_references.py",
        "vntyper/scripts/install_references_config.json",
        "vntyper/scripts/install_references_logging.py",
        "vntyper/scripts/reference_bundle.py",
        "vntyper/scripts/reference_download.py",
        "vntyper/scripts/reference_integrity.py",
        "vntyper/scripts/reference_provenance.py",
        "vntyper/dependencies/advntr/**",
        ".dockerignore",
    )
    expressions = _base_hash_expressions()
    if not expressions:
        pytest.skip("no GitHub Actions workflows present in this tree")
    sample = next(iter(expressions.values()))[0]

    missing = [path for path in required if f"'{path}'" not in sample]
    assert not missing, f"base image hash does not cover: {missing}"


def test_docker_base_push_paths_cover_every_hash_input() -> None:
    """docker-base.yml's push `paths:` filter must list every base-image hash input.

    `test_base_image_hash_inputs_are_identical` only compares the three `hashFiles(...)`
    expressions to each other, so a path all three omitted together would still pass it -
    and nothing else compares the push `paths:` filter to the hash at all. If a hash
    input is missing from `paths:`, a commit touching only that file changes the content
    hash but does not fire `docker-base.yml` on push to `main`; the base then gets built
    lazily by the next `docker-build.yml` run instead of proactively.

    This is a subset check, not an equality check: `paths:` legitimately also lists
    `.github/workflows/docker-base.yml` itself, which is not a hash input.

    Raises:
        AssertionError: If a hash input is absent from the push paths filter.
    """
    expressions = _base_hash_expressions()
    if not expressions:
        pytest.skip("no GitHub Actions workflows present in this tree")

    hash_paths = re.findall(r"'([^']+)'", next(iter(expressions.values()))[0])
    push_paths = _docker_base_push_paths()
    assert push_paths, "docker-base.yml has no push paths: filter to compare against"

    missing = [path for path in hash_paths if path not in push_paths]
    assert not missing, (
        f"docker-base.yml's push paths: filter is missing hash input(s) {missing} - "
        "a commit touching only those files changes the content-addressed base tag "
        "without triggering a proactive rebuild on push to main"
    )


def test_makefile_base_inputs_mirror_every_hash_input() -> None:
    """`BASE_INPUTS` is the fifth mirror of the base-image input list, and the only one
    nothing compared to the hash.

    `make ci-local-docker` uses it to refuse building against the published `:latest`
    base when a base input has changed locally. A path missing from it means that guard
    silently stops firing for that file, and the local Docker check gives false
    assurance about an image CI will build differently.

    Globs are normalised away: the workflows hash `conda/**`, the Makefile passes `conda`
    to `git diff`, and those mean the same set of files.

    Raises:
        AssertionError: If a hash input has no counterpart in `BASE_INPUTS`.
    """
    makefile = WORKFLOWS.parents[1] / "Makefile"
    if not makefile.exists():
        pytest.skip("no Makefile present in this tree")
    expressions = _base_hash_expressions()
    if not expressions:
        pytest.skip("no GitHub Actions workflows present in this tree")

    block = re.search(
        r"^BASE_INPUTS :=(.*?)(?=\n[^\t\s])", makefile.read_text(encoding="utf-8"), re.DOTALL | re.MULTILINE
    )
    assert block, "Makefile has no BASE_INPUTS assignment to compare against"
    base_inputs = set(block.group(1).replace("\\\n", " ").split())

    hash_paths = re.findall(r"'([^']+)'", next(iter(expressions.values()))[0])
    missing = [path for path in hash_paths if path.removesuffix("/**") not in base_inputs]
    assert not missing, (
        f"Makefile BASE_INPUTS is missing base-image hash input(s) {missing} - "
        "`make ci-local-docker` would not notice that you changed them"
    )


def test_unit_coverage_matrix_and_patch_coverage_version_are_fixed() -> None:
    workflow = (WORKFLOWS / "ci-tests.yml").read_text(encoding="utf-8")
    parsed = yaml.safe_load(workflow)
    assert "python-version: ['3.10', '3.11', '3.12', '3.13']" in workflow
    assert "run: make test-unit-cov" in workflow
    assert "matrix.python-version == '3.12'" in workflow
    assert "PATCH_COVERAGE_BASE" in workflow
    assert "`mypy vntyper/ docker/app/ scripts/`" in workflow
    assert parsed["jobs"]["test-unit"]["timeout-minutes"] == 20


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


def _docker_workflow() -> dict[str, Any]:
    return yaml.safe_load((WORKFLOWS / "docker-build.yml").read_text(encoding="utf-8"))


def _docker_test_step() -> dict[str, object]:
    workflow = _docker_workflow()
    steps = workflow["jobs"]["build-and-test"]["steps"]
    matches = [step for step in steps if step.get("name") == "Run Docker tests"]
    assert len(matches) == 1
    return matches[0]


def _run_docker_test_step(
    tmp_path: Path,
    *,
    event_name: str,
    make_status: int = 0,
) -> tuple[subprocess.CompletedProcess[str], str]:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir(parents=True)
    capture = tmp_path / "make-argv.txt"
    fake_make = fake_bin / "make"
    fake_make.write_text(
        """#!/bin/sh
printf '%s\n' "$*" > "$WORKFLOW_CAPTURE"
exit "$WORKFLOW_MAKE_STATUS"
""",
        encoding="utf-8",
    )
    fake_make.chmod(0o755)
    step = _docker_test_step()
    result = subprocess.run(
        ["bash", "-c", str(step["run"])],
        env={
            **os.environ,
            "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
            "EVENT_NAME": event_name,
            "WORKFLOW_CAPTURE": str(capture),
            "WORKFLOW_MAKE_STATUS": str(make_status),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    return result, capture.read_text(encoding="utf-8").strip() if capture.exists() else ""


@pytest.mark.parametrize(
    ("event_name", "expected_make_argv"),
    [
        ("pull_request", "test-docker-quick"),
        ("push", "docker-cram-fixtures test-docker-fast"),
        ("schedule", "docker-cram-fixtures test-docker"),
        ("workflow_dispatch", "docker-cram-fixtures test-docker"),
    ],
)
def test_docker_workflow_selects_exact_event_tier_and_fixture_order(
    tmp_path: Path,
    event_name: str,
    expected_make_argv: str,
) -> None:
    """Swapping a tier or collecting CRAM before candidate fixture generation must fail."""
    result, argv = _run_docker_test_step(tmp_path, event_name=event_name)

    assert result.returncode == 0
    assert argv == expected_make_argv


def test_docker_workflow_rejects_unknown_events_and_propagates_make_failure(tmp_path: Path) -> None:
    """An unknown event or failed fixture/test target must never reach a green Docker check."""
    unknown, unknown_argv = _run_docker_test_step(tmp_path / "unknown", event_name="repository_dispatch")
    failed, failed_argv = _run_docker_test_step(
        tmp_path / "failure",
        event_name="push",
        make_status=17,
    )

    assert unknown.returncode != 0
    assert unknown_argv == ""
    assert failed.returncode == 17
    assert failed_argv == "docker-cram-fixtures test-docker-fast"


def test_docker_workflow_step_is_strict_ordered_and_has_full_matrix_budget() -> None:
    """Soft failure, early collection, or the old one-adVNTR job budget must fail."""
    workflow = _docker_workflow()
    job = workflow["jobs"]["build-and-test"]
    step = _docker_test_step()
    steps = job["steps"]
    names = [item.get("name") for item in steps]

    assert job["timeout-minutes"] == 120
    assert step["env"] == {
        "EVENT_NAME": "${{ github.event_name }}",
        "VNTYPER_TEST_IMAGE": "vntyper:test",
        "VNTYPER_TEST_DATA_SKIP_DOWNLOAD": "1",
    }
    assert "set -euo pipefail" in str(step["run"])
    assert "continue-on-error" not in step
    assert "|| true" not in str(step["run"])
    assert names.index("Build application image") < names.index("Final test data verification")
    assert names.index("Final test data verification") < names.index("Run Docker tests")
