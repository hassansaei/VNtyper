"""Parsed and executable contracts for the release workflows."""

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

from scripts.release_policy import REQUIRED_CHECK_NAMES

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]


def _workflow(name: str) -> dict[str, Any]:
    """Load one workflow without YAML 1.1 coercing the ``on`` key."""
    return yaml.load((ROOT / ".github" / "workflows" / name).read_text(), Loader=yaml.BaseLoader)


def _step_with_id(workflow: dict[str, Any], job: str, step_id: str) -> dict[str, Any]:
    return next(step for step in workflow["jobs"][job]["steps"] if step.get("id") == step_id)


def _step_using(workflow: dict[str, Any], job: str, action: str) -> dict[str, Any]:
    return next(step for step in workflow["jobs"][job]["steps"] if step.get("uses") == action)


def _heredoc_python(run: str, invocation: str = "python - <<'PY'") -> str:
    after_invocation = run.split(invocation + "\n", maxsplit=1)[1]
    return after_invocation.split("\nPY", maxsplit=1)[0]


def test_main_push_forces_every_substantive_ci_and_docker_scope_true() -> None:
    """A path-only main commit must not skip evidence required by a later tag."""
    ci = _workflow("ci-tests.yml")
    docker = _workflow("docker-build.yml")

    assert ci["concurrency"] == {
        "group": "ci-${{ github.event_name == 'push' && github.sha || github.ref }}",
        "cancel-in-progress": "${{ github.event_name == 'pull_request' }}",
    }
    assert docker["concurrency"]["group"] == "docker-${{ github.event_name == 'push' && github.sha || github.ref }}"
    assert docker["concurrency"]["cancel-in-progress"] == "false"

    for output in ("python", "docs"):
        assert "github.event_name == 'push'" in ci["jobs"]["changes"]["outputs"][output]
        assert f"steps.filter.outputs.{output} == 'true'" in ci["jobs"]["changes"]["outputs"][output]
    assert "github.event_name == 'push'" in docker["jobs"]["changes"]["outputs"]["image"]
    assert "steps.filter.outputs.image == 'true'" in docker["jobs"]["changes"]["outputs"]["image"]

    assert {ci["jobs"][job]["if"] for job in ("lint", "typecheck", "test-unit")} == {
        "needs.changes.outputs.python == 'true'"
    }
    assert ci["jobs"]["docs"]["if"] == "needs.changes.outputs.docs == 'true'"
    assert set(ci["jobs"]["ci-success"]["needs"]) == {"changes", "lint", "typecheck", "test-unit", "docs"}
    assert docker["jobs"]["base-status"]["if"] == "needs.changes.outputs.image == 'true'"
    assert "needs.changes.outputs.image == 'true'" in docker["jobs"]["build-and-test"]["if"]
    assert set(docker["jobs"]["docker-success"]["needs"]) == {
        "changes",
        "base-status",
        "build-base",
        "build-and-test",
    }


def test_pr_filters_permissions_and_fork_missing_base_diagnostic_remain_intact() -> None:
    """Main evidence changes must not broaden fork permissions or remove PR filtering."""
    ci = _workflow("ci-tests.yml")
    docker = _workflow("docker-build.yml")

    assert "pull_request" in ci["on"]
    assert "pull_request" in docker["on"]
    assert "pull_request_target" not in ci["on"]
    assert "pull_request_target" not in docker["on"]
    assert ci["permissions"] == {}
    assert docker["permissions"] == {}

    ci_filters = _step_with_id(ci, "changes", "filter")["with"]["filters"]
    docker_filters = _step_with_id(docker, "changes", "filter")["with"]["filters"]
    assert "vntyper/**" in ci_filters
    assert "docs/**" in ci_filters
    assert "vntyper/**" in docker_filters
    assert "docker/**" in docker_filters

    base_check = _step_with_id(docker, "base-status", "check")["run"]
    assert "github.event.pull_request.head.repo.full_name" in base_check
    assert "a fork cannot push to ghcr.io" in base_check
    assert docker["jobs"]["build-and-test"]["permissions"] == {"contents": "read", "packages": "write"}


def test_workflow_job_names_expand_to_exact_release_check_contract() -> None:
    """Renaming or omitting a substantive job must break exact-SHA release gating."""
    ci = _workflow("ci-tests.yml")
    docker = _workflow("docker-build.yml")

    actual = [
        ci["jobs"]["lint"]["name"],
        ci["jobs"]["typecheck"]["name"],
        *(
            ci["jobs"]["test-unit"]["name"].replace("${{ matrix.python-version }}", version)
            for version in ci["jobs"]["test-unit"]["strategy"]["matrix"]["python-version"]
        ),
        ci["jobs"]["docs"]["name"],
        ci["jobs"]["ci-success"]["name"],
        docker["jobs"]["build-and-test"]["name"],
        docker["jobs"]["docker-success"]["name"],
    ]

    assert tuple(actual) == REQUIRED_CHECK_NAMES


def test_docker_build_uses_short_sha_source_and_explicit_oci_identity_labels() -> None:
    """The source tag and full identity labels must survive metadata refactoring."""
    docker = _workflow("docker-build.yml")
    metadata = _step_using(docker, "build-and-test", "docker/metadata-action@v6")
    tags = metadata["with"]["tags"].splitlines()
    build = _step_using(docker, "build-and-test", "docker/build-push-action@v7")
    labels = build["with"]["labels"]
    package = _step_with_id(docker, "build-and-test", "package")

    assert "type=sha" in tags
    assert "type=ref,event=tag" not in tags
    assert "type=ref,event=branch" in tags
    assert "type=ref,event=pr" in tags
    assert "runpy.run_path('vntyper/version.py')" in package["run"]
    assert "org.opencontainers.image.created=" in labels
    assert "org.opencontainers.image.revision=${{ steps.package.outputs.revision }}" in labels
    assert "org.opencontainers.image.version=${{ steps.package.outputs.version }}" in labels
    assert "steps.meta.outputs.labels" not in labels
    assert "revision=%s\\n" in package["run"]
    assert '"$GITHUB_SHA"' in package["run"]


def test_release_evidence_uses_post_push_registry_digest_and_attempt_scoped_artifact() -> None:
    """Evidence must describe the pushed manifest, not a local build-action output."""
    docker = _workflow("docker-build.yml")
    steps = docker["jobs"]["build-and-test"]["steps"]
    push_index = next(index for index, step in enumerate(steps) if step.get("name") == "Push image")
    registry_index = next(index for index, step in enumerate(steps) if step.get("id") == "registry-digest")
    evidence_index = next(index for index, step in enumerate(steps) if step.get("id") == "release-evidence")
    upload_index = next(index for index, step in enumerate(steps) if step.get("uses") == "actions/upload-artifact@v5")

    assert push_index < registry_index < evidence_index < upload_index
    registry = steps[registry_index]
    assert registry["if"] == "${{ github.event_name == 'push' && github.ref == 'refs/heads/main' }}"
    assert "sha-$SHORT_SHA" in registry["run"]
    assert "--format '{{.Manifest.Digest}}'" in registry["run"]
    assert "--format '{{json .Image}}'" in registry["run"]
    assert "org.opencontainers.image.revision" in registry["run"]
    assert "org.opencontainers.image.version" in registry["run"]

    evidence = steps[evidence_index]
    assert evidence["env"] == {
        "DIGEST": "${{ steps.registry-digest.outputs.value }}",
        "REVISION": "${{ github.sha }}",
        "RUN_ATTEMPT": "${{ github.run_attempt }}",
        "RUN_ID": "${{ github.run_id }}",
        "VERSION": "${{ steps.package.outputs.version }}",
    }
    assert '"contract_version": 1' in evidence["run"]
    upload = steps[upload_index]
    assert upload["if"] == "${{ github.event_name == 'push' && github.ref == 'refs/heads/main' }}"
    assert upload["with"] == {
        "name": "docker-release-evidence-${{ github.sha }}-${{ github.run_attempt }}",
        "path": "docker-release-evidence.json",
        "if-no-files-found": "error",
        "retention-days": "90",
    }


def test_release_evidence_serializer_writes_typed_exact_identity_payload(tmp_path: Path) -> None:
    """The embedded serializer must preserve full SHA, digest, IDs, and package version."""
    docker = _workflow("docker-build.yml")
    release_evidence = _step_with_id(docker, "build-and-test", "release-evidence")
    script = _heredoc_python(release_evidence["run"])
    revision = "0123456789abcdef0123456789abcdef01234567"
    digest = "sha256:" + "ab" * 32
    env = os.environ | {
        "DIGEST": digest,
        "REVISION": revision,
        "RUN_ATTEMPT": "3",
        "RUN_ID": "781",
        "VERSION": "2.4.6",
    }

    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert json.loads((tmp_path / "docker-release-evidence.json").read_text()) == {
        "contract_version": 1,
        "sha": revision,
        "digest": digest,
        "run_id": 781,
        "run_attempt": 3,
        "revision": revision,
        "version": "2.4.6",
    }
