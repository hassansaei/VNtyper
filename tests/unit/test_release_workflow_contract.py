"""Parsed and executable contracts for the release workflows."""

import json
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pytest
import yaml

from scripts.release_policy import REQUIRED_CHECK_NAMES

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]

README_RENAMES = {
    "# VNtyper 2.0 - A Pipeline": "# VNtyper 2 - A Pipeline",
    "**VNtyper 2.0** is an advanced": "**VNtyper 2** is an advanced",
    "VNtyper 2.0 uses modern Python packaging": "VNtyper 2 uses modern Python packaging",
    "VNtyper 2.0 offers multiple subcommands": "VNtyper 2 offers multiple subcommands",
    "Docker image for VNtyper 2.0 is provided": "Docker image for VNtyper 2 is provided",
    "VNtyper 2.0 integrates multiple steps": "VNtyper 2 integrates multiple steps",
    "VNtyper 2.0 relies on several tools": "VNtyper 2 relies on several tools",
    "If you use VNtyper 2.0 in your research": "If you use VNtyper 2 in your research",
}


def _workflow(name: str) -> dict[str, Any]:
    """Load one workflow without YAML 1.1 coercing the ``on`` key."""
    return yaml.load((ROOT / ".github" / "workflows" / name).read_text(), Loader=yaml.BaseLoader)


def _step_with_id(workflow: dict[str, Any], job: str, step_id: str) -> dict[str, Any]:
    return next(step for step in workflow["jobs"][job]["steps"] if step.get("id") == step_id)


def _step_using(workflow: dict[str, Any], job: str, action: str) -> dict[str, Any]:
    return next(step for step in workflow["jobs"][job]["steps"] if step.get("uses") == action)


def _run_step(
    tmp_path: Path,
    workflow: dict[str, Any],
    job: str,
    step_id: str,
    env: dict[str, str],
) -> subprocess.CompletedProcess[str]:
    """Execute one checked-in shell step with controlled external observations."""
    step = _step_with_id(workflow, job, step_id)
    return subprocess.run(
        ["bash", "-euo", "pipefail", "-c", step["run"]],
        cwd=tmp_path,
        env={**os.environ, **env},
        text=True,
        capture_output=True,
        check=False,
    )


def _heredoc_python(run: str, invocation: str = "python - <<'PY'") -> str:
    after_invocation = run.split(invocation + "\n", maxsplit=1)[1]
    return after_invocation.split("\nPY", maxsplit=1)[0]


def _heredoc_pythons(run: str) -> list[str]:
    scripts: list[str] = []
    lines = iter(run.splitlines())
    for line in lines:
        if not line.strip().startswith("python") or not line.rstrip().endswith("<<'PY'"):
            continue
        body: list[str] = []
        for body_line in lines:
            if body_line == "PY":
                break
            body.append(body_line)
        scripts.append("\n".join(body))
    return scripts


def _write_executable(path: Path, source: str) -> None:
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)


def _release_runtime(tmp_path: Path, fixture: Mapping[str, object]) -> dict[str, str]:
    """Install deterministic GitHub/registry fakes for executable release-step tests."""
    controller_scripts = tmp_path / "controller" / "scripts"
    controller_scripts.mkdir(parents=True)
    shutil.copy2(ROOT / "scripts" / "release_policy.py", controller_scripts)
    shutil.copy2(ROOT / "scripts" / "release_evidence.py", controller_scripts)
    shutil.copy2(ROOT / "scripts" / "release_registry.py", controller_scripts)
    manifest_script = ROOT / "scripts" / "release_manifest.py"
    if manifest_script.is_file():
        shutil.copy2(manifest_script, controller_scripts)
        local_scripts = tmp_path / "scripts"
        local_scripts.mkdir()
        shutil.copy2(manifest_script, local_scripts)
    fixture_path = tmp_path / "fake-fixture.json"
    fixture_path.write_text(json.dumps(fixture), encoding="utf-8")
    (tmp_path / "mutation.log").write_text("", encoding="utf-8")
    (tmp_path / "create.log").write_text("", encoding="utf-8")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _write_executable(
        fake_bin / "gh",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

fixture_path = pathlib.Path(os.environ["FAKE_FIXTURE"])
fixture = json.loads(fixture_path.read_text(encoding="utf-8"))
counter_path = fixture_path.with_name("fake-counts.json")
counts = json.loads(counter_path.read_text(encoding="utf-8")) if counter_path.is_file() else {}
args = sys.argv[1:]
if args and args[0] == "api":
    endpoint = next(arg for arg in args if arg.startswith("repos/"))
    if "/check-runs" in endpoint:
        key = "checks"
    elif "/actions/workflows/docker-build.yml/runs" in endpoint:
        key = "runs"
    else:
        run_id = endpoint.split("/actions/runs/", 1)[1].split("/", 1)[0]
        key = f"artifacts:{run_id}"
    counts[key] = counts.get(key, 0) + 1
    counter_path.write_text(json.dumps(counts), encoding="utf-8")
    if counts[key] <= fixture.get("failures", {}).get(key, 0):
        raise SystemExit(1)
    print(json.dumps(fixture.get("payloads", {}).get(key, [])))
    raise SystemExit(0)
if args[:2] == ["run", "download"]:
    if "--repo" not in args or args[args.index("--repo") + 1] != os.environ["GITHUB_REPOSITORY"]:
        print("gh run download requires the exact explicit repository outside a checkout", file=sys.stderr)
        raise SystemExit(2)
    run_id = args[2]
    destination = pathlib.Path(args[args.index("-D") + 1])
    evidence = fixture.get("downloads", {}).get(run_id)
    if evidence is None:
        raise SystemExit(1)
    destination.mkdir(parents=True)
    (destination / "docker-release-evidence.json").write_text(json.dumps(evidence), encoding="utf-8")
    raise SystemExit(0)
with pathlib.Path(os.environ["MUTATION_LOG"]).open("a", encoding="utf-8") as handle:
    handle.write("gh " + " ".join(args) + "\\n")
raise SystemExit(97)
""",
    )
    _write_executable(
        fake_bin / "docker",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

fixture = json.loads(pathlib.Path(os.environ["FAKE_FIXTURE"]).read_text(encoding="utf-8"))
args = sys.argv[1:]
if args[:3] == ["buildx", "imagetools", "create"]:
    create_log = pathlib.Path(os.environ["CREATE_LOG"])
    with create_log.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(args) + "\\n")
    if "--prefer-index=false" not in args:
        raise SystemExit(95)
    positional = []
    skip_next = False
    for arg in args[3:]:
        if skip_next:
            skip_next = False
            continue
        if arg == "--tag":
            skip_next = True
            continue
        if not arg.startswith("-"):
            positional.append(arg)
    if "--dry-run" in args:
        if "--tag" in args or len(positional) != 1:
            raise SystemExit(94)
        print(json.dumps({"source": positional[0], "dry_run": True}))
        raise SystemExit(0)
    if "--tag" not in args or len(positional) != 1:
        raise SystemExit(93)
    target = args[args.index("--tag") + 1]
    source = positional[0]
    with pathlib.Path(os.environ["MUTATION_LOG"]).open("a", encoding="utf-8") as handle:
        handle.write("docker " + " ".join(args) + "\\n")
    source_record = fixture.get("docker", {}).get(source)
    if source_record is None:
        raise SystemExit(92)
    if target not in fixture.get("fail_reinspect", []):
        fixture.setdefault("docker", {})[target] = source_record
        pathlib.Path(os.environ["FAKE_FIXTURE"]).write_text(json.dumps(fixture), encoding="utf-8")
    raise SystemExit(0)
if args[:3] != ["buildx", "imagetools", "inspect"]:
    with pathlib.Path(os.environ["MUTATION_LOG"]).open("a", encoding="utf-8") as handle:
        handle.write("docker " + " ".join(args) + "\\n")
    raise SystemExit(97)
counter_path = pathlib.Path(os.environ["FAKE_FIXTURE"]).with_name("fake-counts.json")
counts = json.loads(counter_path.read_text(encoding="utf-8")) if counter_path.is_file() else {}
format_value = args[args.index("--format") + 1]
inspect_key = f"inspect:{args[3]}:{format_value}"
counts[inspect_key] = counts.get(inspect_key, 0) + 1
counter_path.write_text(json.dumps(counts), encoding="utf-8")
failure = fixture.get("inspect_failures", {}).get(args[3])
if failure is not None and counts[inspect_key] <= failure["times"]:
    print(failure["stderr"], file=sys.stderr)
    raise SystemExit(failure.get("exit_code", 1))
record = fixture.get("docker", {}).get(args[3])
if record is None:
    print(f"{args[3]}: manifest unknown", file=sys.stderr)
    raise SystemExit(1)
if format_value == "{{.Manifest.Digest}}":
    print(f"Name:      {args[3]}")
    print("MediaType: application/vnd.docker.distribution.manifest.v2+json")
    print(f"Digest:    {record['digest']}")
elif format_value == "{{json .Manifest}}":
    descriptor = record.get("manifest", {"mediaType": "application/vnd.oci.image.manifest.v1+json",
                                          "digest": record["digest"], "size": 3464})
    print(json.dumps(descriptor, indent=2))
elif format_value == "{{json .Image}}":
    print(json.dumps(record["image"]))
else:
    raise SystemExit(96)
""",
    )
    _write_executable(fake_bin / "sleep", "#!/usr/bin/env bash\nexit 0\n")
    return {
        "COORDINATOR_ROOT": str(tmp_path / "controller"),
        "CREATE_LOG": str(tmp_path / "create.log"),
        "FAKE_FIXTURE": str(fixture_path),
        "GH_TOKEN": "test-token",
        "GITHUB_OUTPUT": str(tmp_path / "github-output"),
        "GITHUB_REPOSITORY": "hassansaei/VNtyper",
        "MUTATION_LOG": str(tmp_path / "mutation.log"),
        "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
        "RELEASE_API_ATTEMPTS": "3",
        "RELEASE_API_RETRY_SECONDS": "0",
        "RELEASE_ALIAS_INSPECT_ATTEMPTS": "3",
        "RELEASE_ALIAS_INSPECT_RETRY_SECONDS": "0",
        "RELEASE_POLL_ATTEMPTS": "1",
        "RELEASE_POLL_INTERVAL_SECONDS": "0",
        "SHA": "a" * 40,
        "SHORT_SHA": "a" * 7,
        "SOURCE_DIGEST": "sha256:" + "b" * 64,
        "VERSION": "2.0.10",
        "IMAGE": "ghcr.io/hassansaei/vntyper",
    }


def _docker_run(
    run_id: int,
    *,
    attempt: int = 1,
    status: str = "completed",
    conclusion: str | None = "success",
    event: str = "push",
    branch: str = "main",
    sha: str = "a" * 40,
) -> dict[str, object]:
    return {
        "id": run_id,
        "run_attempt": attempt,
        "head_sha": sha,
        "head_branch": branch,
        "event": event,
        "status": status,
        "conclusion": conclusion,
        "html_url": f"https://github.test/runs/{run_id}",
    }


def _check_runs(*, conclusion: str = "success", status: str = "completed") -> list[dict[str, object]]:
    return [
        {
            "id": index,
            "name": name,
            "head_sha": "a" * 40,
            "status": status,
            "conclusion": conclusion if status == "completed" else None,
            "html_url": f"https://github.test/checks/{index}",
            "app": {"slug": "github-actions"},
        }
        for index, name in enumerate(REQUIRED_CHECK_NAMES, start=1)
    ]


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
    assert "--format '{{json .Manifest}}'" in registry["run"]
    assert "release_manifest.py digest" in registry["run"]
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
    assert "tempfile.mkstemp" in evidence["run"]
    assert "os.replace" in evidence["run"]
    upload = steps[upload_index]
    assert upload["if"] == "${{ github.event_name == 'push' && github.ref == 'refs/heads/main' }}"
    assert upload["with"] == {
        "name": "docker-release-evidence-${{ github.sha }}-${{ github.run_attempt }}",
        "path": "docker-release-evidence.json",
        "if-no-files-found": "error",
        "retention-days": "90",
    }


def test_every_workflow_manifest_site_uses_descriptor_json_and_one_fail_closed_parser() -> None:
    """All seven registry identity reads must share the canonical descriptor parser."""
    docker_source = (ROOT / ".github" / "workflows" / "docker-build.yml").read_text(encoding="utf-8")
    publish_source = (ROOT / ".github" / "workflows" / "publish-pypi.yml").read_text(encoding="utf-8")
    combined = docker_source + publish_source

    assert "{{.Manifest.Digest}}" not in combined
    assert docker_source.count("{{json .Manifest}}") == 2
    assert publish_source.count("{{json .Manifest}}") == 5
    assert docker_source.count("release_manifest.py digest") == 2
    assert publish_source.count("release_manifest.py digest") == 5


def test_docker_build_digest_steps_execute_descriptor_parser_for_base_and_pushed_image(tmp_path: Path) -> None:
    """Both Docker workflow shell sites must reduce descriptor JSON to one canonical digest."""
    docker = _workflow("docker-build.yml")
    digest = "sha256:" + "b" * 64
    image = "ghcr.io/hassansaei/vntyper"
    base_image = f"{image}-base"
    base_tag = "base-deadbeefdeadbeef"
    fixture = {
        "docker": {
            f"{base_image}:{base_tag}": {"digest": digest, "image": {"config": {"Labels": {}}}},
            f"{image}:sha-aaaaaaa": {
                "digest": digest,
                "image": {
                    "config": {
                        "Labels": {
                            "org.opencontainers.image.revision": "a" * 40,
                            "org.opencontainers.image.version": "2.0.10",
                        }
                    }
                },
            },
        }
    }

    base_root = tmp_path / "base"
    base_env = _release_runtime(base_root, fixture)
    base_output = base_root / "github-output"
    base_run = _step_with_id(docker, "build-and-test", "base")["run"]
    base_run = base_run.replace("${{ steps.img.outputs.base }}", base_image).replace(
        "${{ steps.hash.outputs.tag }}", base_tag
    )
    base_completed = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", base_run],
        cwd=base_root,
        env={**os.environ, **base_env, "GITHUB_OUTPUT": str(base_output)},
        text=True,
        capture_output=True,
        check=False,
    )

    assert base_completed.returncode == 0, base_completed.stderr
    assert base_output.read_text(encoding="utf-8") == f"ref={base_image}@{digest}\n"

    pushed_root = tmp_path / "pushed"
    pushed_env = _release_runtime(pushed_root, fixture)
    pushed_output = pushed_root / "github-output"
    pushed_completed = _run_step(
        pushed_root,
        docker,
        "build-and-test",
        "registry-digest",
        pushed_env
        | {
            "GITHUB_OUTPUT": str(pushed_output),
            "GITHUB_SHA": "a" * 40,
            "IMAGE": image,
            "VERSION": "2.0.10",
        },
    )

    assert pushed_completed.returncode == 0, pushed_completed.stderr
    assert pushed_output.read_text(encoding="utf-8") == f"value={digest}\n"


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


def test_release_trigger_uses_default_branch_repository_dispatch_for_production() -> None:
    """Production policy must come from the default branch, never a tagged historical commit."""
    publish = _workflow("publish-pypi.yml")

    assert publish["name"] == "Publish PyPI and promote GHCR"
    assert set(publish["on"]) == {"repository_dispatch", "workflow_dispatch"}
    assert publish["on"]["repository_dispatch"] == {"types": ["vntyper_release"]}
    assert publish["on"]["workflow_dispatch"]["inputs"] == {
        "tag": {
            "description": "Existing strict vMAJOR.MINOR.PATCH tag to inspect without writes",
            "required": "true",
            "type": "string",
        }
    }
    assert publish["permissions"] == {}
    assert publish["concurrency"] == {
        "group": "release-${{ inputs.tag || github.event.client_payload.tag }}",
        "cancel-in-progress": "false",
    }
    assert set(publish["jobs"]) == {
        "validate-release",
        "wait-for-release-gates",
        "build-package",
        "promote-ghcr",
        "publish-pypi",
        "release-summary",
    }
    assert publish["jobs"]["validate-release"]["permissions"] == {"actions": "read", "contents": "read"}

    commands = "\n".join(step.get("run", "") for job in publish["jobs"].values() for step in job.get("steps", []))
    assert "git tag" not in commands
    assert "git push" not in commands
    assert "gh release create" not in commands


def test_package_build_is_unprivileged_and_uploads_one_exact_artifact() -> None:
    """Package construction must not share publication or registry authority."""
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["build-package"]
    checkout = _step_using(publish, "build-package", "actions/checkout@v7")
    setup_python = _step_using(publish, "build-package", "actions/setup-python@v7")
    package = _step_with_id(publish, "build-package", "package")
    upload = _step_using(publish, "build-package", "actions/upload-artifact@v5")

    assert job["needs"] == ["validate-release", "wait-for-release-gates"]
    assert job["permissions"] == {"contents": "read"}
    assert "environment" not in job
    assert job["outputs"] == {
        "artifact_name": "${{ steps.package.outputs.artifact_name }}",
        "package_summary_json": "${{ steps.package-result.outputs.package_summary_json }}",
    }
    assert checkout["with"] == {
        "ref": "${{ needs.validate-release.outputs.sha }}",
        "persist-credentials": "false",
        "path": "candidate",
    }
    assert setup_python["with"]["python-version"] == "3.12"
    assert package["working-directory"] == "candidate"
    assert "python3 -m venv .release-venv" in package["run"]
    assert ".release-venv/bin/python -m build" in package["run"]
    assert ".release-venv/bin/twine check dist/*" in package["run"]
    assert (
        "python-dist-${{ needs.validate-release.outputs.version }}-${{ github.run_id }}-${{ github.run_attempt }}"
        in package["run"]
    )
    assert upload["with"] == {
        "name": "${{ steps.package.outputs.artifact_name }}",
        "path": "candidate/dist/",
        "if-no-files-found": "error",
        "retention-days": "7",
    }
    assert "secrets." not in json.dumps(job)
    assert "id-token" not in json.dumps(job)
    assert "packages" not in json.dumps(job["permissions"])
    assert publish["jobs"]["promote-ghcr"]["needs"] == [
        "validate-release",
        "wait-for-release-gates",
        "build-package",
    ]


def test_pypi_publish_is_protected_oidc_only_and_rerun_safe() -> None:
    """Only the protected publisher may request OIDC, with no long-lived token path."""
    workflow_path = ROOT / ".github" / "workflows" / "publish-pypi.yml"
    raw = workflow_path.read_text(encoding="utf-8")
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["publish-pypi"]
    download = _step_using(publish, "publish-pypi", "actions/download-artifact@v5")
    publisher = _step_with_id(publish, "publish-pypi", "publish")
    preflight = _step_with_id(publish, "validate-release", "pypi-environment")

    assert job["needs"] == ["validate-release", "wait-for-release-gates", "build-package", "promote-ghcr"]
    assert job["if"] == (
        "${{ github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release' }}"
    )
    assert job["environment"] == {"name": "pypi"}
    assert job["permissions"] == {"id-token": "write"}
    assert [
        name
        for name, candidate in publish["jobs"].items()
        if candidate.get("permissions", {}).get("id-token") == "write"
    ] == ["publish-pypi"]
    assert job["outputs"] == {"publish_summary_json": "${{ steps.result.outputs.publish_summary_json }}"}
    assert download["with"] == {
        "name": "${{ needs.build-package.outputs.artifact_name }}",
        "path": "dist/",
    }
    assert publisher["uses"] == "pypa/gh-action-pypi-publish@dc37677b2e1c63e2034f94d8a5b11f265b73ba33"
    assert publisher["with"] == {"packages-dir": "dist/", "skip-existing": "true"}
    assert "password" not in publisher["with"]
    assert "user" not in publisher["with"]
    assert "# v1.14.2; verified 2026-08-10 against upstream release/v1" in raw
    for forbidden in (
        "PYPI_API_TOKEN",
        "TWINE_PASSWORD",
        "TWINE_USERNAME",
        "DOCKER_USERNAME",
        "DOCKER_PASSWORD",
    ):
        assert forbidden not in raw
    assert "/pending_deployments" not in raw
    assert "review_pending_deployments" not in raw
    assert "--method POST" not in preflight["run"]
    job_source = json.dumps(job)
    assert "actions/checkout" not in job_source
    assert "actions/setup-python" not in job_source
    assert "python -m build" not in job_source
    assert "pip install" not in job_source


def test_package_result_preserves_partial_failure_state(tmp_path: Path) -> None:
    """A failed build still reports the artifact identity and files it produced."""
    publish = _workflow("publish-pypi.yml")
    dist = tmp_path / "candidate" / "dist"
    dist.mkdir(parents=True)
    (dist / "vntyper-2.0.10-py3-none-any.whl").write_bytes(b"partial")
    output = tmp_path / "github-output"

    completed = _run_step(
        tmp_path,
        publish,
        "build-package",
        "package-result",
        {
            "ARTIFACT_NAME": "python-dist-2.0.10-41-2",
            "GITHUB_OUTPUT": str(output),
            "PACKAGE_OUTCOME": "failure",
        },
    )

    assert completed.returncode == 0, completed.stderr
    output_name, encoded = output.read_text(encoding="utf-8").strip().split("=", maxsplit=1)
    assert output_name == "package_summary_json"
    assert json.loads(encoded) == {
        "artifact_name": "python-dist-2.0.10-41-2",
        "files": ["vntyper-2.0.10-py3-none-any.whl"],
        "step_outcome": "failure",
    }


@pytest.mark.parametrize(
    ("existed_before", "publish_outcome", "expected_result"),
    (
        ("false", "success", "published"),
        ("true", "success", "already-existed-skip"),
        ("false", "failure", "failed"),
        ("true", "failure", "failed"),
    ),
)
def test_publish_result_distinguishes_new_existing_and_failed_reruns(
    tmp_path: Path,
    existed_before: str,
    publish_outcome: str,
    expected_result: str,
) -> None:
    """Partial PyPI reruns must remain distinguishable from first publication and failure."""
    publish = _workflow("publish-pypi.yml")
    output = tmp_path / "github-output"

    completed = _run_step(
        tmp_path,
        publish,
        "publish-pypi",
        "result",
        {
            "EXISTED_BEFORE": existed_before,
            "GITHUB_OUTPUT": str(output),
            "PUBLISH_OUTCOME": publish_outcome,
        },
    )

    assert completed.returncode == 0, completed.stderr
    output_name, encoded = output.read_text(encoding="utf-8").strip().split("=", maxsplit=1)
    assert output_name == "publish_summary_json"
    assert json.loads(encoded) == {
        "result": expected_result,
        "step_outcome": publish_outcome,
        "existed_before": existed_before == "true",
    }


def test_current_container_commands_use_a_published_release_alias() -> None:
    """Every runnable install example must name a tag the release workflow publishes.

    The release job promotes a verified digest to `v<version>`, `<version>`, the `X.Y`
    and `X` series tags and `latest`, so those aliases exist and the surfaces must stop
    saying they arrive only after a future release (#214). `main` stays documented as the
    rolling branch tag, but no runnable example may send a user to it.
    """
    surfaces = (
        ROOT / "README.md",
        ROOT / "docker" / "README.md",
        ROOT / "docs" / "getting-started" / "installation.md",
        ROOT / "docs" / "user-guide" / "docker.md",
    )
    active_kinds: set[str] = set()
    unsupported = re.compile(r"(?<![A-Za-z0-9_.-])(?:docker://)?saei/vntyper")
    published = re.compile(r"ghcr\.io/hassansaei/vntyper:(?:latest|v?\d+(?:\.\d+){0,2})(?![A-Za-z0-9._-])")
    for path in surfaces:
        text = path.read_text(encoding="utf-8")
        normalized = " ".join(text.lower().split())
        assert "rolling" in normalized
        # The aliases exist now; claiming otherwise is the #214 defect, inverted.
        assert "unreleased" not in normalized
        assert "first gated release" not in normalized
        assert "`latest`" in text
        assert "`vX.Y.Z`" in text
        assert "`X.Y.Z`" in text
        assert "`main`" in text
        for block in re.findall(r"```bash\n(.*?)```", text, flags=re.DOTALL):
            kinds = {kind for kind in ("docker pull", "docker run", "apptainer pull") if kind in block}
            if not kinds:
                continue
            active_kinds.update(kinds)
            assert published.search(block) is not None
            assert "ghcr.io/hassansaei/vntyper:main" not in block
            assert unsupported.search(block) is None
    assert active_kinds == {"docker pull", "docker run", "apptainer pull"}


def test_generation_renames_grammar_and_protected_identities_are_exact() -> None:
    """Only the nine approved generation-prose targets may drop the minor component."""
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    dockerfile = (ROOT / "docker" / "Dockerfile").read_text(encoding="utf-8")
    for before, after in README_RENAMES.items():
        assert before not in readme
        assert readme.count(after) == 1
    assert "This version is a refactored version of VNtyper v1 integrates" not in readme
    assert "This refactored version of VNtyper v1 integrates" in readme
    assert 'org.opencontainers.image.description="VNtyper 2.0 - MUC1 VNTR genotyping pipeline for ADTKD-MUC1"' not in (
        dockerfile
    )
    assert (
        dockerfile.count(
            'org.opencontainers.image.description="VNtyper 2 - MUC1 VNTR genotyping pipeline for ADTKD-MUC1"'
        )
        == 1
    )

    assert 'title: "VNtyper"' in (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert "site_name: VNtyper" in (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    assert "## VNtyper Version:" in (ROOT / "vntyper" / "scripts" / "kestrel_genotyping.py").read_text(encoding="utf-8")
    assert (ROOT / "snakemake" / "vntyper2.smk").is_file()
    assert (ROOT / "snakemake" / "run_vntyper2.sh").is_file()
    assert "Before VNtyper 2.0.6 all of" in (ROOT / "docs" / "cli" / "online.md").read_text(encoding="utf-8")


def test_validate_release_resolves_in_controller_and_tests_exact_candidate() -> None:
    """An old tag must use current policy while validating metadata from its own commit."""
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["validate-release"]
    steps = job["steps"]
    checkouts = [step for step in steps if step.get("uses") == "actions/checkout@v7"]

    assert job["outputs"] == {
        "mode": "${{ steps.resolve.outputs.mode }}",
        "tag": "${{ steps.resolve.outputs.tag }}",
        "version": "${{ steps.resolve.outputs.version }}",
        "sha": "${{ steps.resolve.outputs.sha }}",
        "short_sha": "${{ steps.resolve.outputs.short_sha }}",
        "summary_json": "${{ steps.validate-result.outputs.summary_json }}",
    }
    assert checkouts == [
        {
            "uses": "actions/checkout@v7",
            "with": {"fetch-depth": "0", "persist-credentials": "false", "path": "controller"},
        },
        {
            "uses": "actions/checkout@v7",
            "with": {
                "ref": "${{ steps.resolve.outputs.sha }}",
                "persist-credentials": "false",
                "path": "candidate",
            },
        },
    ]

    resolve = _step_with_id(publish, "validate-release", "resolve")
    candidate = _step_with_id(publish, "validate-release", "candidate")
    assert resolve["working-directory"] == "controller"
    assert resolve["env"] == {
        "DISPATCH_TAG": "${{ inputs.tag }}",
        "EVENT_ACTION": "${{ github.event.action }}",
        "EVENT_NAME": "${{ github.event_name }}",
        "PRODUCTION_TAG": "${{ github.event.client_payload.tag }}",
    }
    assert candidate["working-directory"] == "candidate"
    assert "from scripts.release_policy import parse_release_tag" in resolve["run"]
    assert resolve["run"].index("parse_release_tag") < resolve["run"].index('git rev-parse "refs/tags/${TAG}^{commit}"')
    assert "git fetch --no-tags origin main" in resolve["run"]
    assert 'git merge-base --is-ancestor "$SHA" origin/main' in resolve["run"]
    assert "MODE=dry-run" in resolve["run"]
    assert "TAG=$DISPATCH_TAG" in resolve["run"]
    assert 'elif [ "$EVENT_NAME" = "repository_dispatch" ] && [ "$EVENT_ACTION" = "vntyper_release" ]' in resolve["run"]
    assert "MODE=production" in resolve["run"]
    assert "TAG=$PRODUCTION_TAG" in resolve["run"]
    assert "PUSH_SHA" not in resolve["run"]
    assert "PUSH_TAG" not in resolve["run"]

    assert 'runpy.run_path("vntyper/version.py")' in candidate["run"]
    assert "python3 -m venv .release-venv" in candidate["run"]
    assert ".release-venv/bin/pip install pytest packaging PyYAML requests" in candidate["run"]
    assert ".release-venv/bin/pytest -m unit tests/unit/test_version_consistency.py -q" in candidate["run"]


def test_dispatch_is_read_only_and_metadata_eligible_before_any_check_polling() -> None:
    """Validation stays read-only while later exact-SHA gates own polling and registry reads."""
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["validate-release"]
    commands = "\n".join(step.get("run", "") for step in job["steps"])

    assert set(publish["jobs"]) == {
        "validate-release",
        "wait-for-release-gates",
        "build-package",
        "promote-ghcr",
        "publish-pypi",
        "release-summary",
    }
    assert job["permissions"] == {"actions": "read", "contents": "read"}
    assert "id-token" not in job["permissions"]
    assert "packages" not in job["permissions"]
    assert "classify_check_runs" not in commands
    assert "/check-runs" not in commands
    assert "--method POST" not in commands


def test_validation_summary_preserves_structured_mismatch_observations(tmp_path: Path) -> None:
    """A failed candidate test must retain each independently observed version and verdict."""
    publish = _workflow("publish-pypi.yml")
    candidate = _step_with_id(publish, "validate-release", "candidate")
    validate_result = _step_with_id(publish, "validate-release", "validate-result")
    candidate_scripts = _heredoc_pythons(candidate["run"])
    assert len(candidate_scripts) == 2

    candidate_root = tmp_path / "candidate"
    (candidate_root / "vntyper").mkdir(parents=True)
    (candidate_root / "docs" / "about").mkdir(parents=True)
    (candidate_root / "vntyper" / "version.py").write_text('__version__ = "7.8.8"\n')
    (candidate_root / "CITATION.cff").write_text('version: "7.8.7"\n')
    (candidate_root / "docs" / "about" / "changelog.md").write_text("## 7.8.6\n")
    candidate_env = os.environ | {"MODE": "dry-run", "SHA": "1" * 40, "TAG": "v7.8.9", "VERSION": "7.8.9"}

    first = subprocess.run(
        [sys.executable, "-c", candidate_scripts[0]],
        cwd=candidate_root,
        env=candidate_env,
        check=False,
        capture_output=True,
        text=True,
    )
    assert first.returncode == 0, first.stderr
    second = subprocess.run(
        [sys.executable, "-c", candidate_scripts[1], "1"],
        cwd=candidate_root,
        env=candidate_env,
        check=False,
        capture_output=True,
        text=True,
    )
    assert second.returncode == 1
    assert "candidate release metadata does not match" in second.stderr

    summary_env = os.environ | {
        "CANDIDATE_OUTCOME": "failure",
        "MODE": "dry-run",
        "RESOLVE_OUTCOME": "success",
        "SHA": "1" * 40,
        "TAG": "v7.8.9",
        "VERSION": "7.8.9",
    }
    summary = subprocess.run(
        [sys.executable, "-c", _heredoc_python(validate_result["run"], "python - <<'PY' >> \"$GITHUB_OUTPUT\"")],
        cwd=tmp_path,
        env=summary_env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert summary.returncode == 0, summary.stderr
    prefix, encoded = summary.stdout.strip().split("=", maxsplit=1)
    assert prefix == "summary_json"
    assert json.loads(encoded) == {
        "mode": "dry-run",
        "tag": "v7.8.9",
        "version": "7.8.9",
        "sha": "1" * 40,
        "main_ancestor": True,
        "resolve_outcome": "success",
        "candidate_outcome": "failure",
        "version_validation": {
            "expected_version": "7.8.9",
            "package": {"observed": "7.8.8", "matches": False},
            "citation": {"observed": "7.8.7", "matches": False},
            "changelog": {"observed": "7.8.6", "matches": False},
            "version_test_exit_code": 1,
            "version_test_passed": False,
        },
    }


@pytest.mark.parametrize("invalid_contents", ["", "{", "[]"])
def test_validation_result_always_emits_structured_fallback_for_invalid_candidate_phase_json(
    tmp_path: Path,
    invalid_contents: str,
) -> None:
    """A torn candidate observation must not suppress the always-run validation output."""
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    (candidate / "candidate-version-observations.json").write_text(invalid_contents, encoding="utf-8")
    output = tmp_path / "github-output"
    publish = _workflow("publish-pypi.yml")

    completed = _run_step(
        tmp_path,
        publish,
        "validate-release",
        "validate-result",
        {
            "CANDIDATE_OUTCOME": "failure",
            "GITHUB_OUTPUT": str(output),
            "MODE": "production",
            "RESOLVE_OUTCOME": "success",
            "SHA": "1" * 40,
            "TAG": "v7.8.9",
            "VERSION": "7.8.9",
        },
    )

    assert completed.returncode == 0, completed.stderr
    name, encoded = output.read_text(encoding="utf-8").strip().split("=", maxsplit=1)
    assert name == "summary_json"
    payload = json.loads(encoded)
    assert payload["candidate_outcome"] == "failure"
    assert payload["version_validation"] == {
        "expected_version": "7.8.9",
        "package": {"observed": None, "matches": False},
        "citation": {"observed": None, "matches": False},
        "changelog": {"observed": None, "matches": False},
        "version_test_exit_code": None,
        "version_test_passed": False,
    }


def test_candidate_phase_json_uses_atomic_replace_for_initial_and_updated_snapshots() -> None:
    """Both candidate writers must publish a complete observation object atomically."""
    candidate = _step_with_id(_workflow("publish-pypi.yml"), "validate-release", "candidate")
    scripts = _heredoc_pythons(candidate["run"])

    assert len(scripts) == 2
    assert all("os.replace" in script for script in scripts)
    assert all("tempfile.mkstemp" in script for script in scripts)


@pytest.mark.parametrize(
    ("package", "citation", "changelog"),
    [
        ("7.8.8", "7.8.9", "7.8.9"),
        ("7.8.9", "7.8.8", "7.8.9"),
        ("7.8.9", "7.8.9", "7.8.8"),
        ("7.8.8", "7.8.8", "7.8.8"),
    ],
)
def test_candidate_step_rejects_each_false_version_verdict_even_when_pytest_passes(
    tmp_path: Path,
    package: str,
    citation: str,
    changelog: str,
) -> None:
    """No tag mismatch may pass merely because candidate metadata is internally consistent."""
    publish = _workflow("publish-pypi.yml")
    candidate = _step_with_id(publish, "validate-release", "candidate")
    candidate_root = tmp_path / "candidate"
    fake_bin = tmp_path / "bin"
    (candidate_root / "vntyper").mkdir(parents=True)
    (candidate_root / "docs" / "about").mkdir(parents=True)
    fake_bin.mkdir()
    (candidate_root / "vntyper" / "version.py").write_text(f'__version__ = "{package}"\n')
    (candidate_root / "CITATION.cff").write_text(f'version: "{citation}"\n')
    (candidate_root / "docs" / "about" / "changelog.md").write_text(f"## {changelog}\n")

    fake_python3 = fake_bin / "python3"
    fake_python3.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'test "$1" = -m\n'
        'test "$2" = venv\n'
        'mkdir -p "$3/bin"\n'
        "printf '#!/usr/bin/env bash\\nexit 0\\n' > \"$3/bin/pip\"\n"
        "printf '#!/usr/bin/env bash\\nexit 0\\n' > \"$3/bin/pytest\"\n"
        'chmod +x "$3/bin/pip" "$3/bin/pytest"\n'
    )
    fake_python3.chmod(0o755)
    env = os.environ | {
        "MODE": "production",
        "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
        "SHA": "1" * 40,
        "TAG": "v7.8.9",
        "VERSION": "7.8.9",
    }

    completed = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", candidate["run"]],
        cwd=candidate_root,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode != 0
    observations = json.loads((candidate_root / "candidate-version-observations.json").read_text())
    assert observations["version_test_exit_code"] == 0
    assert observations["version_test_passed"] is True
    assert False in {
        observations["package"]["matches"],
        observations["citation"]["matches"],
        observations["changelog"]["matches"],
    }


@pytest.mark.parametrize(
    ("event_name", "event_action", "dispatch_tag", "production_tag"),
    [
        ("push", "", "", "v7.8.9"),
        ("repository_dispatch", "not_a_release", "", "v7.8.9"),
        ("repository_dispatch", "vntyper_release", "", "HOSTILE"),
        ("workflow_dispatch", "", "HOSTILE", ""),
    ],
)
def test_resolve_rejects_wrong_events_and_hostile_tags_before_git(
    tmp_path: Path,
    event_name: str,
    event_action: str,
    dispatch_tag: str,
    production_tag: str,
) -> None:
    """Only the named dispatches may reach Git, and tag text must be parsed first."""
    publish = _workflow("publish-pypi.yml")
    resolve = _step_with_id(publish, "validate-release", "resolve")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    git_log = tmp_path / "git-called"
    injection_marker = tmp_path / "injected"
    hostile_tag = f"v7.8.9;touch {injection_marker}"
    dispatch_tag = hostile_tag if dispatch_tag == "HOSTILE" else dispatch_tag
    production_tag = hostile_tag if production_tag == "HOSTILE" else production_tag
    fake_git = fake_bin / "git"
    fake_git.write_text(f"#!/usr/bin/env bash\nprintf called > {git_log}\nexit 97\n")
    fake_git.chmod(0o755)
    env = os.environ | {
        "DISPATCH_TAG": dispatch_tag,
        "EVENT_ACTION": event_action,
        "EVENT_NAME": event_name,
        "GITHUB_OUTPUT": str(tmp_path / "github-output"),
        "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
        "PRODUCTION_TAG": production_tag,
        "PUSH_TAG": production_tag,
    }

    completed = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", resolve["run"]],
        cwd=ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode != 0
    assert not git_log.exists()
    assert not injection_marker.exists()


def test_exact_sha_gate_has_bounded_polling_auth_and_read_only_permissions() -> None:
    """Missing auth, skipped checks, or an unbounded loop must never release a candidate."""
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["wait-for-release-gates"]
    preflight = _step_with_id(publish, "wait-for-release-gates", "evidence-preflight")
    poll = _step_with_id(publish, "wait-for-release-gates", "poll")
    evidence = _step_with_id(publish, "wait-for-release-gates", "evidence")
    dry_aliases = _step_with_id(publish, "wait-for-release-gates", "dry-run-aliases")

    assert job["needs"] == "validate-release"
    poll_attempts = int(job["env"]["RELEASE_POLL_ATTEMPTS"])
    poll_interval = int(job["env"]["RELEASE_POLL_INTERVAL_SECONDS"])
    api_attempts = int(job["env"]["RELEASE_API_ATTEMPTS"])
    api_retry = int(job["env"]["RELEASE_API_RETRY_SECONDS"])
    alias_attempts = int(dry_aliases["env"]["RELEASE_ALIAS_INSPECT_ATTEMPTS"])
    alias_retry = int(dry_aliases["env"]["RELEASE_ALIAS_INSPECT_RETRY_SECONDS"])
    worst_poll_and_retry_seconds = (poll_attempts - 1) * poll_interval + poll_attempts * (api_attempts - 1) * api_retry
    candidate_retry_groups = 5  # preflight runs/artifact plus evidence runs/artifact/download
    candidate_retry_seconds = candidate_retry_groups * (api_attempts - 1) * api_retry
    alias_retry_seconds = 5 * (alias_attempts - 1) * alias_retry
    full_budget_seconds = worst_poll_and_retry_seconds + candidate_retry_seconds + alias_retry_seconds + 10 * 60
    assert job["timeout-minutes"] == "120"
    assert int(job["timeout-minutes"]) * 60 > full_budget_seconds
    assert job["permissions"] == {
        "actions": "read",
        "checks": "read",
        "contents": "read",
        "packages": "read",
    }
    assert job["env"] == {
        "COORDINATOR_ROOT": "${{ github.workspace }}/controller",
        "RELEASE_POLL_ATTEMPTS": "120",
        "RELEASE_POLL_INTERVAL_SECONDS": "30",
        "RELEASE_API_ATTEMPTS": "3",
        "RELEASE_API_RETRY_SECONDS": "5",
    }
    assert preflight["env"]["GH_TOKEN"] == "${{ secrets.GITHUB_TOKEN }}"
    assert poll["env"]["GH_TOKEN"] == "${{ secrets.GITHUB_TOKEN }}"
    assert evidence["env"]["GH_TOKEN"] == "${{ secrets.GITHUB_TOKEN }}"
    assert "while true" not in poll["run"]
    assert 'for attempt in $(seq 1 "$RELEASE_POLL_ATTEMPTS")' in poll["run"]
    assert 'for API_ATTEMPT in $(seq 1 "$RELEASE_API_ATTEMPTS")' in poll["run"]
    assert "max_attempts=int(sys.argv[3])" in poll["run"]
    assert poll["run"].rstrip().endswith("exit 1")
    for name in REQUIRED_CHECK_NAMES:
        assert name in poll["run"]


def test_production_registry_login_precedes_inspection_and_manual_dispatch_cannot_login() -> None:
    """Only the authenticated production event may receive registry credentials."""
    publish = _workflow("publish-pypi.yml")
    steps = publish["jobs"]["wait-for-release-gates"]["steps"]
    login_index = next(index for index, step in enumerate(steps) if step.get("uses") == "docker/login-action@v4")
    evidence_index = next(index for index, step in enumerate(steps) if step.get("id") == "evidence")
    login = steps[login_index]

    assert login_index < evidence_index
    assert (
        login["if"] == "${{ github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release' }}"
    )
    assert login["with"] == {
        "registry": "ghcr.io",
        "username": "${{ github.actor }}",
        "password": "${{ secrets.GITHUB_TOKEN }}",
    }
    assert "docker buildx imagetools inspect" in steps[evidence_index]["run"]
    assert "docker login" not in steps[evidence_index]["run"]


def test_source_recovery_uses_only_exact_sha_attempt_artifacts_and_digest_inspection() -> None:
    """Fallback may change the evidence run, never the candidate SHA or provenance digest."""
    publish = _workflow("publish-pypi.yml")
    preflight = _step_with_id(publish, "wait-for-release-gates", "evidence-preflight")["run"]
    evidence = _step_with_id(publish, "wait-for-release-gates", "evidence")["run"]

    assert "docker-release-evidence-${SHA}-${RUN_ATTEMPT}" in preflight
    assert "docker-release-evidence-${SHA}-${RUN_ATTEMPT}" in evidence
    assert "docker-release-evidence-${SHA}-${CANDIDATE_ATTEMPT}" in evidence
    assert 'gh run download "$RUN_ID" --repo "$GITHUB_REPOSITORY"' in evidence
    assert '"$IMAGE@$DIGEST"' in evidence
    assert '"$IMAGE:sha-${SHORT_SHA}"' in evidence
    assert "sha-${SHORT_SHA}" in evidence
    assert ":main" not in evidence
    assert "rerun this existing Docker Build run" in evidence


def test_evidence_download_requires_exact_repo_outside_checkout_and_workflow_supplies_it(tmp_path: Path) -> None:
    """The selected artifact must remain downloadable from the non-git GitHub workspace."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _successful_evidence_fixture()) | {"GITHUB_WORKSPACE": str(tmp_path)}
    assert not (tmp_path / ".git").exists()

    without_repo = subprocess.run(
        ["gh", "run", "download", "41", "-n", "evidence", "-D", "without-repo"],
        cwd=tmp_path,
        env={**os.environ, **env},
        text=True,
        capture_output=True,
        check=False,
    )
    wrong_repo = subprocess.run(
        ["gh", "run", "download", "41", "--repo", "other/project", "-n", "evidence", "-D", "wrong-repo"],
        cwd=tmp_path,
        env={**os.environ, **env},
        text=True,
        capture_output=True,
        check=False,
    )
    exact_repo = subprocess.run(
        [
            "gh",
            "run",
            "download",
            "41",
            "--repo",
            "hassansaei/VNtyper",
            "-n",
            "evidence",
            "-D",
            "exact-repo",
        ],
        cwd=tmp_path,
        env={**os.environ, **env},
        text=True,
        capture_output=True,
        check=False,
    )
    workflow = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence", env)

    assert without_repo.returncode == 2
    assert wrong_repo.returncode == 2
    assert exact_repo.returncode == 0, exact_repo.stderr
    assert workflow.returncode == 0, workflow.stderr


def test_evidence_provenance_and_recovery_fields_are_job_outputs_for_final_summary_dataflow() -> None:
    """The future final summary must not reconstruct source or recovery identity from logs."""
    publish = _workflow("publish-pypi.yml")

    assert publish["jobs"]["wait-for-release-gates"]["outputs"] == {
        "source_ref": "${{ steps.evidence.outputs.source_ref }}",
        "source_digest": "${{ steps.evidence.outputs.source_digest }}",
        "source_revision": "${{ steps.evidence.outputs.source_revision }}",
        "source_version": "${{ steps.evidence.outputs.source_version }}",
        "evidence_contract_version": "${{ steps.evidence.outputs.evidence_contract_version }}",
        "short_tag_collision": "${{ steps.evidence.outputs.short_tag_collision }}",
        "source_run_id": "${{ steps.evidence.outputs.source_run_id }}",
        "source_run_attempt": "${{ steps.evidence.outputs.source_run_attempt }}",
        "source_run_url": "${{ steps.evidence.outputs.source_run_url }}",
        "source_artifact_name": "${{ steps.evidence.outputs.source_artifact_name }}",
        "source_artifact_url": "${{ steps.evidence.outputs.source_artifact_url }}",
        "source_artifact_download_url": "${{ steps.evidence.outputs.source_artifact_download_url }}",
        "recovery_instruction": "${{ steps.evidence.outputs.recovery_instruction }}",
        "preflight_summary_json": "${{ steps.evidence-preflight-result.outputs.preflight_summary_json }}",
        "check_summary_json": "${{ steps.poll-result.outputs.check_summary_json }}",
        "evidence_summary_json": "${{ steps.evidence-result.outputs.evidence_summary_json }}",
        "dry_run_alias_summary_json": "${{ steps.dry-run-aliases.outputs.dry_run_alias_summary_json }}",
    }


def test_release_phase_json_is_published_only_by_atomic_rename() -> None:
    """A failed producer must never expose a partial JSON file to an always-run serializer."""
    publish = _workflow("publish-pypi.yml")
    phase_files = [
        ("wait-for-release-gates", "evidence-preflight", "preflight-state.json"),
        ("wait-for-release-gates", "evidence-preflight", "eligible-runs.json"),
        ("wait-for-release-gates", "poll", "poll.json"),
        ("wait-for-release-gates", "evidence", "evidence-eligible-runs.json"),
        ("wait-for-release-gates", "evidence", "selected-run.json"),
        ("wait-for-release-gates", "evidence", "selected-artifact.json"),
        ("wait-for-release-gates", "evidence", "validated-evidence.json"),
        ("wait-for-release-gates", "evidence", "image-config.json"),
        ("wait-for-release-gates", "evidence", "short-image-config.json"),
        ("wait-for-release-gates", "evidence", "validated-image.json"),
        ("promote-ghcr", "promote", "plan.json"),
        ("promote-ghcr", "promote", "alias-progress.json"),
    ]

    for job, step_id, filename in phase_files:
        run = _step_with_id(publish, job, step_id)["run"]
        assert f"{filename}.tmp" in run
        assert f"mv -- {filename}.tmp {filename}" in run


def test_preflight_retries_api_and_falls_back_to_older_exact_artifact(tmp_path: Path) -> None:
    """A missing newest artifact must recover through the next eligible run in numeric order."""
    sha = "a" * 40
    fixture: dict[str, object] = {
        "failures": {"runs": 2},
        "payloads": {
            "runs": [{"workflow_runs": [_docker_run(41, attempt=3), _docker_run(42, attempt=2)]}],
            "artifacts:42": [{"artifacts": []}],
            "artifacts:41": [{"artifacts": [{"name": f"docker-release-evidence-{sha}-3", "expired": False}]}],
        },
    }
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence-preflight", env)
    result = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "evidence-preflight-result",
        env | {"PREFLIGHT_OUTCOME": "success"},
    )

    assert completed.returncode == 0, completed.stderr
    assert result.returncode == 0, result.stderr
    summary = json.loads((tmp_path / "preflight-state.json").read_text(encoding="utf-8"))
    assert set(summary) == {"sha", "state", "reason", "eligible_runs", "selected_run", "step_outcome"}
    assert summary["state"] == "eligible"
    assert [run["id"] for run in summary["eligible_runs"]] == [42, 41]
    assert summary["selected_run"] == {
        "id": 41,
        "run_attempt": 3,
        "html_url": "https://github.test/runs/41",
    }
    assert json.loads((tmp_path / "fake-counts.json").read_text())["runs"] == 3
    assert (tmp_path / "mutation.log").read_text() == ""


@pytest.mark.parametrize(
    ("fixture", "expected_state", "reason_fragment"),
    [
        (
            {
                "payloads": {
                    "runs": [{"workflow_runs": [_docker_run(42)]}],
                    "artifacts:42": [{"artifacts": []}],
                }
            },
            "ineligible",
            "none has its exact attempt-qualified",
        ),
        (
            {"failures": {"runs": 3}, "payloads": {}},
            "infrastructure-failure",
            "failed after 3 attempts",
        ),
    ],
)
def test_preflight_terminal_failures_are_structured_and_poll_is_not_run(
    tmp_path: Path,
    fixture: dict[str, object],
    expected_state: str,
    reason_fragment: str,
) -> None:
    """Deterministic ineligibility and API exhaustion need distinct recovery semantics."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    preflight = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence-preflight", env)
    serialized = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "evidence-preflight-result",
        env | {"PREFLIGHT_OUTCOME": "failure"},
    )
    poll = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "poll-result",
        env | {"POLL_OUTCOME": "skipped"},
    )

    assert preflight.returncode != 0
    assert serialized.returncode == 0, serialized.stderr
    assert poll.returncode == 0, poll.stderr
    summary = json.loads((tmp_path / "preflight-state.json").read_text())
    poll_summary = json.loads((tmp_path / "poll.json").read_text())
    assert set(summary) == {"sha", "state", "reason", "eligible_runs", "selected_run", "step_outcome"}
    assert summary["state"] == expected_state
    assert reason_fragment in summary["reason"]
    assert poll_summary["action"] == "not-run"
    assert poll_summary["preflight_state"] == expected_state
    assert poll_summary["preflight_reason"] == summary["reason"]
    assert "infrastructure_error" not in poll_summary


def test_poll_serializer_reserves_missing_snapshot_diagnostic_for_nonterminal_preflight(tmp_path: Path) -> None:
    """An unexplained missing poll after pending preflight must not be mislabeled not-run."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, {"payloads": {}})
    (tmp_path / "preflight-state.json").write_text(
        json.dumps(
            {
                "sha": "a" * 40,
                "state": "pending",
                "reason": "Docker Build is still running",
                "eligible_runs": [],
                "selected_run": None,
                "step_outcome": "success",
            }
        ),
        encoding="utf-8",
    )

    completed = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "poll-result",
        env | {"POLL_OUTCOME": "skipped"},
    )

    assert completed.returncode == 0, completed.stderr
    payload = json.loads((tmp_path / "poll.json").read_text())
    assert payload["action"] == "fail"
    assert payload["infrastructure_error"] == "poll step ended before a Check Runs snapshot"
    assert "preflight_state" not in payload


@pytest.mark.parametrize("invalid_contents", ["", "{", "[]"])
def test_preflight_result_always_replaces_invalid_phase_json_with_structured_failure(
    tmp_path: Path,
    invalid_contents: str,
) -> None:
    """A torn preflight write must still yield one structured job output."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, {"payloads": {}})
    (tmp_path / "preflight-state.json").write_text(invalid_contents, encoding="utf-8")

    completed = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "evidence-preflight-result",
        env | {"PREFLIGHT_OUTCOME": "failure"},
    )

    assert completed.returncode == 0, completed.stderr
    payload = json.loads((tmp_path / "preflight-state.json").read_text(encoding="utf-8"))
    assert payload == {
        "sha": "a" * 40,
        "state": "infrastructure-failure",
        "reason": "preflight ended before structured state was written",
        "eligible_runs": [],
        "selected_run": None,
        "step_outcome": "failure",
    }
    name, encoded = (tmp_path / "github-output").read_text(encoding="utf-8").strip().split("=", maxsplit=1)
    assert name == "preflight_summary_json"
    assert json.loads(encoded) == payload


@pytest.mark.parametrize("invalid_contents", ["{", "[]"])
def test_poll_result_always_replaces_invalid_phase_json_with_structured_failure(
    tmp_path: Path,
    invalid_contents: str,
) -> None:
    """A nonempty malformed poll snapshot must not escape through the job output."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, {"payloads": {}})
    (tmp_path / "poll.json").write_text(invalid_contents, encoding="utf-8")

    completed = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "poll-result",
        env | {"POLL_OUTCOME": "failure"},
    )

    assert completed.returncode == 0, completed.stderr
    payload = json.loads((tmp_path / "poll.json").read_text(encoding="utf-8"))
    assert payload["action"] == "fail"
    assert payload["infrastructure_error"] == "poll step ended before a Check Runs snapshot"
    assert payload["step_outcome"] == "failure"


@pytest.mark.parametrize(
    ("filename", "invalid_contents", "availability_key"),
    [
        ("plan.json", "", "plan_available"),
        ("plan.json", "[]", "plan_available"),
        ("alias-progress.json", "{", "alias_progress_available"),
        ("alias-progress.json", "{}", "alias_progress_available"),
    ],
)
def test_promotion_result_always_serializes_invalid_or_wrong_shape_phase_json(
    tmp_path: Path,
    filename: str,
    invalid_contents: str,
    availability_key: str,
) -> None:
    """Promotion recovery must retain a structured result after a torn phase write."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _promotion_fixture())
    (tmp_path / filename).write_text(invalid_contents, encoding="utf-8")

    completed = _run_step(
        tmp_path,
        publish,
        "promote-ghcr",
        "promote-result",
        env | {"PROMOTE_OUTCOME": "failure"},
    )

    assert completed.returncode == 0, completed.stderr
    payload = json.loads((tmp_path / "alias-summary.json").read_text(encoding="utf-8"))
    assert payload["step_outcome"] == "failure"
    assert payload[availability_key] is False
    assert isinstance(payload["plan"], dict)
    assert isinstance(payload["alias_progress"], list)


def test_ineligible_preflight_prints_exact_run_recovery_context(tmp_path: Path) -> None:
    """Operators need the exact run URL and rerun action, not a generic rebuild suggestion."""
    fixture = {
        "payloads": {
            "runs": [{"workflow_runs": [_docker_run(42)]}],
            "artifacts:42": [{"artifacts": []}],
        }
    }
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence-preflight", env)

    assert completed.returncode != 0
    assert "Docker Build" in completed.stderr
    assert "a" * 40 in completed.stderr
    assert "https://github.test/runs/42" in completed.stderr
    assert "rerun this existing Docker Build run" in completed.stderr


@pytest.mark.parametrize(
    ("checks", "failures", "expected_action", "expected_returncode", "expected_calls"),
    [
        (_check_runs(), 2, "success", 0, 3),
        (_check_runs(conclusion="skipped"), 0, "fail", 1, 1),
        (_check_runs(status="in_progress"), 0, "timeout", 1, 1),
        (_check_runs(), 3, "fail", 1, 3),
    ],
)
def test_poll_executes_retry_skip_timeout_and_api_exhaustion_paths(
    tmp_path: Path,
    checks: list[dict[str, object]],
    failures: int,
    expected_action: str,
    expected_returncode: int,
    expected_calls: int,
) -> None:
    """The actual shell must fail closed for skipped, timed-out, and unavailable checks."""
    fixture = {"failures": {"checks": failures}, "payloads": {"checks": [{"check_runs": checks}]}}
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "poll", env)

    assert completed.returncode == expected_returncode
    payload = json.loads((tmp_path / "poll.json").read_text())
    assert payload["action"] == expected_action
    if failures == 3:
        assert payload["infrastructure_error"] == "Check Runs API failed after 3 attempts"
    assert json.loads((tmp_path / "fake-counts.json").read_text())["checks"] == expected_calls
    assert (tmp_path / "mutation.log").read_text() == ""


def _successful_evidence_fixture(**evidence_overrides: object) -> dict[str, Any]:
    sha = "a" * 40
    digest = "sha256:" + "b" * 64
    evidence = {
        "contract_version": 1,
        "sha": sha,
        "digest": digest,
        "run_id": 41,
        "run_attempt": 3,
        "revision": sha,
        "version": "2.0.10",
    }
    evidence.update(evidence_overrides)
    image = {
        "config": {
            "Labels": {
                "org.opencontainers.image.revision": sha,
                "org.opencontainers.image.version": "2.0.10",
            }
        }
    }
    artifact_name = f"docker-release-evidence-{sha}-3"
    return {
        "payloads": {
            "runs": [{"workflow_runs": [_docker_run(41, attempt=3), _docker_run(42, attempt=2)]}],
            "artifacts:42": [{"artifacts": []}],
            "artifacts:41": [
                {
                    "artifacts": [
                        {
                            "id": 501,
                            "name": artifact_name,
                            "expired": False,
                            "url": "https://api.github.test/artifacts/501",
                            "archive_download_url": "https://api.github.test/artifacts/501/zip",
                        }
                    ]
                }
            ],
        },
        "downloads": {"41": evidence},
        "docker": {
            f"ghcr.io/hassansaei/vntyper@{digest}": {"digest": digest, "image": image},
            "ghcr.io/hassansaei/vntyper:sha-aaaaaaa": {"digest": digest, "image": image},
        },
    }


def test_evidence_step_selects_fallback_run_and_exports_exact_digest_identity(tmp_path: Path) -> None:
    """The executable workflow must bind outputs to the older run that owns exact evidence."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _successful_evidence_fixture())

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence", env)
    serialized = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "evidence-result",
        env | {"EVIDENCE_OUTCOME": "success"},
    )

    assert completed.returncode == 0, completed.stderr
    assert serialized.returncode == 0, serialized.stderr
    assert json.loads((tmp_path / "selected-run.json").read_text()) == {
        "id": 41,
        "run_attempt": 3,
        "html_url": "https://github.test/runs/41",
    }
    assert json.loads((tmp_path / "selected-artifact.json").read_text()) == {
        "id": 501,
        "name": f"docker-release-evidence-{'a' * 40}-3",
        "expired": False,
        "url": "https://api.github.test/artifacts/501",
        "archive_download_url": "https://api.github.test/artifacts/501/zip",
    }
    output = (tmp_path / "github-output").read_text().splitlines()
    digest = "sha256:" + "b" * 64
    assert output[:-1] == [
        f"source_ref=ghcr.io/hassansaei/vntyper@{digest}",
        f"source_digest={digest}",
        f"source_revision={'a' * 40}",
        "source_version=2.0.10",
        "evidence_contract_version=1",
        "short_tag_collision=false",
        "source_run_id=41",
        "source_run_attempt=3",
        "source_run_url=https://github.test/runs/41",
        f"source_artifact_name=docker-release-evidence-{'a' * 40}-3",
        "source_artifact_url=https://api.github.test/artifacts/501",
        "source_artifact_download_url=https://api.github.test/artifacts/501/zip",
        (
            "recovery_instruction=Rerun existing Docker Build run https://github.test/runs/41 "
            f"(run 41, attempt 3) to regenerate exact artifact docker-release-evidence-{'a' * 40}-3; "
            "artifact API URL: https://api.github.test/artifacts/501."
        ),
    ]
    assert output[-1].startswith("evidence_summary_json=")
    assert json.loads((tmp_path / "validated-image.json").read_text()) == {
        "source_ref": f"ghcr.io/hassansaei/vntyper@{digest}",
        "source_digest": digest,
        "source_revision": "a" * 40,
        "source_version": "2.0.10",
        "evidence_contract_version": 1,
        "short_tag_collision": False,
    }
    summary = json.loads((tmp_path / "evidence-summary.json").read_text(encoding="utf-8"))
    assert summary["state"] == "verified"
    assert summary["source_ref"] == f"ghcr.io/hassansaei/vntyper@{digest}"
    assert summary["evidence_contract_version"] == 1
    assert summary["short_tag_collision"] is False
    assert summary["selected_run"]["html_url"] == "https://github.test/runs/41"
    assert summary["selected_artifact"]["url"] == "https://api.github.test/artifacts/501"
    assert summary["selected_artifact"]["archive_download_url"] == "https://api.github.test/artifacts/501/zip"
    assert summary["recovery_instruction"].startswith("Rerun existing Docker Build run https://github.test/runs/41")
    assert (tmp_path / "mutation.log").read_text() == ""


@pytest.mark.parametrize(
    ("overrides", "error_fragment"),
    [
        ({"contract_version": 2}, "contract_version"),
        ({"revision": "c" * 40}, "revision"),
        ({"version": "2.0.9"}, "version"),
    ],
)
def test_evidence_step_rejects_wrong_contract_digest_revision_and_version(
    tmp_path: Path,
    overrides: dict[str, object],
    error_fragment: str,
) -> None:
    """Downloaded JSON cannot redirect release provenance away from validated identity."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _successful_evidence_fixture(**overrides))

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence", env)

    assert completed.returncode != 0
    assert error_fragment in completed.stderr
    assert not (tmp_path / "validated-evidence.json").exists()
    assert (tmp_path / "mutation.log").read_text() == ""


def test_evidence_step_rejects_registry_digest_mismatch(tmp_path: Path) -> None:
    """A registry manifest digest cannot differ from the evidence provenance digest."""
    fixture = _successful_evidence_fixture()
    immutable_ref = f"ghcr.io/hassansaei/vntyper@{'sha256:' + 'b' * 64}"
    fixture["docker"][immutable_ref]["digest"] = "sha256:" + "c" * 64
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence", env)

    assert completed.returncode != 0
    assert "manifest digest does not match evidence digest" in completed.stderr
    assert "https://github.test/runs/41" in completed.stderr
    assert (tmp_path / "mutation.log").read_text() == ""


def test_evidence_step_reports_short_prefix_collision_without_substitution(tmp_path: Path) -> None:
    """A colliding short tag is allowed only after immutable evidence has been validated."""
    fixture = _successful_evidence_fixture()
    other_sha = "a" * 7 + "d" * 33
    short_ref = "ghcr.io/hassansaei/vntyper:sha-aaaaaaa"
    fixture["docker"][short_ref] = {
        "digest": "sha256:" + "d" * 64,
        "image": {
            "config": {
                "Labels": {
                    "org.opencontainers.image.revision": other_sha,
                    "org.opencontainers.image.version": "2.0.9",
                }
            }
        },
    }
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence", env)

    assert completed.returncode == 0, completed.stderr
    assert "short_tag_collision=true" in (tmp_path / "github-output").read_text()
    assert f"source_digest={'sha256:' + 'b' * 64}" in (tmp_path / "github-output").read_text()


def test_evidence_step_fails_when_immutable_manifest_is_missing_and_preserves_run_url(tmp_path: Path) -> None:
    """Missing source images must identify the existing Docker Build run to rerun."""
    fixture = _successful_evidence_fixture()
    fixture["docker"] = {}
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence", env)
    summary = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "evidence-result",
        env | {"EVIDENCE_OUTCOME": "failure"},
    )

    assert completed.returncode != 0
    assert "https://github.test/runs/41" in completed.stderr
    assert "rerun this existing Docker Build run" in completed.stderr
    assert summary.returncode == 0, summary.stderr
    payload = json.loads((tmp_path / "evidence-summary.json").read_text())
    assert payload["state"] == "failed"
    assert payload["selected_run"]["html_url"] == "https://github.test/runs/41"
    assert payload["selected_artifact"] == {
        "available": True,
        "id": 501,
        "name": f"docker-release-evidence-{'a' * 40}-3",
        "url": "https://api.github.test/artifacts/501",
        "archive_download_url": "https://api.github.test/artifacts/501/zip",
    }
    assert payload["evidence"]["digest"] == "sha256:" + "b" * 64
    assert payload["evidence_contract_version"] == 1
    assert payload["source_ref"] is None
    assert payload["short_tag_collision"] is None
    assert payload["recovery_instruction"].startswith("Rerun existing Docker Build run https://github.test/runs/41")
    assert payload["image"]["available"] is False


def test_evidence_step_without_selected_artifact_names_every_candidate_run_url(tmp_path: Path) -> None:
    """Artifact absence must retain all exact candidate URLs in logs and structured recovery data."""
    fixture = _successful_evidence_fixture()
    fixture["payloads"]["artifacts:41"] = [{"artifacts": []}]
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "evidence", env)
    summary = _run_step(
        tmp_path,
        publish,
        "wait-for-release-gates",
        "evidence-result",
        env | {"EVIDENCE_OUTCOME": "failure"},
    )

    assert completed.returncode != 0
    assert "https://github.test/runs/42" in completed.stderr
    assert "https://github.test/runs/41" in completed.stderr
    assert "run URL unavailable" not in completed.stderr
    assert summary.returncode == 0, summary.stderr
    payload = json.loads((tmp_path / "evidence-summary.json").read_text(encoding="utf-8"))
    assert [run["html_url"] for run in payload["candidate_runs"]] == [
        "https://github.test/runs/42",
        "https://github.test/runs/41",
    ]
    assert payload["selected_run"]["available"] is False
    assert payload["selected_artifact"]["available"] is False
    assert payload["recovery_instruction"] == (
        "Rerun one of these existing exact-SHA Docker Build runs to regenerate its exact attempt-qualified "
        "evidence artifact: https://github.test/runs/42, https://github.test/runs/41."
    )


def _promotion_fixture(
    aliases: Mapping[str, tuple[str, str | None]] | None = None,
    *,
    fail_reinspect: tuple[str, ...] = (),
    inspect_failures: Mapping[str, Mapping[str, object]] | None = None,
) -> dict[str, object]:
    """Return registry observations for executable promotion-step tests."""
    image = "ghcr.io/hassansaei/vntyper"
    source_digest = "sha256:" + "b" * 64
    source_record = {
        "digest": source_digest,
        "image": {
            "config": {
                "Labels": {
                    "org.opencontainers.image.revision": "a" * 40,
                    "org.opencontainers.image.version": "2.0.10",
                }
            }
        },
    }
    docker: dict[str, object] = {f"{image}@{source_digest}": source_record}
    for alias, (digest, version) in (aliases or {}).items():
        docker[f"{image}:{alias}"] = {
            "digest": digest,
            "image": {"config": {"Labels": {"org.opencontainers.image.version": version}}},
        }
    return {
        "docker": docker,
        "fail_reinspect": [f"{image}:{alias}" for alias in fail_reinspect],
        "inspect_failures": dict(inspect_failures or {}),
    }


def _create_records(tmp_path: Path) -> list[list[str]]:
    return [json.loads(line) for line in (tmp_path / "create.log").read_text().splitlines()]


def test_ghcr_promotion_has_fixed_global_lock_exact_guard_and_digest_only_source() -> None:
    """A release write must acquire one repository lock and copy only the verified digest."""
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["promote-ghcr"]
    promote = _step_with_id(publish, "promote-ghcr", "promote")
    header = (ROOT / ".github" / "workflows" / "publish-pypi.yml").read_text(encoding="utf-8")

    assert job["needs"] == ["validate-release", "wait-for-release-gates", "build-package"]
    assert job["if"] == "${{ github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release' }}"
    assert job["concurrency"] == {
        "group": "vntyper-ghcr-release-promotion",
        "cancel-in-progress": "false",
    }
    assert job["permissions"] == {"contents": "read", "packages": "write"}
    assert "mutual-exclusion lock" in header
    assert "not an unbounded queue" in header
    assert "superseded pending" in header
    assert "must be rerun" in header
    assert "plan_alias_updates" in promote["run"]
    assert 'for ALIAS in "v$VERSION" "$VERSION" "${VERSION%.*}" "${VERSION%%.*}" latest' in promote["run"]
    assert 'docker buildx imagetools create --prefer-index=false --tag "$IMAGE:$ALIAS"' in promote["run"]
    assert '"$IMAGE@$SOURCE_DIGEST"' in promote["run"]
    assert '"$IMAGE:main"' not in promote["run"]
    assert promote["run"].count("release_registry.py classify-absence") == 1
    assert "grep -Ei" not in promote["run"]


def test_promotion_executes_five_digest_copies_and_rerun_converges_to_noops(tmp_path: Path) -> None:
    """A complete rerun must observe the five prior writes and perform no new mutations."""
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _promotion_fixture())
    source = f"ghcr.io/hassansaei/vntyper@{'sha256:' + 'b' * 64}"

    first = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert first.returncode == 0, first.stderr
    records = _create_records(tmp_path)
    assert len(records) == 5
    assert [record[record.index("--tag") + 1].rsplit(":", maxsplit=1)[1] for record in records] == [
        "v2.0.10",
        "2.0.10",
        "2.0",
        "2",
        "latest",
    ]
    assert all(record[-1] == source for record in records)
    progress = json.loads((tmp_path / "alias-progress.json").read_text())
    assert len(progress) == 5
    assert all(row["attempted"] and row["write_succeeded"] and row["verified"] for row in progress)
    assert {row["final_digest"] for row in progress} == {"sha256:" + "b" * 64}

    second = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert second.returncode == 0, second.stderr
    assert len(_create_records(tmp_path)) == 5
    plan = json.loads((tmp_path / "plan.json").read_text())
    assert [item["decision"] for item in plan["plan"]] == ["no-op"] * 5
    assert (tmp_path / "alias-progress.json").read_text() == "[]"


@pytest.mark.parametrize("conflicting_alias", ["v2.0.10", "latest"])
def test_equal_version_digest_conflict_fails_before_any_alias_write(
    tmp_path: Path,
    conflicting_alias: str,
) -> None:
    """Conflicting immutable or floating identity must abort the whole plan before mutation."""
    other_digest = "sha256:" + "c" * 64
    source_digest = "sha256:" + "b" * 64
    aliases = {
        "v2.0.10": (source_digest, "2.0.10"),
        "2.0.10": (source_digest, "2.0.10"),
        "2.0": (source_digest, "2.0.10"),
        "2": (source_digest, "2.0.10"),
        "latest": (source_digest, "2.0.10"),
    }
    aliases[conflicting_alias] = (other_digest, "2.0.10")
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _promotion_fixture(aliases))

    completed = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert completed.returncode != 0
    assert "release alias conflict" in completed.stderr
    assert (tmp_path / "mutation.log").read_text() == ""
    assert _create_records(tmp_path) == []


@pytest.mark.parametrize(
    ("alias", "existing_version", "failure_stderr", "expected_returncode"),
    [
        ("v2.0.10", "2.0.10", "request timed out", 1),
        ("latest", "2.1.0", "unexpected status from HEAD request: 503 Service Unavailable", 0),
    ],
)
def test_transient_inspection_cannot_hide_and_overwrite_existing_alias(
    tmp_path: Path,
    alias: str,
    existing_version: str,
    failure_stderr: str,
    expected_returncode: int,
) -> None:
    """A retry must reveal immutable conflicts and newer floating aliases before planning writes."""
    image = "ghcr.io/hassansaei/vntyper"
    source_digest = "sha256:" + "b" * 64
    existing_digest = "sha256:" + "c" * 64
    other_aliases = {
        candidate: (source_digest, "2.0.10")
        for candidate in ("v2.0.10", "2.0.10", "2.0", "2", "latest")
        if candidate != alias
    }
    other_aliases[alias] = (existing_digest, existing_version)
    fixture = _promotion_fixture(
        other_aliases,
        inspect_failures={f"{image}:{alias}": {"times": 1, "stderr": failure_stderr}},
    )
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert completed.returncode == expected_returncode
    if alias == "v2.0.10":
        assert "release alias conflict" in completed.stderr
    else:
        assert "newer than candidate" in completed.stderr
    assert _create_records(tmp_path) == []
    registry = json.loads((tmp_path / "fake-fixture.json").read_text())["docker"]
    assert registry[f"{image}:{alias}"]["digest"] == existing_digest


@pytest.mark.parametrize(
    "failure_stderr",
    [
        "unauthorized: authentication required",
        "unexpected status from HEAD request: 401 Unauthorized: manifest not found",
        "dial tcp: network is unreachable",
        "request timed out",
        "unexpected status from HEAD request: 503 Service Unavailable",
        "unexpected status from HEAD request: 503 Service Unavailable: repository not found",
        "unexpected status from credentials endpoint: 404 Not Found",
    ],
)
def test_ambiguous_alias_inspection_failure_aborts_before_every_registry_write(
    tmp_path: Path,
    failure_stderr: str,
) -> None:
    """Only an authoritative registry not-found response may be treated as an absent alias."""
    image = "ghcr.io/hassansaei/vntyper"
    hidden_alias = f"{image}:latest"
    fixture = _promotion_fixture(
        inspect_failures={hidden_alias: {"times": 3, "stderr": failure_stderr}},
    )
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert completed.returncode != 0
    assert "alias inspection failed without authoritative not-found" in completed.stderr
    assert _create_records(tmp_path) == []
    assert (tmp_path / "mutation.log").read_text() == ""


def test_authoritative_registry_404_may_plan_create_after_bounded_retries(tmp_path: Path) -> None:
    """An explicit registry 404 remains the only failed inspection that proves alias absence."""
    image = "ghcr.io/hassansaei/vntyper"
    fixture = _promotion_fixture(
        inspect_failures={
            f"{image}:latest": {
                "times": 3,
                "stderr": (
                    "unexpected status from HEAD request to "
                    "https://ghcr.io/v2/hassansaei/vntyper/manifests/latest: 404 Not Found"
                ),
            }
        },
    )
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert completed.returncode == 0, completed.stderr
    assert len(_create_records(tmp_path)) == 5


@pytest.mark.parametrize("completed_prefix", range(6))
def test_promotion_converges_after_every_completed_alias_prefix(tmp_path: Path, completed_prefix: int) -> None:
    """Each partial prefix must safely converge without rewriting completed aliases."""
    digest = "sha256:" + "b" * 64
    ordered = ("v2.0.10", "2.0.10", "2.0", "2", "latest")
    aliases = dict.fromkeys(ordered[:completed_prefix], (digest, "2.0.10"))
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _promotion_fixture(aliases))

    completed = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert completed.returncode == 0, completed.stderr
    records = _create_records(tmp_path)
    assert [record[record.index("--tag") + 1].rsplit(":", maxsplit=1)[1] for record in records] == list(
        ordered[completed_prefix:]
    )
    registry = json.loads((tmp_path / "fake-fixture.json").read_text())["docker"]
    assert all(registry[f"ghcr.io/hassansaei/vntyper:{alias}"]["digest"] == digest for alias in ordered)


def test_newer_floating_alias_emits_a_notice_without_writes(tmp_path: Path) -> None:
    """An anti-downgrade decision must be observable without moving protected aliases."""
    digest = "sha256:" + "b" * 64
    aliases = {
        "v2.0.10": (digest, "2.0.10"),
        "2.0.10": (digest, "2.0.10"),
        "2.0": (digest, "2.0.10"),
        "2": (digest, "2.0.10"),
        "latest": (digest, "2.1.0"),
    }
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _promotion_fixture(aliases))

    completed = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)

    assert completed.returncode == 0, completed.stderr
    assert "::notice title=GHCR alias latest::" in completed.stderr
    assert "newer than candidate" in completed.stderr
    assert (tmp_path / "mutation.log").read_text() == ""


def test_failed_reinspection_preserves_write_succeeded_without_verified(tmp_path: Path) -> None:
    """A completed registry write must remain visible when post-write verification fails."""
    digest = "sha256:" + "b" * 64
    aliases = {
        "v2.0.10": (digest, "2.0.10"),
        "2.0.10": (digest, "2.0.10"),
        "2": (digest, "2.0.10"),
        "latest": (digest, "2.0.10"),
    }
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, _promotion_fixture(aliases, fail_reinspect=("2.0",)))

    completed = _run_step(tmp_path, publish, "promote-ghcr", "promote", env)
    summary = _run_step(
        tmp_path,
        publish,
        "promote-ghcr",
        "promote-result",
        env | {"PROMOTE_OUTCOME": "failure"},
    )

    assert completed.returncode != 0
    assert summary.returncode == 0, summary.stderr
    progress = json.loads((tmp_path / "alias-progress.json").read_text())
    assert progress == [
        {
            "alias": "2.0",
            "decision": "create",
            "reason": "Alias is absent and will be created.",
            "previous_digest": None,
            "previous_version": None,
            "attempted": True,
            "write_succeeded": True,
            "verified": False,
            "final_digest": None,
        }
    ]
    output_name, encoded = (tmp_path / "github-output").read_text().strip().split("=", maxsplit=1)
    assert output_name == "alias_summary_json"
    assert json.loads(encoded)["alias_progress"] == progress


def test_manual_alias_dry_run_executes_one_untagged_probe_and_zero_writes(tmp_path: Path) -> None:
    """Manual mode may exercise Buildx carbon-copy validation but cannot tag any alias."""
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["wait-for-release-gates"]
    step = _step_with_id(publish, "wait-for-release-gates", "dry-run-aliases")
    env = _release_runtime(tmp_path, _promotion_fixture())

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "dry-run-aliases", env)

    assert step["if"] == "${{ github.event_name == 'workflow_dispatch' }}"
    assert job["outputs"]["dry_run_alias_summary_json"] == (
        "${{ steps.dry-run-aliases.outputs.dry_run_alias_summary_json }}"
    )
    assert step["run"].count("release_registry.py classify-absence") == 1
    assert "grep -Ei" not in step["run"]
    assert completed.returncode == 0, completed.stderr
    records = _create_records(tmp_path)
    assert len(records) == 1
    assert "--dry-run" in records[0]
    assert "--prefer-index=false" in records[0]
    assert "--tag" not in records[0]
    assert records[0][-1] == f"ghcr.io/hassansaei/vntyper@{'sha256:' + 'b' * 64}"
    assert (tmp_path / "mutation.log").read_text() == ""
    output_name, encoded = (tmp_path / "github-output").read_text().strip().split("=", maxsplit=1)
    assert output_name == "dry_run_alias_summary_json"
    summary = json.loads(encoded)
    assert [item["alias"] for item in summary["plan"]] == ["v2.0.10", "2.0.10", "2.0", "2", "latest"]
    assert not any(item["execute"] for item in summary["plan"])


def test_manual_alias_dry_run_fails_when_any_alias_conflicts(tmp_path: Path) -> None:
    """Read-only mode must not report success for a plan production would reject."""
    source_digest = "sha256:" + "b" * 64
    conflicting_digest = "sha256:" + "c" * 64
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(
        tmp_path,
        _promotion_fixture({"latest": (conflicting_digest, "2.0.10")}),
    )

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "dry-run-aliases", env)

    assert completed.returncode != 0
    assert "dry-run alias conflict" in completed.stderr
    assert "latest" in completed.stderr
    assert (tmp_path / "mutation.log").read_text() == ""
    assert not (tmp_path / "github-output").exists()
    records = _create_records(tmp_path)
    assert len(records) == 1
    assert "--dry-run" in records[0]
    assert "--tag" not in records[0]
    assert records[0][-1] == f"ghcr.io/hassansaei/vntyper@{source_digest}"


def test_manual_alias_inspection_aborts_on_non_authoritative_failure_before_dry_run_probe(tmp_path: Path) -> None:
    """Dry-run planning must not describe an uncertain registry observation as absence."""
    image = "ghcr.io/hassansaei/vntyper"
    fixture = _promotion_fixture(
        inspect_failures={
            f"{image}:latest": {
                "times": 3,
                "stderr": "unauthorized: authentication required",
            }
        }
    )
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, "wait-for-release-gates", "dry-run-aliases", env)

    assert completed.returncode != 0
    assert "alias inspection failed without authoritative not-found" in completed.stderr
    assert _create_records(tmp_path) == []
    assert not (tmp_path / "github-output").exists()


@pytest.mark.parametrize(
    ("job", "step"),
    [("promote-ghcr", "promote"), ("wait-for-release-gates", "dry-run-aliases")],
)
def test_missing_credential_executable_aborts_both_modes_before_create_or_probe(
    tmp_path: Path,
    job: str,
    step: str,
) -> None:
    """A local credential failure containing ``not found`` never proves manifest absence."""
    image = "ghcr.io/hassansaei/vntyper"
    fixture = _promotion_fixture(
        inspect_failures={
            f"{image}:latest": {
                "times": 3,
                "stderr": (
                    'error getting credentials - err: exec: "docker-credential-gh": '
                    'executable file not found in $PATH, out: ""'
                ),
            }
        }
    )
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, job, step, env)

    assert completed.returncode != 0
    assert "alias inspection failed without authoritative not-found" in completed.stderr
    assert _create_records(tmp_path) == []
    assert (tmp_path / "mutation.log").read_text() == ""
    assert not (tmp_path / "github-output").exists()


@pytest.mark.parametrize(
    "failure_stderr",
    (
        (
            "HEAD https://ghcr.io/v2/hassansaei/vntyper/manifests/latest failed\n"
            "unexpected status from token request: 404 Not Found"
        ),
        (
            "unexpected status from HEAD request to "
            "https://ghcr.io/v2/hassansaei/vntyper/manifests/unrelated: 404 Not Found"
        ),
        "manifest unknown: not found",
        "ghcr.io/hassansaei/vntyper:latest:\nmanifest unknown",
    ),
)
@pytest.mark.parametrize(
    ("job", "step"),
    (("promote-ghcr", "promote"), ("wait-for-release-gates", "dry-run-aliases")),
)
def test_hostile_unbound_absence_output_retries_then_aborts_before_create_or_probe(
    tmp_path: Path,
    failure_stderr: str,
    job: str,
    step: str,
) -> None:
    """Cross-line, wrong-manifest, and unbound statuses cannot authorize a registry operation."""
    image = "ghcr.io/hassansaei/vntyper"
    reference = f"{image}:latest"
    fixture = _promotion_fixture(
        inspect_failures={reference: {"times": 3, "stderr": failure_stderr}},
    )
    publish = _workflow("publish-pypi.yml")
    env = _release_runtime(tmp_path, fixture)

    completed = _run_step(tmp_path, publish, job, step, env)

    assert completed.returncode != 0
    assert "alias inspection failed without authoritative not-found" in completed.stderr
    assert "after 3 attempts" in completed.stderr
    assert _create_records(tmp_path) == []
    assert (tmp_path / "mutation.log").read_text() == ""
    assert not (tmp_path / "github-output").exists()
    counts = json.loads((tmp_path / "fake-counts.json").read_text(encoding="utf-8"))
    assert counts[f"inspect:{reference}:{{{{json .Manifest}}}}"] == 3


def test_release_guidance_records_current_controller_rollout_and_recovery_contract() -> None:
    """Maintainer guidance must describe the implemented default-branch release controller exactly."""
    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    followups = (ROOT / "docs" / "development" / "ci-followups.md").read_text(encoding="utf-8")
    workflow_source = (ROOT / ".github" / "workflows" / "publish-pypi.yml").read_text(encoding="utf-8")
    publish = _workflow("publish-pypi.yml")
    normalized_agents = " ".join(agents.split())

    assert all(name in normalized_agents for name in REQUIRED_CHECK_NAMES)
    for literal in (
        "vntyper_release",
        "client_payload.tag",
        "repository_dispatch",
        "default branch",
        "sha-<7>",
        "org.opencontainers.image.revision",
        "org.opencontainers.image.version",
        "digest",
        "main",
        "latest",
        "vX.Y.Z",
        "X.Y.Z",
        "fixed `vntyper-ghcr-release-promotion`",
        "existing tag",
        "dry run",
        "`pypi` environment",
        "`id-token: write`",
        "`PYPI_API_TOKEN` has been deleted after the green OIDC releases",
        "historical tagged commits cannot retrieve the obsolete credential",
        "proven short-prefix collision continues from the evidence digest",
        "ambiguous short-tag drift fails closed",
    ):
        assert literal in normalized_agents
    assert "mypy vntyper/ docker/app/ scripts/" in agents
    assert "`scripts/` is linted and formatted but is not type-checked" not in agents
    assert re.search(
        r"```text\nphase \| job \| permissions \| retry/recovery\n"
        r"validation \| validate-release \| actions: read, contents: read \|.*\n"
        r"gates \| wait-for-release-gates \| actions: read, checks: read, contents: read, packages: read \|.*\n"
        r"build \| build-package \| contents: read \|.*\n"
        r"promotion \| promote-ghcr \| contents: read, packages: write \|.*\n"
        r"publish \| publish-pypi \| id-token: write \|.*\n"
        r"summary \| release-summary \| none \|.*\n```",
        agents,
    )

    b4 = followups.split("### B4.", maxsplit=1)[1].split("### B5.", maxsplit=1)[0]
    assert "RESOLVED" in b4
    assert "cleanup complete" in b4
    assert "31465885545" in b4
    assert "31464328451" in b4
    assert "PyPI 2.0.12" in b4
    assert "PYPI_API_TOKEN has been deleted" in b4
    assert "remains intentionally pending" not in b4

    assert "client_payload.tag" in workflow_source
    assert "default branch" in workflow_source
    assert "push:" not in workflow_source.split("permissions:", maxsplit=1)[0]
    production_guard = "${{ github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release' }}"
    assert publish["jobs"]["promote-ghcr"]["if"] == production_guard
    assert publish["jobs"]["publish-pypi"]["if"] == production_guard


def test_release_summary_has_complete_always_run_safe_dataflow() -> None:
    """The diagnostic job must retain every upstream result without obtaining release authority."""
    publish = _workflow("publish-pypi.yml")
    job = publish["jobs"]["release-summary"]
    step = _step_with_id(publish, "release-summary", "summarize")

    assert job["needs"] == [
        "validate-release",
        "wait-for-release-gates",
        "build-package",
        "promote-ghcr",
        "publish-pypi",
    ]
    assert job["if"] == "${{ always() }}"
    assert job["permissions"] == {}
    assert "environment" not in job
    assert step["env"] == {
        "VALIDATE_RESULT": "${{ needs.validate-release.result }}",
        "VALIDATE_JSON": "${{ needs.validate-release.outputs.summary_json }}",
        "GATES_RESULT": "${{ needs.wait-for-release-gates.result }}",
        "PREFLIGHT_JSON": "${{ needs.wait-for-release-gates.outputs.preflight_summary_json }}",
        "CHECKS_JSON": "${{ needs.wait-for-release-gates.outputs.check_summary_json }}",
        "EVIDENCE_JSON": "${{ needs.wait-for-release-gates.outputs.evidence_summary_json }}",
        "BUILD_RESULT": "${{ needs.build-package.result }}",
        "PACKAGE_JSON": "${{ needs.build-package.outputs.package_summary_json }}",
        "PROMOTE_RESULT": "${{ needs.promote-ghcr.result }}",
        "ALIASES_JSON": "${{ needs.promote-ghcr.outputs.alias_summary_json }}",
        "DRY_ALIASES_JSON": "${{ needs.wait-for-release-gates.outputs.dry_run_alias_summary_json }}",
        "PUBLISH_RESULT": "${{ needs.publish-pypi.result }}",
        "PUBLISH_JSON": "${{ needs.publish-pypi.outputs.publish_summary_json }}",
        "EVENT_NAME": "${{ github.event_name }}",
    }
    assert "Dry run performed no production writes." in step["run"]
    assert "raw value intentionally not emitted" in step["run"]
    assert "json.loads(raw)" in step["run"]
    assert not {"GH_TOKEN", "GITHUB_TOKEN", "PYPI_API_TOKEN", "TWINE_PASSWORD"} & set(step["env"])


def test_release_summary_renders_all_provenance_and_handles_partial_invalid_state(tmp_path: Path) -> None:
    """Summary output must expose structured provenance while suppressing malformed raw values."""
    workflow = _workflow("publish-pypi.yml")
    summary_path = tmp_path / "summary.md"
    validation = {
        "mode": "dry-run",
        "tag": "v2.0.10",
        "sha": "a" * 40,
        "main_ancestor": True,
        "observed_versions": {"package": "2.0.10", "citation": "2.0.10", "changelog": "2.0.10"},
        "matches": {"package": True, "citation": True, "changelog": True},
    }
    preflight = {
        "state": "eligible",
        "reason": "checks exist",
        "eligible_runs": [{"id": 41, "url": "https://github.test/runs/41"}],
        "selected_run": {"id": 41, "url": "https://github.test/runs/41"},
    }
    checks = {
        "attempt": 2,
        "elapsed_seconds": 60,
        "verdicts": [
            {"name": name, "verdict": "success", "url": "https://github.test/check/1"} for name in REQUIRED_CHECK_NAMES
        ],
    }
    evidence = {
        "run": {"id": 41, "attempt": 3},
        "artifact": {"id": 501, "name": "docker-release-evidence"},
        "evidence_contract_version": 1,
        "source_ref": "ghcr.io/hassansaei/vntyper@sha256:" + "b" * 64,
        "digest": "sha256:" + "b" * 64,
        "revision": "a" * 40,
        "version": "2.0.10",
    }
    aliases = {
        "plan": [{"decision": decision} for decision in ("create", "advance", "no-op", "skip-newer", "fail-conflict")],
        "alias_progress": [{"attempted": True, "write_succeeded": True, "verified": False}],
    }
    package = {"artifact_name": "python-dist-2.0.10-1-1", "files": ["vntyper.whl", "vntyper.tar.gz"]}
    env = {
        "GITHUB_STEP_SUMMARY": str(summary_path),
        "EVENT_NAME": "workflow_dispatch",
        "VALIDATE_RESULT": "success",
        "VALIDATE_JSON": json.dumps(validation),
        "GATES_RESULT": "failure",
        "PREFLIGHT_JSON": json.dumps(preflight),
        "CHECKS_JSON": json.dumps(checks),
        "EVIDENCE_JSON": json.dumps(evidence),
        "BUILD_RESULT": "success",
        "PACKAGE_JSON": json.dumps(package),
        "PROMOTE_RESULT": "skipped",
        "ALIASES_JSON": json.dumps(aliases),
        "DRY_ALIASES_JSON": "",
        "PUBLISH_RESULT": "skipped",
        "PUBLISH_JSON": "not-json SUPER_SECRET_MUST_NOT_LEAK",
    }

    completed = _run_step(tmp_path, workflow, "release-summary", "summarize", env)

    assert completed.returncode == 0, completed.stderr
    rendered = summary_path.read_text(encoding="utf-8")
    assert "Dry run performed no production writes." in rendered
    for title in (
        "Validation",
        "Docker evidence preflight",
        "Exact-SHA checks",
        "Docker evidence",
        "Package",
        "GHCR aliases",
        "Dry-run aliases",
        "PyPI",
    ):
        assert f"## {title}" in rendered
    for literal in (
        "main_ancestor",
        "observed_versions",
        "elapsed_seconds",
        "evidence_contract_version",
        "source_ref",
        "alias_progress",
        "write_succeeded",
        "fail-conflict",
        "vntyper.whl",
        "Structured output: unavailable",
        "Structured output: invalid JSON",
    ):
        assert literal in rendered
    assert "SUPER_SECRET_MUST_NOT_LEAK" not in rendered
