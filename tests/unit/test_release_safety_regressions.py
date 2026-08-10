"""Hostile regressions for release-authoritative Docker mutations."""

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest
import yaml

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]
IMAGE = "ghcr.io/hassansaei/vntyper"
SOURCE_DIGEST = "sha256:" + "b" * 64
OTHER_DIGEST = "sha256:" + "c" * 64
MAIN_PUSH_GUARD = "${{ github.event_name == 'push' && github.ref == 'refs/heads/main' }}"


def _workflow(name: str) -> dict[str, Any]:
    """Load one workflow without YAML 1.1 coercing the ``on`` key."""
    return yaml.load((ROOT / ".github" / "workflows" / name).read_text(), Loader=yaml.BaseLoader)


def _named_step(workflow: dict[str, Any], job: str, name: str) -> dict[str, Any]:
    """Return one workflow step selected by its display name."""
    return next(step for step in workflow["jobs"][job]["steps"] if step.get("name") == name)


def _id_step(workflow: dict[str, Any], job: str, step_id: str) -> dict[str, Any]:
    """Return one workflow step selected by its identifier."""
    return next(step for step in workflow["jobs"][job]["steps"] if step.get("id") == step_id)


def _write_executable(path: Path, source: str) -> None:
    """Write one executable used only inside an isolated test directory."""
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)


def test_only_exact_main_push_may_publish_release_authoritative_application_tags() -> None:
    """Nightly/manual validation must not overwrite the push-bound short-SHA image."""
    workflow = _workflow("docker-build.yml")
    job = workflow["jobs"]["build-and-test"]
    push = _named_step(workflow, "build-and-test", "Push image")
    registry = _id_step(workflow, "build-and-test", "registry-digest")
    evidence = _id_step(workflow, "build-and-test", "release-evidence")
    upload = next(step for step in job["steps"] if step.get("uses") == "actions/upload-artifact@v5")
    docker_tests = _named_step(workflow, "build-and-test", "Run Docker tests")["run"]

    assert "schedule" in workflow["on"]
    assert "workflow_dispatch" in workflow["on"]
    assert 'github.event_name }}" = "schedule"' in docker_tests
    assert 'github.event_name }}" = "workflow_dispatch"' in docker_tests
    assert "inputs.full" not in docker_tests
    assert {push["if"], registry["if"], evidence["if"], upload["if"]} == {MAIN_PUSH_GUARD}
    assert "docker push" in push["run"]
    assert (
        "type=sha"
        in next(step for step in job["steps"] if step.get("uses") == "docker/metadata-action@v6")["with"]["tags"]
    )


def test_release_contract_documents_nightly_immutability_and_legacy_alias_migration() -> None:
    """The reviewed safety correction must remain in every governing contract."""
    paths = (
        ROOT / "AGENTS.md",
        ROOT / "docs" / "superpowers" / "specs" / "2026-08-10-milestone-6-release-and-naming-design.md",
        ROOT / "docs" / "superpowers" / "plans" / "2026-08-10-milestone-6-release-and-naming-plan.md",
    )
    for path in paths:
        text = " ".join(path.read_text(encoding="utf-8").split())
        assert "scheduled and manual Docker validation never publish application tags" in text
        assert "legacy rolling `main`" in text
        assert "missing or unrecognized version label" in text
        assert "fails closed" in text
        assert "skip-unorderable" not in text


def test_main_push_publication_shell_still_pushes_main_and_short_sha_tags(tmp_path: Path) -> None:
    """Narrowing the event guard must retain the checked-in main publication loop."""
    workflow = _workflow("docker-build.yml")
    step = _named_step(workflow, "build-and-test", "Push image")
    tags = f"{IMAGE}:main\n{IMAGE}:sha-aaaaaaa"
    run = step["run"].replace("${{ steps.meta.outputs.tags }}", tags)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    command_log = tmp_path / "docker-commands.jsonl"
    _write_executable(
        fake_bin / "docker",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

with pathlib.Path(os.environ["DOCKER_COMMAND_LOG"]).open("a", encoding="utf-8") as handle:
    handle.write(json.dumps(sys.argv[1:]) + "\\n")
""",
    )

    completed = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", run],
        cwd=tmp_path,
        env={
            **os.environ,
            "DOCKER_COMMAND_LOG": str(command_log),
            "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
        },
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert [json.loads(line) for line in command_log.read_text(encoding="utf-8").splitlines()] == [
        ["tag", "vntyper:test", f"{IMAGE}:main"],
        ["push", f"{IMAGE}:main"],
        ["tag", "vntyper:test", f"{IMAGE}:sha-aaaaaaa"],
        ["push", f"{IMAGE}:sha-aaaaaaa"],
    ]


def _promotion_environment(tmp_path: Path, floating_version: str | None) -> dict[str, str]:
    """Install deterministic registry observations for the real promotion shell."""
    controller = tmp_path / "controller" / "scripts"
    controller.mkdir(parents=True)
    for script in ("release_policy.py", "release_registry.py", "release_manifest.py"):
        shutil.copy2(ROOT / "scripts" / script, controller / script)
    observations = {
        f"{IMAGE}:v2.0.10": {"digest": SOURCE_DIGEST, "version": "2.0.10"},
        f"{IMAGE}:2.0.10": {"digest": SOURCE_DIGEST, "version": "2.0.10"},
        f"{IMAGE}:2.0": {"digest": OTHER_DIGEST, "version": floating_version},
        f"{IMAGE}:2": {"digest": SOURCE_DIGEST, "version": "2.0.10"},
        f"{IMAGE}:latest": {"digest": SOURCE_DIGEST, "version": "2.0.10"},
    }
    observations_path = tmp_path / "observations.json"
    observations_path.write_text(json.dumps(observations), encoding="utf-8")
    mutation_log = tmp_path / "mutation.log"
    mutation_log.write_text("", encoding="utf-8")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _write_executable(
        fake_bin / "docker",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

args = sys.argv[1:]
if args[:3] == ["buildx", "imagetools", "create"]:
    with pathlib.Path(os.environ["MUTATION_LOG"]).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(args) + "\\n")
    raise SystemExit(0)
if args[:3] != ["buildx", "imagetools", "inspect"]:
    raise SystemExit(97)
observations = json.loads(pathlib.Path(os.environ["OBSERVATIONS"]).read_text(encoding="utf-8"))
record = observations[args[3]]
format_value = args[args.index("--format") + 1]
if format_value == "{{json .Manifest}}":
    print(json.dumps({"mediaType": "application/vnd.oci.image.manifest.v1+json",
                      "digest": record["digest"], "size": 3464}))
elif format_value == "{{json .Image}}":
    print(json.dumps({"config": {"Labels": {"org.opencontainers.image.version": record["version"]}}}))
else:
    raise SystemExit(96)
""",
    )
    _write_executable(fake_bin / "sleep", "#!/usr/bin/env bash\nexit 0\n")
    return {
        **os.environ,
        "IMAGE": IMAGE,
        "MUTATION_LOG": str(mutation_log),
        "OBSERVATIONS": str(observations_path),
        "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}",
        "RELEASE_ALIAS_INSPECT_ATTEMPTS": "3",
        "RELEASE_ALIAS_INSPECT_RETRY_SECONDS": "0",
        "SOURCE_DIGEST": SOURCE_DIGEST,
        "VERSION": "2.0.10",
    }


@pytest.mark.parametrize("floating_version", (None, "nightly", "2.0"))
def test_unrecognized_floating_label_aborts_real_promotion_before_any_write(
    tmp_path: Path,
    floating_version: str | None,
) -> None:
    """An unorderable registry label must not leave a stale alias behind with green status."""
    workflow = _workflow("publish-pypi.yml")
    promote = _id_step(workflow, "promote-ghcr", "promote")
    env = _promotion_environment(tmp_path, floating_version)

    completed = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", promote["run"]],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode != 0
    assert "release alias conflict" in completed.stderr
    assert "cannot be ordered safely" in completed.stderr
    assert (tmp_path / "mutation.log").read_text(encoding="utf-8") == ""
