#!/usr/bin/env python3
"""
tests/docker/test_image_structure.py

Fast structural smoke tests for the VNtyper application image.

Complements tests/docker/test_docker_pipeline.py rather than replacing it. This tier
needs no Zenodo test data and asserts only on the image's structure, so it runs on every
PR in a couple of seconds; the pipeline tier still owns output correctness.

Two fixtures, two costs:
  - `image_metadata` reads image config from the daemon and starts no container at all.
  - `probe` starts exactly ONE container, pipes image_probe.py to the image's own
    interpreter, and returns the parsed JSON. Every runtime assertion below reads that
    one blob, so twenty assertions cost one container start rather than twenty.

Marked `smoke` only, NOT `docker`: `make test-docker` runs `-m docker`, so adding both
markers would make the slow tier re-run all of these.

These tests never build an image - see the note in conftest.py about a 1042s fixture.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import pytest

IMAGE = os.environ.get("VNTYPER_TEST_IMAGE", "vntyper:local")
PROBE = Path(__file__).parent / "image_probe.py"

# Generous on purpose. The image is ~4.8 GiB; this catches a regression of the class
# documented in docker/Dockerfile (a 680 MB duplicated reference layer, a 1.66 GB
# recursive chown), not ordinary drift. A tight budget is a flaky budget.
MAX_IMAGE_BYTES = 6 * 1024**3

# org.opencontainers.image.version and .ref.name are inherited from the Ubuntu base
# ("20.04", "ubuntu") on local builds and only corrected by docker/metadata-action in
# CI, so they are deliberately not asserted.
REQUIRED_LABELS = (
    "org.opencontainers.image.title",
    "org.opencontainers.image.source",
    "org.opencontainers.image.licenses",
)
REQUIRED_BINARIES = ("bwa", "samtools", "fastp", "bcftools", "java")

# Floors, never equality: BWA index sizes shift with the indexer version. These only
# need to catch a truncated download or a half-written layer.
MIN_REFERENCE_BYTES = {
    "bwa_reference_hg19": 200 * 1024**2,
    "bwa_reference_hg38": 200 * 1024**2,
}


def _docker_available() -> bool:
    """Report whether a usable Docker daemon is present.

    This runs at collection time via the module-level skipif below, so it must never
    block: a daemon that is present but wedged (starting up, or stuck on its socket)
    would otherwise stall the whole test session. Treat non-responsive as unavailable.

    Returns:
        bool: True if the `docker` CLI exists and the daemon answers within 10s.
    """
    if shutil.which("docker") is None:
        return False
    try:
        completed = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            check=False,
            timeout=10,
        )
    except (subprocess.TimeoutExpired, OSError):
        return False
    return completed.returncode == 0


pytestmark = [
    pytest.mark.smoke,
    pytest.mark.skipif(not _docker_available(), reason="Docker daemon not available"),
]


@pytest.fixture(scope="session")
def image_metadata() -> dict[str, Any]:
    """Return image configuration straight from the daemon; starts no container.

    Returns:
        dict: Parsed `docker image inspect` output.
    """
    result = subprocess.run(
        ["docker", "image", "inspect", IMAGE],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.skip(f"image {IMAGE} not present locally; set VNTYPER_TEST_IMAGE or build it first")
    return json.loads(result.stdout)[0]


@pytest.fixture(scope="session")
def probe(image_metadata: dict[str, Any]) -> dict[str, Any]:
    """Run image_probe.py inside the image once and return its JSON report.

    `--network none` guarantees no assertion here can silently depend on Zenodo, UCSC
    or conda-forge being reachable.

    Args:
        image_metadata: Ensures the image exists before a container is started.

    Returns:
        dict: The probe's report.
    """
    with PROBE.open("rb") as handle:
        result = subprocess.run(
            [
                "docker",
                "run",
                "--rm",
                "-i",
                "--network",
                "none",
                "--entrypoint",
                "/opt/conda/envs/vntyper/bin/python",
                IMAGE,
                "-",
            ],
            stdin=handle,
            capture_output=True,
            text=True,
            check=False,
        )
    assert result.returncode == 0, f"probe failed in {IMAGE}:\n{result.stderr}"
    return json.loads(result.stdout)


# --- image metadata: no container started ------------------------------------


def test_runs_as_non_root(image_metadata: dict[str, Any]) -> None:
    """The image must declare the unprivileged user."""
    assert image_metadata["Config"]["User"] == "appuser"


def test_healthcheck_declared(image_metadata: dict[str, Any]) -> None:
    """A HEALTHCHECK must survive Dockerfile refactors."""
    assert image_metadata["Config"]["Healthcheck"]["Test"][:1] == ["CMD"]


def test_entrypoint_is_exec_form(image_metadata: dict[str, Any]) -> None:
    """Exec form keeps signal handling intact under docker stop / k8s."""
    assert image_metadata["Config"]["Entrypoint"] == ["/usr/local/bin/entrypoint.sh"]


def test_production_entrypoint_streams_child_output_before_the_child_exits(
    image_metadata: dict[str, Any], tmp_path: Path
) -> None:
    """A live child sentinel must reach Docker logs through the real entrypoint.

    The bind mount replaces only the command exercised by the entrypoint; the
    image, entrypoint, conda environment and ``conda run`` invocation remain the
    production ones. The child stays alive for 30 seconds so observing the
    sentinel and then observing ``State.Running=true`` proves output arrived
    before process exit rather than being flushed during teardown.

    Args:
        image_metadata: Ensures the configured image exists before launch.
        tmp_path: Scratch location for the controlled child executable.
    """
    assert image_metadata["Config"]["Entrypoint"] == ["/usr/local/bin/entrypoint.sh"]
    sentinel = "VNTYPER_ENTRYPOINT_STREAM_SENTINEL"
    child = tmp_path / "vntyper"
    child.write_text(f'#!/bin/bash\nprintf "%s\\n" "{sentinel}"\nsleep 30\n', encoding="utf-8")
    child.chmod(0o755)

    started = subprocess.run(
        [
            "docker",
            "run",
            "--detach",
            "--rm",
            "--network",
            "none",
            "--mount",
            f"type=bind,src={child},dst=/opt/conda/envs/vntyper/bin/vntyper,readonly",
            IMAGE,
            "vntyper",
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=15,
    )
    assert started.returncode == 0, started.stderr
    container_id = started.stdout.strip()
    try:
        deadline = time.monotonic() + 8
        logs = ""
        while sentinel not in logs and time.monotonic() < deadline:
            captured = subprocess.run(
                ["docker", "logs", container_id],
                capture_output=True,
                text=True,
                check=False,
                timeout=5,
            )
            logs = captured.stdout + captured.stderr
            if sentinel not in logs:
                time.sleep(0.05)

        assert sentinel in logs, f"child output was not streamed while container {container_id} was alive:\n{logs}"
        running = subprocess.run(
            ["docker", "inspect", "--format", "{{.State.Running}}", container_id],
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
        assert running.returncode == 0, running.stderr
        assert running.stdout.strip() == "true", "sentinel appeared only after the controlled child exited"
    finally:
        subprocess.run(
            ["docker", "rm", "--force", container_id],
            capture_output=True,
            check=False,
            timeout=15,
        )


def test_port_exposed(image_metadata: dict[str, Any]) -> None:
    """docker-compose and k8s wiring depend on the declared port."""
    assert "8000/tcp" in image_metadata["Config"]["ExposedPorts"]


def test_workdir_is_prefix(image_metadata: dict[str, Any]) -> None:
    """config.json's reference paths are relative and resolve against WORKDIR."""
    assert image_metadata["Config"]["WorkingDir"] == "/opt/vntyper"


@pytest.mark.parametrize("label", REQUIRED_LABELS)
def test_oci_label_present(image_metadata: dict[str, Any], label: str) -> None:
    """Provenance labels for a published scientific tool.

    Args:
        image_metadata: Image config.
        label: Label key that must be set.
    """
    assert image_metadata["Config"]["Labels"].get(label)


def test_image_size_budget(image_metadata: dict[str, Any]) -> None:
    """Fail loudly on a large size regression rather than discovering it in a registry."""
    size = image_metadata["Size"]
    assert size < MAX_IMAGE_BYTES, f"image is {size / 1024**3:.2f} GiB, budget is {MAX_IMAGE_BYTES / 1024**3:.2f} GiB"


# --- runtime probe: one container start, shared by every test below ----------


def test_uid_is_1001(probe: dict[str, Any]) -> None:
    """Volume-mount permissions depend on the effective UID, not just Config.User."""
    assert probe["uid"] == 1001


def test_cli_version_matches_package(probe: dict[str, Any]) -> None:
    """A bare --version proves the entry point runs; equality proves *what* ran.

    Catches a failed `pip install --no-deps`, a stale layer cache serving an old
    checkout, and a broken console-script entry point.
    """
    assert probe["cli_version"]["rc"] == 0, probe["cli_version"]["err"]
    assert probe["cli_version"]["out"].split()[-1] == probe["package_version"]


@pytest.mark.parametrize("env", ["vntyper", "envadvntr", "shark_env"])
def test_conda_env_exists(probe: dict[str, Any], env: str) -> None:
    """All three environments must survive the COPY --from into the runtime stage.

    Args:
        probe: Probe report.
        env: Conda environment name.
    """
    assert probe["conda_env_dirs"][env]


@pytest.mark.parametrize(("env", "expected"), [("vntyper", "3.12"), ("envadvntr", "2.7")])
def test_interpreter_version(probe: dict[str, Any], env: str, expected: str) -> None:
    """adVNTR needs Python 2.7; the pipeline needs 3.12. Neither may drift.

    Args:
        probe: Probe report.
        env: Conda environment name.
        expected: Expected major.minor version.
    """
    assert probe["interpreters"][env]["out"] == expected


def test_shark_binary_present(probe: dict[str, Any]) -> None:
    """shark_env ships a single C++ binary and no interpreter - assert what it is."""
    assert probe["shark_binary"]


def test_advntr_imports_under_py27(probe: dict[str, Any]) -> None:
    """The compiled adVNTR egg is the most fragile artifact in the image.

    It is built in a cloned tree that is deleted afterwards, so a broken .so would
    otherwise only surface under `--extra-modules advntr` in the slow tier.
    """
    assert probe["advntr_import"]["rc"] == 0, probe["advntr_import"]["err"]


@pytest.mark.parametrize("binary", REQUIRED_BINARIES)
def test_binary_resolves(probe: dict[str, Any], binary: str) -> None:
    """A conda solve that drops a tool is invisible until mid-pipeline.

    Args:
        probe: Probe report.
        binary: Executable expected on PATH.
    """
    assert probe["binaries"][binary], f"{binary} is not on PATH in the image"


def test_kestrel_jar_present(probe: dict[str, Any]) -> None:
    """The Kestrel JAR is package-data; a packaging change can drop it silently."""
    assert probe["kestrel_jar"]


def test_no_compiler_leaked(probe: dict[str, Any]) -> None:
    """Verify the multi-stage split held.

    `git` is excluded on purpose - condaforge/miniforge3 installs it in its own layer,
    so asserting on it would be permanently red. See image_probe.py.
    """
    assert probe["leaked_build_tools"] == {}


def test_every_declared_reference_exists(probe: dict[str, Any]) -> None:
    """The regression this whole tier exists for.

    config.json names reference files that a .dockerignore rule, a changed COPY, or a
    stale base image can silently omit - this actually happened: three paths
    (both adVNTR databases and SHARK's muc1_region_hg19.fa) went missing when
    .dockerignore excluded all of reference/. Nothing in the build failed, because
    Dockerfile.base's `test -f` guards only run during a base build.
    """
    missing = {key: entry["resolved"] for key, entry in probe["references"].items() if not entry["exists"]}
    assert not missing, f"config declares paths absent from the image: {missing}"


@pytest.mark.parametrize(("key", "minimum"), sorted(MIN_REFERENCE_BYTES.items()))
def test_reference_not_truncated(probe: dict[str, Any], key: str, minimum: int) -> None:
    """Catch truncated downloads and half-written layers with a floor, never equality.

    Args:
        probe: Probe report.
        key: config.json reference key.
        minimum: Minimum plausible size in bytes.
    """
    assert probe["references"][key]["size"] >= minimum
