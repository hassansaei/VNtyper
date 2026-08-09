"""Every `conda run` that launches a long-running child must stream its output.

`conda run` buffers **all** child stdout and stderr until the child exits unless
``--no-capture-output`` is passed. A pipeline running under the buffering form therefore
produces no output at all in ``docker logs`` no matter how much progress it is making,
which is why #178 was reported as "hangs indefinitely, no output, no error" rather than
as a stage log ending at a specific step (#213).

There are **two** such invocations, and the issue only names one:

* ``docker/entrypoint.sh`` -- the container's own entry point, which is how a CLI user
  runs the image.
* ``docker/app/tasks.py`` -- **two** invocations: ``build_vntyper_command``, which is how
  every *web* job runs the pipeline, and the ``vntyper cohort`` command built inline in
  ``run_cohort_analysis_job``. Fixing only the entrypoint makes the Celery worker's own
  output stream while the pipeline it launches stays buffered, so a web job would still
  run for an hour and log nothing.

An earlier version of this module asserted only on ``build_vntyper_command`` while
claiming to cover "every long-running child", and the cohort launcher was buffered behind
that claim. The source-level scan below exists so the claim is actually true: it finds
every ``conda run`` in the module, however it is built.

This test reads both as **source text**, so it fails on a new invocation added later. It
lives in the unit tier rather than the docker tier deliberately: the docker tier needs a
daemon and does not run on every PR, and the property being checked is a property of the
source, not of a built image. It sits under ``tests/unit/web/`` because that package's
``conftest.py`` is what puts ``docker`` on ``sys.path``, which is how ``app.tasks``
becomes importable.

The build-time invocations in ``docker/Dockerfile`` and ``docker/Dockerfile.base`` are
excluded on purpose. They run during an image build rather than during a job, and both
Dockerfiles are inputs to the content hash that binds the two halves of the image
(AGENTS.md trap 10), so editing them forces a base rebuild for no runtime benefit.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[3]

#: `conda run` invocations in the entrypoint that hand control to a long-running child.
#: The `--version` probe at :108 is excluded: it is a liveness check whose output is sent
#: to /dev/null, it exits immediately, and buffering it changes nothing.
ENTRYPOINT_EXEC_SITES = 4


def _entrypoint_exec_lines() -> list[str]:
    """Return the `exec conda run ...` lines from the container entrypoint.

    Returns:
        list[str]: One entry per invocation, whitespace-normalised. Continuation lines
        are joined so a flag on the next line still counts as present.
    """
    text = (REPO_ROOT / "docker" / "entrypoint.sh").read_text(encoding="utf-8")
    joined = re.sub(r"\\\s*\n\s*", " ", text)
    return [" ".join(line.split()) for line in joined.splitlines() if "exec conda run" in line]


def test_the_entrypoint_has_the_expected_number_of_exec_sites():
    """A new invocation must be classified deliberately, not slip past the check."""
    lines = _entrypoint_exec_lines()
    assert len(lines) == ENTRYPOINT_EXEC_SITES, (
        f"expected {ENTRYPOINT_EXEC_SITES} `exec conda run` sites in docker/entrypoint.sh, "
        f"found {len(lines)}. If you added one, decide whether it streams and update "
        f"ENTRYPOINT_EXEC_SITES. Lines: {lines}"
    )


def test_every_entrypoint_exec_streams_its_child_output():
    """#213: all four buffered, so a running pipeline logged nothing."""
    buffering = [line for line in _entrypoint_exec_lines() if "--no-capture-output" not in line]
    assert not buffering, (
        "these `conda run` invocations buffer their child's output until it exits, so a "
        f"running pipeline emits nothing in `docker logs`: {buffering}"
    )


def test_the_web_worker_launches_the_pipeline_with_streaming_enabled():
    """The path every web job takes, which #213 does not mention.

    Asserted against the built argument vector rather than the source text, because this
    one is a list of tokens and a substring check would pass on a comment.
    """
    from app.tasks import build_vntyper_command

    command = build_vntyper_command(
        alignment_path="/data/sample.cram",
        output_dir="/out",
        thread=4,
        reference_assembly="hg38",
    )
    assert command[:2] == ["conda", "run"], f"unexpected launcher: {command[:2]}"
    assert "--no-capture-output" in command, (
        "the Celery worker launches the pipeline through its own `conda run`; without "
        f"--no-capture-output every web job's output is buffered until it exits: {command}"
    )


def test_every_conda_run_in_the_worker_streams_its_child_output():
    """Catches an invocation built inline rather than through a helper.

    ``run_cohort_analysis_job`` assembles its own argument list, so a test that only
    calls ``build_vntyper_command`` cannot see it. This reads the module as source text
    and checks each ``conda`` / ``run`` pair, so a third launcher added later fails here
    rather than shipping buffered.
    """
    source = (REPO_ROOT / "docker" / "app" / "tasks.py").read_text(encoding="utf-8")
    launchers = re.findall(r'"conda",\s*\n\s*"run",\s*\n(.*?)"-n",', source, re.DOTALL)
    assert launchers, "no `conda run` argument vectors found in docker/app/tasks.py"
    buffering = [block for block in launchers if "--no-capture-output" not in block]
    assert not buffering, (
        f"{len(buffering)} of {len(launchers)} `conda run` launchers in docker/app/tasks.py "
        "buffer their child's output until it exits"
    )


def test_the_streaming_flag_precedes_the_environment_selection():
    """`conda run` parses its own options before `-n`; order is not cosmetic."""
    from app.tasks import build_vntyper_command

    command = build_vntyper_command(
        alignment_path="/data/sample.bam",
        output_dir="/out",
        thread=1,
        reference_assembly="hg19",
    )
    assert command.index("--no-capture-output") < command.index("-n")
