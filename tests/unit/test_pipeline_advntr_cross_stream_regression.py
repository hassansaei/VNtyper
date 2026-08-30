"""End-to-end regression coverage for mixed-stream adVNTR version answers (#289)."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from shlex import quote

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness

pytestmark = pytest.mark.unit


def _write_mixed_stream_advntr(path: Path, invocation_log: Path) -> None:
    """Create an adVNTR probe that separates mamba diagnostics from its answer."""
    path.write_text(
        "#!/bin/sh\n"
        f"printf '%s\\n' \"$*\" >> {quote(str(invocation_log))}\n"
        "printf '%s\\n' \"warning libmamba Cannot lock '/tmp/mamba/proc'\"\n"
        "printf '%s\\n' \"Waiting for other mamba process to finish\"\n"
        "printf '%s\\n' '2.0.4' >&2\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def test_a_stderr_advntr_answer_survives_stdout_process_lock_diagnostics(tmp_path: Path) -> None:
    """Restoring stdout priority would abort a valid cross-stream adVNTR pipeline run."""
    invocation_log = tmp_path / "advntr-invocations.log"
    fake_advntr = tmp_path / "advntr-mixed-stream"
    _write_mixed_stream_advntr(fake_advntr, invocation_log)
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = str(fake_advntr)

    harness = run_pipeline_under_harness(tmp_path / "out", config=config, extra_modules=["advntr"])

    assert harness.error is None
    harness.stages["run_kestrel"].assert_called_once()
    assert invocation_log.read_text(encoding="utf-8") == "--version\n"
    assert harness.kwargs("get_tool_versions")["version_overrides"] == {"advntr": "2.0.4"}
