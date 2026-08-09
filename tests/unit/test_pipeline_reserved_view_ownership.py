"""Full-pipeline ownership guard contract for safely replaced view symlinks."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("link_state", ["stale", "dangling"])
@pytest.mark.parametrize("input_mode", ["bam", "cram", "fastq"])
def test_reserved_alignment_view_symlink_is_admitted_for_safe_builder_replacement(
    tmp_path: Path,
    input_mode: str,
    link_state: str,
) -> None:
    """Only the pre-existing snapshot is checked; the builder rechecks before replacement."""
    output = tmp_path / "output"
    stage = output / "fastq_bam_processing"
    stage.mkdir(parents=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    if input_mode == "fastq":
        selected_input = inputs / "reads.fastq.gz"
        reference = inputs / "reference.fa"
        selected_input.write_bytes(b"reads")
        reference.write_bytes(b"reference")
        reserved_view = stage / "post_alignment.bam"
        input_kwargs: dict[str, Any] = {"fastq1": str(selected_input), "bwa_reference": str(reference)}
    else:
        selected_input = inputs / f"selected.{input_mode}"
        selected_input.write_bytes(b"selected")
        reserved_view = stage / f"input.{input_mode}"
        input_kwargs = {input_mode: str(selected_input)}

    stale_target = inputs / f"different-{input_mode}"
    if link_state == "stale":
        stale_target.write_bytes(b"different patient")
    reserved_view.symlink_to(stale_target)

    harness = run_pipeline_under_harness(output, **input_kwargs)

    assert harness.error is None
    if input_mode == "fastq":
        harness.call("process_fastq")
        harness.call("align_and_sort_fastq")
    else:
        harness.call("run_preflight")
        harness.call("process_bam_to_fastq")


@pytest.mark.parametrize("input_mode", ["bam", "fastq"])
def test_nonreserved_symlink_remains_rejected_before_pipeline_work(tmp_path: Path, input_mode: str) -> None:
    """The safe-builder exception does not admit any other pre-existing symlink."""
    output = tmp_path / "output"
    output.mkdir()
    external = tmp_path / "external"
    external.write_bytes(b"external bytes")
    malicious = output / "pipeline_summary.json"
    malicious.symlink_to(external)
    original_inode = external.stat().st_ino
    if input_mode == "bam":
        alignment = tmp_path / "patient" / "patient.bam"
        alignment.parent.mkdir()
        alignment.write_bytes(b"patient")
        input_kwargs: dict[str, Any] = {"bam": str(alignment)}
    else:
        fastq = tmp_path / "reads.fastq.gz"
        reference = tmp_path / "reference.fa"
        fastq.write_bytes(b"reads")
        reference.write_bytes(b"reference")
        input_kwargs = {"fastq1": str(fastq), "bwa_reference": str(reference)}

    harness = run_pipeline_under_harness(output, expect_failure=True, **input_kwargs)

    assert isinstance(harness.error, SystemExit)
    assert external.stat().st_ino == original_inode
    assert external.read_bytes() == b"external bytes"
    for stage in harness.stages.values():
        stage.assert_not_called()
