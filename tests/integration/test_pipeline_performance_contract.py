"""Real command-count acceptance checks for milestone-four pipeline work."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration


def test_a_kestrel_only_bam_run_does_not_index_the_final_slice(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """A-PERF-2 counts the real index command from the pipeline log.

    A-PERF-2 was about *double* indexing: the pipeline had to index the final slice
    once rather than twice, and this test counted the real command to prove it. Since
    #262 the count in a Kestrel-only run is **zero**, which is a deliberate revision of
    that contract rather than a relaxation of this assertion.

    The reason is a complete consumer trace of ``<name>_sliced.bam.bai``: its only
    readers are ``run_advntr`` and ``downsample_bam_if_needed``, both reachable only
    from ``pipeline.py``'s ``if "advntr" in extra_modules`` block. Coverage reads
    ``plan.view_path`` with ``plan.stable_index_path``, not this file, and
    ``artifact_names.py``'s ``sliced_bam`` entry has no consumer at all. The index is
    not a required artifact of any golden-cohort case or of any row in
    ``tests/compatibility/real_success_baseline.json``.

    The invocation is unchanged, ``--keep-intermediates`` included, so this still
    measures the configuration A-PERF-2 measured. The adVNTR arm -- where the index is
    still emitted exactly once -- is covered by the golden-cohort gate's
    ``executed_commands`` comparison on its adVNTR cases and by
    ``tests/unit/test_fastq_bam_command_wiring.py``, rather than by a second real run
    here: adVNTR costs 15-25 minutes and is deliberately off the merge path.
    """
    repository = Path(__file__).parents[2]
    source = repository / "tests/data/remapped/bwa/hg19/example_b178_hg19_bwa.bam"
    output = tmp_path / "output"
    log_file = tmp_path / "pipeline.log"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "vntyper.cli",
            "--config-path",
            str(repository / "vntyper/config.json"),
            "--log-level",
            "DEBUG",
            "--log-file",
            str(log_file),
            "pipeline",
            "--bam",
            str(source),
            "--reference-assembly",
            "hg19",
            "--threads",
            "2",
            "--keep-intermediates",
            "--output-dir",
            str(output),
        ],
        cwd=repository,
        capture_output=True,
        text=True,
        check=False,
    )

    commands = [
        line
        for line in log_file.read_text(encoding="utf-8").splitlines()
        if "Re-indexing BAM file with command: samtools index" in line and "output_sliced.bam" in line
    ]
    assert completed.returncode == 0, completed.stderr
    assert "Pipeline finished successfully." in completed.stderr
    assert commands == []
