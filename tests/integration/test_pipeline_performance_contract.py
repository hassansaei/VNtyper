"""Real command-count acceptance checks for milestone-four pipeline work."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration


def test_default_bam_run_indexes_the_final_slice_exactly_once(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """A-PERF-2 counts the real index command from the pipeline log."""
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
    assert completed.returncode == 1
    assert "FASTQ layout 'mixed' cannot be consumed without dropping reads" in completed.stderr
    assert len(commands) == 1
