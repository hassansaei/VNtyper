"""Focused tests for threading local CRAM header references into preflight."""

from pathlib import Path

import pytest

from vntyper.scripts.pipeline_alignment import build_alignment_preflight_kwargs

pytestmark = pytest.mark.unit


def test_preflight_kwargs_resolve_local_header_references_against_the_original_cram(tmp_path: Path) -> None:
    """The owned input boundary, not the run-local view directory, anchors relative UR values."""
    input_cram = tmp_path / "patient" / "input.cram"
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t10\t20\n", encoding="utf-8")

    result = build_alignment_preflight_kwargs(
        in_path=input_cram,
        output_dir=tmp_path / "run" / "fastq_bam_processing",
        output_name="input",
        file_format="cram",
        config={},
        threads=2,
        bed_file=bed,
        reference_assembly="hg19",
        fast_mode=True,
        alignment_header="@SQ\tSN:chr1\tLN:100\tUR:reference.fa\n",
    )

    assert result["header_reference_paths"] == (str(input_cram.parent / "reference.fa"),)
