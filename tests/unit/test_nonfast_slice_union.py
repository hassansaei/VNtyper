"""Unit contract for disjoint non-fast target slicing and complete recovery."""

import pytest

from vntyper.scripts.command_builders import build_samtools_slice_command

pytestmark = pytest.mark.unit


def test_the_nonfast_slice_excludes_reads_recovered_by_the_complete_flag_four_scan() -> None:
    """The target half of the merge must be disjoint from whole-stream recovery."""
    command = build_samtools_slice_command(
        samtools_path="samtools",
        in_bam="/data/sample.bam",
        output_bam="/out/output_sliced.bam.partial",
        region="chr1:155158000-155163000",
        exclude_unmapped=True,
        index_output=False,
    )

    assert command == (
        "samtools view -P -b -F 4 /data/sample.bam chr1:155158000-155163000 -o /out/output_sliced.bam.partial"
    )
