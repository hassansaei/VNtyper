"""Resource-lifetime tests for the extracted alignment coverage boundary."""

from pathlib import Path

import pytest

from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.pipeline_coverage import calculate_alignment_coverage

pytestmark = pytest.mark.unit


def test_coverage_failure_still_releases_the_alignment_binding(tmp_path: Path) -> None:
    """A failed final consumer must not leak its parent-held descriptor."""
    alignment = tmp_path / "patient.bam"
    alignment.write_bytes(b"patient alignment")
    binding = AlignmentBinding(str(alignment))
    plan = AlignmentPlan(
        input_path=str(alignment),
        view_path="/run/input.bam",
        file_format="bam",
        index_path="/run/input.bam.bai",
        reference_path=None,
        reference_source="not-required",
        uncovered_contigs=(),
        unmapped_scan="indexed",
        binding=binding,
    )

    def fail_coverage(**kwargs: object) -> None:
        del kwargs
        raise RuntimeError("coverage failed")

    with pytest.raises(RuntimeError, match="coverage failed"):
        calculate_alignment_coverage(
            plan=plan,
            region="chr1:10-20",
            reference_assembly="hg19",
            threads=2,
            config={},
            output_dir=str(tmp_path / "coverage"),
            coverage_calculator=fail_coverage,
            region_resolver=lambda **kwargs: "unused",
        )

    assert binding.is_open is False
