"""Resource-lifetime tests for the extracted alignment coverage boundary."""

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.pipeline_coverage import calculate_alignment_coverage

pytestmark = pytest.mark.unit


def test_coverage_receives_the_descriptor_bound_index_path() -> None:
    """The final htslib consumer must not receive the replaceable public index path."""
    binding = MagicMock(spec=AlignmentBinding)
    binding.index_view_path = "/proc/123/fd/12"
    plan = AlignmentPlan(
        input_path="/input/patient.cram",
        view_path="/run/input.cram",
        file_format="cram",
        index_path="/run/input.cram.crai",
        reference_path="/reference.fa",
        reference_source="test",
        uncovered_contigs=(),
        unmapped_scan="indexed",
        binding=binding,
    )
    received: list[object] = []

    def record_coverage(**kwargs: object) -> None:
        received.append(kwargs["index_path"])

    calculate_alignment_coverage(
        plan=plan,
        region="chr1:10-20",
        reference_assembly="hg19",
        threads=2,
        config={},
        output_dir="/coverage",
        coverage_calculator=record_coverage,
        region_resolver=lambda **kwargs: "unused",
    )

    assert received == ["/proc/123/fd/12"]


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


def test_region_resolution_failure_still_releases_the_alignment_binding(tmp_path: Path) -> None:
    """The lifetime guard must include fallback region resolution through the plan view."""
    alignment = tmp_path / "patient.bam"
    alignment.write_bytes(b"patient alignment")
    view = tmp_path / "run" / "input.bam"
    view.parent.mkdir()
    binding = AlignmentBinding(str(alignment))
    binding.install_view(view)
    plan = AlignmentPlan(
        input_path=str(alignment),
        view_path=str(view),
        file_format="bam",
        index_path=f"{view}.bai",
        reference_path=None,
        reference_source="not-required",
        uncovered_contigs=(),
        unmapped_scan="indexed",
        binding=binding,
    )

    def fail_resolution(**kwargs: object) -> str:
        del kwargs
        raise RuntimeError("region resolution failed")

    try:
        with pytest.raises(RuntimeError, match="region resolution failed"):
            calculate_alignment_coverage(
                plan=plan,
                region=None,
                reference_assembly="hg19",
                threads=2,
                config={},
                output_dir=str(tmp_path / "coverage"),
                coverage_calculator=lambda **kwargs: None,
                region_resolver=fail_resolution,
            )

        assert binding.is_open is False
        assert not view.exists()
    finally:
        plan.close()


def test_coverage_failure_is_not_replaced_by_cleanup_failure(tmp_path: Path) -> None:
    """The stage diagnostic is primary; binding cleanup is secondary while it propagates."""

    class CleanupFailingPlan(AlignmentPlan):
        def close(self) -> None:
            raise RuntimeError("cleanup failure")

    plan = CleanupFailingPlan(
        input_path="/input/patient.bam",
        view_path="/run/input.bam",
        file_format="bam",
        index_path="/run/input.bam.bai",
        reference_path=None,
        reference_source="not-required",
        uncovered_contigs=(),
        unmapped_scan="indexed",
    )

    def fail_coverage(**kwargs: object) -> None:
        del kwargs
        raise ValueError("primary coverage failure")

    with pytest.raises(ValueError, match="primary coverage failure"):
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


def test_the_assembly_alias_is_resolved_before_the_geometry_is_looked_up():
    """`hg19` is the default, and the geometry is keyed `GRCh37` (#222).

    Found in adversarial review of the PR: looking the user-facing name up
    directly missed on the default and on six of the eight supported aliases, so
    the build-comparable columns wrote `NA` on almost every real run while every
    test passed by handing the assembly dict in directly.
    """
    import json
    from pathlib import Path

    from vntyper.scripts.alignment_contract import AlignmentPlan
    from vntyper.scripts.pipeline_coverage import calculate_alignment_coverage

    config = json.loads(Path("vntyper/config.json").read_text(encoding="utf-8"))
    seen = {}

    def recording_calculator(**kwargs):
        seen.update(kwargs)

    for alias in ("hg19", "hg38", "GRCh37", "GRCh38", "hg19_ensembl", "hg38_ensembl"):
        seen.clear()
        calculate_alignment_coverage(
            plan=AlignmentPlan(
                input_path="/data/s.bam",
                view_path="/data/s.bam",
                file_format="bam",
                index_path="/data/s.bam.bai",
                reference_path=None,
                reference_source="none",
                uncovered_contigs=(),
                unmapped_scan="none",
            ),
            region="chr1:155160500-155162000",
            reference_assembly=alias,
            threads=1,
            config=config,
            output_dir="/tmp",
            coverage_calculator=recording_calculator,
            region_resolver=lambda **_: "chr1:155160500-155162000",
        )
        assert seen["assembly_config"] is not None, f"{alias} resolved to no geometry"
        assert "vntr_array_coords" in seen["assembly_config"], alias
