"""Unit tests for focused alignment-to-FASTQ conversion orchestration."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vntyper.scripts import fastq_bam_processing
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.alignment_fastq_conversion import (
    plan_alignment_fastq_conversion,
    run_alignment_fastq_conversion,
)

pytestmark = pytest.mark.unit


def test_conversion_paths_stay_beside_the_output_for_a_hostile_name(tmp_path: Path) -> None:
    output = tmp_path / "patient's output; literal"

    paths = plan_alignment_fastq_conversion(output=output, output_name="sample $(literal)")

    assert paths.fastq_r1 == output / "sample $(literal)_R1.fastq.gz"
    assert paths.fastq_r2 == output / "sample $(literal)_R2.fastq.gz"
    assert paths.fastq_other == output / "sample $(literal)_other.fastq.gz"
    assert paths.fastq_single == output / "sample $(literal)_single.fastq.gz"
    assert paths.sort_tmp_prefix == output / "sample $(literal)_sort_tmp"
    assert paths.log_file == output / "sample $(literal)_sort_fastq.log"


def test_conversion_wires_exact_builder_arguments_and_critical_runner(tmp_path: Path) -> None:
    paths = plan_alignment_fastq_conversion(output=tmp_path, output_name="output")
    runner = MagicMock(return_value=True)

    with patch(
        "vntyper.scripts.alignment_fastq_conversion.build_bam_to_fastq_command",
        return_value="conversion command",
    ) as builder:
        result = run_alignment_fastq_conversion(
            paths=paths,
            final_bam=tmp_path / "output_sliced.bam",
            samtools_path="samtools",
            threads=4,
            command_runner=runner,
        )

    builder.assert_called_once_with(
        samtools_path="samtools",
        in_bam=tmp_path / "output_sliced.bam",
        threads=4,
        fastq_r1=paths.fastq_r1,
        fastq_r2=paths.fastq_r2,
        fastq_other=paths.fastq_other,
        fastq_single=paths.fastq_single,
        sort_tmp_prefix=paths.sort_tmp_prefix,
    )
    runner.assert_called_once_with("conversion command", str(paths.log_file), critical=True)
    assert result == tuple(
        str(path) for path in (paths.fastq_r1, paths.fastq_r2, paths.fastq_other, paths.fastq_single)
    )


def test_conversion_raises_the_existing_message_when_runner_returns_false(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    paths = plan_alignment_fastq_conversion(output=tmp_path, output_name="output")

    with pytest.raises(RuntimeError, match=r"^BAM to FASTQ conversion failed\.$"):
        run_alignment_fastq_conversion(
            paths=paths,
            final_bam=tmp_path / "output_sliced.bam",
            samtools_path="samtools",
            threads=4,
            command_runner=MagicMock(return_value=False),
        )

    assert "BAM to FASTQ conversion failed." in caplog.messages


def test_conversion_propagates_runner_exceptions(tmp_path: Path) -> None:
    paths = plan_alignment_fastq_conversion(output=tmp_path, output_name="output")
    failure = OSError("runner failed")

    with pytest.raises(OSError, match="runner failed"):
        run_alignment_fastq_conversion(
            paths=paths,
            final_bam=tmp_path / "output_sliced.bam",
            samtools_path="samtools",
            threads=4,
            command_runner=MagicMock(side_effect=failure),
        )


def _fast_plan() -> AlignmentPlan:
    return AlignmentPlan(
        input_path="/data/sample.bam",
        view_path="/data/sample.bam",
        file_format="bam",
        index_path="/data/sample.bam.bai",
        reference_path=None,
        reference_source="test",
        uncovered_contigs=(),
        unmapped_scan="indexed",
    )


@pytest.mark.parametrize(("delete_intermediates", "expected_exists"), [(False, True), (True, False)])
def test_process_delegates_with_caller_runner_and_cleans_after_success(
    tmp_path: Path,
    delete_intermediates: bool,
    expected_exists: bool,
) -> None:
    intermediate = tmp_path / "output_unmapped.bam"
    intermediate.write_bytes(b"intermediate")
    runner = MagicMock(return_value=True)
    expected = ("r1", "r2", "other", "single")

    with (
        patch.object(fastq_bam_processing, "run_command", runner),
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "publish_partial"),
        patch.object(fastq_bam_processing, "run_alignment_fastq_conversion", return_value=expected) as conversion,
    ):
        result = fastq_bam_processing.process_bam_to_fastq(
            output=tmp_path,
            output_name="output",
            threads=4,
            config={"tools": {"samtools": "samtools"}},
            plan=_fast_plan(),
            fast_mode=True,
            delete_intermediates=delete_intermediates,
        )

    conversion.assert_called_once()
    assert conversion.call_args.kwargs["command_runner"] is runner
    assert conversion.call_args.kwargs["paths"].sort_tmp_prefix == tmp_path / "output_sort_tmp"
    assert result == expected
    assert intermediate.exists() is expected_exists


def test_conversion_failure_propagates_before_cleanup(tmp_path: Path) -> None:
    intermediate = tmp_path / "output_unmapped.bam"
    intermediate.write_bytes(b"intermediate")
    failure = RuntimeError("conversion failed")

    with (
        patch.object(fastq_bam_processing, "run_command", return_value=True),
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "publish_partial"),
        patch.object(fastq_bam_processing, "run_alignment_fastq_conversion", side_effect=failure),
        pytest.raises(RuntimeError, match="conversion failed"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=tmp_path,
            output_name="output",
            threads=4,
            config={"tools": {"samtools": "samtools"}},
            plan=_fast_plan(),
            fast_mode=True,
            delete_intermediates=True,
        )

    assert intermediate.exists()
