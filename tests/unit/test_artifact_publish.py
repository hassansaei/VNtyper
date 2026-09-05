"""Unit tests for atomic artifact publishing and cleanup (#314)."""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vntyper.scripts import alignment_processing, fastq_bam_processing, kestrel_genotyping
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.artifact_publish import (
    PARTIAL_SUFFIX,
    discard_partial,
    partial_output,
    partial_path,
    publish_partial,
)

pytestmark = pytest.mark.unit


def test_partial_path_is_deterministic_and_sibling(tmp_path: Path) -> None:
    final = tmp_path / "output_sliced.bam"
    partial = partial_path(final)
    assert partial == tmp_path / f"output_sliced.bam{PARTIAL_SUFFIX}"
    assert partial.parent == final.parent
    assert partial_path(str(final)) == partial


def test_publish_partial_replaces_atomically(tmp_path: Path) -> None:
    final = tmp_path / "output_sorted.bam"
    partial = partial_path(final)
    partial.write_bytes(b"BAM_DATA")

    published = publish_partial(partial, final)
    assert published == final
    assert final.read_bytes() == b"BAM_DATA"
    assert not partial.exists()


def test_publish_partial_raises_on_missing_source(tmp_path: Path) -> None:
    final = tmp_path / "target.bam"
    partial = partial_path(final)
    with pytest.raises(FileNotFoundError, match="Partial artifact not found"):
        publish_partial(partial, final)


def test_discard_partial_removes_existing_file_and_symlink(tmp_path: Path) -> None:
    partial = tmp_path / "file.bam.partial"
    partial.write_bytes(b"temp")
    discard_partial(partial)
    assert not partial.exists()

    # Repeated discard is a no-op
    discard_partial(partial)
    discard_partial(None)

    # Discarding a broken symlink unlinks the symlink itself
    target = tmp_path / "nonexistent"
    symlink = tmp_path / "symlink.bam.partial"
    symlink.symlink_to(target)
    assert symlink.is_symlink()
    discard_partial(symlink)
    assert not symlink.is_symlink()


def test_partial_output_publishes_on_clean_exit(tmp_path: Path) -> None:
    final = tmp_path / "result.bam"
    with partial_output(final) as partial:
        assert partial == partial_path(final)
        partial.write_bytes(b"VALID_BAM")

    assert final.exists()
    assert final.read_bytes() == b"VALID_BAM"
    assert not partial.exists()


def test_partial_output_discards_on_exception(tmp_path: Path) -> None:
    final = tmp_path / "result.bam"
    with pytest.raises(RuntimeError, match="stage exploded"), partial_output(final) as partial:
        partial.write_bytes(b"CORRUPT_OR_PARTIAL")
        raise RuntimeError("stage exploded")

    assert not final.exists()
    assert not partial.exists()


def test_partial_output_raises_when_partial_was_not_created(tmp_path: Path) -> None:
    final = tmp_path / "missing.bam"
    with pytest.raises(FileNotFoundError, match="Expected partial artifact not created"), partial_output(final):
        pass

    assert not final.exists()


def test_partial_output_checks_non_empty_when_requested(tmp_path: Path) -> None:
    final = tmp_path / "empty.bai"
    with (
        pytest.raises(OSError, match="Generated artifact is empty"),
        partial_output(final, check_non_empty=True) as partial,
    ):
        partial.touch()  # 0 bytes

    assert not final.exists()
    assert not partial.exists()


# --------------------------------------------------------------------------------------
# Pipeline stage failure and cleanup tests (#314)
# --------------------------------------------------------------------------------------


def _plan(tmp_path: Path, file_format: str = "bam") -> AlignmentPlan:
    in_dir = tmp_path / "in"
    in_dir.mkdir(parents=True, exist_ok=True)
    in_file = in_dir / f"input.{file_format}"
    in_file.write_bytes(b"DATA")
    index = in_dir / f"input.{file_format}.bai"
    index.write_bytes(b"INDEX")
    return AlignmentPlan(
        input_path=str(in_file),
        view_path=str(in_file),
        file_format=file_format,
        index_path=str(index),
        reference_path=str(tmp_path / "ref.fa") if file_format == "cram" else None,
        reference_source="test",
        uncovered_contigs=(),
        unmapped_scan="indexed",
    )


def test_process_bam_to_fastq_slice_failure_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"

    def fail_slice(cmd, log, **kwargs):
        partial_slice.parent.mkdir(parents=True, exist_ok=True)
        partial_slice.write_bytes(b"partial_data")
        return False

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=fail_slice),
        pytest.raises(RuntimeError, match="region slicing failed"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=True,
        )

    assert not partial_slice.exists()


def test_process_bam_to_fastq_slice_exception_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"

    def explode_slice(cmd, log, **kwargs):
        partial_slice.parent.mkdir(parents=True, exist_ok=True)
        partial_slice.write_bytes(b"partial_data")
        raise RuntimeError("boom")

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=explode_slice),
        pytest.raises(RuntimeError, match="boom"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=True,
        )

    assert not partial_slice.exists()


def test_process_bam_to_fastq_unmapped_failure_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"
    partial_unmapped = tmp_path / "run" / "output_unmapped.bam.partial"

    def run_cmd(cmd, log, **kwargs):
        if "-F 4" in cmd:
            partial_slice.parent.mkdir(parents=True, exist_ok=True)
            partial_slice.write_bytes(b"slice_data")
            return True
        if "-f 4" in cmd:
            partial_unmapped.parent.mkdir(parents=True, exist_ok=True)
            partial_unmapped.write_bytes(b"unmapped_data")
            return False
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="filtering failed"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=False,
        )

    assert not partial_unmapped.exists()


def test_process_bam_to_fastq_unmapped_exception_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"
    partial_unmapped = tmp_path / "run" / "output_unmapped.bam.partial"

    def run_cmd(cmd, log, **kwargs):
        if "-F 4" in cmd:
            partial_slice.parent.mkdir(parents=True, exist_ok=True)
            partial_slice.write_bytes(b"slice_data")
            return True
        if "-f 4" in cmd:
            partial_unmapped.parent.mkdir(parents=True, exist_ok=True)
            partial_unmapped.write_bytes(b"unmapped_data")
            raise RuntimeError("filter boom")
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="filter boom"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=False,
        )

    assert not partial_unmapped.exists()


def test_process_bam_to_fastq_merge_failure_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"
    partial_merged = tmp_path / "run" / "output_sliced_unmapped.bam.partial"

    def run_cmd(cmd, log, **kwargs):
        if "-F 4" in cmd:
            partial_slice.parent.mkdir(parents=True, exist_ok=True)
            partial_slice.write_bytes(b"slice_data")
            return True
        if "-f 4" in cmd:
            unmapped = tmp_path / "run" / "output_unmapped.bam.partial"
            unmapped.parent.mkdir(parents=True, exist_ok=True)
            unmapped.write_bytes(b"unmapped_data")
            return True
        if "samtools merge" in cmd:
            partial_merged.parent.mkdir(parents=True, exist_ok=True)
            partial_merged.write_bytes(b"merged_data")
            return False
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="merging failed"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=False,
        )

    assert not partial_merged.exists()
    assert not partial_slice.exists()


def test_process_bam_to_fastq_merge_exception_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"
    partial_merged = tmp_path / "run" / "output_sliced_unmapped.bam.partial"

    def run_cmd(cmd, log, **kwargs):
        if "-F 4" in cmd:
            partial_slice.parent.mkdir(parents=True, exist_ok=True)
            partial_slice.write_bytes(b"slice_data")
            return True
        if "-f 4" in cmd:
            unmapped = tmp_path / "run" / "output_unmapped.bam.partial"
            unmapped.parent.mkdir(parents=True, exist_ok=True)
            unmapped.write_bytes(b"unmapped_data")
            return True
        if "samtools merge" in cmd:
            partial_merged.parent.mkdir(parents=True, exist_ok=True)
            partial_merged.write_bytes(b"merged_data")
            raise RuntimeError("merge boom")
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="merge boom"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=False,
        )

    assert not partial_merged.exists()
    assert not partial_slice.exists()


def test_process_bam_to_fastq_indexing_failure_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"
    partial_bai = tmp_path / "run" / "output_sliced.bam.bai.partial"

    def run_cmd(cmd, log, **kwargs):
        if "-F 4" in cmd or "-P -b" in cmd:
            partial_slice.parent.mkdir(parents=True, exist_ok=True)
            partial_slice.write_bytes(b"slice_data")
            return True
        if "samtools index" in cmd:
            partial_bai.parent.mkdir(parents=True, exist_ok=True)
            partial_bai.write_bytes(b"bad_index")
            return False
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="Re-indexing BAM file failed"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=True,
            needs_advntr=True,
        )

    assert not partial_bai.exists()


def test_process_bam_to_fastq_indexing_empty_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"
    partial_bai = tmp_path / "run" / "output_sliced.bam.bai.partial"

    def run_cmd(cmd, log, **kwargs):
        if "-F 4" in cmd or "-P -b" in cmd:
            partial_slice.parent.mkdir(parents=True, exist_ok=True)
            partial_slice.write_bytes(b"slice_data")
            return True
        if "samtools index" in cmd:
            partial_bai.parent.mkdir(parents=True, exist_ok=True)
            partial_bai.touch()  # 0 bytes
            return True
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="empty"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=True,
            needs_advntr=True,
        )

    assert not partial_bai.exists()


def test_process_bam_to_fastq_indexing_exception_discards_partial(tmp_path: Path) -> None:
    plan = _plan(tmp_path)
    partial_slice = tmp_path / "run" / "output_sliced.bam.partial"
    partial_bai = tmp_path / "run" / "output_sliced.bam.bai.partial"

    def run_cmd(cmd, log, **kwargs):
        if "-F 4" in cmd or "-P -b" in cmd:
            partial_slice.parent.mkdir(parents=True, exist_ok=True)
            partial_slice.write_bytes(b"slice_data")
            return True
        if "samtools index" in cmd:
            partial_bai.parent.mkdir(parents=True, exist_ok=True)
            partial_bai.write_bytes(b"index_data")
            raise RuntimeError("index boom")
        return True

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="index boom"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path / "run"),
            output_name="output",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=True,
            needs_advntr=True,
        )

    assert not partial_bai.exists()


def test_align_and_sort_fastq_alignment_exception_discards_partial(tmp_path: Path) -> None:
    out_dir = tmp_path / "alignment"
    sorted_bam = out_dir / "sample_sorted.bam"
    partial_bam = partial_path(sorted_bam)

    def explode(cmd, log, **kwargs):
        partial_bam.parent.mkdir(parents=True, exist_ok=True)
        partial_bam.write_bytes(b"partial_bam")
        raise RuntimeError("bwa boom")

    with (
        patch.object(alignment_processing, "check_bwa_index", return_value=True),
        patch.object(alignment_processing, "run_command", side_effect=explode),
        pytest.raises(RuntimeError, match="bwa boom"),
    ):
        alignment_processing.align_and_sort_fastq(
            fastq1=Path("r1.fq"),
            fastq2=Path("r2.fq"),
            output_dir=out_dir,
            output_name="sample",
            reference=Path("ref.fa"),
            threads=4,
            config={"tools": {"samtools": "samtools", "bwa": "bwa"}},
        )

    assert not partial_bam.exists()


def test_align_and_sort_fastq_indexing_exception_discards_partial(tmp_path: Path) -> None:
    out_dir = tmp_path / "alignment"
    sorted_bam = out_dir / "sample_sorted.bam"
    final_bai = sorted_bam.with_suffix(".bam.bai")
    partial_bai = partial_path(final_bai)

    def run_cmd(cmd, log, **kwargs):
        if "samtools sort" in cmd:
            partial_path(sorted_bam).parent.mkdir(parents=True, exist_ok=True)
            partial_path(sorted_bam).write_bytes(b"BAM_DATA")
            return True
        if "samtools index" in cmd:
            partial_bai.parent.mkdir(parents=True, exist_ok=True)
            partial_bai.write_bytes(b"partial_index")
            raise RuntimeError("index boom")
        return True

    with (
        patch.object(alignment_processing, "check_bwa_index", return_value=True),
        patch.object(alignment_processing, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="index boom"),
    ):
        alignment_processing.align_and_sort_fastq(
            fastq1=Path("r1.fq"),
            fastq2=Path("r2.fq"),
            output_dir=out_dir,
            output_name="sample",
            reference=Path("ref.fa"),
            threads=4,
            config={"tools": {"samtools": "samtools", "bwa": "bwa"}},
        )

    assert not partial_bai.exists()


def test_align_and_sort_fastq_indexing_empty_discards_partial(tmp_path: Path) -> None:
    out_dir = tmp_path / "alignment"
    sorted_bam = out_dir / "sample_sorted.bam"
    final_bai = sorted_bam.with_suffix(".bam.bai")
    partial_bai = partial_path(final_bai)

    def run_cmd(cmd, log, **kwargs):
        if "samtools sort" in cmd:
            partial_path(sorted_bam).parent.mkdir(parents=True, exist_ok=True)
            partial_path(sorted_bam).write_bytes(b"BAM_DATA")
            return True
        if "samtools index" in cmd:
            partial_bai.parent.mkdir(parents=True, exist_ok=True)
            partial_bai.touch()  # empty
            return True
        return True

    with (
        patch.object(alignment_processing, "check_bwa_index", return_value=True),
        patch.object(alignment_processing, "run_command", side_effect=run_cmd),
    ):
        result = alignment_processing.align_and_sort_fastq(
            fastq1=Path("r1.fq"),
            fastq2=Path("r2.fq"),
            output_dir=out_dir,
            output_name="sample",
            reference=Path("ref.fa"),
            threads=4,
            config={"tools": {"samtools": "samtools", "bwa": "bwa"}},
        )

    assert result is None
    assert not partial_bai.exists()


def test_downsample_bam_view_failure_discards_partial(tmp_path: Path) -> None:
    in_bam = tmp_path / "input.bam"
    in_bam.write_bytes(b"BAM")
    partial_down = tmp_path / "input_downsampled.bam.partial"

    def fail_view(cmd, check=True):
        partial_down.write_bytes(b"corrupt")
        raise subprocess.CalledProcessError(1, cmd)

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "calculate_vntr_coverage", return_value={"mean": 400.0}),
        patch("subprocess.run", side_effect=fail_view),
    ):
        res = fastq_bam_processing.downsample_bam_if_needed(
            bam_path=in_bam,
            max_coverage=100,
            reference_assembly="hg19",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            coverage_dir=str(tmp_path),
            coverage_prefix="cov",
        )

    assert res == in_bam
    assert not partial_down.exists()


def test_downsample_bam_indexing_empty_discards_partial(tmp_path: Path) -> None:
    in_bam = tmp_path / "input.bam"
    in_bam.write_bytes(b"BAM")
    down_bam = tmp_path / "input_downsampled.bam"
    sorted_down_bam = down_bam.with_suffix(".sorted.bam")
    partial_sorted = partial_path(sorted_down_bam)
    final_bai = Path(f"{sorted_down_bam}.bai")
    partial_bai = partial_path(final_bai)

    def run_subp(cmd, check=True):
        if "view" in cmd:
            partial_path(down_bam).write_bytes(b"down")
        elif "sort" in cmd:
            partial_sorted.write_bytes(b"sorted")
        elif "index" in cmd:
            partial_bai.touch()  # empty
        return MagicMock(returncode=0)

    with (
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "calculate_vntr_coverage", return_value={"mean": 400.0}),
        patch("subprocess.run", side_effect=run_subp),
    ):
        res = fastq_bam_processing.downsample_bam_if_needed(
            bam_path=in_bam,
            max_coverage=100,
            reference_assembly="hg19",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            coverage_dir=str(tmp_path),
            coverage_prefix="cov",
        )

    assert res == in_bam
    assert not partial_sorted.exists()
    assert not partial_bai.exists()


def test_kestrel_convert_sam_empty_index_discards_partial(tmp_path: Path) -> None:
    sam = tmp_path / "output.sam"
    sam.write_text("@HD\tVN:1.6\n", encoding="utf-8")
    partial_bam = tmp_path / "output.bam.partial"
    partial_bai = tmp_path / "output.bam.bai.partial"

    def run_cmd(cmd, log_file=None, **kwargs):
        if not cmd.startswith("samtools index"):
            partial_bam.write_bytes(b"BAM")
        else:
            partial_bai.touch()  # empty
        return True

    with (
        patch.object(kestrel_genotyping, "run_command", side_effect=run_cmd),
        pytest.raises(RuntimeError, match="empty"),
    ):
        kestrel_genotyping.convert_sam_to_bam_and_index(str(sam), str(tmp_path))

    assert not partial_bai.exists()
