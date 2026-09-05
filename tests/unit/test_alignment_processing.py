# tests/unit/test_alignment_processing.py

"""
Unit tests for alignment processing scripts.
Tests include functionality for checking BWA index completeness
and aligning/sorting FASTQ files.
"""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts import alignment_processing
from vntyper.scripts.alignment_processing import align_and_sort_fastq, check_bwa_index
from vntyper.scripts.alignment_target_io import BWA_INDEX_EXTENSIONS
from vntyper.scripts.artifact_publish import partial_path

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


def test_check_bwa_index_all_present(tmp_path, test_config, caplog):
    """
    Test that BWA index checking passes when all required files are present.
    """
    ref_path = tmp_path / "ref.fa"
    ref_path.touch()

    # Typically .amb, .ann, .bwt, .pac, .sa
    for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
        (tmp_path / f"ref.fa{ext}").touch()

    with caplog.at_level(logging.WARNING):
        result = check_bwa_index(ref_path)

    assert result is True, "check_bwa_index should return True when all index files exist."
    assert "Missing BWA index files" not in caplog.text


def test_check_bwa_index_uses_the_same_configured_suffixes_as_pipeline_safety(tmp_path):
    ref_path = tmp_path / "ref.fa"
    ref_path.touch()
    (tmp_path / "ref.fa.custom-index").touch()
    config = {"tool_params": {"bwa_index_extensions": [".custom-index"]}}

    assert check_bwa_index(ref_path, config) is False
    for extension in BWA_INDEX_EXTENSIONS:
        (tmp_path / f"ref.fa{extension}").touch()
    assert check_bwa_index(ref_path, config) is True


def test_align_and_sort_fastq_missing_tools(tmp_path, test_config):
    """
    Test alignment fails when required tools are not configured.
    """
    config_for_test = {"tools": {}}  # No 'samtools' or 'bwa' keys

    output_dir = tmp_path / "align_output"
    output_dir.mkdir()

    # We still pick up file paths from the test config if needed:
    fastq1 = test_config["file_resources"][5]["filename"]
    fastq2 = test_config["file_resources"][6]["filename"]
    reference = test_config["file_resources"][0]["filename"]

    result = align_and_sort_fastq(
        fastq1=fastq1,
        fastq2=fastq2,
        reference=Path(reference),
        output_dir=output_dir,
        output_name="test",
        threads=4,
        config=config_for_test,
    )
    assert result is None, (
        "Expected None when tools are missing, because 'samtools'/'bwa' are not configured in config_for_test['tools']."
    )


def test_align_and_sort_fastq_failure_leaves_no_partial_or_final_files(tmp_path: Path) -> None:
    """A failed BWA/sort command leaves neither .partial nor final BAM."""
    out_dir = tmp_path / "alignment"
    sorted_bam = out_dir / "sample_sorted.bam"
    partial_bam = partial_path(sorted_bam)

    with (
        patch.object(alignment_processing, "check_bwa_index", return_value=True),
        patch.object(alignment_processing, "run_command", return_value=False),
    ):
        result = alignment_processing.align_and_sort_fastq(
            fastq1="r1.fq",
            fastq2="r2.fq",
            output_dir=out_dir,
            output_name="sample",
            reference="ref.fa",
            threads=4,
            config={"tools": {"samtools": "samtools", "bwa": "bwa"}},
        )

    assert result is None
    assert not sorted_bam.exists()
    assert not partial_bam.exists()


def test_align_and_sort_fastq_index_failure_leaves_no_partial_index(tmp_path: Path) -> None:
    """A failed index command leaves no partial index file."""
    out_dir = tmp_path / "alignment"
    sorted_bam = out_dir / "sample_sorted.bam"
    final_bai = sorted_bam.with_suffix(".bam.bai")
    partial_bai = partial_path(final_bai)

    def mock_run(cmd, log, critical=True):
        if "samtools sort" in cmd:
            partial_path(sorted_bam).write_bytes(b"BAM_DATA")
            return True
        return False

    with (
        patch.object(alignment_processing, "check_bwa_index", return_value=True),
        patch.object(alignment_processing, "run_command", side_effect=mock_run),
    ):
        result = alignment_processing.align_and_sort_fastq(
            fastq1="r1.fq",
            fastq2="r2.fq",
            output_dir=out_dir,
            output_name="sample",
            reference="ref.fa",
            threads=4,
            config={"tools": {"samtools": "samtools", "bwa": "bwa"}},
        )

    assert result is None
    assert sorted_bam.exists()  # BAM succeeded and was published
    assert not final_bai.exists()
    assert not partial_bai.exists()

