#!/usr/bin/env python3
# tests/unit/test_alignment_processing.py

"""
Unit tests for alignment processing scripts.
Tests include functionality for checking BWA index completeness
and aligning/sorting FASTQ files.
"""

import pytest
import logging
from pathlib import Path
from vntyper.scripts.alignment_processing import check_bwa_index, align_and_sort_fastq


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

    assert (
        result is True
    ), "check_bwa_index should return True when all index files exist."
    assert "Missing BWA index files" not in caplog.text


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
        "Expected None when tools are missing, because 'samtools'/'bwa' "
        "are not configured in config_for_test['tools']."
    )
