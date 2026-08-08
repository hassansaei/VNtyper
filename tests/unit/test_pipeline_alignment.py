"""Pure unit tests for pipeline alignment target and preflight preparation."""

from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts.pipeline_alignment import (
    build_alignment_preflight_kwargs,
    format_regions_as_bed,
    prepare_alignment_target,
)

pytestmark = pytest.mark.unit


def test_format_regions_as_bed_preserves_the_slice_shape() -> None:
    assert format_regions_as_bed("chr1:10-20, chr2:30-40") == "chr1\t10\t20\nchr2\t30\t40\n"


def test_format_regions_as_bed_rejects_a_malformed_region() -> None:
    with pytest.raises(ValueError, match="Invalid region format: chr1-10-20"):
        format_regions_as_bed("chr1-10-20")


def test_prepare_alignment_target_keeps_a_provided_bed_byte_identical(tmp_path: Path) -> None:
    provided = tmp_path / "provided.bed"
    original = b"chr1\t10\t20\n# operator annotation\n"
    provided.write_bytes(original)

    result = prepare_alignment_target(
        input_type="BAM",
        bam="/input/sample.bam",
        cram=None,
        output_dir=tmp_path / "output",
        reference_assembly="hg19",
        config={},
        bed_file=provided,
        custom_regions=None,
    )

    assert result == provided
    assert provided.read_bytes() == original


def test_prepare_alignment_target_rejects_a_missing_provided_bed(tmp_path: Path) -> None:
    missing = tmp_path / "missing.bed"

    with pytest.raises(FileNotFoundError, match=f"BED file not found: {missing}"):
        prepare_alignment_target(
            input_type="CRAM",
            bam=None,
            cram="/input/sample.cram",
            output_dir=tmp_path / "output",
            reference_assembly="hg38",
            config={},
            bed_file=missing,
            custom_regions=None,
        )


def test_prepare_alignment_target_resolves_alignment_region_before_writing(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()

    with mock.patch(
        "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
        return_value="chr1:101-202",
    ) as get_region:
        result = prepare_alignment_target(
            input_type="CRAM",
            bam=None,
            cram="/input/sample.cram",
            output_dir=output,
            reference_assembly="hg38",
            config={},
            bed_file=None,
            custom_regions=None,
        )

    assert result == output / "predefined_regions_hg38.bed"
    assert result.read_text(encoding="utf-8") == "chr1\t101\t202\n"
    get_region.assert_called_once_with(
        bam_file="/input/sample.cram",
        reference_assembly="hg38",
        region_type="bam_region",
        config={},
    )


def test_prepare_alignment_target_rejects_a_missing_selected_alignment_path(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="CRAM input path is required for target resolution"):
        prepare_alignment_target(
            input_type="CRAM",
            bam=None,
            cram=None,
            output_dir=tmp_path,
            reference_assembly="hg38",
            config={},
            bed_file=None,
            custom_regions=None,
        )


def test_prepare_alignment_target_uses_get_defaults_for_fastq_config(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()

    with pytest.raises(ValueError, match="Missing configuration for region: bam_region_hg19"):
        prepare_alignment_target(
            input_type="FASTQ",
            bam=None,
            cram=None,
            output_dir=output,
            reference_assembly="hg19",
            config={},
            bed_file=None,
            custom_regions=None,
        )


@pytest.mark.parametrize("alignment_header", [None, "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n"])
def test_build_alignment_preflight_kwargs_pins_exact_bed_and_header_contigs(
    tmp_path: Path, alignment_header: str | None
) -> None:
    bed = tmp_path / "target.bed"

    result = build_alignment_preflight_kwargs(
        in_path="/input/sample.cram",
        output_dir=tmp_path / "stage",
        output_name="input",
        file_format="cram",
        config={},
        threads=7,
        bed_file=bed,
        reference_assembly="hg38",
        fast_mode=True,
        alignment_header=alignment_header,
    )

    assert result == {
        "in_path": "/input/sample.cram",
        "output_dir": str(tmp_path / "stage"),
        "output_name": "input",
        "file_format": "cram",
        "config": {},
        "threads": 7,
        "region": None,
        "bed_file": bed,
        "reference_assembly": "hg38",
        "reference_fasta": None,
        "header_contigs": () if alignment_header is None else ("chr1", "chr2"),
        "m5": None,
        "fast_mode": True,
    }
