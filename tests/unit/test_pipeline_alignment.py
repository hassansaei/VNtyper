"""Pure unit tests for pipeline alignment target and preflight preparation."""

import hashlib
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts.alignment_preflight import resolve_reference
from vntyper.scripts.pipeline_alignment import (
    build_alignment_preflight_kwargs,
    format_regions_as_bed,
    prepare_alignment_target,
)

pytestmark = pytest.mark.unit


def _tree_digest(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


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


def test_generated_bed_rejects_a_symlink_to_the_patient_alignment_without_writing_it(tmp_path: Path) -> None:
    input_root = tmp_path / "input"
    output = tmp_path / "output"
    input_root.mkdir()
    output.mkdir()
    alignment = input_root / "patient.bam"
    alignment.write_bytes(b"patient-alignment")
    (output / "predefined_regions_hg19.bed").symlink_to(alignment)
    before = _tree_digest(input_root)
    error: ValueError | None = None

    with mock.patch(
        "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
        return_value="chr1:10-20",
    ):
        try:
            prepare_alignment_target(
                input_type="BAM",
                bam=str(alignment),
                cram=None,
                output_dir=output,
                reference_assembly="hg19",
                config={},
                bed_file=None,
                custom_regions=None,
            )
        except ValueError as caught:
            error = caught

    assert _tree_digest(input_root) == before
    assert error is not None and "symlink" in str(error).lower()


def test_generated_bed_rejects_the_patient_input_tree_as_its_output_root(tmp_path: Path) -> None:
    input_root = tmp_path / "input"
    input_root.mkdir()
    alignment = input_root / "patient.cram"
    alignment.write_bytes(b"patient-cram")
    before = _tree_digest(input_root)
    error: ValueError | None = None

    with mock.patch(
        "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
        return_value="chr1:10-20",
    ):
        try:
            prepare_alignment_target(
                input_type="CRAM",
                bam=None,
                cram=str(alignment),
                output_dir=input_root,
                reference_assembly="hg19",
                config={},
                bed_file=None,
                custom_regions=None,
            )
        except ValueError as caught:
            error = caught

    assert _tree_digest(input_root) == before
    assert error is not None and "input tree" in str(error).lower()


def test_generated_bed_rejects_a_hardlink_to_the_patient_source_index(tmp_path: Path) -> None:
    input_root = tmp_path / "input"
    output = tmp_path / "output"
    input_root.mkdir()
    output.mkdir()
    alignment = input_root / "patient.bam"
    alignment.write_bytes(b"patient-alignment")
    source_index = input_root / "patient.bam.bai"
    source_index.write_bytes(b"patient-index")
    destination = output / "predefined_regions_hg19.bed"
    destination.hardlink_to(source_index)
    before = _tree_digest(input_root)
    error: ValueError | None = None

    with mock.patch(
        "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
        return_value="chr1:10-20",
    ):
        try:
            prepare_alignment_target(
                input_type="BAM",
                bam=str(alignment),
                cram=None,
                output_dir=output,
                reference_assembly="hg19",
                config={},
                bed_file=None,
                custom_regions=None,
            )
        except ValueError as caught:
            error = caught

    assert _tree_digest(input_root) == before
    assert error is not None and "protected" in str(error).lower()


def test_generated_bed_rejects_an_existing_non_regular_destination(tmp_path: Path) -> None:
    destination = tmp_path / "custom_regions.bed"
    destination.mkdir()

    with pytest.raises(ValueError, match="regular file"):
        prepare_alignment_target(
            input_type="FASTQ",
            bam=None,
            cram=None,
            output_dir=tmp_path,
            reference_assembly="hg19",
            config={},
            bed_file=None,
            custom_regions="chr1:10-20",
        )


def test_generated_bed_regular_rerun_is_an_atomic_replacement(tmp_path: Path) -> None:
    destination = tmp_path / "custom_regions.bed"
    destination.write_text("stale\n", encoding="utf-8")
    old_inode = destination.stat().st_ino

    result = prepare_alignment_target(
        input_type="FASTQ",
        bam=None,
        cram=None,
        output_dir=tmp_path,
        reference_assembly="hg19",
        config={},
        bed_file=None,
        custom_regions="chr2:30-40",
    )

    assert result == destination
    assert destination.read_text(encoding="utf-8") == "chr2\t30\t40\n"
    assert destination.stat().st_ino != old_inode


def test_generated_bed_failed_atomic_install_preserves_the_regular_rerun(tmp_path: Path) -> None:
    destination = tmp_path / "custom_regions.bed"
    destination.write_bytes(b"previous-good-target\n")

    with (
        mock.patch("os.replace", side_effect=OSError("install failed")),
        pytest.raises(OSError, match="install failed"),
    ):
        prepare_alignment_target(
            input_type="FASTQ",
            bam=None,
            cram=None,
            output_dir=tmp_path,
            reference_assembly="hg19",
            config={},
            bed_file=None,
            custom_regions="chr2:30-40",
        )

    assert destination.read_bytes() == b"previous-good-target\n"


@pytest.mark.parametrize("alignment_header", [None, "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n"])
def test_build_alignment_preflight_kwargs_pins_exact_bed_and_header_contigs(
    tmp_path: Path, alignment_header: str | None
) -> None:
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t10\t20\n", encoding="utf-8")

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


def test_preflight_kwargs_select_the_m5_for_the_first_active_bed_target(tmp_path: Path) -> None:
    bed = tmp_path / "target.bed"
    bed.write_text("\n# first real target follows\nchr2\t30\t40\nchr1\t10\t20\n", encoding="utf-8")
    header = "@SQ\tSN:chr1\tLN:100\tM5:first-checksum\n@SQ\tM5:target-checksum\tLN:200\tSN:chr2\n"

    result = build_alignment_preflight_kwargs(
        in_path="/input/sample.cram",
        output_dir=tmp_path / "stage",
        output_name="input",
        file_format="cram",
        config={},
        threads=2,
        bed_file=bed,
        reference_assembly="hg38",
        fast_mode=False,
        alignment_header=header,
    )

    assert result["header_contigs"] == ("chr1", "chr2")
    assert result["m5"] == "target-checksum"


def test_terminal_reference_diagnostic_uses_the_bed_target_and_its_m5(tmp_path: Path) -> None:
    bed = tmp_path / "target.bed"
    bed.write_text("# ignored\n\nchr2\t30\t40\n", encoding="utf-8")
    kwargs = build_alignment_preflight_kwargs(
        in_path="/input/sample.cram",
        output_dir=tmp_path,
        output_name="input",
        file_format="cram",
        config={},
        threads=2,
        bed_file=bed,
        reference_assembly="hg38",
        fast_mode=False,
        alignment_header=("@SQ\tSN:chr1\tLN:100\tM5:first-checksum\n@SQ\tSN:chr2\tLN:200\tM5:target-checksum\n"),
    )

    with (
        mock.patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(False, "decode failed")),
        pytest.raises(ValueError) as error,
    ):
        resolve_reference(
            kwargs["in_path"],
            (),
            kwargs["region"],
            kwargs["bed_file"],
            kwargs["config"],
            kwargs["threads"],
            kwargs["output_dir"],
            kwargs["output_name"],
            kwargs["header_contigs"],
            kwargs["m5"],
        )

    assert "contig=chr2" in str(error.value)
    assert "M5=target-checksum" in str(error.value)
