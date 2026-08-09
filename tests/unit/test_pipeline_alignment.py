"""Pure unit tests for pipeline alignment target and preflight preparation."""

import hashlib
import json
import os
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.alignment_preflight import resolve_reference
from vntyper.scripts.pipeline_alignment import (
    build_alignment_preflight_kwargs,
    format_regions_as_bed,
    prepare_alignment_target,
    prepare_input_alignment_preflight,
)
from vntyper.scripts.preflight_error_io import PreflightErrorContext
from vntyper.scripts.reference_resolution_environment import restore_reference_resolution

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
    output = tmp_path / "output"
    output.mkdir()
    original = b"chr1\t10\t20\n# operator annotation\n"
    provided.write_bytes(original)

    result = prepare_alignment_target(
        input_type="BAM",
        bam="/input/sample.bam",
        cram=None,
        output_dir=output,
        reference_assembly="hg19",
        config={},
        bed_file=provided,
        custom_regions=None,
    )

    assert result == output / "operator_regions.bed"
    assert result.read_bytes() == original
    assert provided.read_bytes() == original


def test_operator_bed_is_materialized_once_and_rerun_refreshes_only_the_owned_copy(tmp_path: Path) -> None:
    """Replacing the operator path after target preparation must not retarget later consumers."""
    output = tmp_path / "output"
    output.mkdir()
    provided = tmp_path / "provided.bed"
    first_bytes = "chr1\t10\t20\r\n# café\r\n".encode()
    second_bytes = b"chr2\t30\t40\n"
    provided.write_bytes(first_bytes)

    first_target = prepare_alignment_target(
        input_type="BAM",
        bam="/input/sample.bam",
        cram=None,
        output_dir=output,
        reference_assembly="hg19",
        config={},
        bed_file=provided,
        custom_regions=None,
    )
    replacement = tmp_path / "replacement.bed"
    replacement.write_bytes(second_bytes)
    replacement.replace(provided)

    assert first_target != provided
    assert first_target.read_bytes() == first_bytes
    assert provided.read_bytes() == second_bytes

    second_target = prepare_alignment_target(
        input_type="BAM",
        bam="/input/sample.bam",
        cram=None,
        output_dir=output,
        reference_assembly="hg19",
        config={},
        bed_file=provided,
        custom_regions=None,
    )

    assert second_target == first_target
    assert second_target.read_bytes() == second_bytes
    assert not tuple(output.glob(".*.tmp"))


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


def test_generated_bed_rejects_an_output_directory_nested_in_the_patient_tree(tmp_path: Path) -> None:
    """Containment, not only equality, would let pipeline artifacts mutate input."""
    input_root = tmp_path / "input"
    nested_output = input_root / "results"
    input_root.mkdir()
    nested_output.mkdir()
    alignment = input_root / "patient.cram"
    alignment.write_bytes(b"patient-cram")
    before = _tree_digest(input_root)

    with (
        mock.patch(
            "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
            return_value="chr1:10-20",
        ),
        pytest.raises(ValueError, match="input tree"),
    ):
        prepare_alignment_target(
            input_type="CRAM",
            bam=None,
            cram=str(alignment),
            output_dir=nested_output,
            reference_assembly="hg19",
            config={},
            bed_file=None,
            custom_regions=None,
        )

    assert _tree_digest(input_root) == before


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
    reference = tmp_path / "full reference.fa"
    error_output = tmp_path / "run-output"

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
        reference_fasta=reference,
        error_output_dir=error_output,
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
        "coverage_region": None,
        "reference_assembly": "hg38",
        "reference_fasta": str(reference),
        "header_reference_paths": (),
        "header_references": (),
        "has_remote_header_reference": False,
        "header_contigs": () if alignment_header is None else ("chr1", "chr2"),
        "m5": None,
        "header_m5s": (),
        "fast_mode": True,
        "error_output_dir": str(error_output),
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


def test_input_alignment_preflight_rejects_a_non_alignment_type_before_pinning(tmp_path: Path) -> None:
    """The owned seam cannot silently reinterpret a FASTQ path as an alignment."""
    with (
        mock.patch("vntyper.scripts.pipeline_alignment.pin_reference_resolution", return_value=None) as pin,
        pytest.raises(ValueError, match="requires BAM or CRAM input, got: FASTQ"),
    ):
        prepare_input_alignment_preflight(
            in_path=tmp_path / "reads.fastq.gz",
            input_type="FASTQ",
            output_dir=tmp_path / "output",
            config={},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    pin.assert_not_called()


def test_input_alignment_binding_precedes_validation_and_survives_source_replacement(tmp_path: Path) -> None:
    """Every evidence reader must consume the inode opened before quickcheck, not its mutable source name."""
    output = tmp_path / "run-output"
    output.mkdir()
    source = tmp_path / "patient-input" / "sample.bam"
    source.parent.mkdir()
    source.write_bytes(b"alignment-A")
    replacement = tmp_path / "replacement.bam"
    replacement.write_bytes(b"alignment-B")
    seen_view: list[str] = []

    def validate_bound(path: str, *, cwd: str | None, log_dir: str | None) -> None:
        assert cwd == "/project"
        assert log_dir == str(output)
        assert Path(path).read_bytes() == b"alignment-A"
        seen_view.append(path)
        replacement.replace(source)

    def read_bound_header(path: str, config: dict) -> str:
        del config
        assert path == seen_view[0]
        assert Path(path).read_bytes() == b"alignment-A"
        return "@SQ\tSN:chr1\tLN:249250621\n"

    def resolve_bound_region(*args: object, **kwargs: object) -> str:
        del args
        assert kwargs["bam_file"] == seen_view[0]
        assert Path(str(kwargs["bam_file"])).read_bytes() == b"alignment-A"
        return "chr1:10-20"

    def finish_preflight(**kwargs: object) -> AlignmentPlan:
        binding = kwargs["binding"]
        assert isinstance(binding, AlignmentBinding)
        assert binding.is_open
        assert kwargs["in_path"] == str(source)
        assert kwargs["bound_view_path"] == seen_view[0]
        return AlignmentPlan(
            input_path=str(source),
            view_path=seen_view[0],
            file_format="bam",
            index_path=f"{seen_view[0]}.bai",
            reference_path=None,
            reference_source="not-required",
            uncovered_contigs=(),
            unmapped_scan="indexed",
            binding=binding,
        )

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", side_effect=read_bound_header),
        mock.patch(
            "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
            side_effect=resolve_bound_region,
        ),
        mock.patch("vntyper.scripts.pipeline_alignment.run_preflight", side_effect=finish_preflight),
    ):
        prepared = prepare_input_alignment_preflight(
            in_path=source,
            input_type="BAM",
            output_dir=output,
            config={},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions="chr1:10-20",
            reference_fasta=None,
            fast_mode=False,
            alignment_validator=validate_bound,
            validation_cwd="/project",
        )

    try:
        assert source.read_bytes() == b"alignment-B"
        assert Path(prepared.plan.view_path).read_bytes() == b"alignment-A"
        assert prepared.bed_file.read_bytes() == b"chr1\t10\t20\n"
    finally:
        prepared.plan.close()


def test_owned_cram_boundary_drops_an_out_of_tree_header_uri_and_keeps_the_cli_reference(
    tmp_path: Path,
) -> None:
    """An untrusted header path cannot block the explicit reference at the owned pipeline boundary."""
    output = tmp_path / "run-output"
    output.mkdir()
    source = tmp_path / "upload" / "sample.cram"
    source.parent.mkdir()
    source.write_bytes(b"owned CRAM")
    outside_reference = tmp_path / "private" / "header-reference.fa"
    explicit_reference = tmp_path / "operator-reference.fa"
    explicit_reference.write_text(">chr1\nAAAA\n", encoding="ascii")
    header = f"@SQ\tSN:chr1\tM5:digest\tLN:100\tUR:{outside_reference.as_uri()}\n"

    def finish_preflight(**kwargs: object) -> AlignmentPlan:
        binding = kwargs["binding"]
        assert isinstance(binding, AlignmentBinding)
        assert kwargs["reference_fasta"] == str(explicit_reference)
        assert kwargs["header_reference_paths"] == ()
        assert kwargs["header_references"] == ()
        return AlignmentPlan(
            input_path=str(source),
            view_path=str(kwargs["bound_view_path"]),
            file_format="cram",
            index_path=f"{kwargs['bound_view_path']}.crai",
            reference_path=str(explicit_reference),
            reference_source="cli",
            uncovered_contigs=(),
            unmapped_scan="indexed",
            binding=binding,
        )

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", return_value=header),
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_declared_assembly"),
        mock.patch(
            "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
            return_value="chr1:10-20",
        ),
        mock.patch("vntyper.scripts.pipeline_alignment.run_preflight", side_effect=finish_preflight),
    ):
        prepared = prepare_input_alignment_preflight(
            in_path=source,
            input_type="CRAM",
            output_dir=output,
            config={},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions="chr1:10-20",
            reference_fasta=explicit_reference,
            fast_mode=True,
        )

    try:
        assert prepared.plan.reference_source == "cli"
        assert prepared.plan.reference_path == str(explicit_reference)
    finally:
        prepared.plan.close()
        restore_reference_resolution(prepared.previous_ref_path)


def test_cram_quickcheck_failure_uses_the_alignment_validation_phase_and_releases_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A malformed alignment must not be reported as a reference-policy failure."""
    output = tmp_path / "run-output"
    output.mkdir()
    input_path = tmp_path / "patient-input" / "sample.cram"
    input_path.parent.mkdir()
    input_path.write_bytes(b"malformed CRAM")
    original_ref_path = "/operator/original/%s"
    monkeypatch.setenv("REF_PATH", original_ref_path)
    original = ValueError("quickcheck rejected malformed CRAM")

    with pytest.raises(ValueError) as raised:
        prepare_input_alignment_preflight(
            in_path=input_path,
            input_type="CRAM",
            output_dir=output,
            config={"cram": {"local_ref_path": "/local/cache/%s"}},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
            alignment_validator=mock.Mock(side_effect=original),
        )

    assert raised.value is original
    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8")) == {
        "code": "alignment_header_invalid",
        "message": "Alignment preflight rejected the alignment header or declared assembly; verify the input "
        "and --reference-assembly.",
        "candidates": [],
    }
    assert os.environ["REF_PATH"] == original_ref_path
    assert not os.path.lexists(output / "fastq_bam_processing" / "input.cram")


def test_input_alignment_header_failure_has_a_stable_artifact_and_restores_ref_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Assembly/header rejection belongs to a curated outer preflight phase."""
    output = tmp_path / "run-output"
    output.mkdir()
    input_path = tmp_path / "patient-input" / "sample.bam"
    input_path.parent.mkdir()
    input_path.touch()
    monkeypatch.setenv("REF_PATH", "operator-value")

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", return_value="@SQ\tSN:chr1\n"),
        mock.patch(
            "vntyper.scripts.pipeline_alignment.enforce_declared_assembly",
            side_effect=ValueError("private assembly detail"),
        ),
        mock.patch("vntyper.scripts.pipeline_alignment.run_preflight") as preflight,
        pytest.raises(ValueError, match="private assembly detail"),
    ):
        prepare_input_alignment_preflight(
            in_path=input_path,
            input_type="BAM",
            output_dir=output,
            config={"cram": {"local_ref_path": "/local/cache/%s"}},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    assert os.environ["REF_PATH"] == "operator-value"
    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8")) == {
        "code": "alignment_header_invalid",
        "message": "Alignment preflight rejected the alignment header or declared assembly; verify the input "
        "and --reference-assembly.",
        "candidates": [],
    }
    preflight.assert_not_called()


def test_default_cram_remote_header_uri_fails_at_the_owned_policy_boundary_before_work(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A header network target is rejected before assembly, target, probe, or stage work."""
    output = tmp_path / "run-output"
    output.mkdir()
    (output / "preflight_error.json").write_text('{"code":"stale"}\n', encoding="utf-8")
    monkeypatch.setenv("REF_PATH", "operator-value")
    remote_uri = "http://127.0.0.1:8765/private/reference.fa"
    header = f"@SQ\tSN:chr7\tLN:100\tUR:{remote_uri}\n"
    input_path = tmp_path / "patient-input" / "sample.cram"
    input_path.parent.mkdir()
    input_path.touch()

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", return_value=header),
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_declared_assembly") as assembly,
        mock.patch("vntyper.scripts.pipeline_alignment.prepare_alignment_target") as target,
        mock.patch("vntyper.scripts.pipeline_alignment.run_preflight") as preflight,
        mock.patch("vntyper.scripts.alignment_preflight.capture_command") as capture,
        pytest.raises(ValueError) as raised,
    ):
        prepare_input_alignment_preflight(
            in_path=input_path,
            input_type="CRAM",
            output_dir=output,
            config={"cram": {"local_ref_path": "/local/cache/%s"}},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    message = str(raised.value)
    assert "contig=chr7" in message
    assert "scheme=http" in message
    assert remote_uri not in message
    assert os.environ["REF_PATH"] == "operator-value"
    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8")) == {
        "code": "reference_policy_invalid",
        "message": message,
        "candidates": [],
    }
    assembly.assert_not_called()
    target.assert_not_called()
    preflight.assert_not_called()
    capture.assert_not_called()


@pytest.mark.parametrize("invalid", ["false", 0, None])
def test_owned_cram_boundary_rejects_non_boolean_ambient_policy_before_header_work(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    invalid: object,
) -> None:
    """The pipeline boundary cannot reinterpret JSON strings, numbers, or null as a waiver."""
    output = tmp_path / "run-output"
    output.mkdir()
    monkeypatch.setenv("REF_PATH", "operator-value")

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header") as header,
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_declared_assembly") as assembly,
        mock.patch("vntyper.scripts.pipeline_alignment.prepare_alignment_target") as target,
        mock.patch("vntyper.scripts.pipeline_alignment.run_preflight") as preflight,
        pytest.raises(ValueError, match="allow_ambient_reference_resolution must be true or false"),
    ):
        prepare_input_alignment_preflight(
            in_path=tmp_path / "patient-input" / "sample.cram",
            input_type="CRAM",
            output_dir=output,
            config={"cram": {"allow_ambient_reference_resolution": invalid}},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    assert os.environ["REF_PATH"] == "operator-value"
    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8"))["code"] == (
        "reference_policy_invalid"
    )
    header.assert_not_called()
    assembly.assert_not_called()
    target.assert_not_called()
    preflight.assert_not_called()


def test_outer_boundary_preserves_the_inner_structured_reference_payload(tmp_path: Path) -> None:
    """The outer preparation owner must not replace candidate-specific diagnostics."""
    output = tmp_path / "run-output"
    output.mkdir()
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t10\t20\n", encoding="utf-8")
    input_path = tmp_path / "patient-input" / "sample.cram"
    input_path.parent.mkdir()
    input_path.touch()
    payload = {
        "code": "reference_unresolved",
        "message": "curated reference failure",
        "candidates": [["cli", "reference.fa", "probe exited non-zero"]],
    }

    def fail_with_structured_payload(*args: object, **kwargs: object) -> None:
        del args
        context = kwargs["failure_context"]
        assert isinstance(context, PreflightErrorContext)
        context.payload = payload
        raise RuntimeError("private probe detail")

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", return_value="@SQ\tSN:chr1\n"),
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_declared_assembly"),
        mock.patch(
            "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
            return_value="chr1:10-20",
        ),
        mock.patch("vntyper.scripts.pipeline_alignment.run_preflight", side_effect=fail_with_structured_payload),
        pytest.raises(RuntimeError, match="private probe detail"),
    ):
        prepare_input_alignment_preflight(
            in_path=input_path,
            input_type="CRAM",
            output_dir=output,
            config={},
            threads=1,
            reference_assembly="hg19",
            bed_file=bed,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8")) == payload


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


def test_preflight_kwargs_preserve_m5_values_for_every_header_contig(tmp_path: Path) -> None:
    """Terminal stream diagnostics can attribute a later BED failure to its own checksum."""
    bed = tmp_path / "targets.bed"
    bed.write_text("chr1\t10\t20\nchr2\t30\t40\n", encoding="utf-8")

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
        alignment_header=("@SQ\tSN:chr1\tLN:100\tM5:first\n@SQ\tSN:chr2\tLN:200\tM5:second\n"),
    )

    assert kwargs["header_m5s"] == (("chr1", "first"), ("chr2", "second"))
