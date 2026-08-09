"""Alignment-mode pipeline setup must not overwrite operator-owned inputs."""

import json
from copy import deepcopy
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import pipeline as pipeline_module
from vntyper.scripts.alignment_target_io import protect_alignment_inputs
from vntyper.scripts.pipeline_alignment import prepare_input_alignment_preflight

pytestmark = pytest.mark.unit


def test_input_guard_leaves_invalid_reference_policy_to_persisted_preflight_owner(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = inputs / "input.cram"
    alignment.write_bytes(b"operator-owned")
    bed = inputs / "target.bed"
    bed.write_text("chr1\t1\t2\n", encoding="utf-8")
    config = deepcopy(MINIMAL_CONFIG)
    config["cram"] = {"reference_candidate_order": ["cli"]}

    protect_alignment_inputs(output, alignment, "cram", bed, None, config, "hg19")
    with (
        mock.patch("vntyper.scripts.pipeline_alignment.pin_reference_resolution", return_value=None),
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", return_value="@SQ\tSN:chr1\n"),
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_header_reference_policy"),
        mock.patch("vntyper.scripts.pipeline_alignment.enforce_declared_assembly"),
        mock.patch("vntyper.scripts.pipeline_alignment.prepare_alignment_target", return_value=bed),
        mock.patch(
            "vntyper.scripts.pipeline_alignment.get_region_string_with_fallback",
            return_value="chr1:1-2",
        ),
        pytest.raises(ValueError, match="terminal htslib_resolved"),
    ):
        prepare_input_alignment_preflight(
            in_path=alignment,
            input_type="CRAM",
            output_dir=output,
            config=config,
            threads=1,
            reference_assembly="hg19",
            bed_file=bed,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    assert json.loads((output / "preflight_error.json").read_text(encoding="utf-8")) == {
        "code": "reference_policy_invalid",
        "message": "CRAM reference candidate policy is invalid; configure a list ending in exactly one terminal "
        "htslib-resolved candidate.",
        "candidates": [],
    }


def test_exact_null_reference_override_does_not_protect_unused_family_fallback(tmp_path: Path) -> None:
    """The pre-write guard must protect the same candidates preflight can select."""
    output = tmp_path / "output"
    output.mkdir()
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = inputs / "input.cram"
    alignment.write_bytes(b"operator-owned")
    config = deepcopy(MINIMAL_CONFIG)
    config["reference_data"]["cram_reference_hg19"] = str(output / "unused-family.fa")
    config["reference_data"]["cram_reference_hg19_ncbi"] = None

    harness = run_pipeline_under_harness(
        output,
        config=config,
        cram=str(alignment),
        reference_assembly="hg19_ncbi",
    )

    assert harness.error is None
    harness.stages["validate_bam_file"].assert_called_once()


@pytest.mark.parametrize("input_type", ["BAM", "CRAM"])
@pytest.mark.parametrize("protected_role", ["bed", "reference_fasta"])
@pytest.mark.parametrize("destination_name", ["preflight_error.json", "pipeline_summary.json"])
def test_alignment_operator_input_at_late_root_artifact_fails_before_work(
    input_type: str,
    protected_role: str,
    destination_name: str,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    destination = output / destination_name
    destination.write_bytes(b"operator-owned")
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = inputs / f"input.{input_type.lower()}"
    bed = inputs / "target.bed"
    reference = inputs / "reference.fa"
    for source in (alignment, bed, reference):
        source.write_bytes(b"operator-owned")
    if protected_role == "bed":
        bed = destination
    else:
        reference = destination
    selected_input: dict[str, Any] = {input_type.lower(): str(alignment)}
    caplog.clear()

    with mock.patch.object(pipeline_module, "prepare_input_alignment_preflight", autospec=True) as prepare:
        harness = run_pipeline_under_harness(
            output,
            bed_file=str(bed),
            reference_fasta=str(reference),
            expect_failure=True,
            **selected_input,
        )

    assert isinstance(harness.error, SystemExit)
    assert "inside pipeline output root" in caplog.text
    assert destination.read_bytes() == b"operator-owned"
    harness.stages["validate_bam_file"].assert_not_called()
    prepare.assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()


@pytest.mark.parametrize("reference_key", ["cram_reference_hg19", "bwa_reference_hg19"])
@pytest.mark.parametrize("destination_name", ["preflight_error.json", "pipeline_summary.json"])
def test_configured_cram_reference_at_late_root_artifact_fails_before_work(
    reference_key: str,
    destination_name: str,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    destination = output / destination_name
    destination.write_bytes(b"operator-owned")
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = inputs / "input.cram"
    alignment.write_bytes(b"operator-owned")
    config = deepcopy(MINIMAL_CONFIG)
    config["reference_data"][reference_key] = str(destination)

    harness = run_pipeline_under_harness(
        output,
        config=config,
        cram=str(alignment),
        reference_fasta=None,
        expect_failure=True,
    )

    assert isinstance(harness.error, SystemExit)
    assert destination.read_bytes() == b"operator-owned"
    harness.stages["validate_bam_file"].assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()


@pytest.mark.parametrize("reference_key", ["cram_reference_hg19", "bwa_reference_hg19"])
@pytest.mark.parametrize("sidecar_suffix", [".fai", ".gzi", ".dict"])
def test_configured_cram_reference_sidecar_alias_fails_before_work(
    reference_key: str,
    sidecar_suffix: str,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    reference = inputs / "reference.fa"
    reference.write_bytes(b"reference")
    sidecar = (
        reference.with_suffix(sidecar_suffix) if sidecar_suffix == ".dict" else Path(f"{reference}{sidecar_suffix}")
    )
    sidecar.write_bytes(b"operator-owned-sidecar")
    (output / "pipeline_summary.json").hardlink_to(sidecar)
    alignment = inputs / "input.cram"
    alignment.write_bytes(b"operator-owned")
    config = deepcopy(MINIMAL_CONFIG)
    config["reference_data"][reference_key] = str(reference)

    harness = run_pipeline_under_harness(output, config=config, cram=str(alignment), expect_failure=True)

    assert isinstance(harness.error, SystemExit)
    assert sidecar.read_bytes() == b"operator-owned-sidecar"
    harness.stages["validate_bam_file"].assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()


@pytest.mark.parametrize("input_type", ["BAM", "CRAM"])
def test_alignment_mode_rejects_redirected_stage_directory_before_work(
    input_type: str,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    redirected = tmp_path / "redirected"
    redirected.mkdir()
    (output / "fastq_bam_processing").symlink_to(redirected, target_is_directory=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = inputs / f"input.{input_type.lower()}"
    bed = inputs / "target.bed"
    reference = inputs / "reference.fa"
    for source in (alignment, bed, reference):
        source.write_bytes(b"operator-owned")
    selected_input: dict[str, Any] = {input_type.lower(): str(alignment)}
    caplog.clear()

    with mock.patch.object(pipeline_module, "prepare_input_alignment_preflight", autospec=True) as prepare:
        harness = run_pipeline_under_harness(
            output,
            bed_file=str(bed),
            reference_fasta=str(reference),
            expect_failure=True,
            **selected_input,
        )

    assert isinstance(harness.error, SystemExit)
    assert "symlink directory" in caplog.text
    harness.stages["validate_bam_file"].assert_not_called()
    prepare.assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()


@pytest.mark.parametrize("input_type", ["BAM", "CRAM"])
def test_alignment_mode_safe_rerun_accepts_regular_artifacts_and_preflight_view(
    input_type: str,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    stage = output / "fastq_bam_processing"
    stage.mkdir(parents=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = inputs / f"input.{input_type.lower()}"
    bed = inputs / "target.bed"
    reference = inputs / "reference.fa"
    for source in (alignment, bed, reference):
        source.write_bytes(b"operator-owned")
    (stage / f"input.{input_type.lower()}").symlink_to(alignment)
    (stage / f"input.{input_type.lower()}.bai").write_bytes(b"stale-index")
    (output / "preflight_error.json").write_bytes(b"stale-artifact")
    (output / "pipeline_summary.json").write_bytes(b"stale-artifact")

    protect_alignment_inputs(output, alignment, input_type.lower(), bed, reference, MINIMAL_CONFIG, "hg19")


@pytest.mark.parametrize("protected_role", ["alignment", "bed", "reference_fasta"])
def test_alignment_mode_rejects_unknown_output_hardlink_to_operator_input(
    protected_role: str,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = inputs / "input.cram"
    bed = inputs / "target.bed"
    reference = inputs / "reference.fa"
    protected = {"alignment": alignment, "bed": bed, "reference_fasta": reference}
    for source in protected.values():
        source.write_bytes(b"operator-owned")
    (output / "unknown-late-artifact").hardlink_to(protected[protected_role])

    with pytest.raises(ValueError, match="output-tree entry aliases protected input"):
        protect_alignment_inputs(output, alignment, "cram", bed, reference, MINIMAL_CONFIG, "hg19")


@pytest.mark.parametrize("input_mode", ["bam", "fastq"])
@pytest.mark.parametrize("output_route", ["direct", "symlinked"])
@pytest.mark.parametrize("artifact_name", ["pipeline_summary.json", "coverage/coverage_summary.tsv"])
def test_unrelated_output_hardlink_is_rejected_before_mkdir_or_pipeline_work(
    tmp_path: Path,
    input_mode: str,
    output_route: str,
    artifact_name: str,
) -> None:
    """Every later output entry needs exclusive ownership, not only input non-aliasing."""
    actual_output = tmp_path / "actual-output"
    actual_output.mkdir()
    artifact = actual_output / artifact_name
    artifact.parent.mkdir(parents=True, exist_ok=True)
    external = tmp_path / "unrelated-external-artifact"
    external.write_bytes(b"unrelated-external-bytes")
    artifact.hardlink_to(external)
    protected_inode = external.stat().st_ino
    protected_bytes = external.read_bytes()
    output = actual_output
    if output_route == "symlinked":
        output = tmp_path / "output-link"
        output.symlink_to(actual_output, target_is_directory=True)

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    if input_mode == "bam":
        alignment = inputs / "input.bam"
        alignment.write_bytes(b"operator-bam")
        input_kwargs: dict[str, Any] = {"bam": str(alignment)}
    else:
        fastq = inputs / "reads.fastq.gz"
        reference = inputs / "reference.fa"
        fastq.write_bytes(b"operator-fastq")
        reference.write_bytes(b"operator-reference")
        input_kwargs = {"fastq1": str(fastq), "bwa_reference": str(reference)}

    with (
        mock.patch.object(Path, "mkdir", autospec=True) as mkdir,
        mock.patch.object(
            pipeline_module,
            "create_output_directories",
            wraps=pipeline_module.create_output_directories,
        ) as create_dirs,
    ):
        harness = run_pipeline_under_harness(
            output,
            create_output_dir=False,
            expect_failure=True,
            **input_kwargs,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert external.stat().st_ino == protected_inode
    assert artifact.stat().st_ino == protected_inode
    assert external.read_bytes() == protected_bytes
    mkdir.assert_not_called()
    create_dirs.assert_not_called()
    for stage in harness.stages.values():
        stage.assert_not_called()
