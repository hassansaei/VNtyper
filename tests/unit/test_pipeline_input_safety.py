"""Direct-FASTQ pipeline destinations must not overwrite operator inputs."""

from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, PipelineHarness, run_pipeline_under_harness
from vntyper.scripts import pipeline as pipeline_module
from vntyper.scripts.alignment_target_io import BWA_INDEX_EXTENSIONS, validate_fastq_pipeline_destinations

pytestmark = pytest.mark.unit


DESTRUCTIVE_DESTINATIONS = (
    "fastq_bam_processing/output_R1.fastq.gz",
    "fastq_bam_processing/output_R2.fastq.gz",
    "fastq_bam_processing/output.html",
    "fastq_bam_processing/output.json",
    "fastq_bam_processing/output_fastp.log",
    "alignment_processing/output_sorted.bam",
    "alignment_processing/output_sorted.bam.bai",
    "alignment_processing/output_alignment.log",
    "alignment_processing/output_index.log",
    "fastq_bam_processing/post_alignment.bam",
    "fastq_bam_processing/post_alignment.bam.bai",
    "fastq_bam_processing/post_alignment.bai",
    "fastq_bam_processing/post_alignment.bam.csi",
    "fastq_bam_processing/post_alignment.csi",
    "fastq_bam_processing/post_alignment_index.log",
    "fastq_bam_processing/post_alignment_idxstats.log",
    "fastq_bam_processing/post_alignment_alignment_probe.log",
    "fastq_bam_processing/output_sliced.bam",
    "fastq_bam_processing/output_unmapped.bam",
    "fastq_bam_processing/output_sliced_unmapped.bam",
    "fastq_bam_processing/output_other.fastq.gz",
    "fastq_bam_processing/output_single.fastq.gz",
    "fastq_bam_processing/output_slice.log",
    "fastq_bam_processing/output_filter.log",
    "fastq_bam_processing/output_merge.log",
    "fastq_bam_processing/output_index.log",
    "fastq_bam_processing/output_sort_fastq.log",
    "fastq_bam_processing/output_sliced.bam.bai",
    "fastq_bam_processing/output_sliced.bai",
    "fastq_bam_processing/output_sliced.bam.csi",
    "fastq_bam_processing/output_sliced.csi",
)

PROTECTED_INPUT_ROLES = (
    "fastq1",
    "fastq2",
    "bed",
    "bwa_reference",
    *(f"bwa_index{extension}" for extension in BWA_INDEX_EXTENSIONS),
)


def _run_alias_case(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    *,
    protected_role: str,
    destination_name: str,
    alias_kind: str = "direct",
) -> tuple[Path, PipelineHarness, mock.MagicMock]:
    output = tmp_path / "output"
    destination = output / destination_name
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(b"operator-owned")
    input_root = tmp_path / "inputs"
    input_root.mkdir()
    inputs = {
        "fastq1": input_root / "r1.fastq.gz",
        "fastq2": input_root / "r2.fastq.gz",
        "bed": input_root / "target.bed",
        "bwa_reference": input_root / "reference.fa",
    }
    for path in inputs.values():
        path.write_bytes(b"operator-owned")
    bwa_indexes = {extension: Path(f"{inputs['bwa_reference']}{extension}") for extension in BWA_INDEX_EXTENSIONS}
    for bwa_index in bwa_indexes.values():
        bwa_index.write_bytes(b"operator-owned")

    index_extension = protected_role.removeprefix("bwa_index") if protected_role.startswith("bwa_index") else None
    protected_path = bwa_indexes[index_extension] if index_extension is not None else inputs[protected_role]
    if alias_kind == "direct":
        if index_extension is not None:
            destination.unlink()
            destination.hardlink_to(protected_path)
        else:
            protected_path.unlink()
            inputs[protected_role] = destination
    else:
        protected_path.unlink()
        if alias_kind == "hardlink":
            protected_path.hardlink_to(destination)
        else:
            protected_path.symlink_to(destination)

    config = deepcopy(MINIMAL_CONFIG)
    config["tool_params"] = {"bwa_index_extensions": list(BWA_INDEX_EXTENSIONS)}
    caplog.clear()
    with mock.patch.object(pipeline_module, "create_output_directories", autospec=True) as create_dirs:
        harness = run_pipeline_under_harness(
            output,
            config=config,
            fastq1=str(inputs["fastq1"]),
            fastq2=str(inputs["fastq2"]),
            bed_file=str(inputs["bed"]),
            bwa_reference=str(inputs["bwa_reference"]),
            expect_failure=True,
        )
    return destination, harness, create_dirs


@pytest.mark.parametrize("protected_role", PROTECTED_INPUT_ROLES)
@pytest.mark.parametrize("destination_name", DESTRUCTIVE_DESTINATIONS)
def test_every_operator_input_alias_of_every_direct_fastq_destination_fails_before_work(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    protected_role: str,
    destination_name: str,
) -> None:
    destination, harness, create_dirs = _run_alias_case(
        tmp_path,
        caplog,
        protected_role=protected_role,
        destination_name=destination_name,
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert "protected input" in caplog.text or "inside pipeline output root" in caplog.text
    assert destination.read_bytes() == b"operator-owned"
    create_dirs.assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()
    harness.stages["process_fastq"].assert_not_called()
    harness.stages["align_and_sort_fastq"].assert_not_called()
    harness.stages["run_preflight"].assert_not_called()
    harness.stages["process_bam_to_fastq"].assert_not_called()


@pytest.mark.parametrize("protected_role", PROTECTED_INPUT_ROLES)
@pytest.mark.parametrize("alias_kind", ["hardlink", "symlink"])
def test_filesystem_aliases_across_stages_fail_before_work(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    protected_role: str,
    alias_kind: str,
) -> None:
    destination, harness, create_dirs = _run_alias_case(
        tmp_path,
        caplog,
        protected_role=protected_role,
        destination_name="fastq_bam_processing/output_other.fastq.gz",
        alias_kind=alias_kind,
    )

    assert isinstance(harness.error, SystemExit)
    assert destination.read_bytes() == b"operator-owned"
    create_dirs.assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()


def test_existing_single_link_regular_pipeline_artifacts_remain_rerunnable(tmp_path: Path) -> None:
    output = tmp_path / "output"
    for destination_name in DESTRUCTIVE_DESTINATIONS:
        if destination_name == "fastq_bam_processing/post_alignment.bam":
            continue
        destination = output / destination_name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(b"stale-run-artifact")
    (output / "fastq_bam_processing/post_alignment.bam").symlink_to(output / "alignment_processing/output_sorted.bam")
    input_root = tmp_path / "inputs"
    input_root.mkdir()
    fastq1 = input_root / "r1.fastq.gz"
    fastq2 = input_root / "r2.fastq.gz"
    bed = input_root / "target.bed"
    reference = input_root / "reference.fa"
    for source in (fastq1, fastq2, bed, reference, *(Path(f"{reference}{ext}") for ext in BWA_INDEX_EXTENSIONS)):
        source.write_bytes(b"operator-owned")
    config = deepcopy(MINIMAL_CONFIG)
    config["tool_params"] = {"bwa_index_extensions": list(BWA_INDEX_EXTENSIONS)}

    harness = run_pipeline_under_harness(
        output,
        config=config,
        fastq1=str(fastq1),
        fastq2=str(fastq2),
        bed_file=str(bed),
        bwa_reference=str(reference),
    )

    harness.call("process_fastq")


def test_generated_exact_bed_cannot_alias_another_pipeline_destination(tmp_path: Path) -> None:
    output = tmp_path / "output"
    generated_bed = output / "predefined_regions_hg19.bed"
    generated_bed.parent.mkdir()
    generated_bed.write_bytes(b"generated-target")
    destination = output / "fastq_bam_processing/output_single.fastq.gz"
    destination.parent.mkdir()
    destination.hardlink_to(generated_bed)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fastq1 = inputs / "r1.fastq.gz"
    reference = inputs / "reference.fa"
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")

    with pytest.raises(ValueError, match="multiple hard links"):
        validate_fastq_pipeline_destinations(output, fastq1, None, None, reference, MINIMAL_CONFIG)


def test_generated_bed_inside_output_is_pipeline_owned_and_safe(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    generated_bed = output / "predefined_regions_hg19.bed"
    generated_bed.write_bytes(b"generated-target")
    fastq1 = tmp_path / "r1.fastq.gz"
    reference = tmp_path / "reference.fa"
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")

    validate_fastq_pipeline_destinations(output, fastq1, None, None, reference, MINIMAL_CONFIG)


@pytest.mark.parametrize("protected_role", PROTECTED_INPUT_ROLES)
def test_every_operator_input_inside_output_root_is_rejected(protected_role: str, tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fastq1 = inputs / "r1.fastq.gz"
    fastq2 = inputs / "r2.fastq.gz"
    bed = inputs / "target.bed"
    reference = inputs / "reference.fa"
    for source in (fastq1, fastq2, bed, reference):
        source.write_bytes(b"operator-owned")
    indexes = {extension: Path(f"{reference}{extension}") for extension in BWA_INDEX_EXTENSIONS}
    for index in indexes.values():
        index.write_bytes(b"operator-owned")

    if protected_role.startswith("bwa_index"):
        extension = protected_role.removeprefix("bwa_index")
        indexes[extension].unlink()
        indexes[extension].symlink_to(output / "protected-index")
        (output / "protected-index").write_bytes(b"operator-owned")
    else:
        contained = output / protected_role
        contained.write_bytes(b"operator-owned")
        if protected_role == "fastq1":
            fastq1 = contained
        elif protected_role == "fastq2":
            fastq2 = contained
        elif protected_role == "bed":
            bed = contained
        else:
            reference = contained

    config = deepcopy(MINIMAL_CONFIG)
    config["tool_params"] = {"bwa_index_extensions": list(BWA_INDEX_EXTENSIONS)}
    with pytest.raises(ValueError, match="inside pipeline output root"):
        validate_fastq_pipeline_destinations(output, fastq1, fastq2, bed, reference, config)


@pytest.mark.parametrize("protected_role", ["fastq1", "fastq2", "bed", "bwa_reference"])
def test_resolved_operator_input_inside_symlinked_output_root_is_rejected(
    protected_role: str,
    tmp_path: Path,
) -> None:
    actual_output = tmp_path / "actual-output"
    actual_output.mkdir()
    output = tmp_path / "output-link"
    output.symlink_to(actual_output, target_is_directory=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    values = {
        "fastq1": inputs / "r1.fastq.gz",
        "fastq2": inputs / "r2.fastq.gz",
        "bed": inputs / "target.bed",
        "bwa_reference": inputs / "reference.fa",
    }
    for source in values.values():
        source.write_bytes(b"operator-owned")
    contained = actual_output / protected_role
    contained.write_bytes(b"operator-owned")
    values[protected_role] = contained

    with pytest.raises(ValueError, match="inside pipeline output root"):
        validate_fastq_pipeline_destinations(
            output,
            values["fastq1"],
            values["fastq2"],
            values["bed"],
            values["bwa_reference"],
            MINIMAL_CONFIG,
        )


@pytest.mark.parametrize("protected_role", PROTECTED_INPUT_ROLES)
def test_unknown_output_leaf_alias_of_every_external_operator_input_is_rejected(
    protected_role: str,
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fastq1 = inputs / "r1.fastq.gz"
    fastq2 = inputs / "r2.fastq.gz"
    bed = inputs / "target.bed"
    reference = inputs / "reference.fa"
    for source in (fastq1, fastq2, bed, reference):
        source.write_bytes(b"operator-owned")
    indexes = {extension: Path(f"{reference}{extension}") for extension in BWA_INDEX_EXTENSIONS}
    for index in indexes.values():
        index.write_bytes(b"operator-owned")
    protected = {
        "fastq1": fastq1,
        "fastq2": fastq2,
        "bed": bed,
        "bwa_reference": reference,
        **{f"bwa_index{extension}": index for extension, index in indexes.items()},
    }[protected_role]
    (output / "unknown-late-artifact").hardlink_to(protected)
    config = deepcopy(MINIMAL_CONFIG)
    config["tool_params"] = {"bwa_index_extensions": list(BWA_INDEX_EXTENSIONS)}

    with pytest.raises(ValueError, match="output-tree entry aliases protected input"):
        validate_fastq_pipeline_destinations(output, fastq1, fastq2, bed, reference, config)


def test_symlinked_output_subdirectory_is_rejected_before_work(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    redirected = tmp_path / "redirected"
    redirected.mkdir()
    (output / "fastq_bam_processing").symlink_to(redirected, target_is_directory=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fastq1 = inputs / "r1.fastq.gz"
    bed = inputs / "target.bed"
    reference = inputs / "reference.fa"
    for source in (fastq1, bed, reference):
        source.write_bytes(b"operator-owned")

    with pytest.raises(ValueError, match="symlink directory"):
        validate_fastq_pipeline_destinations(output, fastq1, None, bed, reference, MINIMAL_CONFIG)


def test_unknown_output_file_symlink_is_rejected_before_work(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    external = tmp_path / "external"
    external.write_bytes(b"external")
    (output / "unknown-late-artifact").symlink_to(external)
    fastq1 = tmp_path / "r1.fastq.gz"
    bed = tmp_path / "target.bed"
    reference = tmp_path / "reference.fa"
    for source in (fastq1, bed, reference):
        source.write_bytes(b"operator-owned")

    with pytest.raises(ValueError, match="unsafe file symlink"):
        validate_fastq_pipeline_destinations(output, fastq1, None, bed, reference, MINIMAL_CONFIG)


def test_non_directory_output_root_is_rejected(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.write_bytes(b"not-a-directory")
    fastq1 = tmp_path / "r1.fastq.gz"
    reference = tmp_path / "reference.fa"
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")

    with pytest.raises(ValueError, match="must be a directory"):
        validate_fastq_pipeline_destinations(output, fastq1, None, None, reference, MINIMAL_CONFIG)


def test_isolated_symlinked_output_root_remains_safe(tmp_path: Path) -> None:
    actual_output = tmp_path / "actual-output"
    actual_output.mkdir()
    output = tmp_path / "output-link"
    output.symlink_to(actual_output, target_is_directory=True)
    fastq1 = tmp_path / "r1.fastq.gz"
    reference = tmp_path / "reference.fa"
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")

    validate_fastq_pipeline_destinations(output, fastq1, None, None, reference, MINIMAL_CONFIG)


def test_symlinked_output_root_is_scanned_for_protected_hardlinks(tmp_path: Path) -> None:
    actual_output = tmp_path / "actual-output"
    actual_output.mkdir()
    output = tmp_path / "output-link"
    output.symlink_to(actual_output, target_is_directory=True)
    fastq1 = tmp_path / "r1.fastq.gz"
    reference = tmp_path / "reference.fa"
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")
    (actual_output / "unknown-late-artifact").hardlink_to(reference)

    with pytest.raises(ValueError, match="output-tree entry aliases protected input"):
        validate_fastq_pipeline_destinations(output, fastq1, None, None, reference, MINIMAL_CONFIG)


@pytest.mark.parametrize("configured_extensions", [[], [".custom-index"]])
def test_canonical_bwa_indexes_cannot_be_removed_from_pipeline_protection(
    configured_extensions: list[str],
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    fastq1 = tmp_path / "r1.fastq.gz"
    reference = tmp_path / "reference.fa"
    canonical_index = Path(f"{reference}.amb")
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")
    canonical_index.write_bytes(b"canonical-index")
    (output / "unknown-late-artifact").hardlink_to(canonical_index)
    config = {"tool_params": {"bwa_index_extensions": configured_extensions}}

    with pytest.raises(ValueError, match="output-tree entry aliases protected input"):
        validate_fastq_pipeline_destinations(output, fastq1, None, None, reference, config)


def test_non_regular_output_tree_entry_is_rejected(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    os.mkfifo(output / "named-pipe")
    fastq1 = tmp_path / "r1.fastq.gz"
    reference = tmp_path / "reference.fa"
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")

    with pytest.raises(ValueError, match="non-regular entry"):
        validate_fastq_pipeline_destinations(output, fastq1, None, None, reference, MINIMAL_CONFIG)


def test_absent_output_root_is_safe_for_external_operator_inputs(tmp_path: Path) -> None:
    fastq1 = tmp_path / "r1.fastq.gz"
    reference = tmp_path / "reference.fa"
    fastq1.write_bytes(b"fastq")
    reference.write_bytes(b"reference")

    validate_fastq_pipeline_destinations(tmp_path / "new-output", fastq1, None, None, reference, MINIMAL_CONFIG)


@pytest.mark.parametrize(
    ("protected_role", "custom_regions", "generated_name"),
    [
        ("fastq1", None, "predefined_regions_hg19.bed"),
        ("bwa_reference", None, "predefined_regions_hg19.bed"),
        ("fastq1", "chr1:1-2", "custom_regions.bed"),
        ("bwa_reference", "chr1:1-2", "custom_regions.bed"),
    ],
)
def test_generated_bed_cannot_clobber_operator_input_before_global_guard(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    protected_role: str,
    custom_regions: str | None,
    generated_name: str,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    generated_target = output / generated_name
    generated_target.write_bytes(b"operator-owned")
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fastq1 = generated_target if protected_role == "fastq1" else inputs / "r1.fastq.gz"
    reference = generated_target if protected_role == "bwa_reference" else inputs / "reference.fa"
    if fastq1 != generated_target:
        fastq1.write_bytes(b"fastq")
    if reference != generated_target:
        reference.write_bytes(b"reference")
    caplog.clear()

    with mock.patch.object(pipeline_module, "prepare_alignment_target", autospec=True) as prepare_target:
        prepare_target.return_value = generated_target
        harness = run_pipeline_under_harness(
            output,
            fastq1=str(fastq1),
            bwa_reference=str(reference),
            custom_regions=custom_regions,
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert "inside pipeline output root" in caplog.text
    prepare_target.assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()
    assert generated_target.read_bytes() == b"operator-owned"


def test_bwa_reference_at_late_summary_destination_is_rejected_before_work(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    destination, harness, create_dirs = _run_alias_case(
        tmp_path,
        caplog,
        protected_role="bwa_reference",
        destination_name="pipeline_summary.json",
    )

    assert isinstance(harness.error, SystemExit)
    assert "inside pipeline output root" in caplog.text
    create_dirs.assert_not_called()
    harness.stages["get_tool_versions"].assert_not_called()
    assert destination.read_bytes() == b"operator-owned"
