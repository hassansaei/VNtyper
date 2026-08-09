"""Pre-write ownership of optional pipeline archive destinations."""

from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import pipeline as pipeline_module
from vntyper.scripts.alignment_contract import index_candidate_names
from vntyper.scripts.alignment_target_io import BWA_INDEX_EXTENSIONS, reference_index_paths

pytestmark = pytest.mark.unit

OPERATOR_BYTES = b"operator-owned-bytes-must-not-be-replaced-by-an-archive"


def _write_operator_file(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(OPERATOR_BYTES)
    return path


def _reference_roles(prefix: str, reference: Path) -> dict[str, Path]:
    sidecars = reference_index_paths(reference)
    return {
        prefix: reference,
        f"{prefix}_fai": sidecars[0],
        f"{prefix}_gzi": sidecars[1],
        f"{prefix}_dict": sidecars[2],
    }


def _operator_case(tmp_path: Path, mode: str) -> tuple[dict[str, Path], dict[str, Any], dict[str, Any]]:
    inputs = tmp_path / "inputs"
    config = deepcopy(MINIMAL_CONFIG)
    if mode == "fastq":
        fastq1 = inputs / "reads_R1.fastq.gz"
        fastq2 = inputs / "reads_R2.fastq.gz"
        bed = inputs / "target.bed"
        reference = inputs / "bwa.fa"
        roles = {
            "fastq1": fastq1,
            "fastq2": fastq2,
            "bed": bed,
            "bwa_reference": reference,
            **{
                f"bwa_index_{extension.removeprefix('.')}": Path(f"{reference}{extension}")
                for extension in (*BWA_INDEX_EXTENSIONS, ".custom")
            },
        }
        config["tool_params"] = {"bwa_index_extensions": [".custom"]}
        kwargs = {
            "fastq1": str(fastq1),
            "fastq2": str(fastq2),
            "bed_file": str(bed),
            "bwa_reference": str(reference),
        }
        return roles, kwargs, config

    alignment = inputs / f"patient.{mode}"
    bed = inputs / "target.bed"
    roles = {"alignment": alignment, "bed": bed}
    for position, candidate in enumerate(index_candidate_names(str(alignment), mode)):
        roles[f"alignment_index_{position}"] = Path(candidate)
    kwargs = {mode: str(alignment), "bed_file": str(bed)}
    if mode == "cram":
        explicit = inputs / "explicit.fa"
        configured_cram = inputs / "configured-cram.fa"
        configured_bwa = inputs / "configured-bwa.fa"
        roles.update(_reference_roles("explicit_reference", explicit))
        roles.update(_reference_roles("configured_cram_reference", configured_cram))
        roles.update(_reference_roles("configured_bwa_reference", configured_bwa))
        config["reference_data"]["cram_reference_hg19"] = str(configured_cram)
        config["reference_data"]["bwa_reference_hg19"] = str(configured_bwa)
        kwargs["reference_fasta"] = str(explicit)
    return roles, kwargs, config


ARCHIVE_ROLE_CASES = (
    *(pytest.param("fastq", role, id=f"fastq-{role}") for role in _operator_case(Path("/tmp"), "fastq")[0]),
    *(pytest.param("bam", role, id=f"bam-{role}") for role in _operator_case(Path("/tmp"), "bam")[0]),
    *(pytest.param("cram", role, id=f"cram-{role}") for role in _operator_case(Path("/tmp"), "cram")[0]),
)


@pytest.mark.parametrize("mode, protected_role", ARCHIVE_ROLE_CASES)
@pytest.mark.parametrize("archive_format, suffix", [("zip", ".zip"), ("tar.gz", ".tar.gz")])
@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_archive_alias_of_every_operator_role_fails_before_pipeline_work(
    tmp_path: Path,
    mode: str,
    protected_role: str,
    archive_format: str,
    suffix: str,
    alias_kind: str,
) -> None:
    """The pre-existing snapshot check refuses resolved and same-inode aliases."""
    roles, input_kwargs, config = _operator_case(tmp_path, mode)
    for path in roles.values():
        _write_operator_file(path)
    protected = roles[protected_role]
    output = tmp_path / "results"
    archive = Path(f"{output}{suffix}")
    if alias_kind == "symlink":
        archive.symlink_to(protected)
    else:
        archive.hardlink_to(protected)
    protected_inode = protected.stat().st_ino
    archive_inode = archive.stat().st_ino

    with (
        mock.patch.object(Path, "mkdir", autospec=True) as mkdir,
        mock.patch.object(
            pipeline_module, "create_output_directories", wraps=pipeline_module.create_output_directories
        ) as create_dirs,
    ):
        harness = run_pipeline_under_harness(
            output,
            create_output_dir=False,
            config=config,
            archive_results=True,
            archive_format=archive_format,
            expect_failure=True,
            **input_kwargs,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert protected.stat().st_ino == protected_inode
    assert archive.stat().st_ino == archive_inode
    assert protected.read_bytes() == OPERATOR_BYTES
    mkdir.assert_not_called()
    create_dirs.assert_not_called()
    for stage in harness.stages.values():
        stage.assert_not_called()


@pytest.mark.parametrize("archive_format, suffix", [("zip", ".zip"), ("tar.gz", ".tar.gz")])
@pytest.mark.parametrize("mode", ["fastq", "bam"])
def test_lexically_exact_archive_input_collision_fails_before_work(
    tmp_path: Path,
    archive_format: str,
    suffix: str,
    mode: str,
) -> None:
    """A literal operator-input destination is refused without relying on inode evidence."""
    output = tmp_path / "results"
    archive = _write_operator_file(Path(f"{output}{suffix}"))
    reference = _write_operator_file(tmp_path / "reference.fa")
    input_kwargs: dict[str, Any] = (
        {"fastq1": str(archive), "bwa_reference": str(reference)} if mode == "fastq" else {"bam": str(archive)}
    )
    original_inode = archive.stat().st_ino

    with mock.patch.object(Path, "mkdir", autospec=True) as mkdir:
        harness = run_pipeline_under_harness(
            output,
            create_output_dir=False,
            archive_results=True,
            archive_format=archive_format,
            expect_failure=True,
            **input_kwargs,
        )

    assert isinstance(harness.error, SystemExit)
    assert archive.stat().st_ino == original_inode
    assert archive.read_bytes() == OPERATOR_BYTES
    mkdir.assert_not_called()
    for stage in harness.stages.values():
        stage.assert_not_called()


def test_trailing_output_separator_still_guards_sibling_archive_before_work(tmp_path: Path) -> None:
    """Archive ownership uses the normalized output name, never an in-directory dotfile."""
    output = tmp_path / "results"
    archive = _write_operator_file(tmp_path / "results.zip")
    reference = _write_operator_file(tmp_path / "reference.fa")
    original_inode = archive.stat().st_ino

    with mock.patch.object(Path, "mkdir", autospec=True) as mkdir:
        harness = run_pipeline_under_harness(
            output,
            pipeline_output_dir=f"{output}{os.sep}",
            create_output_dir=False,
            fastq1=str(archive),
            bwa_reference=str(reference),
            archive_results=True,
            archive_format="zip",
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert archive.stat().st_ino == original_inode
    assert archive.read_bytes() == OPERATOR_BYTES
    mkdir.assert_not_called()
    for stage in harness.stages.values():
        stage.assert_not_called()


def test_unsupported_archive_format_fails_before_directory_or_stage_work(tmp_path: Path) -> None:
    """Archive-format validation is part of the pre-write pipeline boundary."""
    fastq = _write_operator_file(tmp_path / "reads.fastq.gz")
    reference = _write_operator_file(tmp_path / "reference.fa")
    output = tmp_path / "results"

    with (
        mock.patch.object(Path, "mkdir", autospec=True) as mkdir,
        mock.patch.object(
            pipeline_module, "create_output_directories", wraps=pipeline_module.create_output_directories
        ) as create_dirs,
    ):
        harness = run_pipeline_under_harness(
            output,
            create_output_dir=False,
            fastq1=str(fastq),
            bwa_reference=str(reference),
            archive_results=True,
            archive_format="unsupported",
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not output.exists()
    mkdir.assert_not_called()
    create_dirs.assert_not_called()
    for stage in harness.stages.values():
        stage.assert_not_called()


@pytest.mark.parametrize("archive_format, suffix", [("zip", ".zip"), ("tar.gz", ".tar.gz")])
def test_archive_resolving_into_patient_tree_fails_before_work(
    tmp_path: Path,
    archive_format: str,
    suffix: str,
) -> None:
    """The archive may not resolve to unenumerated state beside a patient alignment."""
    patient_tree = tmp_path / "patient"
    alignment = _write_operator_file(patient_tree / "patient.bam")
    patient_note = _write_operator_file(patient_tree / "do-not-replace.txt")
    output = tmp_path / "results"
    archive = Path(f"{output}{suffix}")
    archive.symlink_to(patient_note)
    original_inode = patient_note.stat().st_ino

    with mock.patch.object(Path, "mkdir", autospec=True) as mkdir:
        harness = run_pipeline_under_harness(
            output,
            create_output_dir=False,
            bam=str(alignment),
            archive_results=True,
            archive_format=archive_format,
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert patient_note.stat().st_ino == original_inode
    assert patient_note.read_bytes() == OPERATOR_BYTES
    mkdir.assert_not_called()
    for stage in harness.stages.values():
        stage.assert_not_called()


@pytest.mark.parametrize("archive_format, suffix", [("zip", ".zip"), ("tar.gz", ".tar.gz")])
def test_fastq_inputs_and_archive_can_share_a_safe_sibling_directory(
    tmp_path: Path,
    archive_format: str,
    suffix: str,
) -> None:
    """Ordinary sibling layout remains supported when no paths alias."""
    fastq = _write_operator_file(tmp_path / "reads.fastq.gz")
    reference = _write_operator_file(tmp_path / "reference.fa")
    output = tmp_path / "results"

    harness = run_pipeline_under_harness(
        output,
        create_output_dir=False,
        fastq1=str(fastq),
        bwa_reference=str(reference),
        archive_results=True,
        archive_format=archive_format,
    )

    assert harness.error is None
    assert Path(f"{output}{suffix}").is_file()
    assert fastq.read_bytes() == OPERATOR_BYTES


def test_trailing_output_separator_writes_sibling_archive_not_dotfile(tmp_path: Path) -> None:
    """The writer and its ownership guard derive the same normalized destination."""
    fastq = _write_operator_file(tmp_path / "reads.fastq.gz")
    reference = _write_operator_file(tmp_path / "reference.fa")
    output = tmp_path / "results"

    harness = run_pipeline_under_harness(
        output,
        pipeline_output_dir=f"{output}{os.sep}",
        fastq1=str(fastq),
        bwa_reference=str(reference),
        archive_results=True,
        archive_format="zip",
    )

    assert harness.error is None
    assert (tmp_path / "results.zip").is_file()
    assert not (output / ".zip").exists()
