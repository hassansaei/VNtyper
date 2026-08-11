"""Pipeline logging must not open any operator-owned input as its log."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from vntyper import cli

pytestmark = pytest.mark.unit

CANONICAL_BWA_INDEXES = (".amb", ".ann", ".bwt", ".pac", ".sa")


@pytest.fixture(autouse=True)
def _restore_root_logger():
    """Restore handlers that the real CLI logging setup replaces."""
    root = logging.getLogger()
    saved_handlers = list(root.handlers)
    saved_level = root.level
    try:
        yield
    finally:
        for handler in root.handlers:
            if handler not in saved_handlers:
                handler.close()
        root.handlers[:] = saved_handlers
        root.setLevel(saved_level)


def _config(reference: Path) -> dict[str, Any]:
    """Build the complete configuration surface used before CLI dispatch."""
    return {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19"},
        "reference_data": {"bwa_reference_hg19": str(reference)},
        "tool_params": {"bwa_index_extensions": [".custom-index"]},
    }


def _pipeline_case(tmp_path: Path, role: str) -> tuple[list[str], dict[str, Any], Path]:
    """Return argv, config and the protected path for one input role."""
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    bam = inputs / "sample.bam"
    cram = inputs / "sample.cram"
    fastq1 = inputs / "reads_R1.fastq.gz"
    fastq2 = inputs / "reads_R2.fastq.gz"
    bed = inputs / "regions.bed"
    cram_reference = inputs / "cram-reference.fa"
    bwa_reference = inputs / "bwa-reference.fa"
    for path in (bam, cram, fastq1, fastq2, bed, cram_reference, bwa_reference):
        path.write_bytes(f"operator-owned:{path.name}".encode())
    bam_indexes = (Path(f"{bam}.bai"), bam.with_suffix(".bai"), Path(f"{bam}.csi"), bam.with_suffix(".csi"))
    cram_indexes = (
        Path(f"{cram}.crai"),
        cram.with_suffix(".crai"),
        Path(f"{cram}.csi"),
        cram.with_suffix(".csi"),
    )
    bwa_indexes = tuple(Path(f"{bwa_reference}{suffix}") for suffix in (*CANONICAL_BWA_INDEXES, ".custom-index"))
    cram_reference_indexes = (
        Path(f"{cram_reference}.fai"),
        Path(f"{cram_reference}.gzi"),
        cram_reference.with_suffix(".dict"),
    )
    for path in (*bam_indexes, *cram_indexes, *cram_reference_indexes, *bwa_indexes):
        path.write_bytes(f"operator-owned:{path.name}".encode())

    config = _config(bwa_reference)
    output = tmp_path / "output"
    cases: dict[str, tuple[list[str], Path]] = {
        "bam": (["pipeline", "-o", str(output), "--bam", str(bam)], bam),
        "bam_index": (["pipeline", "-o", str(output), "--bam", str(bam)], bam_indexes[0]),
        "cram": (["pipeline", "-o", str(output), "--cram", str(cram)], cram),
        "cram_index": (["pipeline", "-o", str(output), "--cram", str(cram)], cram_indexes[0]),
        "fastq1": (["pipeline", "-o", str(output), "--fastq1", str(fastq1)], fastq1),
        "fastq2": (
            ["pipeline", "-o", str(output), "--fastq1", str(fastq1), "--fastq2", str(fastq2)],
            fastq2,
        ),
        "bed": (["pipeline", "-o", str(output), "--bam", str(bam), "--bed-file", str(bed)], bed),
        "cram_reference": (
            ["pipeline", "-o", str(output), "--cram", str(cram), "--reference-fasta", str(cram_reference)],
            cram_reference,
        ),
        "cram_reference_fai": (
            ["pipeline", "-o", str(output), "--cram", str(cram), "--reference-fasta", str(cram_reference)],
            cram_reference_indexes[0],
        ),
        "cram_reference_gzi": (
            ["pipeline", "-o", str(output), "--cram", str(cram), "--reference-fasta", str(cram_reference)],
            cram_reference_indexes[1],
        ),
        "cram_reference_dict": (
            ["pipeline", "-o", str(output), "--cram", str(cram), "--reference-fasta", str(cram_reference)],
            cram_reference_indexes[2],
        ),
        "bwa_reference": (
            ["pipeline", "-o", str(output), "--fastq1", str(fastq1), "--reference-assembly", "hg19"],
            bwa_reference,
        ),
        "canonical_bwa_index": (
            ["pipeline", "-o", str(output), "--fastq1", str(fastq1), "--reference-assembly", "hg19"],
            bwa_indexes[0],
        ),
        "configured_bwa_index": (
            ["pipeline", "-o", str(output), "--fastq1", str(fastq1), "--reference-assembly", "hg19"],
            bwa_indexes[-1],
        ),
    }
    argv, protected = cases[role]
    return argv, config, protected


def _alias_log(tmp_path: Path, protected: Path, alias_kind: str) -> Path:
    """Create an explicit log path with the requested relationship to an input."""
    if alias_kind == "direct":
        return protected
    log_file = tmp_path / f"{alias_kind}.log"
    if alias_kind == "symlink":
        log_file.symlink_to(protected)
    else:
        os.link(protected, log_file)
    return log_file


def test_real_debug_logging_cannot_append_to_an_explicit_bam_log_alias(tmp_path: Path, monkeypatch) -> None:
    """Opening the logger itself must not append records to a BAM input."""
    argv, config, bam = _pipeline_case(tmp_path, "bam")
    original = bam.read_bytes()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: None)

    exit_code = None
    try:
        cli.main(["--log-level", "DEBUG", "--log-file", str(bam), *argv])
    except SystemExit as exc:
        exit_code = exc.code

    assert bam.read_bytes() == original
    assert exit_code == 1


@pytest.mark.parametrize(
    "role",
    [
        "bam",
        "bam_index",
        "cram",
        "cram_index",
        "fastq1",
        "fastq2",
        "bed",
        "cram_reference",
        "cram_reference_fai",
        "cram_reference_gzi",
        "cram_reference_dict",
        "bwa_reference",
        "canonical_bwa_index",
        "configured_bwa_index",
    ],
)
@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_explicit_pipeline_log_alias_is_rejected_before_logging_or_dispatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    role: str,
    alias_kind: str,
) -> None:
    """Every supported input role gets the same pre-open alias protection."""
    argv, config, protected = _pipeline_case(tmp_path, role)
    original = protected.read_bytes()
    log_file = _alias_log(tmp_path, protected, alias_kind)
    setup = mock.Mock()
    handler = mock.Mock()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", setup)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handler)

    with mock.patch.object(Path, "mkdir", autospec=True) as mkdir, pytest.raises(SystemExit) as raised:
        cli.main(["--log-file", str(log_file), *argv])

    assert raised.value.code == 1
    assert protected.read_bytes() == original
    mkdir.assert_not_called()
    setup.assert_not_called()
    handler.assert_not_called()


@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_default_pipeline_log_alias_is_rejected_before_logging(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    alias_kind: str,
) -> None:
    """The implicit ``output/pipeline.log`` receives identical protection."""
    output = tmp_path / "output"
    output.mkdir()
    fastq = tmp_path / "reads.fastq.gz"
    fastq.write_bytes(b"operator-fastq")
    default_log = output / "pipeline.log"
    reference = default_log if alias_kind == "direct" else tmp_path / "reference.fa"
    if alias_kind == "direct":
        reference.write_bytes(b"operator-reference")
    else:
        reference.write_bytes(b"operator-reference")
        if alias_kind == "symlink":
            default_log.symlink_to(reference)
        else:
            os.link(reference, default_log)
    original = reference.read_bytes()
    config = _config(reference)
    setup = mock.Mock()
    handler = mock.Mock()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", setup)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handler)

    with pytest.raises(SystemExit) as raised:
        cli.main(["pipeline", "-o", str(output), "--fastq1", str(fastq), "--reference-assembly", "hg19"])

    assert raised.value.code == 1
    assert reference.read_bytes() == original
    setup.assert_not_called()
    handler.assert_not_called()


def test_safe_pipeline_log_preserves_setup_then_dispatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A non-alias log retains the existing setup-before-handler lifecycle."""
    argv, config, _bam = _pipeline_case(tmp_path, "bam")
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(["--log-file", str(tmp_path / "safe.log"), *argv])

    assert events == ["setup", "handler"]


@pytest.mark.parametrize("protected_role", ["fastq", "bed", "reference"])
def test_log_under_file_input_parent_is_safe_when_alignment_patient_tree_is_separate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    protected_role: str,
) -> None:
    """File inputs forbid aliases, not unrelated sibling output directories."""
    case = tmp_path / "case"
    patient = tmp_path / "patient"
    case.mkdir()
    patient.mkdir()
    output = case / "results"
    fastq = case / "reads.fastq.gz"
    bed = case / "regions.bed"
    reference = case / "reference.fa"
    bam = patient / "sample.bam"
    cram = patient / "sample.cram"
    for path in (fastq, bed, reference, bam, cram):
        path.write_bytes(f"operator-owned:{path.name}".encode())
    config = _config(reference)
    if protected_role == "fastq":
        argv = ["pipeline", "-o", str(output), "--fastq1", str(fastq)]
    elif protected_role == "bed":
        argv = ["pipeline", "-o", str(output), "--bam", str(bam), "--bed-file", str(bed)]
    else:
        argv = [
            "pipeline",
            "-o",
            str(output),
            "--cram",
            str(cram),
            "--reference-fasta",
            str(reference),
        ]
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(argv)

    assert events == ["setup", "handler"]


def test_non_pipeline_logging_is_unchanged(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The new pipeline guard must not change another subcommand's logging."""
    events: list[str] = []
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "cohort", lambda *args, **kwargs: events.append("handler"))

    cli.main(["--log-file", str(tmp_path / "cohort.log"), "cohort", "-i", str(tmp_path), "-o", str(tmp_path)])

    assert events == ["setup", "handler"]


@pytest.mark.parametrize("reference_part", ["reference", "fai", "gzi", "dict"])
@pytest.mark.parametrize("config_source", ["config_cram_reference", "config_bwa_reference"])
def test_selected_config_cram_candidate_cannot_be_opened_as_the_log(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    config_source: str,
    reference_part: str,
) -> None:
    """CRAM's exact ordered candidate policy also owns each FASTA sidecar."""
    input_dir = tmp_path / "inputs"
    reference_dir = tmp_path / "references"
    input_dir.mkdir()
    reference_dir.mkdir()
    cram = input_dir / "sample.cram"
    cram.write_bytes(b"operator-cram")
    family_reference = reference_dir / f"{config_source}.fa"
    family_reference.write_bytes(b"operator-reference")
    protected_paths = {
        "reference": family_reference,
        "fai": Path(f"{family_reference}.fai"),
        "gzi": Path(f"{family_reference}.gzi"),
        "dict": family_reference.with_suffix(".dict"),
    }
    protected = protected_paths[reference_part]
    if protected != family_reference:
        protected.write_bytes(b"operator-reference-index")
        log_file = tmp_path / f"{reference_part}.log"
        os.link(protected, log_file)
    else:
        log_file = protected
    config = {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19_ncbi"},
        "reference_data": {
            "cram_reference_hg19": str(family_reference) if config_source == "config_cram_reference" else None,
            "bwa_reference_hg19": str(family_reference) if config_source == "config_bwa_reference" else None,
        },
    }
    setup = mock.Mock()
    handler = mock.Mock()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", setup)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handler)

    with pytest.raises(SystemExit) as raised:
        cli.main(
            [
                "--log-file",
                str(log_file),
                "pipeline",
                "--cram",
                str(cram),
                "--reference-assembly",
                "hg19_ncbi",
            ]
        )

    assert raised.value.code == 1
    setup.assert_not_called()
    handler.assert_not_called()


def test_null_label_override_does_not_reintroduce_the_family_cram_candidate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An exact null override keeps the family fallback disabled for that label."""
    input_dir = tmp_path / "inputs"
    reference_dir = tmp_path / "references"
    input_dir.mkdir()
    reference_dir.mkdir()
    cram = input_dir / "sample.cram"
    cram.write_bytes(b"operator-cram")
    family_reference = reference_dir / "family-reference.fa"
    family_reference.write_bytes(b"unused-family-reference")
    config = {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19_ncbi"},
        "reference_data": {
            "cram_reference_hg19": str(family_reference),
            "cram_reference_hg19_ncbi": None,
            "bwa_reference_hg19": None,
        },
    }
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(
        [
            "--log-file",
            str(family_reference),
            "pipeline",
            "--cram",
            str(cram),
            "--reference-assembly",
            "hg19_ncbi",
        ]
    )

    assert events == ["setup", "handler"]


def test_configured_reference_is_protected_without_validating_probe_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Input ownership is independent of the later CRAM probe policy."""
    input_dir = tmp_path / "inputs"
    reference_dir = tmp_path / "references"
    input_dir.mkdir()
    reference_dir.mkdir()
    cram = input_dir / "sample.cram"
    cram.write_bytes(b"operator-cram")
    omitted_reference = reference_dir / "omitted-reference.fa"
    omitted_reference.write_bytes(b"not-selected")
    config = {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19"},
        "cram": {"reference_candidate_order": ["cli", "htslib_resolved"]},
        "reference_data": {"cram_reference_hg19": str(omitted_reference)},
    }
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    with pytest.raises(SystemExit) as raised:
        cli.main(["--log-file", str(omitted_reference), "pipeline", "--cram", str(cram)])

    assert raised.value.code == 1
    assert events == []


def test_invalid_cram_candidate_policy_is_left_to_persisted_preflight(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Log safety must not steal malformed-policy ownership from preflight."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    cram = input_dir / "sample.cram"
    cram.write_bytes(b"operator-cram")
    config = {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19"},
        "cram": {"reference_candidate_order": ["cli"]},
        "reference_data": {"cram_reference_hg19": str(tmp_path / "references" / "reference.fa")},
    }
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(
        [
            "--log-file",
            str(tmp_path / "logs" / "pipeline.log"),
            "pipeline",
            "--cram",
            str(cram),
        ]
    )

    assert events == ["setup", "handler"]


@pytest.mark.parametrize("malformed_value", [7, {"path": "reference.fa"}])
def test_malformed_configured_reference_value_is_left_to_preflight(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    malformed_value: object,
) -> None:
    """The pre-open guard only interprets strings as configured paths."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    cram = input_dir / "sample.cram"
    cram.write_bytes(b"operator-cram")
    config = {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19"},
        "reference_data": {"cram_reference_hg19": malformed_value},
    }
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(
        [
            "--log-file",
            str(tmp_path / "logs" / "pipeline.log"),
            "pipeline",
            "--cram",
            str(cram),
        ]
    )

    assert events == ["setup", "handler"]


@pytest.mark.parametrize("malformed_value", [7, {"path": "reference.fa"}])
def test_malformed_fastq_bwa_reference_is_left_to_pipeline_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    malformed_value: object,
) -> None:
    """The logging guard must not interpret a malformed BWA value as a path."""
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    fastq = inputs / "reads.fastq.gz"
    fastq.write_bytes(b"operator-fastq")
    config = {
        "cli_defaults": {"log_level": "INFO", "log_file": None},
        "default_values": {"reference_assembly": "hg19"},
        "reference_data": {"bwa_reference_hg19": malformed_value},
    }
    events: list[str] = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: events.append("setup"))
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda *args, **kwargs: events.append("handler"))

    cli.main(
        [
            "--log-file",
            str(tmp_path / "logs" / "pipeline.log"),
            "pipeline",
            "--fastq1",
            str(fastq),
        ]
    )

    assert events == ["setup", "handler"]


@pytest.mark.parametrize("log_route", ["direct-output", "symlinked-output", "symlinked-parent"])
def test_log_creation_inside_an_input_tree_is_rejected_before_parent_setup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log_route: str,
) -> None:
    """A new log cannot mutate a patient directory through lexical or resolved parents."""
    patient = tmp_path / "patient"
    patient.mkdir()
    bam = patient / "sample.bam"
    bam.write_bytes(b"operator-bam")
    config = _config(tmp_path / "reference.fa")
    setup = mock.Mock()
    handler = mock.Mock()
    monkeypatch.setattr(cli, "load_config", lambda _path=None: config)
    monkeypatch.setattr(cli, "setup_logging", setup)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", handler)

    if log_route == "direct-output":
        output = patient / "run"
        argv = ["pipeline", "-o", str(output), "--bam", str(bam)]
    elif log_route == "symlinked-output":
        actual_output = patient / "run"
        actual_output.mkdir()
        output = tmp_path / "output-link"
        output.symlink_to(actual_output, target_is_directory=True)
        argv = ["pipeline", "-o", str(output), "--bam", str(bam)]
    else:
        log_parent = tmp_path / "log-parent"
        log_parent.symlink_to(patient, target_is_directory=True)
        argv = ["--log-file", str(log_parent / "new.log"), "pipeline", "-o", str(tmp_path / "out"), "--bam", str(bam)]

    with mock.patch.object(Path, "mkdir", autospec=True) as mkdir, pytest.raises(SystemExit) as raised:
        cli.main(argv)

    assert raised.value.code == 1
    mkdir.assert_not_called()
    setup.assert_not_called()
    handler.assert_not_called()
    assert not (patient / "pipeline.log").exists()
    assert not (patient / "new.log").exists()


class TestTheGuardUsesTheSameReferenceTheRunWill:
    """A guard that inspects a different file than the run uses is not a guard."""

    @pytest.mark.parametrize(
        "label,key",
        [
            ("hg38_ensembl", "bwa_reference_hg38_ensembl"),
            ("hg38_ncbi", "bwa_reference_GRCh38"),
            ("GRCh37", "bwa_reference_GRCh37"),
            ("hg19_ensembl", "bwa_reference_hg19_ensembl"),
        ],
    )
    def test_the_exact_reference_cannot_be_used_as_a_log_file(self, tmp_path, label, key):
        from vntyper.scripts.cli_logging_safety import validate_pipeline_log_destination

        reference = tmp_path / "exact.fa"
        reference.write_text(">chr1\nACGT\n")
        config = {"reference_data": {key: str(reference), "bwa_reference_hg38": str(tmp_path / "ucsc.fa")}}
        args = argparse.Namespace(
            fastq1="r1.fq",
            fastq2="r2.fq",
            bam=None,
            cram=None,
            bed_file=None,
            reference_fasta=None,
            reference_assembly=label,
            log_file=str(reference),
        )
        with pytest.raises(ValueError, match="exact.fa"):
            validate_pipeline_log_destination(str(reference), args, config)

    @pytest.mark.parametrize("suffix", [".amb", ".ann", ".bwt", ".pac", ".sa"])
    def test_a_bwa_sidecar_of_the_exact_reference_cannot_be_used_either(self, tmp_path, suffix):
        from vntyper.scripts.cli_logging_safety import validate_pipeline_log_destination

        reference = tmp_path / "exact.fa"
        reference.write_text(">chr1\nACGT\n")
        sidecar = tmp_path / f"exact.fa{suffix}"
        sidecar.write_bytes(b"\x00")
        config = {"reference_data": {"bwa_reference_hg38_ensembl": str(reference)}}
        args = argparse.Namespace(
            fastq1="r1.fq",
            fastq2="r2.fq",
            bam=None,
            cram=None,
            bed_file=None,
            reference_fasta=None,
            reference_assembly="hg38_ensembl",
            log_file=str(sidecar),
        )
        with pytest.raises(ValueError):
            validate_pipeline_log_destination(str(sidecar), args, config)
