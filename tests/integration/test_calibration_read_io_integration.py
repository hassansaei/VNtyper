"""Generated alignment and FASTQ inputs exercise the calibration read adapter."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pysam
import pytest

from vntyper.scripts.calibration_intake_contract import InputArtifact
from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

pytestmark = pytest.mark.integration


def artifact(path: Path, artifact_format: str, *, mate_path: Path | None = None) -> InputArtifact:
    return InputArtifact(
        key=f"synthetic-{artifact_format.lower()}",
        specimen_key="synthetic-specimen",
        path=str(path),
        format=artifact_format,  # type: ignore[arg-type]
        mate_path=str(mate_path) if mate_path is not None else None,
        assembly="synthetic-assembly",
        assay_class="synthetic-assay",
        input_scope="full",
        preprocessing_id="synthetic-preprocessing-v1",
        replicate_group="synthetic-replicate",
        expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )


def write_fastq_pair(first: Path, second: Path) -> None:
    first.write_text(
        "@synthetic-pair/1 1:N:0:SYNTHETIC\nACGA\n+synthetic-pair/1 1:N:0:SYNTHETIC\nIJKL\n",
        encoding="ascii",
    )
    second.write_text(
        "@synthetic-pair/2 2:N:0:SYNTHETIC\nTTGC\n+synthetic-pair/2\nLMNO\n",
        encoding="ascii",
    )


def segment(header: pysam.AlignmentHeader, *, mate: int, flags: int | None = None) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment(header)
    record.query_name = "synthetic-pair"
    record.query_sequence = "ACGA" if mate == 1 else "GCAA"
    record.query_qualities = pysam.qualitystring_to_array("IJKL" if mate == 1 else "ONML")
    record.flag = flags if flags is not None else (97 if mate == 1 else 145)
    record.reference_id = 0
    record.reference_start = 10 if mate == 1 else 30
    record.mapping_quality = 37 if mate == 1 else 29
    record.cigar = ((0, 4),)
    record.next_reference_id = 0
    record.next_reference_start = 30 if mate == 1 else 10
    record.template_length = 24 if mate == 1 else -24
    return record


def generated_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    reference = tmp_path / "reference with spaces.fa"
    reference.write_text(">synthetic-contig\n" + "ACGT" * 100 + "\n", encoding="ascii")
    pysam.faidx(str(reference))
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "synthetic-contig", "LN": 400}]}
    )
    bam = tmp_path / "reads with spaces.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
        output.write(segment(header, mate=1))
        output.write(segment(header, mate=2))
        output.write(segment(header, mate=1, flags=353))
        output.write(segment(header, mate=2, flags=2193))
    cram = tmp_path / "reads with spaces.cram"
    with (
        pysam.AlignmentFile(str(bam), "rb") as source,
        pysam.AlignmentFile(str(cram), "wc", header=source.header, reference_filename=str(reference)) as output,
    ):
        for record in source.fetch(until_eof=True):
            output.write(record)
    first = tmp_path / "first reads.fastq"
    second = tmp_path / "second reads.fastq"
    write_fastq_pair(first, second)
    return reference, bam, cram, first, second


def test_generated_bam_and_fastq_have_the_same_name_neutral_logical_reads(tmp_path: Path) -> None:
    _reference, bam, _cram, first, second = generated_inputs(tmp_path)

    bam_result = fingerprint_input_artifact(artifact(bam, "BAM"), temporary_parent=tmp_path)
    fastq_result = fingerprint_input_artifact(
        artifact(first, "FASTQ_PAIR", mate_path=second), temporary_parent=tmp_path
    )

    assert not Path(f"{bam}.bai").exists()
    assert bam_result.logical.primary_record_count == 2
    assert fastq_result.logical.primary_record_count == 2
    assert bam_result.logical.named_sequence_sha256 == fastq_result.logical.named_sequence_sha256
    assert bam_result.logical.unnamed_sequence_sha256 == fastq_result.logical.unnamed_sequence_sha256
    assert bam_result.logical.alignment_sha256 != fastq_result.logical.alignment_sha256


def test_generated_cram_uses_verified_reference_and_restores_ref_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference, bam, cram, _first, _second = generated_inputs(tmp_path)
    reference_sha256 = hashlib.sha256(reference.read_bytes()).hexdigest()
    monkeypatch.setenv("REF_PATH", "https://synthetic.invalid/reference")

    cram_result = fingerprint_input_artifact(
        artifact(cram, "CRAM"),
        reference_path=reference,
        reference_sha256=reference_sha256,
        temporary_parent=tmp_path,
    )
    bam_result = fingerprint_input_artifact(artifact(bam, "BAM"), temporary_parent=tmp_path)

    assert cram_result.logical == bam_result.logical
    assert os.environ["REF_PATH"] == "https://synthetic.invalid/reference"


@pytest.mark.parametrize("declared_format", ["BAM", "CRAM"])
def test_generated_alignment_content_must_match_declared_format(tmp_path: Path, declared_format: str) -> None:
    reference, bam, cram, _first, _second = generated_inputs(tmp_path)
    path = cram if declared_format == "BAM" else bam
    kwargs = (
        {
            "reference_path": reference,
            "reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
        }
        if declared_format == "CRAM"
        else {}
    )

    with pytest.raises((ValueError, RuntimeError), match="declared format|failed to open"):
        fingerprint_input_artifact(
            artifact(path, declared_format),  # type: ignore[arg-type]
            temporary_parent=tmp_path,
            **kwargs,  # type: ignore[arg-type]
        )


def test_truncated_generated_bam_fails_without_an_index_or_partial_result(tmp_path: Path) -> None:
    _reference, bam, _cram, _first, _second = generated_inputs(tmp_path)
    truncated = tmp_path / "truncated.bam"
    truncated.write_bytes(bam.read_bytes()[:-20])

    with pytest.raises(RuntimeError, match="failed to (open|read) alignment"):
        fingerprint_input_artifact(artifact(truncated, "BAM"), temporary_parent=tmp_path)
