"""Strict local readers connect declared artifacts to duplicate-audit fingerprints."""

from __future__ import annotations

import gzip
import hashlib
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from vntyper.scripts.calibration_intake_contract import InputArtifact

pytestmark = pytest.mark.unit


def artifact(path: Path, *, artifact_format: str = "FASTQ_PAIR", mate_path: Path | None = None) -> InputArtifact:
    return InputArtifact(
        key="synthetic-artifact",
        specimen_key="synthetic-specimen",
        path=str(path),
        format=artifact_format,  # type: ignore[arg-type]
        mate_path=str(mate_path) if mate_path is not None else None,
        assembly="synthetic-assembly",
        assay_class="synthetic-assay",
        input_scope="full",
        preprocessing_id="synthetic-preprocessing-v1",
        replicate_group="synthetic-replicate",
        expected_sha256=None,
    )


def write_pair(
    first: Path,
    second: Path,
    *,
    first_name: str = "synthetic-read/1",
    second_name: str = "synthetic-read/2",
    sequence_1: str = "ACGT",
    sequence_2: str = "TGCA",
    quality_1: str = "IJKL",
    quality_2: str = "LMNO",
    gzip_encoded: bool = False,
) -> None:
    first_bytes = f"@{first_name}\n{sequence_1}\n+\n{quality_1}\n".encode("ascii")
    second_bytes = f"@{second_name}\n{sequence_2}\n+\n{quality_2}\n".encode("ascii")
    if gzip_encoded:
        first.write_bytes(gzip.compress(first_bytes, mtime=0))
        second.write_bytes(gzip.compress(second_bytes, mtime=0))
    else:
        first.write_bytes(first_bytes)
        second.write_bytes(second_bytes)


def test_plain_fastq_pair_hashes_exact_bytes_and_preserves_both_mates(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    directory = tmp_path / "paths with spaces"
    directory.mkdir()
    first = directory / "first reads"
    second = directory / "second reads"
    write_pair(first, second)
    declared = replace(
        artifact(first, mate_path=second), expected_sha256=hashlib.sha256(first.read_bytes()).hexdigest()
    )

    result = fingerprint_input_artifact(declared, temporary_parent=tmp_path)

    assert result.artifact_key == declared.key
    assert result.byte_sha256 == hashlib.sha256(first.read_bytes()).hexdigest()
    assert result.mate_byte_sha256 == hashlib.sha256(second.read_bytes()).hexdigest()
    assert result.logical.primary_record_count == 2
    assert result.logical.sequence_identity_reliable
    assert sorted(path.name for path in tmp_path.iterdir()) == ["paths with spaces"]


def test_typed_read_evidence_preserves_the_public_fingerprint_contract(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import (
        ArtifactReadEvidence,
        fingerprint_input_artifact,
        fingerprint_input_artifact_evidence,
    )

    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    write_pair(first, second)
    declared = artifact(first, mate_path=second)

    evidence = fingerprint_input_artifact_evidence(declared, temporary_parent=tmp_path)

    assert isinstance(evidence, ArtifactReadEvidence)
    assert evidence.verified_reference_sha256 is None
    assert evidence.fingerprint == fingerprint_input_artifact(declared, temporary_parent=tmp_path)


def test_symlink_aba_inputs_cannot_mix_original_hashes_with_alternate_logical_reads(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    original_first = tmp_path / "original-first.fastq"
    original_second = tmp_path / "original-second.fastq"
    alternate_first = tmp_path / "alternate-first.fastq"
    alternate_second = tmp_path / "alternate-second.fastq"
    write_pair(original_first, original_second, sequence_1="AAAA", sequence_2="TTTT")
    write_pair(alternate_first, alternate_second, sequence_1="CCCC", sequence_2="GGGG")
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    first.symlink_to(original_first)
    second.symlink_to(original_second)

    with pytest.raises(ValueError, match="symlink"):
        fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)


def test_gzip_fastq_is_detected_from_bytes_without_filename_guessing(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    first = tmp_path / "first.data"
    second = tmp_path / "second.data"
    write_pair(first, second, gzip_encoded=True)

    result = fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)

    assert result.logical.primary_record_count == 2


def test_bare_fastq_names_and_duplicate_occurrences_are_preserved(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    write_pair(first, second, first_name="synthetic-read", second_name="synthetic-read")
    once = fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)
    first.write_bytes(first.read_bytes() * 2)
    second.write_bytes(second.read_bytes() * 2)

    repeated = fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)

    assert once.logical.primary_record_count == 2
    assert repeated.logical.primary_record_count == 4
    assert once.logical.unnamed_sequence_sha256 != repeated.logical.unnamed_sequence_sha256


def test_casava_comments_and_matching_repeated_plus_identifiers_are_supported(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    first.write_text(
        "@synthetic-read/1 1:N:0:SYNTHETIC\nACGT\n+synthetic-read/1 1:N:0:SYNTHETIC\nIIII\n",
        encoding="ascii",
    )
    second.write_text(
        "@synthetic-read/2 2:N:0:SYNTHETIC\nTGCA\n+synthetic-read/2\nIIII\n",
        encoding="ascii",
    )

    result = fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)

    assert result.logical.primary_record_count == 2


@pytest.mark.parametrize(
    ("first_text", "second_text", "message"),
    [
        ("@synthetic/1\nACGT\n+\n", "@synthetic/2\nTGCA\n+\nIIII\n", "truncated"),
        ("@synthetic/1\nACGT\n+\nIII\n", "@synthetic/2\nTGCA\n+\nIIII\n", "quality"),
        ("@synthetic/1\nACGT\n-\nIIII\n", "@synthetic/2\nTGCA\n+\nIIII\n", "separator"),
        ("synthetic/1\nACGT\n+\nIIII\n", "@synthetic/2\nTGCA\n+\nIIII\n", "header"),
        ("@ synthetic/1\nACGT\n+\nIIII\n", "@synthetic/2\nTGCA\n+\nIIII\n", "header"),
        ("@one/1\nACGT\n+\nIIII\n", "@two/2\nTGCA\n+\nIIII\n", "names"),
        ("@synthetic/1\nACGT\n+\nIIII\n", "@synthetic\nTGCA\n+\nIIII\n", "names"),
        ("@synthetic/2\nACGT\n+\nIIII\n", "@synthetic/1\nTGCA\n+\nIIII\n", "names"),
        ("@synthetic\nACGT\n+\nIIII\n", "", "same number"),
        ("@synthetic\nACGZ\n+\nIIII\n", "@synthetic\nTGCA\n+\nIIII\n", "sequence"),
        ("@synthetic\nACGT\n+\n III\n", "@synthetic\nTGCA\n+\nIIII\n", "Phred"),
        (
            "@synthetic/1 2:N:0:SYNTHETIC\nACGT\n+\nIIII\n",
            "@synthetic/2 2:N:0:SYNTHETIC\nTGCA\n+\nIIII\n",
            "CASAVA",
        ),
        (
            "@synthetic/1 1:N:0:SYNTHETIC\nACGT\n+other/1\nIIII\n",
            "@synthetic/2 2:N:0:SYNTHETIC\nTGCA\n+\nIIII\n",
            "repeated",
        ),
    ],
)
def test_fastq_structure_and_pairing_fail_closed(
    tmp_path: Path, first_text: str, second_text: str, message: str
) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    first.write_text(first_text, encoding="ascii")
    second.write_text(second_text, encoding="ascii")

    with pytest.raises(ValueError, match=message):
        fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)


def test_truncated_gzip_fails_without_leaking_a_partial_fingerprint(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    first = tmp_path / "first.fastq.gz"
    second = tmp_path / "second.fastq.gz"
    write_pair(first, second, gzip_encoded=True)
    first.write_bytes(first.read_bytes()[:-5])

    with pytest.raises(RuntimeError, match="failed to read compressed FASTQ"):
        fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)

    assert sorted(path.name for path in tmp_path.iterdir()) == ["first.fastq.gz", "second.fastq.gz"]


def test_expected_hash_is_checked_before_records_are_opened(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_read_io

    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    write_pair(first, second)
    declared = replace(artifact(first, mate_path=second), expected_sha256="0" * 64)

    with (
        mock.patch.object(calibration_read_io, "_iter_fastq_pair", side_effect=AssertionError("opened")),
        pytest.raises(ValueError, match="expected SHA256"),
    ):
        calibration_read_io.fingerprint_input_artifact(declared, temporary_parent=tmp_path)


def test_source_byte_change_during_fastq_iteration_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from vntyper.scripts import calibration_read_io

    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    write_pair(first, second)
    original = calibration_read_io.fingerprint_read_records

    def changing_fingerprint(records, *, temporary_parent):
        result = original(records, temporary_parent=temporary_parent)
        first.write_bytes(first.read_bytes() + b"\n")
        return result

    monkeypatch.setattr(calibration_read_io, "fingerprint_read_records", changing_fingerprint)

    with pytest.raises(RuntimeError, match="changed during read fingerprinting"):
        calibration_read_io.fingerprint_input_artifact(artifact(first, mate_path=second), temporary_parent=tmp_path)


def alignment_record(**changes: object) -> SimpleNamespace:
    values: dict[str, object] = {
        "query_name": "synthetic-read/1",
        "query_sequence": "ACGT",
        "query_qualities": (40, 41, 42, 43),
        "flag": 65,
        "mapping_quality": 37,
        "reference_id": 0,
        "reference_name": "synthetic-contig",
        "reference_start": 10,
        "cigartuples": ((0, 4),),
        "next_reference_id": 0,
        "next_reference_name": "synthetic-contig",
        "next_reference_start": 30,
        "template_length": 24,
    }
    values.update(changes)
    return SimpleNamespace(**values)


class FakeAlignment:
    def __init__(self, records: list[object], *, is_bam: bool, is_cram: bool) -> None:
        self.records = records
        self.is_bam = is_bam
        self.is_cram = is_cram
        self.closed = False
        self.fetch_calls: list[dict[str, object]] = []

    def __enter__(self):
        return self

    def __exit__(self, *_args: object) -> None:
        self.closed = True

    def fetch(self, **kwargs: object):
        self.fetch_calls.append(kwargs)
        yield from self.records


def test_bam_uses_whole_file_iteration_and_maps_all_identity_fields(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_read_io
    from vntyper.scripts.calibration_read_identity import PrimaryReadRecord

    path = tmp_path / "input with spaces.bam"
    path.write_bytes(b"synthetic-bam-placeholder")
    primary = alignment_record()
    mate_two = alignment_record(query_name="synthetic-read/2", flag=129)
    unpaired = alignment_record(
        query_name="synthetic-unpaired",
        flag=0,
        next_reference_id=-1,
        next_reference_name=None,
        next_reference_start=-1,
        template_length=0,
    )
    secondary = alignment_record(flag=321)
    supplementary = alignment_record(flag=2113)
    alignment = FakeAlignment([primary, mate_two, unpaired, secondary, supplementary], is_bam=True, is_cram=False)
    captured: list[PrimaryReadRecord] = []
    original_fingerprint = calibration_read_io.fingerprint_read_records

    def capture(records, *, temporary_parent):
        captured.extend(records)
        return original_fingerprint(captured, temporary_parent=temporary_parent)

    with (
        mock.patch.object(calibration_read_io.pysam, "AlignmentFile", return_value=alignment) as opener,
        mock.patch.object(calibration_read_io, "fingerprint_read_records", side_effect=capture),
    ):
        result = calibration_read_io.fingerprint_input_artifact(
            artifact(path, artifact_format="BAM"), temporary_parent=tmp_path
        )

    opened_descriptor, mode = opener.call_args.args
    assert isinstance(opened_descriptor, int)
    assert mode == "rb"
    assert opener.call_args.kwargs == {"duplicate_filehandle": False}
    assert alignment.fetch_calls == [{"until_eof": True}]
    assert alignment.closed
    assert result.logical.primary_record_count == 3
    assert captured[0] == PrimaryReadRecord(
        name="synthetic-read/1",
        sequence="ACGT",
        qualities=(40, 41, 42, 43),
        mate=1,
        flags=65,
        mapping_quality=37,
        contig="synthetic-contig",
        position_zero_based=10,
        cigar=((0, 4),),
        mate_contig="synthetic-contig",
        mate_position_zero_based=30,
        template_length=24,
    )
    assert [record.mate for record in captured[:3]] == [1, 2, 0]


def test_declared_alignment_format_must_match_opened_content(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_read_io

    path = tmp_path / "declared.bam"
    path.write_bytes(b"synthetic-cram-placeholder")
    alignment = FakeAlignment([alignment_record()], is_bam=False, is_cram=True)

    with (
        mock.patch.object(calibration_read_io.pysam, "AlignmentFile", return_value=alignment),
        pytest.raises(ValueError, match="declared format"),
    ):
        calibration_read_io.fingerprint_input_artifact(artifact(path, artifact_format="BAM"), temporary_parent=tmp_path)

    assert alignment.closed


def test_alignment_iteration_failure_closes_handle_and_hides_input_name(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_read_io

    path = tmp_path / "private-identifier.bam"
    path.write_bytes(b"synthetic-bam-placeholder")
    alignment = FakeAlignment([alignment_record(), OSError("private-identifier.bam")], is_bam=True, is_cram=False)

    def failing_fetch(**kwargs: object):
        alignment.fetch_calls.append(kwargs)
        yield alignment.records[0]
        error = alignment.records[1]
        assert isinstance(error, OSError)
        raise error

    alignment.fetch = failing_fetch  # type: ignore[method-assign]
    with (
        mock.patch.object(calibration_read_io.pysam, "AlignmentFile", return_value=alignment),
        pytest.raises(RuntimeError, match="failed to read alignment") as caught,
    ):
        calibration_read_io.fingerprint_input_artifact(artifact(path, artifact_format="BAM"), temporary_parent=tmp_path)

    assert "private-identifier" not in str(caught.value)
    assert alignment.closed


def test_cram_requires_verified_reference_and_restores_resolution_environment(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_read_io

    path = tmp_path / "input.cram"
    reference = tmp_path / "reference.fa"
    path.write_bytes(b"synthetic-cram-placeholder")
    reference.write_text(">synthetic-contig\nACGT\n", encoding="ascii")
    reference_sha256 = hashlib.sha256(reference.read_bytes()).hexdigest()
    alignment = FakeAlignment([alignment_record()], is_bam=False, is_cram=True)

    snapshot: Path | None = None

    def open_snapshot(descriptor, mode, **kwargs):
        nonlocal snapshot
        assert isinstance(descriptor, int)
        assert mode == "rc"
        assert kwargs["duplicate_filehandle"] is False
        snapshot = Path(kwargs["reference_filename"])
        assert snapshot != reference
        assert snapshot.read_bytes() == reference.read_bytes()
        assert snapshot.stat().st_mode & 0o777 == 0o600
        assert snapshot.parent.stat().st_mode & 0o777 == 0o700
        return alignment

    with (
        mock.patch.object(calibration_read_io.pysam, "AlignmentFile", side_effect=open_snapshot) as opener,
        mock.patch.object(calibration_read_io, "pin_reference_resolution", return_value="previous") as pin,
        mock.patch.object(calibration_read_io, "restore_reference_resolution") as restore,
    ):
        evidence = calibration_read_io.fingerprint_input_artifact_evidence(
            artifact(path, artifact_format="CRAM"),
            reference_path=reference,
            reference_sha256=reference_sha256,
            temporary_parent=tmp_path,
        )
    assert opener.call_count == 1
    assert snapshot is not None
    pin.assert_called_once_with(
        {"cram": {"allow_ambient_reference_resolution": False, "local_ref_path": str(snapshot)}}
    )
    restore.assert_called_once_with("previous")
    assert evidence.fingerprint.logical.primary_record_count == 1
    assert evidence.verified_reference_sha256 == reference_sha256
    assert not snapshot.parent.exists()


def test_cram_restores_resolution_environment_when_content_is_not_cram(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_read_io

    path = tmp_path / "input.cram"
    reference = tmp_path / "reference.fa"
    path.write_bytes(b"synthetic-bam-placeholder")
    reference.write_text(">synthetic-contig\nACGT\n", encoding="ascii")
    reference_sha256 = hashlib.sha256(reference.read_bytes()).hexdigest()
    alignment = FakeAlignment([alignment_record()], is_bam=True, is_cram=False)

    with (
        mock.patch.object(calibration_read_io.pysam, "AlignmentFile", return_value=alignment),
        mock.patch.object(calibration_read_io, "pin_reference_resolution", return_value=None),
        mock.patch.object(calibration_read_io, "restore_reference_resolution") as restore,
        pytest.raises(ValueError, match="declared format"),
    ):
        calibration_read_io.fingerprint_input_artifact(
            artifact(path, artifact_format="CRAM"),
            reference_path=reference,
            reference_sha256=reference_sha256,
            temporary_parent=tmp_path,
        )

    restore.assert_called_once_with(None)
    assert alignment.closed


@pytest.mark.parametrize(
    "changes",
    [
        {"artifact_format": "CRAM", "reference_path": None, "reference_sha256": None},
        {"artifact_format": "CRAM", "reference_sha256": None},
        {"artifact_format": "CRAM", "reference_sha256": "0" * 64},
        {"artifact_format": "CRAM", "reference_path": Path("relative.fa"), "reference_sha256": "0" * 64},
        {"artifact_format": "BAM", "reference_sha256": "0" * 64},
    ],
)
def test_reference_contract_fails_before_alignment_open(tmp_path: Path, changes: dict[str, object]) -> None:
    from vntyper.scripts import calibration_read_io

    path = tmp_path / "input"
    reference = tmp_path / "reference"
    path.write_bytes(b"input")
    reference.write_bytes(b"reference")
    kwargs: dict[str, object] = {"reference_path": reference, "reference_sha256": None}
    kwargs.update({key: value for key, value in changes.items() if key != "artifact_format"})

    with (
        mock.patch.object(calibration_read_io.pysam, "AlignmentFile", side_effect=AssertionError("opened")),
        pytest.raises(ValueError),
    ):
        calibration_read_io.fingerprint_input_artifact(
            artifact(path, artifact_format=changes["artifact_format"]),  # type: ignore[arg-type]
            temporary_parent=tmp_path,
            **kwargs,  # type: ignore[arg-type]
        )


def test_directly_forged_artifact_and_nonregular_input_fail_closed(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    directory = tmp_path / "directory"
    directory.mkdir()
    malformed = replace(artifact(directory, mate_path=directory), expected_sha256=True)  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="artifact"):
        fingerprint_input_artifact(malformed, temporary_parent=tmp_path)

    valid_shape = replace(malformed, expected_sha256=None)
    with pytest.raises(ValueError, match="regular file"):
        fingerprint_input_artifact(valid_shape, temporary_parent=tmp_path)


@pytest.mark.parametrize(
    "malformed",
    [
        object(),
        replace(artifact(Path("input"), mate_path=Path("mate")), key=""),
        replace(artifact(Path("input"), mate_path=Path("mate")), format=[]),  # type: ignore[arg-type]
        replace(artifact(Path("input"), mate_path=Path("mate")), input_scope="unknown"),  # type: ignore[arg-type]
        replace(artifact(Path("input"), mate_path=Path("mate")), mate_path=None),
        replace(artifact(Path("input"), artifact_format="BAM"), mate_path="unexpected"),
    ],
)
def test_every_directly_forged_artifact_shape_is_revalidated(malformed: object, tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    with pytest.raises(ValueError, match="[Aa]rtifact|FASTQ"):
        fingerprint_input_artifact(malformed, temporary_parent=tmp_path)  # type: ignore[arg-type]


def test_invalid_temporary_parent_and_missing_input_fail_before_record_iteration(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_read_io import fingerprint_input_artifact

    missing = artifact(tmp_path / "missing", mate_path=tmp_path / "also-missing")
    with pytest.raises(ValueError, match="temporary parent"):
        fingerprint_input_artifact(missing, temporary_parent=tmp_path / "missing-parent")
    with pytest.raises(RuntimeError, match="could not be read"):
        fingerprint_input_artifact(missing, temporary_parent=tmp_path)
