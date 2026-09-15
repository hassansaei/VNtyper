"""Strict local input readers for calibration duplicate-audit fingerprints."""

from __future__ import annotations

import gzip
import hashlib
import logging
import os
import re
import stat
from collections.abc import Generator
from contextlib import ExitStack, closing
from dataclasses import dataclass
from pathlib import Path
from typing import IO, NoReturn, cast

import pysam

from vntyper.scripts.calibration_identity import ArtifactFingerprint
from vntyper.scripts.calibration_intake_contract import InputArtifact
from vntyper.scripts.calibration_read_fingerprints import LogicalReadFingerprint, fingerprint_read_records
from vntyper.scripts.calibration_read_identity import PrimaryReadRecord, canonical_sequence_name
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution, restore_reference_resolution

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_GZIP_MAGIC = b"\x1f\x8b"
_BUFFER_SIZE = 1024 * 1024
_FORMATS = frozenset({"BAM", "CRAM", "FASTQ_PAIR"})
_SCOPES = frozenset({"full", "regional"})
_CASAVA_COMMENT = re.compile(r"(?P<mate>[12]):[YN]:[0-9]+:[!-~]+\Z")


@dataclass(frozen=True)
class _FileEvidence:
    device: int
    inode: int
    size: int
    modified_ns: int
    changed_ns: int
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def fingerprint_input_artifact(
    artifact: InputArtifact,
    *,
    reference_path: Path | None = None,
    reference_sha256: str | None = None,
    temporary_parent: Path | None = None,
) -> ArtifactFingerprint:
    """Read and fingerprint one declared BAM, CRAM, or paired FASTQ artifact.

    BAM and CRAM inputs are read sequentially with ``until_eof=True`` and need
    no index. CRAM decoding is restricted to one verified local reference.
    Physical bytes and stable file identity are checked before and after logical
    fingerprinting. Exceptions never include artifact paths or read names.

    Args:
        artifact: Strict typed intake artifact declaration.
        reference_path: Pinned local FASTA required only for CRAM.
        reference_sha256: Expected lowercase SHA256 of the CRAM reference.
        temporary_parent: Optional parent for digest-only external-sort files.

    Returns:
        Physical byte hashes and logical multiset fingerprints.

    Raises:
        ValueError: If declarations, formats, FASTQ structure, or hashes conflict.
        RuntimeError: If files change or an input cannot be read completely.
    """
    _validate_artifact(artifact)
    _validate_temporary_parent(temporary_parent)
    input_path = Path(artifact.path)
    mate_path = Path(artifact.mate_path) if artifact.mate_path is not None else None
    reference = _validate_reference_contract(artifact.format, reference_path, reference_sha256)

    input_evidence = _stable_file_evidence(input_path, "calibration input file")
    if artifact.expected_sha256 is not None and input_evidence.sha256 != artifact.expected_sha256:
        _fail("calibration input file differs from its expected SHA256")
    mate_evidence = _stable_file_evidence(mate_path, "calibration mate file") if mate_path is not None else None
    reference_evidence = (
        _stable_file_evidence(reference, "calibration CRAM reference") if reference is not None else None
    )
    if reference_evidence is not None and reference_evidence.sha256 != reference_sha256:
        _fail("calibration CRAM reference differs from its expected SHA256")

    if artifact.format == "FASTQ_PAIR":
        records = _iter_fastq_pair(input_path, cast(Path, mate_path))
        with closing(records):
            logical = fingerprint_read_records(records, temporary_parent=temporary_parent)
    else:
        logical = _fingerprint_alignment(
            input_path,
            artifact.format,
            reference,
            temporary_parent=temporary_parent,
        )

    _assert_unchanged(input_path, input_evidence, "calibration input file")
    if mate_path is not None and mate_evidence is not None:
        _assert_unchanged(mate_path, mate_evidence, "calibration mate file")
    if reference is not None and reference_evidence is not None:
        _assert_unchanged(reference, reference_evidence, "calibration CRAM reference")
    return ArtifactFingerprint(
        artifact.key,
        input_evidence.sha256,
        mate_evidence.sha256 if mate_evidence is not None else None,
        logical,
    )


def _validate_artifact(artifact: InputArtifact) -> None:
    if not isinstance(artifact, InputArtifact):
        _fail("calibration read adapter requires an InputArtifact")
    for value in (
        artifact.key,
        artifact.specimen_key,
        artifact.path,
        artifact.assembly,
        artifact.assay_class,
        artifact.preprocessing_id,
        artifact.replicate_group,
    ):
        if not isinstance(value, str) or not value:
            _fail("calibration read adapter artifact strings must be non-empty")
    if (
        not isinstance(artifact.format, str)
        or artifact.format not in _FORMATS
        or not isinstance(artifact.input_scope, str)
        or artifact.input_scope not in _SCOPES
    ):
        _fail("calibration read adapter artifact contract is invalid")
    if artifact.format == "FASTQ_PAIR":
        if not isinstance(artifact.mate_path, str) or not artifact.mate_path:
            _fail("calibration FASTQ pair requires a mate file")
    elif artifact.mate_path is not None:
        _fail("calibration alignment artifact forbids a mate file")
    if artifact.expected_sha256 is not None and (
        not isinstance(artifact.expected_sha256, str) or _SHA256.fullmatch(artifact.expected_sha256) is None
    ):
        _fail("calibration read adapter artifact expected SHA256 is invalid")


def _validate_temporary_parent(temporary_parent: Path | None) -> None:
    if temporary_parent is not None and (not isinstance(temporary_parent, Path) or not temporary_parent.is_dir()):
        _fail("calibration read adapter temporary parent must be an existing directory")


def _validate_reference_contract(
    artifact_format: str,
    reference_path: Path | None,
    reference_sha256: str | None,
) -> Path | None:
    if artifact_format != "CRAM":
        if reference_path is not None or reference_sha256 is not None:
            _fail("calibration reference arguments are permitted only for CRAM")
        return None
    if not isinstance(reference_path, Path) or not reference_path.is_absolute():
        _fail("calibration CRAM requires an absolute pinned local reference path")
    if not isinstance(reference_sha256, str) or _SHA256.fullmatch(reference_sha256) is None:
        _fail("calibration CRAM requires an expected lowercase reference SHA256")
    return reference_path


def _stable_file_evidence(path: Path, label: str) -> _FileEvidence:
    try:
        path_before = path.stat()
    except OSError:
        raise RuntimeError(f"{label} could not be read") from None
    if not stat.S_ISREG(path_before.st_mode):
        _fail(f"{label} must be a regular file")
    try:
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode):
                _fail(f"{label} must be a regular file")
            digest = hashlib.sha256()
            for chunk in iter(lambda: handle.read(_BUFFER_SIZE), b""):
                digest.update(chunk)
            after = os.fstat(handle.fileno())
        path_after = path.stat()
    except OSError:
        raise RuntimeError(f"{label} could not be read") from None
    before_identity = _stat_identity(before)
    if (
        before_identity != _stat_identity(path_before)
        or before_identity != _stat_identity(after)
        or before_identity != _stat_identity(path_after)
    ):
        raise RuntimeError(f"{label} changed while its digest was computed")
    return _FileEvidence(*before_identity, digest.hexdigest())


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def _assert_unchanged(path: Path, expected: _FileEvidence, label: str) -> None:
    if _stable_file_evidence(path, label) != expected:
        raise RuntimeError(f"{label} changed during read fingerprinting")


def _open_fastq(stack: ExitStack, path: Path) -> IO[bytes]:
    raw = stack.enter_context(path.open("rb"))
    magic = raw.read(2)
    raw.seek(0)
    if magic == _GZIP_MAGIC:
        return cast(IO[bytes], stack.enter_context(gzip.GzipFile(fileobj=raw, mode="rb")))
    return raw


def _fastq_line(handle: IO[bytes]) -> bytes | None:
    raw = handle.readline()
    if not raw:
        return None
    if not raw.endswith(b"\n"):
        _fail("calibration FASTQ record is truncated")
    line = raw[:-1]
    if line.endswith(b"\r"):
        line = line[:-1]
    return line


def _read_fastq_record(handle: IO[bytes], mate: int) -> tuple[str, str, tuple[int, ...]] | None:
    header = _fastq_line(handle)
    if header is None:
        return None
    lines = [header]
    for _ in range(3):
        line = _fastq_line(handle)
        if line is None:
            _fail("calibration FASTQ record is truncated")
        lines.append(line)
    header, sequence, separator, quality = lines
    if not header.startswith(b"@") or len(header) == 1:
        _fail("calibration FASTQ header is invalid")
    try:
        full_header = header[1:].decode("ascii")
        decoded_sequence = sequence.decode("ascii")
    except UnicodeDecodeError:
        _fail("calibration FASTQ header and sequence must be ASCII")
    fields = full_header.split(maxsplit=1)
    if not fields or full_header[0].isspace():
        _fail("calibration FASTQ header is invalid")
    name = fields[0]
    if len(fields) == 2:
        first_comment_field = fields[1].split(maxsplit=1)[0]
        casava = _CASAVA_COMMENT.fullmatch(first_comment_field)
        if first_comment_field.startswith(("1:", "2:")) and casava is None:
            _fail("calibration FASTQ CASAVA comment is malformed")
        if casava is not None and int(casava.group("mate")) != mate:
            _fail("calibration FASTQ CASAVA mate field conflicts with its stream")
    if not separator.startswith(b"+"):
        _fail("calibration FASTQ separator must begin with '+'")
    repeated = separator[1:]
    if repeated:
        try:
            repeated_text = repeated.decode("ascii")
        except UnicodeDecodeError:
            _fail("calibration FASTQ repeated identifier must be ASCII")
        if repeated_text not in {name, full_header}:
            _fail("calibration FASTQ repeated identifier differs from its header")
    if len(quality) != len(sequence):
        _fail("calibration FASTQ quality length differs from sequence length")
    if any(value < 33 or value > 126 for value in quality):
        _fail("calibration FASTQ quality must use Phred+33 bytes")
    return name, decoded_sequence, tuple(value - 33 for value in quality)


def _validate_paired_names(first: str, second: str) -> None:
    first_has_suffix = first.endswith(("/1", "/2"))
    second_has_suffix = second.endswith(("/1", "/2"))
    if first_has_suffix or second_has_suffix:
        if not first.endswith("/1") or not second.endswith("/2"):
            _fail("calibration FASTQ mate names use inconsistent terminal suffixes")
        if canonical_sequence_name(first, 1) != canonical_sequence_name(second, 2):
            _fail("calibration FASTQ mate names differ")
    elif first != second:
        _fail("calibration FASTQ mate names differ")


def _fastq_record(name: str, sequence: str, qualities: tuple[int, ...], mate: int) -> PrimaryReadRecord:
    return PrimaryReadRecord(
        name=name,
        sequence=sequence,
        qualities=qualities,
        mate=mate,
        flags=69 if mate == 1 else 133,
        mapping_quality=0,
        contig=None,
        position_zero_based=-1,
        cigar=(),
        mate_contig=None,
        mate_position_zero_based=-1,
        template_length=0,
    )


def _iter_fastq_pair(first_path: Path, second_path: Path) -> Generator[PrimaryReadRecord, None, None]:
    try:
        with ExitStack() as stack:
            first_handle = _open_fastq(stack, first_path)
            second_handle = _open_fastq(stack, second_path)
            while True:
                first = _read_fastq_record(first_handle, 1)
                second = _read_fastq_record(second_handle, 2)
                if first is None and second is None:
                    return
                if first is None or second is None:
                    _fail("calibration FASTQ mates must contain the same number of records")
                _validate_paired_names(first[0], second[0])
                yield _fastq_record(*first, mate=1)
                yield _fastq_record(*second, mate=2)
    except (gzip.BadGzipFile, EOFError):
        raise RuntimeError("failed to read compressed FASTQ input completely") from None
    except OSError:
        raise RuntimeError("failed to read calibration FASTQ input") from None


def _mate_ordinal(flags: int) -> int:
    first = bool(flags & 0x40)
    second = bool(flags & 0x80)
    if first and not second:
        return 1
    if second and not first:
        return 2
    return 0


def _primary_record(record: pysam.AlignedSegment) -> PrimaryReadRecord:
    qualities = record.query_qualities
    cigar = record.cigartuples
    return PrimaryReadRecord(
        name=cast(str, record.query_name),
        sequence=record.query_sequence,
        qualities=tuple(qualities) if qualities is not None else None,
        mate=_mate_ordinal(record.flag),
        flags=record.flag,
        mapping_quality=record.mapping_quality,
        contig=record.reference_name if record.reference_id >= 0 else None,
        position_zero_based=record.reference_start,
        cigar=tuple(cigar) if cigar is not None else (),
        mate_contig=record.next_reference_name if record.next_reference_id >= 0 else None,
        mate_position_zero_based=record.next_reference_start,
        template_length=record.template_length,
    )


def _fingerprint_alignment(
    path: Path,
    artifact_format: str,
    reference_path: Path | None,
    *,
    temporary_parent: Path | None,
) -> LogicalReadFingerprint:
    previous_reference_path: str | None = None
    if artifact_format == "CRAM":
        if reference_path is None:
            _fail("calibration CRAM requires a pinned local reference path")
        previous_reference_path = pin_reference_resolution(
            {"cram": {"allow_ambient_reference_resolution": False, "local_ref_path": str(reference_path)}}
        )
    try:
        try:
            if artifact_format == "CRAM":
                alignment = pysam.AlignmentFile(str(path), "rc", reference_filename=str(reference_path))
            else:
                alignment = pysam.AlignmentFile(str(path), "rb")
        except (OSError, ValueError):
            raise RuntimeError("failed to open alignment input") from None
        with alignment:
            if (artifact_format == "BAM" and not alignment.is_bam) or (
                artifact_format == "CRAM" and not alignment.is_cram
            ):
                _fail("calibration alignment content differs from its declared format")
            try:
                records = (_primary_record(record) for record in alignment.fetch(until_eof=True))
                return fingerprint_read_records(records, temporary_parent=temporary_parent)
            except OSError:
                raise RuntimeError("failed to read alignment input completely") from None
    finally:
        if artifact_format == "CRAM":
            restore_reference_resolution(previous_reference_path)
