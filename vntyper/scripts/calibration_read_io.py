"""Strict local input readers for calibration duplicate-audit fingerprints."""

from __future__ import annotations

import gzip
import logging
import os
import re
import stat
from collections.abc import Generator, Iterator
from contextlib import ExitStack, closing, contextmanager, suppress
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import IO, BinaryIO, NoReturn, cast

import pysam

from vntyper.scripts.calibration_identity import ArtifactFingerprint
from vntyper.scripts.calibration_intake_contract import InputArtifact
from vntyper.scripts.calibration_read_fingerprints import LogicalReadFingerprint, fingerprint_read_records
from vntyper.scripts.calibration_read_identity import PrimaryReadRecord, canonical_sequence_name
from vntyper.scripts.calibration_read_sources import PinnedReadSource
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution, restore_reference_resolution

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_GZIP_MAGIC = b"\x1f\x8b"
_FORMATS = frozenset({"BAM", "CRAM", "FASTQ_PAIR"})
_SCOPES = frozenset({"full", "regional"})
_CASAVA_COMMENT = re.compile(r"(?P<mate>[12]):[YN]:[0-9]+:[!-~]+\Z")


@dataclass(frozen=True)
class ArtifactReadEvidence:
    """A read fingerprint plus the verified CRAM reference digest, when used."""

    fingerprint: ArtifactFingerprint
    verified_reference_sha256: str | None


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
    return fingerprint_input_artifact_evidence(
        artifact,
        reference_path=reference_path,
        reference_sha256=reference_sha256,
        temporary_parent=temporary_parent,
    ).fingerprint


def fingerprint_input_artifact_evidence(
    artifact: InputArtifact,
    *,
    reference_path: Path | None = None,
    reference_sha256: str | None = None,
    temporary_parent: Path | None = None,
) -> ArtifactReadEvidence:
    """Fingerprint one artifact and return its verified reference evidence.

    Args:
        artifact: Strict typed intake artifact declaration.
        reference_path: Pinned local FASTA required only for CRAM.
        reference_sha256: Expected lowercase SHA256 of the CRAM reference.
        temporary_parent: Optional parent for private reference and sort files.

    Returns:
        The artifact fingerprint and the observed CRAM reference digest.

    Raises:
        ValueError: If declarations, formats, FASTQ structure, or hashes conflict.
        RuntimeError: If files change or an input cannot be read completely.
    """
    _validate_artifact(artifact)
    _validate_temporary_parent(temporary_parent)
    input_path = Path(artifact.path)
    mate_path = Path(artifact.mate_path) if artifact.mate_path is not None else None
    reference = _validate_reference_contract(artifact.format, reference_path, reference_sha256)

    with ExitStack() as stack:
        input_source = stack.enter_context(PinnedReadSource.open(input_path, "calibration input file"))
        mate_source = (
            stack.enter_context(PinnedReadSource.open(mate_path, "calibration mate file"))
            if mate_path is not None
            else None
        )
        reference_source = (
            stack.enter_context(PinnedReadSource.open(reference, "calibration CRAM reference"))
            if reference is not None
            else None
        )
        input_sha256 = input_source.digest()
        if artifact.expected_sha256 is not None and input_sha256 != artifact.expected_sha256:
            _fail("calibration input file differs from its expected SHA256")
        mate_sha256 = mate_source.digest() if mate_source is not None else None
        observed_reference_sha256: str | None = None
        if reference_source is None:
            logical = _fingerprint_source(
                artifact.format,
                input_source,
                mate_source,
                None,
                temporary_parent=temporary_parent,
            )
        else:
            with _private_reference_snapshot(
                reference_source,
                cast(str, reference_sha256),
                temporary_parent=temporary_parent,
            ) as snapshot:
                reference_snapshot, observed_reference_sha256 = snapshot
                logical = _fingerprint_source(
                    artifact.format,
                    input_source,
                    mate_source,
                    reference_snapshot,
                    temporary_parent=temporary_parent,
                )
        input_source.verify()
        if mate_source is not None:
            mate_source.verify()
        if reference_source is not None:
            reference_source.verify()
        fingerprint = ArtifactFingerprint(artifact.key, input_sha256, mate_sha256, logical)
        return ArtifactReadEvidence(fingerprint, observed_reference_sha256)


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


def _open_fastq(stack: ExitStack, source: PinnedReadSource) -> IO[bytes]:
    raw = stack.enter_context(os.fdopen(source.duplicate_descriptor(), "rb"))
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


def _iter_fastq_pair(
    first_source: PinnedReadSource, second_source: PinnedReadSource
) -> Generator[PrimaryReadRecord, None, None]:
    try:
        with ExitStack() as stack:
            first_handle = _open_fastq(stack, first_source)
            second_handle = _open_fastq(stack, second_source)
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


def _write_private_reference(source: PinnedReadSource, target: Path) -> str:
    descriptor = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0), 0o600)
    try:
        os.fchmod(descriptor, 0o600)
        with os.fdopen(descriptor, "wb", closefd=False) as handle:
            return source.copy_and_digest(cast(BinaryIO, handle))
    finally:
        os.close(descriptor)


def _make_private_fasta_index(reference_path: Path) -> None:
    try:
        pysam.faidx(str(reference_path))  # type: ignore[attr-defined]
    except (OSError, ValueError):
        raise RuntimeError("failed to index the verified calibration CRAM reference snapshot") from None
    for child in reference_path.parent.iterdir():
        metadata = child.lstat()
        if not stat.S_ISREG(metadata.st_mode):
            raise RuntimeError("verified calibration CRAM reference indexing produced an invalid artifact")
        child.chmod(0o600)


@contextmanager
def _private_reference_snapshot(
    source: PinnedReadSource,
    expected_sha256: str,
    *,
    temporary_parent: Path | None,
) -> Iterator[tuple[Path, str]]:
    try:
        with TemporaryDirectory(prefix=".vntyper-cram-reference-", dir=temporary_parent) as raw_directory:
            directory = Path(raw_directory)
            directory.chmod(0o700)
            reference_path = directory / "reference.fa"
            observed_sha256 = _write_private_reference(source, reference_path)
            if observed_sha256 != expected_sha256:
                _fail("calibration CRAM reference differs from its expected SHA256")
            _make_private_fasta_index(reference_path)
            yield reference_path, observed_sha256
    except OSError:
        raise RuntimeError("failed to create a private calibration CRAM reference snapshot") from None


def _fingerprint_source(
    artifact_format: str,
    input_source: PinnedReadSource,
    mate_source: PinnedReadSource | None,
    reference_snapshot: Path | None,
    *,
    temporary_parent: Path | None,
) -> LogicalReadFingerprint:
    if artifact_format == "FASTQ_PAIR":
        if mate_source is None:
            _fail("calibration FASTQ pair requires a mate file")
        records = _iter_fastq_pair(input_source, mate_source)
        with closing(records):
            return fingerprint_read_records(records, temporary_parent=temporary_parent)
    return _fingerprint_alignment(
        input_source,
        artifact_format,
        reference_snapshot,
        temporary_parent=temporary_parent,
    )


def _close_descriptor_if_open(descriptor: int) -> None:
    with suppress(OSError):
        os.fstat(descriptor)
        os.close(descriptor)


def _fingerprint_alignment(
    source: PinnedReadSource,
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
        descriptor = source.duplicate_descriptor()
        try:
            if artifact_format == "CRAM":
                alignment = pysam.AlignmentFile(
                    descriptor,
                    "rc",
                    reference_filename=str(reference_path),
                    duplicate_filehandle=False,
                )
            else:
                alignment = pysam.AlignmentFile(descriptor, "rb", duplicate_filehandle=False)
        except (OSError, ValueError):
            _close_descriptor_if_open(descriptor)
            raise RuntimeError("failed to open alignment input") from None
        try:
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
            _close_descriptor_if_open(descriptor)
    finally:
        if artifact_format == "CRAM":
            restore_reference_resolution(previous_reference_path)
