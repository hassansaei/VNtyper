"""Strict depth and fragment-support I/O for VNTR length measurement."""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Literal, cast

import pysam

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_annotation import (
    Interval,
    LengthAnnotation,
    decode_length_annotation,
    encode_length_annotation,
)
from vntyper.scripts.length_feature_provenance import (
    LengthFeatureContext,
    QueryInterval,
    encode_length_feature_context,
)
from vntyper.scripts.length_features import DepthPosition
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution, restore_reference_resolution

_EXCLUDED_ALIGNMENT_FLAGS = 0x4 | 0x100 | 0x200 | 0x400
_SAMTOOLS_VERSION_PATTERN = re.compile(r"^samtools ([^\s]+)$")
_HTSLIB_VERSION_PATTERN = re.compile(r"^Using htslib ([^\s]+)$")
_UNSIGNED_INTEGER_PATTERN = re.compile(r"^(0|[1-9][0-9]*)$")
_BUFFER_SIZE = 1024 * 1024


@dataclass(frozen=True)
class FragmentPositionEvidence:
    """Selected aligned-base contributions and fragment identities at one position."""

    base_contributions: int
    fragment_ids: frozenset[str]

    def __post_init__(self) -> None:
        if not _is_int(self.base_contributions) or self.base_contributions < 0:
            raise ValueError("fragment base contributions must be a non-negative integer")
        if not isinstance(self.fragment_ids, frozenset) or any(
            not isinstance(value, str) or not value or value.strip() != value for value in self.fragment_ids
        ):
            raise ValueError("fragment identities must be a frozenset of non-empty trimmed strings")
        if len(self.fragment_ids) > self.base_contributions:
            raise ValueError("distinct fragment identities cannot exceed base contributions")


@dataclass(frozen=True)
class AlignmentEvidence:
    """Validated header identity and per-position evidence from alignment records."""

    contig_length: int
    assembly: str | None
    positions: Mapping[int, FragmentPositionEvidence]

    def __post_init__(self) -> None:
        if not _is_int(self.contig_length) or self.contig_length <= 0:
            raise ValueError("alignment header contig length must be a positive integer")
        if self.assembly is not None:
            _text(self.assembly, "alignment header assembly")
        if not isinstance(self.positions, Mapping):
            raise ValueError("alignment evidence positions must be a mapping")
        copied = dict(self.positions)
        if any(not _is_int(position) or position < 0 for position in copied):
            raise ValueError("alignment evidence positions must be non-negative integers")
        validated: dict[int, FragmentPositionEvidence] = {}
        for position, evidence in copied.items():
            if not isinstance(evidence, FragmentPositionEvidence):
                raise ValueError("alignment evidence must contain FragmentPositionEvidence values")
            validated[position] = FragmentPositionEvidence(evidence.base_contributions, evidence.fragment_ids)
        object.__setattr__(self, "positions", MappingProxyType(validated))


@dataclass(frozen=True)
class _FileIdentity:
    device: int
    inode: int
    size: int
    modified_ns: int
    sha256: str


def build_samtools_depth_argv(
    samtools_path: Path,
    input_path: Path,
    reference_path: Path,
    bed_path: Path,
) -> tuple[str, ...]:
    """Build the exact shell-free argv for the frozen base-depth policy.

    Args:
        samtools_path: Absolute path to the pinned samtools executable.
        input_path: Absolute BAM or CRAM path.
        reference_path: Absolute pinned local reference FASTA path.
        bed_path: Absolute path to the exact zero-based interval union BED.

    Returns:
        Immutable argv suitable for :func:`subprocess.run` with ``shell=False``.

    Raises:
        ValueError: If any path is not an absolute ``Path``.
    """
    for value, label in (
        (samtools_path, "samtools executable"),
        (input_path, "length depth input"),
        (reference_path, "length depth reference"),
        (bed_path, "length depth BED"),
    ):
        _absolute_path(value, label)
    return (
        str(samtools_path),
        "depth",
        "-a",
        "-q",
        "0",
        "-Q",
        "0",
        "-b",
        str(bed_path),
        "--reference",
        str(reference_path),
        str(input_path),
    )


def parse_samtools_depth(
    stdout: str,
    original_contig: str,
    expected_positions: set[int],
    fragment_evidence: Mapping[int, FragmentPositionEvidence],
) -> tuple[DepthPosition, ...]:
    """Parse strict samtools depth output and cross-check record-level contributions.

    Args:
        stdout: Complete UTF-8 text emitted by ``samtools depth``.
        original_contig: Exact input contig spelling required in every row.
        expected_positions: Exact zero-based annotation position union.
        fragment_evidence: Independent record-level contribution and fragment evidence.

    Returns:
        Position-sorted immutable depth observations.

    Raises:
        ValueError: If output or independent evidence is incomplete, malformed, or inconsistent.
    """
    if not isinstance(stdout, str):
        raise ValueError("samtools depth output must be text")
    _text(original_contig, "samtools depth contig")
    if not isinstance(expected_positions, set) or any(
        not _is_int(position) or position < 0 for position in expected_positions
    ):
        raise ValueError("expected depth positions must be a set of non-negative integers")
    if not isinstance(fragment_evidence, Mapping):
        raise ValueError("fragment evidence must be a position mapping")

    parsed: dict[int, int] = {}
    for line in stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError("samtools depth rows must contain exactly three tab-separated fields")
        contig, position_text, depth_text = fields
        if contig != original_contig:
            raise ValueError("samtools depth row contig does not match the original contig")
        if (
            _UNSIGNED_INTEGER_PATTERN.fullmatch(position_text) is None
            or _UNSIGNED_INTEGER_PATTERN.fullmatch(depth_text) is None
        ):
            raise ValueError("samtools depth position and depth must be canonical non-negative integers")
        position_one_based = int(position_text)
        depth = int(depth_text)
        if position_one_based <= 0:
            raise ValueError("samtools depth position must be a positive one-based integer")
        position = position_one_based - 1
        if position not in expected_positions:
            raise ValueError("samtools depth contains an out-of-range position")
        if position in parsed:
            raise ValueError("samtools depth contains a duplicate position")
        parsed[position] = depth

    missing = expected_positions - set(parsed)
    if missing:
        raise ValueError("samtools depth output is missing an expected position")
    if set(fragment_evidence) != expected_positions:
        raise ValueError("fragment evidence does not cover the exact expected position union")

    result: list[DepthPosition] = []
    for position in sorted(expected_positions):
        evidence = fragment_evidence[position]
        if not isinstance(evidence, FragmentPositionEvidence):
            raise ValueError("fragment evidence must contain FragmentPositionEvidence values")
        evidence = FragmentPositionEvidence(evidence.base_contributions, evidence.fragment_ids)
        depth = parsed[position]
        if evidence.base_contributions != depth:
            raise ValueError("fragment reader base contribution mismatch with samtools depth")
        result.append(
            DepthPosition(
                contig=original_contig,
                position_zero_based=position,
                depth=depth,
                supporting_fragment_ids=tuple(sorted(evidence.fragment_ids)),
            )
        )
    return tuple(result)


def read_length_depth(
    input_path: Path,
    reference_path: Path,
    annotation: LengthAnnotation,
    context: LengthFeatureContext,
    samtools_path: Path,
) -> tuple[DepthPosition, ...]:
    """Read authoritative base depth and independently cross-checked fragment support.

    Args:
        input_path: Absolute indexed BAM or CRAM path.
        reference_path: Absolute pinned local reference FASTA path.
        annotation: Validated length annotation.
        context: Validated measurement context with exact tool identities and policy.
        samtools_path: Absolute pinned samtools executable path.

    Returns:
        One depth observation for every position in the exact annotation union.

    Raises:
        ValueError: If contracts, files, headers, tools, or output are incompatible.
        RuntimeError: If a tool/read fails or an input changes during measurement.
    """
    for value, label in (
        (input_path, "length depth input"),
        (reference_path, "length depth reference"),
        (samtools_path, "samtools executable"),
    ):
        _absolute_path(value, label)
        if not value.is_file():
            raise ValueError(f"{label} must be an existing regular file")
    if not os.access(samtools_path, os.X_OK):
        raise ValueError("samtools executable must be executable")
    _validate_contract(annotation, context)

    input_identity = _stable_file_identity(input_path, "length depth input")
    reference_identity = _stable_file_identity(reference_path, "length depth reference")
    if input_identity.sha256 != context.input_sha256:
        raise ValueError("length depth input digest mismatch")
    if reference_identity.sha256 != context.reference_fasta_sha256:
        raise ValueError("length depth reference digest mismatch")

    expected_positions = _annotation_positions(annotation)
    if not expected_positions:
        raise ValueError("length depth annotation has no positions to query")
    required_end = max(expected_positions) + 1
    fasta_length = _fasta_contig_length(reference_path, context.original_contig)
    if fasta_length < required_end:
        raise ValueError("reference FASTA contig length does not cover the queried interval union")
    if _samtools_revision(samtools_path) != context.counting_policy.samtools_revision:
        raise ValueError("samtools revision does not match the measurement context")
    _validate_fragment_reader(context)

    is_cram = _has_cram_magic(input_path)
    previous_ref_path: str | None = None
    if is_cram:
        previous_ref_path = pin_reference_resolution(
            {"cram": {"allow_ambient_reference_resolution": False, "local_ref_path": str(reference_path)}}
        )
    try:
        alignment_evidence = _read_alignment_evidence(
            input_path,
            reference_path,
            context.original_contig,
            context.counting_policy.queried_intervals,
            expected_positions,
            is_cram=is_cram,
        )
        _validate_alignment_header(alignment_evidence, context, required_end, fasta_length)
        with tempfile.TemporaryDirectory(prefix=".vntyper-length-depth-") as temporary_name:
            bed_path = Path(temporary_name) / "exact-union.bed"
            _write_bed(bed_path, context.original_contig, context.counting_policy.queried_intervals)
            argv = build_samtools_depth_argv(samtools_path, input_path, reference_path, bed_path)
            completed = subprocess.run(argv, capture_output=True, check=False, text=True, env=None, shell=False)
            if completed.returncode != 0:
                raise RuntimeError("samtools depth failed for length measurement")
            depths = parse_samtools_depth(
                completed.stdout,
                context.original_contig,
                expected_positions,
                alignment_evidence.positions,
            )
    finally:
        if is_cram:
            restore_reference_resolution(previous_ref_path)

    _assert_file_unchanged(input_path, input_identity, "length depth input")
    _assert_file_unchanged(reference_path, reference_identity, "length depth reference")
    return depths


def _validate_contract(annotation: LengthAnnotation, context: LengthFeatureContext) -> None:
    if not isinstance(annotation, LengthAnnotation):
        raise ValueError("length depth annotation must be a LengthAnnotation")
    if not isinstance(context, LengthFeatureContext):
        raise ValueError("length depth context must be a LengthFeatureContext")
    encoded_annotation = encode_length_annotation(annotation)
    decoded_annotation = decode_length_annotation(encoded_annotation)
    if annotation.sha256 != decoded_annotation.sha256:
        raise ValueError("length depth annotation digest does not match canonical content")
    encode_length_feature_context(context)
    if context.assembly != annotation.assembly:
        raise ValueError("length depth assembly mismatch")
    if context.original_contig not in annotation.accepted_contigs:
        raise ValueError("length depth contig mismatch")
    if context.reference_fasta_sha256 != annotation.reference_fasta_sha256:
        raise ValueError("length depth reference digest mismatch")
    if context.annotation_sha256 != annotation.sha256:
        raise ValueError("length depth annotation digest mismatch")
    positions = _annotation_positions(annotation)
    if context.counting_policy.queried_intervals != _positions_to_intervals(positions):
        raise ValueError("length depth queried interval union does not match annotation positions")


def _annotation_positions(annotation: LengthAnnotation) -> set[int]:
    intervals: list[Interval] = []
    for repeated_intervals in (annotation.core, annotation.invariant):
        intervals.extend(repeated_intervals or ())
    intervals.extend(
        interval
        for interval in (annotation.array, annotation.left_flank, annotation.right_flank)
        if interval is not None
    )
    return {position for interval in intervals for position in range(interval.start, interval.end)}


def _positions_to_intervals(positions: set[int]) -> tuple[QueryInterval, ...]:
    if not positions:
        return ()
    result: list[QueryInterval] = []
    ordered = sorted(positions)
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position != previous + 1:
            result.append(QueryInterval(start, previous + 1))
            start = position
        previous = position
    result.append(QueryInterval(start, previous + 1))
    return tuple(result)


def _stable_file_identity(path: Path, label: str) -> _FileIdentity:
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_BUFFER_SIZE), b""):
            digest.update(chunk)
    after = path.stat()
    identity = _stat_tuple(before)
    if identity != _stat_tuple(after):
        raise RuntimeError(f"{label} changed while its digest was computed")
    return _FileIdentity(*identity, digest.hexdigest())


def _assert_file_unchanged(path: Path, expected: _FileIdentity, label: str) -> None:
    observed = _stable_file_identity(path, label)
    if observed != expected:
        raise RuntimeError(f"{label} changed during length measurement")


def _stat_tuple(value: os.stat_result) -> tuple[int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns


def _fasta_contig_length(path: Path, contig: str) -> int:
    lengths: dict[str, int] = {}
    current: str | None = None
    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                name = line[1:].split(maxsplit=1)[0]
                _text(name, "reference FASTA contig")
                if name in lengths:
                    raise ValueError("reference FASTA contains a duplicate contig")
                lengths[name] = 0
                current = name
            elif line:
                if current is None or any(character.isspace() for character in line):
                    raise ValueError("reference FASTA is malformed")
                lengths[current] += len(line)
    length = lengths.get(contig)
    if length is None:
        raise ValueError("reference FASTA is missing the original input contig")
    if length <= 0:
        raise ValueError("reference FASTA contig must contain sequence")
    return length


def _samtools_revision(samtools_path: Path) -> str:
    argv = (str(samtools_path), "--version")
    completed = subprocess.run(argv, capture_output=True, check=False, text=True, env=None, shell=False)
    if completed.returncode != 0:
        raise RuntimeError("samtools version check failed for length measurement")
    lines = completed.stdout.splitlines()
    if len(lines) < 2:
        raise ValueError("samtools version output is truncated")
    samtools_match = _SAMTOOLS_VERSION_PATTERN.fullmatch(lines[0])
    htslib_match = _HTSLIB_VERSION_PATTERN.fullmatch(lines[1])
    if samtools_match is None or htslib_match is None:
        raise ValueError("samtools version output is malformed")
    return f"samtools={samtools_match.group(1)};htslib={htslib_match.group(1)}"


def _validate_fragment_reader(context: LengthFeatureContext) -> None:
    reader = context.fragment_reader
    if reader.name != "pysam" or reader.version != pysam.__version__:
        raise ValueError("pysam version does not match the measurement context")
    htslib_version = getattr(pysam.version, "__htslib_version__", None)
    if reader.htslib_version != htslib_version:
        raise ValueError("pysam htslib version does not match the measurement context")
    if reader.alignment_semantics != "explicit-filtered-aligned-pairs-v1":
        raise ValueError("fragment reader semantics do not match the measurement context")


def _has_cram_magic(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(4) == b"CRAM"


def _read_alignment_evidence(
    input_path: Path,
    reference_path: Path,
    contig: str,
    intervals: tuple[QueryInterval, ...],
    expected_positions: set[int],
    *,
    is_cram: bool,
) -> AlignmentEvidence:
    mutable: dict[int, tuple[int, set[str]]] = {position: (0, set()) for position in expected_positions}
    mode = cast(Literal["rb", "rc"], "rc" if is_cram else "rb")
    try:
        with pysam.AlignmentFile(str(input_path), mode, reference_filename=str(reference_path)) as alignment:
            contig_length, assembly = _header_identity(alignment, contig)
            for interval in intervals:
                for record in alignment.fetch(contig, interval.start, interval.end):
                    if record.flag & _EXCLUDED_ALIGNMENT_FLAGS:
                        continue
                    fragment_id = _fragment_id(record)
                    qualities = record.query_qualities
                    for query_position, reference_position in record.get_aligned_pairs(matches_only=False):
                        if (
                            query_position is None
                            or reference_position is None
                            or reference_position < interval.start
                            or reference_position >= interval.end
                        ):
                            continue
                        if qualities is not None and qualities[query_position] < 0:
                            continue
                        count, fragment_ids = mutable[reference_position]
                        mutable[reference_position] = (count + 1, fragment_ids)
                        fragment_ids.add(fragment_id)
    except (OSError, ValueError) as error:
        raise RuntimeError("failed to read alignment records for length fragment evidence") from error
    return AlignmentEvidence(
        contig_length=contig_length,
        assembly=assembly,
        positions={
            position: FragmentPositionEvidence(count, frozenset(fragment_ids))
            for position, (count, fragment_ids) in mutable.items()
        },
    )


def _header_identity(alignment: pysam.AlignmentFile, contig: str) -> tuple[int, str | None]:
    header = alignment.header.to_dict()
    sequence_records = header.get("SQ")
    if not isinstance(sequence_records, list):
        raise ValueError("alignment header has no sequence dictionary")
    matches = [record for record in sequence_records if isinstance(record, dict) and record.get("SN") == contig]
    if len(matches) != 1:
        raise ValueError("alignment header must contain the original contig exactly once")
    length = matches[0].get("LN")
    if not _is_int(length):
        raise ValueError("alignment header contig length must be a positive integer")
    length_int = cast(int, length)
    if length_int <= 0:
        raise ValueError("alignment header contig length must be a positive integer")
    assembly = matches[0].get("AS")
    if assembly is not None:
        assembly = _text(assembly, "alignment header assembly")
    return length_int, cast(str | None, assembly)


def _fragment_id(record: pysam.AlignedSegment) -> str:
    query_name = _text(record.query_name, "alignment query name")
    read_group: object = record.get_tag("RG") if record.has_tag("RG") else None
    if read_group is not None:
        read_group = _text(read_group, "alignment read group")
    return canonical_sha256({"query_name": query_name, "read_group": read_group})


def _validate_alignment_header(
    evidence: AlignmentEvidence,
    context: LengthFeatureContext,
    required_end: int,
    fasta_length: int,
) -> None:
    if evidence.contig_length != fasta_length:
        raise ValueError("alignment header and reference FASTA contig lengths do not match")
    if evidence.contig_length < required_end:
        raise ValueError("alignment header contig length does not cover the queried interval union")
    if evidence.assembly is not None and evidence.assembly != context.assembly:
        raise ValueError("alignment header assembly does not match the measurement context")


def _write_bed(path: Path, contig: str, intervals: tuple[QueryInterval, ...]) -> None:
    content = "".join(f"{contig}\t{interval.start}\t{interval.end}\n" for interval in intervals)
    path.write_text(content, encoding="utf-8")
    path.chmod(0o600)


def _absolute_path(value: object, label: str) -> Path:
    if not isinstance(value, Path) or not value.is_absolute():
        raise ValueError(f"{label} must be an absolute Path")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a non-empty trimmed string")
    return value


def _is_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)
