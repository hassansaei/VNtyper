"""Efficient indexed locus reader for the standard 13-feature length model."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Literal, cast

import pysam

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_features import DepthPosition
from vntyper.scripts.length_standard_features import (
    STANDARD_ACCEPTED_CONTIGS,
    STANDARD_ASSEMBLY,
    STANDARD_LOCUS_END,
    STANDARD_LOCUS_START,
    STANDARD_REFERENCE_LOCUS_SHA256,
    StandardLengthMeasurement,
    StandardReadSummary,
    extract_standard_length_features,
)
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution, restore_reference_resolution
from vntyper.scripts.region_utils import resolve_assembly_alias

_EXCLUDED_FLAGS = 0x704


def _path(path: object, label: str) -> Path:
    if not isinstance(path, Path) or not path.is_absolute() or not path.is_file():
        raise ValueError(f"standard length {label} must be an existing absolute file Path")
    return path


def _cram(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(4) == b"CRAM"


def _contig(alignment: pysam.AlignmentFile) -> str:
    matches = tuple(contig for contig in alignment.references if contig in STANDARD_ACCEPTED_CONTIGS)
    if len(matches) != 1:
        raise ValueError("standard length alignment must contain exactly one accepted chr1 contig")
    index = alignment.references.index(matches[0])
    if alignment.lengths[index] < STANDARD_LOCUS_END:
        raise ValueError("standard length alignment contig does not cover the locus")
    header = alignment.header.to_dict()
    records = header.get("SQ", [])
    sq = [record for record in records if isinstance(record, dict) and record.get("SN") == matches[0]]
    if len(sq) != 1:
        raise ValueError("standard length alignment sequence header differs")
    assembly = sq[0].get("AS")
    if assembly is not None and assembly != STANDARD_ASSEMBLY:
        raise ValueError("standard length alignment header assembly differs")
    return matches[0]


def _reference_locus(reference_path: Path, contig: str) -> str:
    try:
        with pysam.FastaFile(str(reference_path)) as reference:
            candidates = tuple(name for name in STANDARD_ACCEPTED_CONTIGS if name in reference.references)
            reference_contig = contig if contig in candidates else candidates[0] if len(candidates) == 1 else None
            if reference_contig is None:
                raise ValueError("standard length reference must contain exactly one compatible chr1 contig")
            sequence = reference.fetch(reference_contig, STANDARD_LOCUS_START, STANDARD_LOCUS_END).upper()
    except (OSError, ValueError) as error:
        raise ValueError("standard length reference locus cannot be read through its FASTA index") from error
    if len(sequence) != STANDARD_LOCUS_END - STANDARD_LOCUS_START:
        raise ValueError("standard length reference locus is truncated")
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if digest != STANDARD_REFERENCE_LOCUS_SHA256:
        raise ValueError("standard length reference locus digest differs")
    return digest


def _fragment_id(record: pysam.AlignedSegment) -> str:
    query_name = record.query_name
    if not isinstance(query_name, str) or not query_name:
        raise ValueError("standard length alignment query name is missing")
    read_group = record.get_tag("RG") if record.has_tag("RG") else None
    if read_group is not None and (not isinstance(read_group, str) or not read_group):
        raise ValueError("standard length alignment read group is invalid")
    return canonical_sha256({"query_name": query_name, "read_group": read_group})


def _read(
    input_path: Path,
    reference_path: Path,
    index_path: Path | None,
    *,
    is_cram: bool,
) -> tuple[str, tuple[DepthPosition, ...], StandardReadSummary]:
    counts = [0] * (STANDARD_LOCUS_END - STANDARD_LOCUS_START)
    fragments: list[set[str]] = [set() for _ in counts]
    eligible = mapq_zero = mapq_sum = clipped = gc_bases = query_bases = 0
    mode = cast(Literal["rb", "rc"], "rc" if is_cram else "rb")
    try:
        with pysam.AlignmentFile(
            str(input_path),
            mode,
            reference_filename=str(reference_path),
            index_filename=None if index_path is None else str(index_path),
        ) as alignment:
            contig = _contig(alignment)
            for record in alignment.fetch(contig, STANDARD_LOCUS_START, STANDARD_LOCUS_END):
                if record.flag & _EXCLUDED_FLAGS:
                    continue
                sequence = record.query_sequence
                if sequence:
                    eligible += 1
                    mapq_zero += record.mapping_quality == 0
                    mapq_sum += record.mapping_quality
                    clipped += any(operation == 4 for operation, _ in (record.cigartuples or ()))
                    upper = sequence.upper()
                    gc_bases += upper.count("G") + upper.count("C")
                    query_bases += len(upper)
                fragment = _fragment_id(record)
                qualities = record.query_qualities
                for query_position, reference_position in record.get_aligned_pairs(matches_only=False):
                    if (
                        query_position is None
                        or reference_position is None
                        or reference_position < STANDARD_LOCUS_START
                        or reference_position >= STANDARD_LOCUS_END
                    ):
                        continue
                    if qualities is not None and qualities[query_position] < 0:
                        continue
                    offset = reference_position - STANDARD_LOCUS_START
                    counts[offset] += 1
                    fragments[offset].add(fragment)
    except (OSError, ValueError) as error:
        raise RuntimeError("standard length indexed alignment locus cannot be read") from error
    depths = tuple(
        DepthPosition(
            contig,
            STANDARD_LOCUS_START + offset,
            count,
            tuple(sorted(fragments[offset])),
        )
        for offset, count in enumerate(counts)
    )
    return contig, depths, StandardReadSummary(eligible, mapq_zero, mapq_sum, clipped, gc_bases, query_bases)


def read_standard_length_features(
    input_path: Path,
    reference_path: Path,
    *,
    assembly: str,
    index_path: Path | None = None,
) -> StandardLengthMeasurement:
    """Read the GRCh38 locus once and derive the standard 13-feature measurement.

    This reader intentionally hashes only the public reference locus. Callers must
    keep any descriptor-backed alignment plan alive for this complete call.
    """
    input_path = _path(input_path, "input")
    reference_path = _path(reference_path, "reference")
    if not isinstance(assembly, str) or resolve_assembly_alias(assembly) != STANDARD_ASSEMBLY:
        raise ValueError("standard length reader supports GRCh38 only")
    if index_path is not None:
        index_path = _path(index_path, "index")
    before = os.stat(input_path)
    reference_before = os.stat(reference_path)
    is_cram = _cram(input_path)
    previous_ref_path: str | None = None
    if is_cram:
        previous_ref_path = pin_reference_resolution(
            {"cram": {"allow_ambient_reference_resolution": False, "local_ref_path": str(reference_path)}}
        )
    try:
        # Validate the exact public reference sequence before opening alignment records.
        reference_locus_sha256 = _reference_locus(reference_path, STANDARD_ACCEPTED_CONTIGS[0])
        contig, depths, read_summary = _read(
            input_path,
            reference_path,
            index_path,
            is_cram=is_cram,
        )
    finally:
        if is_cram:
            restore_reference_resolution(previous_ref_path)
    after = os.stat(input_path)
    reference_after = os.stat(reference_path)
    before_identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns)
    after_identity = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns)
    if before_identity != after_identity:
        raise RuntimeError("standard length input changed during indexed locus measurement")
    reference_before_identity = (
        reference_before.st_dev,
        reference_before.st_ino,
        reference_before.st_size,
        reference_before.st_mtime_ns,
        reference_before.st_ctime_ns,
    )
    reference_after_identity = (
        reference_after.st_dev,
        reference_after.st_ino,
        reference_after.st_size,
        reference_after.st_mtime_ns,
        reference_after.st_ctime_ns,
    )
    if reference_before_identity != reference_after_identity:
        raise RuntimeError("standard length reference changed during indexed locus measurement")
    return extract_standard_length_features(
        depths,
        read_summary,
        assembly=STANDARD_ASSEMBLY,
        contig=contig,
        reference_locus_sha256=reference_locus_sha256,
    )
