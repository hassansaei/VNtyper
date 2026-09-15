"""Small development-cohort TSV contract with independent optional truth fields."""

from __future__ import annotations

import csv
import hashlib
import io
import math
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.calibration_secure_io import read_regular_path

_REQUIRED = {"sample_id", "bam", "assembly"}
_OPTIONAL = {"genotype", "allele_1", "allele_2", "group_id", "kestrel_result", "advntr_result"}


@dataclass(frozen=True)
class CohortSample:
    """One declared biological sample; group metadata is never a predictor."""

    sample_id: str
    bam: Path
    assembly: str
    genotype: bool | None
    length_total: float | None
    group_id: str
    kestrel_result: Path | None = None
    advntr_result: Path | None = None
    length_reason: str | None = None


def _path(raw: str, parent: Path) -> Path | None:
    return (parent / raw).resolve() if raw else None


def _length(first: str, second: str) -> float | None:
    if not first and not second:
        return None
    values = tuple(float(value) for value in (first, second) if value)
    if any(not math.isfinite(value) or value <= 0 or not value.is_integer() for value in values):
        raise ValueError("cohort allele counts must be finite positive integers")
    if len(values) != 2:
        return None
    total = sum(values)
    if not math.isfinite(total):
        raise ValueError("cohort total length is not finite")
    return total


def read_cohort_manifest(path: Path) -> tuple[CohortSample, ...]:
    """Read strict TSV rows without opening alignment or caller outcomes.

    Args:
        path: Local TSV with required sample_id, bam and assembly columns.

    Returns:
        Sample-sorted immutable rows, preserving missing targets independently.

    Raises:
        ValueError: If columns, truth, identities, paths or row geometry differ.
    """
    parent = path.resolve().parent
    reader = csv.DictReader(io.StringIO(read_regular_path(path).decode("utf-8-sig")), delimiter="\t")
    fields = reader.fieldnames
    if (
        not fields
        or len(fields) != len(set(fields))
        or not set(fields) >= _REQUIRED
        or set(fields) - (_REQUIRED | _OPTIONAL)
    ):
        raise ValueError("cohort manifest columns must contain required and supported unique fields")
    rows = []
    identities: set[str] = set()
    paths: set[Path] = set()
    inodes: set[tuple[int, int]] = set()
    for raw in reader:
        if None in raw or any(value is None or value != value.strip() for value in raw.values()):
            raise ValueError("cohort manifest rows must match the header and contain trimmed values")
        if any(not raw[key] for key in _REQUIRED):
            raise ValueError("cohort required values must not be empty")
        identity = raw["sample_id"]
        bam = _path(raw["bam"], parent)
        assert bam is not None
        if identity in identities or bam in paths:
            raise ValueError("cohort duplicate sample identity or alignment contribution")
        if bam.exists():
            stat = bam.stat()
            inode = (stat.st_dev, stat.st_ino)
            if inode in inodes:
                raise ValueError("cohort duplicate alignment inode contribution")
            inodes.add(inode)
        genotype = raw.get("genotype", "")
        if genotype not in {"", "unknown", "positive", "negative"}:
            raise ValueError("cohort genotype must be positive, negative, unknown, or empty")
        identities.add(identity)
        paths.add(bam)
        group = "declared:" + raw["group_id"] if raw.get("group_id") else "sample:" + identity
        rows.append(
            CohortSample(
                identity,
                bam,
                raw["assembly"],
                None if genotype in {"", "unknown"} else genotype == "positive",
                _length(raw.get("allele_1", ""), raw.get("allele_2", "")),
                group,
                _path(raw.get("kestrel_result", ""), parent),
                _path(raw.get("advntr_result", ""), parent),
                "incomplete_allele_pair" if bool(raw.get("allele_1")) != bool(raw.get("allele_2")) else None,
            )
        )
    if not rows:
        raise ValueError("cohort manifest must have at least one sample")
    return tuple(sorted(rows, key=lambda row: row.sample_id))


def audit_alignment_duplicates(samples: tuple[CohortSample, ...]) -> dict[str, str | None]:
    """Hash available alignment bytes and refuse repeated physical-file contributions.

    Args:
        samples: Already validated cohort rows.

    Returns:
        Observed content digests, with null for unavailable alignments. Distinct bytes
        do not establish distinct biological identity or semantic readset independence.

    Raises:
        ValueError: For duplicate content or a file changed during hashing.
    """
    observed: dict[str, str | None] = {}
    seen: set[str] = set()
    for sample in samples:
        if not sample.bam.is_file():
            observed[sample.sample_id] = None
            continue
        before = sample.bam.stat()
        digest = hashlib.sha256()
        with sample.bam.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        after = sample.bam.stat()
        if any(
            getattr(before, name) != getattr(after, name)
            for name in ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        ):
            raise ValueError("cohort alignment changed during duplicate audit")
        value = digest.hexdigest()
        if value in seen:
            raise ValueError("cohort duplicate alignment content contribution")
        seen.add(value)
        observed[sample.sample_id] = value
    return observed
