"""Derive a single-end BAM from a declared paired-end cohort fixture."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SingleEndFixture:
    """A validated, portable single-end fixture declaration."""

    name: str
    source_bam: Path
    output_bam: Path


def parse_single_end_fixture(
    entry: Mapping[str, object],
    *,
    root: Path = Path("."),
) -> SingleEndFixture:
    """Validate one ``single_end_bam`` fixture declaration.

    Args:
        entry: Mapping from ``tests/test_data_config.json``.
        root: Repository root used to resolve the declaration's relative paths.

    Returns:
        The validated fixture name and resolved paths.

    Raises:
        ValueError: If the declaration is incomplete, has the wrong kind, uses
            absolute paths, or would overwrite its own source.
    """
    name = entry.get("name")
    kind = entry.get("kind")
    source_value = entry.get("source_bam")
    output_value = entry.get("output_bam")
    values = (name, source_value, output_value)
    if kind != "single_end_bam" or any(not isinstance(value, str) or not value.strip() for value in values):
        raise ValueError(
            "Invalid single-end fixture declaration: kind must be single_end_bam and name/source_bam/output_bam "
            "must be non-empty strings."
        )

    assert isinstance(name, str)
    assert isinstance(source_value, str)
    assert isinstance(output_value, str)
    source_relative = Path(source_value)
    output_relative = Path(output_value)
    paths = (source_relative, output_relative)
    if Path(name).name != name:
        raise ValueError("Invalid single-end fixture declaration: name must be a portable basename.")
    if any(path.is_absolute() or ".." in path.parts or path.suffix != ".bam" for path in paths):
        raise ValueError(
            "Invalid single-end fixture declaration: source_bam and output_bam must be relative BAM paths "
            "contained by the repository."
        )

    source_bam = root / source_relative
    output_bam = root / output_relative
    if source_bam == output_bam:
        raise ValueError("Invalid single-end fixture declaration: output_bam must differ from source_bam.")
    return SingleEndFixture(name=name, source_bam=source_bam, output_bam=output_bam)


def derive_single_end_bam(spec: SingleEndFixture) -> int:
    """Copy every alignment while clearing all paired-end relationship flags.

    The read name, sequence, qualities, alignment and every unrelated flag stay
    unchanged. The output is indexed after the complete BAM has been written.

    Args:
        spec: Validated fixture declaration.

    Returns:
        The number of records copied to the single-end BAM.

    Raises:
        FileNotFoundError: If the declared source BAM is absent.
        ValueError: If source and destination resolve to the same file.
    """
    if not spec.source_bam.is_file():
        raise FileNotFoundError(f"Single-end fixture source BAM does not exist: {spec.source_bam}")
    if spec.source_bam.resolve() == spec.output_bam.resolve():
        raise ValueError("Single-end fixture output BAM must differ from its source BAM.")

    spec.output_bam.parent.mkdir(parents=True, exist_ok=True)
    records = 0
    with (
        pysam.AlignmentFile(str(spec.source_bam), "rb") as source,
        pysam.AlignmentFile(str(spec.output_bam), "wb", template=source) as output,
    ):
        for read in source.fetch(until_eof=True):
            read.is_paired = False
            read.is_read1 = False
            read.is_read2 = False
            read.is_proper_pair = False
            read.mate_is_unmapped = False
            output.write(read)
            records += 1
    pysam.samtools.index(str(spec.output_bam))
    logger.info("Derived %s with %d single-end records at %s", spec.name, records, spec.output_bam)
    return records
