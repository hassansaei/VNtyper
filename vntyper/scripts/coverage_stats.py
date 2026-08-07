"""
coverage_stats.py

Pure coverage statistics for the MUC1 VNTR region.

Extracted from ``fastq_bam_processing.calculate_vntr_coverage`` so the arithmetic
can be tested without samtools, a BAM file or a subprocess. The I/O half - running
``samtools depth`` and writing the summary next to the other artefacts - stays in
``fastq_bam_processing``.

**The output schema is a frozen contract.** ``generate_report.py`` reads
``region_length``, ``uncovered_bases`` and ``percent_uncovered`` out of the TSV by
exact, lowercase key and defaults each to ``0`` when the key is absent, and
``tests/helpers.validate_coverage_output`` reads all nine. Renaming a column
therefore raises nothing anywhere - it makes the HTML report state that a sample
with no coverage at all has 0 bp uncovered. :data:`COVERAGE_COLUMNS` is the single
declaration of that schema, and the TSV writer derives from it.

The schema is split in two because a measurement and a verdict are different things.
:data:`COVERAGE_METRIC_COLUMNS` is what :func:`summarise_coverage` returns;
:data:`COVERAGE_COLUMNS` is that plus ``coverage_qc``, the sample-level QC verdict
``fastq_bam_processing.calculate_vntr_coverage`` merges in from
:func:`~vntyper.scripts.coverage_qc.evaluate_coverage_qc` (#172).

One definition is worth stating explicitly because the opposite was true until #171:
every statistic here is over the **region**, not over the covered positions. ``mean``,
``median``, ``stdev``, ``min`` and ``max`` are computed after uncovered positions are
restored as zeros, so a sample covered at 30x across 10% of the VNTR reports
``mean = 3``, ``min = 0`` and ``percent_uncovered = 90``. The old ``mean = 30`` was the
depth where there happened to be reads, which is not a property of the region and was
systematically too high exactly where coverage was patchy.

``samtools depth`` is called with ``-a`` (:func:`command_builders.build_samtools_depth_command`),
so the depth table already spans the whole region and the zeros are real rows.
``uncovered_bases`` counts them rather than deriving them by subtraction, which is what
makes the two changes one change: under ``-a``, ``region_length - len(values)`` is 0 for
every sample.

Functions:
    parse_region_length: Region string to inclusive base count
    read_depth_values: ``samtools depth`` output file to a list of depths
    summarise_coverage: Depth values plus region length to the frozen schema
    format_coverage_summary: The frozen schema to the exact TSV bytes
"""

from __future__ import annotations

import logging
import statistics
from pathlib import Path

logger = logging.getLogger(__name__)

#: The eight measured statistics. Exactly what :func:`summarise_coverage` returns.
COVERAGE_METRIC_COLUMNS: tuple[str, ...] = (
    "mean",
    "median",
    "stdev",
    "min",
    "max",
    "region_length",
    "uncovered_bases",
    "percent_uncovered",
)

#: The coverage summary schema, in the order the TSV writes it: the measurements plus
#: the QC verdict.
#:
#: Frozen contract (wave-2 C1). ``generate_report.py`` looks these up with
#: ``.get(name, 0)``, so a rename is silent: the report loses coverage rather than
#: failing. Changing this tuple changes the report's input format.
#:
#: The verdict is not a measurement and :func:`summarise_coverage` cannot produce it -
#: it has no thresholds - so ``calculate_vntr_coverage`` fills it in before writing.
#: :func:`format_coverage_summary` still raises ``KeyError`` when it is absent, which is
#: the contract that keeps a caller from writing a summary with no verdict in it (#172).
COVERAGE_COLUMNS: tuple[str, ...] = COVERAGE_METRIC_COLUMNS + ("coverage_qc",)

#: Columns rendered with two decimal places; the rest are written as-is.
_TWO_DECIMAL_COLUMNS = frozenset({"mean", "median", "stdev", "percent_uncovered"})


def parse_region_length(region: str) -> int:
    """
    Compute the inclusive length in bases of a ``contig:start-end`` region string.

    Args:
        region (str): Region string, e.g. ``chr1:155160500-155162000``. Any contig
            naming works - UCSC, Ensembl or an NCBI accession - because only the
            part after the last colon-separated coordinate block is parsed.

    Returns:
        int: ``end - start + 1``, or ``0`` if the string cannot be parsed.

    Warning:
        A malformed region degrades to ``0`` rather than raising. That is the
        pre-existing behaviour and it is preserved, but it means
        ``percent_uncovered`` is then reported as ``0`` for a sample whose region
        could not be read at all. The WARNING log is the only signal.

    Examples:
        >>> parse_region_length("chr1:155160500-155162000")
        1501
        >>> parse_region_length("chr1:100-100")
        1
        >>> parse_region_length("chr1")
        0
    """
    try:
        region_parts = region.split(":")
        if len(region_parts) != 2:
            raise ValueError(f"Invalid region format: {region}")

        pos_range = region_parts[1].split("-")
        if len(pos_range) != 2:
            raise ValueError(f"Invalid position range: {region_parts[1]}")

        start_pos = int(pos_range[0])
        end_pos = int(pos_range[1])
        total_region_length = end_pos - start_pos + 1
        logger.debug(f"VNTR region total length: {total_region_length} bp")
        return total_region_length
    except (ValueError, IndexError) as e:
        logger.warning(f"Could not parse region string: {e}. Setting region length to 0.")
        return 0


def read_depth_values(depth_file: str | Path) -> list[int]:
    """
    Read the depth column out of a ``samtools depth`` output file.

    ``samtools depth -a`` writes one tab-separated row per position in the region:
    contig, 1-based position, depth - zero-depth positions included. Without ``-a``
    the uncovered positions are absent entirely, and
    :func:`summarise_coverage` restores them (#171).

    Args:
        depth_file (str | Path): Path to the depth table.

    Returns:
        list[int]: The third column of every non-blank line, in file order.

    Raises:
        ValueError: If a line has fewer than three fields or a non-integer depth.
        OSError: If the file cannot be read.
    """
    with open(depth_file) as f:
        return [int(line.strip().split("\t")[2]) for line in f if line.strip()]


def summarise_coverage(coverage_values: list[int], total_region_length: int) -> dict:
    """
    Summarise per-base depths into the frozen coverage schema.

    Args:
        coverage_values (list[int]): Depth at each position the depth table
            carries, as returned by :func:`read_depth_values`. Under ``-a`` that
            is every position in the region, zeros included.
        total_region_length (int): Region length in bases, as returned by
            :func:`parse_region_length`. ``0`` means the region could not be
            parsed.

    Returns:
        dict: Exactly the keys in :data:`COVERAGE_METRIC_COLUMNS` - the measurements,
        without ``coverage_qc``, which is the caller's verdict. Value types are
        preserved from the pre-extraction implementation and are deliberately
        mixed: ``median`` is an ``int`` for an odd-length integer sample and a
        ``float`` for an even-length one, ``stdev`` is the integer ``0`` for a
        single observation, and ``min``/``max``/``region_length``/
        ``uncovered_bases`` are ``int``.

    Raises:
        RuntimeError: If ``coverage_values`` is empty. Under ``-a`` an empty depth
            table can only mean the region matched no contig at all - a wrong
            contig name, or a region past the end of the chromosome - because a
            covered-free region would still emit its zero rows. Returning
            ``mean = 0`` instead would let the pipeline carry on and report a
            confident negative genotype.

    Note:
        Every statistic is over the **region**: uncovered positions are restored as
        zeros first, so ``mean`` is the region's mean depth and ``min`` is ``0``
        wherever anything is uncovered (#171). ``uncovered_bases`` counts zero-depth
        positions rather than subtracting the row count from the region length.
    """
    if not coverage_values:
        raise RuntimeError("No coverage data found.")

    # `samtools depth -a` emits one row per position in the region, zeros included, so
    # `coverage_values` *is* the region. A file written without `-a` - a legacy artefact,
    # or a region truncated at a contig end - is short by exactly the positions that had
    # no reads, so restoring them as zeros makes both cases the same base set. `absent` is
    # floored at 0: an unparseable region (length 0) must not subtract padding. See #171.
    absent = max(total_region_length - len(coverage_values), 0)
    base = coverage_values + [0] * absent

    uncovered_bases = sum(1 for depth in base if depth == 0)
    percent_uncovered = 0 if total_region_length <= 0 else uncovered_bases / total_region_length * 100

    return {
        "mean": sum(base) / len(base),
        "median": statistics.median(base),
        "stdev": statistics.stdev(base) if len(base) > 1 else 0,
        "min": min(base),
        "max": max(base),
        "region_length": total_region_length,
        "uncovered_bases": uncovered_bases,
        "percent_uncovered": percent_uncovered,
    }


def format_coverage_summary(stats: dict) -> str:
    """
    Render a coverage summary as the two-line TSV the report reads.

    Args:
        stats (dict): A summary as returned by :func:`summarise_coverage`.

    Returns:
        str: A header line and one data line, both tab-separated and
        newline-terminated. ``mean``, ``median``, ``stdev`` and
        ``percent_uncovered`` are written to two decimal places; the counts are
        written as integers.

    Raises:
        KeyError: If ``stats`` is missing any of :data:`COVERAGE_COLUMNS`.
    """
    header = "\t".join(COVERAGE_COLUMNS)
    values = "\t".join(
        f"{stats[column]:.2f}" if column in _TWO_DECIMAL_COLUMNS else f"{stats[column]}" for column in COVERAGE_COLUMNS
    )
    return f"{header}\n{values}\n"
