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
``tests/helpers.validate_coverage_output`` reads all eight. Renaming a column
therefore raises nothing anywhere - it makes the HTML report state that a sample
with no coverage at all has 0 bp uncovered. :data:`COVERAGE_COLUMNS` is the single
declaration of that schema and both the TSV writer and the returned dict derive
from it.

One definition is worth stating explicitly because it is easy to get wrong when
reading the numbers: ``samtools depth`` emits a row **only** for positions covered
by at least one read. The number of rows is therefore the covered-base count, and
``mean`` is the mean depth over *covered* positions - not over the region. A
sample covered at 30x across 10% of the VNTR reports ``mean = 30`` and
``percent_uncovered = 90``, and it is the second number that tells you the call is
unreliable.

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

#: The coverage summary schema, in the order the TSV writes it.
#:
#: Frozen contract (wave-2 C1). ``generate_report.py`` looks these up with
#: ``.get(name, 0)``, so a rename is silent: the report loses coverage rather than
#: failing. Changing this tuple changes the report's input format.
COVERAGE_COLUMNS: tuple[str, ...] = (
    "mean",
    "median",
    "stdev",
    "min",
    "max",
    "region_length",
    "uncovered_bases",
    "percent_uncovered",
)

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

    ``samtools depth`` writes one tab-separated row per **covered** position:
    contig, 1-based position, depth. Positions with no reads are absent entirely,
    which is why the row count is the covered-base count.

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
        coverage_values (list[int]): Depth at each covered position, as returned
            by :func:`read_depth_values`.
        total_region_length (int): Region length in bases, as returned by
            :func:`parse_region_length`. ``0`` means the region could not be
            parsed.

    Returns:
        dict: Exactly the keys in :data:`COVERAGE_COLUMNS`. Value types are
        preserved from the pre-extraction implementation and are deliberately
        mixed: ``median`` is an ``int`` for an odd-length integer sample and a
        ``float`` for an even-length one, ``stdev`` is the integer ``0`` for a
        single observation, and ``min``/``max``/``region_length``/
        ``uncovered_bases`` are ``int``.

    Raises:
        RuntimeError: If ``coverage_values`` is empty. An empty depth table means
            the region matched nothing - a wrong contig name, a region past the
            end of the chromosome, or a BAM with no reads there. Returning
            ``mean = 0`` instead would let the pipeline carry on and report a
            confident negative genotype.

    Note:
        ``mean`` is the mean over covered positions, not over the region. Read it
        together with ``percent_uncovered``.
    """
    if not coverage_values:
        raise RuntimeError("No coverage data found.")

    covered_bases_count = len(coverage_values)
    zero_coverage_bases = total_region_length - covered_bases_count

    # A region length of 0 means the region string could not be parsed; there is
    # no denominator to divide by, so report 0 rather than crashing.
    percent_uncovered = 0 if total_region_length <= 0 else zero_coverage_bases / total_region_length * 100

    return {
        "mean": sum(coverage_values) / len(coverage_values),
        "median": statistics.median(coverage_values),
        "stdev": statistics.stdev(coverage_values) if len(coverage_values) > 1 else 0,
        "min": min(coverage_values),
        "max": max(coverage_values),
        "region_length": total_region_length,
        "uncovered_bases": zero_coverage_bases,
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
