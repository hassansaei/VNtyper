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
so the depth table normally spans the whole region and the zeros are real rows.
``uncovered_bases`` counts those zero rows rather than deriving them by subtraction, which is
what makes the two changes one change: under ``-a``, ``region_length - len(values)`` is 0 for
every sample.

Normally, not always. A region truncated at a contig end legitimately yields a *short* file,
which is why :func:`summarise_coverage` pads a short input rather than rejecting it. A file
*longer* than its own region is incoherent under either convention, and is refused (#224).

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

#: Largest region span :func:`parse_region_length` accepts, in bases.
#:
#: Not a configuration key, deliberately. ``--config-path`` replaces the whole config rather
#: than merging it (AGENTS.md trap 2), so a new key can go missing - and a sanity bound that
#: aborts a run when its own key is absent is a worse failure than the one it guards. The MUC1
#: VNTR regions are 1501 and 4501 bp, so this carries three orders of magnitude of headroom.
#:
#: The bound exists because ``samtools depth -a`` writes one row per *declared* base, which
#: makes the depth **file** unbounded as well as the list built from it: ``chr1:1-250000000``
#: emitted roughly 48 million rows and exhausted 30 GB of disk before Python allocated
#: anything. Disk is the first failure mode here, not memory - Python's own share is about 16
#: bytes per declared base, or 3.8 GiB for that region. That is why
#: ``fastq_bam_processing.calculate_vntr_coverage`` parses the region *before* it runs
#: samtools rather than after (#224).
MAX_REGION_SPAN_BASES: int = 10_000_000


def parse_region_length(region: str) -> int:
    """
    Compute the inclusive length in bases of a ``contig:start-end`` region string.

    Args:
        region (str): Region string, e.g. ``chr1:155160500-155162000``. Any contig
            naming works - UCSC, Ensembl or an NCBI accession - because only the
            part after the last colon-separated coordinate block is parsed.

    Returns:
        int: ``end - start + 1``, or ``0`` if the string cannot be parsed.

    Raises:
        ValueError: If the region parses but its coordinates are impossible - a start below
            1, an end before the start, or a span above :data:`MAX_REGION_SPAN_BASES`.

    Warning:
        A region that cannot be *parsed at all* still degrades to ``0`` rather than raising.
        That is the pre-existing behaviour and it is preserved, but it means
        ``percent_uncovered`` is then reported as ``0`` for a sample whose region could not be
        read at all. The WARNING log is the only signal.

        That is deliberately distinct from the ``ValueError`` above: "could not be parsed" and
        "parsed to something impossible" are different states, and only the second is a value
        somebody configured (#224).

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
    except (ValueError, IndexError) as e:
        logger.warning(f"Could not parse region string: {e}. Setting region length to 0.")
        return 0

    # This validation sits OUTSIDE the `try` above on purpose. That block catches
    # `(ValueError, IndexError)` and returns 0, so a raise placed inside it would be caught
    # by its own handler and degraded to a silent zero-length region -- which is the failure
    # mode being fixed, not a fallback. "Could not be parsed" and "parsed to something
    # impossible" are different states and must stay different (#224).
    if start_pos < 1:
        msg = f"Invalid region {region!r}: start {start_pos} is below 1. Region coordinates are 1-based."
        logger.error(msg)
        raise ValueError(msg)

    if end_pos < start_pos:
        msg = (
            f"Invalid region {region!r}: end {end_pos} is before start {start_pos}. A reversed region "
            "yields a negative length, which is then reported beside a passing coverage verdict."
        )
        logger.error(msg)
        raise ValueError(msg)

    total_region_length = end_pos - start_pos + 1
    if total_region_length > MAX_REGION_SPAN_BASES:
        msg = (
            f"Configured region {region!r} spans {total_region_length} bases, above the "
            f"{MAX_REGION_SPAN_BASES} limit. `samtools depth -a` writes one row per base in the "
            "region, so this exhausts disk before any coverage summary is produced. The MUC1 VNTR "
            "regions are 1501 and 4501 bp."
        )
        logger.error(msg)
        raise ValueError(msg)

    logger.debug(f"VNTR region total length: {total_region_length} bp")
    return total_region_length


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
        ValueError: If ``coverage_values`` is longer than a positive
            ``total_region_length``. The depth table and the region string then disagree, and
            every statistic would be computed over more bases than the region it describes -
            ``percent_uncovered`` alone can exceed 100 (#224).
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

    # `absent` below pads a SHORT list up to the region. A LONG one has no correct summary:
    # every statistic would be computed over a base set larger than the region it claims to
    # describe, and `percent_uncovered` alone can exceed 100 -- `summarise_coverage([0]*4, 2)`
    # reported 200.0. Clamping would invent a number instead of refusing one.
    #
    # Guarded on a positive length because 0 does not mean "a zero-base region", it means the
    # region string could not be parsed at all, and an unknown region cannot be compared
    # against a row count. `summarise_coverage([10, 20], 0)` is a pinned contract (#224).
    if total_region_length > 0 and len(coverage_values) > total_region_length:
        msg = (
            f"The depth table carries {len(coverage_values)} rows for a region of {total_region_length} "
            "bases. There is no correct summary to return: every statistic would cover more bases than "
            "the region it describes, and percent_uncovered would exceed 100."
        )
        logger.error(msg)
        raise ValueError(msg)

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
