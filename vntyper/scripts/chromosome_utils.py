"""
chromosome_utils.py

Chromosome name detection and resolution utilities.

This module provides functions to detect chromosome naming conventions
from BAM/CRAM headers and resolve chromosome identifiers to support
multiple naming schemes (UCSC, simple numeric, NCBI accessions).

Functions:
    detect_naming_convention: Identify the naming convention used in a BAM file
    get_chromosome_name_from_bam: Get the actual chromosome name for a given number
    validate_chromosome_name: Validate that a chromosome name follows expected patterns
"""

from __future__ import annotations

import logging
import math
import re

logger = logging.getLogger(__name__)

# Chr1 length constants for assembly detection (Phase 1 - Issue #139 fix)
CHR1_LENGTHS = {
    "GRCh37": 249250621,
    "hg19": 249250621,  # Alias for GRCh37
    "GRCh38": 248956422,
    "hg38": 248956422,  # Alias for GRCh38
}

# Chr1 under every naming convention this module already handles elsewhere:
# UCSC ("chr1"), ENSEMBL ("1") and NCBI RefSeq ("NC_000001.10" / ".11").
# The accession version is deliberately not read as evidence of the build --
# the length is the evidence, and a mislabelled accession must not override it.
CHR1_NAMES = frozenset({"chr1", "1"})
NCBI_CHR1_PATTERN = re.compile(r"^nc_0*1(\.\d+)?$")
NAMING_CONVENTION_ERROR_PREFIX = "Ambiguous or unclassifiable chromosome naming convention"
NAMING_POLICY_ERROR_PREFIXES = ("Invalid naming_convention_threshold", "Invalid primary_contig_patterns")


def is_chr1_name(contig_name: str) -> bool:
    """Return whether a contig name denotes chromosome 1.

    Args:
        contig_name (str): Contig name from a BAM/CRAM header.

    Returns:
        bool: True for "chr1", "1" and any ``NC_000001.<version>`` accession,
        in any case; False otherwise.

    Examples:
        >>> is_chr1_name("chr1")
        True
        >>> is_chr1_name("NC_000001.11")
        True
        >>> is_chr1_name("NC_000012.11")
        False
    """
    lowered = contig_name.lower()
    return lowered in CHR1_NAMES or bool(NCBI_CHR1_PATTERN.match(lowered))


def detect_assembly_from_chr1_length(contigs: list[dict]) -> str | None:
    """
    Detect reference assembly from chr1 length - most reliable method.

    This function addresses Issue #139 where assembly detection fails on modern
    reference genomes (hg38) with alternate contigs. Chr1 length is unique per
    assembly and unaffected by the presence of ~170 alternate sequences.

    Args:
        contigs (list[dict]): List of contig dictionaries with 'name' and 'length' keys.
                             Example: [{'name': 'chr1', 'length': 248956422}, ...]

    Returns:
        str | None: Assembly name ('GRCh38', 'hg38', 'GRCh37', 'hg19') or None if:
                   - Chr1 not found in contigs
                   - Chr1 length doesn't match known assemblies
                   - Chr1 length is missing or not an integer
                   - Contigs list is empty

    Examples:
        >>> contigs = [{"name": "chr1", "length": 248956422}]
        >>> detect_assembly_from_chr1_length(contigs)
        'GRCh38'

        >>> contigs = [{"name": "1", "length": 249250621}]  # ENSEMBL naming
        >>> detect_assembly_from_chr1_length(contigs)
        'GRCh37'

    Notes:
        - Handles UCSC ('chr1'), ENSEMBL ('1') and NCBI ('NC_000001.10') naming
        - Returns first matching assembly name (GRCh* names prioritized over hg* aliases)
        - Most reliable detection method (proven in frontend, 95%+ success rate)
        - Unaffected by alternate contigs (chr1_random, chr1_alt, etc.)

    See Also:
        - GitHub Issue #139: Assembly detection fails on hg38 with 195 contigs
        - plan/02_active/2025-10-27_assembly-detection-fix/SENIOR_DEVELOPER_REVIEW.md
    """
    if not contigs:
        logger.debug("Empty contigs list provided to detect_assembly_from_chr1_length")
        return None

    # Find chr1 under any supported naming convention, case-insensitive
    chr1 = next((c for c in contigs if is_chr1_name(c.get("name", ""))), None)

    if not chr1:
        logger.debug("Chr1 not found in contigs (tried 'chr1', '1' and NC_000001.*)")
        return None

    chr1_length = chr1.get("length")
    if chr1_length is None:
        logger.warning("Chr1 found but length is None")
        return None

    if not isinstance(chr1_length, int):
        # A float would compare equal to the integer constants and a string
        # would crash the formatted log line below; neither is a length.
        logger.warning(f"Chr1 found but length is not an integer: {chr1_length!r}")
        return None

    logger.debug(f"Chr1 length: {chr1_length:,} bp")

    # Check against known assemblies (GRCh* prioritized over hg* aliases)
    for assembly in ["GRCh38", "hg38", "GRCh37", "hg19"]:
        if chr1_length == CHR1_LENGTHS[assembly]:
            logger.info(f"Assembly detected from chr1 length: {assembly} ({chr1_length:,} bp)")
            return assembly

    logger.warning(
        f"Unknown chr1 length: {chr1_length:,} bp. "
        f"Expected GRCh38/hg38: {CHR1_LENGTHS['GRCh38']:,} bp or "
        f"GRCh37/hg19: {CHR1_LENGTHS['GRCh37']:,} bp"
    )
    return None


def detect_naming_convention(contig_names: list[str], config: dict | None = None) -> str:
    """
    Detect the chromosome naming convention from a list of contig names.

    This function examines contig names to determine if they follow:
    - UCSC convention: chr1, chr2, ..., chrX, chrY, chrM
    - ENSEMBL simple numeric: 1, 2, ..., X, Y, MT
    - NCBI accessions: NC_000001.XX, NC_000002.XX, ...

    Args:
        contig_names (list[str]): List of contig names from BAM header.
        config (dict | None): Pipeline configuration. ``assembly_detection`` can
            set ``naming_convention_threshold`` and
            ``primary_contig_patterns``. Missing values use the shipped defaults.

    Returns:
        str: Convention identifier ("ucsc", "ensembl", "ncbi", "unknown")

    Raises:
        ValueError: If the configured naming-convention threshold is not a
            finite number from 0 through 1, or if the configured primary
            patterns are not three nonempty, compilable strings in UCSC, NCBI,
            ENSEMBL order.

    Examples:
        >>> detect_naming_convention(["chr1", "chr2", "chrX"])
        'ucsc'
        >>> detect_naming_convention(["1", "2", "X"])
        'ensembl'
        >>> detect_naming_convention(["NC_000001.10", "NC_000002.11"])
        'ncbi'
    """
    if not contig_names:
        logger.warning("Empty contig list provided to detect_naming_convention")
        return "unknown"

    assembly_detection = (config or {}).get("assembly_detection", {})
    threshold = assembly_detection.get("naming_convention_threshold", 0.5)
    if (
        isinstance(threshold, bool)
        or not isinstance(threshold, (int, float))
        or not math.isfinite(threshold)
        or not 0 <= threshold <= 1
    ):
        message = "Invalid naming_convention_threshold: expected a finite number from 0 through 1"
        logger.error(message)
        raise ValueError(message)

    primary_patterns = assembly_detection.get(
        "primary_contig_patterns",
        [r"^chr[0-9XYM]+$", r"^NC_\d{6}\.\d+$", r"^([0-9]+|X|Y|MT?)$"],
    )
    conventions = ("ucsc", "ncbi", "ensembl")
    if not isinstance(primary_patterns, list) or len(primary_patterns) != len(conventions):
        message = "Invalid primary_contig_patterns: expected three patterns in UCSC, NCBI, ENSEMBL order"
        logger.error(message)
        raise ValueError(message)

    compiled_patterns = []
    for convention, pattern in zip(conventions, primary_patterns, strict=True):
        if not isinstance(pattern, str) or not pattern:
            message = f"Invalid primary_contig_patterns: {convention} pattern must be a nonempty string"
            logger.error(message)
            raise ValueError(message)
        try:
            compiled_patterns.append(re.compile(pattern, re.IGNORECASE))
        except re.error as error:
            message = f"Invalid primary_contig_patterns: {convention} pattern cannot be compiled: {error}"
            logger.error(message)
            raise ValueError(message) from error

    counts = dict.fromkeys(conventions, 0)

    for name in contig_names:
        for convention, pattern in zip(conventions, compiled_patterns, strict=True):
            if pattern.match(name):
                counts[convention] += 1
                break

    classified_total = sum(counts.values())
    if classified_total == 0:
        logger.warning("Could not determine naming convention: no contigs matched a configured primary pattern")
        return "unknown"

    winning_count = max(counts.values())
    winners = [convention for convention, count in counts.items() if count == winning_count]
    if len(winners) != 1 or winning_count * 2 <= classified_total or winning_count / classified_total < threshold:
        logger.warning(
            "Could not determine naming convention. "
            f"UCSC: {counts['ucsc']}, ENSEMBL: {counts['ensembl']}, NCBI: {counts['ncbi']}; "
            f"classified contigs: {classified_total}"
        )
        return "unknown"

    convention = winners[0]
    logger.debug(f"Detected {convention.upper()} naming convention ({winning_count}/{classified_total} contigs)")
    return convention


def get_chromosome_name_from_bam(
    bam_file: str, config: dict, chromosome_number: int = 1, reference_assembly: str = "hg19"
) -> str:
    """
    Detect the actual chromosome name used in a BAM file for a given chromosome number.

    This function parses the BAM header, detects the naming convention, and returns
    the appropriate chromosome identifier based on what's actually in the file.

    Args:
        bam_file (str): Path to BAM/CRAM file
        config (dict): Main configuration dictionary
        chromosome_number (int): Chromosome number to resolve (default: 1)
        reference_assembly (str): Reference assembly for NCBI accession lookup

    Returns:
        str: The actual chromosome name (e.g., "chr1", "1", "NC_000001.10")

    Raises:
        ValueError: If the naming convention is ambiguous or unclassifiable,
            or if the chromosome cannot be found in the BAM header.
        FileNotFoundError: If BAM file doesn't exist

    Examples:
        >>> get_chromosome_name_from_bam("sample.bam", config, 1, "GRCh37")
        "NC_000001.10"
        >>> get_chromosome_name_from_bam("sample.bam", config, 1, "hg19")
        "chr1"
    """
    from vntyper.scripts.fastq_bam_processing import (
        extract_bam_header,
        parse_contigs_from_header,
    )

    # Extract BAM header
    try:
        header = extract_bam_header(bam_file, config)
    except Exception as e:
        logger.error(f"Failed to extract BAM header from {bam_file}: {e}")
        raise ValueError(f"Cannot read BAM header: {e}") from e

    # Parse contigs
    contigs = parse_contigs_from_header(header)
    if not contigs:
        raise ValueError(f"No contigs found in BAM header for {bam_file}")

    contig_names = [c["name"] for c in contigs]
    logger.debug(f"Found {len(contig_names)} contigs in BAM header. First 5: {contig_names[:5]}")

    # Detect naming convention
    convention = detect_naming_convention(contig_names, config)
    logger.debug(f"Detected naming convention: {convention}")
    if convention == "unknown":
        message = (
            f"{NAMING_CONVENTION_ERROR_PREFIX} in BAM header for {bam_file}; "
            f"cannot select chromosome {chromosome_number}. Available contigs: {', '.join(contig_names[:10])}..."
        )
        logger.error(message)
        raise ValueError(message)

    # Build expected chromosome name based on convention
    chr_name = _build_chromosome_name(chromosome_number, convention, reference_assembly, config)

    # Verify the chromosome exists in the BAM
    if chr_name not in contig_names:
        # Try to find a match with case-insensitive search
        chr_name_lower = chr_name.lower()
        for name in contig_names:
            if name.lower() == chr_name_lower:
                logger.debug(f"Found case-insensitive match: {name} for {chr_name}")
                return name

        # Chromosome not found - provide helpful error
        raise ValueError(
            f"Chromosome {chromosome_number} not found in BAM file {bam_file}. "
            f"Expected name: '{chr_name}' (convention: {convention}). "
            f"Available contigs: {', '.join(contig_names[:10])}..."
        )

    logger.debug(f"Resolved chromosome {chromosome_number} to '{chr_name}'")
    return chr_name


def _build_chromosome_name(chromosome_number: int, convention: str, reference_assembly: str, config: dict) -> str:
    """
    Build the expected chromosome name based on convention and assembly.

    Args:
        chromosome_number (int): Chromosome number (1-22, 23=X, 24=Y, 25=MT)
        convention (str): Naming convention ("ucsc", "ensembl", "ncbi")
        reference_assembly (str): Reference assembly name
        config (dict): Configuration dictionary

    Returns:
        str: Expected chromosome name

    Raises:
        ValueError: If convention is unknown or chromosome number invalid
    """
    # Validate chromosome number
    if not 1 <= chromosome_number <= 25:
        raise ValueError(f"Invalid chromosome number: {chromosome_number}. Must be 1-25.")

    # Handle special chromosomes
    if chromosome_number == 23:
        chr_suffix = "X"
    elif chromosome_number == 24:
        chr_suffix = "Y"
    elif chromosome_number == 25:
        chr_suffix = "M" if convention == "ucsc" else "MT"
    else:
        chr_suffix = str(chromosome_number)

    # Build name based on convention
    if convention == "ucsc":
        return f"chr{chr_suffix}"

    elif convention == "ensembl":
        return chr_suffix if chromosome_number >= 23 else str(chromosome_number)

    elif convention == "ncbi":
        # Look up NCBI accession from config
        known_naming = config.get("bam_processing", {}).get("known_chromosome_naming", {})

        # Map reference assembly to coordinate set (use registry)
        from vntyper.scripts.reference_registry import get_coordinate_system

        try:
            coord_assembly = get_coordinate_system(reference_assembly)
        except ValueError:
            logger.warning(f"Unknown assembly '{reference_assembly}', defaulting to GRCh37")
            coord_assembly = "GRCh37"

        # Get NCBI accession from config
        naming_info = known_naming.get(coord_assembly, {})
        ncbi_accession = naming_info.get("ncbi")

        if ncbi_accession and chromosome_number == 1:
            # Config provides chr1 accession, use it
            return ncbi_accession
        else:
            # Construct NCBI accession based on assembly and chr number
            return _construct_ncbi_accession(chromosome_number, coord_assembly)

    else:
        message = f"Unknown chromosome naming convention: '{convention}'"
        logger.error(message)
        raise ValueError(message)


def _construct_ncbi_accession(chromosome_number: int, assembly: str) -> str:
    """
    Construct NCBI RefSeq accession for a given chromosome and assembly.

    NCBI accession format: NC_XXXXXX.VERSION
    - Chr 1-22: NC_000001.XX - NC_000022.XX
    - Chr X: NC_000023.XX
    - Chr Y: NC_000024.XX
    - MT: NC_012920.1

    Both spellings of each build are accepted: the UCSC-style alias ("hg19",
    "hg38") and the build name ("GRCh37", "GRCh38"), plus every other alias the
    reference registry knows. `_build_chromosome_name` normalises to the build
    name before calling this, so accepting only "hg19" silently built every
    GRCh37 accession from the GRCh38 table.

    Args:
        chromosome_number (int): Chromosome number (1-25)
        assembly (str): Assembly name or alias (e.g. "hg19", "GRCh37",
            "hg38_ncbi"). An unrecognised name falls back to GRCh37, matching
            `_build_chromosome_name` and `region_utils.resolve_assembly_alias`.

    Returns:
        str: NCBI accession (e.g., "NC_000001.10")

    Raises:
        ValueError: If the chromosome number has no accession.
    """
    # Define NCBI versions for each assembly
    # GRCh37/hg19 versions
    grch37_versions = {
        1: "10",
        2: "11",
        3: "11",
        4: "11",
        5: "9",
        6: "11",
        7: "13",
        8: "10",
        9: "11",
        10: "10",
        11: "9",
        12: "11",
        13: "10",
        14: "8",
        15: "9",
        16: "9",
        17: "10",
        18: "9",
        19: "9",
        20: "10",
        21: "8",
        22: "10",
        23: "10",
        24: "9",
        25: "1",  # 23=X, 24=Y, 25=MT
    }

    # GRCh38/hg38 versions
    grch38_versions = {
        1: "11",
        2: "12",
        3: "12",
        4: "12",
        5: "10",
        6: "12",
        7: "14",
        8: "11",
        9: "12",
        10: "11",
        11: "10",
        12: "12",
        13: "11",
        14: "9",
        15: "10",
        16: "10",
        17: "11",
        18: "10",
        19: "10",
        20: "11",
        21: "9",
        22: "11",
        23: "11",
        24: "10",
        25: "1",  # 23=X, 24=Y, 25=MT
    }

    # Select version table by coordinate system, not by one hardcoded spelling.
    from vntyper.scripts.reference_registry import get_coordinate_system

    try:
        coordinate_system = get_coordinate_system(assembly)
    except ValueError:
        logger.warning(f"Unknown assembly '{assembly}' for NCBI accession lookup, defaulting to GRCh37")
        coordinate_system = "GRCh37"

    versions = grch37_versions if coordinate_system == "GRCh37" else grch38_versions
    version = versions.get(chromosome_number)

    if version is None:
        raise ValueError(f"No NCBI version for chromosome {chromosome_number}")

    # Construct accession
    if chromosome_number == 25:  # Mitochondrial
        return "NC_012920.1"
    elif chromosome_number <= 24:
        # Format: NC_000001 through NC_000024 (zero-padded to 6 digits)
        accession_base = f"NC_{chromosome_number:06d}"
        return f"{accession_base}.{version}"
    else:
        raise ValueError(f"Invalid chromosome number: {chromosome_number}")


def validate_chromosome_name(chromosome_name: str) -> bool:
    """
    Validate that a chromosome name follows expected patterns.

    Args:
        chromosome_name (str): Chromosome name to validate

    Returns:
        bool: True if valid, False otherwise

    Examples:
        >>> validate_chromosome_name("chr1")
        True
        >>> validate_chromosome_name("NC_000001.10")
        True
        >>> validate_chromosome_name("invalid_chr")
        False
    """
    if not chromosome_name:
        return False

    # Valid patterns
    patterns = [
        r"^chr[0-9]+$",  # UCSC: chr1, chr2, ...
        r"^chr[XYM]$",  # UCSC: chrX, chrY, chrM
        r"^[0-9]+$",  # ENSEMBL: 1, 2, ...
        r"^[XYMT]+$",  # ENSEMBL: X, Y, MT
        r"^NC_\d{6}\.\d+$",  # NCBI: NC_000001.10
    ]

    for pattern in patterns:
        if re.match(pattern, chromosome_name, re.IGNORECASE):
            return True

    logger.debug(f"Chromosome name '{chromosome_name}' did not match any pattern")
    return False
