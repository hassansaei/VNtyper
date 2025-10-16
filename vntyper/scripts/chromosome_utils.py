#!/usr/bin/env python3
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

import logging
import re
from typing import List


def detect_naming_convention(contig_names: List[str]) -> str:
    """
    Detect the chromosome naming convention from a list of contig names.

    This function examines contig names to determine if they follow:
    - UCSC convention: chr1, chr2, ..., chrX, chrY, chrM
    - No-prefix numeric: 1, 2, ..., X, Y, MT (used by Ensembl/GRCh)
    - NCBI accessions: NC_000001.XX, NC_000002.XX, ...

    Args:
        contig_names (List[str]): List of contig names from BAM header

    Returns:
        str: Convention identifier ("ucsc", "no_prefix", "ncbi", "unknown")

    Examples:
        >>> detect_naming_convention(["chr1", "chr2", "chrX"])
        'ucsc'
        >>> detect_naming_convention(["1", "2", "X"])
        'no_prefix'
        >>> detect_naming_convention(["NC_000001.10", "NC_000002.11"])
        'ncbi'
    """
    if not contig_names:
        logging.warning("Empty contig list provided to detect_naming_convention")
        return "unknown"

    # Count naming patterns
    ucsc_count = 0
    no_prefix_count = 0
    ncbi_count = 0

    for name in contig_names:
        # UCSC: starts with "chr" followed by number/letter
        if re.match(r'^chr[0-9XYM]+$', name, re.IGNORECASE):
            ucsc_count += 1
        # NCBI: NC_XXXXXX.YY format
        elif re.match(r'^NC_\d{6}\.\d+$', name):
            ncbi_count += 1
        # No-prefix numeric: just digits or X, Y, MT
        elif re.match(r'^([0-9]+|X|Y|MT?)$', name, re.IGNORECASE):
            no_prefix_count += 1

    # Determine convention based on majority
    total = len(contig_names)
    threshold = 0.5  # At least 50% of contigs should match the pattern

    if ucsc_count / total >= threshold:
        logging.debug(
            f"Detected UCSC naming convention ({ucsc_count}/{total} contigs)"
        )
        return "ucsc"
    elif ncbi_count / total >= threshold:
        logging.debug(
            f"Detected NCBI naming convention ({ncbi_count}/{total} contigs)"
        )
        return "ncbi"
    elif no_prefix_count / total >= threshold:
        logging.debug(
            f"Detected no-prefix naming convention ({no_prefix_count}/{total} contigs)"
        )
        return "no_prefix"
    else:
        logging.warning(
            f"Could not determine naming convention. "
            f"UCSC: {ucsc_count}, No-prefix: {no_prefix_count}, NCBI: {ncbi_count}"
        )
        return "unknown"


def get_chromosome_name_from_bam(
    bam_file: str,
    config: dict,
    chromosome_number: int = 1,
    reference_assembly: str = "hg19"
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
        ValueError: If chromosome cannot be found in BAM header
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
        logging.error(f"Failed to extract BAM header from {bam_file}: {e}")
        raise ValueError(f"Cannot read BAM header: {e}")

    # Parse contigs
    contigs = parse_contigs_from_header(header)
    if not contigs:
        raise ValueError(f"No contigs found in BAM header for {bam_file}")

    contig_names = [c["name"] for c in contigs]
    logging.debug(
        f"Found {len(contig_names)} contigs in BAM header. "
        f"First 5: {contig_names[:5]}"
    )

    # Detect naming convention
    convention = detect_naming_convention(contig_names)
    logging.debug(f"Detected naming convention: {convention}")

    # Build expected chromosome name based on convention
    chr_name = _build_chromosome_name(
        chromosome_number, convention, reference_assembly, config
    )

    # Verify the chromosome exists in the BAM
    if chr_name not in contig_names:
        # Try to find a match with case-insensitive search
        chr_name_lower = chr_name.lower()
        for name in contig_names:
            if name.lower() == chr_name_lower:
                logging.debug(
                    f"Found case-insensitive match: {name} for {chr_name}"
                )
                return name

        # Chromosome not found - provide helpful error
        raise ValueError(
            f"Chromosome {chromosome_number} not found in BAM file {bam_file}. "
            f"Expected name: '{chr_name}' (convention: {convention}). "
            f"Available contigs: {', '.join(contig_names[:10])}..."
        )

    logging.debug(f"Resolved chromosome {chromosome_number} to '{chr_name}'")
    return chr_name


def _build_chromosome_name(
    chromosome_number: int,
    convention: str,
    reference_assembly: str,
    config: dict
) -> str:
    """
    Build the expected chromosome name based on convention and assembly.

    Args:
        chromosome_number (int): Chromosome number (1-22, 23=X, 24=Y, 25=MT)
        convention (str): Naming convention ("ucsc", "no_prefix", "ncbi")
        reference_assembly (str): Reference assembly name
        config (dict): Configuration dictionary

    Returns:
        str: Expected chromosome name

    Raises:
        ValueError: If convention is unknown or chromosome number invalid
    """
    # Validate chromosome number
    if not 1 <= chromosome_number <= 25:
        raise ValueError(
            f"Invalid chromosome number: {chromosome_number}. Must be 1-25."
        )

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

    elif convention == "no_prefix":
        return chr_suffix if chromosome_number >= 23 else str(chromosome_number)

    elif convention == "ncbi":
        # Look up NCBI accession from config
        known_naming = config.get("bam_processing", {}).get(
            "known_chromosome_naming", {}
        )

        # Map reference assembly to coordinate set (hg19=GRCh37, hg38=GRCh38)
        assembly_map = {
            "hg19": "hg19",
            "GRCh37": "hg19",
            "hg38": "hg38",
            "GRCh38": "hg38",
            "hg19_nochr": "hg19",
            "hg38_nochr": "hg38"
        }
        coord_assembly = assembly_map.get(reference_assembly, "hg19")

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
        # Unknown convention - return no-prefix format as fallback
        logging.warning(
            f"Unknown naming convention '{convention}', using no-prefix format"
        )
        return chr_suffix if chromosome_number >= 23 else str(chromosome_number)


def _construct_ncbi_accession(chromosome_number: int, assembly: str) -> str:
    """
    Construct NCBI RefSeq accession for a given chromosome and assembly.

    NCBI accession format: NC_XXXXXX.VERSION
    - Chr 1-22: NC_000001.XX - NC_000022.XX
    - Chr X: NC_000023.XX
    - Chr Y: NC_000024.XX
    - MT: NC_012920.1

    Args:
        chromosome_number (int): Chromosome number (1-25)
        assembly (str): Assembly name (hg19 or hg38)

    Returns:
        str: NCBI accession (e.g., "NC_000001.10")
    """
    # Define NCBI versions for each assembly
    # GRCh37/hg19 versions
    grch37_versions = {
        1: "10", 2: "11", 3: "11", 4: "11", 5: "9",
        6: "11", 7: "13", 8: "10", 9: "11", 10: "10",
        11: "9", 12: "11", 13: "10", 14: "8", 15: "9",
        16: "9", 17: "10", 18: "9", 19: "9", 20: "10",
        21: "8", 22: "10", 23: "10", 24: "9", 25: "1"  # 23=X, 24=Y, 25=MT
    }

    # GRCh38/hg38 versions
    grch38_versions = {
        1: "11", 2: "12", 3: "12", 4: "12", 5: "10",
        6: "12", 7: "14", 8: "11", 9: "12", 10: "11",
        11: "10", 12: "12", 13: "11", 14: "9", 15: "10",
        16: "10", 17: "11", 18: "10", 19: "10", 20: "11",
        21: "9", 22: "11", 23: "11", 24: "10", 25: "1"  # 23=X, 24=Y, 25=MT
    }

    # Select version table based on assembly
    versions = grch37_versions if assembly == "hg19" else grch38_versions
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
        r'^chr[0-9]+$',  # UCSC: chr1, chr2, ...
        r'^chr[XYM]$',   # UCSC: chrX, chrY, chrM
        r'^[0-9]+$',     # Simple: 1, 2, ...
        r'^[XYMT]+$',    # Simple: X, Y, MT
        r'^NC_\d{6}\.\d+$'  # NCBI: NC_000001.10
    ]

    for pattern in patterns:
        if re.match(pattern, chromosome_name, re.IGNORECASE):
            return True

    logging.debug(f"Chromosome name '{chromosome_name}' did not match any pattern")
    return False
