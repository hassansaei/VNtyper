#!/usr/bin/env python3
"""
reference_registry.py

Centralized reference assembly registry for VNtyper.

This module provides a single source of truth for:
1. Assembly name aliases and deprecation handling
2. Coordinate systems (biological coordinate spaces)
3. Reference source metadata (UCSC, NCBI, ENSEMBL)

Architecture:
- COORDINATE_SYSTEMS: Define coordinate ranges (GRCh37, GRCh38)
- REFERENCE_SOURCES: Define chromosome naming patterns (ucsc, ncbi, ensembl)
- ASSEMBLY_ALIASES: Map user inputs to canonical names
- ASSEMBLY_METADATA: Combine coordinate systems with reference sources

This design separates coordinate systems from reference sources,
eliminating duplication and improving maintainability.
"""

import logging
from typing import Optional, TypedDict, cast

# =============================================================================
# Type Definitions
# =============================================================================


class AssemblyMetadataDict(TypedDict):
    """Type definition for assembly metadata dictionary."""

    coordinate_system: str
    reference_source: str
    description: str
    deprecated: bool


class CoordinateSystemDict(TypedDict):
    """Type definition for coordinate system dictionary."""

    chromosome: int
    bam_region_coords: str
    vntr_region_coords: str


# =============================================================================
# Coordinate Systems (Biological Truth - Single Source of Truth)
# =============================================================================

COORDINATE_SYSTEMS = {
    "GRCh37": {
        "chromosome": 1,
        "bam_region_coords": "155158000-155163000",
        "vntr_region_coords": "155160500-155162000",
    },
    "GRCh38": {
        "chromosome": 1,
        "bam_region_coords": "155184000-155194000",
        "vntr_region_coords": "155188000-155192500",
    },
}

# =============================================================================
# Reference Sources (Chromosome Naming Patterns)
# =============================================================================

REFERENCE_SOURCES = {
    "ucsc": {
        "description": "UCSC naming convention",
        "chr_format": "chr{n}",  # chr1, chr2, chrX, chrY, chrM
        "example": "chr1",
        "download_base_url": "https://hgdownload.soe.ucsc.edu/goldenPath/",
    },
    "ncbi": {
        "description": "NCBI RefSeq accession naming",
        "chr_format": "NC_{accession}",  # NC_000001.10, NC_000001.11
        "example": "NC_000001.10",
        "download_base_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/",
    },
    "ensembl": {
        "description": "ENSEMBL simple numeric naming",
        "chr_format": "{n}",  # 1, 2, X, Y, MT
        "example": "1",
        "download_base_url": "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/",
    },
}

# =============================================================================
# Known NCBI Accessions (for chr1)
# =============================================================================

KNOWN_NCBI_ACCESSIONS = {
    "GRCh37": "NC_000001.10",
    "GRCh38": "NC_000001.11",
}

# =============================================================================
# Assembly Metadata (Combines Coordinate System + Reference Source)
# =============================================================================

ASSEMBLY_METADATA = {
    # UCSC references (chr1, chr2, etc.)
    "hg19": {
        "coordinate_system": "GRCh37",
        "reference_source": "ucsc",
        "description": "UCSC hg19 (GRCh37 coordinates, chr prefix)",
        "deprecated": False,
    },
    "hg38": {
        "coordinate_system": "GRCh38",
        "reference_source": "ucsc",
        "description": "UCSC hg38 (GRCh38 coordinates, chr prefix)",
        "deprecated": False,
    },
    # Biological assembly names (NCBI)
    "GRCh37": {
        "coordinate_system": "GRCh37",
        "reference_source": "ncbi",
        "description": "GRCh37 biological assembly (NCBI RefSeq accessions)",
        "deprecated": False,
    },
    "GRCh38": {
        "coordinate_system": "GRCh38",
        "reference_source": "ncbi",
        "description": "GRCh38 biological assembly (NCBI RefSeq accessions)",
        "deprecated": False,
    },
    # NCBI references with explicit naming
    "hg19_ncbi": {
        "coordinate_system": "GRCh37",
        "reference_source": "ncbi",
        "description": "NCBI GRCh37 (GRCh37 coordinates, RefSeq accessions)",
        "deprecated": False,
    },
    "hg38_ncbi": {
        "coordinate_system": "GRCh38",
        "reference_source": "ncbi",
        "description": "NCBI GRCh38 (GRCh38 coordinates, RefSeq accessions)",
        "deprecated": False,
    },
    # ENSEMBL references (1, 2, X, Y, MT)
    "hg19_ensembl": {
        "coordinate_system": "GRCh37",
        "reference_source": "ensembl",
        "description": "ENSEMBL GRCh37 (GRCh37 coordinates, simple numeric)",
        "deprecated": False,
    },
    "hg38_ensembl": {
        "coordinate_system": "GRCh38",
        "reference_source": "ensembl",
        "description": "ENSEMBL GRCh38 (GRCh38 coordinates, simple numeric)",
        "deprecated": False,
    },
}

# =============================================================================
# Assembly Aliases (Backward Compatibility + Deprecation)
# =============================================================================

ASSEMBLY_ALIASES = {
    # Canonical UCSC names
    "hg19": "hg19",
    "hg38": "hg38",
    # Canonical biological assembly names
    "GRCh37": "GRCh37",
    "GRCh38": "GRCh38",
    # Canonical explicit source names
    "hg19_ncbi": "hg19_ncbi",
    "hg38_ncbi": "hg38_ncbi",
    "hg19_ensembl": "hg19_ensembl",
    "hg38_ensembl": "hg38_ensembl",
}

# =============================================================================
# Public API Functions
# =============================================================================


def normalize_assembly_name(user_input: str, warn_deprecated: bool = True) -> str:
    """
    Normalize assembly name to canonical form.

    Args:
        user_input (str): User-provided assembly name
        warn_deprecated (bool): Reserved for future use (default: True)

    Returns:
        str: Canonical assembly name (e.g., "hg19", "hg19_ncbi", "hg19_ensembl")

    Raises:
        ValueError: If assembly name is not recognized

    Examples:
        >>> normalize_assembly_name("hg19")
        "hg19"
        >>> normalize_assembly_name("GRCh37")
        "GRCh37"
        >>> normalize_assembly_name("hg19_ensembl")
        "hg19_ensembl"
    """
    if user_input not in ASSEMBLY_ALIASES:
        supported = ", ".join(sorted(ASSEMBLY_ALIASES.keys()))
        raise ValueError(f"Unknown assembly '{user_input}'. Supported assemblies: {supported}")

    return ASSEMBLY_ALIASES[user_input]


def get_coordinate_system(assembly_name: str) -> str:
    """
    Get the coordinate system for a given assembly.

    Args:
        assembly_name (str): Assembly name (canonical or alias)

    Returns:
        str: Coordinate system name ("GRCh37" or "GRCh38")

    Raises:
        ValueError: If assembly not found

    Examples:
        >>> get_coordinate_system("hg19")
        "GRCh37"
        >>> get_coordinate_system("hg38_ensembl")
        "GRCh38"
    """
    canonical = normalize_assembly_name(assembly_name, warn_deprecated=False)
    metadata = ASSEMBLY_METADATA.get(canonical)

    if not metadata:
        raise ValueError(f"No metadata found for assembly '{canonical}'")

    # Cast is safe here because we validate the structure in validate_registry()
    return cast(str, metadata["coordinate_system"])


def get_reference_source(assembly_name: str) -> str:
    """
    Get the reference source for a given assembly.

    Args:
        assembly_name (str): Assembly name (canonical or alias)

    Returns:
        str: Reference source ("ucsc", "ncbi", or "ensembl")

    Raises:
        ValueError: If assembly not found

    Examples:
        >>> get_reference_source("hg19")
        "ucsc"
        >>> get_reference_source("hg19_ncbi")
        "ncbi"
        >>> get_reference_source("hg19_ensembl")
        "ensembl"
    """
    canonical = normalize_assembly_name(assembly_name, warn_deprecated=False)
    metadata = ASSEMBLY_METADATA.get(canonical)

    if not metadata:
        raise ValueError(f"No metadata found for assembly '{canonical}'")

    # Cast is safe here because we validate the structure in validate_registry()
    return cast(str, metadata["reference_source"])


def get_coordinates(assembly_name: str, region_type: str) -> str:
    """
    Get coordinate range for a given assembly and region type.

    Args:
        assembly_name (str): Assembly name (canonical or alias)
        region_type (str): Region type ("bam_region_coords" or "vntr_region_coords")

    Returns:
        str: Coordinate range in "start-end" format (e.g., "155158000-155163000")

    Raises:
        ValueError: If assembly or region type not found

    Examples:
        >>> get_coordinates("hg19", "bam_region_coords")
        "155158000-155163000"
        >>> get_coordinates("hg38_ensembl", "vntr_region_coords")
        "155188000-155192500"
    """
    coord_system = get_coordinate_system(assembly_name)
    coord_data = COORDINATE_SYSTEMS.get(coord_system)

    if not coord_data:
        raise ValueError(f"No coordinates found for coordinate system '{coord_system}'")

    coordinates = coord_data.get(region_type)

    if not coordinates:
        available = ", ".join(coord_data.keys())
        raise ValueError(f"Region type '{region_type}' not found. Available types: {available}")

    # Cast is safe here because we validate the structure in validate_registry()
    return cast(str, coordinates)


def get_assembly_metadata(assembly_name: str) -> dict:
    """
    Get complete metadata for a given assembly.

    Args:
        assembly_name (str): Assembly name (canonical or alias)

    Returns:
        dict: Assembly metadata including coordinate_system, reference_source, description

    Raises:
        ValueError: If assembly not found

    Examples:
        >>> meta = get_assembly_metadata("hg19_ensembl")
        >>> meta["coordinate_system"]
        "GRCh37"
        >>> meta["reference_source"]
        "ensembl"
    """
    canonical = normalize_assembly_name(assembly_name, warn_deprecated=False)
    metadata = ASSEMBLY_METADATA.get(canonical)

    if not metadata:
        raise ValueError(f"No metadata found for assembly '{canonical}'")

    return metadata.copy()


def get_reference_source_info(source_name: str) -> dict:
    """
    Get information about a reference source.

    Args:
        source_name (str): Reference source name ("ucsc", "ncbi", "ensembl")

    Returns:
        dict: Reference source metadata

    Raises:
        ValueError: If source not found

    Examples:
        >>> info = get_reference_source_info("ensembl")
        >>> info["chr_format"]
        "{n}"
        >>> info["example"]
        "1"
    """
    source_info = REFERENCE_SOURCES.get(source_name)

    if not source_info:
        available = ", ".join(REFERENCE_SOURCES.keys())
        raise ValueError(f"Unknown reference source '{source_name}'. Available: {available}")

    return source_info.copy()


def list_assemblies(include_deprecated: bool = False) -> list:
    """
    List all available assembly names.

    Args:
        include_deprecated (bool): Reserved for future use (default: False)

    Returns:
        list: List of assembly names

    Examples:
        >>> list_assemblies()
        ['GRCh37', 'GRCh38', 'hg19', 'hg19_ensembl', 'hg19_ncbi', 'hg38', 'hg38_ensembl', 'hg38_ncbi']
    """
    # All names are canonical now
    return sorted(ASSEMBLY_ALIASES.keys())


def is_deprecated(assembly_name: str) -> bool:
    """
    Check if an assembly name is deprecated.

    Args:
        assembly_name (str): Assembly name to check

    Returns:
        bool: Always False (no deprecated names in this version)

    Examples:
        >>> is_deprecated("hg19")
        False
        >>> is_deprecated("GRCh37")
        False
    """
    # No deprecated names in current version
    return False


def get_all_coordinate_systems() -> list:
    """
    Get all available coordinate systems.

    Returns:
        list: List of coordinate system names

    Examples:
        >>> get_all_coordinate_systems()
        ['GRCh37', 'GRCh38']
    """
    return list(COORDINATE_SYSTEMS.keys())


def get_all_reference_sources() -> list:
    """
    Get all available reference sources.

    Returns:
        list: List of reference source names

    Examples:
        >>> get_all_reference_sources()
        ['ucsc', 'ncbi', 'ensembl']
    """
    return list(REFERENCE_SOURCES.keys())


def resolve_chromosome_name(
    assembly_name: str, chromosome_number: int = 1, detected_convention: Optional[str] = None
) -> str:
    """
    Resolve chromosome name based on assembly and optional detected convention.

    If detected_convention is provided, it overrides the assembly's default source.
    This allows dynamic chromosome naming based on BAM file contents.

    Args:
        assembly_name (str): Assembly name (canonical or alias)
        chromosome_number (int): Chromosome number (1-25, where 23=X, 24=Y, 25=MT)
        detected_convention (str, optional): Detected naming convention from BAM
            ("ucsc", "ncbi", "ensembl")

    Returns:
        str: Chromosome name (e.g., "chr1", "NC_000001.10", "1")

    Raises:
        ValueError: If inputs are invalid

    Examples:
        >>> resolve_chromosome_name("hg19", 1)
        "chr1"
        >>> resolve_chromosome_name("hg19_ncbi", 1)
        "NC_000001.10"
        >>> resolve_chromosome_name("hg19", 1, detected_convention="ensembl")
        "1"
    """
    # Determine which reference source to use
    source = detected_convention or get_reference_source(assembly_name)

    # Get coordinate system for NCBI accession lookup
    coord_system = get_coordinate_system(assembly_name)

    # Handle special chromosomes
    if chromosome_number == 23:
        chr_suffix = "X"
    elif chromosome_number == 24:
        chr_suffix = "Y"
    elif chromosome_number == 25:
        chr_suffix = "M" if source == "ucsc" else "MT"
    elif 1 <= chromosome_number <= 22:
        chr_suffix = str(chromosome_number)
    else:
        raise ValueError(f"Invalid chromosome number: {chromosome_number}. Must be 1-25.")

    # Build chromosome name based on source
    if source == "ucsc":
        return f"chr{chr_suffix}"
    elif source == "ensembl":
        return chr_suffix if chromosome_number >= 23 else str(chromosome_number)
    elif source == "ncbi":
        # Look up NCBI accession
        if chromosome_number == 1:
            return KNOWN_NCBI_ACCESSIONS.get(coord_system, "NC_000001.10")
        else:
            # For other chromosomes, would need full accession table
            # For now, raise NotImplementedError
            raise NotImplementedError(f"NCBI accession for chromosome {chromosome_number} not yet implemented")
    else:
        raise ValueError(f"Unknown reference source: {source}")


# =============================================================================
# Validation Functions
# =============================================================================


def validate_registry() -> tuple[bool, list]:
    """
    Validate the internal consistency of the registry.

    Returns:
        tuple: (is_valid, list_of_errors)

    Examples:
        >>> is_valid, errors = validate_registry()
        >>> if not is_valid:
        ...     print(errors)
    """
    errors = []

    # Check that all canonical names in ASSEMBLY_ALIASES exist in ASSEMBLY_METADATA
    for alias, canonical in ASSEMBLY_ALIASES.items():
        if canonical not in ASSEMBLY_METADATA:
            errors.append(f"Canonical name '{canonical}' from alias '{alias}' not found in ASSEMBLY_METADATA")

    # Check that all coordinate systems in ASSEMBLY_METADATA exist in COORDINATE_SYSTEMS
    for assembly, metadata in ASSEMBLY_METADATA.items():
        coord_system = metadata.get("coordinate_system")
        if coord_system not in COORDINATE_SYSTEMS:
            errors.append(f"Coordinate system '{coord_system}' for assembly '{assembly}' not found")

    # Check that all reference sources in ASSEMBLY_METADATA exist in REFERENCE_SOURCES
    for assembly, metadata in ASSEMBLY_METADATA.items():
        ref_source = metadata.get("reference_source")
        if ref_source not in REFERENCE_SOURCES:
            errors.append(f"Reference source '{ref_source}' for assembly '{assembly}' not found")

    return len(errors) == 0, errors


# Run validation on module import (optional - can be disabled in production)
if __name__ != "__main__":
    is_valid, errors = validate_registry()
    if not is_valid:
        logging.warning(f"Reference registry validation failed with {len(errors)} errors:")
        for error in errors:
            logging.warning(f"  - {error}")
