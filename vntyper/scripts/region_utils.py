"""
region_utils.py

Region string construction utilities.

This module provides functions to build genomic region strings
dynamically based on BAM file chromosome naming conventions.
It enables VNtyper to work with multiple chromosome naming schemes
(UCSC, simple numeric, NCBI accessions) without hardcoding region strings.

Functions:
    get_region_string: Main function to build region strings dynamically
    build_region_string: Combine chromosome name and coordinates
    resolve_assembly_alias: Map assembly aliases to coordinate sets
    get_region_string_with_fallback: Backwards-compatible region resolution
"""

import logging

logger = logging.getLogger(__name__)


def get_region_string(bam_file: str, reference_assembly: str, region_type: str, config: dict) -> str:
    """
    Build a region string dynamically based on BAM chromosome naming.

    This function:
    1. Resolves the assembly alias (hg19=GRCh37, hg38=GRCh38)
    2. Gets coordinates from config for the resolved assembly
    3. Detects the chromosome naming used in the BAM file
    4. Constructs the final region string (e.g., "NC_000001.10:155158000-155163000")

    Args:
        bam_file (str): Path to input BAM/CRAM file
        reference_assembly (str): User-specified assembly
            ("hg19", "GRCh37", "hg38", "GRCh38")
        region_type (str): Type of region
            ("bam_region_coords" or "vntr_region_coords")
        config (dict): Main configuration dictionary

    Returns:
        str: Region string in format "chr:start-end"
            (e.g., "NC_000001.10:155158000-155163000")

    Raises:
        ValueError: If assembly not supported or chromosome not found
        KeyError: If configuration is missing required keys

    Examples:
        >>> get_region_string("sample.bam", "GRCh37", "bam_region_coords", config)
        "NC_000001.10:155158000-155163000"
        >>> get_region_string("sample.bam", "hg38", "vntr_region_coords", config)
        "chr1:155188000-155192500"
    """
    from vntyper.scripts.chromosome_utils import get_chromosome_name_from_bam

    # Resolve assembly alias to coordinate set
    coord_assembly = resolve_assembly_alias(reference_assembly)
    logger.debug(f"Resolved assembly '{reference_assembly}' to coord set '{coord_assembly}'")

    # Get coordinates from config
    assemblies = config.get("bam_processing", {}).get("assemblies", {})
    if not assemblies:
        raise KeyError("Configuration missing 'bam_processing' -> 'assemblies' section")

    assembly_config = assemblies.get(coord_assembly)
    if not assembly_config:
        raise ValueError(
            f"Assembly '{coord_assembly}' not found in configuration. Supported assemblies: {list(assemblies.keys())}"
        )

    coordinates = assembly_config.get(region_type)
    if not coordinates:
        raise KeyError(
            f"Region type '{region_type}' not found for assembly '{coord_assembly}'. "
            f"Available types: {list(assembly_config.keys())}"
        )

    # Get chromosome number (currently hardcoded to 1, could be made configurable)
    chromosome_number = assembly_config.get("chromosome", 1)

    chromosome_name = get_chromosome_name_from_bam(
        bam_file=bam_file, config=config, chromosome_number=chromosome_number, reference_assembly=reference_assembly
    )

    # Build final region string
    region = build_region_string(chromosome_name, coordinates)
    logger.debug(f"Built region string: {region} for {region_type} in assembly {reference_assembly}")

    return region


def build_region_string(chromosome_name: str, coordinates: str) -> str:
    """
    Combine chromosome name and coordinates into a region string.

    Args:
        chromosome_name (str): Chromosome identifier
            (e.g., "chr1", "1", "NC_000001.10")
        coordinates (str): Coordinate range in "start-end" format
            (e.g., "155158000-155163000")

    Returns:
        str: Complete region string (e.g., "chr1:155158000-155163000")

    Raises:
        ValueError: If coordinate format is invalid

    Examples:
        >>> build_region_string("chr1", "155158000-155163000")
        "chr1:155158000-155163000"
        >>> build_region_string("NC_000001.10", "155158000-155163000")
        "NC_000001.10:155158000-155163000"
    """
    if not chromosome_name or not coordinates:
        raise ValueError(f"Invalid inputs: chromosome_name='{chromosome_name}', coordinates='{coordinates}'")

    # Validate coordinate format
    if "-" not in coordinates:
        raise ValueError(f"Invalid coordinate format: '{coordinates}'. Expected format: 'start-end'")

    try:
        start, end = coordinates.split("-")
        int(start)  # Validate as integer
        int(end)  # Validate as integer
    except ValueError as e:
        raise ValueError(f"Invalid coordinate values in '{coordinates}': {e}") from e

    return f"{chromosome_name}:{coordinates}"


def resolve_assembly_alias(reference_assembly: str) -> str:
    """
    Map assembly aliases to their canonical coordinate set names.

    Uses the centralized reference registry to resolve assemblies to
    their coordinate systems (GRCh37 or GRCh38).

    Args:
        reference_assembly (str): User-specified assembly name

    Returns:
        str: Coordinate system name ("GRCh37" or "GRCh38")

    Examples:
        >>> resolve_assembly_alias("hg19")
        "GRCh37"
        >>> resolve_assembly_alias("GRCh37")
        "GRCh37"
        >>> resolve_assembly_alias("hg38")
        "GRCh38"
        >>> resolve_assembly_alias("GRCh38")
        "GRCh38"
    """
    from vntyper.scripts.reference_registry import get_coordinate_system

    try:
        return get_coordinate_system(reference_assembly)
    except ValueError as e:
        logger.warning(f"Unknown assembly '{reference_assembly}', defaulting to 'GRCh37': {e}")
        return "GRCh37"


def get_region_string_with_fallback(bam_file: str, reference_assembly: str, region_type: str, config: dict) -> str:
    """
    Get region string with fallback to legacy config format.

    This function attempts to use the new dynamic resolution method first.
    If that fails (e.g., due to missing config keys), it falls back to
    the old hardcoded region lookup for backwards compatibility.

    Args:
        bam_file (str): Path to BAM/CRAM file
        reference_assembly (str): Reference assembly name
        region_type (str): Region type (with or without "_coords" suffix)
        config (dict): Configuration dictionary

    Returns:
        str: Region string

    Raises:
        ValueError: If region cannot be resolved by either method, or if the
            legacy region names a contig the BAM demonstrably does not contain.
    """
    try:
        # Normalize region_type to include "_coords" suffix
        region_type_with_coords = f"{region_type}_coords" if not region_type.endswith("_coords") else region_type

        # Try new dynamic resolution
        return get_region_string(
            bam_file=bam_file, reference_assembly=reference_assembly, region_type=region_type_with_coords, config=config
        )

    except (KeyError, ValueError) as e:
        from vntyper.scripts.chromosome_utils import NAMING_CONVENTION_ERROR_PREFIX

        if isinstance(e, ValueError) and str(e).startswith(NAMING_CONVENTION_ERROR_PREFIX):
            raise
        logger.warning(f"Dynamic region resolution failed: {e}. Falling back to legacy config lookup.")

        # Fall back to old method: look up hardcoded region in config
        region_key = f"{region_type.replace('_coords', '')}_{reference_assembly}"
        region = config.get("bam_processing", {}).get(region_key)

        if not region:
            raise ValueError(
                f"Region not found in configuration. "
                f"Tried key: '{region_key}'. "
                f"Neither new nor legacy format available."
            ) from e

        _reject_region_absent_from_bam(bam_file=bam_file, region=region, region_key=region_key, config=config)

        logger.info(f"Using legacy region format: {region_key} = {region}")
        return region


def _bam_contig_names(bam_file: str, config: dict) -> list[str]:
    """
    Read the contig names out of a BAM/CRAM header.

    Args:
        bam_file (str): Path to the BAM/CRAM file.
        config (dict): Main configuration dictionary, for the samtools path.

    Returns:
        list[str]: Contig names, or an empty list if the header could not be read
        or carried no ``@SQ`` lines. An empty list means "could not determine",
        not "the file has no contigs" - callers must not treat it as evidence.
    """
    # Imported here rather than at module scope: fastq_bam_processing imports this
    # module, so a top-level import would be circular.
    from vntyper.scripts.fastq_bam_processing import (
        extract_bam_header,
        parse_contigs_from_header,
    )

    try:
        header = extract_bam_header(bam_file, config)
    except Exception as e:
        logger.warning(f"Could not read the header of {bam_file} to verify the fallback region: {e}")
        return []

    return [contig["name"] for contig in parse_contigs_from_header(header)]


def _reject_region_absent_from_bam(bam_file: str, region: str, region_key: str, config: dict) -> None:
    """
    Refuse a legacy fallback region whose contig is not in the BAM header.

    The legacy keys in ``config.json`` are per-spelling, not per-naming-convention:
    ``bam_region_hg38`` is UCSC-named while ``bam_region_GRCh38`` carries an NCBI
    accession. Falling back to the former on an NCBI- or Ensembl-named BAM produces
    a region string that matches nothing. ``samtools view`` returns zero records for
    a contig it does not have **and exits 0**, so the slice is empty, Kestrel sees
    no reads and the report states a confident negative. This turns that into a
    failure.

    Args:
        bam_file (str): Path to the BAM/CRAM file.
        region (str): The legacy region string about to be returned.
        region_key (str): The config key it came from, for the error message.
        config (dict): Main configuration dictionary.

    Raises:
        ValueError: If the header could be read, listed contigs, and none of them
            matches the region's contig (case-insensitively).

    Note:
        Silence is the default when the header cannot be read or lists no contigs.
        Refusing on no evidence would break working runs; this only fires when the
        BAM has positively told us the contig is not there.
    """
    contig = region.split(":")[0]
    contig_names = _bam_contig_names(bam_file, config)

    if not contig_names:
        logger.warning(
            f"Cannot verify that '{contig}' exists in {bam_file}; using legacy region '{region}' unverified."
        )
        return

    if contig.lower() in {name.lower() for name in contig_names}:
        return

    message = (
        f"Legacy fallback region '{region}' (config key '{region_key}') names contig '{contig}', "
        f"which is not in {bam_file}. Available contigs: {', '.join(contig_names[:10])}. "
        f"Slicing this region would return zero reads and produce a false negative. "
        f"Re-run with a --reference-assembly spelling matching the BAM's chromosome naming "
        f"(e.g. 'GRCh38' rather than 'hg38' for NCBI accessions)."
    )
    logger.error(message)
    raise ValueError(message)
