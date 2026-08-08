"""Prepare pipeline alignment targets and exact preflight requests."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TypedDict

from vntyper.scripts.alignment_target_io import install_generated_bed
from vntyper.scripts.fastq_bam_processing import parse_contigs_from_header
from vntyper.scripts.region_utils import get_region_string_with_fallback

logger = logging.getLogger(__name__)


class AlignmentPreflightKwargs(TypedDict):
    """Keyword arguments shared by input and post-alignment preflight calls."""

    in_path: str
    output_dir: str
    output_name: str
    file_format: str
    config: dict
    threads: int
    region: None
    bed_file: Path
    reference_assembly: str
    reference_fasta: None
    header_contigs: tuple[str, ...]
    m5: str | None
    fast_mode: bool


def format_regions_as_bed(regions: str) -> str:
    """Convert comma-separated ``chrom:start-end`` regions to BED rows.

    Args:
        regions: Comma-separated target regions.

    Returns:
        BED-formatted text with one target per line.

    Raises:
        ValueError: If any region does not have the expected shape.
    """
    rows: list[str] = []
    for region in regions.split(","):
        try:
            chrom, positions = region.strip().split(":")
            start, end = positions.strip().split("-")
        except ValueError as error:
            message = f"Invalid region format: {region}. Expected format 'chr:start-end'."
            logger.error(message)
            raise ValueError(message) from error
        rows.append(f"{chrom}\t{start}\t{end}\n")
    return "".join(rows)


def _first_active_bed_contig(bed_text: str) -> str | None:
    for line in bed_text.splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            return stripped.split(maxsplit=1)[0]
    return None


def _target_m5_from_header(header: str, target_contig: str | None) -> str | None:
    if target_contig is None:
        return None
    for line in header.splitlines():
        fields = line.split("\t")
        if not fields or fields[0] != "@SQ":
            continue
        tags: dict[str, str] = {}
        for field in fields[1:]:
            key, separator, value = field.partition(":")
            if separator:
                tags[key] = value
        if tags.get("SN") == target_contig:
            return tags.get("M5")
    return None


def prepare_alignment_target(
    *,
    input_type: str,
    bam: str | None,
    cram: str | None,
    output_dir: str | Path,
    reference_assembly: str,
    config: dict,
    bed_file: str | Path | None,
    custom_regions: str | None,
) -> Path:
    """Resolve and materialize the exact BED target later given to preflight.

    Args:
        input_type: Pipeline input kind (``BAM``, ``CRAM``, or ``FASTQ``).
        bam: BAM input path when selected.
        cram: CRAM input path when selected.
        output_dir: Pipeline output root for generated BED files.
        reference_assembly: Declared input assembly.
        config: Pipeline configuration.
        bed_file: Operator-provided BED path, when present.
        custom_regions: Operator-provided comma-separated regions, when present.

    Returns:
        The existing or generated target BED path.

    Raises:
        FileNotFoundError: If an operator-provided BED does not exist.
        ValueError: If a generated target is absent or malformed.
    """
    if bed_file:
        bed_file_path = Path(bed_file)
        if not bed_file_path.exists():
            logger.error(f"Provided BED file does not exist: {bed_file_path}")
            raise FileNotFoundError(f"BED file not found: {bed_file_path}")
        logger.info(f"Using provided BED file: {bed_file_path}")
        return bed_file_path

    if custom_regions:
        bed_file_path = Path(output_dir) / "custom_regions.bed"
        regions = custom_regions
        description = "Custom regions"
    else:
        if input_type in {"BAM", "CRAM"}:
            input_file = bam if input_type == "BAM" else cram
            if input_file is None:
                message = f"{input_type} input path is required for target resolution."
                logger.error(message)
                raise ValueError(message)
            regions = get_region_string_with_fallback(
                bam_file=input_file,
                reference_assembly=reference_assembly,
                region_type="bam_region",
                config=config,
            )
        else:
            region_key = f"bam_region_{reference_assembly}"
            regions = config.get("bam_processing", {}).get(region_key)
            if not regions:
                message = f"Missing configuration for region: {region_key}"
                logger.error(f"Region key '{region_key}' not found in config.json under 'bam_processing'.")
                raise ValueError(message)
        bed_file_path = Path(output_dir) / f"predefined_regions_{reference_assembly}.bed"
        description = "Predefined regions"

    input_path = bam if input_type == "BAM" else cram if input_type == "CRAM" else None
    file_format = input_type.lower() if input_path is not None else None
    install_generated_bed(
        bed_file_path,
        format_regions_as_bed(regions),
        input_path=input_path,
        file_format=file_format,
    )
    logger.info(f"{description} converted to BED file: {bed_file_path}")
    return bed_file_path


def build_alignment_preflight_kwargs(
    *,
    in_path: str | Path,
    output_dir: str | Path,
    output_name: str,
    file_format: str,
    config: dict,
    threads: int,
    bed_file: Path,
    reference_assembly: str,
    fast_mode: bool,
    alignment_header: str | None = None,
) -> AlignmentPreflightKwargs:
    """Build the exact, BED-shaped request used by a real alignment preflight.

    Args:
        in_path: Alignment path to prove.
        output_dir: Run-local alignment working directory.
        output_name: Basename for the view and preflight logs.
        file_format: ``bam`` or ``cram``.
        config: Pipeline configuration.
        threads: Samtools thread count.
        bed_file: Exact BED later used for the real slice.
        reference_assembly: Declared input assembly.
        fast_mode: Whether downstream processing permits any htslib index.
        alignment_header: Header already read by the assembly guard, if available.

    Returns:
        Typed keyword arguments for :func:`alignment_preflight.run_preflight`.
    """
    header_contigs = (
        tuple(str(contig["name"]) for contig in parse_contigs_from_header(alignment_header) if "name" in contig)
        if alignment_header is not None
        else ()
    )
    target_contig = _first_active_bed_contig(bed_file.read_text(encoding="utf-8"))
    target_m5 = _target_m5_from_header(alignment_header, target_contig) if alignment_header is not None else None
    return {
        "in_path": str(in_path),
        "output_dir": str(output_dir),
        "output_name": output_name,
        "file_format": file_format,
        "config": config,
        "threads": threads,
        "region": None,
        "bed_file": bed_file,
        "reference_assembly": reference_assembly,
        "reference_fasta": None,
        "header_contigs": header_contigs,
        "m5": target_m5,
        "fast_mode": fast_mode,
    }
