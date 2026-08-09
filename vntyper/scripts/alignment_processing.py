# vntyper/scripts/alignment_processing.py

from __future__ import annotations

import logging
from pathlib import Path

from vntyper.scripts.alignment_target_io import bwa_index_paths
from vntyper.scripts.command_builders import (
    build_bwa_align_sort_command,
    build_samtools_index_command,
)
from vntyper.scripts.utils import run_command

logger = logging.getLogger(__name__)


def check_bwa_index(reference: Path, config: dict | None = None) -> bool:
    """
    Check if the BWA index files exist for the given reference genome.

    This function verifies the presence of all required BWA index files
    associated with the provided reference genome. The expected index files
    have the following extensions: .amb, .ann, .bwt, .pac, and .sa.

    Args:
        reference (Path): Path to the reference genome (without extension).
        config (dict | None): Pipeline configuration controlling BWA index suffixes.

    Returns:
        bool: True if all BWA index files exist, False otherwise.
    """
    reference = Path(reference)
    missing_files = [path for path in bwa_index_paths(reference, config or {}) if not path.exists()]

    if missing_files:
        # Log a warning with the list of missing index files
        logger.warning(f"Missing BWA index files for reference {reference}: {[str(f) for f in missing_files]}")
        return False
    return True


def align_and_sort_fastq(
    fastq1: Path,
    fastq2: Path | None,
    reference: Path,
    output_dir: Path,
    output_name: str,
    threads: int,
    config: dict,
) -> str | None:
    """
    Align FASTQ files to the reference genome using BWA, then sort and convert to BAM using Samtools.

    This function performs the following steps:
    1. Checks for the existence of BWA index files for the reference genome.
    2. Aligns one or two FASTQ files to the reference genome using BWA MEM.
    3. Pipes the alignment output to Samtools for conversion to BAM format and sorting.
    4. Indexes the resulting sorted BAM file.

    Args:
        fastq1 (Path): Path to the first FASTQ file.
        fastq2 (Path | None): Optional path to the second FASTQ file.
        reference (Path): Path to the reference genome in FASTA format.
        output_dir (Path): Directory where output files will be saved.
        output_name (str): Base name for the output files.
        threads (int): Number of threads to use for alignment and sorting.
        config (dict): Configuration dictionary with paths to tools and other parameters.

    Returns:
        Optional[str]: Path to the sorted BAM file if successful, or None if the process failed.
    """
    try:
        samtools_path = Path(config["tools"]["samtools"])
        bwa_path = Path(config["tools"]["bwa"])
    except KeyError as e:
        logger.error(f"Missing tool path in configuration: {e}")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sorted_bam_out = output_dir / f"{output_name}_sorted.bam"

    if not check_bwa_index(reference, config):
        logger.error(
            f"BWA index files not found for reference: {reference}. "
            f"Please run 'bwa index {reference}' to generate them."
        )
        return None

    full_command = build_bwa_align_sort_command(
        bwa_path=str(bwa_path),
        samtools_path=str(samtools_path),
        threads=threads,
        reference=reference,
        fastq1=fastq1,
        fastq2=fastq2,
        sorted_bam=sorted_bam_out,
    )
    log_file_alignment = output_dir / f"{output_name}_alignment.log"
    logger.info(f"Executing alignment and sorting with command: {full_command}")

    if not run_command(str(full_command), str(log_file_alignment), critical=True):
        logger.error("BWA alignment and Samtools sorting failed.")
        return None

    if not sorted_bam_out.exists():
        logger.error(
            f"Sorted BAM file {sorted_bam_out} not created. BWA alignment or Samtools sorting might have failed."
        )
        return None

    logger.info("BWA alignment and Samtools sorting completed successfully.")

    logger.info(f"Indexing sorted BAM file: {sorted_bam_out}")
    # No output_bai here, deliberately: sorted_bam_out is the BAM this function just
    # wrote inside output_dir, so samtools' default destination beside it is already
    # inside the run's output directory. `output_bai` exists for the one caller that
    # indexes the user's *input* alignment, whose directory is routinely mounted
    # read-only (#162, #210). The index still receives the caller's configured
    # thread count like every other samtools stage.
    samtools_index_command = build_samtools_index_command(
        samtools_path=str(samtools_path),
        bam_file=sorted_bam_out,
        threads=threads,
    )
    log_file_index = output_dir / f"{output_name}_index.log"

    if not run_command(str(samtools_index_command), str(log_file_index), critical=True):
        logger.error("Samtools indexing failed.")
        return None

    index_file = sorted_bam_out.with_suffix(".bam.bai")
    if not index_file.exists():
        logger.error(f"BAM index file {index_file} not created. Samtools indexing might have failed.")
        return None

    logger.info("Samtools indexing completed successfully.")
    return str(sorted_bam_out)
