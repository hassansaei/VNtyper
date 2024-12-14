# vntyper/scripts/alignment_processing.py

import logging
from pathlib import Path
from typing import Optional

from vntyper.scripts.utils import run_command


def check_bwa_index(reference: Path) -> bool:
    """
    Check if the BWA index files exist for the given reference genome.

    This function verifies the presence of all required BWA index files
    associated with the provided reference genome. The expected index files
    have the following extensions: .amb, .ann, .bwt, .pac, and .sa.

    Args:
        reference (Path): Path to the reference genome (without extension).

    Returns:
        bool: True if all BWA index files exist, False otherwise.
    """
    required_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    reference = Path(reference)
    missing_files = [
        reference.with_suffix(ext)
        for ext in required_extensions
        if not reference.with_suffix(ext).exists()
    ]

    if missing_files:
        # Log a warning with the list of missing index files
        logging.warning(
            f"Missing BWA index files for reference {reference}: "
            f"{[str(f) for f in missing_files]}"
        )
        return False
    return True


def align_and_sort_fastq(
    fastq1: Path,
    fastq2: Path,
    reference: Path,
    output_dir: Path,
    output_name: str,
    threads: int,
    config: dict
) -> Optional[str]:
    """
    Align FASTQ files to the reference genome using BWA, then sort and convert to BAM using Samtools.

    This function performs the following steps:
    1. Checks for the existence of BWA index files for the reference genome.
    2. Aligns paired-end FASTQ files to the reference genome using BWA MEM.
    3. Pipes the alignment output to Samtools for conversion to BAM format and sorting.
    4. Indexes the resulting sorted BAM file.

    Args:
        fastq1 (Path): Path to the first FASTQ file.
        fastq2 (Path): Path to the second FASTQ file.
        reference (Path): Path to the reference genome in FASTA format.
        output_dir (Path): Directory where output files will be saved.
        output_name (str): Base name for the output files.
        threads (int): Number of threads to use for alignment and sorting.
        config (dict): Configuration dictionary with paths to tools and other parameters.

    Returns:
        Optional[str]: Path to the sorted BAM file if successful, or None if the process failed.
    """
    # Retrieve paths to Samtools and BWA executables from the configuration
    try:
        samtools_path = Path(config["tools"]["samtools"])
        bwa_path = Path(config["tools"]["bwa"])
    except KeyError as e:
        logging.error(f"Missing tool path in configuration: {e}")
        return None

    # Create the output directory if it does not exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define the path for the sorted BAM output file
    sorted_bam_out = output_dir / f"{output_name}_sorted.bam"

    # Check if the BWA index files exist for the reference genome
    if not check_bwa_index(reference):
        logging.error(
            f"BWA index files not found for reference: {reference}. "
            f"Please run 'bwa index {reference}' to generate them."
        )
        return None

    # Construct the BWA MEM command for alignment
    bwa_command = (
        f"{bwa_path} mem -t {threads} {reference} {fastq1} {fastq2}"
    )

    # Construct the Samtools command for converting SAM to BAM and sorting
    samtools_view_sort_command = (
        f"{config['tools']['samtools']} view -@ {threads} -b | "
        f"{config['tools']['samtools']} sort -@ {threads} -o {sorted_bam_out}"
    )

    # Combine BWA and Samtools commands using a pipe
    full_command = f"{bwa_command} | {samtools_view_sort_command}"

    # Define the path for the alignment log file
    log_file_alignment = output_dir / f"{output_name}_alignment.log"
    logging.info(f"Executing alignment and sorting with command: {full_command}")

    # Execute the combined BWA and Samtools command
    if not run_command(str(full_command), str(log_file_alignment), critical=True):
        logging.error("BWA alignment and Samtools sorting failed.")
        return None

    # Verify that the sorted BAM file was created successfully
    if not sorted_bam_out.exists():
        logging.error(
            f"Sorted BAM file {sorted_bam_out} not created. "
            f"BWA alignment or Samtools sorting might have failed."
        )
        return None

    logging.info("BWA alignment and Samtools sorting completed successfully.")

    # Index the sorted BAM file using Samtools
    logging.info(f"Indexing sorted BAM file: {sorted_bam_out}")
    samtools_index_command = f"{samtools_path} index {sorted_bam_out}"

    # Define the path for the indexing log file
    log_file_index = output_dir / f"{output_name}_index.log"

    # Execute the Samtools indexing command
    if not run_command(str(samtools_index_command), str(log_file_index), critical=True):
        logging.error("Samtools indexing failed.")
        return None

    # Define the expected path for the BAM index file
    index_file = sorted_bam_out.with_suffix(".bam.bai")

    # Verify that the BAM index file was created successfully
    if not index_file.exists():
        logging.error(
            f"BAM index file {index_file} not created. "
            f"Samtools indexing might have failed."
        )
        return None

    logging.info("Samtools indexing completed successfully.")
    return str(sorted_bam_out)
