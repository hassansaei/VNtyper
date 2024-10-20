# vntyper/scripts/alignment_processing.py

import logging
from pathlib import Path

from vntyper.scripts.utils import run_command


def check_bwa_index(reference):
    """
    Check if the BWA index files exist for the given reference genome.
    The index files should have extensions: .amb, .ann, .bwt, .pac, and .sa.

    Args:
        reference (str or Path): Path to the reference genome (without extension).

    Returns:
        bool: True if all BWA index files exist, False otherwise.
    """
    required_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    reference = Path(reference)
    missing_files = [reference.with_suffix(ext) for ext in required_extensions if not reference.with_suffix(ext).exists()]
    
    if missing_files:
        logging.warning(f"Missing BWA index files for reference {reference}: {[str(f) for f in missing_files]}")
        return False
    return True


def align_and_sort_fastq(
    fastq1,
    fastq2,
    reference,
    output_dir,
    output_name,
    threads,
    config
):
    """
    Align FASTQ files to the reference genome using BWA, sort, and convert to BAM directly using Samtools.

    Args:
        fastq1 (str or Path): Path to the first FASTQ file.
        fastq2 (str or Path): Path to the second FASTQ file.
        reference (str or Path): Path to the reference genome in FASTA format.
        output_dir (str or Path): Directory where output files will be saved.
        output_name (str): Base name for the output files.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary with paths and parameters.

    Returns:
        str or None: Path to the sorted BAM file, or None if the process failed.
    """
    samtools_path = Path(config["tools"]["samtools"])
    bwa_path = Path(config["tools"]["bwa"])
    
    reference = Path(reference)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    sorted_bam_out = output_dir / f"{output_name}_sorted.bam"

    # Check if the BWA index files exist
    if not check_bwa_index(reference):
        logging.error(
            f"BWA index files not found for reference: {reference}. "
            f"Please run 'bwa index {reference}' to generate them."
        )
        return None

    # Construct BWA MEM command
    bwa_command = (
        f"{bwa_path} mem -t {threads} {reference} {fastq1} {fastq2}"
    )
    # Construct Samtools view and sort command
    samtools_view_sort_command = (
        f"{samtools_path} view -@ {threads} -b | "
        f"{samtools_path} sort -@ {threads} -o {sorted_bam_out}"
    )
    
    # Combine commands using pipes
    full_command = f"{bwa_command} | {samtools_view_sort_command}"
    
    log_file_alignment = output_dir / f"{output_name}_alignment.log"
    logging.info(f"Executing alignment and sorting with command: {full_command}")
    
    # Execute the alignment and sorting command
    if not run_command(str(full_command), str(log_file_alignment), critical=True):
        logging.error("BWA alignment and Samtools sorting failed.")
        return None

    if not sorted_bam_out.exists():
        logging.error(
            f"Sorted BAM file {sorted_bam_out} not created. "
            f"BWA alignment or Samtools sorting might have failed."
        )
        return None

    logging.info("BWA alignment and Samtools sorting completed successfully.")

    # Index the sorted BAM file
    logging.info(f"Indexing sorted BAM file: {sorted_bam_out}")
    samtools_index_command = f"{samtools_path} index {sorted_bam_out}"
    
    log_file_index = output_dir / f"{output_name}_index.log"
    if not run_command(str(samtools_index_command), str(log_file_index), critical=True):
        logging.error("Samtools indexing failed.")
        return None

    index_file = sorted_bam_out.with_suffix(".bam.bai")
    if not index_file.exists():
        logging.error(
            f"BAM index file {index_file} not created. "
            f"Samtools indexing might have failed."
        )
        return None

    logging.info("Samtools indexing completed successfully.")
    return str(sorted_bam_out)
