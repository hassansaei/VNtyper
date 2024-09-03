import subprocess as sp
import logging
import os
from vntyper.scripts.utils import run_command


def check_bwa_index(reference):
    """
    Check if the BWA index files exist for the given reference genome.
    The index files should have extensions: .amb, .ann, .bwt, .pac, and .sa.

    Args:
        reference: Path to the reference genome (without extension).

    Returns:
        bool: True if all BWA index files exist, False otherwise.
    """
    required_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    for ext in required_extensions:
        if not os.path.exists(f"{reference}{ext}"):
            return False
    return True


def align_and_sort_fastq(fastq1, fastq2, reference, output_dir, output_name, threads, config):
    """
    Align FASTQ files to the reference genome using BWA, sort, and convert to BAM directly using Samtools.

    Args:
        fastq1: Path to the first FASTQ file.
        fastq2: Path to the second FASTQ file.
        reference: Path to the reference genome in FASTA format.
        output_dir: Directory where output files will be saved.
        output_name: Base name for the output files.
        threads: Number of threads to use.
        config: Configuration dictionary with paths and parameters.

    Returns:
        str: Path to the sorted BAM file, or None if the process failed.
    """
    # Extract paths from the config
    samtools_path = config["tools"]["samtools"]

    sorted_bam_out = os.path.join(output_dir, f"{output_name}_sorted.bam")

    # Check if the BWA index files exist
    if not check_bwa_index(reference):
        logging.error(f"BWA index files not found for reference: {reference}. Please run 'bwa index' on the reference file.")
        return None

    # BWA MEM command piped directly into Samtools for sorting and BAM conversion
    bwa_command = f"bwa mem -t {threads} {reference} {fastq1} {fastq2}"
    samtools_view_sort_command = f"{samtools_path} view -@ {threads} -b | {samtools_path} sort -@ {threads} -o {sorted_bam_out}"

    log_file = os.path.join(output_dir, f"{output_name}_alignment.log")
    if not run_command(f"{bwa_command} | {samtools_view_sort_command}", log_file, critical=True):
        logging.error("BWA alignment and Samtools sorting failed.")
        return None

    if not os.path.exists(sorted_bam_out):
        logging.error(f"Sorted BAM file {sorted_bam_out} not created. BWA alignment or Samtools sorting might have failed.")
        return None

    logging.info("BWA alignment and Samtools sorting completed.")

    # Index the sorted BAM file
    logging.info("Indexing sorted BAM file with Samtools...")
    samtools_index_command = f"{samtools_path} index {sorted_bam_out}"

    log_file = os.path.join(output_dir, f"{output_name}_index.log")
    if not run_command(samtools_index_command, log_file, critical=True):
        logging.error("BAM indexing failed.")
        return None

    if not os.path.exists(f"{sorted_bam_out}.bai"):
        logging.error(f"BAM index file {sorted_bam_out}.bai not created. BAM indexing might have failed.")
        return None

    logging.info("Samtools indexing completed.")

    return sorted_bam_out