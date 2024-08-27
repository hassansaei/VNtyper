import subprocess as sp
import logging

def align_and_sort_fastq(fastq1, fastq2, reference, output_dir, output_name, threads, config):
    """
    Align FASTQ files to the reference genome using BWA, sort the SAM file, and convert to BAM using Samtools.

    Args:
        fastq1: Path to the first FASTQ file.
        fastq2: Path to the second FASTQ file.
        reference: Path to the reference genome in FASTA format.
        output_dir: Directory where output files will be saved.
        output_name: Base name for the output files.
        threads: Number of threads to use.
        config: Configuration dictionary with paths and parameters.
    """
    # Extract paths from the config
    samtools_path = config["tools"]["samtools"]

    sam_out = f"{output_dir}{output_name}.sam"
    bam_out = f"{output_dir}{output_name}.bam"
    sorted_bam_out = f"{output_dir}{output_name}_sorted.bam"
    
    # BWA MEM command
    bwa_command = f"bwa mem -t {threads} {reference} {fastq1} {fastq2} -o {sam_out}"
    
    # Samtools view and sort commands to convert SAM to BAM and sort the BAM file
    samtools_view_command = f"{samtools_path} view -@ {threads} -bS {sam_out} -o {bam_out}"
    samtools_sort_command = f"{samtools_path} sort -@ {threads} -o {sorted_bam_out} {bam_out}"
    
    # Samtools index command
    samtools_index_command = f"{samtools_path} index {sorted_bam_out}"
    
    logging.info("Starting BWA alignment...")
    process = sp.Popen(bwa_command, shell=True)
    process.wait()
    logging.info("BWA alignment completed.")

    logging.info("Converting SAM to BAM with Samtools...")
    process = sp.Popen(samtools_view_command, shell=True)
    process.wait()
    logging.info("SAM to BAM conversion completed.")

    logging.info("Sorting BAM file with Samtools...")
    process = sp.Popen(samtools_sort_command, shell=True)
    process.wait()
    logging.info("BAM sorting completed.")

    logging.info("Indexing sorted BAM file with Samtools...")
    process = sp.Popen(samtools_index_command, shell=True)
    process.wait()
    logging.info("Samtools indexing completed.")
    
    return sorted_bam_out
