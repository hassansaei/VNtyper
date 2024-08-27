#!/usr/bin/env python3

import timeit
import os
import sys
import logging
from pathlib import Path

from vntyper.scripts.utils import setup_logging, create_output_directories, load_config
from vntyper.scripts.file_processing import read_vcf, filter_vcf, filter_indel_vcf
from vntyper.scripts.fastq_bam_processing import process_fastq, process_bam_to_fastq
from vntyper.scripts.kestrel_genotyping import run_kestrel
from vntyper.scripts.motif_processing import process_motifs, preprocessing_insertion, preprocessing_deletion
from vntyper.scripts.advntr_genotyping import run_advntr, process_advntr_output
from vntyper.scripts.alignment_processing import align_and_sort_fastq

def run_pipeline(reference_file, output_dir, ignore_advntr, config, fastq1=None, fastq2=None, bam=None, threads=4, reference_assembly="hg19", fast_mode=False, keep_intermediates=False, delete_intermediates=False):
    """
    Main pipeline function that orchestrates the genotyping process.
    
    Args:
        reference_file: Path to the reference FASTA file.
        output_dir: Path to the output directory.
        ignore_advntr: Boolean indicating whether to skip adVNTR genotyping.
        config: Configuration dictionary.
        fastq1: Path to the first FASTQ file.
        fastq2: Path to the second FASTQ file.
        bam: Path to the BAM file.
        threads: Number of threads to use.
        reference_assembly: Reference assembly to use ("hg19" or "hg38").
        fast_mode: Boolean indicating whether to enable fast mode (skip filtering of unmapped and partially mapped reads).
        keep_intermediates: Boolean indicating whether to keep intermediate files.
        delete_intermediates: Boolean indicating whether to delete intermediate files after processing.
    """
    # Create output directories
    output_dir, temp_dir = create_output_directories(output_dir, config["temp_directory"])

    # Setup logging
    log_file = os.path.join(temp_dir, f"pipeline.log")
    setup_logging(log_file)

    start = timeit.default_timer()

    try:
        # Select the appropriate BAM region and BWA reference file based on the reference assembly
        if reference_assembly == "hg38":
            bam_region = config["bam_processing"]["bam_region_hg38"]
            bwa_reference = config["reference_data"]["bwa_reference_hg38"]
        else:
            bam_region = config["bam_processing"]["bam_region_hg19"]
            bwa_reference = config["reference_data"]["bwa_reference_hg19"]

        # Determine if intermediates should be deleted (delete_intermediates takes precedence over keep_intermediates)
        delete_intermediates = delete_intermediates or not keep_intermediates

        # FASTQ Quality Control or BAM Processing
        if fastq1 and fastq2:
            process_fastq(fastq1, fastq2, threads, output_dir, "output", config)
        elif bam:
            # Convert BAM to FASTQ
            fastq1, fastq2 = process_bam_to_fastq(
                bam, output_dir, "output", threads, config, reference_assembly, fast_mode, delete_intermediates, keep_intermediates
            )

            if not fastq1 or not fastq2:
                raise ValueError("Failed to generate FASTQ files from BAM. Exiting pipeline.")

        # Alignment step
        sorted_bam = None
        if fastq1 and fastq2:
            sorted_bam = align_and_sort_fastq(fastq1, fastq2, bwa_reference, output_dir, "output", threads, config)
        elif bam:
            sorted_bam = bam
        
        # Kestrel Genotyping
        vcf_out = os.path.join(output_dir, "output.vcf")
        vcf_path = Path(vcf_out)
        reference_vntr = config["reference_data"]["muc1_reference_vntr"]  # Extract reference VNTR from config
        kestrel_path = config["tools"]["kestrel"]  # Extract Kestrel tool path from config
        kestrel_settings = config["kestrel_settings"]  # Extract Kestrel-specific settings from config
        
        if fastq1 and fastq2:
            run_kestrel(vcf_path, output_dir, fastq1, fastq2, reference_vntr, kestrel_path, temp_dir, kestrel_settings)  # Pass necessary arguments to run_kestrel
        else:
            raise ValueError("FASTQ files are required for Kestrel genotyping, but none were provided or generated.")

        # Motif and VNTR Processing (Add processing logic if needed)
        
        # adVNTR Genotyping if not skipped
        if not ignore_advntr and sorted_bam:
            run_advntr(reference_file, config["reference_data"]["advntr_reference_vntr"], sorted_bam, output_dir, "output")

        # Final processing and output generation
        # This includes merging results from Kestrel and adVNTR genotyping

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)

    # Log end message
    end_message = """
    ==============================
    Thanks for using VNtyper pipeline!
    Contact: hassan.saei@inserm.fr
    ==============================
    """
    logging.info(end_message)

    stop = timeit.default_timer()
    elapsed_time = (stop - start) / 60
    logging.info(f"Pipeline completed in {elapsed_time:.2f} minutes.")
