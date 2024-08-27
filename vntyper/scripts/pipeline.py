#!/usr/bin/env python3

import timeit
import os
import sys
from pathlib import Path

from vntyper.scripts.utils import setup_logging, create_output_directories
from vntyper.scripts.file_processing import read_vcf, filter_vcf, filter_indel_vcf
from vntyper.scripts.fastq_bam_processing import process_fastq, process_bam_to_fastq
from vntyper.scripts.kestrel_genotyping import run_kestrel, process_kmer
from vntyper.scripts.motif_processing import process_motifs, preprocessing_insertion, preprocessing_deletion
from vntyper.scripts.advntr_genotyping import run_advntr, process_advntr_output

def run_pipeline(reference_file, output_dir, ignore_advntr, fastq1=None, fastq2=None, bam=None, threads=4):
    """
    Main pipeline function that orchestrates the genotyping process.
    
    Args:
        reference_file: Path to the reference FASTA file.
        output_dir: Path to the output directory.
        ignore_advntr: Boolean indicating whether to skip adVNTR genotyping.
        fastq1: Path to the first FASTQ file.
        fastq2: Path to the second FASTQ file.
        bam: Path to the BAM file.
        threads: Number of threads to use.
    """
    # Create output directories
    output_dir, temp_dir = create_output_directories(output_dir, "temp")

    # Setup logging
    log_file = os.path.join(temp_dir, f"pipeline.log")
    setup_logging(log_file)

    # Print welcome message
    welcome_message = """
    ==========================================================================================================
    Given alignment (BAM) or raw file (FASTQ), this tool genotypes MUC1 coding-VNTR 
    -- For rapid genotyping, BAM files are preferred!
    -- User can Skip code-adVNTR genotyping step using --ignore_advntr option (This step will take a while..)
    v. 1.0.0
    This is free non-commercial software. 
    ==========================================================================================================
    """
    print(welcome_message)

    start = timeit.default_timer()

    try:
        # FASTQ Quality Control or BAM Processing
        if fastq1 and fastq2:
            process_fastq(fastq1, fastq2, threads, output_dir, "output")
        elif bam:
            process_bam_to_fastq(bam, output_dir, "output", threads)
        
        # Kestrel Genotyping
        vcf_out = os.path.join(output_dir, "output.vcf")
        vcf_path = Path(vcf_out)
        run_kestrel(vcf_path, output_dir, fastq1, fastq2, reference_file)
        
        # Motif and VNTR Processing (Add processing logic if needed)
        # This is where motif processing would go if applicable
        
        # adVNTR Genotyping if not skipped
        if not ignore_advntr:
            run_advntr(reference_file, reference_file, bam, output_dir, "output")

        # Final processing and output generation
        # This includes merging results from Kestrel and adVNTR genotyping

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

    # Print end message
    end_message = """
    ==============================
    Thanks for using VNtyper pipeline!
    Contact: hassan.saei@inserm.fr
    ==============================
    """
    print(end_message)

    stop = timeit.default_timer()
    elapsed_time = (stop - start) / 60
    print(f"Pipeline completed in {elapsed_time:.2f} minutes.")
