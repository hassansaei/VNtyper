#!/usr/bin/env python3

import timeit
import os
import sys
from pathlib import Path

from arg_parser import parse_arguments
from utils import setup_logging, create_output_directories
from file_processing import read_vcf, filter_vcf, filter_indel_vcf
from fastq_bam_processing import process_fastq, process_bam_to_fastq
from kestrel_genotyping import run_kestrel, process_kmer
from motif_processing import process_motifs, preprocessing_insertion, preprocessing_deletion
from advntr_genotyping import run_advntr, process_advntr_output

def main():
    # Parse arguments
    args = parse_arguments()

    # Create output directories
    output_dir, temp_dir = create_output_directories(args.working_dir, args.output)

    # Setup logging
    log_file = os.path.join(temp_dir, f"{args.output}.log")
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
        if args.fastq:
            process_fastq(args.fastq1, args.fastq2, args.threads, output_dir, args.output)
        elif args.bam:
            process_bam_to_fastq(args.alignment, output_dir, args.output, args.threads)
        
        # Kestrel Genotyping
        vcf_out = os.path.join(output_dir, f"{args.output}.vcf")
        vcf_path = Path(vcf_out)
        run_kestrel(vcf_path, output_dir, args.fastq1, args.fastq2, args.reference_VNTR)
        
        # Motif and VNTR Processing
        # Here you would add motif processing logic as needed

        # adVNTR Genotyping if not skipped
        if not args.ignore_advntr:
            run_advntr(args.reference_file, args.reference_vntr, args.alignment, output_dir, args.output)

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

if __name__ == "__main__":
    main()
