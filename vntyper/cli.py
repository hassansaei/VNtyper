import argparse
import logging
from pathlib import Path
from vntyper.scripts.pipeline import run_pipeline
from vntyper.scripts.fastq_bam_processing import process_fastq, process_bam_to_fastq
from vntyper.scripts.kestrel_genotyping import run_kestrel
from vntyper.scripts.advntr_genotyping import run_advntr, process_advntr_output
from vntyper.scripts.utils import load_config, setup_logging
from vntyper.version import __version__ as VERSION
import sys

def main():
    parser = argparse.ArgumentParser(description="VNtyper CLI: A pipeline for genotyping MUC1-VNTR.", add_help=False)
    
    default_config_path = Path(__file__).parent / "config.json"

    # Adding global flags with short and long options
    parser.add_argument('-h', '--help', action='store_true', help="Show this help message and exit")
    parser.add_argument('-l', '--log-level', help="Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR)", default="INFO")
    parser.add_argument('-f', '--log-file', help="Set the log output file (default is stdout)", default=None)
    parser.add_argument('--config-path', type=Path, help="Path to the config.json file", default=default_config_path)
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {VERSION}')

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subcommand for running the full pipeline
    parser_pipeline = subparsers.add_parser("pipeline", help="Run the full VNtyper pipeline.")
    parser_pipeline.add_argument('--advntr-reference', type=str, choices=["hg19", "hg38"], required=False, help="Override reference assembly for adVNTR genotyping (hg19 or hg38).")
    parser_pipeline.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for the results.")
    parser_pipeline.add_argument('--ignore-advntr', action='store_true', help="Skip adVNTR genotyping of MUC1-VNTR.")
    parser_pipeline.add_argument('--fastq1', type=str, help="Path to the first FASTQ file.")
    parser_pipeline.add_argument('--fastq2', type=str, help="Path to the second FASTQ file.")
    parser_pipeline.add_argument('--bam', type=str, help="Path to the BAM file.")
    parser_pipeline.add_argument('--threads', type=int, default=4, help="Number of threads to use.")
    parser_pipeline.add_argument('--reference-assembly', type=str, choices=["hg19", "hg38"], default="hg19",
                                help="Specify the reference assembly used for the input BAM file alignment.")
    parser_pipeline.add_argument('--fast-mode', action='store_true', help="Enable fast mode (skips filtering for unmapped and partially mapped reads).")
    parser_pipeline.add_argument('--keep-intermediates', action='store_true', help="Keep intermediate files (e.g., BAM slices, temporary files).")
    parser_pipeline.add_argument('--delete-intermediates', action='store_true', help="Delete intermediate files after processing (overrides --keep-intermediates).")

    # Subcommand for FASTQ processing
    parser_fastq = subparsers.add_parser("fastq", help="Process FASTQ files.")
    parser_fastq.add_argument('-r1', '--fastq1', type=str, required=True, help="Path to the first FASTQ file.")
    parser_fastq.add_argument('-r2', '--fastq2', type=str, help="Path to the second FASTQ file.")
    parser_fastq.add_argument('-t', '--threads', type=int, default=4, help="Number of threads to use.")
    parser_fastq.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for processed FASTQ files.")
    
    # Subcommand for BAM processing
    parser_bam = subparsers.add_parser("bam", help="Process BAM files.")
    parser_bam.add_argument('-a', '--alignment', type=str, required=True, help="Path to the BAM file.")
    parser_bam.add_argument('-t', '--threads', type=int, default=4, help="Number of threads to use.")
    parser_bam.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for processed BAM files.")
    parser_bam.add_argument('--reference-assembly', type=str, choices=["hg19", "hg38"], default="hg19",
                            help="Specify the reference assembly to use (hg19 or hg38). Default is hg19.")
    parser_bam.add_argument('--fast-mode', action='store_true', help="Enable fast mode (skips filtering for unmapped and partially mapped reads).")
    parser_bam.add_argument('--keep-intermediates', action='store_true', help="Keep intermediate files (e.g., BAM slices, temporary files).")
    parser_bam.add_argument('--delete-intermediates', action='store_true', help="Delete intermediate files after processing (overrides --keep-intermediates).")

    # Subcommand for Kestrel genotyping
    parser_kestrel = subparsers.add_parser("kestrel", help="Run Kestrel genotyping.")
    parser_kestrel.add_argument('-r', '--reference-vntr', type=str, required=True, help="Path to the MUC1-specific reference VNTR file.")
    parser_kestrel.add_argument('-f1', '--fastq1', type=str, required=True, help="Path to the first FASTQ file.")
    parser_kestrel.add_argument('-f2', '--fastq2', type=str, help="Path to the second FASTQ file.")
    parser_kestrel.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for Kestrel results.")
    
    # Subcommand for adVNTR genotyping
    parser_advntr = subparsers.add_parser("advntr", help="Run adVNTR genotyping.")
    parser_advntr.add_argument('-a', '--alignment', type=str, required=True, help="Path to the BAM file.")
    parser_advntr.add_argument('-m', '--reference-vntr', type=str, required=True, help="Path to the adVNTR reference VNTR database.")
    parser_advntr.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for adVNTR results.")

    # Parse arguments
    args = parser.parse_args()

    # Load config with error handling
    config = None
    try:
        if args.config_path and args.config_path.exists():
            config = load_config(args.config_path)
        else:
            logging.error(f"Configuration file not found at {args.config_path}. Using default values where applicable.")
    except Exception as e:
        logging.critical(f"Failed to load configuration: {e}")
        sys.exit(1)

    # Determine the appropriate BWA reference and adVNTR reference based on the assembly or advntr_reference
    if args.advntr_reference:  # If --advntr-reference is provided, override the references
        if args.advntr_reference == "hg19":
            bwa_reference = config["reference_data"]["bwa_reference_hg19"]
            advntr_reference = config["reference_data"]["advntr_reference_vntr_hg19"]
        elif args.advntr_reference == "hg38":
            bwa_reference = config["reference_data"]["bwa_reference_hg38"]
            advntr_reference = config["reference_data"]["advntr_reference_vntr_hg38"]
    else:  # Use reference-assembly if --advntr-reference is not provided
        if args.reference_assembly == "hg19":
            bwa_reference = config["reference_data"]["bwa_reference_hg19"]
            advntr_reference = config["reference_data"]["advntr_reference_vntr_hg19"]
        else:
            bwa_reference = config["reference_data"]["bwa_reference_hg38"]
            advntr_reference = config["reference_data"]["advntr_reference_vntr_hg38"]

    # Display welcome message if --help is called
    if args.help:
        if config and "welcome_message" in config:
            print(config["welcome_message"])
        parser.print_help()
        sys.exit(0)

    # Setup logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    setup_logging(log_level=log_level, log_file=args.log_file)

    # Handle subcommands
    if args.command == "pipeline":
        run_pipeline(
            bwa_reference=bwa_reference,
            advntr_reference=advntr_reference,
            output_dir=Path(args.output_dir),
            ignore_advntr=args.ignore_advntr,
            config=config,
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            bam=args.bam,
            threads=args.threads,
            reference_assembly=args.reference_assembly,
            fast_mode=args.fast_mode,
            keep_intermediates=args.keep_intermediates,
            delete_intermediates=args.delete_intermediates
        )
    
    elif args.command == "fastq":
        process_fastq(
            fastq1=args.fastq1, 
            fastq2=args.fastq2, 
            threads=args.threads, 
            output_dir=Path(args.output_dir), 
            config=config
        )
    
    elif args.command == "bam":
        process_bam_to_fastq(
            bam=args.alignment, 
            output_dir=Path(args.output_dir), 
            threads=args.threads, 
            config=config, 
            reference_assembly=args.reference_assembly, 
            fast_mode=args.fast_mode, 
            delete_intermediates=args.delete_intermediates, 
            keep_intermediates=args.keep_intermediates
        )
    
    elif args.command == "kestrel":
        run_kestrel(
            vcf_path=Path(args.output_dir) / "output.vcf",
            output_dir=Path(args.output_dir),
            fastq_1=args.fastq1,
            fastq_2=args.fastq2,
            reference_vntr=args.reference_vntr,
            kestrel_path=config["tools"]["kestrel"],
            kestrel_settings=config["kestrel_settings"]
        )
    
    elif args.command == "advntr":
        run_advntr(
            db_file_hg19=args.reference_vntr,
            sorted_bam=args.alignment,
            output=args.output_dir,
            output_name="output",
            config=config
        )

        # After running adVNTR, process the output
        advntr_vcf_path = Path(args.output_dir) / "output_adVNTR.vcf"
        process_advntr_output(
            vcf_path=advntr_vcf_path,
            output=args.output_dir,
            output_name="output",
            config=config
        )

if __name__ == "__main__":
    main()
