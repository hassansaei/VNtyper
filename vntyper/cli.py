import argparse
import logging
from pathlib import Path
from vntyper.scripts.pipeline import run_pipeline
from vntyper.scripts.fastq_bam_processing import process_fastq, process_bam_to_fastq
from vntyper.scripts.kestrel_genotyping import run_kestrel
from vntyper.scripts.advntr_genotyping import run_advntr
from vntyper.scripts.utils import load_config, setup_logging
from vntyper.version import __version__ as VERSION

def main():
    parser = argparse.ArgumentParser(description="VNtyper CLI: A pipeline for genotyping MUC1-VNTR.")
    
    default_config_path = Path(__file__).parent / "config.json"

    # Adding global flags with short and long options
    parser.add_argument('-l', '--log-level', help="Set the logging level", default="INFO")
    parser.add_argument('-f', '--log-file', help="Set the log output file", default=None)
    parser.add_argument('--config-path', type=Path, help="Path to the config.json file", default=default_config_path)
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {VERSION}')

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subcommand for running the full pipeline
    parser_pipeline = subparsers.add_parser("pipeline", help="Run the full VNtyper pipeline.")
    parser_pipeline.add_argument('-r', '--reference-file', type=str, required=False, help="The reference FASTA file.")
    parser_pipeline.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for the results.")
    parser_pipeline.add_argument('--ignore-advntr', action='store_true', help="Skip adVNTR genotyping of MUC1-VNTR.")
    parser_pipeline.add_argument('--fastq1', type=str, help="Path to the first FASTQ file.")
    parser_pipeline.add_argument('--fastq2', type=str, help="Path to the second FASTQ file.")
    parser_pipeline.add_argument('--bam', type=str, help="Path to the BAM file.")
    parser_pipeline.add_argument('--threads', type=int, default=4, help="Number of threads to use.")

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
    
    # Subcommand for Kestrel genotyping
    parser_kestrel = subparsers.add_parser("kestrel", help="Run Kestrel genotyping.")
    parser_kestrel.add_argument('-r', '--reference-vntr', type=str, required=True, help="Path to the MUC1-specific reference VNTR file.")
    parser_kestrel.add_argument('-f1', '--fastq1', type=str, required=True, help="Path to the first FASTQ file.")
    parser_kestrel.add_argument('-f2', '--fastq2', type=str, help="Path to the second FASTQ file.")
    parser_kestrel.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for Kestrel results.")
    
    # Subcommand for adVNTR genotyping
    parser_advntr = subparsers.add_parser("advntr", help="Run adVNTR genotyping.")
    parser_advntr.add_argument('-a', '--alignment', type=str, required=True, help="Path to the BAM file.")
    parser_advntr.add_argument('-r', '--reference-file', type=str, required=True, help="Path to the reference FASTA file.")
    parser_advntr.add_argument('-m', '--reference-vntr', type=str, required=True, help="Path to the adVNTR reference VNTR database.")
    parser_advntr.add_argument('-o', '--output-dir', type=str, default="out", help="Output directory for adVNTR results.")

    # Parse arguments
    args = parser.parse_args()

    # Setup logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    setup_logging(log_level=log_level, log_file=args.log_file)

    # Load config with error handling
    try:
        if args.config_path and args.config_path.exists():
            config = load_config(args.config_path)
        else:
            logging.error(f"Configuration file not found at {args.config_path}. Using default values where applicable.")
            config = None
    except Exception as e:
        logging.critical(f"Failed to load configuration: {e}")
        sys.exit(1)

    # Handle subcommands
    if args.command == "pipeline":
        run_pipeline(
            reference_file=args.reference_file,
            output_dir=Path(args.output_dir),
            ignore_advntr=args.ignore_advntr,
            config=config,
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            bam=args.bam,
            threads=args.threads
        )
    
    elif args.command == "fastq":
        process_fastq(args.fastq1, args.fastq2, args.threads, Path(args.output_dir), config=config)
    
    elif args.command == "bam":
        process_bam_to_fastq(args.alignment, Path(args.output_dir), args.threads, config=config)
    
    elif args.command == "kestrel":
        run_kestrel(args.reference_vntr, args.fastq1, args.fastq2, Path(args.output_dir))
    
    elif args.command == "advntr":
        run_advntr(args.alignment, args.reference_file, args.reference_vntr, Path(args.output_dir))

if __name__ == "__main__":
    main()
