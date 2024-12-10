#!/usr/bin/env python3
# vntyper/cli.py
# VNtyper CLI entry point

import argparse
import importlib.resources as pkg_resources  # For accessing package data
import json
import logging
import sys
from pathlib import Path

from vntyper.scripts.cohort_summary import aggregate_cohort
from vntyper.scripts.fastq_bam_processing import (
    process_bam_to_fastq, process_fastq
)
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.install_references import main as install_references_main
from vntyper.scripts.kestrel_genotyping import run_kestrel
from vntyper.scripts.pipeline import run_pipeline
from vntyper.scripts.utils import setup_logging
from vntyper.version import __version__ as VERSION


def load_config(config_path=None):
    """
    Load the configuration file with fallback to the default package config.

    Args:
        config_path (Path or None): Path to the user-provided config file.

    Returns:
        dict: The loaded configuration dictionary.
    """
    if config_path is not None and Path(config_path).exists():
        # User provided a config path
        with open(config_path, 'r') as f:
            config = json.load(f)
    else:
        # No config path provided or file does not exist; use default config from package data
        try:
            with pkg_resources.open_text('vntyper', 'config.json') as f:
                config = json.load(f)
        except Exception as e:
            logging.error("Error: Default config file not found in package data.")
            logging.error(e)
            sys.exit(1)
    return config


def main():
    """
    Main function to parse arguments and execute corresponding subcommands.
    """
    # Create the main parser and add global arguments
    parser = argparse.ArgumentParser(
        description="VNtyper CLI: A pipeline for genotyping MUC1-VNTR.",
        add_help=True
    )

    # Adding global flags with short and long options
    parser.add_argument(
        '-l', '--log-level',
        help="Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR)",
        default="INFO"
    )
    parser.add_argument(
        '-f', '--log-file',
        help="Set the log output file (default is stdout)",
        default=None
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'%(prog)s {VERSION}'
    )
    parser.add_argument(
        '--config-path',
        type=Path,
        default=None,
        help="Path to the configuration file (config.json). If not provided, the default config will be used.",
        required=False
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subcommand for running the full pipeline
    parser_pipeline = subparsers.add_parser(
        "pipeline",
        help="Run the full VNtyper pipeline."
    )
    parser_pipeline.add_argument(
        '-o', '--output-dir',
        type=str,
        default="out",
        help="Output directory for the results."
    )
    # Added the --extra-modules flag
    parser_pipeline.add_argument(
        '--extra-modules',
        nargs='*',
        default=[],
        help="Optional extra modules to include (e.g., advntr)."
    )
    parser_pipeline.add_argument(
        '--fastq1',
        type=str,
        help="Path to the first FASTQ file."
    )
    parser_pipeline.add_argument(
        '--fastq2',
        type=str,
        help="Path to the second FASTQ file."
    )
    parser_pipeline.add_argument(
        '--bam',
        type=str,
        help="Path to the BAM file."
    )
    parser_pipeline.add_argument(
        '--threads',
        type=int,
        default=4,
        help="Number of threads to use."
    )
    parser_pipeline.add_argument(
        '--reference-assembly',
        type=str,
        choices=["hg19", "hg38"],
        default="hg19",
        help="Specify the reference assembly used for the input BAM file alignment."
    )
    parser_pipeline.add_argument(
        '--fast-mode',
        action='store_true',
        help="Enable fast mode (skips filtering for unmapped and partially mapped reads)."
    )
    parser_pipeline.add_argument(
        '--keep-intermediates',
        action='store_true',
        help="Keep intermediate files (e.g., BAM slices, temporary files)."
    )
    parser_pipeline.add_argument(
        '--delete-intermediates',
        action='store_true',
        help="Delete intermediate files after processing (overrides --keep-intermediates)."
    )
    # New arguments for archiving results
    parser_pipeline.add_argument(
        '--archive-results',
        action='store_true',
        help="Create an archive of the results folder after pipeline completion."
    )
    parser_pipeline.add_argument(
        '--archive-format',
        type=str,
        choices=['zip', 'tar.gz'],
        default='zip',
        help="Format of the archive: 'zip' or 'tar.gz'. Default is 'zip'."
    )
    # Added output_name argument
    parser_pipeline.add_argument(
        '-n', '--output-name',
        type=str,
        default="processed",
        help="Base name for the output files."
    )

    # Mutually exclusive group for custom regions and BED file
    region_group = parser_pipeline.add_mutually_exclusive_group()
    region_group.add_argument(
        '--custom-regions',
        type=str,
        help=(
            "Define custom regions for MUC1 analysis as comma-separated values "
            "(e.g., chr1:1000-2000,chr2:3000-4000)."
        )
    )
    region_group.add_argument(
        '--bed-file',
        type=Path,
        help="Path to a BED file specifying regions for MUC1 analysis."
    )

    # Module-specific argument groups
    module_parsers = {}
    if 'advntr' in sys.argv:
        module_parsers['advntr'] = parser_pipeline.add_argument_group(
            'adVNTR Module Options'
        )
        module_parsers['advntr'].add_argument(
            '--advntr-reference',
            type=str,
            choices=["hg19", "hg38"],
            required=False,
            help="Reference assembly for adVNTR genotyping (hg19 or hg38)."
        )

    # Subcommand for FASTQ processing
    parser_fastq = subparsers.add_parser(
        "fastq",
        help="Process FASTQ files."
    )
    parser_fastq.add_argument(
        '-r1', '--fastq1',
        type=str,
        required=True,
        help="Path to the first FASTQ file."
    )
    parser_fastq.add_argument(
        '-r2', '--fastq2',
        type=str,
        help="Path to the second FASTQ file."
    )
    parser_fastq.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help="Number of threads to use."
    )
    parser_fastq.add_argument(
        '-o', '--output-dir',
        type=str,
        default="out",
        help="Output directory for processed FASTQ files."
    )
    parser_fastq.add_argument(
        '-n', '--output-name',
        type=str,
        default="processed",
        help="Base name for the output FASTQ files."
    )

    # Subcommand for BAM processing
    parser_bam = subparsers.add_parser(
        "bam",
        help="Process BAM files."
    )
    parser_bam.add_argument(
        '-a', '--alignment',
        type=str,
        required=True,
        help="Path to the BAM file."
    )
    parser_bam.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        help="Number of threads to use."
    )
    parser_bam.add_argument(
        '-o', '--output-dir',
        type=str,
        default="out",
        help="Output directory for processed BAM files."
    )
    parser_bam.add_argument(
        '--reference-assembly',
        type=str,
        choices=["hg19", "hg38"],
        default="hg19",
        help="Specify the reference assembly to use (hg19 or hg38). Default is hg19."
    )
    parser_bam.add_argument(
        '--fast-mode',
        action='store_true',
        help="Enable fast mode (skips filtering for unmapped and partially mapped reads)."
    )
    parser_bam.add_argument(
        '--keep-intermediates',
        action='store_true',
        help="Keep intermediate files (e.g., BAM slices, temporary files)."
    )
    parser_bam.add_argument(
        '--delete-intermediates',
        action='store_true',
        help="Delete intermediate files after processing (overrides --keep-intermediates)."
    )
    parser_bam.add_argument(
        '-n', '--output-name',
        type=str,
        default="processed",
        help="Base name for the output FASTQ files."
    )

    # Subcommand for Kestrel genotyping
    parser_kestrel = subparsers.add_parser(
        "kestrel",
        help="Run Kestrel genotyping."
    )
    parser_kestrel.add_argument(
        '-r', '--reference-vntr',
        type=str,
        required=True,
        help="Path to the MUC1-specific reference VNTR file."
    )
    parser_kestrel.add_argument(
        '-f1', '--fastq1',
        type=str,
        required=True,
        help="Path to the first FASTQ file."
    )
    parser_kestrel.add_argument(
        '-f2', '--fastq2',
        type=str,
        help="Path to the second FASTQ file."
    )
    parser_kestrel.add_argument(
        '-o', '--output-dir',
        type=str,
        default="out",
        help="Output directory for Kestrel results."
    )

    # Subcommand for generating reports
    parser_report = subparsers.add_parser(
        "report",
        help="Generate a summary report and visualizations from output data."
    )
    parser_report.add_argument(
        '-o', '--output-dir',
        type=str,
        required=True,
        help="Output directory containing pipeline results."
    )
    parser_report.add_argument(
        '--report-file',
        type=str,
        default="summary_report.html",
        help="Name of the output report file."
    )
    parser_report.add_argument(
        '--log-file',
        type=str,
        default="pipeline.log",
        help="Pipeline log file to include in the report."
    )
    # Subcommand for generating IGV reports
    parser_report.add_argument(
        '--bed-file',
        type=Path,
        help="Path to the BED file for IGV reports."
    )
    parser_report.add_argument(
        '--bam-file',
        type=Path,
        help="Path to the BAM file for IGV reports."
    )
    parser_report.add_argument(
        '--reference-fasta',
        type=Path,
        help="Path to the reference FASTA file for IGV reports."
    )
    parser_report.add_argument(
        '--flanking',
        type=int,
        default=50,
        help="Flanking region size for IGV reports."
    )

    # Subcommand for cohort analysis
    parser_cohort = subparsers.add_parser(
        "cohort",
        help="Aggregate outputs from multiple runs into a single summary file."
    )
    parser_cohort.add_argument(
        '-i', '--input-dirs',
        nargs='+',
        required=True,
        help="List of directories containing output files to aggregate."
    )
    parser_cohort.add_argument(
        '-o', '--output-dir',
        type=str,
        required=True,
        help="Output directory for the aggregated summary."
    )
    parser_cohort.add_argument(
        '--summary-file',
        type=str,
        default="cohort_summary.html",
        help="Name of the cohort summary report file."
    )

    # Subcommand for installing references
    parser_install = subparsers.add_parser(
        "install-references",
        help="Download and set up necessary reference files."
    )
    parser_install.add_argument(
        '-d', '--output-dir',
        type=Path,
        required=True,
        help="Directory where references will be installed."
    )
    parser_install.add_argument(
        '-c', '--config-path',
        type=Path,
        default=None,
        help="Path to the main config.json file to update. If not provided, config update is skipped."
    )
    parser_install.add_argument(
        '--skip-indexing',
        action='store_true',
        help="Skip the bwa indexing step."
    )

    # Parse arguments
    args = parser.parse_args()

    # Display help if no command is provided
    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # Setup logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    if args.log_file:
        log_file_path = Path(args.log_file)
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        setup_logging(log_level=log_level, log_file=str(log_file_path))
    else:
        setup_logging(log_level=log_level, log_file=None)

    # Handle install-references subcommand
    if args.command == "install-references":
        # Execute the install_references_main function with arguments
        install_references_main(
            output_dir=args.output_dir,
            config_path=args.config_path,
            skip_indexing=args.skip_indexing
        )
        sys.exit(0)

    # Load configuration
    try:
        config = load_config(args.config_path)
    except Exception as e:
        logging.critical(f"Failed to load configuration: {e}")
        sys.exit(1)

    # Handle other subcommands
    if args.command != "install-references":
        if args.log_file:
            log_file = args.log_file
        else:
            log_file = Path(args.output_dir) / "pipeline.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        setup_logging(log_level=log_level, log_file=str(log_file))

    # Execute the corresponding subcommand
    if args.command == "pipeline":
        # ----------------------------
        # Custom Validation for Pipeline Inputs
        # ----------------------------
        # Ensure that only one type of input is provided (either BAM or FASTQ)
        if args.bam and (args.fastq1 or args.fastq2):
            parser_pipeline.error("Provide either BAM or FASTQ files, not both.")

        # If BAM is not provided, ensure both FASTQ files are provided
        if not args.bam and (args.fastq1 is None or args.fastq2 is None):
            parser_pipeline.error(
                "When not providing BAM, both --fastq1 and --fastq2 must be specified "
                "for paired-end sequencing."
            )

        # Collect module-specific arguments
        module_args = {}
        if 'advntr' in args.extra_modules:
            module_args['advntr'] = {
                'advntr_reference': args.advntr_reference
            }
            # Remove module-specific args from args to avoid conflicts
            if hasattr(args, 'advntr_reference'):
                delattr(args, 'advntr_reference')
        else:
            module_args['advntr'] = {}

        # Determine bwa_reference based on reference_assembly
        if args.reference_assembly == "hg19":
            bwa_reference = config.get("reference_data", {}).get("bwa_reference_hg19")
        else:
            bwa_reference = config.get("reference_data", {}).get("bwa_reference_hg38")

        # Pass module_args and new arguments to run_pipeline
        run_pipeline(
            bwa_reference=bwa_reference,
            output_dir=Path(args.output_dir),
            extra_modules=args.extra_modules,
            module_args=module_args,
            config=config,
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            bam=args.bam,
            threads=args.threads,
            reference_assembly=args.reference_assembly,
            fast_mode=args.fast_mode,
            keep_intermediates=args.keep_intermediates,
            delete_intermediates=args.delete_intermediates,
            archive_results=args.archive_results,
            archive_format=args.archive_format,
            custom_regions=args.custom_regions,
            bed_file=args.bed_file,
            log_level=log_level
        )

    elif args.command == "fastq":
        process_fastq(
            fastq_1=args.fastq1,
            fastq_2=args.fastq2,
            threads=args.threads,
            output=Path(args.output_dir),
            output_name=args.output_name,
            config=config
        )

    elif args.command == "bam":
        process_bam_to_fastq(
            in_bam=args.alignment,
            output=Path(args.output_dir),
            output_name=args.output_name,
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
            kestrel_path=config.get("tools", {}).get("kestrel"),
            config=config,
            log_level=log_level
        )

    elif args.command == "report":
        generate_summary_report(
            output_dir=Path(args.output_dir),
            template_dir=config.get('paths', {}).get('template_dir', 'vntyper/templates'),
            report_file=args.report_file,
            log_file=args.log_file,
            bed_file=args.bed_file,
            bam_file=args.bam_file,
            fasta_file=args.reference_fasta,
            flanking=args.flanking,
            input_files={},  # You may need to populate this based on your pipeline
            pipeline_version=VERSION,
        )

    elif args.command == "cohort":
        aggregate_cohort(
            input_dirs=args.input_dirs,
            output_dir=Path(args.output_dir),
            summary_file=args.summary_file,
            config=config
        )

    else:
        logging.error(f"Unknown command: {args.command}")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
