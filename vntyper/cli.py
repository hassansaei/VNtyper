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

# Import the online mode function
from vntyper.scripts.online_mode import run_online_mode


def load_config(config_path=None):
    """
    Load the configuration file with fallback to the default package config.

    Args:
        config_path (Path or None): Path to the user-provided config file.

    Returns:
        dict: The loaded configuration dictionary.
    """
    if config_path is not None and Path(config_path).exists():
        with open(config_path, 'r', encoding='utf-8') as f:
            config = json.load(f)
    else:
        # No config path provided or file does not exist; use default config from package data
        try:
            with pkg_resources.open_text('vntyper', 'config.json') as f:
                config = json.load(f)
        except Exception as exc:
            logging.error("Error: Default config file not found in package data.")
            logging.error(exc)
            sys.exit(1)
    return config


def main():
    """
    Main function to parse arguments and execute corresponding subcommands.
    With this setup, global parameters can now be placed before or after
    the subcommand.
    """

    # We'll load an initial config for CLI defaults. We do a try/except because
    # we might overwrite this again if the user specifies --config-path.
    try:
        initial_config = load_config(None)
    except Exception:
        initial_config = {}

    # Fallback lookups for CLI defaults
    default_cli = initial_config.get("cli_defaults", {})
    default_log_level = default_cli.get("log_level", "INFO")
    default_log_file = default_cli.get("log_file", None)

    # Parent parser for global arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        '-l', '--log-level',
        help="Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR)",
        default=default_log_level
    )
    parent_parser.add_argument(
        '-f', '--log-file',
        help="Set the log output file (default is stdout)",
        default=default_log_file
    )
    parent_parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'%(prog)s {VERSION}'
    )
    parent_parser.add_argument(
        '--config-path',
        type=Path,
        default=None,
        help=(
            "Path to the configuration file (config.json). "
            "If not provided, the default config will be used."
        ),
        required=False
    )

    # Main parser that includes the parent parser
    parser = argparse.ArgumentParser(
        description="VNtyper CLI: A pipeline for genotyping MUC1-VNTR.",
        parents=[parent_parser]
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subcommand: pipeline
    # We read further defaults from the loaded config (once we re-load config in "run" step).
    parser_pipeline = subparsers.add_parser(
        "pipeline",
        help="Run the full VNtyper pipeline.",
        parents=[parent_parser]
    )
    parser_pipeline.add_argument(
        '-o', '--output-dir',
        type=str,
        default=None,
        help="Output directory for the results."
    )
    # Changed here (#57): allow --extra-modules multiple times
    parser_pipeline.add_argument(
        '--extra-modules',
        action='append',
        default=[],
        help="Optional extra modules to include (e.g., advntr, shark). "
             "Can be repeated multiple times."
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
        '--cram',
        type=str,
        help="Path to the CRAM file."
    )
    parser_pipeline.add_argument(
        '--threads',
        type=int,
        default=None,
        help="Number of threads to use."
    )
    parser_pipeline.add_argument(
        '--reference-assembly',
        type=str,
        choices=["hg19", "hg38"],
        default=None,
        help="Specify the reference assembly used for the input "
             "BAM/CRAM file alignment."
    )
    parser_pipeline.add_argument(
        '--fast-mode',
        action='store_true',
        help="Enable fast mode (skips filtering for unmapped "
             "and partially mapped reads)."
    )
    parser_pipeline.add_argument(
        '--keep-intermediates',
        action='store_true',
        help="Keep intermediate files (e.g., BAM slices, temporary files)."
    )
    parser_pipeline.add_argument(
        '--delete-intermediates',
        action='store_true',
        help="Delete intermediate files after processing "
             "(overrides --keep-intermediates)."
    )
    parser_pipeline.add_argument(
        '--archive-results',
        action='store_true',
        help="Create an archive of the results folder after "
             "pipeline completion."
    )
    parser_pipeline.add_argument(
        '--archive-format',
        type=str,
        choices=['zip', 'tar.gz'],
        default=None,
        help="Format of the archive: 'zip' or 'tar.gz'."
    )
    parser_pipeline.add_argument(
        '-n', '--output-name',
        type=str,
        default=None,
        help="Base name for the output files."
    )
    parser_pipeline.add_argument(
        '-s', '--sample-name',
        type=str,
        default=None,
        help=("Set the sample name for labeling results. If not provided, "
              "defaults to input BAM or FASTQ name.")
    )
    region_group = parser_pipeline.add_mutually_exclusive_group()
    region_group.add_argument(
        '--custom-regions',
        type=str,
        help="Define custom regions for MUC1 analysis as comma-separated "
             "values (e.g., chr1:1000-2000,chr2:3000-4000)."
    )
    region_group.add_argument(
        '--bed-file',
        type=Path,
        help="Path to a BED file specifying regions for MUC1 analysis."
    )

    # Subcommand: report
    parser_report = subparsers.add_parser(
        "report",
        help="Generate a summary report and visualizations from output data.",
        parents=[parent_parser]
    )
    parser_report.add_argument(
        '-o', '--output-dir',
        type=str,
        required=True,
        help="Directory containing pipeline results (subdirectories, etc.)."
    )
    parser_report.add_argument(
        '--input-dir',
        type=Path,
        default=None,
        help="If provided, search this directory (and its subdirs) for standard pipeline output filenames."
    )
    parser_report.add_argument(
        '--report-file',
        type=str,
        default=None,
        help="Name of the output report file."
    )
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
        default=None,
        help="Flanking region size for IGV reports."
    )

    # Subcommand: cohort
    parser_cohort = subparsers.add_parser(
        "cohort",
        help="Aggregate outputs from multiple runs into a single summary file.",
        parents=[parent_parser]
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
        default=None,
        help="Name of the cohort summary report file."
    )

    # Subcommand: install-references
    parser_install = subparsers.add_parser(
        "install-references",
        help="Download and set up necessary reference files.",
        parents=[parent_parser]
    )
    parser_install.add_argument(
        '-d', '--output-dir',
        type=Path,
        required=True,
        help="Directory where references will be installed."
    )
    parser_install.add_argument(
        '--skip-indexing',
        action='store_true',
        help="Skip the bwa indexing step."
    )

    # Subcommand: online
    parser_online = subparsers.add_parser(
        "online",
        help=("Subset the BAM and submit it to an online vntyper instance, "
              "then retrieve results."),
        parents=[parent_parser]
    )
    parser_online.add_argument(
        '--bam',
        type=str,
        required=True,
        help="Path to the input BAM file."
    )
    parser_online.add_argument(
        '-o', '--output-dir',
        type=str,
        default=None,
        help="Output directory for results."
    )
    parser_online.add_argument(
        '--reference-assembly',
        type=str,
        choices=["hg19", "hg38"],
        default=None,
        help="Reference assembly used."
    )
    parser_online.add_argument(
        '--threads',
        type=int,
        default=None,
        help="Number of threads to use."
    )
    parser_online.add_argument(
        '--email',
        type=str,
        default=None,
        help="Email to receive notifications (optional)."
    )
    parser_online.add_argument(
        '--cohort-id',
        type=str,
        default=None,
        help="Cohort ID to associate the job with (optional)."
    )
    parser_online.add_argument(
        '--passphrase',
        type=str,
        default=None,
        help="Passphrase for the cohort (if required)."
    )
    parser_online.add_argument(
        '--resume',
        action='store_true',
        help="Resume polling a previously submitted job if job_id is found."
    )

    # Parse initial arguments just enough to see what command is used
    args = parser.parse_args()

    # Display help if no command is provided
    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # We create a base logger with user-specified or config-based defaults
    # (we may override or add additional logging in each subcommand).
    try:
        base_config = load_config(args.config_path) if args.config_path else initial_config
    except Exception as exc:
        logging.critical(f"Failed to load configuration: {exc}")
        sys.exit(1)

    def get_conf(key, fallback):
        return base_config.get("default_values", {}).get(key, fallback)

    # Setup base logging now
    log_level_value = getattr(logging, args.log_level.upper(), logging.INFO)
    if args.log_file:
        log_file_path = Path(args.log_file)
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        setup_logging(log_level=log_level_value, log_file=str(log_file_path))
    else:
        fallback_file = get_conf("log_file", None)
        if fallback_file:
            lf_path = Path(fallback_file)
            lf_path.parent.mkdir(parents=True, exist_ok=True)
            setup_logging(log_level=log_level_value, log_file=str(lf_path))
        else:
            setup_logging(log_level=log_level_value, log_file=None)

    # Subcommand: install-references
    if args.command == "install-references":
        install_references_main(
            output_dir=args.output_dir,
            config_path=args.config_path,
            skip_indexing=args.skip_indexing
        )
        sys.exit(0)

    # Load the config again for relevant subcommands
    # (some subcommands, like 'report' or 'cohort', also need config).
    try:
        config = load_config(args.config_path) if args.config_path else initial_config
    except Exception as exc:
        logging.critical(f"Failed to load configuration: {exc}")
        sys.exit(1)

    #
    # -------------------------------------------------------------------------
    # pipeline command
    # -------------------------------------------------------------------------
    #
    if args.command == "pipeline":
        # Possibly unify logs into output_dir/pipeline.log if not specified
        if not args.log_file:
            if args.output_dir:
                log_file = Path(args.output_dir) / "pipeline.log"
                log_file.parent.mkdir(parents=True, exist_ok=True)
                setup_logging(log_level=log_level_value, log_file=str(log_file))
            else:
                # fallback to current dir if no output_dir
                log_file = Path("pipeline.log")
                setup_logging(log_level=log_level_value, log_file=str(log_file))

        if args.output_dir is None:
            args.output_dir = get_conf("output_dir", "out")
        if args.threads is None:
            args.threads = get_conf("threads", 4)
        if args.reference_assembly is None:
            args.reference_assembly = get_conf("reference_assembly", "hg19")
        if args.output_name is None:
            args.output_name = get_conf("output_name", "processed")
        if args.archive_format is None:
            args.archive_format = get_conf("archive_format", "zip")

        import itertools
        flattened_modules = list(itertools.chain.from_iterable(
            m if isinstance(m, list) else [m] for m in args.extra_modules
        ))

        # Validate single input type
        input_types = sum([
            1 if args.bam else 0,
            1 if args.cram else 0,
            1 if (args.fastq1 or args.fastq2) else 0
        ])
        if input_types > 1:
            parser_pipeline.error("Provide either BAM, CRAM, or FASTQ files (not multiples).")

        if (not args.bam and not args.cram
                and (args.fastq1 is None or args.fastq2 is None)):
            parser_pipeline.error(
                "When not providing BAM/CRAM, both --fastq1 and --fastq2 must "
                "be specified for paired-end sequencing."
            )

        # Construct module_args_dict for advntr, etc.
        module_args_dict = {}
        if 'advntr' in flattened_modules:
            if hasattr(args, 'advntr_reference'):
                module_args_dict['advntr'] = {
                    'advntr_reference': args.advntr_reference
                }
                delattr(args, 'advntr_reference')
            else:
                module_args_dict['advntr'] = {}
        else:
            module_args_dict['advntr'] = {}

        # (#62) If user tries to use 'shark' in BAM/CRAM mode, exit with a warning
        if (args.bam or args.cram) and ('shark' in flattened_modules):
            logging.warning(
                "Shark is not supported in BAM mode; please use FASTQ mode or remove the shark flag."
            )
            sys.exit(1)

        # Determine which BWA reference to use from config
        if args.reference_assembly == "hg19":
            bwa_reference = config.get("reference_data", {}).get("bwa_reference_hg19")
        else:
            bwa_reference = config.get("reference_data", {}).get("bwa_reference_hg38")

        sample_name_val = args.sample_name
        if sample_name_val is None:
            if args.bam:
                sample_name_val = Path(args.bam).stem
            elif args.fastq1:
                sample_name_val = Path(args.fastq1).stem
            else:
                sample_name_val = "sample"

        run_pipeline(
            bwa_reference=bwa_reference,
            output_dir=Path(args.output_dir),
            extra_modules=flattened_modules,
            module_args=module_args_dict,
            config=config,
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            bam=args.bam,
            cram=args.cram,
            threads=args.threads,
            reference_assembly=args.reference_assembly,
            fast_mode=args.fast_mode,
            keep_intermediates=args.keep_intermediates,
            delete_intermediates=args.delete_intermediates,
            archive_results=args.archive_results,
            archive_format=args.archive_format,
            custom_regions=args.custom_regions,
            bed_file=args.bed_file,
            log_level=log_level_value,
            sample_name=sample_name_val,
        )

    #
    # -------------------------------------------------------------------------
    # report command
    # -------------------------------------------------------------------------
    #
    elif args.command == "report":
        # Possibly unify logs into output_dir/report.log if not specified
        if not args.log_file:
            report_log_file = Path(args.output_dir) / "report.log"
            report_log_file.parent.mkdir(parents=True, exist_ok=True)
            setup_logging(log_level=log_level_value, log_file=str(report_log_file))
        else:
            # If user gave log_file, we already set it up above
            report_log_file = Path(args.log_file)

        # Fallback from config if not specified on CLI
        # (we do so only if relevant keys exist in config["default_values"])
        if args.report_file is None:
            args.report_file = get_conf("report_file", "summary_report.html")
        if args.flanking is None:
            args.flanking = get_conf("flanking", 50)

        # If user did not provide reference_fasta, fallback to config if present
        if args.reference_fasta is None:
            # for example: config["reference_data"]["muc1_reference_vntr"]
            ref_fasta = config.get("reference_data", {}).get("muc1_reference_vntr")
            if ref_fasta:
                args.reference_fasta = Path(ref_fasta)

        # If user didn't specify --bam-file, we attempt to find a standard location
        # e.g. <input-dir>/kestrel/output.bam
        # (only if --input-dir was provided)
        if args.bam_file is None and args.input_dir:
            candidate_bam = args.input_dir / "kestrel" / "output.bam"
            if candidate_bam.exists():
                args.bam_file = candidate_bam

        # Same approach for bed-file (standard name is "output.bed" in "kestrel")
        if args.bed_file is None and args.input_dir:
            candidate_bed = args.input_dir / "kestrel" / "output.bed"
            if candidate_bed.exists():
                args.bed_file = candidate_bed

        # Now call generate_summary_report
        generate_summary_report(
            output_dir=Path(args.output_dir),
            template_dir=config.get('paths', {}).get('template_dir', 'vntyper/templates'),
            report_file=args.report_file,
            log_file=str(report_log_file),
            bed_file=args.bed_file,
            bam_file=args.bam_file,
            fasta_file=args.reference_fasta,
            flanking=args.flanking,
            input_files={},  # Optionally populate if you want to reference them in the final report
            pipeline_version=VERSION,
            # We pass in the entire config, since we've made it consistent with the pipeline usage.
            config=config
        )

    #
    # -------------------------------------------------------------------------
    # cohort command
    # -------------------------------------------------------------------------
    #
    elif args.command == "cohort":
        # Possibly unify logs into output_dir/pipeline.log if not specified
        if not args.log_file:
            if args.output_dir:
                log_file = Path(args.output_dir) / "pipeline.log"
                log_file.parent.mkdir(parents=True, exist_ok=True)
                setup_logging(log_level=log_level_value, log_file=str(log_file))
            else:
                # fallback to current dir
                log_file = Path("pipeline.log")
                setup_logging(log_level=log_level_value, log_file=str(log_file))

        if args.summary_file is None:
            args.summary_file = get_conf("summary_file", "cohort_summary.html")

        aggregate_cohort(
            input_dirs=args.input_dirs,
            output_dir=Path(args.output_dir),
            summary_file=args.summary_file,
            config=config
        )

    #
    # -------------------------------------------------------------------------
    # online command
    # -------------------------------------------------------------------------
    #
    elif args.command == "online":
        # Possibly unify logs into output_dir/pipeline.log if not specified
        if not args.log_file:
            if args.output_dir:
                log_file = Path(args.output_dir) / "pipeline.log"
                log_file.parent.mkdir(parents=True, exist_ok=True)
                setup_logging(log_level=log_level_value, log_file=str(log_file))
            else:
                # fallback to current dir
                log_file = Path("pipeline.log")
                setup_logging(log_level=log_level_value, log_file=str(log_file))

        if args.output_dir is None:
            args.output_dir = get_conf("output_dir", "out")
        if args.reference_assembly is None:
            args.reference_assembly = get_conf("reference_assembly", "hg19")
        if args.threads is None:
            args.threads = get_conf("threads", 4)

        run_online_mode(
            config=config,
            bam=args.bam,
            output_dir=args.output_dir,
            reference_assembly=args.reference_assembly,
            threads=args.threads,
            email=args.email,
            cohort_id=args.cohort_id,
            passphrase=args.passphrase,
            resume=args.resume
        )

    else:
        logging.error(f"Unknown command: {args.command}")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
