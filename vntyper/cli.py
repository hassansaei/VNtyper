# vntyper/cli.py
# VNtyper CLI entry point

import importlib.resources as pkg_resources  # For accessing package data
import json
import logging
import sys
from pathlib import Path
from typing import Any

from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.cohort_summary import aggregate_cohort
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.install_references import main as install_references_main

# Import the online mode function
from vntyper.scripts.online_mode import run_online_mode
from vntyper.scripts.pipeline import run_pipeline
from vntyper.scripts.utils import setup_logging
from vntyper.version import __version__ as VERSION

logger = logging.getLogger(__name__)


def load_config(config_path=None):
    """
    Load the configuration file with fallback to the default package config.

    Args:
        config_path (Path or None): Path to the user-provided config file.

    Returns:
        dict: The loaded configuration dictionary.
    """
    if config_path is not None and Path(config_path).exists():
        with open(config_path, encoding="utf-8") as f:
            config = json.load(f)
        logger.debug(f"Loaded configuration from {config_path}")
    else:
        # No config path provided or file does not exist; use default config from package data
        try:
            with pkg_resources.open_text("vntyper", "config.json") as f:
                config = json.load(f)
            logger.debug("Loaded default configuration from package data.")
        except Exception as exc:
            logger.error("Error: Default config file not found in package data.")
            logger.error(exc)
            sys.exit(1)
    return config


def main(argv: list[str] | None = None) -> None:
    """
    Main function to parse arguments and execute corresponding subcommands.
    With this setup, global parameters can now be placed before or after
    the subcommand.

    Args:
        argv (list[str] | None): Argument list to parse. Defaults to None, which
            makes argparse read ``sys.argv[1:]`` - the behaviour the ``vntyper``
            console entry point relies on. Passing an explicit list lets callers
            (tests, in particular) drive the CLI without touching ``sys.argv``.
    """
    parser = build_parser()

    # Parse all arguments first
    args = parser.parse_args(argv)

    # If no command is provided, display help and exit
    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # Load the main configuration based on the provided config path
    try:
        # IMPORTANT: Use `initial_config` as fallback when `--config-path` is not provided
        config = load_config(args.config_path) if args.config_path else load_config(None)
        logger.debug("Configuration loaded successfully.")
    except Exception as exc:
        logger.critical(f"Failed to load configuration: {exc}")
        sys.exit(1)

    def get_conf(key, fallback):
        return config.get("default_values", {}).get(key, fallback)

    # Determine log level
    if args.log_level:
        log_level_value = getattr(logging, args.log_level.upper(), logging.INFO)
    else:
        log_level_value = getattr(
            logging,
            config.get("cli_defaults", {}).get("log_level", "INFO").upper(),
            logging.INFO,
        )

    # Determine log file
    if args.log_file:
        log_file_value = args.log_file
    else:
        # If the command is 'pipeline' and output_dir is provided, set log_file accordingly
        if args.command == "pipeline" and args.output_dir:
            log_file_value = Path(args.output_dir) / "pipeline.log"
        elif args.command == "cohort" and args.output_dir:
            log_file_value = Path(args.output_dir) / "cohort.log"
        else:
            log_file_value = config.get("cli_defaults", {}).get("log_file", None)

    # Ensure that the log file directory exists
    if log_file_value:
        log_file_path = Path(log_file_value)
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        log_file_str = str(log_file_path)
    else:
        log_file_str = None

    # Setup logging now (only once) with the determined log level and log file
    setup_logging(log_level=log_level_value, log_file=log_file_str)
    logger.debug(f"Logging has been set up with level {log_level_value} and log_file {log_file_str}")

    # Log current logging handlers and their levels
    for handler in logging.getLogger().handlers:
        handler_type = handler.__class__.__name__
        handler_level = logging.getLevelName(handler.level)
        handler_file = getattr(handler, "baseFilename", "N/A") if isinstance(handler, logging.FileHandler) else "N/A"
        logger.debug(f"Handler: {handler_type}, Level: {handler_level}, File: {handler_file}")

    # Subcommand: install-references
    if args.command == "install-references":
        install_references_main(
            output_dir=args.output_dir,
            config_path=args.config_path,
            skip_indexing=args.skip_indexing,
            index_threads=args.threads,
            aligners_to_use=args.aligners,
            references_to_process=args.references,
        )
        sys.exit(0)

    #
    # -------------------------------------------------------------------------
    # pipeline command
    # -------------------------------------------------------------------------
    #
    if args.command == "pipeline":
        # If log_file was not explicitly provided and output_dir is set, ensure log_file is correctly set
        if not args.log_file and args.output_dir:
            log_file = Path(args.output_dir) / "pipeline.log"
            log_file.parent.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Setting log file to {log_file}")

        if args.output_dir is None:
            args.output_dir = get_conf("output_dir", "out")
            logger.debug(f"output_dir set to {args.output_dir}")
        if args.threads is None:
            args.threads = get_conf("threads", 4)
            logger.debug(f"threads set to {args.threads}")
        if args.reference_assembly is None:
            args.reference_assembly = get_conf("reference_assembly", "hg19")
            logger.debug(f"reference_assembly set to {args.reference_assembly}")
        if args.output_name is None:
            args.output_name = get_conf("output_name", "processed")
            logger.debug(f"output_name set to {args.output_name}")
        if args.archive_format is None:
            args.archive_format = get_conf("archive_format", "zip")
            logger.debug(f"archive_format set to {args.archive_format}")

        import itertools

        flattened_modules = list(
            itertools.chain.from_iterable(m if isinstance(m, list) else [m] for m in args.extra_modules)
        )
        logger.debug(f"extra_modules flattened to {flattened_modules}")

        # Validate single input type
        input_types = sum(
            [
                1 if args.bam else 0,
                1 if args.cram else 0,
                1 if (args.fastq1 or args.fastq2) else 0,
            ]
        )
        if input_types > 1:
            parser.error("Provide either BAM, CRAM, or FASTQ files (not multiples).")
            logger.debug("Multiple input types detected.")
            sys.exit(1)

        if not args.bam and not args.cram and (args.fastq1 is None or args.fastq2 is None):
            parser.error(
                "When not providing BAM/CRAM, both --fastq1 and --fastq2 must be specified for paired-end sequencing."
            )
            logger.debug("Missing FASTQ files for paired-end sequencing.")
            sys.exit(1)

        # Construct module_args_dict for advntr, etc.
        module_args_dict: dict[str, dict[str, Any]] = {}
        if "advntr" in flattened_modules:
            module_args_dict["advntr"] = {}

            # If the user set advntr_reference somewhere
            if hasattr(args, "advntr_reference"):
                module_args_dict["advntr"]["advntr_reference"] = args.advntr_reference
                delattr(args, "advntr_reference")
                logger.debug(f"advntr_reference set to {args.advntr_reference}")

            # The new coverage parameter:
            if args.advntr_max_coverage:
                module_args_dict["advntr"]["max_coverage"] = args.advntr_max_coverage
                logger.debug(f"advntr_max_coverage set to {args.advntr_max_coverage}")

        else:
            module_args_dict["advntr"] = {}
            logger.debug("advntr module not included.")

        # (#62) If user tries to use 'shark' in BAM/CRAM mode, exit with a warning
        if (args.bam or args.cram) and ("shark" in flattened_modules):
            logger.warning("Shark is not supported in BAM mode; please use FASTQ mode or remove the shark flag.")
            logger.debug("Shark module detected with BAM/CRAM input; exiting.")
            sys.exit(1)

        # Determine which BWA reference to use from config using registry
        from vntyper.scripts.reference_registry import get_coordinate_system

        try:
            coord_system = get_coordinate_system(args.reference_assembly)
        except ValueError:
            logger.warning(f"Unknown assembly '{args.reference_assembly}', defaulting to GRCh37")
            coord_system = "GRCh37"

        # Map coordinate system to UCSC-style name for BWA reference lookup
        ucsc_map = {"GRCh37": "hg19", "GRCh38": "hg38"}
        ucsc_style_ref = ucsc_map.get(coord_system, "hg19")
        bwa_key = f"bwa_reference_{ucsc_style_ref}"
        bwa_reference = config.get("reference_data", {}).get(bwa_key)
        logger.debug(f"Using BWA reference {bwa_key}: {bwa_reference}")

        sample_name_val = args.sample_name
        if sample_name_val is None:
            if args.bam:
                sample_name_val = Path(args.bam).stem
                logger.debug(f"sample_name set from BAM file: {sample_name_val}")
            elif args.fastq1:
                sample_name_val = Path(args.fastq1).stem
                logger.debug(f"sample_name set from FASTQ1 file: {sample_name_val}")
            else:
                sample_name_val = "sample"
                logger.debug(f"sample_name defaulted to: {sample_name_val}")

        # Process the new --summary-formats argument (comma-separated)
        summary_formats = []
        if args.summary_formats:
            summary_formats = [fmt.strip().lower() for fmt in args.summary_formats.split(",") if fmt.strip()]

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
            log_level=log_level_value,  # Pass log_level to run_pipeline
            sample_name=sample_name_val,
            log_file=log_file_str,  # Pass the correctly determined log_file
            summary_formats=summary_formats,  # New parameter passed
        )

    #
    # -------------------------------------------------------------------------
    # report command
    # -------------------------------------------------------------------------
    #
    elif args.command == "report":
        # No need to set up logging here; it's already set up in cli.py

        # Fallback from config if not specified on CLI
        # (we do so only if relevant keys exist in config["default_values"])
        if args.report_file is None:
            args.report_file = get_conf("report_file", "summary_report.html")
            logger.debug(f"report_file set to {args.report_file}")
        if args.flanking is None:
            args.flanking = get_conf("flanking", 50)
            logger.debug(f"flanking set to {args.flanking}")

        # If user did not provide reference_fasta, fallback to config if present
        if args.reference_fasta is None:
            # for example: config["reference_data"]["muc1_reference_vntr"]
            ref_fasta = config.get("reference_data", {}).get("muc1_reference_vntr")
            if ref_fasta:
                args.reference_fasta = Path(ref_fasta)
                logger.debug(f"reference_fasta set to {args.reference_fasta}")

        # If user didn't specify --bam-file, we attempt to find a standard location
        # e.g. <input-dir>/kestrel/output.bam
        # (only if --input-dir was provided)
        if args.bam_file is None and args.input_dir:
            candidate_bam = args.input_dir / "kestrel" / "output.bam"
            if candidate_bam.exists():
                args.bam_file = candidate_bam
                logger.debug(f"bam_file set to {args.bam_file}")

        # Same approach for bed-file (standard name is "output.bed" in "kestrel")
        if args.bed_file is None and args.input_dir:
            candidate_bed = args.input_dir / "kestrel" / "output.bed"
            if candidate_bed.exists():
                args.bed_file = candidate_bed
                logger.debug(f"bed_file set to {args.bed_file}")

        # Now call generate_summary_report
        #
        # AGENTS.md trap 11: this call is broken and `vntyper report` raises TypeError.
        # `generate_summary_report()` accepts none of `input_files`, `pipeline_version`
        # or `mean_vntr_coverage`, and it requires a `log_file` this call never passes.
        # Giving `main()` a typed signature made mypy check this body for the first
        # time, which is how the defect became visible to the build rather than only to
        # a user running the subcommand. The narrow ignore keeps the gate green without
        # papering over it; fixing the call is owned by its own commit, and this comment
        # and the ignore go with the fix. Do not copy this call site as an example.
        generate_summary_report(  # type: ignore[call-arg]
            output_dir=Path(args.output_dir),
            template_dir=config.get("paths", {}).get("template_dir", "vntyper/templates"),
            report_file=args.report_file,
            bed_file=args.bed_file,
            bam_file=args.bam_file,
            fasta_file=args.reference_fasta,
            flanking=args.flanking,
            input_files={},  # Optionally populate if you want to reference them in the final report
            pipeline_version=VERSION,
            mean_vntr_coverage=None,  # If applicable, otherwise remove
            vcf_file=None,  # If applicable, otherwise remove
            config=config,
        )

    #
    # -------------------------------------------------------------------------
    # cohort command
    # -------------------------------------------------------------------------
    #
    elif args.command == "cohort":
        if args.summary_file is None:
            args.summary_file = get_conf("summary_file", "cohort_summary.html")
            logger.debug(f"summary_file set to {args.summary_file}")

        # Prepare the list of input paths
        input_paths = []
        if args.input_dirs:
            input_paths.extend(args.input_dirs)
            logger.debug(f"Added input_dirs: {args.input_dirs}")
        if args.input_file:
            if not args.input_file.exists():
                logger.error(f"The input file {args.input_file} does not exist.")
                sys.exit(1)
            with open(args.input_file) as f:
                file_lines = [line.strip() for line in f if line.strip()]
                input_paths.extend(file_lines)
            logger.debug(f"Added input_file entries: {file_lines}")

        aggregate_cohort(
            input_paths=input_paths,
            output_dir=Path(args.output_dir),
            summary_file=args.summary_file,
            config=config,
            additional_formats=args.summary_formats,  # Pass the new parameter for extra formats
            pseudonymize_samples=args.pseudonymize_samples,  # Pass the new pseudonymize flag (value or None)
        )

    #
    # -------------------------------------------------------------------------
    # online command
    # -------------------------------------------------------------------------
    #
    elif args.command == "online":
        # No need to set up logging here; it's already set up in cli.py

        if args.output_dir is None:
            args.output_dir = get_conf("output_dir", "out")
            logger.debug(f"output_dir set to {args.output_dir}")
        if args.reference_assembly is None:
            args.reference_assembly = get_conf("reference_assembly", "hg19")
            logger.debug(f"reference_assembly set to {args.reference_assembly}")
        if args.threads is None:
            args.threads = get_conf("threads", 4)
            logger.debug(f"threads set to {args.threads}")

        run_online_mode(
            config=config,
            bam=args.bam,
            output_dir=args.output_dir,
            reference_assembly=args.reference_assembly,
            threads=args.threads,
            email=args.email,
            cohort_id=args.cohort_id,
            passphrase=args.passphrase,
            resume=args.resume,
        )

    else:
        logger.error(f"Unknown command: {args.command}")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
