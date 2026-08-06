"""
cli_handlers.py

Module Purpose:
---------------
One handler function per ``vntyper`` subcommand, so that ``vntyper/cli.py``'s
``main()`` is nothing but parse-then-dispatch.

Every subcommand body used to live inline in ``main()``, which meant any two
changes to two different subcommands were changes to the same function. Splitting
them apart lets each subcommand be edited, tested and reverted on its own.

The ``report`` handler deliberately lives in its own module
(:mod:`vntyper.scripts.cli_report`) rather than here: it is the one subcommand
known to be broken (AGENTS.md trap 11, fixed in #179) and is owned separately.

Logging is **not** configured here. ``cli.py`` remains the sole place that calls
``setup_logging()``; handlers receive the already-resolved log level and log file
so they can forward them to the stage they invoke.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Protocol

from vntyper.scripts.cohort_summary import aggregate_cohort
from vntyper.scripts.install_references import main as install_references_main

# Import the online mode function
from vntyper.scripts.online_mode import run_online_mode
from vntyper.scripts.pipeline import run_pipeline

logger = logging.getLogger(__name__)


def get_conf(config: dict[str, Any], key: str, fallback: Any) -> Any:
    """Read a CLI default out of ``config["default_values"]``.

    Args:
        config: The loaded configuration mapping.
        key: The ``default_values`` key to read.
        fallback: Value to return when the key is absent.

    Returns:
        Any: The configured value, or ``fallback``.
    """
    return config.get("default_values", {}).get(key, fallback)


class CommandHandler(Protocol):
    """The uniform signature every subcommand handler is dispatched through.

    The signature is uniform rather than minimal so that ``main()`` can dispatch
    through a table and never needs editing when a handler's needs change. A
    handler that does not need a parameter simply ignores it.
    """

    def __call__(
        self,
        args: argparse.Namespace,
        config: dict[str, Any],
        parser: argparse.ArgumentParser,
        log_level_value: int,
        log_file_str: str | None,
    ) -> None:
        """Run the subcommand.

        Args:
            args: The parsed arguments.
            config: The loaded configuration mapping.
            parser: The parser, for ``parser.error()`` on invalid combinations.
            log_level_value: The already-resolved ``logging`` level.
            log_file_str: The already-resolved log file path, or None.
        """
        ...  # pragma: no cover - structural typing only


def handle_install_references(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``install-references`` subcommand.

    Args:
        args: The parsed arguments.
        config: Unused; present for the uniform handler signature.
        parser: Unused; present for the uniform handler signature.
        log_level_value: Unused; present for the uniform handler signature.
        log_file_str: Unused; present for the uniform handler signature.
    """
    install_references_main(
        output_dir=args.output_dir,
        config_path=args.config_path,
        skip_indexing=args.skip_indexing,
        index_threads=args.threads,
        aligners_to_use=args.aligners,
        references_to_process=args.references,
    )
    sys.exit(0)


def handle_pipeline(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``pipeline`` subcommand.

    Args:
        args: The parsed arguments.
        config: The loaded configuration mapping.
        parser: The parser, used to reject invalid input combinations.
        log_level_value: The already-resolved ``logging`` level.
        log_file_str: The already-resolved log file path, or None.
    """
    # If log_file was not explicitly provided and output_dir is set, ensure log_file is correctly set
    if not args.log_file and args.output_dir:
        log_file = Path(args.output_dir) / "pipeline.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Setting log file to {log_file}")

    if args.output_dir is None:
        args.output_dir = get_conf(config, "output_dir", "out")
        logger.debug(f"output_dir set to {args.output_dir}")
    if args.threads is None:
        args.threads = get_conf(config, "threads", 4)
        logger.debug(f"threads set to {args.threads}")
    if args.reference_assembly is None:
        args.reference_assembly = get_conf(config, "reference_assembly", "hg19")
        logger.debug(f"reference_assembly set to {args.reference_assembly}")
    if args.output_name is None:
        args.output_name = get_conf(config, "output_name", "processed")
        logger.debug(f"output_name set to {args.output_name}")
    if args.archive_format is None:
        args.archive_format = get_conf(config, "archive_format", "zip")
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


def handle_cohort(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``cohort`` subcommand.

    Args:
        args: The parsed arguments.
        config: The loaded configuration mapping.
        parser: Unused; present for the uniform handler signature.
        log_level_value: Unused; present for the uniform handler signature.
        log_file_str: Unused; present for the uniform handler signature.
    """
    if args.summary_file is None:
        args.summary_file = get_conf(config, "summary_file", "cohort_summary.html")
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


def handle_online(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``online`` subcommand.

    Args:
        args: The parsed arguments.
        config: The loaded configuration mapping.
        parser: Unused; present for the uniform handler signature.
        log_level_value: Unused; present for the uniform handler signature.
        log_file_str: Unused; present for the uniform handler signature.
    """
    # No need to set up logging here; it's already set up in cli.py

    if args.output_dir is None:
        args.output_dir = get_conf(config, "output_dir", "out")
        logger.debug(f"output_dir set to {args.output_dir}")
    if args.reference_assembly is None:
        args.reference_assembly = get_conf(config, "reference_assembly", "hg19")
        logger.debug(f"reference_assembly set to {args.reference_assembly}")
    if args.threads is None:
        args.threads = get_conf(config, "threads", 4)
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
