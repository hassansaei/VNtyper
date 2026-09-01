# vntyper/cli.py
# VNtyper CLI entry point

import importlib.resources as pkg_resources  # For accessing package data
import json
import logging
import sys
from decimal import Decimal
from pathlib import Path

from vntyper.scripts.artifact_names import validate_output_name
from vntyper.scripts.cli_handlers import (
    CommandHandler,
    handle_cohort,
    handle_install_references,
    handle_online,
    handle_pipeline,
)
from vntyper.scripts.cli_logging_safety import validate_pipeline_log_destination
from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.cli_report import handle_report
from vntyper.scripts.fastp_cutoffs import preserve_exact_fastp_thresholds
from vntyper.scripts.run_configuration import resolve_run_configuration
from vntyper.scripts.utils import setup_logging

logger = logging.getLogger(__name__)

#: Subcommand name -> handler. ``main()`` does nothing but resolve this table, so a
#: change to one subcommand is a change to one function in one module.
HANDLERS: dict[str, CommandHandler] = {
    "install-references": handle_install_references,
    "pipeline": handle_pipeline,
    "report": handle_report,
    "cohort": handle_cohort,
    "online": handle_online,
}


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
            config_text = f.read()
        config = json.loads(config_text)
        exact_config = json.loads(config_text, parse_float=Decimal)
        config = preserve_exact_fastp_thresholds(config, exact_config)
        logger.debug(f"Loaded configuration from {config_path}")
    else:
        # No config path provided or file does not exist; use default config from package data
        try:
            with pkg_resources.open_text("vntyper", "config.json") as f:
                config_text = f.read()
            config = json.loads(config_text)
            exact_config = json.loads(config_text, parse_float=Decimal)
            config = preserve_exact_fastp_thresholds(config, exact_config)
            logger.debug("Loaded default configuration from package data.")
        except Exception as exc:
            logger.error("Error: Default config file not found in package data.")
            logger.error(exc)
            sys.exit(1)
    return config


def main(argv: list[str] | None = None) -> None:
    """
    Parse arguments and dispatch to the subcommand's handler.

    The global options (``-l/--log-level``, ``-f/--log-file``, ``-v/--version``,
    ``--config-path``) are top-level only: they must appear **before** the
    subcommand. ``build_parser`` registers them on the shared parent parser but
    does not pass that parent to ``add_parser``, so ``vntyper pipeline
    --log-level DEBUG`` is a usage error and always has been. This docstring
    previously claimed either position worked;
    ``tests/unit/test_cli_parser_contract.py`` now pins the real behaviour on
    every subcommand, so making both positions work would be a visible change
    rather than an accidental one.

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

    # Resolve the complete profile before loading stage configuration and, most
    # importantly, before creating the output directory or opening a log file. Only
    # `pipeline` executes genotyping decisions; standalone report/cohort commands
    # verify the snapshot recorded by an existing run and never infer current defaults.
    if args.command == "pipeline":
        try:
            args.run_configuration = resolve_run_configuration(args.decision_profile)
        except ValueError as exc:
            logger.critical(str(exc))
            sys.exit(1)

    # Load the main configuration based on the provided config path
    try:
        # IMPORTANT: Use `initial_config` as fallback when `--config-path` is not provided
        config = load_config(args.config_path) if args.config_path else load_config(None)
        logger.debug("Configuration loaded successfully.")
    except Exception as exc:
        logger.critical(f"Failed to load configuration: {exc}")
        sys.exit(1)

    # Refuse a flag the pipeline cannot honour *before* anything is created on
    # disk. The log file for a `pipeline` run is `<output_dir>/pipeline.log`, so
    # validating this in the handler meant the output directory was created and
    # written to before the run died -- and it died with an unhandled ValueError,
    # which is a traceback for what is really a usage error. Exit 1 rather than
    # `parser.error`'s 2: the parser is already built and dispatched by this point,
    # so this reads as a completed run that failed, not as a usage error caught at
    # parse time. `cli_handlers` does use `parser.error` -- and therefore exit 2 --
    # for the input combinations it rejects while parsing.
    if args.command == "pipeline" and getattr(args, "output_name", None) is not None:
        try:
            validate_output_name(args.output_name)
        except ValueError as exc:
            logger.critical(str(exc))
            sys.exit(1)

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

    log_file_path = None
    if log_file_value:
        log_file_path = Path(log_file_value)
        log_file_str = str(log_file_path)
    else:
        log_file_str = None

    # FileHandler opens its destination immediately. Refuse a pipeline log that
    # aliases input state before even its parent-directory setup can run.
    if args.command == "pipeline" and log_file_path is not None:
        try:
            validate_pipeline_log_destination(log_file_path, args, config)
        except ValueError as exc:
            logger.critical(str(exc))
            sys.exit(1)

    if log_file_path is not None:
        log_file_path.parent.mkdir(parents=True, exist_ok=True)

    # Setup logging now (only once) with the determined log level and log file.
    # This is the sole place logging is configured; handlers never reconfigure it.
    setup_logging(log_level=log_level_value, log_file=log_file_str)
    logger.debug(f"Logging has been set up with level {log_level_value} and log_file {log_file_str}")

    # Log current logging handlers and their levels
    for handler in logging.getLogger().handlers:
        handler_type = handler.__class__.__name__
        handler_level = logging.getLevelName(handler.level)
        handler_file = getattr(handler, "baseFilename", "N/A") if isinstance(handler, logging.FileHandler) else "N/A"
        logger.debug(f"Handler: {handler_type}, Level: {handler_level}, File: {handler_file}")

    # Dispatch to the subcommand's handler
    command_handler = HANDLERS.get(args.command)
    if command_handler is None:
        logger.error(f"Unknown command: {args.command}")
        parser.print_help()
        sys.exit(1)

    try:
        command_handler(
            args,
            config=config,
            parser=parser,
            log_level_value=log_level_value,
            log_file_str=log_file_str,
        )
    except ValueError as exc:
        # AGENTS.md's convention for a handler-level validation failure is
        # `logger.error(msg)` followed by `raise ValueError(msg)` - no custom exception
        # classes - so the diagnostic has already been logged before this runs. Without
        # this catch, that clean convention is undone one frame up: the exception
        # escapes `main()` as an uncaught traceback instead of the `sys.exit(1)` every
        # other CLI-level validation error already gets here (see the `--output-name`
        # and `--log-file` checks above). `_resolve_bwa_reference`'s new fail-closed
        # ValueError (a present-but-uninstalled reference) is what surfaced this gap,
        # but the fix is general: any handler that raises `ValueError` gets the same
        # clean presentation, not just `pipeline`.
        logger.critical(str(exc))
        sys.exit(1)


if __name__ == "__main__":
    main()
