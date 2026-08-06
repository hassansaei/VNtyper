# vntyper/cli.py
# VNtyper CLI entry point

import importlib.resources as pkg_resources  # For accessing package data
import json
import logging
import sys
from pathlib import Path

from vntyper.scripts.cli_handlers import (
    CommandHandler,
    handle_cohort,
    handle_install_references,
    handle_online,
    handle_pipeline,
)
from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.cli_report import handle_report
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

    # Load the main configuration based on the provided config path
    try:
        # IMPORTANT: Use `initial_config` as fallback when `--config-path` is not provided
        config = load_config(args.config_path) if args.config_path else load_config(None)
        logger.debug("Configuration loaded successfully.")
    except Exception as exc:
        logger.critical(f"Failed to load configuration: {exc}")
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

    # Ensure that the log file directory exists
    if log_file_value:
        log_file_path = Path(log_file_value)
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        log_file_str = str(log_file_path)
    else:
        log_file_str = None

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

    command_handler(
        args,
        config=config,
        parser=parser,
        log_level_value=log_level_value,
        log_file_str=log_file_str,
    )


if __name__ == "__main__":
    main()
