"""Calibration CLI dispatch with atomic directory installation."""

from __future__ import annotations

import argparse
import logging
from collections.abc import Callable
from pathlib import Path

from vntyper.scripts.calibration_atomic_io import atomic_output

logger = logging.getLogger(__name__)

CalibrationOperation = Callable[[argparse.Namespace, Path], bool]


def handle_calibrate(
    args: argparse.Namespace,
    config: dict[str, object],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run one closed calibration operation and atomically install its output.

    Args:
        args: Parsed calibration arguments.
        config: Unused shared CLI configuration.
        parser: Unused top-level argument parser.
        log_level_value: Already configured logging level.
        log_file_str: Already configured logging destination.

    Raises:
        ValueError: If the operation or output path is invalid, or processing fails.
    """
    del config, parser, log_level_value, log_file_str
    operation = getattr(args, "calibration_operation", None)
    if not isinstance(operation, str):
        raise ValueError("calibration operation must be a string")
    if operation == "intake":
        from vntyper.scripts.cli_calibration_intake import run_calibration_intake

        run_calibration_intake(args)
        return
    producer = OPERATIONS.get(operation)
    if producer is None:
        message = f"unsupported calibration operation: {operation!r}"
        logger.error(message)
        raise ValueError(message)
    successful = _atomic_output(args.output, lambda staging: producer(args, staging))
    if not successful:
        logger.error(
            f"calibration {operation} completed with a failed outcome; "
            f"its complete failed output is installed at {args.output}"
        )
        raise SystemExit(1)


def _atomic_output(output: Path, producer: Callable[[Path], bool]) -> bool:
    """Delegate CLI publication to the shared no-clobber directory primitive."""
    return atomic_output(output, producer)


def _extract(args: argparse.Namespace, output: Path) -> bool:
    """Extract a validated calibration evidence bundle."""
    from vntyper.scripts.calibration_artifacts import extract_artifact_bundle

    return extract_artifact_bundle(args.truth, args.partitions, args.runs, output)


def _fit(args: argparse.Namespace, output: Path) -> bool:
    """Fit one candidate from an extracted evidence bundle."""
    from vntyper.scripts.calibration_artifacts import fit_artifact_bundle

    return fit_artifact_bundle(args.evidence, args.objective, output)


def _validate(args: argparse.Namespace, output: Path) -> bool:
    """Validate a fixed generated profile without selecting another."""
    from vntyper.scripts.calibration_artifacts import validate_artifact_bundle

    return validate_artifact_bundle(args.profile, args.evidence, output)


def _evaluate(args: argparse.Namespace, output: Path) -> bool:
    """Evaluate a fixed profile against one-use locked evidence."""
    from vntyper.scripts.calibration_artifacts import evaluate_artifact_bundle

    return evaluate_artifact_bundle(args.profile, args.evidence, output)


OPERATIONS: dict[str, CalibrationOperation] = {
    "extract": _extract,
    "fit": _fit,
    "validate": _validate,
    "evaluate": _evaluate,
}
