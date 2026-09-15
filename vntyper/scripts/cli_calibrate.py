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
        parser: Top-level argument parser for conditional usage errors.
        log_level_value: Already configured logging level.
        log_file_str: Already configured logging destination.

    Raises:
        ValueError: If the operation or output path is invalid, or processing fails.
    """
    del config, log_level_value, log_file_str
    operation = getattr(args, "calibration_operation", None)
    if not isinstance(operation, str):
        raise ValueError("calibration operation must be a string")
    if operation == "intake":
        from vntyper.scripts.cli_calibration_intake import run_calibration_intake

        try:
            run_calibration_intake(args)
        except RuntimeError as error:
            message = str(error)
            logger.error(message)
            raise ValueError(message) from error
        return
    producer: Callable[[argparse.Namespace, Path], bool] | None
    if operation == "cohort":
        from vntyper.scripts.calibration_cohort import run_cohort_calibration

        producer = run_cohort_calibration
    elif operation == "optimize":
        from vntyper.scripts.calibration_cutoff_optimize import run_cutoff_optimization

        _optimize_arguments(args, parser)
        producer = run_cutoff_optimization
    else:
        target = getattr(args, "target", "dominance")
        _target_arguments(args, parser, target, operation)
        producer = OPERATIONS.get(operation) if target == "dominance" else TARGET_OPERATIONS.get((target, operation))
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


def _optimize_arguments(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Reject the two cross-argument combinations argparse alone cannot express.

    Args:
        args: Parsed ``calibrate optimize`` arguments.
        parser: Top-level parser, which owns the exit code for a usage error.
    """
    if getattr(args, "objective", None) == "max-sensitivity-at-specificity" and args.min_specificity is None:
        parser.error("--objective max-sensitivity-at-specificity requires --min-specificity")
    if getattr(args, "caller", "kestrel") in {"advntr", "both"} and getattr(args, "advntr_executable", None) is None:
        parser.error("--caller advntr and --caller both require --advntr-executable")


def _target_arguments(args: argparse.Namespace, parser: argparse.ArgumentParser, target: str, operation: str) -> None:
    """Reject incompatible target arguments before output or scientific I/O."""
    objectives = {"dominance": "lexicographic-safety-v1", "callers": "caller-safety-v1", "length": "length-total-v1"}
    if target not in objectives:
        parser.error("unsupported calibration target")
    if operation == "fit" and getattr(args, "objective", None) != objectives[target]:
        parser.error(f"calibration target {target} requires --objective {objectives[target]}")
    if operation == "extract":
        _extraction_arguments(args, parser, target)
    prior = any(getattr(args, name, None) is not None for name in ("validation", "authority"))
    if prior and (target == "dominance" or operation == "validate"):
        parser.error("prior validation and authority are accepted only by target locked evaluation or export")
    ledger = getattr(args, "exposure_ledger", None)
    custody = getattr(args, "custody", None)
    if target == "dominance":
        if ledger is not None or custody is not None:
            parser.error("target-aware exposure and custody arguments require --target callers or length")
    else:
        if operation not in {"extract", "export"} and ledger is None:
            parser.error("target calibration requires --exposure-ledger")
        if operation in {"validate", "evaluate"} and custody is None:
            parser.error("target confirmation requires --custody")
        if operation == "evaluate" and (
            getattr(args, "validation", None) is None or getattr(args, "authority", None) is None
        ):
            parser.error("locked target evaluation requires --validation and --authority")


def _extraction_arguments(args: argparse.Namespace, parser: argparse.ArgumentParser, target: str) -> None:
    if target == "dominance":
        if getattr(args, "truth", None) is None or getattr(args, "partitions", None) is None:
            parser.error("dominance extraction requires --truth and --partitions")
        if any(getattr(args, key, None) is not None for key in ("study", "sources", "length_annotation")):
            parser.error("target study/source arguments require --target callers or length")
    else:
        if getattr(args, "study", None) is None or getattr(args, "sources", None) is None:
            parser.error("target extraction requires --study and --sources")
        if getattr(args, "truth", None) is not None or getattr(args, "partitions", None) is not None:
            parser.error("target extraction uses sealed source metadata instead of dominance truth/partitions")
        annotation = getattr(args, "length_annotation", None)
        if (target == "length") != (annotation is not None):
            parser.error("--length-annotation is required exactly for length extraction")


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


def _fit_length(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_length_controller import fit_length_bundle

    return fit_length_bundle(args, output)


def _assess_length(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_length_controller import assess_length_bundle

    return assess_length_bundle(args, output)


def _fit_callers(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_caller_controller import fit_caller_bundle

    return fit_caller_bundle(args, output)


def _assess_callers(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_caller_controller import assess_caller_bundle

    return assess_caller_bundle(args, output)


def _extract_target(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_target_evidence import extract_target_evidence

    return extract_target_evidence(
        args.study,
        args.runs,
        args.sources,
        output,
        expected_target=args.target,
        length_annotation_path=args.length_annotation,
    )


def _validate_target(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_confirmation_controller import confirm_calibration_bundle

    return confirm_calibration_bundle(args, output, role="validation")


def _evaluate_target(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_confirmation_controller import confirm_calibration_bundle

    return confirm_calibration_bundle(args, output, role="locked-heldout")


def _export_target(args: argparse.Namespace, output: Path) -> bool:
    from vntyper.scripts.calibration_export import export_calibration_bundle

    return export_calibration_bundle(args, output)


TARGET_OPERATIONS: dict[tuple[str, str], CalibrationOperation] = {
    ("length", "fit"): _fit_length,
    ("length", "assess"): _assess_length,
    ("callers", "fit"): _fit_callers,
    ("callers", "assess"): _assess_callers,
    **{
        (target, operation): handler
        for target in ("callers", "length")
        for operation, handler in (
            ("extract", _extract_target),
            ("validate", _validate_target),
            ("evaluate", _evaluate_target),
            ("export", _export_target),
        )
    },
}
