"""Argument-only calibration command and runtime bundle registration."""

import argparse
from pathlib import Path


def add_calibrate_subparser(subparsers: argparse._SubParsersAction) -> None:
    """Register the closed calibration command tree.

    Args:
        subparsers: Top-level VNtyper subparser collection.
    """
    calibrate = subparsers.add_parser(
        "calibrate",
        help="Extract, fit, validate, or evaluate an opt-in calibration profile.",
        conflict_handler="resolve",
    )
    operations = calibrate.add_subparsers(dest="calibration_operation", required=True)

    intake = operations.add_parser("intake", help="Audit and normalize declared local calibration inputs.")
    intake.add_argument("--manifest", type=Path, required=True)
    intake.add_argument("--output", type=Path, required=True)
    intake.add_argument("--preprocessing-priority", nargs="+", required=True, metavar="ID")
    intake.add_argument("--cram-references", type=Path, default=None)
    intake.add_argument("--temporary-directory", type=Path, default=None)

    extract = operations.add_parser("extract", help="Extract immutable replay evidence.")
    extract.add_argument("--truth", type=Path, default=None)
    extract.add_argument("--partitions", type=Path, default=None)
    extract.add_argument("--runs", type=Path, required=True)
    extract.add_argument("--output", type=Path, required=True)
    extract.add_argument("--target", choices=["dominance", "callers", "length"], default="dominance")
    extract.add_argument("--study", type=Path, default=None)
    extract.add_argument("--sources", type=Path, default=None)
    extract.add_argument("--length-annotation", type=Path, default=None)

    fit = operations.add_parser("fit", help="Fit the frozen safety-first candidate family.")
    fit.add_argument("--evidence", type=Path, required=True)
    fit.add_argument("--target", choices=["dominance", "callers", "length"], default="dominance")
    fit.add_argument(
        "--objective", required=True, choices=["lexicographic-safety-v1", "caller-safety-v1", "length-total-v1"]
    )
    fit.add_argument("--exposure-ledger", type=Path, default=None)
    fit.add_argument(
        "--advntr-executable",
        type=Path,
        default=None,
        help="Pinned installed adVNTR executable for training an exact-caller background.",
    )
    fit.add_argument("--output", type=Path, required=True)

    for name, help_text in (
        ("validate", "Validate one fixed profile on validation evidence."),
        ("evaluate", "Evaluate one fixed profile on locked held-out evidence."),
    ):
        operation = operations.add_parser(name, help=help_text)
        operation.add_argument("--profile", type=Path, required=True)
        operation.add_argument("--evidence", type=Path, required=True)
        operation.add_argument("--output", type=Path, required=True)
        operation.add_argument("--target", choices=["dominance", "callers", "length"], default="dominance")
        operation.add_argument("--exposure-ledger", type=Path, default=None)
        operation.add_argument("--custody", type=Path, default=None)
        operation.add_argument("--validation", type=Path, default=None)
        operation.add_argument("--authority", type=Path, default=None)

    assess = operations.add_parser("assess", help="Assess one fixed research candidate on development evidence.")
    assess.add_argument("--target", choices=["callers", "length"], required=True)
    assess.add_argument("--profile", type=Path, required=True)
    assess.add_argument("--intake", type=Path, required=True)
    assess.add_argument("--runs", type=Path, required=True)
    assess.add_argument("--exposure-ledger", type=Path, required=True)
    assess.add_argument("--output", type=Path, required=True)

    export = operations.add_parser("export", help="Export the approved aggregate-free portable runtime bundle.")
    export.add_argument("--target", choices=["callers", "length"], required=True)
    for option in ("profile", "validation", "evaluation", "authority", "completion", "output"):
        export.add_argument("--" + option, type=Path, required=True)


def add_pipeline_calibration_arguments(parser: argparse.ArgumentParser) -> None:
    """Register explicit research measurement and approved runtime bundle inputs."""
    parser.add_argument(
        "--calibration-bundle",
        type=Path,
        default=None,
        help="Approved portable caller calibration bundle directory.",
    )
    parser.add_argument(
        "--calibration-context",
        type=Path,
        default=None,
        help="Explicit caller applicability and frozen native capture context JSON.",
    )
    parser.add_argument(
        "--measure-vntr-length-features",
        action="store_true",
        help="Measure research VNTR length features using an explicit annotation and provenance context.",
    )
    parser.add_argument(
        "--length-model",
        type=Path,
        default=None,
        help="Approved portable length model bundle directory; implies feature measurement.",
    )
    parser.add_argument(
        "--length-annotation",
        type=Path,
        default=None,
        help="Explicit length annotation JSON for measurement-only mode.",
    )
    parser.add_argument(
        "--length-context",
        type=Path,
        default=None,
        help="Explicit length measurement provenance context JSON.",
    )
