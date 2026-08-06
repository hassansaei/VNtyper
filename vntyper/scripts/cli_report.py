"""
cli_report.py

Module Purpose:
---------------
The ``vntyper report`` subcommand handler, extracted out of ``cli.py``'s
``main()``.

It lives in its own module rather than beside the other handlers in
:mod:`vntyper.scripts.cli_handlers` because it was the one subcommand known to be
broken -- ``generate_summary_report()`` did not accept the arguments this call
site passed (AGENTS.md trap 11) -- and fixing it was owned separately from every
other subcommand. Keeping it here means that fix touched no other handler.

The handler's job is to fill in what the CLI left unset and then make one call.
Everything it fills in comes from one of three places, in order: an explicit
command-line value, a standard path underneath ``--input-dir``, or
``config["default_values"]``. It deliberately does **not** compute anything the
report generator can read for itself -- the input file list, the pipeline version
and the VNTR coverage all live in ``pipeline_summary.json`` and are read there.

Logging is not configured here; ``cli.py`` remains the sole place that calls
``setup_logging()``.
"""

import argparse
import logging
from pathlib import Path
from typing import Any

from vntyper.scripts.cli_handlers import get_conf
from vntyper.scripts.generate_report import generate_summary_report

logger = logging.getLogger(__name__)

#: Where a finished pipeline run leaves its log, relative to its output directory.
PIPELINE_LOG_NAME = "pipeline.log"


def resolve_log_file(output_dir: str, log_file_str: str | None, explicit: bool) -> str | None:
    """Choose the pipeline log the report should embed.

    ``vntyper report -o out/`` regenerates a report from a run that has already
    finished, and that run's log is ``out/pipeline.log``. The log file ``cli.py``
    resolves is the one *this* process writes to, which without an explicit
    ``--log-file`` is a bare relative default from ``config["cli_defaults"]`` --
    almost never the finished run's log.

    Args:
        output_dir: The ``--output-dir`` value.
        log_file_str: The log file ``cli.py`` resolved, or None.
        explicit: Whether the user passed ``--log-file`` themselves.

    Returns:
        str | None: Path to the log to embed, or None when there is none.
    """
    if not explicit:
        candidate = Path(output_dir) / PIPELINE_LOG_NAME
        if candidate.exists():
            logger.debug(f"log_file set to {candidate} (found in the output directory)")
            return str(candidate)
    return log_file_str


def handle_report(
    args: argparse.Namespace,
    config: dict[str, Any],
    parser: argparse.ArgumentParser,
    log_level_value: int,
    log_file_str: str | None,
) -> None:
    """Run the ``report`` subcommand.

    Args:
        args: The parsed arguments.
        config: The loaded configuration mapping.
        parser: Unused; present for the uniform handler signature.
        log_level_value: Unused; present for the uniform handler signature.
        log_file_str: The log file ``cli.py`` resolved, forwarded to the report
            generator so the report can embed it.
    """
    # No need to set up logging here; it's already set up in cli.py

    # Fallback from config if not specified on CLI
    # (we do so only if relevant keys exist in config["default_values"])
    if args.report_file is None:
        args.report_file = get_conf(config, "report_file", "summary_report.html")
        logger.debug(f"report_file set to {args.report_file}")
    if args.flanking is None:
        args.flanking = get_conf(config, "flanking", 50)
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

    # The pipeline log to embed. `--log-file` wins; otherwise prefer the finished
    # run's own log inside the output directory.
    log_file = resolve_log_file(args.output_dir, log_file_str, explicit=bool(getattr(args, "log_file", None)))

    # Now call generate_summary_report.
    #
    # `input_files`, `pipeline_version` and `mean_vntr_coverage` are deliberately
    # absent: the generator reads all three out of `pipeline_summary.json`, so
    # passing them here would have been a second, divergent source for the same
    # facts even if the signature had accepted them. `log_file` is required and is
    # what makes the Pipeline Log section of the report render.
    generate_summary_report(
        output_dir=Path(args.output_dir),
        template_dir=config.get("paths", {}).get("template_dir", "vntyper/templates"),
        report_file=args.report_file,
        log_file=log_file,
        bed_file=args.bed_file,
        bam_file=args.bam_file,
        fasta_file=args.reference_fasta,
        flanking=args.flanking,
        vcf_file=None,
        config=config,
    )
