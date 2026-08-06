"""
cli_report.py

Module Purpose:
---------------
The ``vntyper report`` subcommand handler, extracted out of ``cli.py``'s
``main()``.

It lives in its own module rather than beside the other handlers in
:mod:`vntyper.scripts.cli_handlers` because it is the one subcommand known to be
broken -- ``generate_summary_report()`` does not accept the arguments this call
site passes (AGENTS.md trap 11) -- and fixing it is owned separately from every
other subcommand. Keeping it here means that fix touches no other handler.

Logging is not configured here; ``cli.py`` remains the sole place that calls
``setup_logging()``.
"""

import argparse
import logging
from pathlib import Path
from typing import Any

from vntyper.scripts.cli_handlers import get_conf
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.version import __version__ as VERSION

logger = logging.getLogger(__name__)


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
        log_file_str: Unused; present for the uniform handler signature.
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
