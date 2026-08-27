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
command-line value, a standard path underneath the effective run directory (see
:func:`resolve_bam_file`, :func:`resolve_bed_file` and :func:`resolve_vcf_file`),
or ``config["default_values"]``. It deliberately does **not** compute anything the
report generator can read for itself -- the input file list, the pipeline version
and the VNTR coverage all live in ``pipeline_summary.json`` and are read there.

Logging is not configured here; ``cli.py`` remains the sole place that calls
``setup_logging()``.
"""

import argparse
import logging
from pathlib import Path
from typing import Any

from vntyper.scripts.artifact_names import select_best_vcf_file
from vntyper.scripts.cli_handlers import get_conf
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.report_assets import DEFAULT_REPORT_IGV

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


def _effective_run_dir(output_dir: str, input_dir: Path | None) -> Path:
    """The run directory whose ``kestrel/`` subdirectory the IGV inputs live under.

    ``--input-dir`` wins when given; otherwise ``--output-dir``. This is the
    directory :func:`resolve_bam_file`, :func:`resolve_bed_file` and
    :func:`resolve_vcf_file` all search under, so the fallback to
    ``--output-dir`` lives in exactly one place.

    Args:
        output_dir: The ``--output-dir`` value.
        input_dir: The ``--input-dir`` value, or None.

    Returns:
        Path: The directory to search under ``kestrel/``.
    """
    return input_dir if input_dir is not None else Path(output_dir)


def resolve_bam_file(output_dir: str, input_dir: Path | None, explicit: Path | None) -> Path | None:
    """Choose the BAM the report's IGV panel should show.

    Resolution order: an explicit ``--bam-file`` wins; otherwise the standard
    ``kestrel/output.bam`` under the effective run directory (see
    :func:`_effective_run_dir`) -- ``--input-dir`` if given, else
    ``--output-dir``.

    Widened from ``--input-dir``-only discovery (#167, REVIEW.md finding 1):
    ``generate_report.py`` only runs IGV generation when ``bed_file`` is both
    set and exists, so without the ``--output-dir`` fallback,
    ``vntyper report -o <run>`` -- the invocation ``--output-dir``'s own help
    text describes ("Directory containing pipeline results") -- resolved no
    BED at all and skipped IGV generation entirely, regardless of what the VCF
    resolved to. ``--output-dir`` is already the basis :func:`resolve_log_file`
    reads ``pipeline.log`` from, so this makes the handler consistent rather
    than introducing a new convention.

    Args:
        output_dir: The ``--output-dir`` value.
        input_dir: The ``--input-dir`` value, or None.
        explicit: The ``--bam-file`` value, or None.

    Returns:
        Path | None: Path to the BAM, or None when there is none.
    """
    if explicit is not None:
        return explicit
    candidate = _effective_run_dir(output_dir, input_dir) / "kestrel" / "output.bam"
    if candidate.exists():
        logger.debug(f"bam_file set to {candidate} (found under the run directory)")
        return candidate
    return None


def resolve_bed_file(output_dir: str, input_dir: Path | None, explicit: Path | None) -> Path | None:
    """Choose the BED the report's IGV panel should show.

    Same resolution order and rationale as :func:`resolve_bam_file`.

    Args:
        output_dir: The ``--output-dir`` value.
        input_dir: The ``--input-dir`` value, or None.
        explicit: The ``--bed-file`` value, or None.

    Returns:
        Path | None: Path to the BED, or None when there is none.
    """
    if explicit is not None:
        return explicit
    candidate = _effective_run_dir(output_dir, input_dir) / "kestrel" / "output.bed"
    if candidate.exists():
        logger.debug(f"bed_file set to {candidate} (found under the run directory)")
        return candidate
    return None


def resolve_vcf_file(output_dir: str, input_dir: Path | None, explicit: Path | None) -> str | None:
    """Choose the VCF the report's IGV panel should show.

    ``handle_report`` used to pass ``vcf_file=None`` unconditionally, so a
    regenerated report never had a variant track even though the run being
    re-reported had written the VCF (#167).
    :func:`~vntyper.scripts.artifact_names.select_best_vcf_file` is the
    resolver ``pipeline.py`` already uses -- it prefers the compressed VCF and
    returns None rather than raising when neither exists, because bcftools is
    optional. It is reused here rather than reimplemented.

    An explicit ``--vcf-file`` wins and is passed through even when it does
    not exist: the user named it, and ``run_igv_report`` warns and skips the
    track.

    Same effective-run-directory resolution as :func:`resolve_bam_file` and
    :func:`resolve_bed_file` -- see :func:`_effective_run_dir`.

    Args:
        output_dir: The ``--output-dir`` value.
        input_dir: The ``--input-dir`` value, or None.
        explicit: The ``--vcf-file`` value, or None.

    Returns:
        str | None: Path to the VCF, or None when there is none.
    """
    if explicit is not None:
        logger.debug(f"vcf_file set to {explicit} (given on the command line)")
        return str(explicit)
    kestrel_dir = _effective_run_dir(output_dir, input_dir) / "kestrel"
    return select_best_vcf_file(kestrel_dir)


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

    # BAM, BED and VCF for the IGV panel are resolved from one effective run
    # directory (`--input-dir` if given, else `--output-dir`), in the same
    # precedence for each: the explicit flag first, then that directory's
    # standard `kestrel/output.*` layout. See resolve_bam_file/resolve_bed_file/
    # resolve_vcf_file's docstrings for why `--output-dir` is consulted at all
    # (#167, REVIEW.md finding 1) -- without it, `vntyper report -o <run>` with
    # no `--input-dir` resolved no BED and therefore generated no IGV panel.
    args.bam_file = resolve_bam_file(args.output_dir, args.input_dir, args.bam_file)
    args.bed_file = resolve_bed_file(args.output_dir, args.input_dir, args.bed_file)

    # The pipeline log to embed. `--log-file` wins; otherwise prefer the finished
    # run's own log inside the output directory.
    log_file = resolve_log_file(args.output_dir, log_file_str, explicit=bool(getattr(args, "log_file", None)))

    # Now call generate_summary_report.
    #
    # `input_files`, `pipeline_version` and `mean_vntr_coverage` are deliberately
    # absent: the generator reads all three out of `pipeline_summary.json`, so
    # passing them here would have been a second, divergent source for the same
    # facts even if the signature had accepted them. `log_file` is required and is
    # what makes the Pipeline Log section of the report render. `vcf_file` is
    # resolved the same way as bam/bed above, reusing the `select_best_vcf_file`
    # resolver `pipeline.py` already uses (#167); a missing VCF stays a warning,
    # never an error, because bcftools is optional. `getattr` guards a direct
    # `handle_report()` call whose namespace predates `--vcf-file`.
    #
    # `sample_name` is the one identity value that genuinely cannot be read out of
    # the summary: it is what the operator wants this report *called*, and the
    # summary records only the input file names. Passing None -- which is what the
    # pipeline does -- makes the generator derive it from those names instead, so
    # this is an override rather than a second source. `getattr` guards a direct
    # `handle_report()` call whose namespace predates `--sample-name`.
    generate_summary_report(
        output_dir=Path(args.output_dir),
        template_dir=config.get("paths", {}).get("template_dir"),
        report_file=args.report_file,
        log_file=log_file,
        bed_file=args.bed_file,
        bam_file=args.bam_file,
        fasta_file=args.reference_fasta,
        flanking=args.flanking,
        vcf_file=resolve_vcf_file(args.output_dir, args.input_dir, getattr(args, "vcf_file", None)),
        config=config,
        sample_name=getattr(args, "sample_name", None),
        # How the report carries its alignment browser. `getattr` guards a direct
        # `handle_report()` call whose namespace predates `--report-igv`, the same way
        # `--vcf-file` and `--sample-name` are guarded above.
        report_igv=getattr(args, "report_igv", DEFAULT_REPORT_IGV),
    )
