"""
cli_parser.py

Argument-parser construction for the VNtyper CLI.

This module holds nothing but the ``argparse`` tree. It was extracted verbatim from
``vntyper/cli.py``, where it lived inline inside ``main()``: that made the parser
unobtainable without running the program, and left ``cli.py`` at 0% unit coverage while
sitting over the repository's 650-line limit (see "Changing existing code" in
AGENTS.md). ``cli.py`` now imports ``build_parser`` and keeps the dispatch logic.

Nothing here executes a subcommand, loads configuration or configures logging.
Logging setup stays in ``cli.py``, which remains the sole place it happens.

Functions:
    positive_int: Parse a positive integer command-line token.
    build_parser: Build the full VNtyper argument parser, subcommands included.
"""

import argparse
import logging
import math
from pathlib import Path

from vntyper.scripts.reference_registry import list_assemblies
from vntyper.scripts.report_assets import DEFAULT_REPORT_IGV, REPORT_IGV_MODES
from vntyper.version import __version__ as VERSION

logger = logging.getLogger(__name__)


def positive_int(value: str) -> int:
    """Parse a positive integer command-line token.

    Args:
        value: Raw command-line token.

    Returns:
        int: The parsed value, guaranteed to be at least one.

    Raises:
        argparse.ArgumentTypeError: If the token is not an integer or is less
            than one. argparse reports this as a usage error with exit code 2.
    """
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"expected an integer, got {value!r}") from error
    if number < 1:
        raise argparse.ArgumentTypeError(f"expected a thread count of at least 1, got {number}")
    return number


def unit_fraction(value: str) -> float:
    """Parse a unit fraction command-line token between 0 and 1 inclusive.

    Args:
        value: Raw command-line token.

    Returns:
        float: The parsed fraction in [0.0, 1.0].

    Raises:
        argparse.ArgumentTypeError: If the token is not numeric or outside [0, 1].
            argparse reports this as a usage error with exit code 2.
    """
    try:
        number = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"expected a numeric fraction from 0 to 1, got {value!r}") from error
    if not (0.0 <= number <= 1.0) or not math.isfinite(number):
        raise argparse.ArgumentTypeError(f"expected a fraction between 0 and 1, got {number}")
    return number


#: ``--report-igv``'s choices and default, taken from the module that implements them
#: rather than restated here. Two spellings of the same list is how a mode ends up
#: accepted by the parser and rejected by the generator.
REPORT_IGV_CHOICES = REPORT_IGV_MODES
REPORT_IGV_DEFAULT = DEFAULT_REPORT_IGV

#: One help string, used by both subcommands. `pipeline` and `report` produce the same
#: artifact by different routes, so an option that read differently on the two would be
#: describing a difference that does not exist.
REPORT_IGV_HELP = (
    "How the report carries its alignment browser. 'embedded' (default) writes the "
    "vendored, gzipped igv.js into the report, so the archived file is a complete "
    "alignment browser needing neither a second file nor a network - about 500 KB. "
    "'sidecar' leaves it out and points the reader at the self-contained "
    "igv_report.html written beside it. 'off' produces no alignment browser at all."
)


def add_calibrate_subparser(subparsers: argparse._SubParsersAction) -> None:
    """Register the closed four-operation calibration command tree.

    Args:
        subparsers: Top-level VNtyper subparser collection.
    """
    calibrate = subparsers.add_parser(
        "calibrate",
        help="Extract, fit, validate, or evaluate an opt-in calibration profile.",
        conflict_handler="resolve",
    )
    operations = calibrate.add_subparsers(dest="calibration_operation", required=True)

    extract = operations.add_parser("extract", help="Extract immutable replay evidence.")
    extract.add_argument("--truth", type=Path, required=True)
    extract.add_argument("--partitions", type=Path, required=True)
    extract.add_argument("--runs", type=Path, required=True)
    extract.add_argument("--output", type=Path, required=True)

    fit = operations.add_parser("fit", help="Fit the frozen safety-first candidate family.")
    fit.add_argument("--evidence", type=Path, required=True)
    fit.add_argument("--objective", required=True, choices=["lexicographic-safety-v1"])
    fit.add_argument("--output", type=Path, required=True)

    for name, help_text in (
        ("validate", "Validate one fixed profile on validation evidence."),
        ("evaluate", "Evaluate one fixed profile on locked held-out evidence."),
    ):
        operation = operations.add_parser(name, help=help_text)
        operation.add_argument("--profile", type=Path, required=True)
        operation.add_argument("--evidence", type=Path, required=True)
        operation.add_argument("--output", type=Path, required=True)


def build_parser() -> argparse.ArgumentParser:
    """
    Build the VNtyper argument parser.

    A new parser is constructed on every call rather than cached at module level:
    ``argparse`` parsers carry mutable state, so handing out a shared instance would
    let one caller's parse leak into the next.

    Returns:
        argparse.ArgumentParser: The top-level parser, with the ``pipeline``,
        ``report``, ``cohort``, ``install-references``, ``online`` and ``calibrate`` subparsers
        registered under ``dest="command"``.
    """
    # Parent parser for global arguments
    parent_parser = argparse.ArgumentParser(add_help=False, conflict_handler="resolve")
    parent_parser.add_argument(
        "-l",
        "--log-level",
        help="Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR)",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    parent_parser.add_argument("-f", "--log-file", help="Set the log output file (default is stdout)")
    parent_parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    parent_parser.add_argument(
        "--config-path",
        type=Path,
        default=None,
        help=("Path to the configuration file (config.json). If not provided, the default config will be used."),
        required=False,
    )

    # Main parser that includes the parent parser
    parser = argparse.ArgumentParser(
        description="VNtyper CLI: A pipeline for genotyping MUC1-VNTR.",
        parents=[parent_parser],
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    add_calibrate_subparser(subparsers)

    # Subcommand: pipeline
    parser_pipeline = subparsers.add_parser(
        "pipeline", help="Run the full VNtyper pipeline.", conflict_handler="resolve"
    )
    parser_pipeline.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for the results.",
    )
    parser_pipeline.add_argument(
        "--decision-profile",
        type=Path,
        default=None,
        help=(
            "Path to one complete explicit decision profile. Omit it to use the verified packaged profile; "
            "profiles are never overlaid, discovered, or selected automatically."
        ),
    )
    # Changed here (#57): allow --extra-modules multiple times
    parser_pipeline.add_argument(
        "--extra-modules",
        action="append",
        default=[],
        help=(
            "Optional extra modules to include: advntr, shark. Can be repeated, or given as a "
            "comma-separated list ('--extra-modules advntr,shark'). Unknown names are rejected."
        ),
    )
    parser_pipeline.add_argument("--fastq1", type=str, help="Path to the first FASTQ file.")
    parser_pipeline.add_argument("--fastq2", type=str, help="Path to the second FASTQ file.")
    parser_pipeline.add_argument("--bam", type=str, help="Path to the BAM file.")
    parser_pipeline.add_argument("--cram", type=str, help="Path to the CRAM file.")
    parser_pipeline.add_argument("--reference-fasta", type=Path, help="Path to the reference FASTA for CRAM decoding.")
    parser_pipeline.add_argument(
        "--threads",
        type=positive_int,
        default=None,
        help="Number of threads to use. Applies to alignment, samtools, fastp and adVNTR. "
        "adVNTR honours it from 2.0.0 onward, where -t parallelises read decoding; "
        "before that its -t flag was a no-op on the short-read frameshift path.",
    )
    parser_pipeline.add_argument(
        "--reference-assembly",
        type=str,
        choices=list_assemblies(include_deprecated=False),  # Only canonical names
        default=None,
        help="Reference assembly for BAM/CRAM alignment. "
        "Options: hg19, hg38 (UCSC), hg19_ncbi, hg38_ncbi (NCBI RefSeq), hg19_ensembl, hg38_ensembl (ENSEMBL).",
    )
    parser_pipeline.add_argument(
        "--fast-mode",
        action="store_true",
        help="Enable fast mode (skips filtering for unmapped and partially mapped reads).",
    )
    parser_pipeline.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Compatibility flag: intermediate files (BAM slices, temporary files) are already kept by default, "
        "so this flag changes nothing. Use --delete-intermediates to remove them.",
    )
    parser_pipeline.add_argument(
        "--delete-intermediates",
        action="store_true",
        help="Delete intermediate files after processing (wins when --keep-intermediates is also given).",
    )
    parser_pipeline.add_argument(
        "--archive-results",
        action="store_true",
        help="Create an archive of the results folder after pipeline completion.",
    )
    parser_pipeline.add_argument(
        "--archive-format",
        type=str,
        choices=["zip", "tar.gz"],
        default=None,
        help="Format of the archive: 'zip' or 'tar.gz'.",
    )
    parser_pipeline.add_argument(
        "-n",
        "--output-name",
        type=str,
        default=None,
        help=(
            "Base name for the output files. Currently fixed at 'output': the report "
            "generator, the report subcommand and the Kestrel stage each name their "
            "files from that literal and take no basename argument, so any other value "
            "is rejected rather than silently ignored. Use --output-dir to separate runs."
        ),
    )
    parser_pipeline.add_argument(
        "-s",
        "--sample-name",
        type=str,
        default=None,
        help=("Set the sample name for labeling results. If not provided, defaults to input BAM or FASTQ name."),
    )
    region_group = parser_pipeline.add_mutually_exclusive_group()
    region_group.add_argument(
        "--custom-regions",
        type=str,
        help="Define custom regions for MUC1 analysis as comma-separated values (e.g., chr1:1000-2000,chr2:3000-4000).",
    )
    region_group.add_argument(
        "--bed-file",
        type=Path,
        help="Path to a BED file specifying regions for MUC1 analysis.",
    )
    parser_pipeline.add_argument(
        "--advntr-max-coverage",
        type=int,
        default=None,
        help="Max coverage (e.g. 300) for quick adVNTR mode.",
    )
    # New argument: additional summary output formats (comma-separated list)
    parser_pipeline.add_argument(
        "--summary-formats",
        type=str,
        default="",
        help=(
            "Comma-separated list of additional summary output formats to generate "
            "(supported: csv, tsv). JSON is always generated."
        ),
    )
    # Declared on `pipeline` as well as on `report` because an ordinary run never goes
    # through the `report` subcommand: `handle_pipeline` calls `run_pipeline`, which
    # calls `generate_summary_report` itself. An option only `report` accepted would be
    # unreachable for every run that produces a report in the first place (#242).
    parser_pipeline.add_argument(
        "--report-igv",
        type=str,
        default=REPORT_IGV_DEFAULT,
        choices=list(REPORT_IGV_CHOICES),
        help=REPORT_IGV_HELP,
    )

    # Subcommand: report
    parser_report = subparsers.add_parser(
        "report",
        help="Generate a summary report and visualizations from output data.",
        conflict_handler="resolve",
    )
    parser_report.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory containing pipeline results (subdirectories, etc.).",
    )
    parser_report.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="If provided, search this directory (and its subdirs) for standard pipeline output filenames.",
    )
    parser_report.add_argument(
        "--report-file",
        type=str,
        default=None,
        help="Name of the output report file: a single file name, written inside --output-dir. A path is refused.",
    )
    parser_report.add_argument("--bed-file", type=Path, help="Path to the BED file for IGV reports.")
    parser_report.add_argument("--bam-file", type=Path, help="Path to the BAM file for IGV reports.")
    parser_report.add_argument(
        "--vcf-file",
        type=Path,
        default=None,
        help="Path to the VCF file for the IGV variant track. Discovered under the run directory when omitted.",
    )
    parser_report.add_argument(
        "--reference-fasta",
        type=Path,
        help="Path to the reference FASTA file for IGV reports.",
    )
    parser_report.add_argument(
        "--flanking",
        type=int,
        default=None,
        help="Flanking region size for IGV reports.",
    )
    parser_report.add_argument(
        "-s",
        "--sample-name",
        type=str,
        default=None,
        help=(
            "What the report calls its sample, in the title, the heading and the header block. "
            "Derived from the run's own input file names when omitted."
        ),
    )
    parser_report.add_argument(
        "--report-igv",
        type=str,
        default=REPORT_IGV_DEFAULT,
        choices=list(REPORT_IGV_CHOICES),
        help=REPORT_IGV_HELP,
    )

    # Subcommand: cohort
    parser_cohort = subparsers.add_parser(
        "cohort",
        help="Aggregate outputs from multiple runs into a single summary file.",
        conflict_handler="resolve",
    )
    cohort_group = parser_cohort.add_mutually_exclusive_group(required=True)
    cohort_group.add_argument(
        "-i",
        "--input-dirs",
        nargs="+",
        help="List of directories containing output files to aggregate.",
    )
    cohort_group.add_argument(
        "--input-file",
        type=Path,
        help="Path to a newline-separated text file listing directories or zip files to aggregate.",
    )
    parser_cohort.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for the aggregated summary.",
    )
    parser_cohort.add_argument(
        "--summary-file",
        type=str,
        default=None,
        help=(
            "Name of the cohort summary report file: a single file name, written inside "
            "--output-dir. A path is refused."
        ),
    )
    # New argument for cohort: additional output formats (comma-separated)
    parser_cohort.add_argument(
        "--summary-formats",
        type=str,
        default="",
        help=(
            "Comma-separated list of additional summary output formats to generate "
            "(supported: csv, tsv, json). HTML is always generated."
        ),
    )
    # New flag: pseudonymize sample names in cohort mode.
    # This argument optionally accepts a basename; if not provided, the default "sample_" is used.
    parser_cohort.add_argument(
        "--pseudonymize-samples",
        nargs="?",
        const="sample_",
        default=None,
        help=(
            "Replace sample names with a prefix plus a truncated SHA-256 digest of the name. "
            "Optionally provide a basename for pseudonyms (default is 'sample_'). This is "
            "obfuscation for readability, not a privacy control: the digest is unsalted and "
            "unkeyed, so where sample names are drawn from a guessable set anyone holding the "
            "report can recover them by testing candidates at one hash each. Treat a shared "
            "report as identifying."
        ),
    )
    parser_cohort.add_argument(
        "--rare-allele-max-frequency",
        type=unit_fraction,
        default=None,
        help=(
            "Maximum frequency threshold (0 to 1) for marking calls in the cohort call "
            "frequency table as Below_Cutoff (default: 0.05). Calls with frequency above "
            "this threshold are not omitted."
        ),
    )

    # Subcommand: install-references
    parser_install = subparsers.add_parser(
        "install-references",
        help="Download and set up necessary reference files.",
        conflict_handler="resolve",
    )
    parser_install.add_argument(
        "-d",
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where references will be installed.",
    )
    parser_install.add_argument("--skip-indexing", action="store_true", help="Skip the indexing step.")
    parser_install.add_argument(
        "-t",
        "--threads",
        type=positive_int,
        default=4,
        help="Number of threads to use for indexing (default: 4).",
    )
    parser_install.add_argument(
        "--aligners",
        nargs="+",
        default=None,
        metavar="ALIGNER",
        help=(
            "Specific aligners to use (e.g., bwa bwa-mem2 minimap2). If not specified, only BWA will be used (default)."
        ),
    )
    parser_install.add_argument(
        "--references",
        nargs="+",
        default=None,
        metavar="REFERENCE",
        help=(
            "Specific references to process (e.g., hg19 hg38 GRCh37 GRCh38). Default: hg19 hg38 (UCSC references only)."
        ),
    )
    parser_install.add_argument(
        "--derive-only",
        action="store_true",
        help=(
            "Rebuild only the derived reference files from what is already installed, without "
            "downloading anything. Three files are derived rather than downloaded: the two MUC1 "
            "region FASTAs, cut from an installed chromosome with samtools faidx, and the merged "
            "MUC1 motif FASTA built from two seeds. Each is verified against its committed "
            "checksum. Use this when a tree has its genomes and seeds but is missing a derived "
            "file: a full --from-source run would re-download and re-index six chromosomes to "
            "rebuild three small ones."
        ),
    )
    parser_install.add_argument(
        "--from-source",
        action="store_true",
        help=(
            "Build references from their upstream sources instead of downloading the "
            "published bundle. Slower (downloads and BWA-indexes six chromosome FASTAs) "
            "and it needs four seed files (MUC1_motifs_Rev_com.fa, code-adVNTR_RUs.fa, "
            "vntr_db_advntr.zip, filter_config.json), fetched automatically into the "
            "output directory from berntpopp/vntyper-data unless already staged there - "
            "which is what the bundle build workflow does before running this same path."
        ),
    )
    parser_install.add_argument(
        "--release-spec",
        type=Path,
        default=None,
        metavar="PATH",
        help=("Release builds only: take every source URL and digest from this file. Requires --from-source."),
    )

    # Subcommand: online
    parser_online = subparsers.add_parser(
        "online",
        help=("Subset the BAM and submit it to an online vntyper instance, then retrieve results."),
        conflict_handler="resolve",
    )
    parser_online.add_argument("--bam", type=str, required=True, help="Path to the input BAM file.")
    parser_online.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for results.",
    )
    parser_online.add_argument(
        "--reference-assembly",
        type=str,
        choices=list_assemblies(include_deprecated=False),  # Only canonical names
        default=None,
        help="Reference assembly used (hg19, hg38, GRCh37, GRCh38, hg19_ncbi, hg38_ncbi, hg19_ensembl, hg38_ensembl).",
    )
    parser_online.add_argument(
        "--threads",
        type=positive_int,
        default=None,
        help="Number of threads to use. Applies to alignment, samtools, fastp and adVNTR. "
        "adVNTR honours it from 2.0.0 onward, where -t parallelises read decoding; "
        "before that its -t flag was a no-op on the short-read frameshift path.",
    )
    parser_online.add_argument(
        "--email",
        type=str,
        default=None,
        help="Email to receive notifications (optional).",
    )
    parser_online.add_argument(
        "--cohort-id",
        type=str,
        default=None,
        help="Cohort ID to associate the job with (optional).",
    )
    parser_online.add_argument(
        "--passphrase",
        type=str,
        default=None,
        help="Passphrase for the cohort (if required).",
    )
    parser_online.add_argument(
        "--resume",
        action="store_true",
        help="Resume polling a previously submitted job if job_id is found.",
    )

    return parser
