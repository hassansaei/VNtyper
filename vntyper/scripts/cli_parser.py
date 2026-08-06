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
    build_parser: Build the full VNtyper argument parser, subcommands included.
"""

import argparse
import logging
from pathlib import Path

from vntyper.scripts.reference_registry import list_assemblies
from vntyper.version import __version__ as VERSION

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """
    Build the VNtyper argument parser.

    A new parser is constructed on every call rather than cached at module level:
    ``argparse`` parsers carry mutable state, so handing out a shared instance would
    let one caller's parse leak into the next.

    Returns:
        argparse.ArgumentParser: The top-level parser, with the ``pipeline``,
        ``report``, ``cohort``, ``install-references`` and ``online`` subparsers
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
    # Changed here (#57): allow --extra-modules multiple times
    parser_pipeline.add_argument(
        "--extra-modules",
        action="append",
        default=[],
        help="Optional extra modules to include (e.g., advntr, shark). Can be repeated multiple times.",
    )
    parser_pipeline.add_argument("--fastq1", type=str, help="Path to the first FASTQ file.")
    parser_pipeline.add_argument("--fastq2", type=str, help="Path to the second FASTQ file.")
    parser_pipeline.add_argument("--bam", type=str, help="Path to the BAM file.")
    parser_pipeline.add_argument("--cram", type=str, help="Path to the CRAM file.")
    parser_pipeline.add_argument("--threads", type=int, default=None, help="Number of threads to use.")
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
        help="Keep intermediate files (e.g., BAM slices, temporary files).",
    )
    parser_pipeline.add_argument(
        "--delete-intermediates",
        action="store_true",
        help="Delete intermediate files after processing (overrides --keep-intermediates).",
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
        help="Base name for the output files.",
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
    parser_report.add_argument("--report-file", type=str, default=None, help="Name of the output report file.")
    parser_report.add_argument("--bed-file", type=Path, help="Path to the BED file for IGV reports.")
    parser_report.add_argument("--bam-file", type=Path, help="Path to the BAM file for IGV reports.")
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
        help="Name of the cohort summary report file.",
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
            "Pseudonymize sample names to protect sensitive information. "
            "Optionally provide a basename for pseudonyms (default is 'sample_')."
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
        type=int,
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
    parser_online.add_argument("--threads", type=int, default=None, help="Number of threads to use.")
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
