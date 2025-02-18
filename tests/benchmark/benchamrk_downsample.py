#!/usr/bin/env python3
"""
Downsampling and Benchmarking Script for BAM Files

This script subsets a BAM file to the MUC1 region, calculates coverage in the VNTR region,
and downsamples the BAM to specified fractions or absolute coverage levels. Optionally,
it can run `vntyper` on the downsampled BAM files and summarize the results in a tabular format.
When using the new adVNTR benchmarking feature (via --run-advntr), the vntyper command
will be modified to include the extra module and the resulting adVNTR output will be parsed
and added to the summary.

Usage:
    python downsample_bam.py \
        --input-bam path/to/input.bam \
        --output-dir path/to/output_directory \
        --fractions 0.1 0.2 0.5 \
        --coverages 10 20 30 \
        --seed 42 \
        --threads 4 \
        --muc1-region 'chr1:155158000-155163000' \
        --vntr-region 'chr1:155160500-155162000' \
        [--run-vntyper] \
        [--run-advntr] \
        [--vntyper-path path/to/vntyper] \
        [--vntyper-options "additional vntyper options"] \
        [--reference-assembly hg19] \
        [--keep-intermediates] \
        [--archive-results] \
        [--fast-mode] \
        [--summary-output summary.csv]
"""

import argparse
import hashlib
import logging
import os
import re
import subprocess
from pathlib import Path
from typing import List, Optional, Dict
import csv
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def run_command(command: List[str], description: str, critical: bool = True) -> None:
    """
    Run a shell command and handle errors.

    Args:
        command (List[str]): Command and arguments to execute.
        description (str): Description of the command for logging.
        critical (bool): If True, exit on failure. Otherwise, log error.
    """
    logging.info(f"Running: {description}")
    logging.debug(f"Command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to {description}. Error: {e}")
        if critical:
            raise


def compute_md5(file_path: Path, chunk_size: int = 1_048_576) -> str:
    """
    Compute MD5 checksum of a file.

    Args:
        file_path (Path): Path to the file.
        chunk_size (int): Size of chunks to read.

    Returns:
        str: MD5 checksum.
    """
    hash_md5 = hashlib.md5()
    logging.debug(f"Computing MD5 for {file_path}")
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def subset_bam(
    samtools: str,
    input_bam: Path,
    output_bam: Path,
    region: str,
    fraction: float,
    seed: int,
    threads: int,
) -> None:
    """
    Subset and optionally downsample BAM file to a specific region and fraction.

    Args:
        samtools (str): Path to samtools executable.
        input_bam (Path): Input BAM file.
        output_bam (Path): Output subset BAM file.
        region (str): Genomic region (e.g., 'chr1:155158000-155163000').
        fraction (float): Fraction of reads to retain (1.0 for full subset).
        seed (int): Seed for reproducibility.
        threads (int): Number of threads to use.
    """
    if fraction < 1.0:
        logging.info(f"Subsampling BAM with fraction {fraction}")
        # Samtools view uses -s seed.fraction where seed is integer and fraction is float between 0-1
        # Combine seed and fraction into a single float: seed.fraction
        subsample_param = float(f"{seed}.{int(fraction * 1000):03d}")  # e.g., 42.100
        view_command = [
            samtools,
            "view",
            "-s",
            f"{subsample_param}",
            "-@",
            str(threads),
            "-b",
            "-o",
            str(output_bam),
            str(input_bam),
            region,
        ]
        run_command(view_command, description=f"subsample BAM to fraction {fraction}")
    else:
        logging.info(f"Subsetting BAM to region {region} without subsampling")
        view_command = [
            samtools,
            "view",
            "-@",
            str(threads),
            "-b",
            "-o",
            str(output_bam),
            str(input_bam),
            region,
        ]
        run_command(view_command, description=f"subset BAM to {region}")

    # Sort the subset/subsampled BAM
    sorted_bam = output_bam.with_suffix(".sorted.bam")
    run_command(
        [
            samtools,
            "sort",
            "-@",
            str(threads),
            "-o",
            str(sorted_bam),
            str(output_bam),
        ],
        description="sort subset/subsampled BAM",
    )
    # Replace unsorted BAM with sorted BAM
    output_bam.unlink()
    sorted_bam.rename(output_bam)
    # Index the BAM
    run_command(
        [samtools, "index", str(output_bam)],
        description="index subset/subsampled BAM",
    )


def calculate_vntr_coverage(
    samtools: str,
    bam_file: Path,
    region: str,
    threads: int,
    output_dir: Path,
    output_name: str,
) -> float:
    """
    Calculate mean coverage over the VNTR region using samtools depth.

    Args:
        samtools (str): Path to samtools executable.
        bam_file (Path): Path to the BAM file.
        region (str): Genomic region in 'chr:start-end' format.
        threads (int): Number of threads to use.
        output_dir (Path): Directory to store the coverage log.
        output_name (str): Base name for the coverage log file.

    Returns:
        float: Mean coverage over the VNTR region.
    """
    coverage_output = output_dir / f"{output_name}_vntr_coverage.txt"
    depth_command = [
        samtools,
        "depth",
        "-@",
        str(threads),
        "-r",
        region,
        str(bam_file),
    ]
    logging.info(f"Calculating VNTR coverage for {bam_file} in region {region}")
    with coverage_output.open("w") as fout:
        subprocess.run(depth_command, stdout=fout, check=True)
    # Calculate mean coverage
    coverage_values = []
    with coverage_output.open("r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                try:
                    coverage = int(parts[2])
                    coverage_values.append(coverage)
                except ValueError:
                    continue
    mean_coverage = (
        sum(coverage_values) / len(coverage_values) if coverage_values else 0
    )
    logging.info(f"Mean VNTR coverage: {mean_coverage:.2f}")
    return mean_coverage


def run_vntyper(
    vntyper_path: str,
    bam_file: Path,
    reference_assembly: str,
    threads: int,
    output_dir: Path,
    keep_intermediates: bool,
    archive_results: bool,
    fast_mode: bool,
    additional_options: Optional[str] = None,
) -> Path:
    """
    Run the vntyper command on a given BAM file.

    Args:
        vntyper_path (str): Path to the vntyper executable.
        bam_file (Path): Path to the input BAM file.
        reference_assembly (str): Reference assembly (e.g., hg19, hg38).
        threads (int): Number of threads to use.
        output_dir (Path): Directory to store vntyper outputs.
        keep_intermediates (bool): Whether to keep intermediate files.
        archive_results (bool): Whether to archive results.
        fast_mode (bool): Whether to run in fast mode.
        additional_options (str, optional): Additional command-line options for vntyper.

    Returns:
        Path: Path to the vntyper output directory.
    """
    output_subdir = output_dir / f"{bam_file.stem}_vntyper_output"
    output_subdir.mkdir(parents=True, exist_ok=True)

    command = [
        vntyper_path,
        "-l",
        "DEBUG",
        "pipeline",
        "--bam",
        str(bam_file),
        "--threads",
        str(threads),
        "--reference-assembly",
        reference_assembly,
        "-o",
        str(output_subdir),
    ]

    if keep_intermediates:
        command.append("--keep-intermediates")
    if archive_results:
        command.append("--archive-results")
    if fast_mode:
        command.append("--fast-mode")
    if additional_options:
        # Split the additional_options string into individual arguments
        command.extend(additional_options.split())

    logging.info(f"Executing vntyper for {bam_file.name}")
    run_command(command, description=f"vntyper on {bam_file.name}")

    return output_subdir


def summarize_vntr_results(vntyper_output_dir: Path) -> Optional[Dict]:
    """
    Summarize vntyper results from the output directory.

    Args:
        vntyper_output_dir (Path): Path to the vntyper output directory.

    Returns:
        dict or None: Dictionary with summary fields or None if not found.
    """
    # Search for kestrel_result.tsv anywhere under vntyper_output_dir
    result_files = list(vntyper_output_dir.rglob("kestrel_result.tsv"))
    if not result_files:
        logging.warning(f"No kestrel_result.tsv found in {vntyper_output_dir}")
        return None

    # Take the first match
    result_file = result_files[0]
    logging.info(f"Parsing Kestrel results from {result_file}")

    try:
        # Add comment='#' to skip lines starting with '#' or '##'
        df = pd.read_csv(result_file, sep="\t", comment="#")

        # Verify that required columns exist
        required_columns = [
            "Estimated_Depth_AlternateVariant",
            "Estimated_Depth_Variant_ActiveRegion",
            "Depth_Score",
            "Confidence",
        ]
        for col in required_columns:
            if col not in df.columns:
                logging.error(f"Missing expected column '{col}' in {result_file}")
                return None

        summary = {
            # Use the vntyper output directory name as a fallback
            "file_analyzed": vntyper_output_dir.name,
            "Estimated_Depth_AlternateVariant": df[
                "Estimated_Depth_AlternateVariant"
            ].tolist(),
            "Estimated_Depth_Variant_ActiveRegion": df[
                "Estimated_Depth_Variant_ActiveRegion"
            ].tolist(),
            "Depth_Score": df["Depth_Score"].tolist(),
            "Confidence": df["Confidence"].tolist(),
        }
        return summary
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse vntyper results from {result_file}: {e}")
        return None
    except Exception as e:
        logging.error(
            f"Unexpected error while parsing vntyper results from {result_file}: {e}"
        )
        return None


def summarize_advntr_results(vntyper_output_dir: Path) -> Optional[Dict]:
    """
    Summarize adVNTR results from the vntyper output directory.

    Args:
        vntyper_output_dir (Path): Path to the vntyper output directory.

    Returns:
        dict or None: Dictionary with adVNTR summary, or None if not found.
    """
    advntr_file = vntyper_output_dir / "advntr" / "output_adVNTR.vcf"
    if not advntr_file.is_file():
        logging.warning(f"adVNTR result file not found in {vntyper_output_dir}")
        return None
    try:
        with advntr_file.open("r") as f:
            lines = [line.rstrip("\n") for line in f if line.strip()]
        data_lines = [line for line in lines if not line.startswith("#")]
        if len(data_lines) == 1:
            fields = data_lines[0].split("\t")
            if fields[0].strip().lower() == "negative":
                logging.debug("adVNTR result indicates a negative outcome.")
                return {"advntr_result": "Negative"}
        header = None
        for line in lines:
            if line.startswith("#VID"):
                header = line.lstrip("#").strip().split("\t")
                break
        if header is not None:
            from io import StringIO

            csv_data = "\n".join(data_lines)
            df = pd.read_csv(StringIO(csv_data), sep="\t", header=None)
            df.columns = header
        else:
            df = pd.read_csv(advntr_file, sep="\t", comment="#")
        summary = {"advntr_result": df.to_dict(orient="records")}
        return summary
    except Exception as e:
        logging.error(f"Error parsing adVNTR results: {e}")
        return None


def parse_pipeline_log(vntyper_output_dir: Path) -> Optional[float]:
    """
    Parse the pipeline.log file to extract the analysis time in minutes.

    Args:
        vntyper_output_dir (Path): Directory where the vntyper output is located.

    Returns:
        float or None: Analysis time in minutes if found, else None.
    """
    pipeline_log = vntyper_output_dir / "pipeline.log"
    if not pipeline_log.is_file():
        logging.warning(f"pipeline.log not found in {vntyper_output_dir}")
        return None

    try:
        with pipeline_log.open("r") as f:
            lines = f.readlines()
            if not lines:
                logging.warning("pipeline.log is empty.")
                return None
            last_line = lines[-1].strip()
            logging.debug(f"Last line of pipeline.log: {last_line}")
            # Example line:
            # 2025-01-15 13:37:21,126 - root - INFO - Pipeline completed in 5.28 minutes.
            match = re.search(r"Pipeline completed in (\d+\.?\d*) minutes\.", last_line)
            if match:
                analysis_time = float(match.group(1))
                logging.info(f"Analysis time extracted: {analysis_time} minutes")
                return analysis_time
            else:
                logging.warning("Could not parse analysis time from pipeline.log.")
                return None
    except Exception as e:
        logging.error(f"Error reading pipeline.log: {e}")
        return None


def calculate_required_fraction(
    desired_coverage: int, current_coverage: float
) -> float:
    """
    Calculate the fraction of reads to keep to achieve desired coverage.

    Args:
        desired_coverage (int): Desired coverage (e.g., 10x).
        current_coverage (float): Current mean coverage.

    Returns:
        float: Fraction of reads to retain.
    """
    if current_coverage == 0:
        logging.warning("Current coverage is 0. Cannot calculate required fraction.")
        return 1.0
    fraction = desired_coverage / current_coverage
    fraction = min(max(fraction, 0.0), 1.0)  # Clamp between 0 and 1
    logging.info(
        f"Calculating fraction to achieve {desired_coverage}x coverage: {fraction:.4f}"
    )
    return fraction


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Downsample BAM to specified fractions or absolute coverage levels. Optionally run vntyper (and adVNTR via --run-advntr) on downsampled BAMs and summarize results."
    )
    parser.add_argument(
        "--input-bam",
        required=True,
        type=Path,
        help="Path to the input BAM file.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory to store downsampled BAM files.",
    )
    parser.add_argument(
        "--fractions",
        nargs="+",
        type=float,
        default=[],
        help="List of fractions to downsample (e.g., 0.1 0.2 0.5).",
    )
    parser.add_argument(
        "--coverages",
        nargs="+",
        type=int,
        default=[],
        help="List of absolute coverage levels to downsample to (e.g., 10 20 30).",
    )
    parser.add_argument(
        "--muc1-region",
        type=str,
        required=True,
        help="Genomic region for MUC1 (e.g., 'chr1:155158000-155163000').",
    )
    parser.add_argument(
        "--vntr-region",
        type=str,
        required=True,
        help="Genomic region for VNTR coverage calculation (e.g., 'chr1:155160500-155162000').",
    )
    parser.add_argument(
        "--samtools",
        type=str,
        default="samtools",
        help="Path to samtools executable. Defaults to 'samtools' in PATH.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for subsampling reproducibility.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use for samtools operations.",
    )
    # Arguments for vntyper
    parser.add_argument(
        "--run-vntyper",
        action="store_true",
        help="Flag to indicate whether to run vntyper on the downsampled BAM files.",
    )
    parser.add_argument(
        "--run-advntr",
        action="store_true",
        help="Flag to indicate whether to include adVNTR mode (adds '--extra-modules advntr' to vntyper and parses its results).",
    )
    parser.add_argument(
        "--vntyper-path",
        type=str,
        default="vntyper",
        help="Path to the vntyper executable. Defaults to 'vntyper' in PATH.",
    )
    parser.add_argument(
        "--vntyper-options",
        type=str,
        default="",
        help="Additional command-line options for vntyper.",
    )
    parser.add_argument(
        "--reference-assembly",
        type=str,
        default="hg19",
        help="Reference assembly for vntyper (e.g., hg19, hg38).",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate files generated by vntyper.",
    )
    parser.add_argument(
        "--archive-results",
        action="store_true",
        help="Archive results generated by vntyper.",
    )
    parser.add_argument(
        "--fast-mode",
        action="store_true",
        help="Run vntyper in fast mode.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default="vntyper_summary.csv",
        help="Path to the summary CSV file for vntyper results.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Validate input BAM
    if not args.input_bam.is_file():
        logging.error(f"Input BAM file does not exist: {args.input_bam}")
        exit(1)

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize summary data list
    vntyper_summary = []

    # Step 1: Subset BAM to MUC1 region without subsampling
    subset_bam_path = args.output_dir / f"{args.input_bam.stem}_MUC1_subset.bam"
    subset_bam(
        samtools=args.samtools,
        input_bam=args.input_bam,
        output_bam=subset_bam_path,
        region=args.muc1_region,
        fraction=1.0,  # Keep all reads in the MUC1 region
        seed=args.seed,
        threads=args.threads,
    )

    # Step 2: Calculate current coverage in VNTR region
    current_coverage = calculate_vntr_coverage(
        samtools=args.samtools,
        bam_file=subset_bam_path,
        region=args.vntr_region,
        threads=args.threads,
        output_dir=args.output_dir,
        output_name=subset_bam_path.stem,
    )

    # Step 3: Downsample by fractions
    for fraction in args.fractions:
        if not (0 < fraction <= 1):
            logging.warning(f"Invalid fraction {fraction}. Skipping.")
            continue
        output_bam = (
            args.output_dir
            / f"{args.input_bam.stem}_downsampled_{int(fraction*100)}p.bam"
        )
        subset_bam(
            samtools=args.samtools,
            input_bam=subset_bam_path,
            output_bam=output_bam,
            region=args.muc1_region,
            fraction=fraction,
            seed=args.seed,
            threads=args.threads,
        )
        # Calculate MD5
        md5sum = compute_md5(output_bam)
        md5_file = output_bam.with_suffix(".bam.md5")
        with md5_file.open("w") as f:
            f.write(f"{md5sum}  {output_bam.name}\n")
        logging.info(f"Created downsampled BAM: {output_bam}")
        logging.info(f"MD5 checksum: {md5sum}")

        # If vntyper is to be run, execute it
        if args.run_vntyper:
            # If adVNTR mode is requested, ensure '--extra-modules advntr' is included
            vntyper_options = args.vntyper_options
            if args.run_advntr and "--extra-modules" not in vntyper_options:
                vntyper_options = (vntyper_options + " --extra-modules advntr").strip()
            vntyper_output_dir = run_vntyper(
                vntyper_path=args.vntyper_path,
                bam_file=output_bam,
                reference_assembly=args.reference_assembly,
                threads=args.threads,
                output_dir=args.output_dir,
                keep_intermediates=args.keep_intermediates,
                archive_results=args.archive_results,
                fast_mode=args.fast_mode,
                additional_options=vntyper_options,
            )
            # Parse pipeline.log for analysis time
            analysis_time = parse_pipeline_log(vntyper_output_dir)
            # Summarize vntyper results
            summary = summarize_vntr_results(vntyper_output_dir)
            # Summarize adVNTR results if requested
            advntr_summary_dict = {}
            if args.run_advntr:
                advntr_summary_dict = summarize_advntr_results(vntyper_output_dir)
            if summary:
                vntyper_summary.append(
                    {
                        "file_analyzed": summary.get("file_analyzed", output_bam.name),
                        "method": "fraction",
                        "value": fraction,
                        "confidence": ", ".join(
                            map(str, summary.get("Confidence", []))
                        ),
                        "Estimated_Depth_AlternateVariant": ", ".join(
                            map(
                                str, summary.get("Estimated_Depth_AlternateVariant", [])
                            )
                        ),
                        "Estimated_Depth_Variant_ActiveRegion": ", ".join(
                            map(
                                str,
                                summary.get("Estimated_Depth_Variant_ActiveRegion", []),
                            )
                        ),
                        "Depth_Score": ", ".join(
                            map(str, summary.get("Depth_Score", []))
                        ),
                        "analysis_time_minutes": (
                            analysis_time if analysis_time is not None else ""
                        ),
                        "advntr_result": (
                            str(advntr_summary_dict.get("advntr_result"))
                            if advntr_summary_dict
                            else ""
                        ),
                    }
                )
            else:
                vntyper_summary.append(
                    {
                        "file_analyzed": output_bam.name,
                        "method": "fraction",
                        "value": fraction,
                        "confidence": "No Results",
                        "Estimated_Depth_AlternateVariant": "",
                        "Estimated_Depth_Variant_ActiveRegion": "",
                        "Depth_Score": "",
                        "analysis_time_minutes": "",
                        "advntr_result": (
                            str(advntr_summary_dict.get("advntr_result"))
                            if advntr_summary_dict
                            else ""
                        ),
                    }
                )

    # Step 4: Downsample to absolute coverages
    for coverage in args.coverages:
        fraction = calculate_required_fraction(
            desired_coverage=coverage, current_coverage=current_coverage
        )
        output_bam = (
            args.output_dir / f"{args.input_bam.stem}_downsampled_{coverage}x.bam"
        )
        subset_bam(
            samtools=args.samtools,
            input_bam=subset_bam_path,
            output_bam=output_bam,
            region=args.muc1_region,
            fraction=fraction,
            seed=args.seed,
            threads=args.threads,
        )
        # Calculate MD5
        md5sum = compute_md5(output_bam)
        md5_file = output_bam.with_suffix(".bam.md5")
        with md5_file.open("w") as f:
            f.write(f"{md5sum}  {output_bam.name}\n")
        logging.info(f"Created downsampled BAM: {output_bam}")
        logging.info(f"MD5 checksum: {md5sum}")

        # If vntyper is to be run, execute it
        if args.run_vntyper:
            # If adVNTR mode is requested, ensure '--extra-modules advntr' is included
            vntyper_options = args.vntyper_options
            if args.run_advntr and "--extra-modules" not in vntyper_options:
                vntyper_options = (vntyper_options + " --extra-modules advntr").strip()
            vntyper_output_dir = run_vntyper(
                vntyper_path=args.vntyper_path,
                bam_file=output_bam,
                reference_assembly=args.reference_assembly,
                threads=args.threads,
                output_dir=args.output_dir,
                keep_intermediates=args.keep_intermediates,
                archive_results=args.archive_results,
                fast_mode=args.fast_mode,
                additional_options=vntyper_options,
            )
            # Parse pipeline.log for analysis time
            analysis_time = parse_pipeline_log(vntyper_output_dir)
            # Summarize vntyper results
            summary = summarize_vntr_results(vntyper_output_dir)
            # Summarize adVNTR results if requested
            advntr_summary_dict = {}
            if args.run_advntr:
                advntr_summary_dict = summarize_advntr_results(vntyper_output_dir)
            if summary:
                vntyper_summary.append(
                    {
                        "file_analyzed": summary.get("file_analyzed", output_bam.name),
                        "method": "coverage",
                        "value": coverage,
                        "confidence": ", ".join(
                            map(str, summary.get("Confidence", []))
                        ),
                        "Estimated_Depth_AlternateVariant": ", ".join(
                            map(
                                str, summary.get("Estimated_Depth_AlternateVariant", [])
                            )
                        ),
                        "Estimated_Depth_Variant_ActiveRegion": ", ".join(
                            map(
                                str,
                                summary.get("Estimated_Depth_Variant_ActiveRegion", []),
                            )
                        ),
                        "Depth_Score": ", ".join(
                            map(str, summary.get("Depth_Score", []))
                        ),
                        "analysis_time_minutes": (
                            analysis_time if analysis_time is not None else ""
                        ),
                        "advntr_result": (
                            str(advntr_summary_dict.get("advntr_result"))
                            if advntr_summary_dict
                            else ""
                        ),
                    }
                )
            else:
                vntyper_summary.append(
                    {
                        "file_analyzed": output_bam.name,
                        "method": "coverage",
                        "value": coverage,
                        "confidence": "No Results",
                        "Estimated_Depth_AlternateVariant": "",
                        "Estimated_Depth_Variant_ActiveRegion": "",
                        "Depth_Score": "",
                        "analysis_time_minutes": "",
                        "advntr_result": (
                            str(advntr_summary_dict.get("advntr_result"))
                            if advntr_summary_dict
                            else ""
                        ),
                    }
                )

    # Step 5: Write summary table if vntyper was run
    if args.run_vntyper:
        summary_csv_path = args.summary_output
        logging.info(f"Writing vntyper summary to {summary_csv_path}")
        with summary_csv_path.open("w", newline="") as csvfile:
            fieldnames = [
                "file_analyzed",
                "method",
                "value",
                "confidence",
                "Estimated_Depth_AlternateVariant",
                "Estimated_Depth_Variant_ActiveRegion",
                "Depth_Score",
                "analysis_time_minutes",
                "advntr_result",
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in vntyper_summary:
                writer.writerow(row)
        logging.info(f"Vntyper summary successfully written to {summary_csv_path}")


if __name__ == "__main__":
    main()
