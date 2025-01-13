#!/usr/bin/env python3
"""
Downsampling and Benchmarking Script for BAM Files

This script subsets a BAM file to the MUC1 region, calculates coverage in the VNTR region,
and downsamples the BAM to specified fractions or absolute coverage levels. The downsampled
BAM files are sorted, indexed, and accompanied by MD5 checksums for verification.

Usage:
    python downsample_bam.py \
        --input-bam path/to/input.bam \
        --output-dir path/to/output_directory \
        --fractions 0.1 0.2 0.5 \
        --coverages 10 20 30 \
        --seed 42 \
        --threads 4 \
        --muc1-region 'chr1:155158000-155163000' \
        --vntr-region 'chr1:155160500-155162000'
"""

import argparse
import hashlib
import logging
import os
import subprocess
from pathlib import Path
from typing import List, Optional


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
    mean_coverage = sum(coverage_values) / len(coverage_values) if coverage_values else 0
    logging.info(f"Mean VNTR coverage: {mean_coverage:.2f}")
    return mean_coverage


def subsample_bam(
    samtools: str,
    input_bam: Path,
    output_bam: Path,
    fraction: float,
    seed: int,
    threads: int,
) -> None:
    """
    Subsample BAM file to a specific fraction.

    Args:
        samtools (str): Path to samtools executable.
        input_bam (Path): Input BAM file.
        output_bam (Path): Output subsampled BAM file.
        fraction (float): Fraction of reads to retain.
        seed (int): Seed for reproducibility.
        threads (int): Number of threads to use.
    """
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
    ]
    run_command(view_command, description=f"subsample BAM to fraction {fraction}")

    # Sort the subsampled BAM
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
        description="sort subsampled BAM",
    )
    # Replace unsorted BAM with sorted BAM
    output_bam.unlink()
    sorted_bam.rename(output_bam)
    # Index the BAM
    run_command(
        [samtools, "index", str(output_bam)],
        description="index subsampled BAM",
    )


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
        description="Downsample BAM to specified fractions or absolute coverage levels."
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
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Validate input BAM
    if not args.input_bam.is_file():
        logging.error(f"Input BAM file does not exist: {args.input_bam}")
        exit(1)

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

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
        output_bam = args.output_dir / f"{args.input_bam.stem}_downsampled_{int(fraction*100)}p.bam"
        subsample_bam(
            samtools=args.samtools,
            input_bam=subset_bam_path,
            output_bam=output_bam,
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

    # Step 4: Downsample to absolute coverages
    for coverage in args.coverages:
        fraction = calculate_required_fraction(desired_coverage=coverage, current_coverage=current_coverage)
        output_bam = args.output_dir / f"{args.input_bam.stem}_downsampled_{coverage}x.bam"
        subsample_bam(
            samtools=args.samtools,
            input_bam=subset_bam_path,
            output_bam=output_bam,
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


if __name__ == "__main__":
    main()
