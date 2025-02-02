#!/usr/bin/env python3
"""
Benchmarking vntyper on simulated BAM files

This script reads a CSV/TSV file containing (at least) two columns:
  - one with the path to each BAM file, and 
  - one with the known mutation status ("positive" or "negative").

For each sample, the script will:
  1. Check that the BAM file exists.
  2. If the vntyper results do not already exist (or if --recompute is set),
     run vntyper on the BAM file.
  3. Parse the vntyper output (specifically the "Confidence" field) to decide
     whether the sample is positive or negative.
     The logic here is: if the reported confidence is not "negative" (e.g.
     it is one of "Low_Precision", "High_Precision", or "High_Precision*"),
     then the sample is predicted to be positive.
  4. Compare the prediction with the expected status.
  5. Compute test statistics (confusion matrix, sensitivity, specificity,
     precision, NPV, and accuracy) and write a detailed per-sample summary as well
     as these aggregated statistics (including standard errors and 95% confidence intervals)
     to output CSV files.

Usage:
    python benchmark_vntyper.py \
        --sample-info samples.csv \
        --delimiter , \
        --bam-col bam \
        --status-col status \
        --vntyper-path /path/to/vntyper \
        --reference-assembly hg19 \
        --threads 4 \
        --output-dir ./vntyper_outputs \
        --summary-output benchmark_summary.csv \
        --stats-output benchmark_stats.csv \
        [--recompute] \
        [--keep-intermediates] [--archive-results] [--fast-mode] [--vntyper-options "additional options"]

Note:
  - The input CSV/TSV file must have a header row containing the column names.
  - By default, the script assumes the BAM file column is named "bam" and the
    mutation status column is named "status". Use --bam-col and --status-col to override.
  - This script uses the categorical nature of the Confidence field output by vntyper.
"""

import argparse
import csv
import hashlib
import logging
import math
import subprocess
from pathlib import Path
from typing import List, Optional, Dict

import pandas as pd

# -----------------------
# Logging configuration
# -----------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# -----------------------
# Utility functions
# -----------------------
def run_command(command: List[str], description: str, critical: bool = True) -> None:
    """
    Run a shell command and handle errors.
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
    Compute the MD5 checksum of a file.
    """
    hash_md5 = hashlib.md5()
    logging.debug(f"Computing MD5 for {file_path}")
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

# -----------------------
# vntyper-specific functions
# -----------------------
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
    Run vntyper on the given BAM file.

    Args:
        vntyper_path (str): Path to the vntyper executable.
        bam_file (Path): Input BAM file.
        reference_assembly (str): Reference assembly (e.g., hg19).
        threads (int): Number of threads to use.
        output_dir (Path): Directory to store vntyper outputs.
        keep_intermediates (bool): Whether to keep intermediate files.
        archive_results (bool): Whether to archive results.
        fast_mode (bool): Whether to run in fast mode.
        additional_options (str, optional): Extra options for vntyper.

    Returns:
        Path: Directory where vntyper outputs are stored.
    """
    sample_output = output_dir / f"{bam_file.stem}_vntyper_output"
    sample_output.mkdir(parents=True, exist_ok=True)

    command = [
        vntyper_path,
        "-l", "DEBUG",
        "pipeline",
        "--bam", str(bam_file),
        "--threads", str(threads),
        "--reference-assembly", reference_assembly,
        "-o", str(sample_output),
    ]
    if keep_intermediates:
        command.append("--keep-intermediates")
    if archive_results:
        command.append("--archive-results")
    if fast_mode:
        command.append("--fast-mode")
    if additional_options:
        command.extend(additional_options.split())

    logging.info(f"Executing vntyper for {bam_file.name}")
    run_command(command, description=f"vntyper run on {bam_file.name}")

    return sample_output

def summarize_vntyper_results(vntyper_output_dir: Path) -> Optional[Dict]:
    """
    Parse vntyper results to extract key metrics from the 'kestrel_result.tsv' file.
    Expects the file to contain a column "Confidence" whose value is either
    "negative" or one of:
       "Low_Precision", "High_Precision", "High_Precision*".
    
    Returns:
        A dictionary with the key "Confidence" (a list of strings), or None if parsing fails.
    """
    result_files = list(vntyper_output_dir.rglob("kestrel_result.tsv"))
    if not result_files:
        logging.warning(f"No kestrel_result.tsv found in {vntyper_output_dir}")
        return None

    result_file = result_files[0]
    logging.info(f"Parsing Kestrel results from {result_file}")

    try:
        df = pd.read_csv(result_file, sep="\t", comment="#")
        if "Confidence" not in df.columns:
            logging.error(f"Missing expected column 'Confidence' in {result_file}")
            return None

        summary = {
            "Confidence": df["Confidence"].dropna().astype(str).tolist()
        }
        return summary
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse {result_file}: {e}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error while parsing {result_file}: {e}")
        return None

def determine_call(summary: Optional[Dict]) -> str:
    """
    Determine the predicted mutation status based on the vntyper summary.
    
    If no summary is available or if the only reported confidence is "negative",
    the sample is called "negative". Otherwise (if any reported confidence is not
    "negative"), the sample is called "positive".
    """
    if summary is None:
        return "negative"
    confidences = summary.get("Confidence", [])
    if not confidences:
        return "negative"
    for conf in confidences:
        if conf.strip().lower() != "negative":
            return "positive"
    return "negative"

# -----------------------
# Confidence interval helper
# -----------------------
def calc_ci(p: float, n: int, z: float = 1.96):
    """
    Calculate the standard error and 95% confidence interval for a proportion.
    
    Args:
        p (float): Estimated proportion.
        n (int): Denominator (sample size).
        z (float): z-score for desired confidence level (default 1.96 for ~95%).
    
    Returns:
        Tuple: (standard error, lower bound, upper bound) with bounds clamped to [0, 1].
    """
    if n == 0:
        return 0, 0, 0
    se = math.sqrt(p * (1 - p) / n)
    lower = max(0, p - z * se)
    upper = min(1, p + z * se)
    return se, lower, upper

# -----------------------
# Argument parsing
# -----------------------
def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark vntyper on simulated BAM files using known mutation status."
    )
    parser.add_argument(
        "--sample-info",
        type=Path,
        required=True,
        help="Path to CSV/TSV file with sample info. Must contain columns for BAM path and mutation status.",
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        default=",",
        help="Delimiter for the sample info file (default: ','). Use '\\t' for TSV.",
    )
    parser.add_argument(
        "--bam-col",
        type=str,
        default="bam",
        help="Column name for BAM file paths (default: 'bam').",
    )
    parser.add_argument(
        "--status-col",
        type=str,
        default="status",
        help="Column name for known mutation status (default: 'status').",
    )
    parser.add_argument(
        "--vntyper-path",
        type=str,
        default="vntyper",
        help="Path to the vntyper executable (default: 'vntyper' from PATH).",
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
        help="Reference assembly for vntyper (default: hg19).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use (default: 4).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("./vntyper_outputs"),
        help="Directory to store vntyper outputs (default: ./vntyper_outputs).",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path("benchmark_summary.csv"),
        help="Path to write detailed per-sample results (default: benchmark_summary.csv).",
    )
    parser.add_argument(
        "--stats-output",
        type=Path,
        default=Path("benchmark_stats.csv"),
        help="Path to write the overall confusion matrix and test statistics (default: benchmark_stats.csv).",
    )
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="Force recomputation of vntyper results even if they already exist.",
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
    return parser.parse_args()

# -----------------------
# Main processing
# -----------------------
def main():
    args = parse_arguments()

    # Ensure output directory exists.
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # List to hold per-sample results.
    results = []

    # Counters for confusion matrix.
    tp = fp = tn = fn = 0

    # Read sample information from the CSV/TSV file.
    with args.sample_info.open("r") as infile:
        reader = csv.DictReader(infile, delimiter=args.delimiter)
        for row in reader:
            bam_path = Path(row[args.bam_col].strip())
            truth_status = row[args.status_col].strip().lower()
            sample_id = bam_path.stem

            logging.info(f"Processing sample {sample_id} with expected status '{truth_status}'")

            if not bam_path.is_file():
                logging.error(f"BAM file not found: {bam_path}. Skipping sample.")
                continue

            # Determine the expected vntyper output directory.
            vntyper_out_dir = args.output_dir / f"{bam_path.stem}_vntyper_output"
            # Check if results already exist (i.e. a kestrel_result.tsv file is present)
            if not args.recompute and vntyper_out_dir.exists() and any(vntyper_out_dir.rglob("kestrel_result.tsv")):
                logging.info(f"vntyper output already exists for {bam_path.name}, skipping recomputation")
            else:
                try:
                    vntyper_out_dir = run_vntyper(
                        vntyper_path=args.vntyper_path,
                        bam_file=bam_path,
                        reference_assembly=args.reference_assembly,
                        threads=args.threads,
                        output_dir=args.output_dir,
                        keep_intermediates=args.keep_intermediates,
                        archive_results=args.archive_results,
                        fast_mode=args.fast_mode,
                        additional_options=args.vntyper_options,
                    )
                except Exception as e:
                    logging.error(f"vntyper run failed for {bam_path}: {e}")
                    predicted_status = "error"
                    summary = None
                else:
                    # Optionally compute MD5 of the BAM (for record keeping).
                    md5sum = compute_md5(bam_path)
                    md5_file = bam_path.with_suffix(bam_path.suffix + ".md5")
                    with md5_file.open("w") as mf:
                        mf.write(f"{md5sum}  {bam_path.name}\n")
            # In any case, try to parse existing results.
            summary = summarize_vntyper_results(vntyper_out_dir)
            predicted_status = determine_call(summary)

            # Update confusion matrix counts (skip if error).
            if predicted_status != "error":
                if truth_status == "positive":
                    if predicted_status == "positive":
                        tp += 1
                    else:
                        fn += 1
                elif truth_status == "negative":
                    if predicted_status == "negative":
                        tn += 1
                    else:
                        fp += 1

            # Save per-sample result.
            results.append({
                "sample_id": sample_id,
                "bam_path": str(bam_path),
                "expected_status": truth_status,
                "predicted_status": predicted_status,
                "vntyper_output": str(vntyper_out_dir),
            })

    total = tp + tn + fp + fn
    logging.info("=== Confusion Matrix ===")
    logging.info(f"True Positives:  {tp}")
    logging.info(f"False Negatives: {fn}")
    logging.info(f"False Positives: {fp}")
    logging.info(f"True Negatives:  {tn}")
    logging.info(f"Total evaluated: {total}")

    # Calculate test statistics.
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0  # Recall
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    accuracy = (tp + tn) / total if total > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0   # Positive Predictive Value (PPV)
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0           # Negative Predictive Value (NPV)

    # Compute standard errors and 95% confidence intervals.
    sens_n = tp + fn
    spec_n = tn + fp
    prec_n = tp + fp
    npv_n = tn + fn
    acc_n = total

    sens_se, sens_lower, sens_upper = calc_ci(sensitivity, sens_n)
    spec_se, spec_lower, spec_upper = calc_ci(specificity, spec_n)
    prec_se, prec_lower, prec_upper = calc_ci(precision, prec_n)
    npv_se, npv_lower, npv_upper = calc_ci(npv, npv_n)
    acc_se, acc_lower, acc_upper = calc_ci(accuracy, acc_n)

    logging.info("=== Test Statistics ===")
    logging.info(f"Sensitivity (Recall): {sensitivity:.3f} (SE: {sens_se:.3f}, CI: [{sens_lower:.3f}, {sens_upper:.3f}])")
    logging.info(f"Specificity: {specificity:.3f} (SE: {spec_se:.3f}, CI: [{spec_lower:.3f}, {spec_upper:.3f}])")
    logging.info(f"Precision (PPV): {precision:.3f} (SE: {prec_se:.3f}, CI: [{prec_lower:.3f}, {prec_upper:.3f}])")
    logging.info(f"NPV: {npv:.3f} (SE: {npv_se:.3f}, CI: [{npv_lower:.3f}, {npv_upper:.3f}])")
    logging.info(f"Accuracy: {accuracy:.3f} (SE: {acc_se:.3f}, CI: [{acc_lower:.3f}, {acc_upper:.3f}])")

    # Write detailed per-sample results to summary CSV.
    with args.summary_output.open("w", newline="") as csvfile:
        fieldnames = ["sample_id", "bam_path", "expected_status", "predicted_status", "vntyper_output"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for rec in results:
            writer.writerow(rec)
    logging.info(f"Per-sample benchmark summary written to {args.summary_output}")

    # Write the overall test statistics and confusion matrix with SE and CIs to the stats output file.
    with args.stats_output.open("w", newline="") as csvfile:
        fieldnames = [
            "TP", "FN", "FP", "TN", "Total",
            "Sensitivity", "Sensitivity_SE", "Sensitivity_CI_lower", "Sensitivity_CI_upper",
            "Specificity", "Specificity_SE", "Specificity_CI_lower", "Specificity_CI_upper",
            "Precision", "Precision_SE", "Precision_CI_lower", "Precision_CI_upper",
            "NPV", "NPV_SE", "NPV_CI_lower", "NPV_CI_upper",
            "Accuracy", "Accuracy_SE", "Accuracy_CI_lower", "Accuracy_CI_upper"
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({
            "TP": tp,
            "FN": fn,
            "FP": fp,
            "TN": tn,
            "Total": total,
            "Sensitivity": f"{sensitivity:.3f}",
            "Sensitivity_SE": f"{sens_se:.3f}",
            "Sensitivity_CI_lower": f"{sens_lower:.3f}",
            "Sensitivity_CI_upper": f"{sens_upper:.3f}",
            "Specificity": f"{specificity:.3f}",
            "Specificity_SE": f"{spec_se:.3f}",
            "Specificity_CI_lower": f"{spec_lower:.3f}",
            "Specificity_CI_upper": f"{spec_upper:.3f}",
            "Precision": f"{precision:.3f}",
            "Precision_SE": f"{prec_se:.3f}",
            "Precision_CI_lower": f"{prec_lower:.3f}",
            "Precision_CI_upper": f"{prec_upper:.3f}",
            "NPV": f"{npv:.3f}",
            "NPV_SE": f"{npv_se:.3f}",
            "NPV_CI_lower": f"{npv_lower:.3f}",
            "NPV_CI_upper": f"{npv_upper:.3f}",
            "Accuracy": f"{accuracy:.3f}",
            "Accuracy_SE": f"{acc_se:.3f}",
            "Accuracy_CI_lower": f"{acc_lower:.3f}",
            "Accuracy_CI_upper": f"{acc_upper:.3f}",
        })
    logging.info(f"Overall test statistics written to {args.stats_output}")

if __name__ == "__main__":
    main()
