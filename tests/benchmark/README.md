# VNtyper Benchmarking and Plotting

This repository provides three main scripts to help you:

1. **Downsample and benchmark** BAM files for MUC1 VNTR analysis, optionally running `vntyper`.  
2. **Benchmark** vntyper results against known mutation status using a dedicated benchmarking script.  
3. **Plot** the resulting `vntyper_summary.csv` as scatter plots (one subplot per metric).

## Table of Contents
1. [Requirements](#requirements)
2. [Scripts Overview](#scripts-overview)
   - [benchamrk_downsample.py](#benchamrk_downsamplepy)
   - [benchmark_vntyper.py](#benchmark_vntyperpy)
   - [plot_vntyper_summary.py](#plot_vntyper_summarypy)
3. [Usage](#usage)
   - [Example: downsample and run vntyper](#example-downsample-and-run-vntyper)
   - [Example: benchmark vntyper results](#example-benchmark-vntyper-results)
   - [Example: generate scatter plots](#example-generate-scatter-plots)
4. [Contact](#contact)

---

## Requirements

- **Python 3.7+**  
- **Samtools** (installed and in PATH, if you're downsampling)  
- **vntyper** (installed and in PATH if using the `--run-vntyper` or benchmarking functionality)  
- **Matplotlib** and **pandas** (for the plotting script). Install with:
  ```bash
  pip install matplotlib pandas
  ```

---

## Scripts Overview

### benchamrk_downsample.py

> **Location**: `tests/benchmark/benchamrk_downsample.py`

**Purpose**  
Downsamples a BAM file to a specific **region** (MUC1) and optionally **further downsamples** to desired fractions or absolute coverage levels. If `--run-vntyper` is set, each downsampled BAM is processed with `vntyper` (with an optional adVNTR mode via `--run-advntr`), and a summary CSV is generated.

**Usage**  
```bash
python tests/benchmark/benchamrk_downsample.py \
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
    [--reference-assembly hg19|hg38] \
    [--keep-intermediates] \
    [--archive-results] \
    [--fast-mode] \
    [--summary-output summary.csv]
```

**Key Arguments**  
- `--input-bam`: Path to your input BAM file.  
- `--output-dir`: Directory where downsampled BAMs (and optional `vntyper` outputs) are stored.  
- `--fractions`: List of fractions to downsample (range 0–1).  
- `--coverages`: List of absolute coverage targets to downsample to (integer).  
- `--muc1-region`: MUC1 region in `chr:start-end` format.  
- `--vntr-region`: Sub-region for coverage calculation.  
- `--run-vntyper`: Run `vntyper` on downsampled BAMs if set.  
- `--run-advntr`: Include adVNTR mode (appends `--extra-modules advntr` to the vntyper command and parses its results).  
- `--vntyper-path`: Path to the `vntyper` executable (defaults to `vntyper` in PATH).  
- `--summary-output`: Path to the final CSV containing the summary of results.  

**What It Produces**  
- Downsampled BAM files (e.g., `_downsampled_50x.bam`) plus `*.bam.md5` checksums.  
- If `--run-vntyper` is set: a summary CSV (default: `vntyper_summary.csv`) containing columns like:
  ```
  file_analyzed,method,value,confidence,
  Estimated_Depth_AlternateVariant,Estimated_Depth_Variant_ActiveRegion,Depth_Score,analysis_time_minutes,advntr_result
  ```
- For each downsampled BAM, if `vntyper` is run, additional subdirectories with `kestrel_result.tsv` (and possibly adVNTR results) may be created.

---

### benchmark_vntyper.py

> **Location**: `tests/benchmark/benchmark_vntyper.py`

**Purpose**  
Benchmarks `vntyper` results on (simulated) BAM files by comparing the predicted mutation status (derived from the categorical `Confidence` field) with the known status provided in a CSV/TSV file. The script computes a confusion matrix and test statistics, including sensitivity (recall), specificity, precision (PPV), negative predictive value (NPV), and accuracy.

**New Feature**  
- **Result Caching:**  
  The script checks if a `vntyper` output directory already exists (with a `kestrel_result.tsv` file) and skips recomputation unless the `--recompute` flag is set.

**Usage**  
```bash
python tests/benchmark/benchmark_vntyper.py \
    --sample-info path/to/sample_info.csv \
    --delimiter , \
    --bam-col bam \
    --status-col status \
    --vntyper-path path/to/vntyper \
    --reference-assembly hg38 \
    --threads 4 \
    --output-dir path/to/vntyper_outputs \
    --summary-output benchmark_summary.csv \
    --stats-output benchmark_stats.csv \
    [--recompute] \
    [--keep-intermediates] [--archive-results] [--fast-mode] \
    [--vntyper-options "additional options"]
```

**Key Arguments**  
- `--sample-info`: Path to the CSV/TSV file containing the BAM file paths and known mutation status.  
- `--delimiter`: Field delimiter for the sample info file (e.g., `,` for CSV or `\t` for TSV).  
- `--bam-col` and `--status-col`: Column names for the BAM file paths and mutation status (default: `bam` and `status`).  
- `--vntyper-path`: Path to the `vntyper` executable (defaults to `vntyper` in PATH).  
- `--reference-assembly`: Reference assembly (e.g., hg19 or hg38).  
- `--threads`: Number of threads for processing.  
- `--output-dir`: Directory to store `vntyper` outputs.  
- `--summary-output`: Path to the per-sample summary CSV file.  
- `--stats-output`: Path to the overall test statistics CSV file.  
- `--recompute`: If set, forces recomputation of `vntyper` results even if they already exist.

**What It Produces**  
- **Per-Sample Summary CSV:** Contains details per sample (sample ID, BAM path, expected and predicted statuses, output directory).  
- **Overall Test Statistics CSV:** Contains the confusion matrix (TP, FN, FP, TN, Total) and computed test statistics (Sensitivity, Specificity, Precision, NPV, Accuracy).

---

### plot_vntyper_summary.py

> **Location**: `tests/benchmark/plot_vntyper_summary.py`

**Purpose**  
Generates **scatter plots** (no lines) from a `vntyper_summary.csv`. Plots 3 subplots for:
1. **Estimated_Depth_AlternateVariant** (X-axis) vs. `value` (Y-axis),
2. **Estimated_Depth_Variant_ActiveRegion** (X-axis) vs. `value` (Y-axis),
3. **Depth_Score** (X-axis) vs. `value` (Y-axis).

The points are **colored** by `confidence`:
- **High_Precision** → red
- **Low_Precision** → orange
- **Negative** → black
- **Else** → gray

**Usage**  
```bash
python tests/benchmark/plot_vntyper_summary.py \
    --input-csv path/to/vntyper_summary.csv \
    --output-png path/to/vntyper_summary.png
```

**Key Arguments**  
- `--input-csv`: The summary CSV file created by `benchamrk_downsample.py` when `--run-vntyper` is used.  
- `--output-png`: File path for the generated PNG plot (default: `vntyper_summary_plot.png`).

**What It Produces**  
- A **PNG** with 3 subplots (side by side):
  - **X-axis**: each metric (`Estimated_Depth_AlternateVariant`, etc.)
  - **Y-axis**: the `value` column in the CSV (i.e., fraction or coverage).  
  - Points color-coded by `confidence`.  
  - Axes starting at 0.

---

## Usage

### Example: downsample and run vntyper

```bash
python tests/benchmark/benchamrk_downsample.py \
    --input-bam tests/data/example_c495.bam \
    --output-dir out/benchmark/example_c495 \
    --coverages 25 50 100 125 150 175 200 225 250 300 500 \
    --seed 42 \
    --threads 16 \
    --muc1-region 'chr1:155158000-155163000' \
    --vntr-region 'chr1:155160500-155162000' \
    --run-vntyper \
    --vntyper-path vntyper \
    --reference-assembly hg19 \
    --keep-intermediates \
    --archive-results \
    --fast-mode \
    --summary-output out/benchmark/example_c495/vntyper_summary.csv
```

This command will:

1. Create `out/benchmark/example_c495/` (if not existing).  
2. Subset the BAM to the MUC1 region.  
3. Downsample to various coverage levels.  
4. Run `vntyper` (and optionally adVNTR if `--run-advntr` is added) on each downsampled BAM.  
5. Write the final summary CSV to `out/benchmark/example_c495/vntyper_summary.csv`.

Afterwards, check `out/benchmark/example_c495` for:
- Downsampled BAMs  
- `kestrel_result.tsv` subdirectories (if any)  
- The CSV summary

---

### Example: benchmark vntyper results

```bash
python tests/benchmark/benchmark_vntyper.py \
    --sample-info tests/data/sample_info.csv \
    --delimiter , \
    --bam-col bam \
    --status-col status \
    --vntyper-path vntyper \
    --reference-assembly hg38 \
    --threads 8 \
    --output-dir out/benchmark/benchmark_results \
    --summary-output out/benchmark/benchmark_summary.csv \
    --stats-output out/benchmark/benchmark_stats.csv
```

This command will:

1. Read the sample info CSV, which contains the BAM file paths and known mutation statuses.
2. For each sample, check if a corresponding `vntyper` output exists.  
   - If results already exist, they will be reused unless the `--recompute` flag is specified.
3. Parse the vntyper results and compare the predicted status with the expected status.
4. Write a per-sample summary CSV (`benchmark_summary.csv`) and an overall statistics CSV (`benchmark_stats.csv`) containing the confusion matrix and test statistics (Sensitivity, Specificity, Precision, NPV, Accuracy).

---

### Example: generate scatter plots

```bash
python tests/benchmark/plot_vntyper_summary.py \
    --input-csv out/benchmark/example_c495/vntyper_summary.csv \
    --output-png out/benchmark/example_c495/vntyper_summary.png
```

- Reads the summary CSV.
- Produces a `vntyper_summary.png` with 3 subplots for the metrics versus the `value` column.

Check `out/benchmark/example_c495/vntyper_summary.png` to see your scatter plots.
