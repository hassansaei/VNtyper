# VNtyper Benchmarking and Plotting

This repository provides two main scripts to help you:

1. **Downsample and benchmark** BAM files for MUC1 VNTR analysis, optionally running `vntyper`.  
2. **Plot** the resulting `vntyper_summary.csv` as scatter plots (one subplot per metric).

## Table of Contents
1. [Requirements](#requirements)
2. [Scripts Overview](#scripts-overview)
   - [downsample_bam.py](#downsample_bampy)
   - [plot_vntyper_summary.py](#plot_vntyper_summarypy)
3. [Usage](#usage)
   - [Example: downsample and run vntyper](#example-downsample-and-run-vntyper)
   - [Example: generate scatter plots](#example-generate-scatter-plots)
4. [Contact](#contact)

---

## Requirements

- **Python 3.7+**  
- **Samtools** (installed and in PATH, if you're downsampling)  
- **vntyper** (installed and in PATH if using the `--run-vntyper` functionality)  
- **Matplotlib** and **pandas** (for the plotting script). Install with:
  ```bash
  pip install matplotlib pandas
  ```

---

## Scripts Overview

### downsample_bam.py

> **Location**: `tests/benchmark/downsample_bam.py`

**Purpose**  
Downsamples a BAM file to a specific **region** (MUC1) and optionally **downsamples** further to desired fractions or absolute coverage levels. If `--run-vntyper` is set, each downsampled BAM is processed with `vntyper`, and a summary CSV is generated.

**Usage**  
```bash
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
- `--vntyper-path`: Path to the `vntyper` executable (defaults to `vntyper` in PATH).  
- `--summary-output`: Path to the final CSV containing the summary of results.  

**What It Produces**  
- Downsampled BAM files (e.g., `_downsampled_50x.bam`) plus `*.bam.md5` checksums.  
- If `--run-vntyper` is set: a summary CSV (default: `vntyper_summary.csv`) containing columns like:
  ```
  file_analyzed,method,value,confidence,
  Estimated_Depth_AlternateVariant,Estimated_Depth_Variant_ActiveRegion,Depth_Score
  ```
- For each downsampled BAM, if `vntyper` is run, additional subdirectories with `kestrel_result.tsv` may be created.

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
python plot_vntyper_summary.py \
    --input-csv path/to/vntyper_summary.csv \
    --output-png path/to/vntyper_summary.png
```

**Key Arguments**  
- `--input-csv`: The summary CSV file created by `downsample_bam.py` when `--run-vntyper` is used.  
- `--output-png`: File path for the generated PNG plot (default `vntyper_summary_plot.png`).  

**What It Produces**  
- A **PNG** with 3 subplots (side by side):
  - **X-axis**: each metric (`Estimated_Depth_AlternateVariant`, etc.)
  - **Y-axis**: the `value` column in the CSV (i.e., fraction or coverage).  
  - Points color-coded by `confidence`.  
  - Axes start at 0.  

---

## Usage

### Example: downsample and run vntyper

```bash
python tests/benchmark/downsample_bam.py \
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
3. Downsample to coverage levels 25, 50, 100, etc.  
4. Run `vntyper` on each downsampled BAM.  
5. Write the final summary CSV to `out/benchmark/example_c495/vntyper_summary.csv`.

Afterwards, check `out/benchmark/example_c495` for:
- Downsampled BAMs  
- `kestrel_result.tsv` subdirectories (if any)  
- The CSV summary

### Example: generate scatter plots

```bash
python tests/benchmark/plot_vntyper_summary.py \
    --input-csv out/benchmark/example_c495/vntyper_summary.csv \
    --output-png out/benchmark/example_c495/vntyper_summary.png
```

- Reads the summary CSV
- Produces a `vntyper_summary.png` with 3 subplots for the 3 metrics (X-axis) vs. the `value` column (Y-axis).

Check `out/benchmark/example_c495/vntyper_summary.png` to see your scatter plots.
