# Cohort Analysis

The `cohort` subcommand aggregates results from multiple VNtyper pipeline runs into a single summary report with visualizations.

## Basic Usage

```bash
vntyper cohort \
    -i results/sample1/ results/sample2/ results/sample3/ \
    -o cohort_output/
```

This scans each directory for `pipeline_summary.json`, extracts Kestrel and adVNTR results, and generates an HTML cohort report.

## Input Methods

=== "Directories"

    Pass directories directly with `-i` / `--input-dirs`:

    ```bash
    vntyper cohort -i results/sample1/ results/sample2/ -o cohort_output/
    ```

    VNtyper searches recursively for `pipeline_summary.json` in each directory.

=== "Input file"

    List directories or zip files in a text file (one per line):

    ```bash
    vntyper cohort --input-file sample_dirs.txt -o cohort_output/
    ```

    ```text title="sample_dirs.txt"
    results/sample1/
    results/sample2/
    /data/archived/sample3.zip
    ```

    Zip files are automatically extracted to temporary directories for processing.

## Pseudonymization

Protect sample identities by replacing directory names with pseudonyms:

```bash
vntyper cohort -i results/sample1/ results/sample2/ -o cohort_output/ \
    --pseudonymize-samples
```

This uses the default prefix `sample_` followed by the first 5 characters of an MD5 hash of the original name. Specify a custom prefix:

```bash
vntyper cohort -i results/* -o cohort_output/ \
    --pseudonymize-samples "patient_"
```

A `pseudonymization_table.tsv` mapping pseudonyms to original names is saved in the output directory.

## Output Formats

HTML is always generated. Request additional machine-readable formats:

```bash
vntyper cohort -i results/* -o cohort_output/ \
    --summary-formats csv,tsv,json
```

This produces:

| File | Content |
|------|---------|
| `cohort_summary.html` | Interactive HTML report (always generated) |
| `cohort_kestrel.csv` | Kestrel results in CSV |
| `cohort_kestrel.tsv` | Kestrel results in TSV |
| `cohort_kestrel.json` | Kestrel results in JSON |
| `cohort_advntr.csv` | adVNTR results in CSV (if adVNTR was run) |
| `cohort_advntr.tsv` | adVNTR results in TSV (if adVNTR was run) |
| `cohort_advntr.json` | adVNTR results in JSON (if adVNTR was run) |

## HTML Report Contents

The cohort summary report includes:

- **Donut charts** showing the distribution of Positive, Positive (Flagged), and Negative results for both Kestrel and adVNTR
- **Kestrel results table** with per-sample variant calls, confidence levels, and flags
- **adVNTR results table** (if adVNTR data is present)
- **Additional statistics** including runtime, coverage metrics, pipeline version, reference assembly, and alignment pipeline for each sample
