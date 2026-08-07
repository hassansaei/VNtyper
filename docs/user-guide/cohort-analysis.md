# Cohort Analysis

The `cohort` subcommand aggregates results from multiple VNtyper 2 pipeline runs into a single summary report with visualizations.

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

    VNtyper 2 searches recursively for `pipeline_summary.json` in each directory.

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

This uses the default prefix `sample_` followed by the first 12 hex characters of a SHA-256 digest of the original name. Specify a custom prefix:

```bash
vntyper cohort -i results/* -o cohort_output/ \
    --pseudonymize-samples "patient_"
```

The digest and its width are configuration, not code. They live under `cohort.pseudonym` in `vntyper/config.json`:

```json
"cohort": {
  "pseudonym": {
    "algorithm": "sha256",
    "digest_characters": 12
  }
}
```

Override them with `--config-path`, remembering that a config file **replaces** the shipped one rather than merging into it. An algorithm your Python's `hashlib` does not offer is refused by name rather than silently substituted.

A `pseudonymization_table.tsv` mapping pseudonyms to original names is saved in the output directory.

## Sample Identity

A sample's name is only half of its identity. The other half is the **namespace** it came from: the name of the input that produced it -- `job_a` for an archive `job_a.zip`, `run1` for `-i /data/run1`. Uniqueness belongs to the pair, following HL7 FHIR's `Identifier` datatype, where a value is only ever "unique within the context of the system" that issued it. Two web jobs that each uploaded a file called `sample.bam` are therefore the same *value* in two different systems -- two patients, not a collision.

So when two or more discovered samples share a name, **those samples only** are reported as `namespace/name` -- `job_a/sample` and `job_b/sample` -- and the cohort runs to completion with both patients in it. A sample whose name is already unique keeps it untouched, so neither its row nor its pseudonym moves because some other pair of samples in the cohort happened to collide.

The run stops only where the inputs leave an ambiguity that cannot be reduced:

- **Two inputs with the same name.** Two archives both called `job.zip`, or two directory inputs both called `sample`, are one namespace, so qualification yields `job/sample` twice and nothing is left to tell the two apart. The error names both input paths; rename one of them, or give the two runs distinct recorded input files.
- **A digest collision.** Two distinct identities whose pseudonyms come out equal. Here the names are perfectly good and only the digest is too narrow, so the fix is to widen `cohort.pseudonym.digest_characters` rather than to rename anything. At 12 hex characters this is a tripwire rather than a practical risk (about 1.8e-9 for a 1,000-sample cohort).

Either way, refusing to run beats merging two patients' genotypes into one row.

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
| `cohort_stats.csv` | Per-sample statistics in CSV |
| `cohort_stats.tsv` | Per-sample statistics in TSV |
| `cohort_stats.json` | Per-sample statistics in JSON |

`cohort_stats_*` carries the same rows as the report's additional statistics table: runtime, version, assembly and alignment pipeline per sample, plus every coverage metric under a `cov_` prefix -- `cov_mean`, `cov_percent_uncovered` and `cov_coverage_qc` among them. It is the only machine-readable cohort output carrying a coverage figure; before VNtyper 2.0.8 that table reached the HTML report and nothing else.

## HTML Report Contents

The cohort summary report includes:

- **Donut charts** showing the distribution of Positive, Positive (Flagged), and Negative results for both Kestrel and adVNTR
- **Kestrel results table** with per-sample variant calls, confidence levels, and flags
- **adVNTR results table** (if adVNTR data is present)
- **Additional statistics** including runtime, coverage metrics, pipeline version, reference assembly, and alignment pipeline for each sample
