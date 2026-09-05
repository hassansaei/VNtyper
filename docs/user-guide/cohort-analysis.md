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

    Zip files extract to temporary staging directories during processing.

## Pseudonymization

Replace directory names with pseudonyms to mask identifiers in report tables:

```bash
vntyper cohort -i results/sample1/ results/sample2/ -o cohort_output/ \
    --pseudonymize-samples
```

This prefixes the first 12 hex characters of a SHA-256 digest of the original name with `sample_`. Specify a custom prefix if needed:

```bash
vntyper cohort -i results/* -o cohort_output/ \
    --pseudonymize-samples "patient_"
```

Digest parameters are configured under `cohort.pseudonym` in `vntyper/config.json`:

```json
"cohort": {
  "pseudonym": {
    "algorithm": "sha256",
    "digest_characters": 12
  }
}
```

Override them with `--config-path`. A custom config file replaces default configuration rather than merging. Unsupported hash algorithms raise configuration errors immediately.

A `pseudonymization_table.tsv` file mapping pseudonyms to original names is written to the output directory.

### What pseudonymization does not protect against

The digest is unsalted and unkeyed: only the sample identity enters the calculation. Anyone possessing the report without `pseudonymization_table.tsv` can recover guessable names (such as `sample_1` or sequential barcodes) by dictionary attack at one hash per candidate.

Pseudonymization provides visual obfuscation for presentations and collaborative sharing, not cryptographic de-identification. Treat pseudonymized reports as identifiable patient data. Restrict access to `pseudonymization_table.tsv`.

Pseudonyms remain stable for an invariant sample identity. If sample names collide, identities expand to `namespace/name` format (see [Sample Identity](#sample-identity)), altering the resulting digest.

## Sample Identity

Sample identity comprises both the sample name and its input namespace: the archive stem for `job_a.zip`, or the directory name for `-i /data/run1`. Identity pairs reflect HL7 FHIR `Identifier` conventions, where uniqueness is evaluated within the issuing system. Two web jobs each processing `sample.bam` represent distinct patients from separate namespaces.

When discovered samples share identical names, those specific records expand to `namespace/name` format (such as `job_a/sample` and `job_b/sample`). Unique sample names remain unmodified.

Execution halts when ambiguity cannot be resolved:

- **Identical input names:** Two archives named `job.zip` or two input directories named `sample` create identical `job/sample` namespaces. Rename one input source to proceed.
- **Digest collisions:** Distinct identities yielding identical truncated hashes. Resolve by increasing `cohort.pseudonym.digest_characters`. At 12 hex characters, collision probability remains negligible (approximately 1.8e-9 across 1,000 samples).

## Output Formats

HTML reports generate automatically. Additional formats are requested via `--summary-formats`:

```bash
vntyper cohort -i results/* -o cohort_output/ \
    --summary-formats csv,tsv,json
```

Outputs written to the destination directory:

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
| `cohort_call_frequency.csv` | Cohort call frequency summary in CSV |
| `cohort_call_frequency.tsv` | Cohort call frequency summary in TSV |
| `cohort_call_frequency.json` | Cohort call frequency summary in JSON |

`cohort_stats_*` exports sample runtime, version, assembly, alignment details, and coverage metrics prefixed with `cov_` (including `cov_mean`, `cov_percent_uncovered`, and `cov_coverage_qc`).

## HTML Report Contents

The cohort report includes:

- **Donut charts** displaying Positive, Positive (Flagged), and Negative proportions for Kestrel and adVNTR.
- **Kestrel results table** containing per-sample calls, confidence assignments, and flags.
- **adVNTR results table** (when adVNTR results exist).
- **Additional statistics** recording runtime, coverage QC, pipeline versions, reference assemblies, and aligners.
- **Call frequency table** aggregating variant calls across the cohort, sorted ascending by frequency, marking calls at or below threshold in the `Below_Cutoff` column.

## Cohort Call Frequency

The call frequency table aggregates variant calls across all cohort samples:

- **Grouping key:** For `Molecular_Identity_Status` values of `unique` or `legacy-selected-among-multiple` with defined molecular identities, calls group by `Molecular_Identity`. Otherwise, grouping falls back to caller format `<Motifs>:<POS>:<REF>:<ALT>`.
- **Grouping key kind:** The `Grouping_Key_Kind` column records `molecular-identity` or `caller-representation` to separate distinct identity tiers.
- **Cohort denominator:** Frequencies are calculated relative to total cohort size (`cohort_size`), including negative and uncalled samples.
- **Placeholder exclusion:** Negative placeholder records are excluded from call frequency groups.
- **Threshold marking:** Calls with `Frequency <= rare_allele_max_frequency` receive `Below_Cutoff = "yes"` (`"no"` otherwise). All calls remain visible in exports.
- **Configurable cutoff:** The cutoff defaults to `0.05` via `cohort.rare_allele_max_frequency` in `config.json` and can be overridden with `--rare-allele-max-frequency`.
