# cohort

Aggregate outputs from multiple pipeline runs into a single cohort summary.

## Synopsis

```
vntyper [global-options] cohort
    (-i <dir> [<dir> ...] | --input-file <file>)
    -o <dir>
    [--summary-file <name>]
    [--summary-formats <formats>]
    [--pseudonymize-samples [<basename>]]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-i, --input-dirs` | list | — | List of directories containing output files to aggregate. Mutually exclusive with `--input-file` |
| `--input-file` | path | — | Path to a newline-separated text file listing directories or zip files to aggregate. Mutually exclusive with `-i` |
| `-o, --output-dir` | path | (required) | Output directory for the aggregated summary |
| `--summary-file` | string | `cohort_summary.html` | Name of the cohort summary report file |
| `--summary-formats` | string | `""` | Comma-separated list of additional summary output formats to generate (supported: `csv`, `tsv`, `json`). HTML is always generated |
| `--pseudonymize-samples` | string (optional value) | — | Pseudonymize sample names to protect sensitive information. When used without a value, defaults to the basename `sample_`. Optionally provide a custom basename |

One of `-i/--input-dirs` or `--input-file` is required.

## Pseudonymization

The `--pseudonymize-samples` flag supports two modes:

- **Default basename:** `--pseudonymize-samples` (no value) uses `sample_` as the prefix, producing names like `sample_a3f8c`, `sample_1b2d0` (prefix + first 5 characters of the MD5 hash of the original sample name).
- **Custom basename:** `--pseudonymize-samples patient_` uses the provided prefix, producing names like `patient_a3f8c`, `patient_1b2d0`.

## Examples

Aggregate results from multiple directories:

```bash
vntyper cohort -i results/sample1/ results/sample2/ results/sample3/ -o cohort_output/
```

Aggregate from a file listing directories:

```bash
vntyper cohort --input-file sample_dirs.txt -o cohort_output/
```

Generate additional output formats with pseudonymized sample names:

```bash
vntyper cohort -i results/*/ -o cohort_output/ \
    --summary-formats csv,tsv,json --pseudonymize-samples
```

Use a custom pseudonymization prefix:

```bash
vntyper cohort -i results/*/ -o cohort_output/ --pseudonymize-samples patient_
```
