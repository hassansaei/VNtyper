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

The `--pseudonymize-samples` flag supports two modes. In both, the pseudonym is the prefix followed by the first 12 hex characters of the SHA-256 digest of the original sample name, so it is stable across runs and machines:

- **Default basename:** `--pseudonymize-samples` (no value) uses `sample_` as the prefix, so `sample1` and `sample2` are reported as `sample_e85130791f31` and `sample_5a9392784e07`.
- **Custom basename:** `--pseudonymize-samples patient_` uses the provided prefix, giving `patient_e85130791f31` and `patient_5a9392784e07` for the same two samples.

The digest algorithm and the number of characters kept are read from `cohort.pseudonym` in `config.json` (`{"algorithm": "sha256", "digest_characters": 12}`), so `--config-path` can change either -- bearing in mind that a supplied config replaces the shipped one rather than merging into it. An algorithm that is not available in the running interpreter's `hashlib` is refused by name.

The pseudonym map is one-to-one or the run stops: two input directories that share a basename, or two sample names that collide on the digest, raise an error naming both instead of silently reporting the two patients as one sample.

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
