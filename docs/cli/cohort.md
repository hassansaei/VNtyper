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
| `--pseudonymize-samples` | string (optional value) | — | Replace sample names with a prefix plus a truncated SHA-256 digest. Obfuscation for readability, **not** a privacy control -- the digest is unsalted and unkeyed, so guessable names are recoverable from the report alone; see [Cohort analysis](../user-guide/cohort-analysis.md#what-pseudonymization-does-not-protect-against). Defaults to the basename `sample_`; optionally provide a custom basename |

One of `-i/--input-dirs` or `--input-file` is required.

## Pseudonymization

The `--pseudonymize-samples` flag supports two modes. In both, the pseudonym is the prefix followed by the first 12 hex characters of the SHA-256 digest of the original sample name, so it is stable across runs and machines:

- **Default basename:** `--pseudonymize-samples` (no value) uses `sample_` as the prefix, so `sample1` and `sample2` are reported as `sample_e85130791f31` and `sample_5a9392784e07`.
- **Custom basename:** `--pseudonymize-samples patient_` uses the provided prefix, giving `patient_e85130791f31` and `patient_5a9392784e07` for the same two samples.

The digest algorithm and the number of characters kept are read from `cohort.pseudonym` in `config.json` (`{"algorithm": "sha256", "digest_characters": 12}`), so `--config-path` can change either -- bearing in mind that a supplied config replaces the shipped one rather than merging into it. An algorithm that is not available in the running interpreter's `hashlib` is refused by name.

## Sample identity

A sample's identity is a namespace and a value. The value is the sample's own name -- its directory name, or the stem of the input file that a zipped run recorded. The namespace is the name of the input it was reached through: the archive stem for `job_a.zip`, the directory name for `-i /data/run1`. The input's *name* is used rather than its path, so `-i /data/run1` and `-i ./run1` are the same namespace and a pseudonym survives the cohort being relocated.

Uniqueness is a property of the pair, not of the value alone. When two or more discovered samples share a value -- two web jobs that both uploaded `sample.bam`, for instance -- those samples are reported as `namespace/value` (`job_a/sample` and `job_b/sample`) and the cohort completes. Samples whose value is already unique are left exactly as they are, so a collision elsewhere in the cohort never moves their identity or their pseudonym.

The run stops only on ambiguity qualification cannot reduce: two inputs sharing a name (two archives both called `job.zip`, or two directory inputs both called `sample`), which qualify to one identity twice; and a digest collision between two identities that are already distinct, whose fix is a wider `cohort.pseudonym.digest_characters`. Both raise an error naming the two samples and the inputs they came from, instead of silently reporting two patients as one.

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
