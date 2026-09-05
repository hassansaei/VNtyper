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
    [--rare-allele-max-frequency <float>]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-i, --input-dirs` | list | None | Space-separated list of run directories to aggregate (mutually exclusive with `--input-file`) |
| `--input-file` | path | None | Text file listing paths to run directories or zip archives, one per line (mutually exclusive with `-i`) |
| `-o, --output-dir` | path | (required) | Target directory for aggregated cohort files |
| `--summary-file` | string | `cohort_summary.html` | File name for HTML cohort report |
| `--summary-formats` | string | `""` | Comma-separated summary formats (`csv`, `tsv`, `json`). HTML generates automatically |
| `--pseudonymize-samples` | string (optional value) | None | Replace sample identifiers with a prefix and truncated SHA-256 digest. Obfuscation for readability, not a privacy control (guessable names are vulnerable to dictionary matching; see [Cohort analysis](../user-guide/cohort-analysis.md#what-pseudonymization-does-not-protect-against)). Defaults to prefix `sample_` |
| `--rare-allele-max-frequency` | float | `0.05` | Frequency cutoff (0 to 1) for marking calls as `Below_Cutoff` in call frequency tables. Calls above threshold remain visible |

Supply either `-i/--input-dirs` or `--input-file`.

## Decision-profile grouping

For schema-3 inputs, cohort processing validates each run-local decision-profile snapshot, exporting `Decision_Profile_ID`, `Decision_Profile_Revision`, and `Decision_Profile_SHA256` alongside variant rows. Samples group by matching profile hash. If multiple profile hashes appear, pooled performance plots split into distinct profile cohorts to prevent discordant policies from mixing.

Legacy runs without profile records are designated `decision profile not recorded by legacy run`.

## Pseudonymization

The `--pseudonymize-samples` option operates in two modes, generating prefixes followed by 12 hex characters of the SHA-256 digest of the sample identifier:

- **Default prefix:** `--pseudonymize-samples` applies `sample_`, producing identifiers such as `sample_e85130791f31`.
- **Custom prefix:** `--pseudonymize-samples patient_` applies the specified prefix, producing `patient_e85130791f31`.

Digest algorithm and character length configure under `cohort.pseudonym` in `config.json`.

## Sample identity

Sample identity combines the sample name and its input namespace. The value represents the sample name; the namespace reflects the input source (the archive stem for `job_a.zip`, or directory name for `-i /data/run1`).

When sample names collide across runs, conflicting rows format as `namespace/value` (e.g. `job_a/sample` and `job_b/sample`). Unique names remain untouched.

Processing halts when identical namespaces collide (such as two input archives named `job.zip` containing `sample.bam`) or upon hash collisions.

## Examples

Aggregate multiple run directories:

```bash
vntyper cohort -i results/sample1/ results/sample2/ results/sample3/ -o cohort_output/
```

Aggregate directories listed in a manifest file:

```bash
vntyper cohort --input-file sample_dirs.txt -o cohort_output/
```

Export tabular formats with pseudonymized sample names:

```bash
vntyper cohort -i results/*/ -o cohort_output/ \
    --summary-formats csv,tsv,json --pseudonymize-samples
```

Apply a custom pseudonym prefix:

```bash
vntyper cohort -i results/*/ -o cohort_output/ --pseudonymize-samples patient_
```

Specify custom allele frequency cutoff:

```bash
vntyper cohort -i results/*/ -o cohort_output/ --rare-allele-max-frequency 0.01
```
