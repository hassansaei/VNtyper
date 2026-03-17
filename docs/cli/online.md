# online

Subset a BAM file and submit it to an online VNtyper 2 instance (vntyper.org) for analysis, then poll for completion and download results.

## Synopsis

```
vntyper [global-options] online
    --bam <file>
    [-o <dir>]
    [--reference-assembly <assembly>]
    [--threads <int>]
    [--email <address>]
    [--cohort-id <id>] [--passphrase <string>]
    [--resume]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--bam` | path | (required) | Path to the input BAM file |
| `-o, --output-dir` | path | `out` | Output directory for results |
| `--reference-assembly` | choice | `hg19` | Reference assembly used for alignment. Options: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ncbi`, `hg38_ncbi`, `hg19_ensembl`, `hg38_ensembl` |
| `--threads` | int | `4` | Number of threads to use |
| `--email` | string | — | Email address to receive notifications (optional) |
| `--cohort-id` | string | — | Cohort ID to associate the job with (optional) |
| `--passphrase` | string | — | Passphrase for the cohort, if required |
| `--resume` | flag | off | Resume polling a previously submitted job if a `job_id` is found in the output directory |

## Workflow

The `online` command follows a submit-poll-download workflow:

1. **Subset:** The input BAM is subsetted to the MUC1 VNTR region based on the specified reference assembly.
2. **Submit:** The subsetted BAM is uploaded to the vntyper.org server for analysis.
3. **Poll:** The command periodically checks the server for job completion.
4. **Download:** Once complete, results are downloaded to the output directory.

Use `--resume` to skip submission and resume polling a previously submitted job.

## Examples

Submit a BAM for online analysis:

```bash
vntyper online --bam sample.bam -o results/ --reference-assembly hg38
```

Submit with email notifications and cohort association:

```bash
vntyper online --bam sample.bam -o results/ \
    --email user@example.com --cohort-id my_cohort --passphrase secret123
```

Resume polling a previously submitted job:

```bash
vntyper online --bam sample.bam -o results/ --resume
```
