# online

Subset a BAM file to the MUC1 region, upload to vntyper.org, poll for completion, and retrieve output files.

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
| `--bam` | path | (required) | Path to input BAM file |
| `-o, --output-dir` | path | `out` | Directory for downloaded results |
| `--reference-assembly` | choice | `hg19` | Reference assembly: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ncbi`, `hg38_ncbi`, `hg19_ensembl`, `hg38_ensembl` |
| `--threads` | int | `4` | Thread count for local BAM extraction |
| `--email` | string | None | Email address for job completion notice |
| `--cohort-id` | string | None | Cohort UUID provided upon cohort creation (not its alias) |
| `--passphrase` | string | None | Cohort passphrase (mandatory when `--cohort-id` is specified; maximum 72 UTF-8 bytes) |
| `--resume` | flag | off | Resume polling for an existing job using local `job_id.txt` |

## Workflow

The command coordinates four stages:

1. **Subset:** Extracts the MUC1 region using local samtools according to the chosen reference assembly.
2. **Submit:** Posts the sliced BAM to vntyper.org.
3. **Poll:** Queries the status endpoint periodically.
4. **Download:** Retrieves results archive upon completion.

Use `--resume` to reconnect to an in-flight job without repeating upload.

## Exit codes and polling limits

`vntyper online` exits **1** when a remote job fails, terminates unexpectedly, returns no job ID, or exceeds the polling limit. Before VNtyper 2.0.6 all of those logged a message and exited **0**, so a wrapping `subprocess.run(..., check=True)` treated a failed genotyping run as a success.

Polling limits are configured in `config.json`:

| Key | Default | Description |
|-----|---------|-------------|
| `poll_interval_seconds` | `10` | Interval between status requests |
| `poll_timeout_seconds` | `14400` (4 h) | Polling cutoff before timeout |

A timeout indicates that the remote job is still executing. Re-issue the command with `--resume` to continue monitoring.

Job failure details are kept in server logs; reference the job ID when requesting administrative review.

## Examples

Submit a BAM file:

```bash
vntyper online --bam sample.bam -o results/ --reference-assembly hg38
```

Submit with email notification and cohort association:

```bash
vntyper online --bam sample.bam -o results/ \
    --email user@example.com \
    --cohort-id 4f9c1a72-5e30-4b8d-9a61-7c2e0d5b83fa --passphrase secret123
```

Resume polling for an active submission:

```bash
vntyper online --bam sample.bam -o results/ --resume
```
