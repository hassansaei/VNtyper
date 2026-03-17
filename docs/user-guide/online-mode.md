# Online Mode

The `online` subcommand submits a BAM file to [vntyper.org](https://vntyper.org) for cloud-based analysis. VNtyper 2 automatically subsets your BAM to the MUC1 region before uploading, so only a small file is transferred.

## Basic Usage

```bash
vntyper online --bam sample.bam -o results/
```

This will:

1. Extract the MUC1 region from your BAM using samtools
2. Upload the subset BAM to the vntyper.org API
3. Poll for job completion
4. Download results to the output directory

## Email Notifications

Receive an email when your job completes:

```bash
vntyper online --bam sample.bam -o results/ \
    --email user@example.com
```

## Cohort Submission

Associate your submission with a cohort on vntyper.org:

```bash
vntyper online --bam sample.bam -o results/ \
    --cohort-id my_study_2024 \
    --passphrase secret123
```

## Resume Polling

If the CLI is interrupted during polling, resume without re-uploading:

```bash
vntyper online --bam sample.bam -o results/ --resume
```

This reads the saved `job_id.txt` from the output directory and resumes status polling.

## Options Reference

| Option | Description |
|--------|-------------|
| `--bam` | Path to input BAM file (required) |
| `-o`, `--output-dir` | Output directory for results (default: `out`) |
| `--reference-assembly` | Assembly used for alignment (default: hg19) |
| `--threads` | Number of threads for BAM subsetting (default: 4) |
| `--email` | Email for job completion notification |
| `--cohort-id` | Cohort ID to associate the job with |
| `--passphrase` | Passphrase for cohort access |
| `--resume` | Resume polling a previously submitted job |

!!! note
    The online mode API endpoint is configured in `config.json` under `api.base_url` (default: `https://vntyper.org/api`).
