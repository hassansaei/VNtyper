# Online Mode

The `online` subcommand submits a BAM file to [vntyper.org](https://vntyper.org) for cloud analysis. VNtyper 2 extracts the MUC1 region locally before upload, transmitting only the sliced target.

## Basic Usage

```bash
vntyper online --bam sample.bam -o results/
```

Workflow steps:

1. Extract the MUC1 region from the BAM using samtools.
2. Upload the subset BAM to the vntyper.org API.
3. Poll the server for job status.
4. Download completed results into the output directory.

## Email Notifications

Receive notification upon job completion:

```bash
vntyper online --bam sample.bam -o results/ \
    --email user@example.com
```

## Cohort Submission

Associate an upload with an existing cohort on vntyper.org:

```bash
vntyper online --bam sample.bam -o results/ \
    --cohort-id 4f9c1a72-5e30-4b8d-9a61-7c2e0d5b83fa \
    --passphrase secret123
```

`--cohort-id` requires the UUID returned during cohort creation (not its human-readable alias). The UUID and passphrase together authenticate reading, downloading, and joint analysis.

!!! note "Request constraints"
    Invalid submissions trigger immediate rejection:

    - **Passphrases are mandatory credentials:** Both identifier and passphrase must accompany every cohort request.
    - **Cohorts lacking passphrases remain locked:** Missing passphrase hashes prevent access; such records expire automatically by TTL.
    - **Identifiers must originate from the service:** Cohort aliases are limited to 64 ASCII characters without control codes.
    - **Passphrases are limited to 72 UTF-8 bytes:** Exceeding 72 bytes triggers an error.
    - **Filenames require plain ASCII:** Permitted characters comprise alphanumeric symbols, dots, hyphens, underscores, and plus signs, ending with `.bam`, `.cram`, `.bai`, or `.crai`.
    - **Server errors return generic messages:** Failure details write to server logs; reference the job ID when reporting issues.
    - **Server deployments require `REDIS_PASSWORD`:** The backend service halts without this environment setting.

## Resume Polling

Resume interrupted polling without re-uploading the BAM:

```bash
vntyper online --bam sample.bam -o results/ --resume
```

This reads `job_id.txt` from the output directory and resumes status checks.

## Options Reference

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--bam` | path | (required) | Path to input BAM file |
| `-o, --output-dir` | path | `out` | Directory for downloaded results |
| `--reference-assembly` | choice | `hg19` | Reference assembly used for alignment |
| `--threads` | int | `4` | Thread count for local BAM extraction |
| `--email` | string | None | Email address for completion notice |
| `--cohort-id` | string | None | Target cohort UUID |
| `--passphrase` | string | None | Cohort passphrase |
| `--resume` | flag | off | Resume polling using local `job_id.txt` |

!!! note
    API endpoints configure under `api.base_url` in `config.json` (default: `https://vntyper.org/api`).
