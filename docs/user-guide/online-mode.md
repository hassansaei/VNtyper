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
    --cohort-id 4f9c1a72-5e30-4b8d-9a61-7c2e0d5b83fa \
    --passphrase secret123
```

`--cohort-id` is the identifier the service returned when the cohort was
created — not the alias. An alias is a label: it is unique, and it is checked
against the cohort when you supply it, but it never selects a cohort on its own.
Keep the identifier the creation response gave you; it and the passphrase
together are what open the cohort for reading, downloading and joint analysis.

!!! note "The rules a request has to satisfy"
    These are the constraints on cohorts and uploads. None of them affects a
    submission that is well-formed; each is listed because a request that breaks
    one is refused rather than repaired.

    - **A passphrase is required, and it is the only credential.** Both the
      identifier and the passphrase must be supplied on every cohort request.
      An alias is a label, not a credential.
    - **A cohort with no stored passphrase cannot be opened.** An absent hash is
      read as "cannot be opened", never as "no passphrase needed"; such records
      age out on their own retention TTL.
    - **A cohort identifier must be one the service issued**, and an alias may be
      at most 64 characters with no control characters. Both become part of a
      storage key, so their shape is fixed.
    - **A passphrase may be at most 72 bytes once UTF-8 encoded.** A longer one
      is refused with an explicit message rather than shortened, and the alias
      you asked for stays free for another attempt.
    - **Upload filenames must be plain ASCII.** Letters, digits, dot, dash,
      underscore and plus, starting with a letter or a digit, ending in
      `.bam`/`.cram` (or `.bai`/`.crai` for an index). The extension may be in
      any case, and the name is stored exactly as you sent it. A name outside
      that set is refused rather than repaired.
    - **A failed job reports that it failed, and nothing more.** The command
      line and container paths behind the failure are written to the server log;
      quote the job ID when reporting one.
    - **Operators:** the service will not start without `REDIS_PASSWORD`, and
      there is no default.

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
| `--cohort-id` | Cohort ID to associate the job with, as returned when the cohort was created |
| `--passphrase` | Passphrase for cohort access, required whenever a cohort is named |
| `--resume` | Resume polling a previously submitted job |

!!! note
    The online mode API endpoint is configured in `config.json` under `api.base_url` (default: `https://vntyper.org/api`).
