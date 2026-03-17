# Snakemake Batch Processing

VNtyper 2 includes a Snakemake workflow for processing multiple BAM files in parallel, with support for cluster execution.

## Setup

The workflow file is located at `snakemake/vntyper2.smk`.

### Create the Input File

List your BAM file paths in `bams.txt`, one per line:

```text title="bams.txt"
/data/samples/sample1.bam
/data/samples/sample2.bam
/data/samples/sample3.bam
```

Each BAM filename (without `.bam`) becomes the sample name and output subdirectory name.

### Create Output Directories

The workflow writes results to `results/<sample_name>/` and logs to `logs/<sample_name>.log`.

## Running Locally

```bash
snakemake -s snakemake/vntyper2.smk --use-conda -j 10
```

| Flag | Purpose |
|------|---------|
| `-s` | Path to the Snakemake workflow file |
| `--use-conda` | Activate the `vntyper` conda environment |
| `-j 10` | Run up to 10 jobs in parallel |

## Cluster Execution

For SLURM or other cluster schedulers, use a Snakemake profile:

```bash
snakemake -s snakemake/vntyper2.smk --use-conda --profile slurm/
```

### Resource Allocation

The workflow defines resource requests per job:

| Resource | Value |
|----------|-------|
| Threads | 8 per job |
| Memory | 2.2 GB per thread (17.6 GB total) |
| Time limit | 72 hours |

These are set in the `run_vntyper_pipeline` rule and can be adjusted by editing `snakemake/vntyper2.smk`.

## Workflow Details

The workflow runs one rule per BAM file:

```bash
vntyper --config-path vntyper/config.json pipeline \
    --bam {input.bam} --threads 8 --reference-assembly hg38 \
    --fast-mode --keep-intermediates -o results/{sample}/
```

!!! tip "Customizing the Pipeline Command"
    Edit the `shell` directive in `snakemake/vntyper2.smk` to change the reference assembly, enable adVNTR, or modify other pipeline options.

## Output Structure

```
results/
├── sample1/
│   ├── kestrel/
│   │   └── kestrel_result.tsv
│   ├── pipeline_summary.json
│   └── ...
├── sample2/
│   └── ...
logs/
├── sample1.log
├── sample2.log
```

After batch processing, use the [cohort analysis](cohort-analysis.md) command to aggregate results:

```bash
vntyper cohort -i results/* -o cohort_output/
```
