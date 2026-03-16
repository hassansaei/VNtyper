# Quick Start

Run your first VNtyper analysis in 5 minutes.

## 1. Install VNtyper

Follow the [Installation guide](installation.md) to install VNtyper and its dependencies.

## 2. Download Reference Files

VNtyper needs reference sequences and motif databases before it can run. Download them to a local directory:

```bash
vntyper install-references -d ./references
```

This downloads chromosome 1 references (hg19/hg38), MUC1 motif databases, and adVNTR VNTR databases, then builds BWA indices. See [Reference Setup](reference-setup.md) for details.

## 3. Run the Pipeline

Analyze a BAM file with the default Kestrel genotyping engine:

```bash
vntyper pipeline \
    --bam sample.bam \
    -o results/ \
    --threads 4 \
    --reference-assembly hg38
```

!!! tip "Don't have a BAM file?"
    Download the VNtyper test dataset (~1.1 GB) to try it out:

    ```bash
    make download-test-data
    ```

    Then run the pipeline on the test BAM file located in the test data directory.

For paired-end FASTQ input:

```bash
vntyper pipeline \
    --fastq1 sample_R1.fastq.gz \
    --fastq2 sample_R2.fastq.gz \
    -o results/ \
    --threads 4
```

Add `--fast-mode` to skip optional steps and speed up the analysis.

## 4. View Results

Once the pipeline completes, the output directory contains:

```
results/
  pipeline.log              # Full pipeline log
  kestrel_result.tsv        # Genotyping results (main output)
  alignment/                # Aligned BAM files
  fastq/                    # Extracted FASTQ reads
  kestrel/                  # Kestrel intermediate files
```

The primary output is `kestrel_result.tsv`, which contains detected MUC1 VNTR variants with confidence scores, frame shift analysis, and depth metrics.

## 5. Generate an HTML Report

Create a visual summary report with IGV integration:

```bash
vntyper report \
    -o results/ \
    --input-dir results/
```

Open the generated HTML file in your browser to review:

- VNTR region coverage statistics
- Genotyping calls from Kestrel
- Quality metrics (duplication rate, Q20/Q30 rates)
- Pipeline execution log

## What's Next?

- **[Reference Setup](reference-setup.md)** --- Configure references for different genome assemblies
- **[User Guide](../index.md)** --- Explore advanced pipeline options, optional modules (adVNTR, SHARK), and Docker usage
