# Quick Start

Run your first VNtyper 2 analysis in 5 minutes.

## 1. Install VNtyper 2

Follow the [Installation guide](installation.md) to install VNtyper 2 and its dependencies.

## 2. Download Reference Files

VNtyper 2 needs reference sequences and motif databases before it can run. Download them to a local directory:

```bash
vntyper install-references -d ./references
```

This downloads chromosome 1 references (hg19/hg38) and MUC1 motif databases, then builds BWA indices. See [Reference Setup](reference-setup.md) for details.

## 3. Run the Pipeline

Analyze a BAM file with the default Kestrel genotyping engine:

```bash
vntyper pipeline \
    --bam inputs/sample.bam \
    -o results/sample/ \
    --threads 4 \
    --reference-assembly hg38
```

For BAM and CRAM input, the output root must be outside the directory containing the
alignment. The separate `inputs/` and `results/` trees above meet that requirement.

!!! tip "Test Dataset"
    Download the VNtyper 2 test dataset (~1.1 GB) to test the pipeline:

    ```bash
    make download-test-data
    ```

    Then run the pipeline on the test BAM file located in `tests/test_data/`.

For paired-end FASTQ input:

```bash
vntyper pipeline \
    --fastq1 sample_R1.fastq.gz \
    --fastq2 sample_R2.fastq.gz \
    -o results/ \
    --threads 4
```

Add `--fast-mode` to skip filtering for unmapped and partially mapped reads, speeding up alignment extraction.

## 4. View Results

The pipeline automatically generates an HTML summary report and tabular outputs in the output directory:

```
results/sample/
  pipeline.log                  # Full pipeline log
  pipeline_summary.json         # Machine-readable summary
  summary_report.html           # Interactive HTML report with embedded IGV
  kestrel/
    kestrel_result.tsv          # Genotyping results (main output)
    output_indel.vcf            # Filtered INDEL VCF
    output.bam                  # Kestrel resolved haplotype records
  fastq_bam_processing/         # Extracted FASTQ reads
  alignment_processing/         # BWA-aligned BAM (FASTQ input)
  coverage/                     # Coverage statistics
```

The primary output is `kestrel/kestrel_result.tsv`, which contains detected MUC1 VNTR variants with confidence scores, frameshift analysis, and depth metrics.

Kestrel's `output.bam` is an alignment of resolved haplotype records; it is not a second
copy of the input sequencing reads. The HTML report explains the evidence units beside
its nomenclature reading key.

## 5. Regenerate HTML Reports

The pipeline writes `summary_report.html` by default. To regenerate or customize the report without rerunning the pipeline:

```bash
vntyper report \
    -o results/sample/ \
    --input-dir results/sample/
```

Open the generated HTML file in a web browser to inspect:

- VNTR region coverage statistics
- Genotyping calls from Kestrel
- Quality metrics (duplication rate, Q20/Q30 rates)
- Pipeline execution log

## What's Next?

- **[Reference Setup](reference-setup.md)**: Configure references for different genome assemblies
- **[User Guide](../user-guide/index.md)**: Explore advanced pipeline options, optional modules (adVNTR, SHARK), and Docker usage
