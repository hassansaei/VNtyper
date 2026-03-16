# Optional Modules

VNtyper supports two optional modules that complement the core Kestrel genotyping: **adVNTR** for independent validation using a different algorithmic approach, and **SHARK** for rapid read extraction from large FASTQ datasets.

## adVNTR

### Overview

[adVNTR](https://github.com/mehrdadbakhtiari/adVNTR) uses profile Hidden Markov Models (profile-HMMs) to genotype VNTRs. Unlike Kestrel's k-mer-based approach, adVNTR models the repeat structure probabilistically and identifies variants through alignment of reads against trained VNTR models.

!!! info "Complementary approaches"
    Kestrel and adVNTR use fundamentally different algorithms (k-mer graph vs. profile-HMM). Concordance between the two methods provides strong evidence for a true positive call. Discordance warrants further investigation with orthogonal methods such as SNaPshot or long-read sequencing.

### Configuration

adVNTR targets MUC1 VNTR using **VNTR ID 25561**, which corresponds to the MUC1 coding VNTR locus. Key settings from `advntr_config.json`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `vid` | 25561 | VNTR database ID for MUC1 |
| `threads` | 1 | Parallel threads |
| `additional_commands` | `-aln` | Extra flags (alignment mode) |
| `output_format` | `vcf` | Output format (tsv or vcf) |
| `max_frameshift` | 100 | Maximum frameshift multiplier for filtering |
| `frameshift_multiplier` | 3 | Base multiplier for valid frame patterns |

### Requirements

- Conda environment `envadvntr` with adVNTR installed
- adVNTR reference database for the target assembly (hg19 or hg38)

### Processing

adVNTR output is processed through frameshift filtering analogous to Kestrel's:

- **Deletion frameshifts**: frame values matching `3n + 2` (e.g., 2, 5, 8, 11, ...)
- **Insertion frameshifts**: frame values matching `3n + 1` (e.g., 1, 4, 7, 10, ...)

Variants are annotated with repeat unit (RU) identity, position, REF, and ALT using the MUC1 RU FASTA reference. adVNTR-specific flagging rules can be configured independently.

### Cross-Matching

When both Kestrel and adVNTR results are available, VNtyper performs a cross-match comparison. For each pair of variants (one from each caller), the pipeline:

1. Determines variant type (Insertion, Deletion, or Other) based on REF/ALT lengths
2. Computes the **allele change** -- the net inserted or deleted sequence after removing the shared prefix
3. Evaluates a configurable match logic expression (default: allele change and variant type must both match)

The cross-match result (`cross_match_results.tsv`) records all pairwise comparisons and an overall concordance flag ("Yes" if at least one pair matches).

### Runtime

adVNTR genotyping typically requires approximately 9 minutes per sample, significantly longer than Kestrel. Optional BAM downsampling (`--max-coverage`) can reduce runtime for high-coverage samples.

## SHARK

### Overview

[SHARK](https://github.com/AlgoLab/shark) is a rapid read extraction tool that identifies reads likely originating from a target region using k-mer matching against a reference sequence. It operates directly on FASTQ files without requiring alignment.

### When to Use SHARK

SHARK is designed for scenarios where MUC1 coverage is sparse relative to total sequencing volume:

- **Whole-exome sequencing (WES)** -- where MUC1 VNTR reads are a small fraction of total reads
- **Whole-genome sequencing (WGS)** -- where extracting MUC1 reads from a large FASTQ is computationally expensive via alignment

By pre-filtering reads with SHARK, the downstream alignment and genotyping steps operate on a much smaller read set, reducing both runtime and memory usage.

### Requirements

- Conda environment `shark_env` with SHARK installed
- MUC1 region FASTA reference (configured in `shark_config.json`)

### Limitations

- **FASTQ input only** -- SHARK cannot process BAM/CRAM files. For aligned input, the pipeline uses samtools region extraction instead.
- SHARK filtering runs **before** fastp QC, so filtered reads still undergo quality control downstream.

### Execution

SHARK is invoked with paired-end FASTQ input and produces filtered FASTQ files containing only reads matching the MUC1 region:

```
shark -r <muc1_region.fa> -1 R1.fastq -2 R2.fastq \
  -o filtered_R1.fastq -p filtered_R2.fastq -t <threads>
```

The filtered FASTQs replace the original inputs for all subsequent pipeline steps.
