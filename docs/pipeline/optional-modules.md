# Optional Modules

VNtyper supports two optional modules to complement Kestrel genotyping: **adVNTR** for independent validation using profile Hidden Markov Models, and **SHARK** for rapid read extraction from large FASTQ datasets.

## adVNTR

### Overview

[adVNTR](https://github.com/mehrdadbakhtiari/adVNTR) uses profile Hidden Markov Models (profile-HMMs) to model repeat structures, identifying variants through read-to-model alignment.

!!! info "Complementary algorithms"
    Kestrel (k-mer de Bruijn graph) and adVNTR (profile-HMM) rely on distinct algorithmic principles. Concordance between the two methods provides high confidence. Discordance indicates complex repeat structures that warrant orthogonal investigation (such as SNaPshot or long-read sequencing).

### Configuration

adVNTR targets the MUC1 coding VNTR using **VNTR ID 25561**. Key settings in `advntr_config.json`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `vid` | 25561 | Database ID for MUC1 VNTR |
| `threads` | 1 | CPU thread allocation |
| `additional_commands` | `""` | Optional adVNTR flags; pass `-aln` to generate alignment sidecars |
| `output_format` | `vcf` | Format of output files (`tsv` or `vcf`) |
| `max_frameshift` | 100 | Maximum frameshift multiplier for filtering |
| `frameshift_multiplier` | 3 | Base multiplier for frame patterns |

### Requirements

- Dedicated conda environment `envadvntr` with adVNTR installed.
- Reference database for the target assembly (hg19 or hg38).

### Processing

VNtyper does not consume adVNTR's optional `.aln` sidecar by default. Users can request it by configuring `additional_commands: "-aln"`.

Frameshift filtering follows MUC1 biological rules:

- **Deletion frameshifts**: Net base change matching `3n + 2` (2, 5, 8, 11, ...).
- **Insertion frameshifts**: Net base change matching `3n + 1` (1, 4, 7, 10, ...).

Variants are annotated with repeat unit (RU) identity and position from adVNTR's state string. Reference and alternate alleles are resolved using the MUC1 RU FASTA reference. If that FASTA fails to resolve, RU and POS remain available while REF and ALT report `Not applicable`.

### Cross-Matching

When both Kestrel and adVNTR run, VNtyper performs pairwise variant cross-matching:

1. Classifies variant type (Insertion, Deletion, or Other) based on allele lengths.
2. Computes the **allele change**: the net inserted or deleted sequence after trimming common prefixes.
3. Evaluates concordance logic (default: variant type and net sequence change must match).

Results are written to `cross_match_results.tsv`.

### Runtime

adVNTR requires approximately 9 minutes per sample. BAM downsampling (`--advntr-max-coverage`) reduces runtime on high-coverage libraries.

## SHARK

### Overview

SHARK performs k-mer-based read filtering to extract MUC1-originating reads directly from paired FASTQ inputs without requiring prior whole-genome or whole-exome alignment.

### Operational Use

Use SHARK when processing whole-genome or whole-exome FASTQ files without an aligned BAM. SHARK filters raw FASTQs directly, avoiding the overhead of aligning full sequencing datasets prior to MUC1 extraction.

!!! tip "Aligned inputs are faster"
    When BAM or CRAM files are available, use `--bam` or `--cram`. Region slicing via `samtools view` is substantially faster than k-mer filtering across raw FASTQ files.

### Requirements

- Conda environment `shark_env` containing the SHARK binary.
- Assembly-specific MUC1 region FASTA files.

### Constraints

- **FASTQ input only**: Aligned BAM/CRAM inputs use samtools slicing instead.
- **Coordinate-specific FASTA selection**: Assembly flags (`--reference-assembly`) select between two coordinate models (GRCh37 vs. GRCh38). Assemblies `hg19`, `GRCh37`, `hg19_ncbi`, and `hg19_ensembl` all map to the GRCh37 FASTA.
- Filtered reads continue through standard fastp QC and BWA alignment.

### Assembly-Specific Reference Selection

SHARK matches k-mers against assembly-specific reference FASTAs (Issue #152). Because 40.6% of canonical 17-mers in the hg38 region are absent in hg19, filtering hg38 reads against an hg19 reference loses reads (retaining 3.2% to 34.7% fewer reads across test cohorts).

Reference resolution in `vntyper/modules/shark/shark_filtering.py` evaluates three tiers:

1. **`config["reference_data"]`**: Paths installed via `vntyper install-references`.
2. **`shark_config.json`**: Shipped defaults (`reference/muc1_region_hg19.fa`, `reference/muc1_region_hg38.fa`).
3. **Legacy flat key `muc1_region_fasta`**: Consulted only if assembly-specific keys are absent.

Explicit `null` entries are treated as intentionally disabled and raise an error rather than falling back.

### Execution

SHARK filters paired FASTQ reads against the region FASTA:

```
shark -r <muc1_region.fa> -1 R1.fastq -2 R2.fastq \
  -o filtered_R1.fastq -p filtered_R2.fastq -t <threads> -k 17 -c 0.6
```

### Verification and Read Pairing

`pipeline_summary.json` logs execution parameters:

- `shark_version`: Binary package version string from conda metadata.
- `shark_k`: K-mer size (default: `17`).
- `shark_c`: Confidence threshold (default: `0.6`).
- `kept_reads_r1` and `kept_reads_r2`: Retained read counts per mate.

The pipeline enforces symmetric read pairing. If `kept_reads_r1 != kept_reads_r2`, the stage fails closed with a `ValueError` before generating summary records.
