# Optional Modules

VNtyper 2 supports two optional modules that complement the core Kestrel genotyping: **adVNTR** for independent validation using a different algorithmic approach, and **SHARK** for rapid read extraction from large FASTQ datasets.

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

When both Kestrel and adVNTR results are available, VNtyper 2 performs a cross-match comparison. For each pair of variants (one from each caller), the pipeline:

1. Determines variant type (Insertion, Deletion, or Other) based on REF/ALT lengths
2. Computes the **allele change** -- the net inserted or deleted sequence after removing the shared prefix
3. Evaluates a configurable match logic expression (default: allele change and variant type must both match)

The cross-match result (`cross_match_results.tsv`) records all pairwise comparisons and an overall concordance flag ("Yes" if at least one pair matches).

### Runtime

adVNTR genotyping typically requires approximately 9 minutes per sample, significantly longer than Kestrel. Optional BAM downsampling (`--advntr-max-coverage`) can reduce runtime for high-coverage samples.

## SHARK

### Overview

The SHARK module in VNtyper 2 is a re-implementation of the SHARK concept for MUC1-targeted read extraction. The original SHARK article (Denti et al., *Bioinformatics* 2021) does not provide a publicly available code repository. VNtyper 2's SHARK implementation identifies reads likely originating from the MUC1 region using k-mer matching against a reference sequence, operating directly on FASTQ files without requiring alignment.

### When to Use SHARK

SHARK is for when you only have FASTQ files (no BAM/CRAM) and want to avoid processing entire exome or genome FASTQs through the pipeline. Instead of aligning all reads and then extracting the MUC1 region, SHARK extracts MUC1-relevant reads directly from the raw FASTQs before any alignment occurs.

This is the typical scenario:

- You have **whole-exome or whole-genome FASTQ files** and no aligned BAM
- You want to **skip aligning the full dataset** just to extract MUC1 reads

!!! tip "BAM input is always faster"
    If you have an aligned BAM file, use `--bam` instead. BAM mode extracts MUC1 reads via samtools region slicing, which is much faster than SHARK. SHARK is only useful when BAM files are not available.

### Requirements

- Conda environment `shark_env` with SHARK installed
- MUC1 region FASTA references (configured in `shark_config.json`, one per coordinate system)

### Limitations

- **FASTQ input only** -- SHARK cannot process BAM/CRAM files. For aligned input, the pipeline uses samtools region extraction instead.
- **`--reference-assembly` selects between two region FASTAs, not eight** -- see [Reference assembly selects the region](#reference-assembly-selects-the-region) below. It distinguishes coordinate system (GRCh37 vs. GRCh38), not chromosome-naming source, so `hg19`, `GRCh37`, `hg19_ncbi` and `hg19_ensembl` all select the same GRCh37-coordinate FASTA.
- SHARK filtering runs **before** fastp QC, so filtered reads still undergo quality control downstream.
- After SHARK filtering, the pipeline still performs BWA alignment and full postprocessing on the filtered reads.

### Reference assembly selects the region

SHARK filters reads by matching k-mers against a MUC1 region FASTA, and which FASTA it uses
now follows `--reference-assembly`'s coordinate system (issue #152). This replaces an
earlier decision ([#187](https://github.com/hassansaei/VNtyper/issues/187)) to keep a single
hg19-based region FASTA for both assemblies -- reopened once the cost of that shortcut was
measured rather than assumed: **40.6% of the hg38 region's canonical 17-mers are absent from
the hg19 region**, and filtering hg38-appropriate reads against the hg19 region measurably
loses reads -- across the seven `tests/data/` cohort samples, the hg38 region retains
**3.2--34.7% more reads** than the hg19 region on the same input.

`select_muc1_region_fasta()` in `vntyper/modules/shark/shark_filtering.py` resolves the
region FASTA in three tiers, most authoritative first:

1. **`config["reference_data"]`**, keyed `muc1_region_fasta_hg19` / `muc1_region_fasta_hg38`
   -- what `vntyper install-references` writes. This is consulted first so that references
   installed into a custom `--output-dir` are honoured: `--config-path` replaces the main
   `config.json` but never touches `shark_config.json`, so without this layer an installed
   tree would be invisible to SHARK.
2. **`shark_config.json`'s `shark_settings`**, keyed the same way -- the shipped default,
   pointing at `reference/muc1_region_hg19.fa` and `reference/muc1_region_hg38.fa`.
3. **The legacy flat `muc1_region_fasta` key**, used only when `shark_settings` carries
   **no** `muc1_region_fasta_*` entry at all. A config with one keyed entry but not the
   other is treated as an incomplete keyed config, not as pre-#152 legacy, and is not
   silently patched from the flat key.

Resolution at every tier is by key **membership**: a key present with value `null` is a
deliberate "disabled" for that assembly and raises rather than falling through to the next
tier or the legacy key. See [Configuration](../user-guide/configuration.md) and
[Reference Assemblies](../user-guide/reference-assemblies.md#the-fallback-and-why-a-complete-installation-never-uses-it)
for the same membership rule as it applies to BWA and adVNTR reference selection.

### Execution

SHARK is invoked with paired-end FASTQ input and produces filtered FASTQ files containing only reads matching the MUC1 region:

```
shark -r <muc1_region.fa> -1 R1.fastq -2 R2.fastq \
  -o filtered_R1.fastq -p filtered_R2.fastq -t <threads>
```

The filtered FASTQs replace the original inputs for all subsequent pipeline steps.
