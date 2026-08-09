# Reference Assemblies

VNtyper 2 supports multiple reference genome assemblies with automatic chromosome naming detection.

## Supported Assemblies

| Assembly Name | Coordinate System | Chromosome Naming | Example Chr1 |
|--------------|-------------------|-------------------|--------------|
| `hg19` | GRCh37 | UCSC | `chr1` |
| `hg38` | GRCh38 | UCSC | `chr1` |
| `GRCh37` | GRCh37 | NCBI RefSeq | `NC_000001.10` |
| `GRCh38` | GRCh38 | NCBI RefSeq | `NC_000001.11` |
| `hg19_ncbi` | GRCh37 | NCBI RefSeq | `NC_000001.10` |
| `hg38_ncbi` | GRCh38 | NCBI RefSeq | `NC_000001.11` |
| `hg19_ensembl` | GRCh37 | Ensembl | `1` |
| `hg38_ensembl` | GRCh38 | Ensembl | `1` |

## MUC1 VNTR Region Coordinates

All assemblies within the same coordinate system use identical coordinates:

| Coordinate System | BAM Extraction Region | VNTR Region |
|-------------------|-----------------------|-------------|
| GRCh37 | 155158000--155163000 | 155160500--155162000 |
| GRCh38 | 155184000--155194000 | 155188000--155192500 |

The BAM extraction region is wider to capture flanking reads. The VNTR region is the precise target used for coverage calculation.

## Auto-Detection from BAM Headers

When processing BAM or CRAM input, VNtyper 2 detects the chromosome naming convention from the file header and constructs the correct region string automatically. For example, a BAM aligned to an Ensembl reference (chromosomes named `1`, `2`, ...) with `--reference-assembly hg19` will produce the region `1:155158000-155163000`.

## When to Use `--reference-assembly`

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ --reference-assembly hg38
```

Specify this option when:

- Your BAM is aligned to **hg38/GRCh38** (VNtyper 2 defaults to hg19)
- You want to explicitly select the chromosome naming convention (e.g., `hg19_ensembl` for Ensembl-style names)
- You are using FASTQ input (no BAM header available for auto-detection)

!!! note
    The assembly choice affects which BWA reference index is used for alignment and which genomic coordinates are used for read extraction. Using the wrong assembly will produce incorrect results.

## The declared-assembly check

BAM and CRAM input is checked against `--reference-assembly` before any region is
resolved. The chr1 length in the alignment header decides which build the file actually
describes; if that disagrees with what you declared, the run **stops with an error**
naming both builds and the value you should have used.

A header can also disagree with *itself*: a hybrid header naming chromosome 1 twice,
once as `1` and once as `chr1`, at two lengths belonging to two different builds. That
also **stops the run**, because whichever alias the region string resolves to, one of
the two contigs does not carry MUC1 where the declared build says it does. No value of
`--reference-assembly` reconciles such a file; re-generate it against a single
reference.

This matters because the failure it prevents is invisible. The MUC1 VNTR sits about
30 kb apart between the two builds, so declaring the wrong one extracts a region that
does not contain the VNTR. Kestrel then finds no supporting reads, and the report says
the sample is negative — a confident, wrong answer with exit code 0.

Three things the check deliberately does **not** do:

- **Aliases are not disagreements.** `hg19`, `GRCh37`, `hg19_ncbi` and `hg19_ensembl`
  all name the GRCh37 coordinate system, and any of them agrees with a GRCh37 header.
  The contig *naming* convention is reported but never used to decide.
- **It does not guess.** If the header cannot be read, carries no chr1, or carries a
  chr1 whose length matches neither build, the result is *undetermined*: VNtyper logs a
  warning and continues. Undetermined is neither a pass nor a failure — it means the
  question could not be answered, which is not the same as the header answering it
  contradictorily.
- **It does not apply to FASTQ input.** A FASTQ has no header of its own, and the header
  produced after alignment describes the BWA reference rather than the sample.

## Installing References

Before running the pipeline, install reference files:

```bash
vntyper install-references -d /path/to/references --threads 4
```

To install specific assemblies or aligners:

```bash
vntyper install-references -d /path/to/references \
    --references hg19 hg38 \
    --aligners bwa bwa-mem2
```
