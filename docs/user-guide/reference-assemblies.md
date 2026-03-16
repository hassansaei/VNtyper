# Reference Assemblies

VNtyper supports multiple reference genome assemblies with automatic chromosome naming detection.

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

When processing BAM or CRAM input, VNtyper detects the chromosome naming convention from the file header and constructs the correct region string automatically. For example, a BAM aligned to an Ensembl reference (chromosomes named `1`, `2`, ...) with `--reference-assembly hg19` will produce the region `1:155158000-155163000`.

## When to Use `--reference-assembly`

```bash
vntyper pipeline --bam sample.bam -o results/ --reference-assembly hg38
```

Specify this option when:

- Your BAM is aligned to **hg38/GRCh38** (VNtyper defaults to hg19)
- You want to explicitly select the chromosome naming convention (e.g., `hg19_ensembl` for Ensembl-style names)
- You are using FASTQ input (no BAM header available for auto-detection)

!!! note
    The assembly choice affects which BWA reference index is used for alignment and which genomic coordinates are used for read extraction. Using the wrong assembly will produce incorrect results.

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
