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

## Physical identity: eight labels, six files

The eight assembly labels above map to six physical reference files. Contig naming differs by reference source (UCSC `chr1` vs. NCBI `NC_000001.10` vs. Ensembl `1`), so BWA and CRAM reference files are keyed by `(coordinate system, source)`. `GRCh37` and `hg19_ncbi` name the same NCBI RefSeq FASTA for GRCh37; `GRCh38` and `hg38_ncbi` name the identical NCBI RefSeq FASTA for GRCh38.

| Label(s) | Physical file id | Config key VNtyper writes and reads |
|----------|-------------------|--------------------------------------|
| `hg19` | `hg19` | `bwa_reference_hg19` |
| `GRCh37`, `hg19_ncbi` | `GRCh37` | `bwa_reference_GRCh37` |
| `hg19_ensembl` | `hg19_ensembl` | `bwa_reference_hg19_ensembl` |
| `hg38` | `hg38` | `bwa_reference_hg38` |
| `GRCh38`, `hg38_ncbi` | `GRCh38` | `bwa_reference_GRCh38` |
| `hg38_ensembl` | `hg38_ensembl` | `bwa_reference_hg38_ensembl` |

`vntyper install-references` writes the `bwa_reference_*` key for each physical file installed (see [install-references](../cli/install-references.md)). Default installation provisions `bwa_reference_hg19` and `bwa_reference_hg38`. Adding `--references GRCh38` installs the physical genome shared by `GRCh38` and `hg38_ncbi`.

adVNTR and SHARK require only two target files: one adVNTR database and one MUC1 region FASTA per coordinate system. Contig naming is irrelevant once VNtyper slices out the MUC1 region. These are keyed strictly by coordinate system (`advntr_reference_vntr_hg19`/`_hg38`, `muc1_region_fasta_hg19`/`_hg38`). All four labels sharing a coordinate system (`hg19`, `GRCh37`, `hg19_ncbi`, `hg19_ensembl`) resolve to identical adVNTR and SHARK references.

Keying by coordinate system prevents cross-build assignment: every label resolves to its coordinate system first, preventing unmapped targets from defaulting to hg19.

### The fallback, and why a complete installation never uses it

Reference resolution checks keys in order of specificity: exact assembly label (if specialized, such as `bwa_reference_hg19_ncbi`), physical file key from the table above, then UCSC coordinate key (`bwa_reference_hg19` or `bwa_reference_hg38`). A key explicitly defined as `null` indicates a deliberately disabled reference and halts fallback.

If no matching key exists, VNtyper falls back to the UCSC key and logs a warning identifying the substituted assembly. A complete `vntyper install-references` setup never triggers this fallback. The warning indicates a missing reference index that requires installation.

## MUC1 VNTR Region Coordinates

Assemblies within the same coordinate system share identical coordinates:

| Coordinate System | BAM Extraction Region | VNTR Region |
|-------------------|-----------------------|-------------|
| GRCh37 | 155158000-155163000 | 155160500-155162000 |
| GRCh38 | 155184000-155194000 | 155188000-155192500 |

The BAM extraction region spans flanking sequence to retain boundary reads. The VNTR region defines the target evaluated during coverage calculations.

## Auto-Detection from BAM Headers

For BAM or CRAM input, VNtyper 2 inspects the file header to determine contig naming and constructs the extraction region dynamically. An Ensembl BAM (contigs `1`, `2`) run with `--reference-assembly hg19` targets `1:155158000-155163000`.

## When to Use `--reference-assembly`

```bash
vntyper pipeline --bam inputs/sample.bam -o results/sample/ --reference-assembly hg38
```

Specify this option when:

- Input alignments map to hg38/GRCh38 (VNtyper 2 defaults to hg19).
- Explicit chromosome naming is required (such as `hg19_ensembl`).
- Running with FASTQ input, where no BAM header exists for auto-detection.

!!! note
    The assembly choice dictates the BWA reference index and coordinate targets. Declaring an incompatible assembly invalidates genotyping results.

## The declared-assembly check

BAM and CRAM inputs are checked against `--reference-assembly` before region slicing. The header chr1 length determines the true build. Disagreement with the declared assembly halts execution with an error indicating expected values.

Hybrid headers containing duplicate conflicting contigs (for example, contig `1` and `chr1` possessing discordant lengths from two builds) also halt execution.

This validation prevents false-negative calls: because MUC1 positions differ by approximately 30 kb between GRCh37 and GRCh38, an incorrect assembly extracts off-target sequence. Kestrel finds no variant k-mers, erroneously generating an exit-code-0 negative call.

Specific check boundaries:

- **Aliases do not disagree:** `hg19`, `GRCh37`, `hg19_ncbi`, and `hg19_ensembl` all map to GRCh37 coordinates and validate against GRCh37 headers.
- **Undetermined headers log warnings:** Headers missing chr1 or carrying unknown lengths log warnings and proceed without halting.
- **FASTQ input bypasses the check:** Raw FASTQ files lack alignment headers; alignment coordinates follow the selected reference directly.

## Installing References

Install reference bundles before executing the pipeline. Standard runs download pre-indexed, checksummed releases from [`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data):

```bash
vntyper install-references -d /path/to/references
```

To install specific physical targets:

```bash
vntyper install-references -d /path/to/references --references hg19 hg38 GRCh37 GRCh38
```

Selecting `--references hg19` installs common MUC1 motif FASTAs and both adVNTR databases alongside the assembly index.
