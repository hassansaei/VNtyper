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

The eight assembly labels above do not each name a separate genome file. Contig naming
differs by *reference source* (UCSC `chr1` vs. NCBI `NC_000001.10` vs. Ensembl `1`), so the
BWA and CRAM reference files are keyed by `(coordinate system, source)` — but `GRCh37` and
`hg19_ncbi` are two labels for the **same** NCBI file, as are `GRCh38` and `hg38_ncbi`: both
name the NCBI RefSeq FASTA for their coordinate system, just spelled differently. That
leaves six physical files behind the eight labels:

| Label(s) | Physical file id | Config key VNtyper writes and reads |
|----------|-------------------|--------------------------------------|
| `hg19` | `hg19` | `bwa_reference_hg19` |
| `GRCh37`, `hg19_ncbi` | `GRCh37` | `bwa_reference_GRCh37` |
| `hg19_ensembl` | `hg19_ensembl` | `bwa_reference_hg19_ensembl` |
| `hg38` | `hg38` | `bwa_reference_hg38` |
| `GRCh38`, `hg38_ncbi` | `GRCh38` | `bwa_reference_GRCh38` |
| `hg38_ensembl` | `hg38_ensembl` | `bwa_reference_hg38_ensembl` |

`vntyper install-references` writes the `bwa_reference_*` key for each physical file it
actually installed, out of these six possible keys (see
[install-references](../cli/install-references.md)) -- the default run installs only
`bwa_reference_hg19` and `bwa_reference_hg38`; `--references GRCh38` additionally installs
the physical file that both `GRCh38` and `hg38_ncbi` runs read.

adVNTR and SHARK need no such distinction: only two files of each exist in total — one
adVNTR database and one MUC1 region FASTA per **coordinate system** — because contig naming
is irrelevant once VNtyper has already sliced out the MUC1 region. Those are keyed by
coordinate system alone (`advntr_reference_vntr_hg19`/`_hg38`,
`muc1_region_fasta_hg19`/`_hg38`), so all four labels sharing a coordinate system — `hg19`,
`GRCh37`, `hg19_ncbi`, `hg19_ensembl` — resolve to the same adVNTR database and the same
SHARK region FASTA.

Before this was fixed, `pipeline.py` held its own four-entry map from assembly to adVNTR
database and defaulted every label it did not recognise to hg19 — so `--reference-assembly
hg38_ncbi` or `hg38_ensembl` silently validated against the **GRCh37** adVNTR database
instead of GRCh38's. Keying by coordinate system instead of by the raw label closes that:
every label maps to its coordinate system first, so there is no "unrecognised label" case
left to default anywhere.

### The fallback, and why a complete installation never uses it

Reading a reference walks a fixed key order, most specific first: the exact label (when the
config specialises one, e.g. `bwa_reference_hg19_ncbi`), then the physical-file key from the
table above, then the UCSC-family key for the coordinate system (`bwa_reference_hg19` /
`bwa_reference_hg38`). A key that is *present* wins even when its value is `null` — that is
a deliberate "disabled", not a miss to fall through past.

If none of those keys is present, VNtyper falls back to the UCSC-family key and logs a
warning naming the assembly and the key it used instead, because that run is now using
`hg19`/`hg38` **UCSC** sequence for a request that named a different source. This is not a
normal outcome: `vntyper install-references` writes the physical-file key for every
reference it installs, so a complete installation never triggers this fallback — it only
fires against a `config.json` that was hand-edited or only partially installed. Seeing the
warning is the signal to run `vntyper install-references` for the missing assembly rather
than to trust the fallback file.

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

Before running the pipeline, install reference files. By default this fetches the
published, checksummed reference release from
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data) rather than building
anything locally — see [Reference Setup](../getting-started/reference-setup.md) for the
trust model and [install-references](../cli/install-references.md) for every flag.

```bash
vntyper install-references -d /path/to/references
```

To install specific physical files (the six ids from the table above, not the eight
labels — `GRCh37` also covers `hg19_ncbi`, `GRCh38` also covers `hg38_ncbi`):

```bash
vntyper install-references -d /path/to/references --references hg19 hg38 GRCh37 GRCh38
```

`--references hg19` alone still installs the common MUC1 motif FASTAs and both adVNTR
databases — those are not selectable per assembly, because only one of each exists.
