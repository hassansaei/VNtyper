# Reference Setup

VNtyper 2 requires reference files before running the pipeline. The `install-references` command downloads and prepares everything automatically.

## What Gets Downloaded

The reference installation includes:

| Component | Description |
|-----------|-------------|
| **Chromosome 1 sequences** | UCSC hg19/hg38 (default), NCBI GRCh37/GRCh38, or Ensembl references |
| **BWA indices** | Alignment indices built locally during installation |
| **MUC1 motif databases** | Pairwise and self-merged MUC1 motif FASTA files with samtools indices |
| **adVNTR VNTR database** | Database files for the optional adVNTR module (requires adding to `--references`) |

## Basic Installation

Download references with default settings (hg19 + hg38, BWA aligner):

```bash
vntyper install-references -d /path/to/references
```

This command downloads reference sequences, builds BWA indices, and indexes motif files. It may take 10--30 minutes depending on your internet connection and CPU.

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `-d`, `--output-dir` | *(required)* | Directory where references will be installed |
| `--references` | `hg19 hg38` | Reference assemblies to download (e.g., `hg19 hg38 GRCh37 GRCh38`) |
| `--aligners` | `bwa` | Aligners to build indices for (e.g., `bwa bwa-mem2 minimap2`) |
| `--skip-indexing` | `false` | Skip BWA index building (useful if indices already exist) |
| `-t`, `--threads` | `4` | Number of threads for indexing |

### Examples

Download only hg38 references:

```bash
vntyper install-references -d ./references --references hg38
```

Download references for both BWA and BWA-MEM2:

```bash
vntyper install-references -d ./references --aligners bwa bwa-mem2
```

Download all supported assemblies including NCBI naming:

```bash
vntyper install-references -d ./references --references hg19 hg38 GRCh37 GRCh38
```

Skip indexing (download sequences only):

```bash
vntyper install-references -d ./references --skip-indexing
```

## Output Directory Structure

After installation, the reference directory looks like this:

```
references/
  alignment/
    chr1.hg19.fa             # hg19 chromosome 1 sequence (extracted)
    chr1.hg19.fa.amb         # BWA index files
    chr1.hg19.fa.ann
    chr1.hg19.fa.bwt
    chr1.hg19.fa.pac
    chr1.hg19.fa.sa
    chr1.hg38.fa             # hg38 chromosome 1 sequence (extracted)
    chr1.hg38.fa.*           # BWA index files
  vntr_db_advntr/            # adVNTR database (only if included in --references)
  All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa
  All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa.fai
  MUC1_motifs_Rev_com.fa
  MUC1_motifs_Rev_com.fa.fai
  code-adVNTR_RUs.fa
  code-adVNTR_RUs.fa.fai
```

## Storage Requirements

Expect approximately 2 GB of disk space for the default installation (hg19 + hg38 references with BWA indices). Adding NCBI or Ensembl references increases this proportionally.

!!! note "MD5 verification"
    VNtyper 2 automatically generates MD5 checksums for all downloaded files, writing them to `md5_checksums.txt` in the output directory. Check `install_references.log` for detailed installation progress.
