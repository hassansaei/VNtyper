# install-references

Download and set up the reference files required by VNtyper 2.

By default, this fetches the published, checksummed reference release from
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data) rather than building
each reference from its original upstream host. See
[Reference Setup](../getting-started/reference-setup.md) for the trust model, what the
release publishes, and what `--from-source` costs.

## Synopsis

```
vntyper [global-options] install-references
    -d <dir>
    [--skip-indexing]
    [-t <int>]
    [--aligners <aligner> ...]
    [--references <reference> ...]
    [--from-source]
    [--release-spec <path>]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-d, --output-dir` | path | (required) | Directory where references will be installed |
| `--skip-indexing` | flag | off | Skip the indexing step. Only affects `--from-source` — the default bundle path re-indexes only in reaction to a detected BWA version mismatch |
| `-t, --threads` | int | `4` | Threads for indexing. Only affects `--from-source` |
| `--aligners` | list | `bwa` only | `--from-source` only: aligners to build indices for (e.g., `bwa`, `bwa-mem2`, `minimap2`) |
| `--derive-only` | flag | off | Rebuild the derived reference files (`muc1_region_hg19.fa`, `muc1_region_hg38.fa`, `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`) from genomes and seeds already installed. Downloads nothing. A derivation this tree cannot rebuild — its source genome or its seeds absent — is **skipped rather than failed**, and whatever is already at its path is verified instead, so every one of the three is checked against its committed digest either way. On a bundle-installed tree that means the two region FASTAs are rebuilt and the motif FASTA is verified in place, since the bundle does not stage `filter_config.json`. Cannot be combined with `--from-source` |
| `--references` | list | `hg19 hg38` | Physical references to install: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`, `hg38_ensembl`. `GRCh37` and `GRCh38` are the files a `hg19_ncbi`/`hg38_ncbi` run also uses — see [Reference Assemblies](../user-guide/reference-assemblies.md) |
| `--from-source` | flag | off | Build every selected reference from its upstream source instead of fetching the published bundle. Slower — downloads and BWA-indexes each selected chromosome FASTA itself — and needs four seed files (`MUC1_motifs_Rev_com.fa`, `code-adVNTR_RUs.fa`, `vntr_db_advntr.zip`, `filter_config.json`). All four are fetched automatically from [`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data)'s `seeds/` directory when not already staged in the output directory — a staged local copy always wins over a download |
| `--release-spec` | path | (none) | `--from-source` only: take every source URL and digest from this file instead of the shipped `install_references_config.json`. Used for release builds |

Whatever `--references` you select, the run also installs the MUC1 motif FASTAs and both
adVNTR databases — those are common assets shared by every assembly, not selectable per
reference, because only one of each exists.

## Examples

Install references with default settings (hg19 + hg38, from the published bundle):

```bash
vntyper install-references -d /path/to/refs
```

Install every physical reference, including NCBI and Ensembl naming:

```bash
vntyper install-references -d /path/to/refs --references hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl
```

Install a single assembly — this still installs the common MUC1/adVNTR assets:

```bash
vntyper install-references -d /path/to/refs --references hg19
```

Rebuild from upstream sources with multiple aligners and 8 threads, instead of the published bundle:

```bash
vntyper install-references -d /path/to/refs --from-source --aligners bwa bwa-mem2 -t 8
```

Point the main configuration at the installed references:

```bash
vntyper --config-path /path/to/config.json install-references -d /path/to/refs
```
