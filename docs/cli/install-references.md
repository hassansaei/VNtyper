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
| `--references` | list | `hg19 hg38` | Physical references to install: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`, `hg38_ensembl`. `GRCh37` and `GRCh38` are the files a `hg19_ncbi`/`hg38_ncbi` run also uses — see [Reference Assemblies](../user-guide/reference-assemblies.md) |
| `--from-source` | flag | off | Build every selected reference from its upstream source instead of fetching the published bundle. Slower — downloads and BWA-indexes each selected chromosome FASTA itself — and requires the MUC1 seed files (`MUC1_motifs_Rev_com.fa`, `filter_config.json`) already present in the output directory |
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
