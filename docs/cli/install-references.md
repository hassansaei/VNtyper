# install-references

Download and configure reference files required by VNtyper 2.

By default, this command downloads verified pre-built bundles from [`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data). See [Reference Setup](../getting-started/reference-setup.md) for verification details and `--from-source` requirements.

## Synopsis

```
vntyper [global-options] install-references
    -d <dir>
    [--skip-indexing]
    [-t <int>]
    [--aligners <aligner> ...]
    [--references <reference> ...]
    [--derive-only]
    [--from-source]
    [--release-spec <path>]
```

## Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-d, --output-dir` | path | (required) | Destination path for reference installation |
| `--skip-indexing` | flag | off | Skip BWA indexing (applies to `--from-source`; pre-built bundles re-index only on BWA version mismatches) |
| `-t, --threads` | int | `4` | Indexing thread count (used during `--from-source` operations) |
| `--aligners` | list | `bwa` | Aligners for which indices are built (`bwa`, `bwa-mem2`, `minimap2`). Only used with `--from-source` |
| `--derive-only` | flag | off | Rebuild derived files (`muc1_region_hg19.fa`, `muc1_region_hg38.fa`, `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`) from local files without downloading. Missing sources skip derivation and verify existing files in place. Cannot combine with `--from-source` |
| `--references` | list | `hg19 hg38` | Physical reference genomes to install: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`, `hg38_ensembl` |
| `--from-source` | flag | off | Download and index chromosome FASTAs directly from primary sources instead of the release bundle. Requires seeds (`MUC1_motifs_Rev_com.fa`, `code-adVNTR_RUs.fa`, `vntr_db_advntr_v2.zip`, `filter_config.json`) |
| `--release-spec` | path | None | File providing explicit source URLs and digests for release builds (requires `--from-source`) |

Every run installs shared assets (MUC1 motif FASTAs and adVNTR databases) regardless of selected reference assemblies.

## Examples

Install default reference set (hg19 and hg38 bundle):

```bash
vntyper install-references -d /path/to/refs
```

Install all reference genomes:

```bash
vntyper install-references -d /path/to/refs --references hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl
```

Install a single assembly:

```bash
vntyper install-references -d /path/to/refs --references hg19
```

Rebuild from upstream sources with multiple aligner indices:

```bash
vntyper install-references -d /path/to/refs --from-source --aligners bwa bwa-mem2 -t 8
```

Update custom configuration to point to installed references:

```bash
vntyper --config-path /path/to/config.json install-references -d /path/to/refs
```
