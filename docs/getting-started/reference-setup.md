# Reference Setup

VNtyper 2 requires reference files before running the pipeline. The `install-references` command downloads and prepares everything automatically.

## Where the references come from

By default, `install-references` fetches a published, checksummed reference release from
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data) — currently
[`refs-v1`](https://github.com/berntpopp/vntyper-data/releases/tag/refs-v1) — instead of
downloading and building each reference from six different upstream hosts. The release
publishes six genome bundles (`hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`,
`hg38_ensembl` — the physical files the eight `--reference-assembly` labels resolve to; see
[Reference Assemblies](../user-guide/reference-assemblies.md)) plus one `muc1` bundle
carrying the MUC1 motif FASTAs and both adVNTR databases, alongside `SHA256SUMS`, a
`release-manifest.json`, a `BUILD_INFO.json` recording the exact builder commit and tool
versions, and a `verification-report.json` in which all 16 checks passed.

### Trust model

Each asset's SHA-256 digest is committed in `vntyper/scripts/install_references_config.json`
— that committed value, not anything fetched at install time, is the trust anchor.
`install-references` downloads an asset, verifies it against that digest **before**
extracting a single byte, and only unpacks it once the digest matches; the whole install
activates atomically, so a failed or interrupted run never leaves the output directory
half-populated. The release's own top-level `SHA256SUMS` file exists for a maintainer
reviewing a draft release by hand — it corroborates the committed digests, but it is never
itself the thing `install-references` trusts, since a release asset could not be tampered
with in a way a loose `SHA256SUMS` sitting beside it would also fail to reflect.

### Building from source instead

`--from-source` rebuilds every requested reference from its original upstream source (UCSC,
NCBI, Ensembl) rather than fetching the bundle. It is slower — it downloads and BWA-indexes
up to six chromosome FASTAs itself — and it needs four seed files that VNtyper itself does
not carry: `MUC1_motifs_Rev_com.fa`, `code-adVNTR_RUs.fa`, `vntr_db_advntr.zip` and
`filter_config.json`. All four are fetched automatically from
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data)'s `seeds/` directory
(pinned to an immutable commit in `install_references_config.json`, never a mutable branch)
when they are not already present in the output directory — a seed staged there ahead of
time always wins over a download, which is what lets the workflow that builds the published
release stage all four first and run `--from-source` against them without any network access
for the seeds themselves. Most users installing references for a pipeline run do not need
this flag — it exists for reproducing or auditing a release, not as the default install path.

## What Gets Downloaded

The reference installation includes:

| Component | Description |
|-----------|-------------|
| **Chromosome 1 sequences** | UCSC hg19/hg38 (default), NCBI GRCh37/GRCh38, or Ensembl references, each with a pre-built BWA index |
| **MUC1 motif databases** | Pairwise and self-merged MUC1 motif FASTA files with samtools indices |
| **MUC1 region FASTAs** | Per-coordinate-system region FASTAs SHARK filters against (`muc1_region_hg19.fa`, `muc1_region_hg38.fa`) |
| **adVNTR VNTR databases** | Both `hg19_muc1.db` and `hg38_muc1.db`, installed unconditionally — not gated by `--references` |

## Basic Installation

Download references with default settings (hg19 + hg38, BWA aligner):

```bash
vntyper install-references -d /path/to/references
```

Fetching and verifying the bundle is fast; most of the wall time is BWA re-indexing, which
only runs if the locally installed `bwa` differs in version from the one the bundle's index
was built with. `--skip-indexing` and `-t`/`--threads` therefore only matter for
`--from-source`.

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `-d`, `--output-dir` | *(required)* | Directory where references will be installed |
| `--references` | `hg19 hg38` | Physical references to install: `hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`, `hg38_ensembl` |
| `--aligners` | `bwa` | `--from-source` only: aligners to build indices for (e.g., `bwa bwa-mem2 minimap2`) |
| `--skip-indexing` | `false` | `--from-source` only: skip BWA index building |
| `-t`, `--threads` | `4` | `--from-source` only: threads for indexing |
| `--from-source` | `false` | Build every reference from its upstream source instead of fetching the published bundle |
| `--release-spec` | *(none)* | `--from-source` only: take every source URL and digest from this file instead of the shipped config |

### Examples

Download only hg38 references:

```bash
vntyper install-references -d ./references --references hg38
```

Download all supported physical references, including NCBI and Ensembl naming:

```bash
vntyper install-references -d ./references --references hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl
```

Rebuild from upstream sources instead of the published bundle (slower; fetches the four MUC1/adVNTR seed files from `berntpopp/vntyper-data` unless already staged):

```bash
vntyper install-references -d ./references --from-source
```

## Output Directory Structure

After installation, the reference directory looks like this:

```
references/
  alignment/
    chr1.hg19.fa             # hg19 chromosome 1 sequence (extracted)
    chr1.hg19.fa.amb         # BWA index files (shipped pre-built in the bundle)
    chr1.hg19.fa.ann
    chr1.hg19.fa.bwt
    chr1.hg19.fa.pac
    chr1.hg19.fa.sa
    chr1.hg38.fa             # hg38 chromosome 1 sequence (extracted)
    chr1.hg38.fa.*           # BWA index files
  vntr_db_advntr/            # both adVNTR databases, always installed
    hg19_muc1.db
    hg38_muc1.db
  muc1_region_hg19.fa        # SHARK's GRCh37-coordinate region FASTA
  muc1_region_hg19.fa.fai
  muc1_region_hg38.fa        # SHARK's GRCh38-coordinate region FASTA
  muc1_region_hg38.fa.fai
  All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa
  All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa.fai
  MUC1_motifs_Rev_com.fa
  MUC1_motifs_Rev_com.fa.fai
  code-adVNTR_RUs.fa
  code-adVNTR_RUs.fa.fai
  install_references.log
```

If `--config-path` names a `config.json`, `install-references` writes each installed file's
**absolute** path into its `reference_data` section, using exactly the config keys the
pipeline reads (see [Reference Assemblies](../user-guide/reference-assemblies.md)). Writing
absolute rather than relative paths is deliberate: installing into a directory other than
the default with `--output-dir` used to leave `config.json` holding the shipped relative
paths, and the run died as soon as the process's working directory did not happen to line
up with them (issue #163). Two of the old keys also pointed at files nobody could open — the
downloaded `.gz` archive rather than the extracted FASTA, and the `.zip` archive rather than
the two extracted `.db` files — which is fixed the same way, by writing the key that names
the file the pipeline actually opens.

## Storage Requirements

Expect approximately 2 GB of disk space for the default installation (hg19 + hg38 references with BWA indices). Adding NCBI or Ensembl references increases this proportionally.

!!! note "Verification, not MD5"
    `install-references` no longer writes an `md5_checksums.txt` — the trust model is the
    SHA-256 digest committed in `install_references_config.json`, checked before any byte
    of a downloaded asset is extracted (see [Trust model](#trust-model) above). Check
    `install_references.log` in the output directory for the full record of what was
    fetched, verified and installed.
