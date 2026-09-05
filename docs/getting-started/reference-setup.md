# Reference Setup

VNtyper 2 requires reference files before running the pipeline. The `install-references` command downloads and prepares everything automatically.

## Where the references come from

By default, `install-references` fetches a published, checksummed reference release from
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data) (currently
[`refs-v1`](https://github.com/berntpopp/vntyper-data/releases/tag/refs-v1)) instead of
downloading and building each reference from six different upstream hosts. The release
publishes six genome bundles (`hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`,
`hg38_ensembl`, representing the physical files the eight `--reference-assembly` labels resolve to; see
[Reference Assemblies](../user-guide/reference-assemblies.md)) plus one `muc1` bundle
carrying the MUC1 motif FASTAs and both adVNTR databases, alongside `SHA256SUMS`, a
`release-manifest.json`, a `BUILD_INFO.json` recording the exact builder commit and tool
versions, and a `verification-report.json` in which all 16 checks passed.

### Trust model

Each asset's SHA-256 digest is committed in `vntyper/scripts/install_references_config.json`.
That committed value, not anything fetched at install time, is the trust anchor.
`install-references` downloads an asset, verifies it against that digest **before**
extracting a single byte, and only unpacks it once the digest matches. The whole install
activates atomically, so a failed or interrupted run never leaves the output directory
half-populated. The release's top-level `SHA256SUMS` file exists for a maintainer
reviewing a draft release by hand: it corroborates the committed digests, but it is never
itself the authority `install-references` trusts, since a release asset could not be tampered
with in a way a loose `SHA256SUMS` file sitting beside it would fail to reflect.

### Building from source instead

`--from-source` rebuilds every requested reference from its original upstream source (UCSC,
NCBI, Ensembl) rather than fetching the bundle. It is slower because it downloads and BWA-indexes
up to six chromosome FASTAs itself, and it needs four seed files that VNtyper itself does
not carry: `MUC1_motifs_Rev_com.fa`, `code-adVNTR_RUs.fa`, `vntr_db_advntr_v2.zip`, and
`filter_config.json`. All four are fetched automatically from
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data)'s `seeds/` directory
(pinned to an immutable commit in `install_references_config.json`, never a mutable branch)
when they are not already present in the output directory. A seed staged there ahead of
time wins over a download, allowing the release build workflow to stage all four first
and run `--from-source` without network access for the seeds themselves. Most users installing
references for a pipeline run do not need this flag: it exists for reproducing or auditing a release,
not as the default install path.

## What Gets Downloaded

The reference installation includes:

| Component | Description |
|-----------|-------------|
| **Chromosome 1 sequences** | UCSC hg19/hg38 (default), NCBI GRCh37/GRCh38, or Ensembl references, each with a pre-built BWA index |
| **MUC1 motif databases** | Pairwise and self-merged MUC1 motif FASTA files with samtools indices |
| **MUC1 region FASTAs** | Per-coordinate-system region FASTAs SHARK filters against (`muc1_region_hg19.fa`, `muc1_region_hg38.fa`) |
| **adVNTR VNTR databases** | Both `hg19_muc1.db` and `hg38_muc1.db`, installed unconditionally (not gated by `--references`) |

## Basic Installation

Download references with default settings (hg19 + hg38, BWA aligner):

```bash
vntyper install-references -d /path/to/references
```

Fetching and verifying the bundle is fast; most wall time is BWA re-indexing, which
only runs if the locally installed `bwa` differs in version from the one the bundle index
was built with. `--skip-indexing` and `-t`/`--threads` therefore only apply to
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
| `--derive-only` | `false` | Rebuild only derived files from installed references without downloading |
| `--release-spec` | *(none)* | `--from-source` only: take every source URL and digest from this file instead of shipped config |

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

## Downloaded files and derived files

Not every file in the reference tree is downloaded. Three are **derived** (built locally from
files that were downloaded), and knowing which is which matters when a tree is incomplete:

| File | How it is produced | From |
| --- | --- | --- |
| `alignment/chr1.<assembly>.fa` (×6) | downloaded | UCSC, NCBI RefSeq, Ensembl |
| `MUC1_motifs_Rev_com.fa`, `code-adVNTR_RUs.fa` | downloaded (seeds) | `berntpopp/vntyper-data`, pinned by commit |
| `filter_config.json` | downloaded (seed): **`--from-source` only**, see below | `berntpopp/vntyper-data`, pinned by commit |
| `vntr_db_advntr_v2/*.db` | downloaded | `berntpopp/vntyper-data` |
| **`muc1_region_hg19.fa`** | **derived** | `samtools faidx chr1.hg19.fa chr1:155158000-155163000` |
| **`muc1_region_hg38.fa`** | **derived** | `samtools faidx chr1.hg38.fa chr1:155184000-155194000` |
| **`All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`** | **derived** | merged from `MUC1_motifs_Rev_com.fa` + `filter_config.json` |

`filter_config.json` is the one seed a **bundle-installed tree does not have**. The bundle
ships the merged motif FASTA pre-built, so nothing on that path needs the filter rules,
and they are not staged beside it. Only `--from-source`, which builds that FASTA itself,
downloads it. This distinction governs `--derive-only` below.

Reference **bytes** are hosted outside this repository: VNtyper seeds in
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data), and genomes at UCSC,
NCBI, and Ensembl. Verification remains local: every source URL, expected
checksum, derivation rule, and exact coordinate boundary is defined in
`vntyper/scripts/install_references_config.json`. Hosting moved; the digests that
verify fetched or derived files did not.

**The two standard paths handle this automatically.** The published bundle ships derived files
pre-built, and `--from-source` builds them at the end of its run.

### Rebuilding just the derived files

`--derive-only` rebuilds derived files from what is already on disk. It **downloads
nothing**:

```bash
vntyper install-references -d /path/to/references --derive-only
```

Below is output from a run against the Docker reference tree (message text only):

```text
Deriving reference files from 6 installed genome(s) in /opt/vntyper/reference: GRCh37, GRCh38, hg19, hg19_ensembl, hg38, hg38_ensembl
Skipping All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa: seed(s) filter_config.json are not in this tree. The published bundle ships this file pre-built and does not stage its seeds beside it, so this is the normal shape of a bundle-installed tree rather than a fault.
  ✓ verified muc1_region_hg19.fa
Derived muc1_region_hg19.fa from chr1.hg19.fa at chr1:155158000-155163000
  ✓ verified muc1_region_hg38.fa
Derived muc1_region_hg38.fa from chr1.hg38.fa at chr1:155184000-155194000
Derived and verified 2 of 3 reference file(s) in /opt/vntyper/reference. Not rebuilt in this tree: All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa.
  ✓ verified All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa
Of those, already present and matching their committed digests: All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa
```

On a tree that contains the seeds (such as one built by `--from-source`), the
motif FASTA is rebuilt as well, concluding with:
`Derived and verified all 3 reference file(s) against their committed digests.`

Each output is verified against its committed checksum, exactly as on the other two paths. A
derived file that does not match is deleted immediately, and the run exits with an error.

Use this option if a derived file is missing or corrupted. Without it, the only alternative is a
full `--from-source` run, which re-downloads and BWA-indexes six chromosome FASTAs to rebuild
three small files. Do not slice regions with `samtools faidx` manually: manual slicing is
unverified, and coordinate errors generate subtly incorrect reference models rather than clean failures.

**A derivation this tree cannot rebuild is skipped, not failed.** A tree
holding only hg19 derives only the hg19 region, and a bundle-installed
tree has no `filter_config.json`. In both cases, the skip is explicitly logged.

**A file that could not be rebuilt is still verified.** Its digest is committed and the file is
small. If such a file is missing rather than outdated, the
command reports `Of those, missing from the tree: ...` and identifies what to install.

If a file is present but does **not** match its digest, it is deleted and the run fails,
protecting downstream runs from corrupted references.

`--derive-only` cannot be combined with `--from-source`; `--from-source` already derives all files.

### Checking a tree is complete

Every path the pipeline reads is declared in `vntyper/config.json` under `reference_data`.
Check a tree directly with Python:

```bash
python -c "
import json, os
paths = json.load(open('vntyper/config.json'))['reference_data']
missing = [(k, v) for k, v in paths.items() if v and not os.path.exists(v)]
print(f'{len(paths) - len(missing)}/{len(paths)} present')
for key, path in missing:
    print(f'  MISSING {key}: {path}')
"
```

If `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa` is missing, Kestrel exits 0
without writing a VCF, and the pipeline refuses to report that as a negative result. Check
the reference tree before debugging pipeline code.

## Output Directory Structure

After installation, the reference directory layout is:

```
references/
  alignment/
    chr1.hg19.fa             # hg19 chromosome 1 sequence (extracted)
    chr1.hg19.fa.amb         # BWA index files (shipped pre-built in bundle)
    chr1.hg19.fa.ann
    chr1.hg19.fa.bwt
    chr1.hg19.fa.pac
    chr1.hg19.fa.sa
    chr1.hg38.fa             # hg38 chromosome 1 sequence (extracted)
    chr1.hg38.fa.*           # BWA index files
  vntr_db_advntr_v2/         # Both adVNTR databases, always installed
    hg19_muc1.db
    hg38_muc1.db
  muc1_region_hg19.fa        # SHARK GRCh37-coordinate region FASTA
  muc1_region_hg19.fa.fai
  muc1_region_hg38.fa        # SHARK GRCh38-coordinate region FASTA
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
**absolute** path into its `reference_data` section, matching the config keys the
pipeline reads (see [Reference Assemblies](../user-guide/reference-assemblies.md)). Writing
absolute paths prevents working-directory alignment failures when installing into non-default locations.

## Storage Requirements

Expect approximately 2 GB of disk space for default installation (hg19 + hg38 references with BWA indices). Adding NCBI or Ensembl references increases storage requirements proportionally.

!!! note "Verification, not MD5"
    `install-references` does not write an `md5_checksums.txt`. The security model relies on
    the SHA-256 digest committed in `install_references_config.json`, validated before any byte
    is extracted. Check `install_references.log` in the output directory for the execution audit trail.
