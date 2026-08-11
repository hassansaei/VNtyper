# Milestone 5 — References and configuration

**Design document.** Covers GitHub milestone 5 (“4. References and configuration”):
[#163](https://github.com/hassansaei/VNtyper/issues/163),
[#152](https://github.com/hassansaei/VNtyper/issues/152),
[#193](https://github.com/hassansaei/VNtyper/issues/193),
[#217](https://github.com/hassansaei/VNtyper/issues/217),
plus the `shlex.quote` fix reverted on #198 and the disposition of PR
[#164](https://github.com/hassansaei/VNtyper/pull/164).

Written against `main` @ `dc8ec22` (v2.0.11). Every claim below was re-verified at that
commit; the audit comments on the issues were written at `39b4a54` (v2.0.7) and one
material thing has changed since — see F1.

**Revision 3.** This document and its implementation plan were each put through an
adversarial review (Codex `gpt-5.6-sol`, high reasoning effort, read-only against this
worktree). Every finding was verified against the code before being accepted.

*Round 1, on the design* — 4 blockers, 5 majors, 3 minors. The blockers: a fifth, unlisted
BWA reference reader inside the log-file safety guard (§4.2); a key scheme that could never
resolve `hg19_ncbi`/`hg38_ncbi`, because they are aliases of `GRCh37`/`GRCh38` (§4.2); a
bundle with no consumer-side trust anchor (§4.1); and a bootstrap ordering in which the
first bundle could not be built through the shared `--from-source` path (§4.6).

*Round 2, on the plan* — 8 blockers, 13 majors, 4 minors. The ones that changed this
document: the two-tier key scheme would have removed an existing tested capability, so the
resolver is now three-tier (§4.2); `Dockerfile.base` runs the installer without the VNtyper
package, so every new import would have failed the base build (§4.5b); installation must
*merge* rather than replace, or a second `--references` run deletes the first (§4.2); and
the release step is operator-only — `AGENTS.md` forbids an agent from tagging, and
`publish-pypi.yml` has no `push.tags` trigger at all (§9).

*Owner decisions after that review, 2026-08-11:* **two PRs** rather than one pre-merge-pinned
PR (§4.6, D9), and **keep `install_references.py` self-contained** rather than copying four
modules into the Docker refs stage (§4.5b).

---

## 1. What this milestone is actually about

Three of the four issues are symptoms of one disease: **VNtyper has no single source of
truth mapping `(reference kind, assembly) → file path`.** Five independent, mutually
inconsistent mappings exist today, and a fifth namespace is written by
`install-references` that nothing reads. The fourth issue, #217, is about where the
reference *bytes* come from — currently third-party hosts, at image-build time, with no
content pinning.

The milestone therefore has two halves:

- **Resolution** — one resolver, five readers, one writer. Fixes #163 and #152.
- **Provenance** — reference data moves to `berntpopp/vntyper-data`, published as
  versioned, checksummed, verified release bundles. Fixes #217, and dissolves the
  `raw.githubusercontent.com/.../main/...` coupling that PR #164 tripped over.

#193 and the `shlex.quote` revert ride along because they are parked on the same base-image
rebuild.

---

## 2. Verified findings

### F1 — The repo already chose a key scheme, and half the code does not know

Milestone 4 (PR #230) introduced `vntyper/scripts/reference_resolution.py:18`,
`configured_reference_candidates()`. It resolves reference paths as
`{prefix}_{exact_assembly}` first, then falls back to `{prefix}_{ucsc_family_assembly}`,
computing the family through `reference_registry`. It already applies this to
`bwa_reference_*`.

But `vntyper/scripts/cli_handlers.py:264-275` still collapses by hand:

```python
coord_system = get_coordinate_system(args.reference_assembly)
ucsc_map = {"GRCh37": "hg19", "GRCh38": "hg38"}
ucsc_style_ref = ucsc_map.get(coord_system, "hg19")
bwa_key = f"bwa_reference_{ucsc_style_ref}"
```

So for `--reference-assembly hg38_ensembl`, the CRAM preflight looks for
`bwa_reference_hg38_ensembl` and falls back to `bwa_reference_hg38`, while the realignment
step looks *only* at `bwa_reference_hg38`. **The two halves of the same run can disagree
about which FASTA is meant.**

This reframes #163: the fix is to extend a decision the repo already made, not to invent
one. It also sets the shape — exact key first, documented fallback second, registry-driven.

### F2 — adVNTR database selection is a raw dict with a wrong default

`vntyper/scripts/pipeline.py:485-491`:

```python
ref_map = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}
ucsc_style_ref = ref_map.get(reference_assembly, "hg19")
```

Four of the eight registry assemblies — `hg19_ncbi`, `hg38_ncbi`, `hg19_ensembl`,
`hg38_ensembl` — are absent from the dict. Two of them (`hg38_ncbi`, `hg38_ensembl`)
therefore silently select the **wrong coordinate system’s** adVNTR database, with no
warning. Genotype-affecting.

### F3 — `install-references --config-path` writes into a namespace nothing reads

`install_references.py:781-820` builds only `ucsc_*`, `ncbi_*`, `ensembl_*`, `vntyper_*`,
`own_repo_*` keys, and `update_config()` (`:328-365`) is a pure insert. The nine keys
`reference_data` actually holds are never touched. `docker/Dockerfile:50-56` states the
same thing independently.

### F4 — The written paths name files that cannot be used

- `ucsc_hg19` → `<out>/alignment/chr1.hg19.fa.gz`, while indexing runs on
  `target_path.with_suffix("")` (`:436`). The key names a gzip BWA cannot read.
- `vntyper_vntr_db_advntr` → `<out>/vntr_db_advntr.zip`, while the pipeline reads
  `advntr_reference_vntr_hg{19,38}` → `reference/vntr_db_advntr/*.db` (`:809` vs
  `pipeline.py:492-498`). Nothing bridges the two.

### F5 — `reference/md5_checksums.txt` is a build artefact, not a manifest

`write_md5_checksums()` (`:551-568`) writes `output_dir/md5_checksums.txt`. The committed
copy is therefore whatever a past `install-references -d reference/` produced: it lists
two Ensembl `.gz` files and three MUC1 FASTAs, and does **not** list
`muc1_region_hg19.fa`. Any future install into `reference/` overwrites it. The “two-tier
convention with an md5 for every Tier-1 file” described on #152 is real as a *convention*
but is not actually implemented, and appending to this file would not survive.

### F6 — Both SHARK references are pure derivations, and this was measured

```
$ samtools faidx reference/alignment/chr1.hg19.fa chr1:155158000-155163000
c9737129069d4855b433b178ebb21e1c  reference/muc1_region_hg19.fa   (committed)
c9737129069d4855b433b178ebb21e1c  derived                          (identical, 50-char lines)
```

The committed SHARK reference reproduces **byte-for-byte** from UCSC hg19 chr1. The hg38
sibling derives just as cleanly: `chr1:155184000-155194000`, 10,001 bp, header
`>chr1:155184000-155194000`.

Soft-masking differs between the two: hg19 region is 11.9 % lowercase (596/5001), hg38
region is 2.1 % (206/10001). This raised the question of whether SHARK is case-sensitive —
**measured and closed**, see F6b.

### F6b — SHARK is case-insensitive, and #152’s missing measurement now exists

Both were run against the locally installed `shark_env` binary (k=17 default, confidence
0.6) on the seven `tests/data/` FASTQ pairs.

**Case-sensitivity — closed.** Running SHARK against the soft-masked region FASTA and
against an upper-cased copy retains byte-identical read counts (57,698 vs 57,698 for hg19;
62,560 vs 62,560 for hg38). No case-normalisation step is needed in the derivations.

**The hg38 defect, quantified.** Retained R1 reads per sample, hg19 region reference versus
hg38 region reference, on the same input:

| sample | input | hg19 ref | hg38 ref | delta |
|---|---:|---:|---:|---:|
| example_6449 | 82,523 | 57,698 | 62,560 | +8.4 % |
| example_66bf | 19,877 | 17,642 | 18,206 | +3.2 % |
| example_6c28 | 59,440 | 21,425 | 28,864 | **+34.7 %** |
| example_7a61 | 481,761 | 2,090 | 2,459 | +17.7 % |
| example_a5c1 | 21,372 | 15,323 | 16,362 | +6.8 % |
| example_b178 | 16,929 | 12,248 | 12,830 | +4.8 % |
| example_dfc3 | 33,819 | 19,242 | 22,304 | +15.9 % |

Positive on every sample, and this is the *conservative* direction — these FASTQs were
derived from an hg19 5 kb extraction, so genuinely hg38-aligned input should show at least
as large a gap.

**The mechanism, directly.** Canonical 17-mer content of the two regions:

```
hg19 region   5,001 bp   4,434 canonical 17-mers
hg38 region  10,001 bp   7,446 canonical 17-mers
shared                   4,422
hg38-only                3,024   40.6 % of hg38   <- what the hg19 filter cannot see
hg19-only                   12    0.3 % of hg19
```

**A user running `--reference-assembly hg38 --extra-modules shark` today is filtering
against a Bloom filter that is missing 40.6 % of the relevant k-mer content.** This is the
“reopening artefact” the #187 close asked for and nobody had gathered.

It also settles a design question the numbers invite: since the hg38 region is a 99.7 %
superset of the hg19 region, one *could* ship a single hg38 FASTA for both assemblies and
satisfy #187’s “keep one FASTA” literally. **Do not.** That would add 3,024 k-mers to the
hg19 filter and change what the currently-validated hg19 path retains — a genotype-affecting
change to the one path that works. Two references, selected by coordinate system, fixes
hg38 while leaving hg19 byte-identical.

### F7 — The writer has zero test coverage

`update_config`, `process_ucsc_references`, `process_vntyper_references`,
`process_own_repository_references` and `write_md5_checksums` appear nowhere under
`tests/`. Only `execute_aligner_index` is covered, and only as *characterisation of the
unquoted behaviour* in `tests/unit/test_shell_quoting.py`.

### F8 — The Docker image is immune to #163 by accident

`docker/Dockerfile.base:148` copies references to `/opt/vntyper/reference` and `:150` sets
`WORKDIR`, and `docker/Dockerfile:50-64` deliberately skips `install-references`. So the
shipped relative paths in `config.json` resolve. #163 is exclusively a source-install /
custom `--output-dir` bug — which is the path the documentation tells users to take.

**Constraint this imposes:** the in-image reference layout must not change.
`tests/docker/test_image_structure.py:56-57,244,345,359` pin the relative paths, the
existence of every declared reference, and a 200 MB floor on the two BWA references.

### F9 — Batching is forced by the content hash

`conda/**`, `docker/Dockerfile.base`, `docker/requirements-web.txt`,
`vntyper/scripts/install_references.py`, `vntyper/scripts/install_references_config.json`,
`vntyper/dependencies/advntr/**`, `reference/**`, `.dockerignore` — hashed identically in
three places (`docker-base.yml:66`, `docker-build.yml:96`, `docker-build.yml:168`),
mirrored by `Makefile:511-513` (`BASE_INPUTS`) and guarded by
`tests/unit/test_workflow_consistency.py:82-85`. Everything in this milestone touches that
set. One rebuild, one PR, from a maintainer branch — a fork’s `GITHUB_TOKEN` cannot push
to GHCR (AGENTS.md trap 10).

### F10 — PR #164 cannot land as it stands

Fork branch (`ElenaPianfetti:fix/correct_reference_paths`), `mergeable=dirty` against a
`main` that has moved 2.0.4 → 2.0.11. Correct diagnosis and a partly-correct patch, but it
deletes `reference/muc1_region_hg19.fa` and adds a download entry pointing at
`raw.githubusercontent.com/hassansaei/VNtyper/main/reference/muc1_region_hg19.fa` — a URL
that 404s the moment the PR merges. It also writes config key `shark_reference_hg19`,
which SHARK does not read.

### F11 — The Ensembl source URLs are unpinned

`install_references_config.json:72` uses `.../pub/grch37/current/fasta/...` and `:78` uses
`.../pub/current_fasta/...`. Both track “current” and change content on every Ensembl
release. Nothing records which release a given base image was built from. This is #217’s
reproducibility argument in its most concrete form.

### F12 — The CLI and the web service disagree about which assemblies exist

`cli_parser.py:333-338` lets `vntyper online` accept all eight registry labels and
`online_mode.py:86-113` submits the string verbatim, but the server enum
(`docker/app/main.py:58-63`) contains only `hg19`, `hg38`, `GRCh37`, `GRCh38`. So
`vntyper online --reference-assembly hg38_ensembl` subsets the BAM, uploads it, and *then*
gets a 422 reported as a generic submission failure. All four source-qualified labels are
affected. Pre-existing, not caused by this milestone, but squarely inside “references and
configuration” and cheap to fix alongside the registry work.

---

## 3. Decisions

| # | Decision | Rationale |
|---|---|---|
| D1 | **Unify on the milestone-4 key scheme.** One resolver in `reference_registry`, five readers routed through it, and the writer emitting the same namespace. | F1: the scheme already exists and is already load-bearing for CRAM. Extending it removes a live inconsistency instead of adding a fifth mapping. |
| D2 | **Build `berntpopp/vntyper-data`**, following the proven `berntpopp/phentrieve-data` pattern. | #217’s four open questions are all answered by an existing, working precedent in the same hands. |
| D3 | **BWA indexes ship in the bundle**, with the BWA version recorded and verified. | This is the whole point: it converts the 25–120 min base build into a network fetch. Version mismatch degrades to a local re-index, so it can never be silently wrong. |
| D4 | **All reference data moves to `vntyper-data`**, including the small MUC1-scoped FASTAs. `reference/` in VNtyper reduces to provenance scripts. | F6 proves the SHARK references are derivable, so nothing is lost. It removes `reference/**` from the base-image hash set and kills the `raw.githubusercontent.com/.../main/...` coupling permanently. |
| D5 | **The three non-derivable seeds live in `vntyper-data/seeds/`.** | `MUC1_motifs_Rev_com.fa`, `code-adVNTR_RUs.fa` and `vntr_db_advntr.zip` have no upstream and no derivation script (`reference/README.md` §3: non-MUC1 entries stripped by hand, hg38 start found by BLAT). They are data, not software. Diverges from phentrieve-data’s metadata-only rule only because HPO is fully upstream-fetchable and MUC1 is not. |
| D6 | **Bundle is the default; `--from-source` is a supported fallback.** | The `--from-source` path *is* what the bundle build workflow runs, so the derivation logic stays exercised by CI rather than rotting. Air-gapped users keep a route. |
| D7 | **GitHub releases only**, no Zenodo. | Matches phentrieve-data exactly. One publishing path, scriptable, proven. |
| D8 | **Tier-1 checksums move into the release manifest**, not a hand-edited `md5_checksums.txt`. | F5: that file is machine-written and gets overwritten. After D4 it is deleted from the repo entirely. |
| D9 | **Two VNtyper PRs: a builder-only PR, then the bundle is published, then the consumer PR.** The data repo pins the builder at PR-1’s **merge commit**. | The bundle build runs VNtyper’s own `--from-source` path (D6), which PR-1 introduces — so nothing can be built until it has merged. Pinning an unmerged branch SHA would work mechanically but would publish an immutable release built from unreviewed code, and published releases cannot be re-cut. Owner decision 2026-08-11: accept two PRs and two base rebuilds in exchange for no re-publish risk. See §4.6. |
| D10 | **Rework PR #164 on a maintainer branch with `Co-Authored-By: ElenaPianfetti`**, then close it with a pointer to the superseding PR. | F10: it can never go green from a fork, and its direction is right. Credit the diagnosis, land it correctly. |

---

## 4. Architecture

> **`reference/` means two different things, and the distinction matters throughout.**
> *In git*, after D4, `reference/` holds only provenance scripts — no FASTAs, no `.fai`, no
> `.zip`, no `md5_checksums.txt`. *At runtime*, `reference/` (or whatever `--output-dir`
> names) is the directory `install-references` populates from the bundle. The shipped
> relative paths in `config.json` continue to name `reference/…` because that is the
> runtime default and the in-image layout; they are not references to tracked files.
>
> Placeholders below marked `<…>` and `"sha256": "..."` are values determined during
> implementation — the reusable workflow’s commit SHA, the upstream digests, and the
> `bcrypt` version the published base resolves. They are not unresolved design questions.

### 4.1 Part A — `berntpopp/vntyper-data`

Mirrors `berntpopp/phentrieve-data` structurally.

```
berntpopp/vntyper-data
├── LICENSE                              MIT (the metadata and scripts)
├── ATTRIBUTION.md                       upstream terms: UCSC, NCBI, Ensembl, adVNTR
├── README.md                            release lifecycle, install-and-verify, policy
├── .gitignore                           build outputs never enter git
├── seeds/                               the three non-derivable artefacts  (D5)
│   ├── MUC1_motifs_Rev_com.fa  (+ .fai)
│   ├── code-adVNTR_RUs.fa      (+ .fai)
│   ├── vntr_db_advntr.zip
│   └── filter_config.json
├── releases/
│   └── refs-v1.json                     the committed release spec
└── .github/workflows/
    └── release-data.yml                 ~15 lines; dispatches a pinned reusable workflow
```

**The build workflow lives in the software repo**, exactly as phentrieve does it:
`hassansaei/VNtyper/.github/workflows/build-reference-bundles.yml`, called via
`workflow_call`, referenced from the data repo at a **full commit SHA**:

```yaml
# berntpopp/vntyper-data/.github/workflows/release-data.yml
on:
  workflow_dispatch:
    inputs:
      release_tag: { type: choice, options: [refs-v1] }   # only committed specs
permissions: { contents: write }
jobs:
  build:
    uses: hassansaei/VNtyper/.github/workflows/build-reference-bundles.yml@<full-sha>
    with:
      release_tag:       ${{ inputs.release_tag }}
      release_spec_path: releases/${{ inputs.release_tag }}.json
      source_commit:     <full-sha>
```

**`releases/refs-v1.json`** pins every input by digest, and — critically — pins Ensembl to
an explicit release rather than `current` (F11):

```jsonc
{
  "release_tag": "refs-v1",
  "source_commit": "<VNtyper full sha>",
  "bwa_version": "0.7.18-r1243",
  "samtools_version": "1.21",
  "sources": {
    "ucsc_hg19":     { "url": "...goldenPath/hg19/chromosomes/chr1.fa.gz", "sha256": "..." },
    "ucsc_hg38":     { "url": "...goldenPath/hg38/chromosomes/chr1.fa.gz", "sha256": "..." },
    "ncbi_GRCh37":   { "url": "...GCF_000001405.25_GRCh37.p13...",         "sha256": "..." },
    "ncbi_GRCh38":   { "url": "...GCF_000001405.40_GRCh38.p14...",         "sha256": "..." },
    "ensembl_hg19":  { "url": ".../pub/release-115/fasta/...",             "sha256": "..." },
    "ensembl_hg38":  { "url": ".../pub/release-115/fasta/...",             "sha256": "..." }
  },
  "seeds": {
    "MUC1_motifs_Rev_com.fa": { "sha256": "..." },
    "code-adVNTR_RUs.fa":     { "sha256": "..." },
    "vntr_db_advntr.zip":     { "sha256": "..." },
    "filter_config.json":     { "sha256": "..." }
  },
  "derivations": [
    { "output": "muc1_region_hg19.fa",
      "command": "samtools faidx chr1.hg19.fa chr1:155158000-155163000",
      "expected_sha256": "..." },
    { "output": "muc1_region_hg38.fa",
      "command": "samtools faidx chr1.hg38.fa chr1:155184000-155194000",
      "expected_sha256": "..." },
    { "output": "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
      "command": "python generate_vntr_reference.py",
      "expected_sha256": "..." }
  ]
}
```

`expected_sha256` on each derivation is what makes the build **self-verifying**: F6 already
established that `muc1_region_hg19.fa` reproduces byte-for-byte, so the workflow asserts it
rather than trusting it. A silent change in a UCSC chromosome file turns the build red
instead of publishing different sequence under the same name.

**Release assets**, per bundle version:

```
vntyper-references-refs-v1-ucsc-hg19.tar.gz      chr1.hg19.fa + .amb/.ann/.bwt/.pac/.sa + .fai
vntyper-references-refs-v1-ucsc-hg38.tar.gz
vntyper-references-refs-v1-ncbi-GRCh37.tar.gz
vntyper-references-refs-v1-ncbi-GRCh38.tar.gz
vntyper-references-refs-v1-ensembl-hg19.tar.gz
vntyper-references-refs-v1-ensembl-hg38.tar.gz
vntyper-references-refs-v1-muc1.tar.gz           the MUC1-scoped FASTAs + adVNTR DBs (~3 MB)
SHA256SUMS
release-manifest.json                            path, size, sha256, source URL, derivation
verification-report.json                         what was checked and what it produced
BUILD_INFO.json                                  bwa/samtools versions, runner, timestamps
```

Each genome asset is a few hundred MB — comfortably under GitHub’s 2 GB per-asset limit,
and split per assembly so `--references hg38` fetches one genome asset rather than all six.

**`vntyper-references-<tag>-muc1.tar.gz` is not selectable — it is an implicit dependency of
every installation.** The MUC1-scoped FASTAs and both adVNTR databases are needed by every
run regardless of assembly. Today this is half-true by accident and half-broken:
`own_repository_references` is processed unconditionally (`install_references.py:771`), but
`vntyper_references` *is* filtered by `--references` (`:670`), so
`vntyper install-references --references hg19` — which is exactly what
`.github/workflows/scheduled-tests.yml:103-110` runs — installs the MUC1 FASTAs and **no
adVNTR database at all**. Making `muc1` implicit fixes that pre-existing inconsistency in
the same change, and a test must assert that `--references hg19` alone yields a complete
reference tree.

#### The trust anchor must live in VNtyper, not next to the assets

Co-hosting the assets and their `SHA256SUMS` in the same release means a single actor can
replace both, and the consumer would accept it. Pinning only `repository` + `release_tag`
is therefore not a content pin.

**Each asset’s SHA-256 is committed in `vntyper/scripts/install_references_config.json`.**
The release’s own `SHA256SUMS` may corroborate, but it cannot be the root of trust. This
also closes a staleness hole in the Docker path: because that file is a base-image
content-hash input, changing any reference byte now *requires* a commit that changes the
base tag. Without the committed digests, a replaced asset under an unchanged tag would leave
`base-status` finding the old base published (`docker-base.yml:63-98`) and source installs
silently disagreeing with container installs under the same nominal reference version.

**Publication policy** (stated in the data repo’s README, following phentrieve-data):
draft release only; a maintainer publishes after reviewing `verification-report.json` and
checking `SHA256SUMS` against the remote assets. Releases are immutable once published.

**Licensing** (#217 open question 4): `ATTRIBUTION.md` records the upstream terms for UCSC,
NCBI, Ensembl and the adVNTR database, and each release spec records the exact source URL
and digest used. MIT covers the repo’s own scripts and metadata, not the sequence.

### 4.2 Part B — reference resolution (#163)

#### The key must be derived from the *physical identity*, not from the label

The obvious design — key on the assembly label, falling back to the UCSC family — **cannot
work**, and finding this is what the adversarial review earned. `ASSEMBLY_METADATA`
(`reference_registry.py:119-144`) contains **eight labels for six physical files**:
`GRCh37` and `hg19_ncbi` are two names for the same NCBI file, as are `GRCh38` and
`hg38_ncbi`. `install_references_config.json:55-67` only ever produces a `GRCh37` and a
`GRCh38` entry, so a label-keyed writer never emits `bwa_reference_hg38_ncbi`, and a
`--reference-assembly hg38_ncbi` run misses its exact key and silently falls back to the
**UCSC** FASTA — aligning against `chr1` when the user asked for `NC_000001.11`.

The fix is to collapse label → physical identity first. Physical identity is exactly
`(coordinate_system, reference_source)`, and the six resulting ids are *already* the entry
names in `install_references_config.json`:

| label | coordinate system | source | physical id | `bwa_reference_*` keys, in order |
|---|---|---|---|---|
| `hg19` | GRCh37 | ucsc | `hg19` | `_hg19` |
| `hg38` | GRCh38 | ucsc | `hg38` | `_hg38` |
| `GRCh37` | GRCh37 | ncbi | `GRCh37` | `_GRCh37`, `_hg19` |
| `hg19_ncbi` | GRCh37 | ncbi | `GRCh37` | `_hg19_ncbi`, `_GRCh37`, `_hg19` |
| `GRCh38` | GRCh38 | ncbi | `GRCh38` | `_GRCh38`, `_hg38` |
| `hg38_ncbi` | GRCh38 | ncbi | `GRCh38` | `_hg38_ncbi`, `_GRCh38`, `_hg38` |
| `hg19_ensembl` | GRCh37 | ensembl | `hg19_ensembl` | `_hg19_ensembl`, `_hg19` |
| `hg38_ensembl` | GRCh38 | ensembl | `hg38_ensembl` | `_hg38_ensembl`, `_hg38` |

The writer emits only the **physical** key, because that is the file it installed. The
label tier exists so a hand-written replacement config can still specialise or disable one
accepted label, which the current suite already requires. The family tier is the
source-degrading compatibility fallback of §4.2.

**New in `vntyper/scripts/reference_registry.py`:**

```python
REFERENCE_KINDS = {
    "bwa":    {"prefix": "bwa_reference",         "keyed_by": "physical"},
    "cram":   {"prefix": "cram_reference",        "keyed_by": "physical"},
    "advntr": {"prefix": "advntr_reference_vntr", "keyed_by": "coordinate_system"},
    "shark":  {"prefix": "muc1_region_fasta",     "keyed_by": "coordinate_system"},
}

def physical_reference_id(assembly: str) -> str: ...   # label -> one of the six ids
def ucsc_family(assembly: str) -> str: ...             # -> "hg19" | "hg38"
def reference_keys(kind: str, assembly: str) -> tuple[str, ...]: ...
```

The `keyed_by` distinction is biological, not incidental:

- **`bwa` / `cram` are physical-identity-keyed**, with a **label key ahead of the physical
  key**. Contig naming differs between sources (`chr1` vs `1` vs `NC_000001.11`), so
  `hg38_ensembl` genuinely wants a different FASTA from `hg38`. Up to three keys are
  returned: label, physical, UCSC family.

  The label tier is not decoration — it is an existing, tested capability.
  `tests/unit/test_reference_resolution.py`'s
  `test_label_specific_reference_keys_override_the_family_fallback_including_null` pins
  that a replacement config may set `bwa_reference_hg19_ncbi` to specialise one accepted
  label, or `cram_reference_hg19_ncbi: null` to disable it. A purely physical scheme would
  route `hg19_ncbi` through `*_GRCh37` and silently ignore both. Dropping that would be a
  regression disguised as a simplification.

  So for `hg19_ncbi`: `("bwa_reference_hg19_ncbi", "bwa_reference_GRCh37", "bwa_reference_hg19")`.
  For `hg19` all three collapse to one key. Duplicates are removed, order preserved.
- **`advntr` / `shark` are coordinate-system-keyed.** Only two adVNTR databases exist
  (`hg19_muc1.db`, `hg38_muc1.db`) and only two MUC1 regions exist (GRCh37 and GRCh38
  coordinates). Source naming is irrelevant to both. One key, and because the UCSC family
  names *are* `hg19`/`hg38`, the existing key names are unchanged.

#### Lookup semantics: membership, never truthiness

`reference_resolution.py:50` already resolves with `reference_data[exact_key] if exact_key
in reference_data else reference_data.get(family_key)` — **membership**, so an explicit
`null` is authoritative and disables the fallback. `tests/unit/test_reference_resolution.py:68-84`
pins that. Every resolver added here must use the same rule, or a config that deliberately
nulls `cram_reference_hg19_ncbi` would silently be re-enabled by a family fallback.

#### The family fallback is *source-degrading*, and must say so

The earlier claim that this is “fully backwards compatible” was too strong. It is
value-compatible — the nine shipped keys keep working — but for any non-UCSC label under a
config that lacks the exact key, the run **silently uses UCSC sequence**. Downstream
detection hides this rather than catching it: `pipeline.py:337-345` aligns FASTQ against
whatever was selected, and `chromosome_utils.py:288-319` then detects the convention from
the BAM *it just produced*, so the run succeeds and looks correct.

Three mitigations, all required:

1. **Ship all six physical keys** in `vntyper/config.json` as relative paths, so a complete
   installation never takes the fallback at all.
2. **Warn on every fallback**, naming the requested assembly, the key that missed, the key
   used, and the resulting reference source.
3. **Record requested assembly, effective key, path and source in the run summary**, so the
   substitution is visible in the report rather than only in a log.

**Readers routed through it — there are five, not four:**

| Site | Today | After |
|---|---|---|
| `cli_handlers.py:264-275` | hand-collapsed `ucsc_map`, silent `None` on miss | `reference_keys("bwa", assembly)`, **fail closed** with a named error |
| `cli_logging_safety.py:43-51` `_selected_bwa_reference` | **a fifth, independent hand-collapse** | same resolver — see the security note below |
| `pipeline.py:485-501` | raw dict, wrong default for 4 labels (F2) | `reference_keys("advntr", assembly)` |
| `shark_filtering.py:76` | one hardcoded key | `reference_keys("shark", assembly)` via the layered lookup of §4.3 |
| `reference_resolution.py:38-51` | computes the family inline | calls `reference_keys()`; duplication removed |

**`cli_logging_safety` is a safety boundary, and missing it is exploitable.**
`_selected_bwa_reference()` re-implements the collapse independently, and its result feeds
the check that refuses to let `--log-file` point at an operator input. It runs *before*
`setup_logging()` (`cli.py:147`), whose `FileHandler` opens the target in append mode
(`utils.py:95-100`). So with `bwa_reference_hg38_ensembl` configured and
`--reference-assembly hg38_ensembl --log-file <that Ensembl FASTA>`, the guard inspects the
**UCSC** file, sees no conflict, and logging then appends text into the Ensembl reference.
The same holds for its `.amb/.ann/.bwt/.pac/.sa` sidecars. This reader is not optional
cleanup — it is part of the fix.

Failing closed on an unresolvable BWA reference is both what the Sourcery review on #164
asked for and the direction milestone 4 established (PR #230, “make CRAM and alignment
inputs fail closed”). Today a missing key yields `bwa_reference = None` and the run dies
later at `pipeline.py:154-156` with a less obvious message.

**Writer.** `install_references_config.json` entries declare `assembly` and `kind`; the
key itself is *derived* from the registry, never hand-written, so the two sides cannot
drift:

```jsonc
"ucsc_references": {
  "hg38": {                                          // the physical id, per the table above
    "kind": "bwa",
    "asset": "vntyper-references-{tag}-ucsc-hg38.tar.gz",
    "asset_sha256": "...",                           // the trust anchor, committed here
    "installed_path": "alignment/chr1.hg38.fa",      // extracted, not .gz  (fixes F4)
    "source_url": "...", "source_sha256": "..."      // --from-source path only
  }
}
```

The entry key *is* the physical id, so `update_config()` derives
`bwa_reference_{physical_id}` from the registry rather than from a hand-written
`config_key` — the writer and the readers cannot drift because neither owns the name.
The written value is the **extracted** file, closing F3 and F4 together.

#### Installation must be transactional, and extraction must be constrained

The current installer is not safe to extend as-is: `download_file` skips any existing
destination without revalidating it (`:50-72`), archives are unpacked with an unrestricted
`tar.extractall()` (`:423-429`), and `update_config` overwrites the config in place
(`:353-365`). Introducing bundle extraction without fixing this converts a network blip
into a half-populated reference tree that the next run treats as complete.

Required contract:

0. **Merge, do not replace.** Seed the staging directory from the existing tree, so
   `--references hg38` after `--references hg19` does not delete hg19, and installing into
   `reference/` does not delete the retained `README.md` and `pseudonymize*`.
1. Download and extract into a **fresh staging directory on the same filesystem** as the
   target.
2. Verify the asset against `asset_sha256` **before** extraction.
3. Reject any archive member with an absolute path, a `..` component, a device node, or a
   link escaping the root. Python's `tarfile` `filter="data"` is the baseline; the check is
   explicit and tested, not implied by the runtime version.
4. Verify every extracted file against `release-manifest.json` **after** extraction.
5. **Atomically activate** the completed tree, then replace `config.json` last, via
   write-temp-and-rename.
6. On any failure, remove the staging directory and leave the previous tree and config
   untouched. The previous tree is deleted only *after* activation succeeds; if it cannot
   be restored it is preserved under a named path and reported, never discarded.
7. **No custom exception classes** (`AGENTS.md:140`): `logger.error(message)` then
   `ValueError` for invalid contents, `RuntimeError` for download or install failure.

Tested against: interruption mid-extract, path traversal, a corrupt sidecar, a pre-existing
partially-populated destination, and rollback.

### 4.3 Part C — SHARK assembly selection (#152)

`vntyper/modules/shark/shark_config.json` becomes assembly-keyed, as specified on #152:

```json
{"shark_settings": {
    "muc1_region_fasta_hg19": "reference/muc1_region_hg19.fa",
    "muc1_region_fasta_hg38": "reference/muc1_region_hg38.fa"}}
```

**But `shark_config.json` alone is not enough, and this is a trap worth naming.** It is a
*separate* config file from `config.json`, so `install-references --config-path
vntyper/config.json` can never update it. Shipping the SHARK paths only there would
reintroduce #163’s exact bug for SHARK: after an install into a custom `--output-dir`, the
region FASTA path would still point at a CWD-relative `reference/…` that does not exist.

`run_shark_filter()` therefore resolves through a **layered lookup**, most authoritative
first:

1. `main_config["reference_data"][key]` — what `install-references` writes, so a custom
   `--output-dir` is honoured;
2. `shark_config["shark_settings"][key]` — the shipped default, in the assembly-keyed shape
   #152 specified;
3. `shark_config["shark_settings"]["muc1_region_fasta"]` — the legacy flat key, **taken only
   when the config is structurally legacy**, i.e. it carries the flat key and *no*
   `muc1_region_fasta_*` key at all;
4. otherwise `ValueError` naming the assembly and every key tried.

**Membership, not truthiness, at every layer** — the same rule as `reference_resolution.py:50`.
A key that is present and `null` is an authoritative “disabled”, not a miss. Without this,
a site that points hg19 SHARK at `/mnt/refs` and deliberately nulls hg38 would silently
fall through to a shipped CWD-relative `reference/muc1_region_hg38.fa` and either die from
an unrelated working directory or, worse, mix references from two different installations.

The layer-3 guard matters for the same reason: once a config has declared the keyed schema,
falling back to the flat key would let an *incomplete* managed installation masquerade as a
legacy one.

`main_config` is already a parameter of `run_shark_filter` (`shark_filtering.py:35`), so no
new plumbing is required. `update_config()` writes `muc1_region_fasta_hg{19,38}` into
`reference_data` alongside the other canonical keys, which is what makes layer 1 work.

Note the AGENTS.md trap-1 interaction: `shark_config` is read into a module global at
import (`shark_filtering.py:26-27`). The lookup must use the `config` argument, not the
global — and the tests must patch accordingly.

The warning added by `39fe0eb` at `shark_filtering.py:64-70` is deleted along with the
docstring paragraph at `:50-54` that promises the parameter is ignored.

`muc1_region_hg38.fa` is **produced by the bundle build**, not committed — it is a
derivation with a pinned expected digest (§4.1). This is strictly better than the Tier-1
plan on #152: no binary blob in git, and the provenance is executable rather than prose.

### 4.4 Part D — the base-rebuild batch (#193)

1. `docker/requirements-web.txt`: add a pinned `bcrypt==<version the published base
   resolves>`, drop `passlib[bcrypt]==1.7.4`. The version must be read off the **published
   base image**, not off a dev machine, so the pin does not change runtime behaviour.
2. **In the same commit**, delete `"bcrypt": "passlib"` from `UNDECLARED_IMPORT_ALLOWANCES`
   (`tests/unit/test_version_consistency.py:317`) and the explanatory block at `:304-313`.
   `test_undeclared_import_allowances_are_all_still_needed` fires on “module now declared
   directly”, so splitting these across commits leaves the tree red in between.
3. Reconcile `pyproject.toml`’s `web` and `dev` extras — both currently carry
   `passlib[bcrypt]==1.7.4`, and `test_pyproject_pins_the_web_versions_the_image_installs`
   plus `test_web_extra_declares_nothing_the_image_does_not_install` compare against
   `requirements-web.txt`, which is the declared authority.
4. Restore the `shlex.quote` fix in `execute_aligner_index` (`install_references.py:239-291`)
   via `command_builders.quote_path`, and **invert** the three characterisation cases in
   `tests/unit/test_shell_quoting.py` — the module docstring already promises they will be
   inverted at this rebuild.

### 4.5 Part E — Docker

`docker/Dockerfile.base`, `refs` stage:

- delete `COPY reference/ /opt/refs/` (there is nothing left to copy after D4);
- `install-references` fetches and unpacks the pinned bundle instead of downloading six
  genomes and BWA-indexing them;
- keep the existing existence guards and extend them to `muc1_region_hg38.fa`;
- **the layout under `/opt/vntyper/reference` is unchanged**, so `config.json`’s relative
  paths, `WORKDIR`, and `tests/docker/test_image_structure.py` are all untouched.

Content-hash bookkeeping: `reference/**` leaves the `hashFiles()` set at all three call
sites — one in `docker-base.yml`, two in `docker-build.yml`, i.e. **two files, three
sites** — plus `Makefile:511-513` (`BASE_INPUTS`) and
`tests/unit/test_workflow_consistency.py:82-85`. The bundle **pin** lives in
`install_references_config.json`, which is already in the set — so bumping the reference
set still triggers a base rebuild, automatically and for the right reason.

### 4.5b Keep `install_references.py` self-contained

`docker/Dockerfile.base:85-99` copies **only** `install_references.py` and its JSON into
`/opt/ir` and runs the file directly, without installing the package. That works because
the file imports **nothing** from `vntyper` today — verified, zero `from vntyper...` lines.
That is not an accident, and it is worth preserving.

Anything the stage needs must also join the base-image content-hash set, or the base could
be built against a stale copy of a module the installer imports. So every module added
here is a module whose every future edit costs a 25–120 min base rebuild. Measured churn
over the last year:

| module | commits/yr | verdict |
|---|---:|---|
| `command_builders.py` | 11 | **do not import** — and it pulls in `samtools_command_fragments` too. Inline the three-line `quote_path` as a private helper instead. |
| `utils.py` | 8 | **do not import** — it imports `command_builders`, so it drags the same chain in. Use `subprocess` directly. |
| `reference_registry.py` | 3 | **import lazily**, inside `canonical_reference_keys`. Docker deliberately passes no `--config-path` (`Dockerfile:50-56`), so that path never executes there. Function-local registry imports are already the house style — `cli_handlers.py:262`, `pipeline.py`. |
| `reference_bundle.py` | new | **copy it.** It is the bundle-fetch path the stage actually runs. Self-contained: `hashlib`, `tarfile`, `shutil`, `tempfile` only. |

So the stage copies **two** Python files, and the hash set grows by **one** new, stable
file — not five churning ones:

```dockerfile
COPY vntyper/__init__.py                            /opt/ir/vntyper/__init__.py
COPY vntyper/scripts/__init__.py                    /opt/ir/vntyper/scripts/__init__.py
COPY vntyper/scripts/install_references.py          /opt/ir/vntyper/scripts/
COPY vntyper/scripts/install_references_config.json /opt/ir/vntyper/scripts/
COPY vntyper/scripts/reference_bundle.py            /opt/ir/vntyper/scripts/
ENV PYTHONPATH=/opt/ir
RUN conda run -n vntyper python -m vntyper.scripts.install_references --output-dir /opt/refs ...
```

Assert the property rather than trusting it: a unit test that parses
`install_references.py` and `reference_bundle.py` and fails if either grows a
module-scope `from vntyper...` import naming anything outside the copied set.

### 4.6 The bootstrap ordering

There is a genuine chicken-and-egg here. D6 says the bundle build workflow runs VNtyper's
own `--from-source` path, so the derivation logic cannot rot — but `--from-source` does not
exist yet. A data repo pinned at any commit predating that work has no builder to call, and
the only way out would be a second, independent derivation implementation in CI: exactly
the drift D6 exists to prevent.

**Resolution (owner decision, 2026-08-11): two VNtyper PRs.**

1. **PR-1, builder only.** `reference_bundle.py`, the derivations, `install_from_source`,
   `--from-source`, `--release-spec`, and the reusable
   `build-reference-bundles.yml`. **All tracked reference data stays in place.** Merge.
   → base rebuild #1.
2. **Publish `refs-v1`.** Create `berntpopp/vntyper-data`, materialise the seeds from
   PR-1's **merge commit** with `git show`, pin `release-data.yml` at that same merge SHA,
   dispatch, verify the draft, publish.
3. **PR-2, consumers.** The registry resolver and all five readers, the canonical writer,
   SHARK selection, the bundle pin and its digests, deletion of the tracked data, the
   `Dockerfile.base` switch, #193, the server enum, docs, and the golden gate.
   → base rebuild #2.

The alternative — pinning PR-1's *pre-merge* branch SHA to collapse this into one PR — was
considered and rejected. It works mechanically (a non-fork branch SHA is a valid `uses:`
target, and VNtyper merges with merge commits so the SHA survives), but it publishes an
**immutable** release built from unreviewed code. Published releases are never re-cut; a
reviewer-requested builder change would mean `refs-v2`, a full rebuild and reverify, and
new digests in every pin. Two rebuilds are the cheaper certainty.

**Never start the deletion step while the pinned tag is in doubt.** PR-2's deletion task
opens by proving the published bundle already carries every file about to be removed.

**The reusable workflow needs two checkouts.** The seeds and the release spec live in
`berntpopp/vntyper-data`; the builder lives in `hassansaei/VNtyper`. `workflow_call` does
not check out the caller's repository automatically. The workflow must check out both
explicitly, log both resolved SHAs, and record them in `verification-report.json`.

---

## 5. Data flow

```
berntpopp/vntyper-data                    hassansaei/VNtyper
──────────────────────                    ──────────────────
seeds/ + releases/refs-v1.json
        │
        │ workflow_dispatch
        ▼
release-data.yml ────────────────────────▶ build-reference-bundles.yml @<full-sha>
                                                    │
                                           download sources, verify sha256
                                           derive muc1_region_hg{19,38}.fa   ─┐
                                           derive All_Pairwise_...filtered.fa │ assert
                                           bwa index each genome              │ expected
                                           tar + SHA256SUMS + manifests      ─┘ sha256
                                                    │
        ◀───────────────────────────────────────────┘  draft release
   maintainer verifies → publishes refs-v1
        │
        │   install_references_config.json pins repository + release_tag
        ▼
   vntyper install-references -d <dir>
        │  fetch asset → verify SHA256SUMS → check BUILD_INFO.bwa_version
        │  (mismatch → re-index locally, logged, never silent)
        ▼
   <dir>/alignment/…, <dir>/vntr_db_advntr/…, <dir>/muc1_region_hg{19,38}.fa
        │
        │  --config-path → update_config() writes canonical keys (registry-derived)
        ▼
   config.json reference_data
        │
        ├── reference_keys("bwa",    assembly) → cli_handlers   → realignment
        ├── reference_keys("cram",   assembly) → reference_resolution → CRAM preflight
        ├── reference_keys("advntr", assembly) → pipeline       → adVNTR DB
        └── reference_keys("shark",  assembly) → shark_filtering→ MUC1 region FASTA
```

---

## 6. Error handling

| Condition | Behaviour |
|---|---|
| Bundle asset unreachable | Fail with the pinned repo/tag/asset named, and point at `--from-source`. Never fall back silently. |
| `SHA256SUMS` mismatch | Fail closed, name the file and both digests. A wrong reference is worse than no reference. |
| `BUILD_INFO.bwa_version` ≠ local BWA | **Warn and re-index locally.** Degraded, not fatal — but never accept an index built by a different BWA without saying so. |
| Derivation digest ≠ `expected_sha256` (bundle build) | Fail the workflow. Upstream sequence changed under a stable URL; that needs a human. |
| `reference_keys("bwa", …)` resolves nothing | Fail closed at `cli_handlers`, naming the assembly and both keys tried. Today: silent `None`, obscure death at `pipeline.py:154-156`. |
| `reference_keys("advntr", …)` resolves nothing | Existing `ValueError` at `pipeline.py:503-505` retained. |
| SHARK key missing for the requested assembly | Walk the four layers of §4.3; if all miss, `ValueError` naming the assembly and every key tried. |
| Unknown `--reference-assembly` | `normalize_assembly_name` already raises with the supported list. `cli_handlers.py:266-268` and `region_utils.py:155-161` currently downgrade this to a warning + GRCh37 default — **tighten to fail closed**, matching milestone 4’s direction. |
| Family fallback taken | Warn, naming requested assembly, missed key, used key, and effective reference source; record all four in the run summary. Never silent. |
| Any key present but `null` | Authoritative “disabled”. Do not fall back. Matches `reference_resolution.py:50`. |
| Bundle extraction fails part-way | Staging directory removed; previous tree and `config.json` untouched. A retry must not treat a half-written tree as complete. |
| Archive member escapes the root | Reject the member and fail the install, naming it. |

**Note on the tightening.** `--reference-assembly` on the command line is already
constrained by argparse `choices`, so no valid CLI invocation changes. What *does* change is
the **config-supplied default**: `cli_handlers.py:194-195` reads
`default_values.reference_assembly` unvalidated, so a replacement config carrying e.g.
`"b37"` runs today as GRCh37-with-warnings and will now stop. That is a deliberate breaking
migration and must be released as one — validate replacement configs at load time and print
the supported list. It does **not** surface as an HTTP 500: web inputs are enum-validated,
and worker CLI failures are recorded as failed jobs (`docker/app/tasks.py:313-342`).

---

## 7. Testing strategy

**Unit — the resolver (new).** Table-driven over all 8 labels × 4 kinds, asserting the exact
key tuple — including that `GRCh38` and `hg38_ncbi` resolve to the *same* key, and `GRCh37`
and `hg19_ncbi` likewise. Regression barrier for F1 and F2.

**End-to-end key agreement, not just tuples.** Assert the writer’s actual output is
resolvable by every reader: run the install path, then for all 8 labels × 5 readers check
that a key is found and that no run silently takes the family fallback. A test that only
compares expected key *tuples* would have passed the label-keyed design that could never
resolve `hg38_ncbi`.

**Unit — the writer (F7: currently 0 %).**
- `update_config` writes canonical keys, not the `ucsc_*`/`own_repo_*` namespace;
- written paths name **extracted** files, never `.gz` or `.zip` (F4);
- bundle fetch: happy path, `asset_sha256` mismatch → fail closed, manifest mismatch → fail
  closed, BWA-version mismatch → re-index with a warning, asset 404 → named error;
- transactional behaviour: interruption mid-extract, path traversal, corrupt sidecar,
  pre-existing partial destination, rollback leaves prior tree and config intact;
- **`--references hg19` alone installs a complete tree** — MUC1 FASTAs *and* both adVNTR
  databases (the pre-existing `:670` filtering bug);
- `--from-source` produces a tree identical to the bundle path.

**Unit — the readers.**
- BWA key for all 8 labels, plus fail-closed on an unresolvable key;
- **`cli_logging_safety` refuses `--log-file` pointing at the exact reference for every
  physical key and every `.amb/.ann/.bwt/.pac/.sa` sidecar**, including hard-link and
  symlink aliases (the blocker in §4.2);
- adVNTR key for all 8 labels — direct regression for F2;
- present-but-`null` at every layer is authoritative, for every kind;
- SHARK selection per assembly; **invert** `TestReferenceAssemblyIsAccceptedAndIgnored` and
  `TestNonHg19AssemblyLogsAWarning` (`tests/unit/test_shark_filtering.py:179,236`) from
  ignore-and-warn into selection tests, and update `SHARK_CONFIG` at `:31` and the command
  assertions at `:68`, `:151`, `:214`;
- **SHARK case-insensitivity guard.** F6b holds for the pinned SHARK 1.2.0
  (`shark 1.2.0 h077b44d_5`; its `to_int` table maps `a/c/g/t` and `A/C/G/T` identically). A
  future base rebuild could change that silently, because cohort reads need not contain a
  discriminating lowercase-only k-mer. Add a deterministic test on a synthetic reference
  whose only possible match crosses lowercase bases, asserting identical output against an
  upper-cased copy.

**Unit — CLI/server contract.** Parity-test `docker/app/main.py`’s `ReferenceAssembly` enum
against `list_assemblies()`, and add an `online`-submission contract test, so F12 cannot
recur.

**Unit — the batch.** Invert the three `execute_aligner_index` cases in
`tests/unit/test_shell_quoting.py`; delete the `bcrypt` allowance in the same commit as the
requirements change.

**Integration — the actual #163 repro.** `install-references --output-dir <tmp>
--config-path <tmp cfg>`, then a pipeline run against that config, asserting it resolves
every reference from `<tmp>`. This is the scenario the issue reports and nothing tests it
today.

**Docker.** Extend `test_every_declared_reference_exists` and the `MINIMUM_SIZES` table
(`tests/docker/test_image_structure.py:56-57,345`) to cover `muc1_region_hg38.fa`. The
relative-path assertion at `:244` must keep passing unchanged — that is the check that the
in-image layout did not move.

**Verification is four commands, not one.** `make check-all` is
`format-check lint type-check-all test-unit` (`Makefile:409`) — it runs neither
`lint-actions` nor anything Docker, and this milestone changes both. Run `make check-all`,
`make ci-local`, `make ci-local-docker` and `pytest tests/docker/`.

**Golden-cohort gate.** Genotype-affecting: adVNTR database selection changes for
`hg38_ncbi` and `hg38_ensembl`, and SHARK’s reference changes for every GRCh38 run. Hand-run
the 60-case matrix (~14 min on 32 cores; generate the CRAM fixtures rather than passing
`--allow-matrix-drift`).

Two prerequisites the gate procedure does not currently state, and must:
- both trees need a **populated** `reference/`, which after D4 means installing and
  verifying the published bundle *first* — `docs/development/golden-cohort-gate.md:897-905`
  assumes tracked data is simply present;
- `scripts/golden_cohort/matrix.py:109` omits `hg19_ncbi` and `hg38_ncbi`, the two labels
  most at risk from the physical-identity change. Adding them to `ASSEMBLIES` achieves
  nothing — that tuple has no executable use; cases are derived from fixture filenames at
  `:230-283` with counts pinned at `:172-190`. Generate explicit alias cases and update the
  pinned counts, the matrix tests, and the gate documentation.

---

## 8. Risks and open technical checks

| | Risk | Mitigation |
|---|---|---|
| ~~R1~~ | ~~SHARK may be case-sensitive.~~ | **Closed by measurement (F6b).** Soft-masked and upper-cased references retain byte-identical read counts. No case-normalisation step in the derivations; expected digests are those of the plain `samtools faidx` output. |
| **R2** | The `bcrypt` pin could change runtime behaviour if taken from a dev machine (locally 5.0.0). | Read the version off the **published base image** (`docker run … pip show bcrypt`) and pin that. |
| **R3** | Ensembl `current` URLs (F11) mean the first bundle cannot reproduce whatever the current base image contains. | Accepted and deliberate: `refs-v1` pins an explicit Ensembl release and records its digest. The discontinuity is recorded in the release spec, which is the point of #217. |
| **R4** | Bootstrapping: the seeds exist only in VNtyper, are copied to `vntyper-data`, then deleted from VNtyper. A mistake here is data loss. | Copy and publish `refs-v1` **first**; delete from VNtyper only in the consolidated PR, only after the release is published and its digests verified against the files being deleted. Git history retains them regardless. |
| **R5** | The base build gains a dependency on GitHub releases. | Net reduction: it currently depends on UCSC *and* NCBI *and* Ensembl *and* `raw.githubusercontent.com` staying up and byte-stable. One pinned, checksummed host replaces four unpinned ones. |
| **R6** | `--from-source` rots if only the bundle path is exercised. | It *is* the bundle build workflow’s code path (D6), so every reference release exercises it end to end. |
| **R7** | The golden-cohort gate is not in CI and must be hand-run. | Explicit step in the plan, with the result recorded on the PR. |
| **R8** | The consolidated PR is large and touches Docker, CI, config and five call sites. | Ordered, individually-green commits; the resolver and its tests land before any reader is switched, so each reader change has a passing barrier already in place. |
| **R9** | `refs-v1` is built from unmerged code (§4.6), so a late review change to the builder invalidates the published bundle. | Build only once the builder half has settled; verify the published bundle against every `expected_sha256` before Pass 2. Fallback is two PRs and a second base build. |
| **R10** | The family fallback silently degrades source for any non-UCSC label whose exact key is absent, and downstream naming detection *hides* it because it reads the BAM the run just produced. | Ship all six physical keys in `config.json`, warn on every fallback, record it in the run summary, and test that a complete installation and the Docker image never take it. |
| **R11** | A future base rebuild could upgrade SHARK to a case-sensitive version, and ordinary cohort reads would not necessarily catch it. | Record SHARK 1.2.0 (`h077b44d_5`) as the verified implementation; add the synthetic lowercase-discriminating regression test in §7. |
| **R12** | Fail-closed on an unknown config-supplied assembly is a genuine breaking change for anyone using a replacement `config.json`. | Release it as a breaking migration: validate replacement configs at load time, print the supported list, and test the config-default path specifically. Not a web 500 — enum-validated there. |

---

## 9. Deliverables

Per §4.6 this is **one VNtyper PR**, opened early and completed in two passes, with the data
repo built against its branch SHA in between.

**Pass 1 — on the PR branch, all tracked reference data still present**
1. `install_references.py`: the source builder and `--from-source`, plus the transactional
   staging/verify/activate contract of §4.2, with tests. Derivations for
   `muc1_region_hg{19,38}.fa` and `All_Pairwise_…_filtered.fa`.
2. `.github/workflows/build-reference-bundles.yml` (`workflow_call`), checking out **both**
   repositories explicitly and recording both SHAs.
3. Create `berntpopp/vntyper-data`: MIT `LICENSE`, `ATTRIBUTION.md`, `README.md`,
   `.gitignore`, `seeds/` copied from this branch, `releases/refs-v1.json`,
   `.github/workflows/release-data.yml` pinned at the PR branch head SHA.
4. Dispatch; verify the draft (`verification-report.json`, `SHA256SUMS`, and each
   derivation against its `expected_sha256`); publish `refs-v1`.

**Pass 2 — same PR branch, now consuming the published bundle**
5. `reference_registry`: `REFERENCE_KINDS`, `physical_reference_id`, `ucsc_family`,
   `reference_keys` + tests.
6. Route all **five** readers through it; fail closed on unresolvable keys; warn and record
   in the run summary on every family fallback. Ship all six physical keys in
   `vntyper/config.json`. (#163, F1, F2)
7. Finish `install_references_config.json`: bundle `repository` + `release_tag` + per-asset
   `asset_sha256` (the trust anchor), physical-id entry names, `kind`, `installed_path`, and
   the `source_url`/`source_sha256` pair for `--from-source`. Make `muc1` an implicit
   dependency of every install. `update_config` writes registry-derived canonical keys
   pointing at extracted files. (#163, #217, F3, F4)
8. `shark_config.json` assembly-keyed; layered selection in `shark_filtering.py` with
   `reference_data` as layer 1 so `--output-dir` is honoured; `update_config` writes
   `muc1_region_fasta_hg{19,38}`; invert the two test classes; drop the #187 warning and
   docstring. (#152)
9. Delete `reference/` data files and `md5_checksums.txt`; keep the provenance scripts.
   Update `.gitignore`. (D4, D8, F5)
10. `Dockerfile.base`: drop the `COPY`, fetch the bundle, extend the guards. Update
   `hashFiles()` at all three call sites (two workflow files), `Makefile` `BASE_INPUTS`, and
   `test_workflow_consistency.py`.
11. `bcrypt` pin + allowance deletion + `pyproject.toml` extras; restore `shlex.quote` and
    invert its three tests. (#193)
12. Route `cli_logging_safety._selected_bwa_reference` (`:43-51`) through the resolver — the
    safety blocker in §4.2. Add `hg19_ncbi`/`hg38_ncbi` to
    `scripts/golden_cohort/matrix.py:108-109`. Add all eight labels to
    `docker/app/main.py:58-63` and parity-test the enum. (#163, F12)
13. Tests per §7. Docs and executable contracts — the migration touches more than prose:
    - `docs/cli/install-references.md`, `docs/getting-started/reference-setup.md`,
      `docs/user-guide/reference-assemblies.md`, `docs/user-guide/configuration.md`,
      `docs/pipeline/optional-modules.md`, `README.md`, `reference/README.md`;
    - `.dockerignore:89-100` — its comments still declare tracked reference files as image
      inputs;
    - `pyproject.toml:166-206` — Ruff excludes `reference/` *because* it is a base-hash
      input; that reason expires;
    - `SPEC.md:23-29` — its proof command executes against a file scheduled for deletion;
    - `AGENTS.md` trap 10 — the base-input list changes;
    - `docs/development/golden-cohort-gate.md` — add the bundle-install prerequisite;
    - `reference/generate_vntr_reference.py:73-79` hardcodes seed filenames that will no
      longer sit beside it. Either give it explicit input paths or move it into
      `vntyper-data` next to `seeds/`. Leaving it in place, unrunnable, is not an option.
14. `make check-all`; hand-run the golden-cohort gate; merge; bump the data repo’s workflow
    pin from the branch SHA to the merge commit (§4.6).

    **The release itself is operator-only.** `AGENTS.md` "Never": *never create, move, or
    push a release tag as an agent*, and `.github/workflows/publish-pypi.yml` has **no**
    `push.tags` trigger — production starts only from an authenticated
    `repository_dispatch` of type `vntyper_release` naming an **already existing** reviewed
    tag. So: prepare the version bump in `vntyper/version.py`, `CITATION.cff` and
    `docs/about/changelog.md` (all three must agree — trap 12), merge, wait for the ten
    exact-full-SHA checks to go green, then **stop and hand over**. The operator creates the
    tag and sends the dispatch. `gh workflow run publish-pypi.yml -f tag=vX.Y.Z` inspects
    without writing.
15. Close #163, #152, #193, #217; close PR #164 with a pointer and
    `Co-Authored-By: ElenaPianfetti` on the relevant commits.

---

## 10. Out of scope

- Migrating `tests/data/` off Zenodo (record 17479292) — unrelated, and it works.
- Adding aligners beyond BWA to the bundle. `bwa-mem2` and `minimap2` remain
  `enabled: false` in `install_references_config.json`; the asset layout leaves room but
  ships nothing for them.
- Changing the MUC1 region coordinates themselves. `COORDINATE_SYSTEMS`
  (`reference_registry.py:54-65`) is unchanged.
- Milestone 6 issues (#173, #175, #202, #189, #33, #20, #215) — parked, needs-decision.
