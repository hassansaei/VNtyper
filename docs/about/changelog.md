# Changelog

All notable changes to VNtyper 2 are documented on this page.

## 2.0.10 (Current)

Milestone 4, "CRAM and input robustness": closes #213, #225, #209, #178, #165 and
#161.

**A converted alignment that contains more than one read layout now stops instead of
discarding one layout.** In the golden cohort, 32 of 50 base BAM cases produce paired
FASTQs plus at least one non-empty single FASTQ. VNtyper 2.0.9 silently ignored those
single reads; 2.0.10 names every produced FASTQ and record count and refuses the mixed
layout. Pure single-end BAMs and `--fastq1` without `--fastq2` are now supported through
fastp, BWA and Kestrel.

**BAM and CRAM inputs are preflighted through a run-local alignment view.** Missing
indexes are built beside that view, where htslib can resolve them, without writing into a
read-only patient input tree. The alignment and freshly generated index are opened and
retained for the plan lifetime; post-preflight random-access commands use that exact index
inode, so replacing the public BAI/CRAI name cannot retarget a slice or depth calculation.
Existing BAI, CSI and CRAI spellings are protected by format, and every derived
BAM/FASTQ/log destination is checked before a processing stage starts.

**Reference-dependent CRAM fails before processing unless one candidate can decode the
real `-P` slice shape.** `--reference-fasta`, configured CRAM/BWA references and local
htslib resolution are tried in configured order; coverage receives the same selected
reference. Ambient network refget is disabled by default, `REF_PATH` is restored on every
exit, overlapping in-process CRAM scopes are serialized, and web jobs transport a curated
preflight code/message without exposing worker paths.

**Unmapped CRAM extraction is fail-closed.** `auto` uses indexed `'*'` extraction only
when idxstats proves no placed-unmapped reads and the exact literal-`'*'` flag-4 count
equals the terminal idxstats count; malformed or unequal evidence selects the stream
path, and forcing an unsafe indexed path raises instead of losing reads. On a zero-placed
`7a61` CRAM, whole-file decoding and idxstats each reported 622,690 records while the raw
indexed fetch returned only 2,690; the count proof refused that 620,000-read loss.

**Archives and web delivery are ownership-bound.** Result archives are built and served
through descriptor-anchored paths, reject unsafe source aliases, and are quarantined or
revoked when a later worker status step fails. Cohort inputs are snapshotted through held
descriptors, so replacing a public path cannot change the bytes a child consumes.

Chromosome naming now votes only across configured primary-contig patterns, requires a
strict majority and returns `unknown` for ties or zero classified contigs. Threshold and
pattern configuration is validated rather than accepted partially.

On the 18 golden-cohort BAM cases that both 2.0.9 and 2.0.10 can complete without
discarding reads, three final alternating runs measured medians of 90.92 s and 88.09 s,
respectively. The non-overlapping ranges (90.58–91.12 vs 87.64–88.27 s) record a 3.1%
median improvement under the predeclared rule.

## 2.0.9

Milestone 2, "Web service and cohort integrity": closes #216, #207, #206, #205, #201,
#167 and #162, plus #210 from the CRAM milestone because #162's acceptance criterion
depends on it.

**Cohort pseudonyms change format, and existing pseudonymization tables do not carry
over.** A pseudonym was the prefix plus five hex characters of MD5 -- 20 bits, which
collides with probability ~37.9% at 1,000 samples, and `sample_mapping` is keyed on it,
so a collision silently merged two patients into one row and `sample_categories` counted
them once. It is now 12 hex characters of SHA-256 (~1.8e-9 at 1,000 samples), with the
algorithm and length in `config.json` under `cohort.pseudonym`. Re-run any cohort whose
table you need to match new output.

**Zip cohort inputs report a different sample name.** A zip whose `pipeline_summary.json`
sits at its root -- the layout the web worker produces -- took its reported sample name
and its sort position from `tempfile.mkdtemp`, so both differed on every run. The identity
now comes from the input file the run itself recorded, and the ordering key contains no
temporary component, so two runs of one cohort produce byte-identical CSV, TSV and JSON
exports. Directory inputs keep the order they had.

**Two samples that cannot be told apart are now reported separately or refused, never
merged.** Samples sharing a name are reported under `<input>/<name>`; two inputs that
share a name as well abort with a message naming both. A cohort reached through two
spellings of the same path -- an absolute parent and a relative child -- is one sample
again rather than two.

**A read-only input directory works.** Nothing in the `vntyper` package writes beside the
input any more: the `samtools quickcheck` log goes to the run's output directory, and a
BAM index is resolved in htslib's order and built into the output directory when absent.
On the web service this also makes the worker's input-directory cleanup reachable for the
first time -- every job used to leave one file and one directory behind, permanently, on
the volume every job shares.

**Two XSS fixes in the reports.** The flag tooltip decoded the server's escaping and
reparsed it, so a `Flag` value from a supplied `pipeline_summary.json` could execute; both
templates now build that markup through DOM APIs. The two igv-reports fragments were
interpolated into a `<script>` block behind a non-empty check and are now re-serialised
through `json.dumps` with every `<` escaped.

**`vntyper report` produces the IGV panel again.** It passed no VCF at all, and resolved
the BAM and BED only under `--input-dir`; all three now resolve from one effective run
directory, and `--vcf-file` was added.

No genotype field changed: the golden-cohort gate was run over the full 60-case matrix
against 2.0.8 and every `Confidence`, `Flag`, `POS`, `REF`/`ALT`, `Motif*`, `Depth_Score`
and `VID` is byte-identical.

## 2.0.8

Milestone 2, "Correctness of reported numbers": closes #171, #172, #174, #203 and #212.

**Coverage numbers change for every sample with an uncovered VNTR base.** Historical
reports are not directly comparable with 2.0.8 output. For any sample whose region span
is unchanged, the corrected mean is recoverable in closed form from the old output:

```
mean_corrected == cov_mean_old * (1 - cov_percent_uncovered_old / 100)
```

### Reported-output changes

- **Every VNTR coverage statistic is now computed over the region, not over the covered
  positions** (#171). `samtools depth` was called without `-a`, so it emitted only
  positions with at least one read, and `mean`, `median`, `stdev` and `min` divided by
  the covered-base count. The bias was largest exactly where coverage is patchy: across
  8215 sample rows from seven cohorts, 1585 carried an inflated mean and 61 were reported
  as meeting the 100x QC threshold while falling below it. `min` is now `0` wherever any
  base is uncovered, which it always should have been. `percent_uncovered` keeps its
  value but changes derivation — it counts zero-depth positions rather than subtracting
  the row count, because with `-a` that subtraction is always zero.
- **The quality gate uses both coverage metrics** (#172). `percent_vntr_uncovered` was
  configured with a threshold of 50.0 and compared to nothing; it drove a report icon and
  no decision, so a sample with acceptable mean coverage and half the VNTR uncovered
  passed QC. That is the wrong failure mode — a patchy VNTR is where a frameshift call is
  missed. `quality_metrics_pass` now requires both. The verdict is emitted as an explicit
  `coverage_qc` column in the coverage summary, in `pipeline_summary.json`, in the report,
  and in a new `cohort_stats.{csv,tsv,json}` cohort export.
  `docs/pipeline/reports.md` already documented this threshold as enforced; it was not,
  and now is.
- **Known 4 bp insertion artifacts are excluded from calls** (#174). A row flagged
  `False_Positive_4bp_Insertion` reached `kestrel_result.tsv` and mapped to
  `High_Precision_flagged`, which the report treats as a **positive finding** — a known
  technical artifact presented as a positive MUC1 call. A new `flag_filter_pass` gate
  excludes rows carrying a flag declared in `kestrel_config.json`'s new `artifact_flags`
  list. Advisory flags such as `Low_Depth_Conserved_Motifs` are unaffected and still only
  deprioritise. The excluded row remains in `kestrel_pre_result.tsv` with the gate
  `False`, so the evidence is never destroyed, and emptying `artifact_flags` restores the
  previous behaviour with no code change.

### Fixes with no effect on genotypes

- **The `POS_fasta` rebase that wrote to a discarded column is removed** (#203).
  `motif_correction_and_annotation` subtracted `position_threshold` from `POS` and, because
  `Series.mask` keeps its name, `DataFrame.update` wrote the result back to `POS` — a
  column nothing reads afterwards. The rebase was also *wrong*: `Motif_fasta` is a verbatim
  copy of the VCF `#CHROM`, which names a **120 bp pair record** of
  `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`, so `POS` is already a coordinate
  in that record's space. `POS_fasta` values are unchanged.
- **`output.bed` writes 0-based half-open intervals** (#203, rider). It emitted
  `[POS, POS+1)`, naming the base *after* the variant, so IGV highlighted the wrong
  position on every row. `generate_bed_file` had no test of its output content; it does now.

### Safety

- **The Kestrel stage is never skipped because `output.vcf` exists** (#212). The skip was
  unconditional, undocumented, untested and gated on no flag, and because it `return`ed it
  also skipped the post-processing that produces `kestrel_result.tsv`. Re-running into a
  directory left by an interrupted run — after a crash, an OOM kill or a cancelled job —
  could therefore re-report a stale result, or produce no result file at all, which both
  the report and cohort mode render as a **negative**. A stale VCF is now removed with a
  warning and Kestrel runs. If no k-mer size produces a result, the run raises instead of
  returning silently.
- **A step whose result file is missing is recorded loudly** (#212). `record_step` logged
  nothing and stored an empty `data` list, indistinguishable from a step that legitimately
  found nothing. It now logs an ERROR and sets `result_file_missing`. The key is added only
  when the file is absent, so a normal run's `pipeline_summary.json` is byte-identical.

### Not included

**#171's region harmonisation.** The hg19 span is 1501 bp and the hg38 span is 4501 bp, so
mean coverage is not comparable across builds. Choosing one homologous span and a flank
convention is a domain decision that restates every hg38 number ever reported, and it is
tracked separately so that each gate diff stays attributable to one change.

### Coverage

The floor moves 80 → 81 (branch-inclusive, 81.70% measured) and the target with it.
`coverage_qc.py` and `coverage_stats.py` are both at 100%; patch coverage of this branch
is 100%.

## 2.0.7

Follow-ups to #179: closes #181–#188, #192 and #194–#197, fixes the defects found while
doing so, commits the golden-cohort gate for the first time, and raises the coverage
floor from 70 to 80 on a strictly harder measurement.

**Branch coverage is now on.** `branch = true` means an `if` that is entered but never
taken no longer counts as covered — the failure mode this codebase actually has. The
floor moved with each measurement (70 → 74 → 79 → 80) and was never lowered to admit a
run. Statement-only coverage of the same suite is 80.77%, so the two numbers are not
interchangeable: deleting `branch = true` would *raise* the reported total while covering
strictly less.

### Reported-output changes

Every one of these is a fix to a rule that was already meant to hold. They are listed
separately because each can change what a report says.

- **adVNTR compound states are judged on the signed net indel, not its magnitude**
  (#182). `abs(Insertion_len - Deletion_length)` discarded the direction of the change,
  and a mixed state satisfies both presence guards, so either arm could admit it —
  `I9_2_A_LEN3&D50_2` (Δ = +2, the non-pathogenic frame) was reported through the
  deletion arm. Both arms now test the sign, which makes the pair equivalent to the
  single signed rule Δ ≡ +1 (mod 3) that Kestrel already applied.
- **A compound state's inserted length is the sum over its `LEN` tokens** (#192).
  `I9_2_A_LEN9&I50_2_A_LEN3` is 12. The previous parser coerced any multi-part remainder
  to **zero**, so a pure-insertion compound collapsed to a net change of 0 and was
  dropped in silence.
- **The Depth_Score mid-band is a closed interval** (#184). A score of exactly
  `high_threshold` is `Low_Precision`, not `High_Precision`.
- **A malformed motif ID drops its own row, not the sample** (#179 W8). The check was an
  aggregate `.max()` over the whole column, so one row with two dashes suppressed every
  call in the sample — a Negative report and exit 0 — while a row with no dash rode
  through unnoticed whenever a sibling row had exactly one.
- **`Potential_Duplicate` sorts on `Depth_Score` alone, stably** (#197).
- **The adVNTR "positive flagged" report rule is guarded on `VID`** (#188/D10), matching
  the "positive" rule beside it; it could previously fire on a Negative call.
- **Indels are classified by length difference** (D2), so a multi-base-REF insertion is
  no longer routed into the deletion file.

### Fail-closed, where a silent wrong answer was possible

Each of these used to log and continue. A run that would previously have produced a
quietly wrong report now aborts.

- `cross_match` raises when a `match_logic` rule cannot be evaluated — "no match" was
  indistinguishable from genuine discordance.
- `filter_final_dataframe` raises when a required gate column is missing from a non-empty
  frame; a missing gate is not a permit (#185).
- `confidence_assignment` requires its calibration keys instead of defaulting them. The
  old defaults were not neutral: `high` fell back to 0.4 where the shipped config sets
  0.00515, a factor of 78.
- `filter_vcf` asserts the pinned Kestrel 1.0.1 allele contract, so a jar swap is loud.
- `extract_unmapped_from_offset` raises on a truncated or non-BAI index instead of
  scanning from offset 0.
- `validate_bam_file` honours its documented `ValueError` — the branch that raises it was
  unreachable.
- `motif_processing` refuses `use_uniform_filtering` with an empty GG allowlist, which
  would delete every GG alternate including canonical dupC (#186). Cannot fire on the
  shipped config.

### Web service and CRAM

- **An accepted CRAM is now run as a CRAM** (#188). The endpoint has always accepted
  `.cram`, but the worker hardcoded `--bam`, and the index fallback named a `.bai` that
  is never written beside a CRAM — so the existence check never found the index the
  worker had just built, and cleanup left the real `.crai` on the shared volume.
- **Patient files come off the shared volume first.** The cleanup block used to open with
  two unguarded Redis calls, so a broker that became unreachable as the task exited left
  the alignment and its index behind — and replaced the pipeline's own exception with a
  connection error.
- **CRAM unmapped-read extraction waits for its writer.** The `tee >(...)` process
  substitution let the shell return while samtools was still flushing: measured on a
  600k-read CRAM, `samtools merge` ran against 199,797 of 200,000 reads and accepted the
  short file with a warning.
- Alignment extensions are matched case-insensitively across all three layers, so the
  `SAMPLE.CRAM` the upload allowlist accepts is no longer rejected by the validator.
  FASTQ stays case-sensitive on purpose — fastp reads `reads.FASTQ.GZ` as text and
  reports 0 reads with exit status 0.

### Evidence

Gated by golden-cohort **run 5** (`4fd638a` → `9816f86`), and no pipeline source changed
after the gated commit. The harness itself is committed for the first time, with
cohort-mode and CRAM cases; see `docs/development/golden-cohort-gate.md`. The branch was
additionally put through an adversarial read-only review, which found the signed-frame
merge blocker above — introduced by the #192 fix that made it reachable — before it could
ship.

## 2.0.6

**Unit coverage went from 25.68% to 70.57%, across 1,744 tests.** The number was not the
point: a mutation experiment had shown it overstated protection badly —
`confidence_assignment.py` had 100% line coverage and a 21% mutation score, meaning its
central condition could be inverted with the whole suite still green. That module now
scores 57.9% raw / 84.6% adjusted at the same 100% line coverage.

**No reported genotype or confidence label changes.** That is measured, not asserted: the
golden-cohort gate ran 122 pipeline runs before and after, and the Kestrel variant set,
`Confidence`, `Flag`, adVNTR genotype and `Flag` are identical on every sample and every
assembly. `vntyper/scripts/kestrel_config.json` is byte-identical to its pre-branch state.
See `docs/development/golden-cohort-gate.md`.

### ⚠️ Required before deploying the image

- **Rotate the two Redis passwords.** They are in this repository's public history;
  removing them from source does not invalidate them.
- **Set `REDIS_PASSWORD`.** The web service now refuses to start without it.
- **Deploy the worker before the API.** `run_vntyper_job` gained an `index_path`
  parameter defaulting to `None`; old-API-to-new-worker is safe, new-API-to-old-worker is
  not.

### Breaking changes (web service only; the CLI is unaffected)

- Alias-only cohort access no longer works — `cohort_id` **and** passphrase are both
  required. Cohorts with no stored passphrase hash are unopenable; treating a missing
  hash as "no passphrase" would be the bypass, and no cohort ever had a working hash
  (see below).
- `/job-status/` returns a generic message for a failed job; the detail is logged
  server-side rather than returned.
- `--output-name` is now refused rather than silently ignored.

### Web service hardening

`docker/app` — 1,679 LOC, internet-facing, processing patient-derived genomic data — had
no tests, no lint, no type checking and no coverage. No automated check had ever read it.
It is now formatted, linted, type-checked and 87% covered, and its path is in CI's filter
so those gates fire.

- **The passphrase check never worked.** passlib 1.7.4's bcrypt backend probe is
  incompatible with bcrypt ≥ 4.1 and raised `ValueError` for every input, so no cohort
  could ever have been protected. Now on `bcrypt` directly.
- **Compose set `CELERY_BROKER_URL`**, which Celery prefers over the constructed URL, so
  the worker's real broker was the unauthenticated one. Removed, with a test blocking its
  return. Two *different* hardcoded Redis defaults were also committed, so with the
  variable unset the API and the worker authenticated differently against the same
  instance.
- Upload filenames are constrained by an allowlist before they build a path; uploads and
  request bodies are bounded; job ids are validated as UUIDs; cohort analysis no longer
  deletes its members' archives; `/cohort-download/` cleans up its temporary archive.

### Pipeline fixes

Each with a test written first and observed failing.

- **`vntyper report` was dead** — an uncaught `TypeError`.
- **The `Motif` column was missing from every per-sample report**; the pipeline emits
  `Motif` and the report keyed on `Motifs`, and unmatched columns are silently dropped.
- **Every negative report was styled as positive**, and screening states with no
  configured message fell through to "the screening was negative" — including
  Kestrel-positive-and-flagged.
- **`--extra-modules advntr,shark` silently disabled both modules** (exact-string
  membership, no comma split), so a typo yielded a silent Kestrel-only run.
- **`_construct_ncbi_accession` compared against `"hg19"`** while its only caller passes
  `"GRCh37"`/`"GRCh38"`, so every non-chr1 GRCh37 NCBI name was built from the GRCh38
  table.
- **`online_mode` exited 0 on a failed job**, so a wrapping `subprocess.run(check=True)`
  saw success.
- **`set -o pipefail` added to three multi-stage pipes** — without it a failing
  `samtools sort` upstream of `samtools fastq` still exits 0, genotyping a truncated
  FASTQ as if complete.
- **The `RU == 7` adVNTR flag now fires.** It was added three hours and 45 minutes
  *before* the `RU` column existed, so it had never once evaluated. RU-7 calls move from
  `positive` to `positive flagged`; no row loses a flag and no genotype field changes.
- Jinja2 autoescaping on both the per-sample and cohort reports; `shlex.quote` on
  interpolated paths, regions and sample names.

### Test infrastructure

- **A patch-coverage gate at 80% of changed lines**, which is what makes "touch a file,
  add tests for it" enforceable — demonstrated failing on an untested one-line change.
- **A mutation-testing harness** with the results and the equivalent-mutant
  classifications committed (`docs/development/mutation-testing.md`).
- **A golden-cohort gate** (`docs/development/golden-cohort-gate.md`).
- `tests/builders.py` domain-object builders, `summary_steps.py` for the five step-name
  literals four modules matched by exact string, and `calibration.json` recording the
  provenance of the ten calibrated constants.

### Characterised, not fixed

Seventeen findings needing a domain decision or exceeding a hardening PR's scope were
each pinned by a test and filed as #181–#197. Three worth reading: **#186**
(`motifs_for_alt_gg` ships empty, inert in the active branch but destructive in the
uniform one), **#181** (every (3n+1)-bp deletion discarded), and **#185**
(`filter_final_dataframe` fails open).

Closes #179.

## 2.0.5

**Docker image now runs Python 3.12.** Python 3.9 reached end of life on 2025-10-31.
Every package on the numerical path is unchanged - bwa, samtools, fastp, bcftools,
openjdk (Kestrel), pysam, numpy, pandas and biopython all resolve to the same versions -
and the pipeline reproduces the expected genotype exactly. Only the reporting stack
moved (matplotlib, plotly, igv-reports, jinja2). The package now requires Python 3.10+
and is tested on 3.10-3.13.

**CI/CD pipeline rebuilt.** The Docker image is split into a rarely-rebuilt base
(conda environments, adVNTR, reference genomes) and a per-commit application image, so a
commit now rebuilds in about 3 minutes instead of 40-70. Images are published only after
their tests pass, and the build uses the repository's own source rather than cloning
GitHub, so a pull request now tests the code in that pull request.

- Fixed `pyproject.toml` declaring `numpy<2.0.0` while the conda environment installed
  numpy 2.0.2: `pip install .` inside the image downgraded conda's numpy, so the
  published image ran a different numerical stack than its own recipe declared.
- Fixed the TSV/CSV summary parsers silently truncating rows whose field count did not
  match the header; a malformed row is now logged and skipped without discarding the
  rest of the file.
- Recovered 30 unit tests that carried no `unit` marker and had therefore never run in
  CI, covering the Issue #136/#145 genotype tie-breaking logic.
- Added a fast image smoke tier that verifies every reference path `config.json`
  declares actually exists in the image.
- Added coverage reporting with a ratcheting floor and an 80% target.
- `make check-all` now runs in seconds and works on a fresh clone; it previously
  required a 1.1 GB download and a Docker daemon.

## 2.0.4

- Migrated all logging to per-module loggers (`logging.getLogger(__name__)`); log
  records now carry the emitting module name instead of `root`.
- Replaced deprecated `datetime.utcnow()` calls; emitted timestamps are unchanged.
- Ruff lint rules are now selected explicitly instead of extending ruff's defaults, so
  lint results no longer change when ruff widens its default rule set.
- mypy configuration moved from Makefile flags into `pyproject.toml`.
- Fixed `.gitignore` ignoring the tracked `docs/` directory, which silently hid newly
  added documentation pages from Git.
- Fixed `make test-quiet`, which used a `--log-cli=false` flag pytest does not accept.
- Added `AGENTS.md` and `CLAUDE.md` with repository instructions for coding agents.

## 2.0.3

- Raised Kestrel max align and haplotype states to 40, fixing GDP inflation
  ([#156](https://github.com/hassansaei/VNtyper/issues/156)).
- Handled `pd.NA` and `None` in flagging condition evaluation to prevent `TypeError`s.
- Added a GDP inflation guard test with an anonymized hg38 sample.
- Added `CITATION.cff` with the Zenodo concept DOI.

## 2.0.2

- Removed specific motifs from `exclude_motifs_right` to prevent unwanted motif filtering.
- Clarified depth score conditions in confidence assignment
  ([#154](https://github.com/hassansaei/VNtyper/issues/154)).
- Added a dynamic version variable to the docs via `mkdocs-macros-plugin`.
- Modernized docs deployment to use GitHub Actions Pages.
- Added a CFLAGS workaround for adVNTR compilation errors on GCC 14+.

## 2.0.1

- Disabled duplicate flagging in Kestrel configuration.
- Cleared `motifs_for_alt_gg` array to prevent unintended variant filtering.
- Fixed flagging to occur before variant selection, preventing selection of flagged variants ([#145](https://github.com/hassansaei/VNtyper/issues/145)).
- Fixed confidence assignment to filter sub-threshold variants at the root rather than via override.
- Raised `Low_Depth_Conserved_Motifs` threshold from 0.2 to 0.4.
- Updated variant selection to prioritize depth score before haplotype count.
- Extracted frameshift deduplication into a helper function for clarity.
- Added comprehensive unit tests for the flagging module.

## 2.0.0

A complete rewrite and modernization of the VNtyper 2 pipeline.

- **Modern Python packaging** using `pyproject.toml` (PEP 517/518/621).
- **Refactored Kestrel postprocessing** with configurable thresholds via `kestrel_config.json`.
- **HTML reports** with embedded IGV integration for interactive variant inspection.
- **Cohort analysis** command with built-in pseudonymization for multi-sample studies.
- **Docker multi-stage build** for reproducible, lightweight container images.
- **Comprehensive test suite** with unit and integration tests, including downloadable test data.
- **Multiple reference assemblies** supported: hg19, hg38, GRCh37, and GRCh38.
- **Modular architecture** separating variant parsing, motif processing, scoring, flagging, and confidence assignment into individual modules.
- **Snakemake workflow** for scalable batch processing.

## 1.x

The original VNtyper release as described in Saei et al., *iScience* 2023. This version provided the initial alignment-free genotyping approach for MUC1 VNTR using Kestrel and adVNTR.
