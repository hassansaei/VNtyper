# Changelog

All notable changes to VNtyper 2 are documented on this page.

## Unreleased

No unreleased changes.

## 2.0.16 (Current)

**Pin what we build on, and let no attempt inherit another's leftovers.**

### Added

- **`vntyper install-references --derive-only`** rebuilds the derived reference files from
  genomes and seeds already on disk, downloading nothing. Recovering a tree missing one of
  them previously meant a full `--from-source` run — re-downloading and BWA-indexing six
  chromosome FASTAs to rebuild three small ones — so in practice it was done by hand, with
  `samtools faidx` retyped from the config, and by hand is where an unverified reference file
  comes from. Every output is checked against its committed digest. A derivation this tree
  cannot rebuild is skipped rather than failed, and whatever is already at its path is
  verified instead, so all three are checked either way.

### Fixed

- **A discarded Kestrel attempt no longer leaves its SAM behind for a later one to convert
  (#255).** Every k-mer size writes to the same `output.vcf` *and* the same `output.sam`. A
  discarded attempt removed the VCF but left the SAM, and a later successful attempt converts
  whatever occupies that path into `output.bam` — the alignment the report's IGV track shows.
  The genotype comes from the VCF and is unaffected, which is why it was easy to miss.
- **`convert_sam_to_bam_and_index` checks samtools' exit status instead of inferring it
  (#255).** Success was inferred from `os.path.exists` on the outputs, which conflates
  "samtools succeeded" with "a file is present". A failed `view` still leaves a truncated BAM,
  so the SAM was deleted and the truncated BAM kept. Both results are checked now, and the SAM
  survives a failure because it is the only copy of the data.
- **A reference entry pointing at a 404 is gone (#253).** It was unreachable — the installer
  drops any `raw_files` entry that is also a derivation output — so the issue's claim that
  `--from-source` is broken by it is false and is corrected on the issue. A 404 was the
  mildest way it could have failed: a URL resolving to a stale but live file would have
  installed silently wrong reference data. Tests pin the rule rather than the entry.
- **The unauthenticated download window is 3 days rather than 7 (#189 — a mitigation).**
  `/download/{job_id}/` takes no credential, so the retention window *is* the exposure window:
  the link is the capability. This reduces exposure duration, not exposure. The completion
  email now names the window and says the link carries no password.
- **`.env.example`'s retention knobs are now actually configurable.** `MAX_RESULT_AGE_DAYS`
  and `COHORT_RETENTION_DAYS` were presented as operator settings, but compose passed a
  literal for the first and never passed the second, so uncommenting either configured
  nothing. Default behaviour is unchanged; the documentation becomes true.

### Changed

- **adVNTR is pinned to a commit, not a branch (#254).** `advntr_genotyping.py` cites exact
  upstream line numbers as evidence for why adVNTR runs at `-t 1`, and a branch can move under
  those citations. `GIT_COMMIT` now names `05fd98a4`, and the installer checks it out and
  aborts if it is absent rather than silently building the branch tip.
- **`install_advntr.sh` resolves its config beside the script**, not in the working directory,
  and resolves symlinks to get there — the image runs it from elsewhere, so a bare relative
  name was never sourced and the pin had no effect on the image it exists to pin.
  `-c/--config` now replaces the shipped config rather than layering over it.

## 2.0.15

**Stop claiming more than we know.** Five defects in which VNtyper stated something it
had not established.

### Fixed

- **A Kestrel VCF that cannot be parsed no longer produces a `Negative` (#223).** A run
  that exited 0 having lost its `#CHROM` header discarded every indel record it carried
  and reported a confident negative genotype, with no error logged at all. The zero-byte
  case is the mild one; the severe one converts a positive into a negative. Each
  configured k-mer size is now tried, and a run where none produces a usable VCF fails
  loudly. A valid header with no records is still a legitimate negative and is unchanged.
- **`summarise_coverage` refuses a depth table longer than its region (#224).** It
  previously reported percentages above 100 -- four uncovered positions in a two-base
  region gave `percent_uncovered: 200.0`. Region coordinates are also validated *before*
  `samtools depth` runs, so a reversed or oversized region fails immediately instead of
  writing an unbounded depth file; a reversed region used to reach the report as a
  negative length beside a passing coverage verdict.
- **The adVNTR thread count and output format come from configuration alone (#247).**
  Both code fallbacks contradicted the shipped `advntr_config.json`, so dropping either
  key changed the emitted command silently.

### Changed

- The golden-cohort harness option `--advntr-threads` is now `--advntr-case-threads`,
  and the `side.json` key with it (#247). The value is unchanged. The old name read as
  the thread count adVNTR received, which was never true: adVNTR runs at `-t 1` from its
  own configuration and never sees the CLI value.
- `docs/pipeline/scoring-and-confidence.md` describes what a confidence tier says about
  the evidence rather than what a reader should do with it (#250). VNtyper is research
  use only, which that page previously contradicted five times.
- `--pseudonymize-samples` and the cohort documentation state what the pseudonym does
  **not** protect against (#227): the digest is unsalted and unkeyed, so a shared report
  must be treated as identifying. Documentation only; the digest is unchanged.
- Planning documents are no longer tracked. `docs/` is now strictly the published site.

## 2.0.14

**CRAM runs fail under Docker with `Too many levels of symbolic links` (#238).**

Every Docker user whose CRAM reference has no `.fai` beside it is affected on 2.0.10
through 2.0.13, including `ghcr.io/hassansaei/vntyper:latest`.

CRAM preflight retained the reference index htslib generates by opening it and then
replacing that same pathname with a symlink to `/proc/<pid>/fd/<n>` — a descriptor whose
own recorded path *is* the pathname just replaced. The entry pointed at itself:

```
reference.fa.fai -> /proc/22/fd/7
/proc/22/fd/7    -> …/reference.fa.fai (deleted)
```

Linux tolerates that because it jumps straight to the descriptor's directory entry; a
container runtime that resolves procfd links by re-walking the link text closes the cycle
and every consumer that opens the index gets `ELOOP`. It is the only binding whose source
and destination are the same path, which is why a reference *with* a shipped `.fai`
decoded fine in the same run. The entry is now retained under its own name rather than
replaced, and both binding classes additionally prove an installed run-local view can be
opened through its published pathname before handing it to samtools.

**Input and output must be two different host directories.**

The same report exposed a second defect. Three containment guards compared pathnames only
and never `os.path.samefile`, so one host directory mounted at two container paths —
`-v .:/opt/vntyper/input -v .:/opt/vntyper/output` — defeated all three: the run wrote its
output tree inside the directory holding the patient alignment, operator BEDs and
references inside the output root went unprotected, and the log-path guard, which runs
before the log directory is created, left a directory in the patient tree even when a
later guard refused the run. All three now compare by inode and name both pathnames in
the rejection.

This is a bypass, not a new policy: the identical layout with a *single* mount was already
rejected. Runs using the documented separate-directory layout are unaffected. `README.md`
and the Docker user guide now state the requirement and show both forms.

## 2.0.13

**One reference resolver, and reference data from verified bundles (#217, #163, #152, #193).**

Three genotype-affecting bugs are fixed. `--reference-assembly hg38_ncbi` and `hg38_ensembl`
silently loaded the **GRCh37** adVNTR database, because an incomplete map defaulted every
unrecognised label to hg19. SHARK filtered against the hg19 MUC1 region regardless of
`--reference-assembly` (#152): 40.6% of the hg38 region's canonical 17-mers are absent from the
hg19 region, and the correct reference retains 3.2-34.7% more reads across seven cohort samples.
And after `install-references --output-dir`, `config.json` kept its shipped relative paths, so the
run died (#163) - the installer wrote keys nothing reads while the seven keys the pipeline does
read were never updated.

Eight assembly labels now resolve through one registry to six physical genome files, plus two
adVNTR databases and two MUC1 region FASTAs keyed by coordinate system. All five readers and the
one writer go through it, so the writer and the readers cannot disagree again.

Reference data moved out of this repository into versioned, checksummed release bundles published
from `berntpopp/vntyper-data`. `install-references` downloads each asset, verifies it against the
digest committed here, extracts it with a hardened extractor and activates it atomically; a
per-file provenance ledger records what was verified, and only files with a record are written
into `config.json`. `--from-source` still rebuilds everything from upstream.

**A fresh clone has no reference data** and must run `vntyper install-references`. The Docker image
installs the bundle at build time and now carries all six physical references, so every assembly
label works out of the box; it is correspondingly larger. A configured-but-missing reference now
fails immediately, naming the assembly, the key, the missing path and the command that fixes it.

## 2.0.12

**Valid mixed alignment conversions are now routed losslessly (#233).** When R1 and R2
contain equal non-zero record counts, every non-empty R1, R2, `other`, and singleton
FASTQ is passed exactly once, in that order, to one Kestrel 1.0.1 sample. The files are
not concatenated or recompressed. Unequal or one-sided mate outputs still fail before
Kestrel as invalid. This restores BAM, CRAM, and adVNTR runs rejected by 2.0.10 and
2.0.11 while retaining genuine singleton and unpaired reads instead of silently
discarding them.

## 2.0.11

When both `--keep-intermediates` and `--delete-intermediates` are supplied,
`--delete-intermediates` now wins as the CLI help has always documented. Earlier
versions kept the intermediate BAM in this conflicting form.

**The integration test matrix covers more outcome and input states.** Successful cases now
assert both present and absent intermediate/archive artifacts and run an independent paired
b178 FASTQ case without SHARK (#71). Issue #226 closes from pre-existing
reference-dependent CRAM fixture implementation and tests; 2.0.11 adds no duplicate fixture
code.

**Runtime quality gates now measure the root scripts.** `make type-check` includes
`scripts/` in mypy (#204), and all 35 Python files there are part of branch-inclusive
coverage: the final measurements were 93.94% scripts-only and 89.17% combined across
5,072 unit tests (#211).

**Mutation sweeps are crash-safe for production source.** Each mutant runs in a
disposable detached worktree over the current non-ignored working state, with import
provenance and a known-killed canary checked before atomic reports are installed in the
real checkout. Even an unhandled interruption can leave only the disposable worktree,
not a live mutant in production source (#208).

**Broad-exception behavior is measured rather than globally rewritten.** The reviewed
BLE001 inventory is 103 diagnostics normally and 108 including explicit suppressions;
fail-open handlers are classified by symbol and linked to behavior tests. BLE001 remains
outside Ruff's global selection, with no global fallback behavior rewrite (#219).

**Release automation is exact-commit and recoverable.** The controller verifies
exact-SHA CI and Docker evidence, promotes the tested GHCR digest to semantic-version and
`latest` aliases with idempotent retry/recovery, and separates unprivileged package building
from PyPI Trusted Publishing through OIDC (#214, #218). Publication and stable aliases
remain pending the first gated release.

**Presentation uses the generation name `VNtyper 2`.** Nine targeted `VNtyper 2.0`
presentation labels now say `VNtyper 2`, while semantic versions, historical statements,
identifiers and machine-readable names remain unchanged (#220).

## 2.0.10

Milestone 4, "CRAM and input robustness": closes #213, #225, #209, #178, #165 and
#161.

**A converted alignment that contains more than one read layout now stops instead of
discarding one layout.** In the measured golden cohort, 32 of 50 base BAM cases failed
rather than discard the non-empty single FASTQ produced beside their paired FASTQs.
VNtyper 2.0.9 silently ignored those reads; 2.0.10 names every produced FASTQ and record
count and refuses the mixed layout. This no-read-drop rule is a hard invariant, not a
configuration option. Regenerate the alignment with consistent pairedness and mate flags,
or provide a homogeneous paired-end or single-end FASTQ input; do not suppress the error
by discarding the named reads. Pure single-end BAMs and `--fastq1` without `--fastq2` are
now supported through fastp, BWA and Kestrel.

**Superseded on 2026-08-11 by #233.** The paragraph above accurately records the policy
shipped in 2.0.10. Current behavior retains the same no-read-drop invariant but routes
equal R1/R2 plus singleton/`other` reads losslessly; only unequal or one-sided mates fail
as invalid.

**BAM and CRAM inputs are preflighted through a run-local alignment view.** Missing
indexes are built beside that view, where htslib can resolve them, without writing into a
read-only patient input tree. The alignment and freshly generated index are opened and
retained for the plan lifetime; post-preflight random-access commands use that exact index
inode, so replacing the public BAI/CRAI name cannot retarget a slice or depth calculation.
Existing BAI, CSI and CRAI spellings are protected by format, and every derived
BAM/FASTQ/log destination is checked before a processing stage starts.

For BAM or CRAM input, the output root must stay outside the directory containing the
alignment and all of that directory's descendants. A layout such as `sample.bam` beside
`results/` is therefore rejected: `results/` is still inside the protected patient input
tree. Put alignments and results under separate roots, for example
`--bam inputs/sample.bam --output-dir results/sample/`, or choose an output path elsewhere.

**Reference-dependent CRAM fails before processing unless one candidate can decode the
real `-P` slice shape.** `--reference-fasta`, configured CRAM/BWA references and local
htslib resolution are tried in configured order; coverage receives the same selected
reference. Ambient network refget is disabled by default, `REF_PATH` is restored on every
exit, overlapping in-process CRAM scopes are serialized, and web jobs transport a curated
preflight code/message without exposing worker paths. Header-derived local references are
confined to the CRAM's own directory; use `--reference-fasta` or a configured reference
for an external path.

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

On the 18 golden-cohort BAM cases that both 2.0.9 and the measured milestone candidate can
complete without discarding reads, three alternating runs measured medians of 90.03 s and
87.52 s, respectively. The non-overlapping ranges (90.00–91.13 vs 87.49–87.83 s) record a
2.8% median improvement under the predeclared rule. The harness recorded candidate
`388f157`; its published-history runtime-tree equivalent is `470fdd6`. Later CRAM-reference
and single-end fixture changes are outside this measured BAM arm. Subsequent
alignment-binding identity hardening adds only a fixed, tiny number of metadata checks and
is timing-neutral by inspection; archive identity hardening runs only under
`--archive-results`, which the timing harness did not enable. This does not relabel the
measured hash as final `HEAD`.

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
