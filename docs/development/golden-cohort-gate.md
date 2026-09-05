# Golden-cohort gate (#179)

Before-versus-after comparison of genotyping output across the local test cohort, executed to verify that behaviour-changing commits in #179 preserve genotype stability.

**Current routing policy (2026-08-11, #233).** The mixed-layout refusals recorded in milestone-4 runs are historical. Equal R1/R2 plus singleton/`other` reads now consume every non-empty FASTQ exactly once under one Kestrel sample; unequal or one-sided mates remain invalid. Run 7's issue #233 `comparison.json` records 32 base cases plus 10 repeat/derived cases, 42 total, without rewriting earlier run evidence. This page also records the final gated candidate SHA and retained comparison artifact digests.

Every run selected as an attestation of record is registered in the table below. Completed superseded executions are not separate attestations: the preliminary issue #233 run at `49b0cc6` was replaced by run 7 at the final executable candidate. **A verdict is intended to attest one candidate commit and nothing after it**, so read the run whose candidate matches the tree you are judging.

**How candidate commits are recorded.** For runs 1 to 4, candidate revisions are operator assertions: harness `1.0.0` did not record `git rev-parse`, check clean working trees, or accept expected commit SHAs. Each side stored filesystem paths that change across branch checkouts. From harness `1.1.0` forward, starting with run 5, each execution records its `HEAD`, active branch, and working-tree status across `vntyper/`, `docker/`, and `scripts/`. The comparator rejects recorded revisions that diverge from `--expect-before-sha` or `--expect-after-sha`, and `--require-clean` aborts if uncommitted edits exist. **Run 5's candidate commit is a recorded cryptographic fact; runs 1 to 4 rely on operator logs.** Run 5's baseline (`4fd638a`) predates this harness change, so only its candidate side records verified revision metadata.

| Run | Candidate ("after") | Baseline | Verdict | Attests |
| --- | --- | --- | --- | --- |
| 1 | `7344c62` | `2fcc6e3` | PASS | The branch as it stood at `7344c62`, and nothing after it |
| 2 | `1792345` | `2fcc6e3` | PASS | The branch at `1792345`, including `d144505`, `c51052c`, `b4059ce`, `7e58eb8`, `2c92096`, and nothing after it |
| 3 | `8537a61` | `2fcc6e3` | PASS | The branch at `8537a61`, including `2ae28c5`, `2aa095a`, `50d7968`, `42c976a`, `4ce5639`, `97033d3`, `22e3d17`, `52a0ec9`, and nothing after it |
| 4 | `ec67fff` | `4fd638a` | PASS with two attributed deltas, neither genotype-affecting | The `fix/issue-181-197-followups` branch at `ec67fff`, against the 2.0.6 release, and nothing after it |
| 5 | `9816f86` | `4fd638a` | DELTAS, both classes fully attributed, every genotype artefact unchanged | The `fix/issue-181-197-followups` branch at `9816f86`, against the 2.0.6 release, and nothing after it |
| 6 | `48f97fe` | `cb593b6` | DELTAS, every one attributed, **no genotype field changed anywhere** | The `fix/milestone-2-correctness-of-reported-numbers` branch at `48f97fe` (milestone 2, #171/#172/#174/#203/#212), against the 2.0.7 release, and nothing after it |
| 7 | `19c8acd` | `4678851` | CANDIDATE PASS; comparison BLOCKED only by baseline-refused successes | Issue #233 at `19c8acd`, against its required regression baseline `4678851`, and nothing after it |
| 8 | `74fcbe0` | `c74e9e5` | DELTAS, every one attributed, **no genotype field changed anywhere** | adVNTR 2.0.x and real `--threads` (#259), against `c74e9e5`, and nothing after it |
| 9 | `edaf44a` | `80ac6be` | IDENTICAL (waived command deltas), **no genotype field changed anywhere** | Atomic BAM and BAI installation (#314), against `80ac6be`, and nothing after it |
| 10 | `936f11e` | `a632aa1` | DELTAS, every one attributed, **no genotype field changed anywhere** | Reporting floor split and profile revision 2 (#311), against `a632aa1`, and nothing after it |
| 11 | `2cf4946` | `a0d27b5` | IDENTICAL (waived command deltas), **no genotype field changed anywhere** | Derived confidence grade and report masthead chip (#173), against `a0d27b5`, and nothing after it |

Runs 1 to 3 evaluate branch `#179` against baseline `2fcc6e3`. Runs 4 and 5 evaluate branch `fix/issue-181-197-followups` against baseline `4fd638a` (the merge-base with `main` and release 2.0.6). **Run 5 supersedes run 4 as the attestation of record for that branch**: candidate `ec67fff` differs from the tip because Phase-5 fixes and review commits landed subsequently. Runs 1 to 3 remain authoritative records for their respective candidate commits.

No run attests branch tips as permanent properties: tips advance while commit SHAs are immutable. To verify whether a working tree matches a recorded run:

```bash
git diff --stat <candidate>..HEAD -- vntyper/ docker/
```

If the diff is empty, the checkout matches the tested codebase. Against run 5's candidate at release 2.0.7, only `vntyper/version.py` differs (a version string reported in logs that does not alter branching logic). `scripts/` is excluded because it contains test harnesses and validation scripts unimported by the pipeline.

**The genotype verdict is PASS in runs 1 to 6.** Every genotype field, `Confidence` label, and `Flag` matches byte-for-byte between baseline and candidate across all samples and assemblies. Run 6 confirms this quantitatively: exactly two columns were added, the only cells modified were within a new column, and **no genotype field was touched anywhere**. Run 7 introduces deliberate reachability updates: it demonstrates that 42 previously rejected input combinations execute successfully and losslessly, meaning no successful baseline outputs exist for direct comparison.

Runs 4 and 5 document non-genotype differences explicitly. In runs 1 to 3, all differences were restricted to report presentation. Run 4 includes two delta classes (a version-probe command string and cohort export file structure). Run 5 includes two deltas (the version-probe string and pruned internal working columns). Each difference maps to a specific commit and leaves genotype calls untouched.

| Metric | Run 1 | Run 2 | Run 3 | Run 4 |
| --- | --- | --- | --- | --- |
| Baseline ("before") | `2fcc6e3` (merge-base with `main`) | `2fcc6e3` | `2fcc6e3` | `4fd638a` (merge-base with `main`, release 2.0.6) |
| Candidate ("after") | `7344c62` | `1792345` | `8537a61` | `ec67fff` |
| Cases per side | 58 (50 BAM x assembly, 5 non-fast-mode, 3 adVNTR) | 58 (identical matrix) | 58 (identical matrix) | 58 (derived from `tests/data`), plus 4 cohort cases |
| Total executions | 116, plus 6 deliberate-mismatch probes | 116, plus 6 probes | 116, plus 6 probes | 130 (65 per side), probes and cohort cases included |
| Non-zero exits | 0 before, 0 after | 0 before, 0 after | 0 before, 0 after | 0 before, 0 after |
| Executed shell commands compared | No | No | **Yes** (480 per side) | **Yes** (1,111 per side across 61 cases) |
| Cohort mode covered | No | No | No | **Yes** (4 cases) |

Run 5 shares run 4's baseline and case count, but incorporates verified commit metadata and updated export deltas as detailed in the Run 5 section.

## Method

### Prerequisite: install reference assets

Every gate execution reads a shared `reference/` directory across both worktrees. Install and verify the bundle once from the published release:

```bash
vntyper install-references -d reference --references hg19 hg38 GRCh37 GRCh38 hg19_ensembl hg38_ensembl
```

Symlink the populated `reference/` directory into both baseline and candidate trees. The cohort contains all BAM fixtures in `tests/data/` across supported assemblies: 7 multi-reference samples across six assemblies (`hg19`, `hg38`, `GRCh37`, `GRCh38`, `hg19_ensembl`, `hg38_ensembl`), their original hg19 subsets, and the hg38 regression case `example_40cf`. Five cases run without `--fast-mode` to exercise unmapped read extraction, three run `--extra-modules advntr`, two run cohort CRAM inputs, and one purpose-built CRAM runs both scan strategies. **Run 6 is the first to include CRAM fixtures.** Fourteen additional cases alias `GRCh37`/`GRCh38` BAMs under `hg19_ncbi`/`hg38_ncbi` identifiers introduced in milestone 5.

**The evaluation matrix contains 78 cases (64 for run 7, 60 for run 6, and 58 for runs 1 to 5).** Figures across early runs reflect the 58-case matrix.

CRAM fixtures are generated dynamically (`scripts/make_cram_fixtures.py`). Fresh checkouts derive 72 cases (50 base, 5 non-fast, 3 adVNTR, 14 alias repeats) until CRAM files are synthesized. Two selected cohort CRAMs and the indexed-safe fixture run in both indexed and stream modes, adding six cases to reach 78.

Outputs compared per case: full `kestrel_result.tsv` rows and headers keyed on `Motifs`/`POS`/`REF`/`ALT`/`Variant`, pre-filter `kestrel_pre_result.tsv`, `output_adVNTR_result.tsv` when adVNTR executes, `coverage_summary.tsv`, HTML screening summaries, recorded pipeline steps, exit codes, and executed subprocess command strings.

Run 4 derives the 50 base cases from `tests/data` dynamically. Five non-fast cases, three adVNTR cases, and three probes follow declared policy (`scripts/golden_cohort/matrix.py`). Run 4 also introduces four `vntyper cohort` evaluations (`cohort_multi`, `cohort_multi_pseudonymized`, `cohort_single`, `cohort_empty`), diffing exports (`cohort_kestrel_{csv,tsv,json}`, `cohort_advntr_{csv,tsv,json}`), rendered tables, and pseudonymization outputs.

### The CRAM group (#188)

`make cram-fixtures` (`scripts/make_cram_fixtures.py`) generates lossless CRAM files matching cohort BAMs in `tests/data/cram/`. Decoded records match source BAMs exactly. Fixtures are built with `no_ref=1` because cohort BAM headers lack `M5` checksums, exercising container decoding, `.crai` indexing, and unmapped read extraction without external reference dependencies. Dedicated fixtures in unit and integration suites validate explicit reference resolution.

| Case | Source BAM | Records | Flag-4 reads (stream) | Guard count | Evaluation purpose |
| --- | --- | --- | --- | --- | --- |
| `b178_hg19_{indexed,stream}_cram` | `b178_hg19_subset` | 34,214 | 4,807 | 329 | Positive control with mixed layout. Raw `'*'` yields 4,478 reads; gate verifies 329 unplaced reads |
| `7a61_hg38_ensembl_{indexed,stream}_cram` | `7a61_hg38_ensembl_bwa` | 985,731 | 634,261 | 11,571 | High unmapped load; verifies stream safety against race conditions and indexed-scan data loss |
| `indexed_safe_{indexed,stream}_cram` | Purpose fixture | 40 | 20 | 0 | Authorized paired fixture: idxstats reports zero unmapped placed reads; indexed `'*'` and flag-4 stream match at count 20 |

All CRAM cases execute without `--fast-mode` to force unmapped read extraction.

### Package resolution verification

To prevent editable install collisions where `sys.meta_path` resolves to the developer checkout regardless of process CWD, each test launches through a wrapper that configures `PYTHONPATH` and verifies that `vntyper.__file__` resides inside its own tree. The wrapper asserts marker module presence: `vntyper.scripts.pipeline_guards` (runs 1 to 3) or `vntyper.scripts.cohort_rules` (run 4) must be absent on baseline and present on candidate checkouts. All gate runs pass this verification.

## Result: run 1, candidate `7344c62`

| Metric | Cases with delta |
| --- | --- |
| Exit code | 0 / 58 |
| Kestrel variant calls (rows added or removed) | 0 / 58 |
| Kestrel `Confidence` | 0 / 58 |
| Kestrel `Flag` | 0 / 58 |
| Other Kestrel columns | 0 / 58 |
| adVNTR variant calls, fields, `Flag` | 0 / 3 |
| Coverage summary | 0 / 58 |
| Recorded pipeline steps | 0 / 58 |
| Screening summary text | 1 / 58 |
| Screening summary emphasis | 11 / 58 |

`kestrel_result.tsv`, `kestrel_pre_result.tsv`, `output_adVNTR_result.tsv`, and `coverage_summary.tsv` match identically after removing timestamped `##` headers.

Calls match across all assemblies and both checkout sides:

| Sample | Motifs | Variant | Confidence | Flag |
| --- | --- | --- | --- | --- |
| `example_6449` | `4-5` | Insertion | `High_Precision*` | Not flagged |
| `example_66bf` | `5C-Q` | Insertion | `High_Precision*` | Not flagged |
| `example_6c28` | `S-Q` | Insertion | `High_Precision*` | Not flagged |
| `example_7a61` | None | None | `Negative` | None |
| `example_a5c1` | `L-6p` | Insertion | `High_Precision` | Not flagged |
| `example_b178` | `D-C` | Insertion | `High_Precision*` | Not flagged |
| `example_dfc3` | `5-E` | Deletion | `High_Precision*` | Not flagged |
| `example_40cf` (hg38) | None | None | `Negative` | None |

## Result: run 2, candidate `1792345`

Same baseline and 58-case matrix. Re-executed to incorporate production commits `d144505` (adVNTR compound repair), `c51052c`/`b4059ce` (assembly validation), `7e58eb8` (shell quoting), and `2c92096` (`--output-name` checking).

| Metric | Cases with delta |
| --- | --- |
| Exit code | 0 / 58 |
| Kestrel variants, confidence, flags, columns | 0 / 58 |
| `kestrel_pre_result.tsv` | 0 / 58 |
| adVNTR variants, fields, flags (excluding intermediate `Insertion_len`) | 0 / 3 |
| Coverage summary | 0 / 58 |
| Pipeline step records | 0 / 58 |
| Screening summary text | 1 / 58 |
| Screening summary emphasis | 11 / 58 |
| Rendered `Motif` table cell (display only) | 48 / 58 |
| HTML entity escaping (display only) | 58 / 58 |

adVNTR calls reproduce run 1 exactly across all cases (VID 25561): `a5c1` (RU 2, P=6.78e-07), `b178` (RU 4, P=3.83e-56), and `dfc3` (`D17_2&D18_2&D19_2&D20_2&D21_2`, RU 2,2,2,2,2, flagged `Polymorphic_Call`).

Assembly verification: the assembly guard returned `agree` across all 58 candidate cases (20 `hg19`, 9 `hg38`, 8 `GRCh38`, 7 each `GRCh37`, `hg19_ensembl`, `hg38_ensembl`) with zero conflicts.

## Result: run 3, candidate `8537a61`

Incorporates `2ae28c5` (`shlex.quote` applied across subprocess boundaries).

| Metric | Cases with delta |
| --- | --- |
| Exit code | 0 / 58 |
| Kestrel variant calls, confidence, flags | 0 / 58 |
| Quoted command strings (quickcheck, SAM->BAM, index, bcftools sort) | **0 / 58** |
| Screening summary text / emphasis | 1 / 11 of 58 |
| Cross-match emphasis (`dfc3_hg19_advntr`, display only) | 1 / 3 adVNTR cases |
| Rendered `Motif` cell / HTML entity escaping / IGV fragments | 48 / 58 / 58 of 58 |

Command string comparisons: 480 commands diffed per side. Quoting was byte-identical across all reached commands because filenames contained no shell metacharacters. Differences were restricted to `set -o pipefail; ` prefixes on sort pipes and header inspection order. Cross-match styling on `dfc3_hg19_advntr` changed from `summary-positive` to standard `summary-box` styling to reflect non-matching results (`2aa095a`).

## Result: run 4, candidate `ec67fff`

Branch `fix/issue-181-197-followups` compared against 2.0.6 release baseline `4fd638a`. Covers 58 per-sample cases, 3 probes, and 4 cohort cases (65 runs per side, 130 total).

| Artifact | Cases with delta | Total compared |
| --- | --- | --- |
| `kestrel_result`, `kestrel_pre_result` | **0** | 59 |
| `advntr_result`, `cross_match_summary` | **0** | 3 |
| `coverage_summary` | **0** | 59 |
| `exit_code` | **0** | 65 |
| `pipeline_step_records` | **0** | 61 |
| `cohort_category_counts`, totals, tables | **0** | 3 |
| `cohort_output_files` | **0** | 4 |
| `executed_commands` | 61 | 61 |
| `cohort_kestrel_{csv,tsv,json}` | 3 each | 3 each |
| `cohort_advntr_{csv,tsv,json}` | 3 each | 3 each |

Attributed deltas:
1. `executed_commands` differed across 61 cases due to removing a duplicated `-h` flag in Kestrel help probes (`2873ad3`). Total commands matched at 1,111 per side.
2. `cohort_*` exports removed leaked internal tracking columns `__row_result` and `__unified` (`90f61fa`). Category counts and summary tables remained identical.

## The evaluation instrument: harness `1.0.0` versus `1.1.0`

Harness `1.1.0` introduced cryptographic verification of revisions, strict matrix bounds, and md5 checks on unparsed sidecars (`pipeline_info.json`, `output_R1.fastq.gz`, `cross_match_results.tsv`). Re-evaluating run 4 outputs with harness `1.1.0` confirms identical outcomes while preventing silent skips on reduced matrices.

## Correction: adVNTR output columns in runs 2 and 3

Early drafts stated that `Insertion_len` was checked in `output_adVNTR_result.tsv`. `Insertion_len` is an internal variable used in frameshift calculation that is not written to disk. The 0/3 delta metric applies to public table columns (`VID`, `Variant`, `NumberOfSupportingReads`, `MeanCoverage`, `Pvalue`, `RU`, `POS`, `REF`, `ALT`, `Flag`).

## Adjudication of differences

- **D1 (Screening emphasis, `5527a49`):** Negative cases lost erroneous `summary-positive` CSS styling. Wording remained unchanged.
- **D2 (Screening sentence on `dfc3_hg19_advntr`, `77d590b`):** Added 15 missing rule mappings in `report_config.json`, replacing a negative fallback with an accurate description of high-precision Kestrel detection and flagged adVNTR calls.
- **D3 (Assembly validation, `078a6c4`):** Correctly identified assemblies across all 58 cases without raising false rejections.
- **D4 (adVNTR `Repeat_Unit_7`, `52f822e`):** Inactive in cohort because no sample carries RU 7.
- **D5 (Compound adVNTR parsing, `a7c3d9e`, `d144505`):** Retains greedy token parsing while constraining string splitting, preserving genotype output.
- **D6 (Pipefail enforcement, `331ea95`):** Subprocess pipes safely exit on upstream errors.

## Run 5: `4fd638a` -> `9816f86`

Evaluates nine commits resolving adversarial review findings for PR #199 using harness `1.1.0`.

- Baseline: `4fd638a` (v2.0.6)
- Candidate: `9816f867c28f` on `fix/issue-181-197-followups`, verified clean
- Executions: 65 per side; package resolution verified
- Verdict: **DELTAS** (all deltas attributed; zero genotype changes)

| Metric | Cases with delta | Total compared |
| --- | --- | --- |
| `kestrel_result`, `kestrel_pre_result` | 0 | 59 |
| `advntr_result`, `cross_match_summary` | 0 | 3 |
| `coverage_summary`, `report_tables` | 0 | 59 |
| `screening_summary`, `exit_code` | 0 | 59 / 65 |
| `cohort_category_counts`, totals, tables | 0 | 3 |
| `cohort_output_files` | 0 | 4 |

Deltas: single command replacement on Kestrel version probes (D6) across 61 cases; cohort export columns pruned of internal fields (29 columns vs 31 in Kestrel, 13 vs 15 in adVNTR).

## Run 6: `cb593b6` -> `48f97fe`, milestone 2

Re-attested against 60-case matrix (including 50 derived CRAM fixtures).

| Metric | Cases with delta | Total compared |
| --- | --- | --- |
| `coverage_summary`, `kestrel_pre_result`, `report_tables` | 61 | 61 |
| `executed_commands` | 61 | 63 |
| `kestrel_result` | 50 (negative sentinels unchanged) | 61 |
| `cohort_kestrel_*`, `cohort_stats_*`, `cohort_tables` | 3 | 3 |
| `advntr_result`, `screening_summary`, `exit_code` | **0** | 3 / 61 / 67 |

Quantitative verdict: exactly two columns added across runs (`flag_filter_pass` and `coverage_qc`). The only cell values modified resided in `flag_filter_pass`. No genotype field changed. Corrected mean coverage satisfied the exact mathematical relation `mean_new == mean_old * (1 - pct_old / 100)` across all 61 cases (#171).

## Run 7: `4678851` -> `19c8acd`, issue #233

Harness `1.4.0`. Clean detached worktrees on both sides. Evaluated 64 pipeline cases, three probes, and four cohort sets. The candidate met **67/67 pipeline/probe outcomes** and **4/4 cohort outcomes** without timeouts, aborted runs, or unverified launches.

The comparison verdict is `BLOCKED` because baseline `4678851` intentionally refused 42 mixed-layout FASTQ declarations that candidate `19c8acd` processes losslessly. Exactly 25 cases ran successfully on both sides, with identical genotype calls and exit codes. Command deltas across 23 cases were restricted to process-substitution FIFO descriptors (`/proc/<pid>/fd/5`) and temporary path noise.

Candidate execution confirmed exact lossless routing across 40cf's 93 singleton reads, b178's singletons, all three adVNTR successes, CRAM streams, and forced-indexed fail-closed guards. This page also records the final gated candidate SHA and retained comparison artifact digests:

| Artifact | SHA-256 |
| --- | --- |
| `comparison.json` | `5b8dc9199cd19fc1142e0a6ba7bd2740d4c0a97b0cdd9e5f8f4b08e51330e88e` |
| `comparison.txt` | `6808936b98be8b8d79decd17c76f89f5f4519a6e1fa9acc3f96c0c9eb6d14cbd` |
| Baseline `side.json` | `d3b17029f55c4a610d708764bf4b9c5298f2caad3f0f3114ce532b79b43b41a3` |
| Candidate `side.json` | `8a0c1a0460934cecf9db19b659c7f219f964bf685bf5f818bf12b2b3b69bac10` |
| Matrix snapshots | `6f09f9350d152ab1b69aa07cf2096aad895a01cca08f23c297608cf772029dd0` |

## Run 8: `c74e9e5` -> `74fcbe0`, adVNTR 2.0.x and real `--threads` (#259)

Harness `1.4.0`. 78 pipeline cases, 3 probes, 4 cohort sets (85 runs per side). Package resolution verified. Distinguishing marker: `vntyper.modules.advntr.advntr_genotyping:resolve_advntr_threads`.

Both sides executed installed adVNTR 2.0.2 with identical references. The gate isolated thread inheritance: baseline ran adVNTR at `-t 1` while candidate ran at `-t 8`.

Genotype outputs: 0 deltas across `advntr_result`, `kestrel_result`, `coverage_summary`, `screening_summary`, and all cohort artifacts. Executed command deltas in 76 cases were limited to `/proc/<pid>/fd/N` process noise; the 3 adVNTR cases demonstrated `-t 8` parameter delivery without genotype deviation.

Artifact digests:
- `result.json`: `3979e774ccaed31e8877f7c2441980849a59a1b44b5fb95bcf6a4082cb992751`
- `result.md`: `a62a0fa5d6064ba6fc4c64c1719c3d3f0b55a0ed30540a7f50f7f769eff943a3`
- Baseline `side.json`: `5321f9206405b77ccf4a198603b1ad1ca1b2f628105aba7a38f8cadaeaddb075`
- Candidate `side.json`: `c8f583664a48615116573697ce1d3114ec45a13482fdf183d7c5d2e757f79335`
- Matrix snapshots: `241aeaf5aa64b8f684848f271e9d1a45e422ba41825adda9ee67b6a32cd7ac68`

## Run 9: `80ac6be` -> `edaf44a`, atomic BAM and BAI installation (#314)

Harness `1.5.0`. 85 runs per side over 78 pipeline cases. Distinguishing marker: `vntyper.scripts.artifact_publish`.

Baseline wrote output BAM and BAI files directly to public target paths. Candidate writes to deterministic sibling `.partial` paths (`<path>.partial`), validating zero exit status and non-empty file creation before atomic `os.replace` publication. On error, `.partial` files are unlinked, preventing corrupt intermediates from populating results.

Verdict: **IDENTICAL** under `--expect-command-delta`. All genotype and QC artifacts match. Command string deltas in 79 of 81 cases reflect `.partial` target names and process descriptor noise.

Artifact digests:
- `result.json`: `3c84f5970a0a145e13116bc44ec46951e225cdbee283b4ae6dc76f0878f54b7b`
- `result.md`: `13344703ac581b5ca75aca6eaa8d8ef1d98786f5b6f8e5b43e31e2d3b2305ef6`
- Baseline `side.json`: `424e33c894b032987cd1971a1b79dfb44c5916c1cf4d4603857cd0b0346d1148`
- Candidate `side.json`: `00bedbf8d8fb1f24b78d91d44b29915dd018c91a7897ca7d9a7627127eb676f9`
- Matrix snapshots: `e3a0509d30d646f836eb129b79edded4960b5890673b2e31cc0f89835933e940`

## Run 10: `a632aa1` -> `936f11e`, reporting floor split and profile revision 2 (#311)

Harness `1.5.0`. 85 runs per side over 78 pipeline cases. Distinguishing marker: `vntyper.scripts.confidence_assignment:HAS_REPORTING_FLOOR_SPLIT`.

Candidate splits the numerical reporting floor from the `Low_Precision` boundary (`0.00469`), mandates `gg_depth_score_threshold`, and updates the packaged decision profile to revision 2 (`0b13d07370491b3ea773e65144891cb30caebcae70b0ef98feb0f2c5ccd2f4a1`).

Because the numerical floor remains 0.00469, all genotype calls, scores, confidence tiers, and summary tables match baseline outputs. Observed deltas are limited to metadata fields in cohort exports (`Decision_Profile_Revision`: `'1'` -> `'2'`) and command string process noise.

Artifact digests:
- `result.json`: `51fe52a3f6b49e7b45e82ad8e09952571b9ef763c445c0113e2a177308f17ed1`
- `result.md`: `235122e7c5c8defcac1e74496a06efc75f71b12873f1449b86e74ea71217f0b1`
- Baseline `side.json`: `5ae4bdc4c3607ff7f08815c16db8beb78f877ca528fa8124b624ea58f4dcaf5e`
- Candidate `side.json`: `8596731beb38755b33470bf3def84e6e899e26833444b9699b9769f7ca3fca8f`
- Matrix snapshots: `d218de77d9db2a7802015a1c76aedd51fd4c61c6dddaef87dc489cde3b738339`

## Run 11: `a0d27b5` -> `2cf4946`, derived confidence grade and report masthead chip (#173)

Harness `1.5.0`. 85 runs per side over 78 pipeline cases. Distinguishing marker: `vntyper.scripts.screening_summary:supports_confidence_grade`.

Candidate calculates derived sample-level confidence grades and displays a corresponding chip in the HTML report masthead (`summary_report.html`) under `state_chips`.

All variant calls, depth scores, confidence tiers, flags, TSVs, BED files, and summary sentences match baseline outputs. Command string differences reflect ephemeral process IDs and temporary path noise.

Verdict: **IDENTICAL** under `--expect-command-delta`.

Artifact digests:
- `result.json`: `fb9565a222b7ebdd6175ef58216537572c0032ab46afb534a4a60c4a67044100`
- `result.md`: `da72d53de325f4516e63f7b85b92aaa3813bbd3c5b1d398b3d292539afd871c6`
- Baseline `side.json`: `5852493a1ddb374e5dce2f059eec4aa3555c80595a71f916c8e717d30798a83b`
- Candidate `side.json`: `4b897bd131fc82197d45110dd198fcb33f2e6568e5a2f014e61616fbdb8bca95`
- Matrix snapshots: `bf8555cd7e734cf8e723e80423ecdbd30e9f83459d8d5e4a42313b63ec82881f`
