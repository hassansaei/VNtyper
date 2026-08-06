# Close comments for #181–#197

One section per issue, ready to paste as a comment before closing. **Do not close any of
these until the golden-cohort run at `ec67fff` has reported**: the six gated issues
(#184, #185, #188, #192, #195, and #196's ratchet) each carry a
`GOLDEN-COHORT: <fill in>` line that must be replaced with the real verdict and run
reference before the comment is posted. The four non-gated issues can go as they stand.

Branch: `fix/issue-181-197-followups`, merge base `4fd638a` (v2.0.6).

Two issues also need their **body text corrected** before closing — #188 and #192. The
correction is included in each close comment rather than left for the next reader to
inherit from a closed issue.

---

## #181 — Frameshift validity rule discards every (3n+1)-bp deletion

Closed by **`88cb88e`** (`docs(scoring): record the +1-frame rule as specification`).

Resolved as **specification, not a code change.** @hassansaei decided the rule is correct
as written: the insertion arm accepts 3n+1 and the deletion arm accepts 3n+2, and the
apparent asymmetry is the intended convention rather than a bug.

**What now pins it**

- `docs/pipeline/scoring-and-confidence.md` records the +1-frame rule as specification,
  with the decision quoted.
- `tests/unit/test_scoring.py::test_a_one_bp_deletion_is_not_a_valid_frameshift_by_design`
  — the one new test, upgraded to a specification docstring that quotes the decision and
  cites this issue.
- The neighbouring tests
  (`test_a_3n_plus_1_bp_deletion_is_a_frameshift_that_the_rule_discards`,
  `test_a_3n_plus_2_bp_deletion_is_the_only_deletion_the_rule_accepts`,
  `test_the_insertion_arm_accepts_3n_plus_1_only`) stay **characterisation** deliberately:
  they also assert direction, `is_frameshift` and row retention, none of which was decided
  here, and upgrading them would ratify behaviour nobody signed off.

**Evidence the test can fail:** mutating `scoring.py:160` produced the expected failure;
reverted, `git diff --quiet -- vntyper/` clean.

No executable line of `scoring.py` changed, so no cohort run is required.

---

## #182 — adVNTR frame filter drops deletions ≡1 (mod 3) and insertions ≡2 (mod 3)

Closed by **`513890e`** (`test(advntr): assert Kestrel and adVNTR share one frameshift
convention`).

Resolved as **specification**: the adVNTR filter and the Kestrel scorer implement the
same convention, and that identity is now asserted rather than assumed.

**What now pins it**

- `tests/unit/test_frameshift_convention_parity.py::test_kestrel_accepts_insertions_at_3n_plus_1_and_deletions_at_3n_plus_2`
- `tests/unit/test_frameshift_convention_parity.py::test_advntr_real_filtering_behaviour_matches_kestrel_for_every_residue`
  — parametrised over every residue class, driving the **real** `advntr_processing_ins` /
  `advntr_processing_del` against `scoring.extract_frameshifts`.
- Two docstrings in `tests/unit/test_advntr_frameshift_filter.py` upgraded to
  specification.

**A correction worth recording.** The drafted parity test was **tautological**: it rebuilt
`ins_frame` / `del_frame` locally and then asserted properties of its own construction.
Those arrays are not module constants — they are built inside
`advntr_processing_ins`/`_del` from `advntr_settings`, which is loaded at import time
(AGENTS.md trap 1). The committed test derives both from the real functions, so it can
actually fail on drift.

**Evidence it can fail:** moving the insertion offset +1 → +2 failed 6 of 9 parametrised
cases (lengths 1, 2, 4, 5, 7, 8 flip; 3, 6, 9 do not), then reverted.

No executable line of `advntr_genotyping.py` changed in this commit — only two comment
blocks were added, verified byte-identical against `513890e^`.

---

## #183 — Decide whether confidence-tier precedence should be specified

Closed by **`451a4b7`** (`test(confidence): specify last-wins tier precedence, per #183`).

Resolved as **specification**: the region-depth demotion being overwritten by a later
`High_Precision` tier is the intended last-wins ordering, not an accident.

**What now pins it**

- `tests/unit/test_confidence_boundaries.py::test_region_depth_demotion_is_overwritten_by_a_later_high_precision_tier`
- `test_the_region_depth_clause_is_only_reachable_through_the_alt_threshold_gap`
- `test_a_low_region_depth_row_is_deliberately_not_capped_at_low_precision`
- The 54-cell boundary matrix (`test_confidence_boundary_matrix`), whose docstring names
  the exact future collision — **#184 changing one of its cells** — so a reader cannot
  mistake the matrix for specification.
- `docs/pipeline/scoring-and-confidence.md` records the decision.

`confidence_assignment.py` was temporarily mutated to prove the test can fail, then
reverted (`git diff --quiet` clean). The last-wins arithmetic was independently
re-derived by hand against `confidence_assignment.py:149-154`.

Note for the record: #184 subsequently *did* change cells of that matrix — five of them,
not the one the plan predicted. The matrix docstring anticipated exactly that.

---

## #184 — cond3 in confidence_assignment.py never changes an outcome

Closed by **`6f03ddf`** (`fix(confidence): keep mid-band Depth_Score at Low_Precision, and
require the calibration constants`).

@hassansaei: *"Any variant with Depth_Score between 0.00469 and 0.00515 (inclusive) must
be Low_Precision, even when alternate depth is >= 21."*

`cond6` already demoted the **open** interval at every alt depth because it was applied
last, so the only divergence from that intent was `Depth_Score == 0.00515` exactly, which
`cond2`/`cond5` promoted. `cond6` becomes `cond_midband`, inclusive at both ends, still
applied last. `cond3` is retained as the expression that names the intent — *"do not
remove this intent"*.

**The second suggested implementation was not used, and the reason is measured rather than
argued.** Changing `cond2`/`cond5` to `> high` alone leaves `alt_depth >= 100,
DS == 0.00515` matching **no condition at all**, so it falls through to the negative
label. Induced: it produced `Negative` on all ten boundary pairs — a reported call
becoming a non-call, strictly worse than the bug.

**Correcting this branch's own analysis: the spec said one matrix cell moves. Five do.**
The 54-cell matrix constructs *fractional* region depths, so the boundary is reachable at
every alt depth there, not only at multiples of 103. The sixth cell survives because its
region depth cannot land exactly on the threshold and sits one ULP above.

**What now pins it** (`tests/unit/test_confidence_boundaries.py`)

- `test_the_top_of_the_mid_band_is_low_precision`, `test_the_bottom_of_the_mid_band_is_low_precision`
- `test_one_ulp_above_the_mid_band_is_still_promoted`
- `test_no_row_at_the_boundary_falls_through_to_negative` — the guard against the
  alternative implementation above
- `test_a_boundary_demotion_changes_which_variant_is_reported` — measured in both
  directions: before POS 67 / `High_Precision*`, after POS 68 / `High_Precision`
- `test_cond5_is_unreachable_at_the_boundary_on_integer_depths`,
  `test_a_fractional_region_depth_reaches_the_boundary_inside_cond5s_window`,
  `test_production_depths_arrive_as_whole_numbers`

Also in this commit, and **not** part of #184: the six `.get()` calibration fallbacks are
now direct subscripts. They encoded a second, wrong calibration — 0.4 where the shipped
config says 0.00515, a factor of 78, reachable by dropping one key. Measured on the old
code: dropping `depth_score_thresholds.high` turned `High_Precision*` into `Low_Precision`
on an ample-depth row. Behaviour is identical under the shipped config, which supplies all
six. Their `EQUIVALENT_MUTANTS` entries are removed from `scripts/mutation_test.py` — the
mutants are deleted rather than hand-excused.

`confidence_assignment.py` is at 100% line and branch coverage.

**GOLDEN-COHORT: `<fill in — run at ec67fff>`.** Read the result with this caveat, which
is stated in the commit and repeated here so the closed issue carries it: **a cohort PASS
is weak evidence for this change.** Exact float equality with 0.00515 requires alt depth
to be an exact multiple of 103, essentially unreachable on real data, so the cohort very
likely contains no boundary row at all. The boundary table and the selection test are the
load-bearing evidence.

Noted, not fixed: `cond3` is now dead in effect and its operands will surface as
equivalent-mutant survivors in a future sweep; and the module docstring's "Convert depth
columns to integers" has always been false — it converts to float, which is precisely the
fact the integer-depth scoping above depends on.

---

## #185 — filter_final_dataframe fails open when a gate column is missing

Closed by **`c489d49`** (`fix(kestrel): raise when a required gate column is missing, per
#185`).

@hassansaei: *"a missing required gate column should raise (abort the run), not be
skipped [...] That is not acceptable for this pipeline."*

AGENTS.md trap 4: the Kestrel stages **mark**, they do not filter. Each appends a boolean
column and `filter_final_dataframe` ANDs them at the end — so a missing column turned a
safety gate into a permit and let variants through that should have been filtered.

All five gates are parametrised so none is left unguarded. `kestrel_pre_result.tsv` is
still written **before** the raise, so a failing run keeps its evidence, and the empty
frame short-circuits after that write — the explicit empty-result path the decision carves
out.

**The named risk was that the raise would be swallowed.** It is not: `kestrel_genotyping.py`
contains zero `try`/`except`, so the chain reaches `pipeline.py:715`, which logs and exits
1, and `SystemExit` is a `BaseException` so it survives `cli.py:97`.

**What now pins it** (`tests/unit/test_kestrel_filtering.py`)

- `test_a_gate_that_is_absent_from_a_non_empty_frame_raises` (parametrised over all five)
- `test_every_gate_missing_at_once_raises_on_the_first_one`
- `test_an_empty_frame_is_the_documented_empty_result_path`
- `test_the_pre_result_tsv_is_still_written_before_the_raise`
- `test_the_only_caller_does_not_swallow_the_raise`
- `test_a_missing_gate_column_aborts_the_run_rather_than_reporting_nothing` — an
  **end-to-end** test driving the real `run_pipeline` with a stage made to stop emitting
  `motif_filter_pass`, asserting `SystemExit`, code 1, and the message at ERROR. Before
  the fix it failed with *"the run did not exit; it returned None"*: the pipeline
  completed normally with a gate silently skipped.

**The check that mattered most was not in the plan:** whether any stage can *legitimately*
omit its gate, in which case the raise would abort healthy runs. All five were read; every
assignment is unconditional on each non-empty return path.
`motif_correction_and_annotation` looks risky because its `keep_cols` drops most columns,
but it returns the untouched deep copy.

**GOLDEN-COHORT: `<fill in — run at ec67fff>`.** Expected no-op — the raise should never
fire on the cohort. A cohort *difference* here would mean a stage really is omitting a
gate, which is a finding, not a reason to weaken the raise. A PASS proves only that no
shipped stage omits a gate on the cohort's inputs.

Also visible to the web service: a job that previously produced a wrong result can now
surface as failed.

---

## #186 — Do not enable use_uniform_filtering

Closed by **`99affb2`** (`fix(motif): refuse uniform filtering with an empty GG
allowlist`).

@hassansaei's comment on this issue was pulled with `gh api` and verified quoted exactly
rather than paraphrased.

The uniform branch would delete the canonical dupC call when the GG allowlist is empty.
Rather than document "do not turn this on", the code now refuses the combination: 17 added
lines in `motif_processing.py`, **no deletions**, and the guard cannot fire on the shipped
config.

**What now pins it** (`tests/unit/test_motif_config_guard.py`)

- `test_the_shipped_config_leaves_uniform_filtering_off`
- `test_uniform_filtering_with_an_empty_gg_allowlist_is_refused`
- `test_uniform_filtering_with_a_populated_allowlist_is_permitted`

`tests/unit/test_motif_filtering_issue_136.py` is unaffected (20 pass). The diff is
additive-only — no `-` lines in `motif_processing.py` — so no genotype can move on the
shipped config, by construction.

**A brief error found while implementing, recorded because it would have masked the
guard:** the drafted `_one_row()` fixture omitted the `Variant` column that
`motif_correction_and_annotation`'s `keep_cols` requires, so both non-trivial tests would
have died on an unrelated `KeyError` instead of exercising the guard.

---

## #187 — SHARK accepts reference_assembly and ignores it

Closed by **`39fe0eb`** (`docs(shark): state that reference_assembly does not select a
region`).

Resolved as **documented behaviour plus a warning**, not a region-selection feature: the
shipped config offers only an hg19 region, so there is nothing for a non-hg19 value to
select. Passing `hg38` now warns; the filtering itself is unchanged.

**What now pins it** (`tests/unit/test_shark_filtering.py`)

- `TestReferenceAssemblyIsAccceptedAndIgnored` — including
  `test_hg38_produces_a_byte_identical_command_to_hg19`,
  `test_an_hg38_run_still_filters_against_the_hg19_region`, and
  `test_the_shipped_config_offers_only_an_hg19_region`
- `TestNonHg19AssemblyLogsAWarning` — `test_a_non_hg19_assembly_is_warned_about`
  (hg38/GRCh38/nonsense) and `test_hg19_and_grch37_do_not_warn`
- `docs/pipeline/optional-modules.md` states it.

**Both directions were verified empirically rather than by reading:** deleting the warning
block fails exactly the 4 tests that assert it and nothing else; making the warning
unconditional fails exactly the 3 negative (hg19/GRCh37/grch37) tests. So
`pytest.ini`'s `log_level = DEBUG` does not produce a false pass here.

A brief error found: it assumed `run_command` returns `(ok, stdout, stderr)`. It returns a
plain `bool` (`utils.py:17-74`), so the file's existing command-capture spy was used
instead.

---

## #188 — Web service passes --bam for CRAM alignments

Closed by **`5f5222e`** (`fix(web): pass --cram for CRAM uploads, and fix the two
regressions that alone would cause`).

### Correcting this issue's diagnosis

**This issue describes a three-production-file change. That is wrong, and leaving it wrong
in a closed issue would mislead the next reader.**

`--cram` has existed since `cli_parser.py:92`; `cli_handlers.py` threads it and
`pipeline.py` branches on it. **Only the task layer was wrong.** But it is not "one line"
either — it is **three lines across two files, and fixes 2 and 3 exist only because of
fix 1**:

1. `docker/app/tasks.py` hardcoded `--bam` for every alignment upload, so every accepted
   CRAM took the BAM code path.
2. `cli_handlers.py` derives the sample name from `--bam` or `--fastq1` and **never from
   `--cram`**. Today a CRAM arrives via `--bam` and gets its file stem; **after the flag
   fix alone**, every CRAM run without an explicit `--sample-name` would be called the
   literal string `"sample"`, in the report and in the filenames.
3. `tasks.py` fell back to `f"{bam_path}.bai"` for every format, but `samtools index` on a
   CRAM writes `.crai` — so cleanup removed a file that never existed and left the real
   index on the shared volume.

So the one-line fix this issue implies would have shipped a regression. That is the part
worth carrying forward.

`build_vntyper_command`, `resolve_index_path` and `is_cram` are extracted to module level
first — verbatim apart from the flag — because the command was built inline inside
`run_vntyper_job` and the flag choice could not otherwise be tested without Celery.
Extension matching is case-insensitive because `uploads.py` accepts `SAMPLE.CRAM`.

**What now pins it**

- `tests/unit/web/test_cram_alignment_handoff.py` (12 tests):
  `test_a_cram_upload_is_passed_to_the_cli_as_cram_not_bam`,
  `test_a_bam_upload_still_uses_the_bam_flag`,
  `test_the_flag_is_chosen_case_insensitively`,
  `test_every_other_flag_is_unchanged_by_the_extraction`,
  `test_a_cram_upload_with_no_index_falls_back_to_the_crai_name`,
  `test_a_cram_run_derives_its_sample_name_from_the_file_stem`, and the rest.
- `tests/unit/web/test_index_handoff.py` is **corrected** in the same commit: its fake
  `samtools index` wrote `.bai` for a CRAM — the stub encoded the defect production had.
  With the corrected stub and the old production code, that test now fails, which is what
  it should have been asserting all along.

**GOLDEN-COHORT: cannot attest this.** The cohort contains no CRAM input at all (#191's
trigger). The evidence is instead an end-to-end run through the real worker with only
Redis faked, on a CRAM built at test time into scratch (never into `tests/data/`): the
command carries `--cram`; `pipeline_summary.json` records `input_files={"cram": …}`; the
VCF sample column is the file stem; the `.crai` is produced and cleaned up; and
`kestrel_result.tsv` is **byte-identical to the BAM run of the same sample** apart from
the timestamp (93 / 9883 / 0.009410098148335525 / High_Precision).

Two induced reverts make it conclusive rather than plausible: reverting to the pre-#188
state leaves `example_a5c1_hg19_subset.cram.crai` on the volume after the CRAM is deleted
(the leak is measured, not predicted), and shipping fix 1 alone yields a VCF sample column
of `sample`.

**Stated rather than glossed:** the genotype-equality assertion is **hand-run, not in
CI** — it needs samtools, reference data and a real pipeline subprocess, all barred from
the unit tier. #188 does not park on the #191 trigger, but it stays unattested in CI until
an integration case exists.

A separate defect found while gathering this evidence, filed separately and **not** fixed
here: every job leaves `<alignment>.quickcheck.log` in its per-job input directory, so the
directory is never reclaimed. It affects BAM and CRAM alike and is independent of #188.

---

## #192 — Define adVNTR Insertion_len semantics for compound variants

Closed by **`2a267aa`** (`fix(advntr): sum every LEN token in a compound state, per
#192`).

@hassansaei: *"use the sum of inserted lengths [...] `I9_2_A_LEN9&I50_2_A_LEN3` to
Insertion_len = 9 + 3 = 12 (not first-LEN-only). [...] Do not keep first-LEN-wins as the
defined semantics."*

### Correcting this issue's description of the old behaviour

**This issue, the #198 PR body and this branch's own spec all describe the previous
behaviour as "first-LEN-wins". It was not, and "first-LEN-wins" describes no input at
all.**

The greedy extract took everything from the first `LEN` to the end and the bounded split
then returned the **remainder**, so the value **collapsed to 0 whenever material followed
the first `LEN` token**, while a compound whose only `LEN` is terminal parsed correctly.
Measured against the real transform at the merge base:

| State | Insertion_len (old) | Deletion_length | frame |
| --- | --- | --- | --- |
| `I22_2_G_LEN1` | 1 | 0 | 1 |
| `I9_2_A_LEN9&I50_2_A_LEN3` | **0** | 0 | 0 |
| `I9_2_A_LEN2&D50_2` | **0** | 1 | 1 |
| `I50_2&I9_2_A_LEN3` | 3 | 0 | 3 |
| `I50_2&D9_2&I80_2_A_LEN7` | 7 | 1 | 6 |
| `I9_2_A_LEN3&D50_2&D51_2` | **0** | 2 | 2 |

Where it collapsed, `frame` became 0 for a pure-insertion compound — and 0 is in neither
`ins_frame` (3n+1) nor `del_frame` (3n+2), **so those states were dropped in silence.**
That is a materially worse defect than this issue describes: not a wrong length, a lost
call.

### And correcting the lost-call example this branch had been carrying

`I9_2_A_LEN3&D50_2&D51_2` does **not** lose a call. It is dropped by the deletion half but
**kept** by the insertion half — frame |3 − 2| = 1, which is 3n+1 — so the reported output
is unchanged and the row merely moves halves. The state that genuinely stops being
reported is **`I9_2_A_LEN2&D50_2&D51_2`**, whose net length is 0, which is therefore in
frame and dropped by both halves. Both are now asserted.

**What now pins it**

- `tests/unit/test_advntr_output_parsing.py::TestSummedInsertionLength`
- `tests/unit/test_advntr_output_parsing.py::TestInsertionLenAndFrameAfterFiltering`
- `scripts/advntr_len_differential.py` — the committed differential sweep

**GOLDEN-COHORT: provably cannot attest this, and the reason is measured.** The cohort's
only compound state, dfc3's `D17_2&D18_2&D19_2&D20_2&D21_2`, carries no `LEN` token, so
`Insertion_len` is 0 under both semantics — confirmed byte-identical on both sides. The
evidence is the differential sweep instead:

```
52,511 probes
13,563 / 13,563 previously-parsing inputs byte-identical
38,943 differing, every one predicted by the oracle
0 violations in the three unchanged classes (no LEN, single terminal LEN,
  LEN only in the last part)
reported-call delta +16,659 / -1,649 across the swept grammar
```

The oracle was stated **before** the sweep ran: a state differs iff material follows its
first `LEN` token. Induced-failure evidence for the rows that pass unchanged: tightening
`LEN(\d+)` to `^LEN(\d+)` produced 21 failures spanning both changed and unchanged rows
and exited the sweep non-zero.

Output changes in both directions. The insertion path is monotone-additive — affected
states were 0 and always dropped, so summation can only make real compound calls appear.
The deletion path can lose a call.

**Left open and posted here for @hassansaei:** `Deletion_length` remains
`Variant.str.count("D")` (`advntr_genotyping.py:177`, `:216`), a count of deletion
**events** used arithmetically as a length in bp. That equals "the sum of deleted lengths"
only if every adVNTR `D` event is a single base. Consistent with the State grammar and
with the example above, but not confirmed. Filed as its own issue rather than assumed.

---

## #194 — Gate tests/ under mypy

Closed by **`65e033c`** (`chore(types): gate docker/app under mypy and fix the annotations
it found`).

23 real errors → 0, with **zero `# type: ignore`, zero relaxed settings, and zero new or
widened overrides**. `mypy_path = "docker"` in `pyproject.toml` puts the `app` package on
the search path for both mypy runs, so the `from app import …` lines in `tests/unit/web/`
resolve to the real signatures either way.

**What was actually wrong, and it was not what it looked like.** The first resolution
spelled `ASGIApp = Callable[..., Awaitable[None]]` — dropping the wrapped app's parameter
types in *production* code so that dict-shaped test doubles would still type-check. Once
the doubles were widened to the contract ASGI and Starlette publish,
`docker/app/request_limits.py` ended as a **comment-only diff**: the defect was in the
doubles the whole way through, and no production type precision was spent at all.

Two genuine production over-declarations remain fixed at the source: `uploads.ByteSource`
(protocol) and `cohorts.CohortRecord` (TypedDict).

**What now pins it**

- `make type-check` runs `mypy vntyper/ docker/app/`; `make type-check-all` adds
  `mypy vntyper/ tests/`. Both are in `make check-all` and `make ci-local`, and CI's
  `typecheck` job runs `type-check-all`.
- The baseline was re-measured independently rather than taken on report: extracting
  `451a4b7` to a scratch tree gives 23 errors (18 in `test_request_size_limit.py`, 1 in
  `test_upload_limits.py`, 4 in `test_cohort_alias.py`); dropping only the new test file
  onto unmodified production code gives 5, exactly the two that `ByteSource` and
  `CohortRecord` address.

**Not done, and filed separately:** `scripts/` is in `RUFF_PATHS` but in **no** mypy
target, so ~3200 lines added by this branch are type-checked by nothing. Adding it
surfaces one pre-existing error in a file this branch does not touch
(`scripts/download_test_data.py:147`), so it is filed rather than folded in. Recorded as
AGENTS.md trap 16.

---

## #195 — motif_processing.py scores 30.9% on mutation testing

Closed by **`11e2300`** (`refactor(motif): extract the decision layer, merge the twin
preprocessors, and fix the column-wide dash gate`).

`motif_processing.py` held **46 of the repository's surviving mutants**, the largest gap
the sweep found. Three changes:

**1. The twin preprocessors are merged.** `preprocessing_insertion` and
`preprocessing_deletion` were identical except the `Variant` literal and held 12 mutants
between them, all column plumbing that no test observed because the existing tests assert
on merged output rather than on the column contract. Both keep their signatures and
delegate to `_preprocess_vcf_frame`. Five of the six per-twin mutants are killed with
pasted induced failures. The sixth, `axis=1` → `axis=2`, is genuinely equivalent — pandas
ignores `axis` when `columns=` is passed, verified on 2.2.2 — so **the inert argument is
deleted rather than the mutant excused.** A mutant that cannot be killed because the code
it mutates does nothing is a signal to remove the code.

**2. The pure decision layer moves to `motif_decisions.py`.** The merge and annotate I/O
stay behind. `motif_decisions.py` is registered in `scripts/mutation_test.py` TARGETS — a
module absent from TARGETS is not mutated at all, so extracting the hard decisions without
registering them would have raised `motif_processing.py`'s score for the worst possible
reason. **Both modules must be quoted together; a `motif_processing.py` figure alone is no
longer comparable to the 30.9% baseline.**

One interface deviation from the plan: `apply_gg_alt_rule` cannot have the specified
signature and stay byte-identical, because the GG gate is evaluated pre-exclusion while
the allowlist is applied post-exclusion. The gate is split out as `has_gg_alternate`.

**3. `motif_processing.py:332` is fixed per row.** `.max()` aggregated the dash check over
the **whole column**, so one malformed `Motifs` value set `motif_filter_pass=False` for
every row of the sample and the run reported Negative while exiting 0. This is W8 from the
#179 audit, never filed until this branch.

**The fix uncovered a second, opposite face the audit did not record**, and the plan's own
fixture did not reproduce the defect it was written for: with `["X-5", "BROKEN"]` **both
rows pass**, because suppression needs a value with two or more dashes. So the same gate
has a false-negative face (a two-dash row suppresses its siblings) and a false-positive
face — **a dashless row rides through on a sibling's dash and is reported as a passing
call on a motif that is not in the reference.** The second is the more dangerous direction
and nothing recorded it.

Malformed rows are now dropped from the decision set but kept in the output marked
`motif_filter_pass=False` with NA annotations, which is what every other rejected row
looks like: this stage marks and `filter_final_dataframe` filters (AGENTS.md trap 4).

**What now pins it**

- `tests/unit/test_motif_decisions.py` — including `TestTheEndToEndOracle` (hashed against
  the pre-extraction code and passed unmodified) and
  `TestMalformedMotifIdsAreContainedToTheirOwnRow`
- `tests/unit/test_motif_preprocessing.py` — the merged preprocessor's column contract

The oracle was initially too weak — a `<` → `<=` mutation survived it — so a POS-60
boundary test was added **before any code moved** and the hash re-derived.

**GOLDEN-COHORT: `<fill in — run at ec67fff>`.** The expected delta is none, since real
motif IDs carry exactly one dash — **but that is an expectation, not evidence, and this
comment says so deliberately.** The oracle is an equivalence proof for the extraction, not
for the dash-gate fix.

Noted, not fixed, and filed separately: the `POS_fasta` threshold rebase updates a column
nothing reads afterwards. Pinned by a test.

---

## #196 — Enable branch coverage

Closed by **`debee23`** (enable + floor 70 → 74), **`c6d967a`** (→ 79) and **`ec67fff`**
(→ 80).

`branch = true` is enabled in `[tool.coverage.run]`, and the floor is now the
**branch-inclusive** figure at every step. Progression, each value taken verbatim from
what `scripts/coverage_gate.py` printed and never from the rounded `TOTAL` column:

| Commit | Measured | Floor |
| --- | --- | --- |
| `debee23` | 74.22% branch-inclusive | 70 → 74 |
| `c6d967a` | 79.02% (after the `cohort_summary.py` split) | 74 → 79 |
| `ec67fff` | **80.24%** (after the gated phase) | 79 → **80** |

80.24% **meets the 80% project target** for the first time, on a strictly harder
measurement than the target was originally stated against. The target stays a warning
rather than becoming the floor: they measure the same number for different purposes, and
collapsing them would leave no headroom above the gate.

**What now pins it**

- `tests/unit/test_coverage_gate.py::test_branch_coverage_is_enabled` — proved by
  commenting out `branch = true`, observing the failure, and reverting.
- `test_read_floor_returns_the_configured_value` — the existing floor pin, which failed
  **naturally** once the floor moved (`assert 74.0 == 70.0`).
- `test_the_ratchet_names_the_new_floor_once_coverage_crosses_an_integer`,
  `test_the_gate_fails_below_the_floor`.

**A finding worth keeping, and the reason a config assertion earns its own test:** deleting
`branch = true` **raises** the reported total (statement-only is 80.77% against 80.24%
branch-inclusive), which clears the floor and passes every other gate while covering
strictly less. The ratchet therefore *cannot* catch that regression — the number moves the
wrong way — and `test_branch_coverage_is_enabled` is the only thing that can. The
measurement is in its docstring.

Also in `debee23`: `docs/development/mutation-testing.md` is **machine-generated** from
`scripts/mutation_test.py`, and `tests/unit/test_mutation_test.py` asserted the superseded
numbers against the *generator*. A correction written into the Markdown would have been
silently reverted by the next `make mutation-render` while the test kept passing. The
correction was ported into the generator instead, and the round trip is proved:
`make mutation-render` leaves the committed file byte-identical. Superseded figures are
still asserted, by a separate test that labels them historical, so the correction cannot
be quietly dropped later.

**GOLDEN-COHORT: not applicable** — no production code changed in any of these three
commits.

---

## #197 — Enabling duplicate_flagging raises KeyError: sort column "Motifs"

Closed by **`6e7cda2`** (`fix(config): order duplicate flagging by Depth_Score alone, per
#197`).

`duplicate_flagging.sort_by` reduced to `Depth_Score` alone, and `flagging.py` now sorts
with `kind="stable"` so tied rows keep input order. The config diff touches **exactly one
key**, leaving `enabled`, `flag_name` and `group_by` alone.

`duplicate_flagging` remains **disabled** in the shipped config, so no genotype can move.
The argument is structural rather than an assertion: `select_single_best_variant` runs
strictly after `add_flags`, and the flag this fix affects is not read before selection.

**What now pins it** (`tests/unit/test_flagging.py`)

- `test_duplicate_ordering_uses_depth_score_only_and_no_motif_column`
- `test_duplicate_flagging_stays_disabled_in_the_shipped_config`
- `test_the_flagged_row_is_deterministic_when_depth_scores_tie`

**Two corrections to the plan, recorded because both are easy to repeat.**

(a) The drafted determinism test **could not have failed**: pandas' quicksort only
reorders ties above ~260 rows on this build, so a 3-row case passes with or without
`kind="stable"`. Rewritten at 300 rows through `add_flags` with the toggle enabled, where
it fails naturally before the fix. Verified by reverting `kind="stable"` and running it
three times.

(b) The claim that `kestrel_genotyping.py:588` "only calls `add_flags` with a duplicates
config when that flag is true" is imprecise. Line 588 reads
`if flagging_rules or duplicates_config.get("enabled", False)`, and `flagging_rules` is
non-empty in the shipped config, so **`add_flags` is always called there**. The real
inertness boundary is one level deeper, inside `add_flags` at `flagging.py:150`. The "no
genotype can move" conclusion still holds, but for a different reason than the plan gave.

`docs/pipeline/flagging.md` records the sort key and why it is a single column.
