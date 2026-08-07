# SPEC — Milestone 2, "Correctness of reported numbers"

Branch: `fix/milestone-2-correctness-of-reported-numbers`
Baseline: `cb593b6` (tip of `main`)
Issues: #171, #172, #174, #203, #212

This document is the contract the implementation is written against. Every
behaviour change below is stated as *before → after*, with the artefact it moves, the
test that will catch it, and the golden-cohort diff it is expected to produce. A diff
the gate reports that is not predicted here is a defect in this spec, not a surprise to
be waved through.

## Scope decisions taken before writing code

Three forks were settled against the code rather than by preference. They are recorded
here because two of them contradict the issue bodies.

### D1 — #203 is settled by the reference FASTA: delete the rebase

The issue offers two readings and calls the choice the maintainer's. The reference
FASTA removes the choice.

```
$ awk '/^>/{if(n)print n,l; n=substr($0,2); l=0; next}{l+=length($0)}END{print n,l}' \
      reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa | head -3
1-2 120
2-3 120
3-4 120
```

Every record in the Kestrel reference is a **pair** of 60 bp motifs, 120 bp long, named
`<left>-<right>`. The chain is:

| Step | Location | Value |
| --- | --- | --- |
| Kestrel writes the VCF against that FASTA | — | `#CHROM` = a pair name, e.g. `X-5` |
| `#CHROM` is renamed | `motif_processing.py:104` | → `Motifs` |
| `Motif_fasta` is assigned | `motif_processing.py:371` | `Motif_fasta = Motifs`, **verbatim** |
| `Motif` is assigned | `motif_processing.py:418`, `:446` | the *half*-motif name (`2`, `X`) |
| BED is written | `kestrel_genotyping.py:632-636` | `{Motif_fasta}\t{POS_fasta}\t{POS_fasta+1}` |

So the name written into column 1 of `output.bed` is the **120 bp pair record**, and
`POS` from the VCF is already a coordinate in that record's space. `POS_fasta == POS` is
correct. Subtracting `position_threshold` (60) would make every row at or above the
threshold wrong by exactly 60 bp against the FASTA the BED names.

The issue body's premise — "`Motif_fasta` is a repeat-unit identifier, so the coordinate
beside it is being interpreted in that coordinate space" — is true of `Motif`, not of
`Motif_fasta`. **Decision: the code is right, the comment is wrong. Delete the dead
rebase, correct the comment, rename the test whose name asserts the opposite.**

Confirmed with the user, 2026-08-07.

### D2 — #171 lands the arithmetic; the region harmonisation is deferred

`#171` bundles two independent changes: the depth arithmetic, and harmonising the
hg19 (1501 bp) and hg38 (4501 bp) region spans. The second moves every hg38 coverage
number ever reported and requires naming a flank convention — a MUC1-domain decision,
and the audit comment on the issue recommends splitting them so each gate diff is
attributable to one change.

**Decision: this PR lands the arithmetic only.** #171 is closed with a comment stating
plainly which of its four suggested items are done and which is not, and a follow-up
issue is filed for the harmonisation. Confirmed with the user, 2026-08-07.

### D3 — the `output.bed` off-by-one is an in-scope, separately-committed rider to #203

`generate_bed_file` writes `{POS}\t{POS+1}`. BED intervals are 0-based half-open; a
1-based VCF `POS` of *p* denotes the interval `[p-1, p)`. The current output names the
base *after* the variant. #203 rewrites this area, and the audit comment notes that no
test covers `generate_bed_file`'s output content at all — so the new test this spec
requires would otherwise have to pin a known-wrong coordinate as if it were intended.

**Decision: fix it, in its own commit, labelled as a rider.** `output.bed` is not a
golden-cohort-compared artefact (`scripts/golden_cohort/compare.py:39-50` lists
`kestrel_result`, `kestrel_pre_result`, `coverage_summary`, … — no BED), so this cannot
move the gate. Confirmed with the user, 2026-08-07.

### D4 — the rule table is 40 rules, not 44

The audit comment on #172 says `report_config.json` carries 44 screening rules. It
carries **40**: 5 `kestrel_result` values × 4 `advntr_result` values × 2
`quality_metrics_pass` values. The cost argument is unchanged — a fourth axis doubles it
to 80 — but the number quoted in the plan is the measured one.

---

## #171 — `samtools depth` without `-a`, and a mean over covered bases

**Root issue. #172 depends on it. Genotype-affecting.**

### The trap that makes this one change and not two

Adding `-a` in isolation silently destroys `percent_uncovered`. It is currently derived
by *subtraction*:

```python
# coverage_stats.py:163-168
covered_bases_count = len(coverage_values)
zero_coverage_bases = total_region_length - covered_bases_count
percent_uncovered = 0 if total_region_length <= 0 else zero_coverage_bases / total_region_length * 100
```

With `-a`, `samtools depth` emits one row per region position including zeros, so
`len(coverage_values) == total_region_length` and `zero_coverage_bases` becomes **0 for
every sample**. The one coverage metric that is correct today would become the one that
is always wrong. The divisor fix and the zero-counting fix must land in the same commit.

### Behaviour change

| | Before | After |
| --- | --- | --- |
| depth command | `samtools depth -@ N -r REGION BAM > OUT` | `samtools depth -a -@ N -r REGION BAM > OUT` |
| `coverage_values` | covered positions only | every position in the region, zeros included |
| `mean` | `sum(v) / len(v)` — divides by covered bases | `sum(base) / len(base)` where `base` spans the region |
| `median`, `stdev`, `min`, `max` | over covered positions | over the region; `min` is 0 wherever anything is uncovered |
| `uncovered_bases` | `region_length - len(v)` | count of zero-depth positions in `base` |
| `percent_uncovered` | unchanged in value | unchanged in value |
| `region_length` | from the region string | from the region string (unchanged) |

`base` is defined once, and defends both the `-a` file and a legacy file:

```python
absent = max(total_region_length - len(coverage_values), 0)
base = coverage_values + [0] * absent
```

* Production (`-a`): `absent == 0`, `base is` the region, zeros are already present.
* A legacy no-`-a` file, or a region truncated at a contig end: the missing positions are
  restored as zeros, and every statistic is over the region regardless.
* `total_region_length <= 0` (unparseable region): `absent == 0`, `base ==
  coverage_values`, `percent_uncovered` stays 0 and the WARNING remains the only signal.

`summarise_coverage([])` still raises `RuntimeError`. Under `-a` an empty depth table can
now only mean the region matched no contig at all, which sharpens the guard rather than
weakening it.

### The free regression test

For any legacy depth list with no zero rows and `region_length > 0`:

```
mean_new == mean_old · (1 − percent_uncovered_old / 100)
```

because `(S/c) · (c/T) = S/T`. This is asserted directly as a unit test — it needs no
fixture, and it is the identity the issue offers for reconciling historical output.
Guarded on `region_length > 0`, per the audit note.

### Both claims confirmed empirically on a shipped test BAM

```
$ B=tests/data/example_40cf_hg38_subset.bam; R=chr1:155188000-155192500   # span 4501
$ samtools depth    -r $R $B | wc -l                          # 4047   (covered only)
$ samtools depth -a -r $R $B | wc -l                          # 4501   (the whole region)
$ samtools depth -a -r $R $B | awk '$3==0' | wc -l            #  454   (zero-depth rows)
$ samtools depth    -r $R $B | awk '{s+=$3;n++}END{print s/n}'  # 164.5325  old mean
$ samtools depth -a -r $R $B | awk '{s+=$3;n++}END{print s/n}'  # 147.9367  new mean
```

* The subtraction trap is real: `region_length − len(values)` is `4501 − 4047 = 454`
  today and would be `4501 − 4501 = 0` after `-a`. The zero-count is 454 — exactly the
  same number, which is why counting zeros is the correct replacement and not merely a
  different one.
* The identity holds exactly: `164.5325 · (1 − 454/4501) = 147.9366`.
* This one sample's reported mean falls by 10.1%.

### Artefact delta

* `<prefix>_summary.tsv` — `mean`, `median`, `stdev`, `min` change for every sample with
  any uncovered base. `max`, `region_length`, `uncovered_bases`, `percent_uncovered` are
  numerically unchanged.
* `<prefix>_vntr_coverage.txt` — grows to one row per region base (1501 or 4501 rows).
* `pipeline_summary.json`, `Coverage Calculation` step — same values as the TSV; the
  `command` string gains `-a`.
* **Genotype-affecting path**: `downsample_bam_if_needed`
  (`fastq_bam_processing.py:389-404`) reads `["mean"]` and computes
  `fraction = max_coverage / current_coverage`. The corrected mean is **≤** the old mean,
  so fewer samples cross `--advntr-max-coverage 300`, and those that do are downsampled
  *less*. adVNTR then genotypes a different read set.

  The audit names the golden-cohort gate (`scripts/golden_cohort_gate.py:131`) as the
  configuration that passes 300. It is not the only one: **`docker/app/tasks.py:215`
  passes `--advntr-max-coverage 300` too**, so the deployed web service runs in the
  genotype-affecting configuration as well. The CLI default remains `None`, which is
  unaffected.

### Tests to add

| Test | Asserts |
| --- | --- |
| `test_the_mean_is_over_the_region_not_over_covered_positions` | `summarise_coverage([20]*10, 100)["mean"] == 2.0` |
| `test_uncovered_bases_counts_zero_depth_rows_when_the_file_is_complete` | `summarise_coverage([10,0,0,10], 4)` → `uncovered_bases == 2`, `percent_uncovered == 50.0` |
| `test_the_minimum_is_zero_wherever_the_region_is_uncovered` | `min == 0` for the above |
| `test_a_legacy_depth_file_is_padded_to_the_region` | `summarise_coverage([20]*10, 100)` → `uncovered_bases == 90` (identical to today) |
| `test_the_closed_form_identity_reconciles_old_and_new_means` | the identity above, parametrised over several (values, T) |
| `test_the_depth_command_requests_every_position` | `-a` present, and `shlex.split` order pinned |
| `test_a_region_longer_than_the_depth_file_still_reports_the_region_length` | truncated-region path |
| `test_a_zero_length_region_reports_zero_uncovered_rather_than_a_negative_count` | replaces the `-2` pin below |

### Tests that must change (deliberately, never silenced)

| File:line | Today | After |
| --- | --- | --- |
| `test_coverage_stats.py:136` | `assert stats["mean"] == 20.0, "the mean is over covered positions"` | `== 2.0`, docstring rewritten to state the region definition |
| `test_coverage_stats.py:127` | `mean == 10.0` on a fully covered region | unchanged in value; docstring clarified |
| `test_coverage_stats.py:152` | single position, `mean == 42.0`, `T=1` | unchanged |
| `test_coverage_stats.py:173` | `uncovered_bases == -2` for `T=0` | `== 0`; a negative count was nonsense, and `absent` is now floored at 0 |
| `test_coverage_stats.py:182` | even-sample median over 4 values, `T=4` | unchanged |
| `test_command_builders.py:241` | pins the depth command string without `-a` | pins it **with** `-a` |
| `test_fastq_bam_command_wiring.py:378` | pins exact TSV bytes `20.00\t20.00\t10.00\t10\t30\t1501\t1498\t99.80` | recomputed under the new arithmetic **and** the #172 column |
| `coverage_stats.py:20-26` module docstring | states covered-position semantics as intentional | states region semantics, and why |

### Predicted golden-cohort diff

* `coverage_summary` — **delta on every case with any uncovered base.** Four columns move.
* `executed_commands` — delta on every case: the depth command gains `-a`. This is
  `commands.jsonl`, the shell commands the run actually executed
  (`scripts/golden_cohort/artifacts.py:252-259`).
* `pipeline_step_records` — **no delta.** Checked, because the obvious prediction is
  wrong twice over: the `Coverage Calculation` record's `command` is the literal string
  `"calculate_vntr_coverage(...)"` (`pipeline.py:489`), not the depth command, so `-a`
  never reaches it; and `coverage_summary.tsv` is in `DIRECTLY_COMPARED_RESULT_FILES`
  (`scripts/golden_cohort/normalise.py:30-36`), so `strip_step_record` drops its `md5sum`
  before comparison. The same holds for `kestrel_result.tsv` under #174.
* `advntr_result` — **possible genotype delta** on cases whose old mean exceeded 300 and
  whose new mean does not, or whose downsampling fraction changes. This is the diff the
  audit warned is not presentational. Each such case must be attributed by name, with the
  old and new mean quoted, before the gate is called green.
* `screening_summary` — delta only where a corrected mean crosses `mean_vntr_coverage:
  100`. The audit measured 61 such samples in a 8215-row cohort; the local 58-case matrix
  may show none.

---

## #172 — QC gates on mean coverage only

**Depends on #171.** Gating on the biased mean would freeze the bias into the rule.

### Shape: widen `quality_metrics_pass`, do not add a fourth axis

```python
# screening_summary.py:272, before
quality_metrics_pass = not (mean_vntr_coverage is not None and mean_vntr_coverage < mean_vntr_cov_threshold)
```

`report_config.json` enumerates **40** screening rules as
`kestrel_result × advntr_result × quality_metrics_pass`, and
`tests/unit/test_screening_summary.py:230/259/272` assert that matrix is complete,
message-covered and non-overlapping. A fourth axis makes it 80 rules and rewrites those
three tests. Widening the *meaning* of the existing axis leaves all 40 rules and all
three tests untouched, and — decisive for a research-use tool whose clinical-sounding text
is config-driven — **requires no new or reworded message**. The 40 shipped messages
already say "quality metrics are below threshold" without naming which metric.

### New pure module: `vntyper/scripts/coverage_qc.py`

The QC verdict is a decision, and AGENTS.md rule 3 puts decisions in a focused pure
module rather than back into `generate_report.py` (574 LOC) or `fastq_bam_processing.py`
(612 LOC).

```python
COVERAGE_QC_PASS = "PASS"
COVERAGE_QC_FAIL = "FAIL"

@dataclass(frozen=True)
class CoverageQC:
    passed: bool
    status: str                 # COVERAGE_QC_PASS | COVERAGE_QC_FAIL
    reasons: tuple[str, ...]    # ("mean_vntr_coverage", "percent_vntr_uncovered"), in that order

def evaluate_coverage_qc(
    mean_vntr_coverage: float | None,
    percent_vntr_uncovered: float | None,
    mean_threshold: float,
    percent_threshold: float,
) -> CoverageQC: ...
```

Rule, stated exactly:

* fails on mean when `mean is not None and mean < mean_threshold`
* fails on uncovered when `pct is not None and pct > percent_threshold`
* `passed` is the conjunction of the two negations; a `None` metric never fails.

The `None`-passes rule is not new — it is the existing `quality_metrics_pass` semantics
(`test_screening_summary.py:388` pins it: "an unmeasured sample is reported as passing")
and the existing `MISSING_AS_OK` icon semantics. It is preserved so that this issue
changes exactly one thing.

Note the asymmetry, which is deliberate and matches `threshold_icon`'s
`higher_better` argument: mean uses `<` (strictly below fails), percent uses `>`
(strictly above fails). A sample at exactly 50.0% uncovered passes.

### Emitting `coverage_qc`

The issue asks for the outcome as an explicit status in the result output, not only as a
template icon. It is written into the coverage summary, which is the one place that
reaches **both** consumers without new plumbing:

* `COVERAGE_COLUMNS` gains a ninth entry, `"coverage_qc"`, so
  `format_coverage_summary` writes it and the TSV carries it.
* `COVERAGE_FIELD_TYPES` (`report_formatting.py:117-126`) gains `"coverage_qc": str`, so
  `parse_coverage_stats` reads it back.
* The cohort path needs **no whitelist edit**: `additional_stats_frame`
  (`cohort_tables.py:189-193`) flattens the whole coverage mapping with a `cov_` prefix
  and applies no projection, so the column appears as `cov_coverage_qc` in the cohort
  statistics table automatically. The audit's warning about `KESTREL_DISPLAY_COLUMNS`
  applies to the *Kestrel* table, which is not where a sample-level coverage verdict
  belongs. **`KESTREL_DISPLAY_COLUMNS` is not modified.**

`summarise_coverage` keeps its two-argument signature and stays free of thresholds.
`calculate_vntr_coverage` reads the two thresholds from `config["thresholds"]`, calls
`evaluate_coverage_qc`, and merges `status` into the stats dict before writing.

Thresholds are read with `.get` defaults matching `config.json` (`100` and `50.0`), because
`--config-path` replaces the whole config rather than merging (AGENTS.md trap 2) and a
`KeyError` deep in the coverage stage is a worse failure than a documented default.

### Wiring into the report

`build_screening_summary`'s signature changes from
`(kestrel_df, advntr_df, advntr_available, mean, mean_threshold, report_config)` to
`(kestrel_df, advntr_df, advntr_available, coverage_qc: CoverageQC, report_config)`, and
sets `quality_metrics_pass = coverage_qc.passed`.

Passing the evaluated verdict rather than four raw numbers means the axis and the emitted
column cannot disagree — they are the same value. `generate_report.py` calls
`evaluate_coverage_qc` once, from the `mean`/`percent_uncovered` it already parses at
`:412-414` plus the two thresholds it already loads, and uses the result for the screening
axis, the `coverage_qc` context key and the existing icons.

It **recomputes** rather than trusting the TSV's `coverage_qc` string, so that
`vntyper report` still works against a `pipeline_summary.json` written before this change.
Same pure function, same inputs, so the two can only agree.

### Artefact delta

* `<prefix>_summary.tsv` — one new trailing column, `coverage_qc`, `PASS` or `FAIL`.
* Cohort statistics table and cohort exports — new `cov_coverage_qc` column.
* HTML report — a new `Coverage QC` row appended to the coverage table
  (`vntyper/templates/report_template.html:211-244`, after `Percentage Uncovered`);
  `screening_summary` text changes for any sample that passed on mean and fails on
  uncovered fraction.

### Tests to add

| Test | Asserts |
| --- | --- |
| `test_a_sample_above_both_thresholds_passes` | `PASS`, empty reasons |
| `test_a_low_mean_fails_and_names_the_mean` | `FAIL`, `reasons == ("mean_vntr_coverage",)` |
| `test_a_patchy_vntr_fails_even_with_acceptable_mean` | the issue's headline case: mean 250, 80% uncovered → `FAIL` |
| `test_both_failures_are_reported_in_declaration_order` | `reasons == ("mean_vntr_coverage", "percent_vntr_uncovered")` |
| `test_a_metric_that_was_never_measured_does_not_fail_the_gate` | `None` inputs → `PASS` |
| `test_the_boundaries_are_inclusive_on_pass` | mean == threshold passes; pct == threshold passes |
| `test_the_coverage_summary_tsv_carries_the_qc_verdict` | ninth column present, correct value |
| `test_a_patchy_sample_fails_the_screening_quality_axis` | through `build_screening_summary` |
| `test_the_cohort_statistics_table_shows_cov_coverage_qc` | end-to-end through `additional_stats_frame` |

### Tests that must change

The audit named three of these. Enumerated against the code, there are **six files and
ten `build_screening_summary` call sites**, not three:

| File:line | Change |
| --- | --- |
| `test_screening_summary.py:329, 351, 357, 364, 371, 378, 384, 392, 404, 417` | **ten** call sites take the current 6-argument signature, not the two the audit lists. All become `CoverageQC` calls. `:392`'s "an unmeasured sample passes" pin is preserved in meaning, not deleted. `:404`/`:417` drive the exception path with an `Exploding()` config and must keep doing so. |
| `test_fastq_bam_command_wiring.py:378` | exact TSV bytes gain a ninth field (also changed by #171) |
| **`tests/helpers.py:405`** | carries its **own** copy of `COVERAGE_COLUMNS`, and `test_helpers.py:83` asserts it equals production's tuple. It must gain `coverage_qc` or the integration tier's validator fails. Missed by the audit. |
| `tests/unit/test_helpers.py:112` | parametrises a dropped-column test over that tuple — gains a ninth case automatically |
| `tests/unit/test_coverage_stats.py:60` | asserts `COVERAGE_COLUMNS == (...)` as an exact literal tuple |
| `tests/unit/test_report_formatting.py:71, 87, 96, 102` | `set(COVERAGE_FIELD_TYPES) == set(COVERAGE_COLUMNS)` and three round-trip assertions |

`validate_coverage_output` only coerces `mean`, `median` and `percent_uncovered` to
float (`tests/helpers.py:452-460`), so a **string** ninth column passes through it
untouched. Checked, because a `float("PASS")` would have failed the whole integration
tier.

Two further consumers were checked and need no change: `tests/support/pipeline_harness.py:238`
mocks `calculate_vntr_coverage` to `{"mean": 100.0}` and never touches the schema, and
`tests/benchmark/benchamrk_downsample.py:172` carries a private copy of the function that
no gate runs.

### Predicted golden-cohort diff

* `coverage_summary` — one new column on every case (on top of #171's value changes).
* `screening_summary` — delta only for samples that pass on mean and fail on uncovered
  fraction. Expected to be rare in the 58-case matrix, which uses subset BAMs over the
  VNTR; if it fires, the case must be named with both metrics quoted.
* `cohort_kestrel_*` — unchanged. `cohort_tables` — new `cov_coverage_qc` column.
* **`is_positive` cannot move.** `screening_summary.py:297-299` derives it as
  `is_finding(kestrel_result, …) or is_finding(advntr_result, …)` — the quality axis is
  not an operand. So widening the QC rule changes the sentence and the report's QC row,
  never whether a sample reads as a finding. Stated because it is the one way this issue
  could have become genotype-affecting, and it is not.

---

## #174 — 4 bp insertion artifacts reported as positive calls

**Independent. Genotype-affecting.**

### The severity, restated from the audit

`report_config.json` maps `Confidence ∈ {High_Precision, High_Precision*}` +
`Flag != "Not flagged"` → `High_Precision_flagged`, and `is_finding`
(`screening_summary.py:201`) returns True for anything that is not the block's `default`.
So an artifact-only sample sets `is_positive = True` and the report styles a known
technical artifact as a positive MUC1 finding.

### Shape: a config-driven artifact gate, not `Flag` in `filter_cols`

Adding `"Flag"` to `filter_cols` is one line and the wrong shape, for two reasons:

1. The gate contract (`kestrel_genotyping.py:848-863`) *raises* when a required column is
   absent, and a negative run's frame legitimately carries no `Flag`
   (`report_formatting.py:68`, `cohort_tables.py:41-42`).
2. It would exclude **every** flag, including `Low_Depth_Conserved_Motifs`, which is an
   advisory annotation rather than an artifact. Suppressing those rows would delete real
   low-depth calls — the exact failure mode this milestone exists to prevent.

Instead: `kestrel_config.json` declares which flags are artifacts, and a derived boolean
gate column is written unconditionally.

```json
"artifact_flags": ["False_Positive_4bp_Insertion"]
```

```python
# flagging.py — new pure function
def add_artifact_gate(df: pd.DataFrame, artifact_flags: Sequence[str]) -> pd.DataFrame:
    """Write ``flag_filter_pass``: False iff the row carries a declared artifact flag."""
```

* Always writes the column, including when `Flag` is absent (then all True) and when the
  frame is empty. So the gate contract can never raise for a missing `flag_filter_pass`.
* `Flag` is a comma-joined string, so membership is tested per element after splitting on
  `", "`, not by substring — a flag named `X` must not match a flag named `XY`.
* `Flag` is coerced with `.fillna("")` before splitting. `add_flags` always writes a
  string (`flagging.py:145`), but the column is also written by `mark_potential_duplicates`
  and can be reindexed by upstream merges, so a NaN must degrade to "no artifact" rather
  than raising inside the gate.
* Called in `process_kmer_results` step 6.5, **unconditionally**, after the existing
  conditional `add_flags` block. Placing it outside that `if` is what guarantees presence.
* `filter_cols` gains `"flag_filter_pass"` → six gates.

Rows are excluded at the filter, which is the repo's established contract ("stages mark,
they do not filter", AGENTS.md trap 4). The evidence survives in
`kestrel_pre_result.tsv`, which carries every row with `flag_filter_pass=False` — so
#131's symptom, a flagged variant vanishing without trace, is answered by the pre-result
rather than by improvising a representation in each output path.

An artifact-only sample therefore reports as negative. That is the correct answer: the
4 bp insertion is a known technical artifact, so the sample has no call.

The `*_flagged` screening states stay reachable through `Low_Depth_Conserved_Motifs` and
`Potential_Duplicate`, so **`report_config.json` and the 40 rules are unchanged**, and
`test_screening_summary.py:189` (asserting `High_Precision_flagged` is a finding) stays
true and stays as-is.

### Artefact delta

* `kestrel_result.tsv` — artifact rows removed. An artifact-only sample yields an empty
  result where it previously yielded one flagged row.
* `kestrel_pre_result.tsv` — one new column, `flag_filter_pass`, on every row.
* HTML report and cohort — an artifact-only sample flips from positive to negative.

### Tests to add

| Test | Asserts |
| --- | --- |
| `test_an_artifact_only_sample_is_not_reported_as_a_call` | the regression test the issue asks for: a `REF=C/ALT=CGGCA` frame through `process_kmer_results` → empty result |
| `test_the_artifact_row_survives_in_the_pre_result_with_its_gate_false` | evidence is not destroyed |
| `test_an_advisory_flag_does_not_exclude_the_row` | `Low_Depth_Conserved_Motifs` still passes the gate |
| `test_a_frame_without_a_flag_column_still_gains_the_gate` | column present, all True |
| `test_an_empty_frame_gains_the_gate_column` | no raise downstream |
| `test_a_flag_name_is_matched_whole_not_as_a_substring` | `"XY"` must not be excluded by an artifact named `"X"` |
| `test_a_row_carrying_both_an_artifact_and_an_advisory_flag_is_excluded` | comma-joined parsing |
| `test_the_shipped_config_declares_the_4bp_insertion_as_an_artifact` | ties the code to `kestrel_config.json` |

### Tests that must change

| File:line | Change |
| --- | --- |
| `test_kestrel_filtering.py:72` | reads `kestrel_genotyping.py` **as source text** and asserts exactly 5 gate columns. Changed **once**, to 6, with `GATE_COLUMNS` extended. The audit notes this tripwire is shared with #173; #173 is not in this milestone, so it is changed here and #173 will find it already correct. |
| `test_haplo_count_and_selection.py:368` (`test_all_flagged_selects_best`) | tests `select_single_best_variant` directly, which still deprioritises rather than excludes. Survives unchanged; a docstring note records that exclusion now happens upstream at the gate. |
| `test_haplo_count_and_selection.py:427` (`test_issue_145_scenario`) | survives; extended with an artifact-flagged competitor to pin the new interaction, per the audit. |

### Predicted golden-cohort diff

* `kestrel_pre_result` — **delta on every case**: one new column.
* `kestrel_result` — delta only on cases whose selected variant is a 4 bp insertion
  artifact. Each must be named.
* `screening_summary`, `report_tables`, `cohort_*` — follow from the above.
* `pipeline_step_records` — **no delta**, for the reason given under #171:
  `kestrel_result.tsv` is directly compared, so its `md5sum` is stripped.

---

## #203 — the `POS_fasta` rebase writes a column nothing reads

**Independent. Not genotype-affecting. Zero predicted gate diff.**

### Behaviour change

Per D1, `POS_fasta` is already correct and the rebase is both dead and semantically
wrong. `motif_processing.py:490-498` becomes:

```python
# POS is a 1-based offset into the 120 bp pair record named by Motif_fasta
# (= the VCF #CHROM, = a record of All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa),
# so it is already the coordinate output.bed needs. No rebase. See #203.
combined_df["POS"] = pd.to_numeric(combined_df["POS"], errors="coerce").fillna(-1).astype(int)
combined_df["POS_fasta"] = combined_df["POS"]
```

`position_threshold` keeps its other, live use: `split_left_right` still uses it to decide
which half of the pair name a row's `Motif` is. Only the subtraction goes.

### Rider (D3): the BED coordinate

```python
# kestrel_genotyping.py:632-636
- bed_file.write(f"{motif_fasta}\t{pos}\t{pos + 1}\n")
+ bed_file.write(f"{motif_fasta}\t{pos - 1}\t{pos}\n")
```

Separate commit, so the gate can attribute it independently even though `output.bed` is
not a compared artefact.

### Artefact delta

* `kestrel_result.tsv` / `kestrel_pre_result.tsv` — **none.** `POS_fasta` values are
  byte-identical; only dead code is removed.
* `output.bed` — start and end each shift down by 1 (the rider). IGV highlights the
  variant base rather than the one after it.

### Tests to add

| Test | Asserts |
| --- | --- |
| `test_generate_bed_file_writes_zero_based_half_open_intervals` | the first content test this function has ever had: `POS 67` → `X-5\t66\t67` |
| `test_generate_bed_file_names_the_pair_record_not_the_half_motif` | column 1 is `Motif_fasta`, and it is a record ID present in the shipped reference FASTA |
| `test_generate_bed_file_returns_none_without_the_two_columns` | existing guard, now covered |
| `test_generate_bed_file_returns_none_for_an_empty_frame` | existing guard, now covered |
| `test_pos_fasta_is_the_vcf_position_because_motif_fasta_names_the_120bp_pair` | replaces the mis-named pin; the name states *why* |

### Tests that must change

| File:line | Change |
| --- | --- |
| `test_motif_decisions.py:140` (`test_the_original_pos_column_is_never_rebased_only_pos_fasta_is_derived`) | The name asserts the opposite of the truth: `POS` was "never rebased" only because the rebase landed on a discarded copy, and `POS_fasta` is a verbatim copy, not "derived". Renamed to `test_pos_and_pos_fasta_are_both_the_vcf_position_in_the_pair_record`, docstring rewritten to state D1's reasoning so the next reader cannot re-derive the wrong conclusion. Assertions on values are unchanged — they were already right. |
| `test_motif_decisions.py:125` (`POS_fasta` == 67 assertion) | unchanged in value; docstring gains the reason. |

### Predicted golden-cohort diff

**None.** If the gate reports a `kestrel_result` or `kestrel_pre_result` delta
attributable to this change, the change is wrong and must be reverted.

---

## #212 — Kestrel skipped when `output.vcf` exists

**Independent. Clinical-safety.**

### Behaviour change

```python
# kestrel_genotyping.py:294-297, deleted
if vcf_path.is_file():
    logger.info("VCF file already exists, skipping Kestrel run...")
    return
```

Replaced by: log a WARNING naming the stale file, unlink it, and run Kestrel. The skip is
unconditional, ungated, undocumented and untested today; #20's `--resume` is the place for
deliberate reuse, and is out of scope here. Deleting the ad-hoc version now is what #212
asks for and is a precondition for #20 rather than a competitor to it.

Unlinking rather than relying on the JAR to overwrite is deliberate: if Kestrel 1.0.1
refused to overwrite, the loop would fail and the run would abort — correct, but the
warning would be about Java rather than about a stale artefact. Removing it first makes
the reason legible in the log.

### The second half — a missing result must be loud

Two guards, at two levels:

1. **At the source.** `run_kestrel` tracks whether any k-mer size completed
   post-processing. If none did, it raises `RuntimeError`. Today the loop can fall out
   silently: `break` is only reached inside `if vcf_path.is_file()`, so a Kestrel run
   that exits 0 without writing a VCF leaves `run_kestrel` returning None. This closes
   outcome 2 — "a confident negative is manufactured" — where it originates.

2. **At the recorder.** `record_step` logs an ERROR and sets `record["result_file_missing"]
   = True` when `result_file` does not exist. **Set only when absent**, so a normal run's
   `pipeline_summary.json` is byte-identical and the gate's `pipeline_step_records`
   artefact does not diff.

This deliberately does **not** add a fourth screening axis for "result missing". The
`RuntimeError` means the pipeline no longer reaches the report in that state, so a state
for it would be unreachable — and the audit's cost warning on #172 applies equally here.

### Artefact delta

* None on a clean run. Both changes are on paths a fresh output directory never takes.
* On a re-run into a populated directory: Kestrel now runs (previously skipped), so
  `output.vcf`, `output.bam` and `kestrel_result.tsv` are regenerated.

### Tests to add

| Test | Asserts |
| --- | --- |
| `test_a_stale_vcf_does_not_skip_the_kestrel_run` | pre-create `output.vcf`, run the stage, assert the Kestrel command was executed — the regression test the issue asks for |
| `test_a_stale_vcf_is_removed_before_the_run` | the file is unlinked, and a WARNING names it |
| `test_a_run_that_produces_no_vcf_raises_rather_than_returning_silently` | `RuntimeError`, not `None` |
| `test_post_processing_runs_after_a_successful_kestrel_run` | `process_kestrel_output` is reached — the two statements the old `return` skipped |
| `test_record_step_flags_a_missing_result_file` | `result_file_missing is True`, ERROR logged |
| `test_record_step_does_not_add_the_flag_when_the_file_exists` | the key is absent — pins the gate-neutrality above |

`grep -rn "VCF file already exists" tests/ docs/` returns nothing today, so there is no
test to update — only tests to add.

### Tests that must change

None.

### Predicted golden-cohort diff

**None.** Every gate case runs into a fresh output directory.

---

## Cross-cutting

### Execution lanes

| Lane | Issues | Depends on |
| --- | --- | --- |
| A | #171 → #172 | sequential; #172 gates on #171's corrected mean |
| B | #174 | — |
| C | #203 (+ BED rider) | — |
| D | #212 | — |

B, C, D are independent of A and of each other. The only file two lanes both touch is
`kestrel_genotyping.py` (B: `filter_cols` and step 6.5; C: `generate_bed_file`; D:
`run_kestrel`), in three disjoint regions.

### Files touched

| File | #171 | #172 | #174 | #203 | #212 |
| --- | --- | --- | --- | --- | --- |
| `command_builders.py` | ✓ | | | | |
| `coverage_stats.py` | ✓ | ✓ | | | |
| `coverage_qc.py` (new) | | ✓ | | | |
| `fastq_bam_processing.py` | | ✓ | | | |
| `report_formatting.py` | | ✓ | | | |
| `screening_summary.py` | | ✓ | | | |
| `generate_report.py` | | ✓ | | | |
| `report_template.html` | | ✓ | | | |
| `config.json` | | (read only) | | | |
| `kestrel_config.json` | | | ✓ | | |
| `flagging.py` | | | ✓ | | |
| `kestrel_genotyping.py` | | | ✓ | ✓ | ✓ |
| `motif_processing.py` | | | | ✓ | |
| `summary.py` | | | | | ✓ |

### Docs and changelog

* `docs/about/changelog.md` — one entry per issue under a new version heading.
* `AGENTS.md` — trap 4 ("stages mark, they do not filter") gains the sixth gate column;
  the pinned-test list gains nothing (the traps section does not enumerate them).
* `docs/development/golden-cohort-gate.md` — a new run row and a result section
  attributing every delta.
* No new page, so no `mkdocs.yml` `nav:` edit is needed. Verified before claiming.

### Version

`vntyper/version.py`, `CITATION.cff` and `docs/about/changelog.md` must agree
(AGENTS.md trap 12). A minor bump to `2.1.0` is proposed: #171, #172 and #174 change
reported values and the QC verdict, which is more than a patch. **No `v*.*.*` tag is
pushed** — tagging publishes to PyPI irreversibly (AGENTS.md "Never", and the user's
standing instruction).

### Definition of done

1. `make check-all` green, output pasted.
2. Golden-cohort gate re-run, every delta attributed to a named change from this spec.
3. One PR, not stacked, closing #171, #172, #174, #203, #212.
4. A comment on #171 stating which of its four items are done and which is deferred, and
   a follow-up issue for the region harmonisation.
5. Changelog and docs updated.
