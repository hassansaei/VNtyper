# Golden-cohort gate (#179)

Before-versus-after comparison of genotyping output across the whole local test cohort,
run to decide whether the deliberately behaviour-changing commits in #179 may ship.

**Verdict: PASS.** Every genotype field, every `Confidence` label and every `Flag` is
byte-identical between the two commits, on every sample and every assembly. The only
differences are in how the HTML report *presents* an unchanged result, and each one is
attributable to a named commit and is a correction of a documented reporting defect.

| | |
| --- | --- |
| Baseline ("before") | `2fcc6e3` — merge-base with `main` |
| Candidate ("after") | `7344c62` — tip of `test/issue-179-test-strategy` |
| Cases per side | 58 (50 BAM x assembly, 5 non-fast-mode, 3 adVNTR) |
| Runs total | 116, plus 6 deliberate-mismatch probes |
| Non-zero exits | 0 before, 0 after |

## Method

The cohort is every BAM under `tests/data/`, run at each assembly it is provided for:
the 7 multi-reference samples at all six assemblies (`hg19`, `hg38`, `GRCh37`, `GRCh38`,
`hg19_ensembl`, `hg38_ensembl`) plus their original hg19 subsets, and the hg38 regression
guard `example_40cf`. Five cases repeat without `--fast-mode` so the unmapped-read pipes
are exercised, and three run `--extra-modules advntr`.

Compared per case: the full `kestrel_result.tsv` row set (all 28 columns, keyed on
`Motifs`/`POS`/`REF`/`ALT`/`Variant`), the pre-filter `kestrel_pre_result.tsv`,
`output_adVNTR_result.tsv` where adVNTR ran, `coverage_summary.tsv`, the report's
screening-summary sentence and its computed emphasis, the recorded pipeline steps, and
the exit code.

### Verifying which code actually ran

This is the part that could have made the whole exercise worthless. The `vntyper` console
script resolves the package through setuptools' editable finder, which is *appended* to
`sys.meta_path` and points at whichever worktree the editable install was made from —
irrespective of the current directory. Demonstrated: with the process CWD set to the
`2fcc6e3` worktree, a plain import still resolved to the candidate worktree and reported
`vntyper.scripts.pipeline_guards` as present, i.e. it would have run the new code.

Both sides were therefore launched through a wrapper that sets `PYTHONPATH` to its own
tree and refuses to start unless (a) `vntyper.__file__` lies under that tree and (b) the
marker module `vntyper.scripts.pipeline_guards`, which exists only after `078a6c4`, is
absent on the before side and present on the after side. Each run logs the resolved path
as its first line; all 58 before-side logs report the baseline tree with the marker
**absent** and all 58 after-side logs report the candidate tree with the marker
**present**. The before-side worktree carried its own `vntyper/config.json`,
`advntr_config.json` and `report_config.json`, so the baseline configuration was in force
as well as the baseline code.

## Result

| Compared | Cases with a delta |
| --- | --- |
| Exit code | 0 / 58 |
| Kestrel variant set (rows added or removed) | 0 / 58 |
| Kestrel `Confidence` | 0 / 58 |
| Kestrel `Flag` | 0 / 58 |
| Any other Kestrel column | 0 / 58 |
| adVNTR variant set, genotype fields, `Flag` | 0 / 3 |
| Coverage summary | 0 / 58 |
| Recorded pipeline steps | 0 / 58 |
| Screening summary **text** | 1 / 58 |
| Screening summary **emphasis** | 11 / 58 |

`kestrel_result.tsv`, `kestrel_pre_result.tsv`, `output_adVNTR_result.tsv` and
`coverage_summary.tsv` hash identically on both sides once the `##` provenance banner
(which carries the analysis timestamp) is removed.

The call each sample yields, identical on all six assemblies and identical on both sides:

| Sample | Motifs | Variant | Confidence | Flag |
| --- | --- | --- | --- | --- |
| `example_6449` | `4-5` | Insertion | `High_Precision*` | Not flagged |
| `example_66bf` | `5C-Q` | Insertion | `High_Precision*` | Not flagged |
| `example_6c28` | `S-Q` | Insertion | `High_Precision*` | Not flagged |
| `example_7a61` | — | — | `Negative` | — |
| `example_a5c1` | `L-6p` | Insertion | `High_Precision` | Not flagged |
| `example_b178` | `D-C` | Insertion | `High_Precision*` | Not flagged |
| `example_dfc3` | `5-E` | Deletion | `High_Precision*` | Not flagged |
| `example_40cf` (hg38 only) | — | — | `Negative` | — |

## Adjudication of every difference

### D1 — screening emphasis on 10 negative cases (`5527a49`, expected)

`example_7a61` at all six assemblies plus its non-fast rerun, and `example_40cf` at hg38
plus its non-fast rerun: the box lost its `summary-positive` styling. The sentence is
unchanged in both — "No variant detected. Note: adVNTR genotyping was not performed."

Emphasis now comes from the computed `summary_is_positive` state rather than from
searching the rendered sentence, so a message whose wording does not match what a text
search expects no longer carries the emphasis of the opposite result. The genotype was
`Negative` on both sides; only the styling of the correct result changed, and it changed
in the safe direction.

### D2 — screening sentence on `dfc3_hg19_advntr` (`77d590b`, expected)

| | |
| --- | --- |
| Before | "The screening was negative (no valid Kestrel or adVNTR data)." |
| After | "Pathogenic frameshift variant detected by Kestrel with high precision, and adVNTR detected the variant with a flagged result. Quality metrics are acceptable. Review the flagged adVNTR result and validate using orthogonal methods…" |

This is the most consequential finding of the exercise, and it is a **fix, not a
regression**. Kestrel called `High_Precision*` and adVNTR called VID 25561
`D17_2&D18_2&D19_2&D20_2&D21_2` flagged `Polymorphic_Call`, on both sides identically.
The state (`kestrel_result = High_Precision`, `advntr_result = positive flagged`,
`quality_metrics_pass = true`) had no rule in `report_config.json` and therefore fell
through to `screening_summary_default`, which is the negative sentence. `77d590b` added
the 15 rules that were missing, so every reachable state now has one and no state falls
through to the default.

Attribution note: the Kestrel and adVNTR result files are byte-identical here, so this is
**not** caused by `52f822e`. The flag was already firing before the change.

### D3 — assembly guard (`078a6c4`, no effect on the cohort)

The guard reached a decided verdict on all 58 after-side cases and passed every one
(20 `hg19`, 9 `hg38`, 8 `GRCh38`, 7 each `GRCh37` / `hg19_ensembl` / `hg38_ensembl`). No
warnings, no undetermined verdicts, **no rejections** — confirming the earlier finding
that it rejects no cohort sample.

A deliberate mismatch (hg19 BAM declared `hg38`, and the reverse) was probed on both
sides. Both sides exit 1; only the failure point moves — from a downstream
`RuntimeError: No coverage data found` to the guard's own message naming the detected
build and the flag to retry with. No run that succeeded before fails now.

### D4 — adVNTR `Repeat_Unit_7` rule (`52f822e`, not exercised)

The revived rule did not fire anywhere: the three adVNTR runs produced `RU` values `2`,
`4` and `2,2,2,2,2`. No cohort sample carries an RU-7 call, so the rule's effect remains
unobserved by this gate. Its blast radius is bounded by construction — `add_flags`
appends a column and filters nothing — and the adVNTR `Flag` column is unchanged on all
three runs.

### D5 — compound adVNTR parsing (`a7c3d9e`, byte-identical)

`example_dfc3` produces the compound call `D17_2&D18_2&D19_2&D20_2&D21_2`, which parsed
to an identical 10-column row on both sides. The crashing input class (a compound call
containing `LEN`) does not occur in this cohort, so the fix is confirmed
non-regressive here but not exercised on the defect itself.

**Follow-up.** This gate could not see it, but `a7c3d9e` was *not* byte-identical off the
cohort: replacing the greedy `(LEN.*)` with a bounded `(LEN\d+)` changed `Insertion_len`
for compound states whose `LEN` is followed by a further `&` part, and therefore changed
which rows survive the frameshift filter — a reported-genotype change for inputs that
never crashed. Restored to crash-only in a later commit on this branch by keeping the
greedy pattern and bounding the *split* instead; a differential sweep of 2380 probes
against `a7c3d9e^` shows 572/572 previously non-crashing inputs identical. The underlying
"what should `Insertion_len` be for a compound call" question is filed for the domain
owner as **B8** in [CI/CD follow-ups](ci-followups.md).

### D6 — pipefail and CRAM samtools (`331ea95`, no effect)

Five cases ran without `--fast-mode`, taking the unmapped and partially-mapped read path
that carries the three newly `pipefail`-guarded pipes. All five produced identical output
and exit 0 on both sides. No pipe stage was failing silently in this cohort. The CRAM
process-substitution change is not exercised: the cohort has no CRAM input.

### D7 — region-string fallback (`5486c84`, not exercised)

Probed with an NCBI-named BAM declared as `hg38`. The dynamic path resolved the region
correctly on both sides, so the legacy fallback was never taken and the new refusal never
triggered. Output identical, exit 0 both sides.

## What this gate does not cover

* CRAM input, and therefore the CRAM branch of `331ea95`.
* FASTQ input (the shark case), which the assembly guard deliberately does not guard.
* An adVNTR call with `RU == 7`, an adVNTR compound call containing `LEN`, and a BAM whose
  header lacks the declared assembly's legacy contig. All three are covered by unit tests
  only; no sample in `tests/data/` reaches them.
* Cohort mode (`vntyper cohort`) and the report subcommand.
