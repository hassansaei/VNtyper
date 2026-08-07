# Golden-cohort gate (#179)

Before-versus-after comparison of genotyping output across the whole local test cohort,
run to decide whether the deliberately behaviour-changing commits in #179 may ship.

Every run the gate has taken is registered in the table below. **A verdict is intended to
attest one candidate commit and nothing after it**, so read the run whose candidate
matches the tree you are judging.

**How the candidate commit is known, and by whom.** For runs 1–4 it is the operator's
record, not the instrument's: the harness that produced them (version `1.0.0`) never ran
`git rev-parse`, never looked at whether the working tree was clean, and never accepted an
expected SHA. Each side's `side.json` recorded a *path*, and a path is a different commit
ten minutes later. So those four candidate SHAs are asserted, and nothing in the run
artefacts can confirm or contradict them. From harness `1.1.0` that changes, and run 5 is
the first to have it: every side records its `HEAD`, its branch and whether `vntyper/`,
`docker/` or `scripts/` had uncommitted changes when it ran, `compare` refuses a recorded
revision that disagrees with `--expect-before-sha` / `--expect-after-sha`, and
`--require-clean` refuses a side that ran over uncommitted edits. **Run 5's candidate is a
recorded fact; runs 1–4's are not.** Run 5's *baseline* still is not — `4fd638a` predates
the harness change, so only the "after" side can prove its revision.

| Run | Candidate ("after") | Baseline | Verdict | Attests |
| --- | --- | --- | --- | --- |
| 1 | `7344c62` | `2fcc6e3` | PASS | the branch as it stood at `7344c62`, and nothing after it |
| 2 | `1792345` | `2fcc6e3` | PASS | the branch at `1792345`, including `d144505`, `c51052c`, `b4059ce`, `7e58eb8`, `2c92096`, and nothing after it |
| 3 | `8537a61` | `2fcc6e3` | PASS | the branch at `8537a61`, including `2ae28c5`, `2aa095a`, `50d7968`, `42c976a`, `4ce5639`, `97033d3`, `22e3d17`, `52a0ec9`, and nothing after it |
| 4 | `ec67fff` | `4fd638a` | PASS with two attributed deltas, neither genotype-affecting | the `fix/issue-181-197-followups` branch at `ec67fff`, against the 2.0.6 release, and nothing after it |
| 5 | `9816f86` | `4fd638a` | DELTAS, both classes fully attributed, every genotype artefact unchanged | the `fix/issue-181-197-followups` branch at `9816f86`, against the 2.0.6 release, and nothing after it |
| 6 | `b27ff9c` | `cb593b6` | DELTAS, every one attributed, **no genotype field changed anywhere** | the `fix/milestone-2-correctness-of-reported-numbers` branch at `b27ff9c` (milestone 2, #171/#172/#174/#203/#212), against the 2.0.7 release, and nothing after it |

Runs 1–3 measure the `#179` branch against the baseline `2fcc6e3`. Runs 4 and 5 measure a
*different* branch — `fix/issue-181-197-followups` — against a *different* baseline,
`4fd638a`, the merge-base with `main` and the 2.0.6 release. **Run 5 supersedes run 4 as
the attestation of record for that branch**: run 4's candidate `ec67fff` is no longer
code-identical to the tip, because the Phase-5 fixes and the commits answering the
adversarial review landed after it. Runs 1–3 remain the measurements the adjudications
below were written against, and each still attests its own candidate exactly.

No run attests "the branch tip" as a standing property — a tip moves, a commit does not.
Each run's candidate is named above and in its own result section. Whether a given tip is
covered by a given run is a question to answer with a command rather than with prose:

```
git diff --stat <candidate>..HEAD -- vntyper/ docker/
```

If that is empty, the tip runs the code the run measured; if it is not, read what it
names before trusting the run. Against run 5's candidate at the 2.0.7 release commit it
names `vntyper/version.py` and nothing else — a version string the pipeline reports but
does not branch on.

(An earlier revision of this section asserted the *result* of exactly that command for
run 4 — "empty", plus a count of the commits after `ec67fff` — and both went stale within
the day. The command does not go stale; its transcribed output does. `scripts/` is left
out deliberately: it holds the gate harness and the mutation tooling, neither of which the
pipeline imports.)

**The genotype verdict is PASS in every run.** Every genotype field, every `Confidence`
label and every `Flag` is byte-identical between baseline and candidate, on every sample
and every assembly, in all **six** runs. Run 6 states this as a measurement rather than a
reading: across its whole matrix, exactly two columns were added and the only column whose
cells changed was one of them, with **no genotype field touched anywhere**.

Runs 4 and 5 are not *only* that, and the difference is named rather than smoothed over.
In runs 1–3 the only differences were in how the HTML report *presents* an unchanged
result. Run 4 carries **two delta classes**, one in a version-probe command string and one
in the cohort export files; run 5 carries two of its own, again a version-probe command
string and a set of leaked working columns. Each is attributable to a named commit and
none touches a genotype — but none is a report-presentation delta either, so those two
runs are recorded with their deltas rather than folded into the sentence above.

| | Run 1 | Run 2 | Run 3 | Run 4 |
| --- | --- | --- | --- | --- |
| Baseline ("before") | `2fcc6e3` — merge-base with `main` | `2fcc6e3` | `2fcc6e3` | `4fd638a` — merge-base with `main`, the 2.0.6 release |
| Candidate ("after") | `7344c62` | `1792345` | `8537a61` | `ec67fff` |
| Cases per side | 58 (50 BAM x assembly, 5 non-fast-mode, 3 adVNTR) | same 58, same matrix | same 58, same matrix | the same 58, **derived from `tests/data` rather than hardcoded**, plus 4 cohort-mode cases |
| Runs total | 116, plus 6 deliberate-mismatch probes | 116, plus 6 probes | 116, plus 6 probes | 130 (65 per side), probes and cohort cases included |
| Non-zero exits | 0 before, 0 after | 0 before, 0 after | 0 before, 0 after | 0 before, 0 after |
| Executed shell commands compared | no | no | **yes** — 480 per side | **yes** — 1,111 per side across 61 cases |
| Cohort mode covered | no | no | no | **yes** — 4 cases |

This side-by-side stops at run 4 on purpose. Run 5 shares run 4's baseline, matrix and
case count, and differs in the two things that are worth stating rather than tabulating —
its candidate is recorded by the instrument instead of asserted by the operator, and its
deltas are a different pair. The "Run 5" section below has them.

## Method

The cohort is every BAM under `tests/data/`, run at each assembly it is provided for:
the 7 multi-reference samples at all six assemblies (`hg19`, `hg38`, `GRCh37`, `GRCh38`,
`hg19_ensembl`, `hg38_ensembl`) plus their original hg19 subsets, and the hg38 regression
guard `example_40cf`. Five cases repeat without `--fast-mode` so the unmapped-read pipes
are exercised, three run `--extra-modules advntr`, and two repeat from a derived CRAM. See
[The CRAM group](#the-cram-group-188) below. **Run 6 is the first to take the CRAM group**;
runs 1–5 predate the fixtures.

**The matrix is 60 cases and was 58 for runs 1–5.** Every `x / 58` figure in the run
sections below is that earlier matrix and is left as measured; run 6's tables are over 60
and are not comparable cell-for-cell with them.

The CRAM fixtures are **derived, not committed** (`scripts/make_cram_fixtures.py`), so a
fresh clone has 58 cases until they are generated. The harness refuses to launch over the
reduced matrix rather than running it silently — run 6 hit exactly that and generated the
50 fixtures instead of passing `--allow-matrix-drift`.

Compared per case: the complete `kestrel_result.tsv` header and row set, keyed on
`Motifs`/`POS`/`REF`/`ALT`/`Variant` — every column that is present, without asserting a
column count. (An earlier version of this sentence said "all 28 columns". There is no
28-column `kestrel_result.tsv` in run 4: 49 of the 59 files carry **27** columns and the
ten negative-call sentinels carry **10**. The comparator never had a count to be right or
wrong about — it diffs `columns_added` / `columns_removed` and the keyed rows — so the
number was decoration, and wrong decoration.) Also the pre-filter `kestrel_pre_result.tsv`,
`output_adVNTR_result.tsv` where adVNTR ran, `coverage_summary.tsv`, the report's
screening-summary sentence and its computed emphasis, the recorded pipeline steps, and
the exit code. Run 3 adds the executed shell command strings, recorded at the
`subprocess` boundary.

Run 4 derives the **50 base cases** from `tests/data` at run time rather than reproducing a
hardcoded list. The rest of the matrix is *not* derived and this page has said otherwise:
the five non-fast ids, the three adVNTR ids and the three probes are declared policy
(`NON_FAST_CASE_IDS`, `ADVNTR_CASE_IDS`, `PROBE_SPECS` in
`scripts/golden_cohort/matrix.py`), resolved against the derived set so that a policy
naming a case the data no longer provides is an error rather than a silent shrink. Only the
adVNTR selection is recoverable from this page; the non-fast one is a reconstruction, and
`matrix.py`'s docstring says so. Run 4 also adds four `vntyper cohort` cases — `cohort_multi`,
`cohort_multi_pseudonymized`, `cohort_single` and `cohort_empty` — comparing each cohort
export (`cohort_kestrel_{csv,tsv,json}`, `cohort_advntr_{csv,tsv,json}`), the rendered
cohort tables, the category counts and totals, the set of cohort output files, and the
pseudonymization table.

### The CRAM group (#188)

VNtyper accepts CRAM, and up to and including run 5 no gate run had ever given it one —
so the CRAM branch of `process_bam_to_fastq` and the process-substitution write race
`175011e` fixed in `build_cram_unmapped_filter_command` were attested by unit tests and one
hand-run equivalence comparison, never by this gate.

`make cram-fixtures` (`scripts/make_cram_fixtures.py`) closes that. It derives a CRAM beside
every cohort BAM under `tests/data/cram/`, mirroring the source layout with `.bam` →
`.cram`, and proves each one lossless: the decoded record stream digests identically to its
source. The fixtures are derived rather than committed because `tests/data/` is
git-ignored and ships as a Zenodo archive. They are written `no_ref=1` — the cohort's BAM
headers carry no `M5` tags, so no reference can be resolved by digest, and
`process_bam_to_fastq` passes an empty `cram_ref_option` unconditionally. **A `no_ref` CRAM
exercises the container format, the CRAM decoder, `.crai` indexing and the unmapped-read
scan. It does not exercise reference resolution, because it needs none.** The ordinary
externally-referenced CRAM a diagnostic lab would send is a different fixture and remains
uncovered.

Two cases are declared, in `CRAM_CASE_IDS` in `scripts/golden_cohort/matrix.py`. Like the
non-fast and adVNTR selections this is policy, not derivation — only the *fixture paths* are
derived, mirrored from each base case's BAM path so they cannot disagree with what
`make_cram_fixtures.py` wrote.

| Case | Repeat of | Records | Unmapped pairs | Why this one |
| --- | --- | --- | --- | --- |
| `b178_hg19_cram` | `b178_hg19_subset` | 34,214 | 4,478 | A known positive (`D-C` insertion, `High_Precision*`, not flagged — run 1's table below). It is what shows a CRAM run still **calls**, not merely that it exits 0. |
| `7a61_hg38_ensembl_cram` | `7a61_hg38_ensembl_bwa` | 985,731 | 622,690 | A heavy unmapped load, and so exposed to the write race `175011e` fixed — measured there at 199,797 of 200,000 unmapped reads present at the instant the shell returned. |

Counts are from `tests/data/cram/manifest.json`. `7a61_hg38_ensembl_bwa` is **not** the
single heaviest case in the cohort — `7a61_hg19_subset` carries 958,804 unmapped pairs and
the six remapped `7a61` cases tie at 623,792 / 622,690. This pair is the one already proven
end to end, and it also covers both layouts the fixture tree mirrors: one top-level subset
BAM and one `remapped/<aligner>/<assembly>/` BAM.

Both run **without `--fast-mode`, and that is the whole point.** `--fast-mode` skips the
unmapped-read extraction entirely, and the CRAM-specific extraction lives inside that
branch — so a fast-mode CRAM case would exercise the slice and the FASTQ conversion and
none of the code the fixtures exist for.

**A declared CRAM case whose fixture has not been derived is skipped and logged at error
level, and the group then comes out short.** That is an ordinary drift mismatch: a strict
build refuses it, and `--allow-matrix-drift` runs it knowingly as a non-attestation run.
There is deliberately no "0 or 2 CRAM cases are both fine" rule — a run without them covers
strictly less than this contract records, which is exactly what the `REDUCED` verdict is
for.

### Verifying which code actually ran

This is the part that could have made the whole exercise worthless. The `vntyper` console
script resolves the package through setuptools' editable finder, which is *appended* to
`sys.meta_path` and points at whichever worktree the editable install was made from —
irrespective of the current directory. Demonstrated: with the process CWD set to the
`2fcc6e3` worktree, a plain import still resolved to the candidate worktree and reported
`vntyper.scripts.pipeline_guards` as present, i.e. it would have run the new code.

Both sides are therefore launched through a wrapper that sets `PYTHONPATH` to its own
tree and refuses to start unless (a) `vntyper.__file__` lies under that tree and (b) the
marker module `vntyper.scripts.pipeline_guards`, which exists only after `078a6c4`, is
absent on the before side and present on the after side. Each run logs the resolved path
as its first line. Every run of every gate has passed that check on its own side, with no
run reaching the pipeline through the wrong package: run 1, all 58 logs per side; runs 2
and 3, all 61 logs per side; run 4, all 65 runs per side, zero aborts. Run 4's marker is
`vntyper.scripts.cohort_rules` rather than `pipeline_guards` — that module does not exist
at `4fd638a` and does at `ec67fff`, which is what the baseline moving to the 2.0.6 release
requires. The before-side worktree carries its own tracked
configuration (`vntyper/config.json` and the report/adVNTR configuration files of that
commit), so the baseline configuration is in force as well as the baseline code.

## Result — run 1, candidate `7344c62`

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

## Result — run 2, candidate `1792345`

Same baseline, same 58-case matrix, same comparison. Re-run because run 1's candidate is
24 commits behind `1792345` and five of those commits change production code: `d144505`
(adVNTR compound-variant repair restored to crash-only), `c51052c` and `b4059ce`
(assembly guard), `7e58eb8` (shell quoting at three call sites) and `2c92096`
(`--output-name` rejection).

| Compared | Cases with a delta |
| --- | --- |
| Exit code | 0 / 58 |
| Kestrel variant set (rows added or removed) | 0 / 58 |
| Kestrel `Confidence` | 0 / 58 |
| Kestrel `Flag` | 0 / 58 |
| Any other Kestrel column, including `Motifs` | 0 / 58 |
| `kestrel_pre_result.tsv` | 0 / 58 |
| adVNTR variant set, genotype fields, `Flag` (**not** `Insertion_len` — corrected below) | 0 / 3 |
| Coverage summary | 0 / 58 |
| Recorded pipeline steps | 0 / 58 |
| Screening summary **text** | 1 / 58 |
| Screening summary **emphasis** | 11 / 58 |
| Rendered `Motif` **cell** (display only) | 48 / 58 |
| HTML entity escaping (display only) | 58 / 58 |

The per-sample call table is unchanged from run 1 — same six motif pairs, same
`Confidence`, same `Flag`, on all six assemblies. The three adVNTR runs reproduce run 1
exactly, VID `25561` throughout:

| Case | RU | Pvalue | Flag |
| --- | --- | --- | --- |
| `a5c1_hg19_advntr` | `2` | `6.78296229901e-07` | Not flagged |
| `b178_hg19_advntr` | `4` | `3.82652062679e-56` | Not flagged |
| `dfc3_hg19_advntr` | `2,2,2,2,2` | `1.5504014332800002e-18` | `Polymorphic_Call` |

`dfc3` is the compound call `D17_2&D18_2&D19_2&D20_2&D21_2`, the input class `d144505`
rewrote. Its row is identical on both sides in every column of the adVNTR output schema.
(This sentence originally read "in every column, `Insertion_len` included". That was
wrong — `Insertion_len` is not a column of that file. See the correction below.)

### Run 2's proof of which code ran

Same failure mode, same defence, re-demonstrated on this tree. With the process CWD set
to the `2fcc6e3` worktree and no `PYTHONPATH`, a script run from outside the tree resolved
`vntyper.__file__` to `…/issue-179-impl/vntyper/__init__.py` and reported
`vntyper.scripts.pipeline_guards` **present** — i.e. it would have run the candidate code
while appearing to run the baseline. With `PYTHONPATH` pinned to the baseline worktree the
same script resolved to `…/before-2fcc6e3/vntyper/__init__.py` with the marker **absent**.

Every one of the 122 runs was launched through a wrapper that prints its resolved
`vntyper.__file__` and marker state as its first line and exits before doing any work
unless both agree with its side. All 61 before-side runs report the baseline tree with
`pipeline_guards=False`; all 61 after-side runs report the candidate tree with
`pipeline_guards=True`. No run reached the pipeline with the wrong package.

### Run 2's assembly-guard verdicts — `b4059ce` rejects nothing

The guard reached `agree` on **58 of 58** after-side cases: 20 `hg19`, 9 `hg38`,
8 `GRCh38`, 7 each `GRCh37` / `hg19_ensembl` / `hg38_ensembl`. Zero `mismatch`, zero
`undetermined`, and — the question `b4059ce` raises — **zero `conflict`**. No cohort
header names two builds at once, so making that state fatal rejects no sample that
previously succeeded. `c51052c` is confirmed the same way: the verdict no longer depends
on contig order and every cohort header still decides.

The two deliberate mismatch probes behave as in run 1: exit 1 on both sides, only the
failure point moves. One consequence is worth naming because it is the sole non-genotype
difference outside the report — on the after side those two runs write **no**
`pipeline_summary.json`, because the guard refuses before the first step is recorded,
whereas the baseline got as far as recording `BAM Header Parsing` and
`BAM to FASTQ Conversion` before failing downstream. Both sides still exit 1; no run that
succeeded before fails now. The naming probe (a `GRCh38` BAM declared `hg38`) exits 0 on
both sides and the guard agrees.

### Run 2's presentation deltas

The four in the run-1 adjudication below reproduce identically (D1, ten negative cases;
D2, `dfc3_hg19_advntr`). Run 2 additionally quantified two the run-1 record mentions but
did not count, both confined to the rendered HTML:

* **Rendered `Motif` cell**, 48 of 58 — `8a76512`. The column now shows the annotated
  motif rather than the raw pair: `4-5`→`5`, `5C-Q`→`5C`, `S-Q`→`Q`, `L-6p`→`L`,
  `D-C`→`D`, `5-E`→`E`. The `Motifs` column of `kestrel_result.tsv` is unchanged on all
  58 cases, asserted directly — this is display only.
* **HTML entity escaping**, 58 of 58 — `bda7e05`. The `&&` inside the recorded samtools
  command is now emitted as `&amp;&amp;`. Two occurrences per report; the command string
  itself is unchanged.

No presentation delta in run 2 is unattributed.

### How run 2 differs procedurally from run 1

Three deviations, none of which touch what is compared: both sides read the BAMs from the
`issue-179-impl` worktree's `tests/data` by absolute path (run 1 used a per-side relative
path), every run used `--threads 2`, and eight cases ran concurrently. The case matrix
itself was reproduced from run 1's `make_matrix.py` verbatim — same 58 ids, same
assemblies, same five non-fast samples, same three adVNTR samples with
`--advntr-max-coverage 300`, same three probes.

## Result — run 3, candidate `8537a61`

Same baseline `2fcc6e3`, same 58-case matrix, same three probes, same comparison. Re-run
because run 2's candidate `1792345` is eight commits behind `8537a61` and one of those
commits changes how production commands are *built*: `2ae28c5` applies `shlex.quote` at
five previously-unquoted shell interpolation sites (`utils.py::validate_bam_file`, the
SAM→BAM conversion, its index and the bcftools sort in `kestrel_genotyping.py`, and the
aligner `index_command` templates in `install_references.py`). Quoting a path that needs
no quoting should be a byte-level no-op in the executed command; run 3 checks that
empirically rather than by argument.

| Compared | Cases with a delta |
| --- | --- |
| Exit code | 0 / 58 |
| Kestrel variant set (rows added or removed) | 0 / 58 |
| Kestrel `Confidence` | 0 / 58 |
| Kestrel `Flag` | 0 / 58 |
| Any other Kestrel column, including `Motifs` | 0 / 58 |
| `kestrel_pre_result.tsv` | 0 / 58 |
| adVNTR variant set, genotype fields, `Flag` (**not** `Insertion_len` — corrected below) | 0 / 3 |
| Coverage summary | 0 / 58 |
| Recorded pipeline steps | 0 / 58 |
| **Quoted-site command strings** (quickcheck, SAM→BAM, index, bcftools sort) | **0 / 58** |
| Screening summary **text** | 1 / 58 |
| Screening summary **emphasis** | 11 / 58 |
| Cross-match **emphasis** (display only) | 1 / 3 adVNTR cases |
| Rendered `Motif` **cell** / column (display only) | 58 / 58 |
| HTML entity escaping in the embedded log (display only) | 58 / 58 |
| IGV script fragment (display only) | 58 / 58 |

The per-sample call table is unchanged from runs 1 and 2 — same six motif pairs, same
`Confidence`, same `Flag`, on all six assemblies, and `example_7a61` / `example_40cf`
`Negative` throughout. The three adVNTR runs reproduce run 2 exactly: VID `25561`, RU `2`
/ `4` / `2,2,2,2,2`, P-values `6.78296229901e-07` / `3.82652062679e-56` /
`1.5504014332800002e-18`, flags `Not flagged` / `Not flagged` / `Polymorphic_Call`.

The assembly guard reached `agree` on **58 of 58** after-side cases again — zero
`mismatch`, zero `conflict`, zero `undetermined` — and the probes behave exactly as in
run 2: both mismatch probes exit 1 on both sides with only the failure point moving (from
`Error calculating coverage summary: No coverage data found` to the guard naming the
detected build and the flag to retry with), and the naming probe exits 0 on both sides.

### Run 3's proof of which code ran

Same failure mode, same defence, re-demonstrated on this tree before the runs started.
With the process CWD set to the `2fcc6e3` worktree and no `PYTHONPATH`, a plain
`import vntyper` resolved to
`…/issue-179-impl/vntyper/__init__.py` with `vntyper.scripts.pipeline_guards`
**present** — the baseline checkout would have executed candidate code. With `PYTHONPATH`
pinned to the baseline worktree the same probe resolved to
`…/gate3/before/vntyper/__init__.py` with the marker **absent**.

Every one of the 122 runs was launched through a wrapper that prints its resolved
`vntyper.__file__` and marker state as its first line and `sys.exit`s before dispatch
unless both agree with its side. **All 61 before-side runs report the baseline tree with
`pipeline_guards=absent`; all 61 after-side runs report the candidate tree with
`pipeline_guards=present`; 0 aborts.** The baseline side ran from its own `git worktree`
at `2fcc6e3` (removed afterwards) carrying its own `vntyper/config.json`, with only the
untracked data and reference directories symlinked in, so the baseline configuration was
in force as well as the baseline code.

### Run 3's command-string comparison — `2ae28c5` is a no-op here

New in run 3: the launcher wraps `subprocess.Popen`/`subprocess.run` and records every
command each run actually executed, so the two sides' command streams are diffed as well
as their outputs. 480 shell commands per side across the 58 cases, an identical count.
After normalising the two things that differ by construction (the per-side output
directory and the source root), the command streams differ in exactly two ways:

* **`set -o pipefail; ` prefixed to the `samtools sort | samtools fastq` pipe**, 58 / 58 —
  `331ea95`, already adjudicated as D6.
* **Extra and reordered `samtools view -H` header reads on the after side** — the assembly
  guard reads the header before slicing. Two cases in most runs against the baseline's
  one-early-one-late ordering; a strictly higher count on five cases. `078a6c4` /
  `c51052c`.

**No command differs in quoting.** The four sites `2ae28c5` touches that this cohort
reaches are byte-identical on every case: `samtools quickcheck` 58/58, `samtools view -Sb`
58/58, `samtools index` 58/58, `bcftools sort` 58/58 (and 3/3, 1/1, 1/1, 1/1 on the
probes). `shlex.quote` leaves anything matching `[\w@%+=:,./-]*` alone and no cohort path
carries a shell metacharacter, so the quoting is invisible in the executed command — which
is what the commit claims and what this run now measures rather than asserts. The fifth
site, the aligner `index_command` templates in `install_references.py`, is **not
exercised**: the references are already installed, so no run built an index.

### Run 3's presentation deltas — every one attributed

| Delta | Cases | Commit |
| --- | --- | --- |
| Screening emphasis lost on negative cases (D1) | 10 | `5527a49` |
| Screening sentence on `dfc3_hg19_advntr` (D2) | 1 | `77d590b` |
| Rendered `Motif` cell shows the annotated motif (`4-5`→`5`, `5C-Q`→`5C`, `S-Q`→`Q`, `L-6p`→`L`, `D-C`→`D`, `5-E`→`E`) | 48 | `8a76512` |
| `Motif` column header plus a `None` cell added to the negative-case table | 10 | `8a76512` |
| HTML entity escaping of the embedded pipeline log (`'`→`&#39;`, `>`→`&gt;`, `&&`→`&amp;&amp;`) | 58 | `bda7e05` |
| Empty IGV fragments plus the `initIGV()` guard (`const tableJson = ;` → `{"headers": [], "rows": []}`) | 58 | `2180de6` |
| New `pipeline_guards` / `chromosome_utils` lines inside the embedded log | 58 | `078a6c4` |
| `set -o pipefail; ` inside the logged BAM→FASTQ command | 58 | `331ea95` |
| **Cross-match emphasis** on `dfc3_hg19_advntr`: `summary-box summary-positive` → `summary-box`, sentence unchanged ("No matches were found between Kestrel and adVNTR results.") | 1 of 3 | **`2aa095a`** |

The last row is the new one this run adds, and it is the delta `2aa095a` exists to
produce: the only cohort case whose cross-match state is negative is `dfc3_hg19_advntr`,
and its "No matches were found" sentence was being rendered in the positive style. The
other two adVNTR cases genuinely match and keep `summary-positive` on both sides. No
genotype field moves with it.

`50d7968` produced no delta here because it changes the **cohort** report
(`cohort_summary.py`) and this gate runs no cohort-mode case; the per-sample report is
untouched by it. No presentation delta in run 3 is unattributed.

### How run 3 differs procedurally from run 2

The matrix, probes and comparison are run 2's, reproduced. Three deviations, none of which
touches what is compared: `--threads 4` (8 for the adVNTR cases, as run 2's driver
specified) and six cases concurrent; the baseline worktree lived under the run's own
scratch directory rather than beside run 2's; and the launcher additionally records
executed commands, which is an observation, not a change to the run.

## Result — run 4, candidate `ec67fff`

New branch, new baseline. Run 4 compares `fix/issue-181-197-followups` at `ec67fff`
against `4fd638a`, the merge-base with `main` and the 2.0.6 release. It is **the first run
in this project's history to cover cohort mode at all** — runs 1–3 did not, and this
page's own "What this gate does not cover" section said so.

The 50 base cases are derived from `tests/data` at run time rather than hardcoded; the
non-fast, adVNTR and probe selections are declared policy resolved against them (see
[Method](#method)). 58 cases (50 base, 5 non-fast, 3 adVNTR), plus the 3 probes, plus
**4 cohort-mode cases** —
`cohort_multi`, `cohort_multi_pseudonymized`, `cohort_single`, `cohort_empty`. 65 runs per
side, 130 in total. Marker module `vntyper.scripts.cohort_rules`, absent at `4fd638a` and
present at `ec67fff`: **all 65 runs on each side verified their package resolution before
doing any work, and there were zero aborts.**

| Compared | Cases with a delta | Cases compared |
| --- | --- | --- |
| `kestrel_result` | **0** | 59 |
| `kestrel_pre_result` | **0** | 59 |
| `advntr_result` | **0** | 3 |
| `coverage_summary` | **0** | 59 |
| `cross_match_summary` | **0** | 3 |
| `exit_code` | **0** | 65 |
| `pipeline_step_records` | **0** | 61 |
| `cohort_category_counts` | **0** | 3 |
| `cohort_category_totals` | **0** | 3 |
| `cohort_tables` | **0** | 3 |
| `cohort_output_files` | **0** | 4 |
| `executed_commands` | 61 | 61 |
| `cohort_kestrel_csv` / `_tsv` / `_json` | 3 each | 3 each |
| `cohort_advntr_csv` / `_tsv` / `_json` | 3 each | 3 each |

Four further artefacts the harness compares are omitted from the table only because they
are also 0: `pipeline_steps` 0/61, `report_tables` 0/59, `screening_summary` 0/59 and
`pseudonymization_table` 0/1. So the screening-summary sentence and emphasis — the source
of every presentation delta in runs 1–3 — do not move at all in run 4.

**Zero deltas on every genotype artefact.** Two delta classes, both intended, both
attributable to a named commit.

The harness's own verdict line reads `DELTAS`, not `PASS`. That is mechanical: it reports
whether anything differed, and something did. The `PASS` on this page is the adjudication
of those differences, made below and open to disagreement — the two are not the same claim
and the raw result file should not be read as endorsing this one.

### Run 4's delta 1 — the duplicate kestrel help flag (`2873ad3`)

`executed_commands` differs on 61 of 61 cases. Exactly one command changed:

| | |
| --- | --- |
| Before | `java -jar vntyper/dependencies/kestrel/kestrel.jar -h -jar vntyper/dependencies/kestrel/kestrel.jar -h` |
| After | `java -jar vntyper/dependencies/kestrel/kestrel.jar -h` |

That is the `get_tool_versions` duplicate-help-flag fix in `2873ad3`. It is a version probe
run once per invocation: it reads no sample data and feeds no genotype. That it appears on
every case is a property of running once per invocation, not evidence of breadth of
effect — and every genotype artefact in the table above is 0.

**The command counts.** 1,111 recorded commands per side across the 61 cases, and the
count matched between the two sides on **every one of the 61 cases** — which is the
statement that matters, because a changed *count* is what a new or dropped subprocess looks
like. The per-case count is not uniform: 42 cases record 18, and the rest run 9, 17, 19,
20, 21, 22, 24 or 28. An earlier version of this section said "the command count is
identical at 18 per side", which took the mode for the whole distribution and multiplied
out to 1,098 rather than the 1,111 actually recorded.

### Run 4's delta 2 — leaked working columns in the cohort exports (`90f61fa`)

`cohort_kestrel_{csv,tsv,json}` and `cohort_advntr_{csv,tsv,json}` each differ on 3 of 3
cohort cases that write exports (`cohort_empty` writes none — see D12 in
[CI/CD follow-ups](ci-followups.md)). That is the removal of the `__row_result` and
`__unified` working columns, which were leaking out of the cohort summary's internals into
every CSV, TSV and JSON export.

The corroborating evidence that only the exports moved is in the same table: the rendered
`cohort_tables` are **unchanged**, 0 of 3, as are `cohort_category_counts`,
`cohort_category_totals` and `cohort_output_files`. A change that altered what the cohort
*reports* — rather than which internal columns it spills into the export files — would have
moved those as well.

### What run 4 does not attest

This page has over-claimed before — see the correction below — so run 4's limits are named
here rather than left to inference. A PASS is only worth what its scope is.

**This is ONE run at the branch tip, not one run per genotype-affecting commit.** The plan
required a run per genotype-affecting commit, so that a failure could be attributed to a
single commit rather than bisected for afterwards. That was traded for a single tip run
plus bisect-on-failure: identical worst case, far cheaper expected case. The consequence
must be stated, because the trade does not buy it back — **a single run cannot detect two
changes producing offsetting deltas that cancel on the same field of the same sample.**
The only mitigation is that the comparison is per-sample and per-field across 58 cases, so
such a cancellation would have to be exact, on the same field, wherever both changes act.
That makes it unlikely. It does not make it impossible, and this run does not exclude it.

**The gate cannot attest #192 at all, and that too is measured rather than assumed.** The
cohort's only compound adVNTR state is `example_dfc3`'s `D17_2&D18_2&D19_2&D20_2&D21_2`,
which carries no `LEN` token — so `Insertion_len` is 0 under both the old and the new
semantics, and `advntr_result` being 0/3 here is silence, not confirmation. The evidence
for #192 is a differential sweep instead: 52,511 probes, 13,563 of 13,563 previously
parsing inputs byte-identical, and all 38,943 differences oracle-predicted.

**#184's PASS is weak evidence.** Exact float equality against `0.00515` requires an
alternate depth that is an exact multiple of 103, so the cohort very likely contains no
row on the boundary at all. The load-bearing evidence for #184 is its boundary table and a
multi-candidate selection test, not this gate.

**#185 is exercised only in the negative.** No cohort case is missing a gate column, so
what the run shows is that the new raise does not fire on healthy input — which is what it
is for, and is not the same as showing that it fires when it should.

**#188 is not exercised.** The cohort had no CRAM input when run 4 was taken, and run 5's
did not either. Its evidence is a hand-run end-to-end CRAM comparison, which is not in CI.
The fixtures and the two CRAM cases arrived afterwards — see
[The CRAM group](#the-cram-group-188) — and no run on this page has taken them.

**The Kestrel allele-shape guard is not exercised, and this is counted rather than
argued.** `102c46f` added `_assert_kestrel_allele_contract` in
`vntyper/scripts/file_processing.py`, which raises on a VCF record whose REF *and* ALT are
both longer than one base. The pinned Kestrel 1.0.1 cannot emit such a record — it anchors
every indel on a single reference base — and the run confirms it: across the **236 Kestrel
VCFs per side** (`output.vcf`, `output_indel.vcf`, `output_insertion.vcf`,
`output_deletion.vcf` for each of the 59 cases that reach Kestrel), **460,849 data records
per side, zero** carry two multi-base alleles. So the new raise never fired, and the
`filter_indel_vcf` re-routing it protects was never taken either. `advntr_result` and
`kestrel_result` both being 0 here is silence about that guard, not confirmation of it; its
evidence is `tests/unit/test_file_processing.py`.

**#195's per-row malformed-motif containment is not exercised either, for the same kind of
reason.** `11e2300` replaced a column-wide `str.count("-").max() != 1` gate — which let one
malformed motif ID suppress an entire sample's call — with a per-row drop. Firing it needs
a `Motifs` value that is not exactly two half-motif names joined by one dash, and the
cohort contains none: across the after side's 118 Kestrel tables, **44,227 non-empty
`Motifs` values (44,178 in `kestrel_pre_result.tsv`, 49 in `kestrel_result.tsv`) all
contain exactly one dash**, and the before side is identical. The containment therefore had
nothing to contain, and what run 4 shows is that the rewrite changes no call on well-formed
input — not that it contains a malformed one. That is
`tests/unit/test_motif_decisions.py`'s job.

**Cohort sample ordering is normalised away and is therefore not attested by this gate.**
The harness sorts cohort rows before comparing them. Two reasons survive scrutiny, and only
the first applies to run 4's pair:

1. **The baseline predates the determinism fix.** At `4fd638a`, `cohort_summary.py`
   iterates the discovery set directly (`for sample_dir in processed_dirs:`). `Path.__hash__`
   is the hash of the path string and Python randomises string hashing per process, so that
   side's row order differs between two runs *of itself*. Comparing order across such a pair
   measures the interpreter's hash seed.
2. **ZIP inputs, on any version.** Each ZIP extracts to
   `tempfile.mkdtemp(prefix="cohort_zip_")`, whose random suffix is part of the path and so
   part of the sort key.

The candidate sorts fixed input directories deterministically (`90f61fa`:
`return sorted(processed_dirs), temp_dirs`), so the ordering fix is attested by
`tests/unit/test_cohort_inputs.py::test_the_discovered_directories_come_back_sorted`,
`::test_the_order_is_lexicographic_by_path_part_rather_than_by_raw_string` and
`::test_processes_with_different_hash_seeds_discover_the_same_order` — not by this run. A
normalisation is a claim that a difference does not matter, and here it also means the fix
to that difference is invisible.

The harness's normalisation note used to cite
`tests/unit/test_cohort_inputs.py::test_discovery_returns_an_unordered_set_today` — a test
`90f61fa` renamed, so the citation named nothing in the repository — and to say discovery
"returns a set, so order is not reproducible", which is false for fixed directories on the
candidate. Both are corrected in `scripts/golden_cohort/compare.py`
(`COHORT_ORDER_WHY`), and `tests/unit/test_golden_cohort_compare.py` now reads
`test_cohort_inputs.py` and fails if the note cites a test that is not defined there.

## The instrument itself — harness `1.0.0` versus `1.1.0`

Runs 1–4 were produced by harness `1.0.0`. A review of that harness found that it could
return `IDENTICAL` over two runs that had both failed producing nothing, and `1.1.0` is
the response — the version run 5 was taken with. What changed, and what each change would
have done to run 4:

| Change | Effect on run 4, measured against its raw artefacts |
| --- | --- |
| Every case's declared `expect_exit` is enforced. It was written seven times in `matrix.py` and read nowhere, so two sides that both exited 1 without writing a genotype artefact compared `absent_both` on every field and earned `IDENTICAL`. | **None.** 0 expectation violations across all 65 cases on both sides; the two mismatch probes exit 1 as declared and the other 63 exit 0. |
| A case expected to exit zero must also have written its declared artefacts (`pipeline_summary.json`, both Kestrel tables, the coverage summary, the report; adVNTR stays optional; `cohort_empty` declares none, since it writes only its log by design). | **None.** All 59 zero-expected pipeline cases wrote all five on both sides; the three exporting cohort cases wrote all seven. |
| `compare` refuses two sides that share a run root, a source tree, a commit or a marker expectation, that are mislabelled, or that recorded no case results. | **None.** Run 4's two sides are properly opposed. |
| Each side records its `HEAD`, branch and working-tree state; `compare` can verify them. | **Not retroactive.** Run 4's sides have no `revision` key, and `compare` warns rather than refuses so existing run roots stay readable. |
| An unfiltered matrix that deviates from the documented per-group contract is refused before launching, a zero-case matrix always, and a clean result over a reduced matrix reads `REDUCED` rather than `IDENTICAL`. (The contract was 50 base / 5 non-fast / 3 adVNTR plus 3 probes when run 4 ran; it is now 50 / 5 / 3 / **2 CRAM** plus 3 probes.) | **None.** Run 4's `matrix.json` records zero mismatches and no filter, against the contract as it stood then. |
| `md5sum` is kept for step result files with no direct comparator — `pipeline_info.json` (which carries the assembly guard's verdict), `output_R1.fastq.gz` and `cross_match_results.tsv` — and dropped only for the three the harness parses row by row. | **None.** Those three checksums are identical between the two sides on 59/59, 59/59 and 3/3 cases. (`kestrel_result.tsv`'s differs on 59 of 59, which is what the original justification was written for and why it stays dropped.) |
| A changed `##` provenance banner now makes a table `differ` instead of being computed and discarded. | **None.** 0 provenance changes across all 180 compared tables. |

So `1.1.0` measures strictly more than `1.0.0` and, on run 4's artefacts, finds exactly
what `1.0.0` found. That is a check on the change, not a defence of the old harness: the
point of the fixes is the runs where the two would *not* agree.

## Correction — runs 2 and 3 overstated what the adVNTR comparison covered

Runs 2 and 3 each recorded a table row reading "adVNTR variant set, genotype fields,
`Insertion_len`, `Flag` — 0 / 3", and run 2's prose said `dfc3`'s row was "identical in
every column, `Insertion_len` included".

**`Insertion_len` is not a column of the adVNTR output.** `final_columns` in
`vntyper/modules/advntr/advntr_genotyping.py` is `VID`, `Variant`,
`NumberOfSupportingReads`, `MeanCoverage`, `Pvalue`, `RU`, `POS`, `REF`, `ALT`, `Flag`.
`Insertion_len` is an intermediate used by the frameshift filter and is dropped before the
file is written, so no comparison of `output_adVNTR_result.tsv` could have observed it.

The row-set comparison those runs made remains valid and their 0/3 result stands — it is
the `Insertion_len` claim specifically that was unsupported. The rows above have been
narrowed and this correction recorded rather than the claim quietly rewritten, because the
whole value of this page is that a reader can tell what was measured from what was
asserted. It also matters downstream: the D5 follow-up below turns on `Insertion_len`
changing for compound states, and this gate never had sight of that field.

## Adjudication of every difference

The deltas below were adjudicated against run 1 and re-observed unchanged in run 2.

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
never crashed. Restored to crash-only by `d144505` — later on this branch, and therefore
inside run 2 but not run 1 — by keeping the greedy pattern and bounding the *split*
instead; a differential sweep of 2380 probes
against `a7c3d9e^` shows 572/572 previously non-crashing inputs identical. The underlying
"what should `Insertion_len` be for a compound call" question is filed for the domain
owner as **B8** in [CI/CD follow-ups](ci-followups.md).

### D6 — pipefail and CRAM samtools (`331ea95`, no effect)

Five cases ran without `--fast-mode`, taking the unmapped and partially-mapped read path
that carries the three newly `pipefail`-guarded pipes. All five produced identical output
and exit 0 on both sides. No pipe stage was failing silently in this cohort. The CRAM
process-substitution change is not exercised: no run up to and including run 5 had a CRAM
input. (Two CRAM cases are now in the matrix — see [The CRAM group](#the-cram-group-188) —
but this adjudication is over runs 1 and 2 and is left as it was measured.)

### D7 — region-string fallback (`5486c84`, not exercised)

Probed with an NCBI-named BAM declared as `hg38`. The dynamic path resolved the region
correctly on both sides, so the legacy fallback was never taken and the new refusal never
triggered. Output identical, exit 0 both sides.

## What this gate does not cover

Still true of every run, run 4 included:

* FASTQ input (the shark case), which the assembly guard deliberately does not guard.
* The `vntyper report` subcommand.
* An adVNTR call with `RU == 7`, an adVNTR compound call containing `LEN`, and a BAM whose
  header lacks the declared assembly's legacy contig. All three are covered by unit tests
  only; no sample in `tests/data/` reaches them. The `LEN` gap is why the gate cannot
  attest #192.
* Reference installation, and therefore the aligner `index_command` quoting site
  `2ae28c5` touches in `install_references.py`: the references are already installed, so
  no gate run builds an index.
* A path containing a shell metacharacter. Run 3 shows `2ae28c5` changes no command in
  this cohort precisely because no cohort path needs quoting; the quoting itself is
  pinned by `tests/unit/test_shell_quoting.py`, not by this gate.
* `Insertion_len`, on any run. It is not a column of `output_adVNTR_result.tsv` — see the
  correction above.
* A Kestrel VCF record carrying two multi-base alleles, and therefore `102c46f`'s
  allele-shape guard: 0 of 460,849 records per side. See run 4's limits above.
* A `Motifs` value that is not two half-motif names joined by one dash, and therefore
  `11e2300`'s per-row containment: 0 of 44,227 non-empty values per side.
* **Which commit each side ran**, on runs 1–4. Harness `1.0.0` recorded a path; the SHAs
  above are the operator's record. Harness `1.1.0` records `HEAD` and the working-tree
  state per side and can be told to verify them.

Still true of runs 1–5, and addressed in the matrix but **not yet by any run**:

* CRAM input, and therefore the CRAM branch of `331ea95`, `175011e` and the whole of #188.
  No run up to and including run 5 fed VNtyper a CRAM; `tests/data` held eight BAMs and no
  CRAM when those runs were taken, and #188's evidence was a hand-run end-to-end CRAM
  comparison that is not in CI. `make cram-fixtures` now derives the fixtures and the matrix
  now declares two CRAM cases (see [The CRAM group](#the-cram-group-188)), so this moves to
  the list below **when a run has actually taken them** — not before.
* A CRAM whose reference is unavailable, which is the ordinary externally-referenced CRAM a
  diagnostic lab sends. The derived fixtures are `no_ref=1` and need no reference, so they
  cannot exercise that failure mode at all.

No longer true, and the change is run 4's:

* Cohort mode (`vntyper cohort`) was uncovered in runs 1–3, which is why `50d7968`
  produced no delta in run 3. Run 4 covers it with four cases. **Cohort *sample ordering*
  remains uncovered** even in run 4, because the harness normalises it away before
  comparing.

Covered only in one direction, run 4:

* #185's raise is exercised only negatively — no cohort case is missing a gate column.
* #184's boundary is almost certainly not reached by any cohort row, so its 0/59 on
  `kestrel_result` is weak evidence rather than confirmation.
* No run detects two changes whose deltas cancel exactly on the same field of the same
  sample. Run 4 is a single tip run, not one run per genotype-affecting commit.

---

## Run 5 — `4fd638a` → `9816f86`, after the adversarial review

Run 5 gates the nine commits that answer an adversarial review of PR #199 (Codex
`gpt-5.6-sol`, `xhigh`, five scoped read-only lanes; 47 findings, 6 Critical). It is the
first run taken with harness **1.1.0**, and the first on this project whose candidate side
can prove which revision it executed.

- Before: `4fd638a` (v2.0.6), **revision not recorded** — see the caveat below
- After: `9816f867c28f` on `fix/issue-181-197-followups`, **clean**, recorded by the
  harness and verified by `compare --expect-after-sha --require-clean`
- 65 runs per side, package resolution verified on every one, both sides
- Verdict: **DELTAS**, in two classes, both fully attributed

### Every genotype artefact is unchanged

| Compared | Cases with a delta | Cases compared |
| --- | --- | --- |
| `kestrel_result` | **0** | 59 |
| `kestrel_pre_result` | **0** | 59 |
| `advntr_result` | **0** | 3 |
| `coverage_summary` | **0** | 59 |
| `report_tables` | **0** | 59 |
| `screening_summary` | **0** | 59 |
| `cross_match_summary` | **0** | 3 |
| `exit_code` | **0** | 65 |
| `pipeline_steps` / `pipeline_step_records` | **0** | 61 |
| `cohort_category_counts` / `cohort_category_totals` / `cohort_tables` | **0** | 3 |
| `pseudonymization_table` | **0** | 1 |
| `cohort_output_files` | **0** | 4 |

### Run 5's delta 1 — the duplicate kestrel help flag, again (`2873ad3`)

`executed_commands` differs on 61 of 61 per-sample cases. The whole of it is one
substitution, measured rather than assumed: normalising the two per-side roots the way the
gate does and diffing all 65 cases gives **122 differing lines across 61 cases — exactly
two per case, one removed and one added — and nothing else**:

```
- java -jar vntyper/dependencies/kestrel/kestrel.jar -h -jar vntyper/dependencies/kestrel/kestrel.jar -h
+ java -jar vntyper/dependencies/kestrel/kestrel.jar -h
```

That is D6. The four cases without the delta are the cohort cases, which never invoke the
version probe. Zero unattributed command lines.

### Run 5's delta 2 — the leaked working columns (`90f61fa`)

The six cohort export artefacts differ on 3 of 3 cohort cases. Set-wise, on both algorithms:

| Export | Before | After | Removed | Added |
| --- | --- | --- | --- | --- |
| `cohort_kestrel.csv` | 31 columns | 29 | `__row_result`, `__unified` | none |
| `cohort_advntr.csv` | 15 columns | 13 | `__row_result`, `__unified` | none |

Exactly the two renderer-created working columns, on both. **No legitimate column was
dropped and none was added**; the projection also reorders columns to lead with the
display set. `cohort_category_counts`, `cohort_category_totals` and `cohort_tables` are
unchanged, so the categorisation those exports feed is unaffected.

### What run 5 does not attest

Everything under *What run 4 does not attest* still applies — it is one run at the tip,
not one per commit, and cannot exclude two changes producing offsetting deltas that cancel
on the same field of the same sample. In addition:

**It cannot attest the adVNTR signed-frame fix (`ad515c6`), which is the most
consequential commit it covers.** `advntr_result` shows 0 deltas on 3 of 3 adVNTR cases —
and that is *silence, not evidence*. A verdict changes only for a **mixed** state
(`Insertion_len >= 1` and `Deletion_length >= 1`) with Δ % 3 == 2, and the cohort contains
no mixed adVNTR state at all: `dfc3` is `D17_2&D18_2&D19_2&D20_2&D21_2`, a pure 5-base
deletion. The evidence for that commit is its 52,511-probe differential sweep — 9,782
states change verdict, every one with Δ % 3 == 2, nothing newly reported, 0 of 23,064 pure
states moved, all now hard failure conditions of the sweep — plus three states in
`advntr_config.json`'s `Polymorphic_Call` list. Not this run.

**The baseline's revision is not recorded.** The before side is a `git archive` extraction
of `4fd638a`, because the shared main worktree belongs to another checkout. An extraction
has no `.git`, so harness 1.1.0 logged a warning and the baseline SHA remains operator
record, exactly as described under *attestation* above. The candidate side **is** recorded
and was verified. The next run taken from two real checkouts can pin both.

**The baseline shares the candidate's `reference/` tree.** The extraction carries tracked
files only, and the adVNTR VNTR database is an installed artefact resolved at the relative
path `reference/vntr_db_advntr/<assembly>_muc1.db` (trap 7: paths are relative to the
process CWD, and each side runs with cwd set to its own tree). The first attempt failed
because of this, and the new expectation check named all six affected cases rather than
comparing them as `absent_both` — which is what the pre-1.1.0 harness would have done.
Sharing one reference tree is sound rather than expedient: `git diff 4fd638a..HEAD --
reference/` is empty, and `reference/**` is a base-image content-hash input that must be
identical on both sides.

**Run 5 ran no CRAM case.** That is a fact about run 5 and does not change: its 65 runs per
side are the 58-case BAM matrix plus 3 probes plus 4 cohort cases, and `175011e` is
therefore attested by the measurements in its own commit message and by a BAM-versus-CRAM
equivalence run, not by this run of this gate.

What *has* changed since run 5 is the matrix, not run 5. #188's fixtures exist —
`make cram-fixtures` derives a verified CRAM beside every cohort BAM — and the matrix now
declares two CRAM cases, `b178_hg19_cram` and `7a61_hg38_ensembl_cram`, both non-fast, both
counted by the drift check (see [The CRAM group](#the-cram-group-188)). **The next run will
cover the CRAM path; run 5 did not.** Until that run is taken and written up here, nothing
on this page attests the CRAM branch, and a reader should treat the two cases as declared
rather than as measured.

## Run 6 — `cb593b6` → `b27ff9c`, milestone 2

Harness `1.2.0`. Both sides clean, package resolution verified on every run, marker
`vntyper.scripts.coverage_qc` expected absent on `before` and present on `after`.

**The first run over the full 60-case matrix.** The harness refused the first launch —
`tests/data` derived 58 cases against the 60 the contract records, because the CRAM
fixtures were absent. Rather than pass `--allow-matrix-drift`, the fixtures were generated
(`scripts/make_cram_fixtures.py`, 50 derived and verified lossless) and the run relaunched.
The refusal is the harness working: *"a reduced run earns the same IDENTICAL verdict as a
full one, which is how a shrinking gate stays invisible."*

| Compared | Cases with a delta | Cases compared |
| --- | --- | --- |
| `coverage_summary` | 61 | 61 |
| `kestrel_pre_result` | 61 | 61 |
| `report_tables` | 61 | 61 |
| `executed_commands` | 61 | 63 |
| `kestrel_result` | 50 | 61 |
| `cohort_kestrel_{csv,tsv,json}` | 3 | 3 |
| `cohort_stats_{csv,tsv,json}` | 3 | 3 |
| `cohort_tables` | 3 | 3 |
| `cohort_output_files` | 3 | 4 |
| **`advntr_result`** | **0** | 3 |
| **`screening_summary`** | **0** | 61 |
| **`exit_code`** | **0** | 67 |
| **`pipeline_step_records`** | **0** | 63 |
| `cross_match_summary`, `cohort_advntr_*`, `cohort_category_*`, `pseudonymization_table`, `pipeline_steps` | 0 | — |

### The verdict, measured rather than argued

Across every case and every artefact, exactly **two** columns were added — `flag_filter_pass`
(117 occurrences) and `coverage_qc` (61) — and the **only** column whose cells changed was
`flag_filter_pass`, the new one:

```
=== columns ADDED across the whole run ===       === columns whose CELLS changed ===
   117x  flag_filter_pass                          46063x  flag_filter_pass
    61x  coverage_qc

=== cells changed on a column that was NOT newly added ===
  NONE - every changed cell is in a column this PR adds

=== any genotype field touched? ===
  NONE
```

`coverage_summary` is compared without a key, so its value changes surface as a row
removed plus a row added rather than as `cells_changed`. One case, in full:

```
BEFORE: mean 164.53  median 117.00  min 1  stdev 174.66  uncovered 454  pct 10.09
AFTER : mean 147.94  median 107.00  min 0  stdev 172.87  uncovered 454  pct 10.09  coverage_qc PASS
```

`uncovered_bases` and `percent_uncovered` are **unchanged**, which is the point: the new
zero-count reproduces the old subtraction exactly, so #171 corrects the four statistics that
were wrong without disturbing the one that was right. `min` moves 1 → 0, the true minimum of
a region with 454 uncovered bases.

**The closed-form identity `mean_new == mean_old · (1 − pct_old/100)` holds on 61 of 61
cases.** That is a complete independent check of #171 over the whole cohort, and it is the
identity #171 itself proposed for reconciling historical output.

### Attribution

| Delta | Attributed to |
| --- | --- |
| `coverage_summary` values (4 fields) + `coverage_qc` column | #171, #172 |
| `executed_commands` (`-a` on the depth command) | #171 |
| `kestrel_result` / `kestrel_pre_result` `columns_added: flag_filter_pass` | #174 |
| `report_tables` (new Coverage QC row) | #172 |
| `cohort_kestrel_*` gaining `flag_filter_pass` | #174 — `cohort_exports.py:14`, "nothing here strips columns" |
| `cohort_stats_*`, `cohort_output_files` | #172, the new export |
| `cohort_tables` (`cov_coverage_qc`, and the `2.0.7 → 2.0.8` version string) | #172 and the release. The normaliser is anchored to `VNtyper Version: ` and does not touch the bare version the statistics table carries. **Not** normalised away: attributing a difference is honest, teaching the gate to stop seeing it is not. |

`kestrel_result` is 50 of 61 rather than 61 of 61 because a negative call writes the
10-column sentinel, which carries no gate columns at all; only the 50 frames that reach the
final filter gain one.

### What this run does NOT attest

Three things, stated because a clean row above invites the opposite reading.

1. **It does not show #171 is genotype-neutral.** The mechanism is real — the corrected mean
   feeds `downsample_bam_if_needed` at `--advntr-max-coverage 300`, which both the gate and
   `docker/app/tasks.py:215` use. `advntr_result` shows 0 deltas because all three configured
   adVNTR cases are **fully covered**, so the corrected mean equals the old one exactly:

   ```
   a5c1_hg19  covered=1501/1501 zeros=0  old=1258.7215  new=1258.7215  IDENTICAL
   b178_hg19  covered=1501/1501 zeros=0  old= 878.2065  new= 878.2065  IDENTICAL
   dfc3_hg19  covered=1501/1501 zeros=0  old=2889.0286  new=2889.0286  IDENTICAL
   ```

   That is a property of `tests/data`, not of the change. The audit cohort behind #171 found
   1585 of 8215 samples with an inflated mean; none is in the local data.

2. **It does not exercise #174's exclusion.** Every positive call in this cohort is
   unflagged, so no artifact row was ever removed — `kestrel_result` shows a column addition
   and no row removals. The exclusion is covered by unit tests, not by this run.

3. **It does not see `output.bed`.** No BED artefact is collected, and the report reader
   extracts summary boxes and literal tables rather than the IGV payload. #203's coordinate
   fix is attested by `tests/unit/test_generate_bed_file.py` alone.
