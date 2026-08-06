# Golden-cohort gate (#179)

Before-versus-after comparison of genotyping output across the whole local test cohort,
run to decide whether the deliberately behaviour-changing commits in #179 may ship.

The gate has been run three times. **A verdict attests one candidate commit and nothing
after it**, so read the run whose candidate matches the tree you are judging.

| Run | Candidate ("after") | Verdict | Attests |
| --- | --- | --- | --- |
| 1 | `7344c62` | PASS | the branch as it stood at `7344c62`, and nothing after it |
| 2 | `1792345` | PASS | the branch at `1792345`, including `d144505`, `c51052c`, `b4059ce`, `7e58eb8`, `2c92096`, and nothing after it |
| 3 | `8537a61` | PASS | the branch at `8537a61`, including `2ae28c5`, `2aa095a`, `50d7968`, `42c976a`, `4ce5639`, `97033d3`, `22e3d17`, `52a0ec9`, and nothing after it |

Run 3 is the attestation of record. It is the run to cite for a tree at `8537a61`; runs 1
and 2 are kept because they are the measurements the adjudications below were written
against, and because each one still attests its own candidate exactly.

None of the three runs attests "the branch tip" as a standing property — a tip moves, a
commit does not. Each run's candidate is named above and in its own result section; if the
branch has commits after `8537a61`, this page does not cover them and says so rather than
implying otherwise.

**Verdict: PASS, all three runs.** Every genotype field, every `Confidence` label and every
`Flag` is byte-identical between baseline and candidate, on every sample and every
assembly. The only differences are in how the HTML report *presents* an unchanged result,
and each one is attributable to a named commit and is a correction of a documented
reporting defect.

| | Run 1 | Run 2 | Run 3 |
| --- | --- | --- | --- |
| Baseline ("before") | `2fcc6e3` — merge-base with `main` | `2fcc6e3` | `2fcc6e3` |
| Candidate ("after") | `7344c62` | `1792345` | `8537a61` |
| Cases per side | 58 (50 BAM x assembly, 5 non-fast-mode, 3 adVNTR) | same 58, same matrix | same 58, same matrix |
| Runs total | 116, plus 6 deliberate-mismatch probes | 116, plus 6 probes | 116, plus 6 probes |
| Non-zero exits | 0 before, 0 after | 0 before, 0 after | 0 before, 0 after |
| Executed shell commands compared | no | no | **yes** — 480 per side |

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
the exit code. Run 3 adds the executed shell command strings, recorded at the
`subprocess` boundary.

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
and 3, all 61 logs per side. The before-side worktree carries its own tracked
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
| adVNTR variant set, genotype fields, `Insertion_len`, `Flag` | 0 / 3 |
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
rewrote. Its row is identical on both sides in every column, `Insertion_len` included.

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
| adVNTR variant set, genotype fields, `Insertion_len`, `Flag` | 0 / 3 |
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
* Cohort mode (`vntyper cohort`) and the report subcommand — and therefore `50d7968`,
  which changes only the cohort report.
* Reference installation, and therefore the aligner `index_command` quoting site
  `2ae28c5` touches in `install_references.py`: the references are already installed, so
  no gate run builds an index.
* A path containing a shell metacharacter. Run 3 shows `2ae28c5` changes no command in
  this cohort precisely because no cohort path needs quoting; the quoting itself is
  pinned by `tests/unit/test_shell_quoting.py`, not by this gate.
