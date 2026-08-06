# Spec: #181–#197 follow-ups, branch coverage, and the measured quality gaps

Status: draft for review
Date: 2026-08-06
Base: `test/issue-179-test-strategy` @ `9104f64` (PR #198, open, green)
Successor branch: `fix/issue-181-197-followups`, from `main` after #198 merges

---

## 1. What this is

PR #198 closed #179 by raising unit coverage from 25.68% to 70.57%, hardening
`docker/app`, and fixing ~20 defects. Seventeen findings it deliberately did **not**
fix were filed as #181–#197, each pinned by a characterisation test so the behaviour
could not drift while a decision was pending.

@hassansaei has now answered nine of them. This spec turns those answers into work,
closes the eight that needed no domain input, and attacks the three measured quality
gaps that #198 recorded but left open.

**Governing constraint.** No value in `vntyper/scripts/kestrel_config.json` changes,
and no confidence, motif or frameshift decision logic changes, without an explicit
decision from @hassansaei recorded in the issue. Every such change below quotes him.
Anything he has not decided stays characterised.

---

## 2. Evidence base

All figures measured at `9104f64` in this session, not carried over from the issues.
Where an issue's own number is stale, both are given.

### 2.1 Coverage

```
pytest -m unit tests/unit --cov --cov-branch
5198 statements, 1530 missing      -> line coverage 70.5656%
1502 branch exits, 640 never taken -> branch-only    57.39%
4530 of 6700 units covered         -> branch-inclusive 67.6119%
fail_under = 70
```

| | #196 says | measured at `9104f64` |
| --- | ---: | ---: |
| branch-inclusive total | 63.80% | **67.6119%** |
| against a floor of | 66 | **70** |
| units needed to clear the floor | 144 | **160** |

### 2.2 Where the uncovered units are

Missing units = missing statements + untaken branch exits.

| Module | Missing units | Branch-incl. | LOC |
| --- | ---: | ---: | ---: |
| `install_references.py` | 426 | 26.0% | 901 |
| `cohort_summary.py` | 346 | 37.8% | 911 |
| `kestrel_genotyping.py` | 187 | 41.9% | 843 |
| `utils.py` | 147 | 38.8% | 377 |
| `fastq_bam_processing.py` | 131 | 60.7% | 612 |
| `pipeline.py` | 122 | 68.5% | 721 |
| `docker/app/tasks.py` | 98 | 61.4% | 447 |
| `generate_report.py` | 96 | 64.8% | 574 |
| **`cross_match.py`** | **96** | **9.4%** | **190** |
| `motif_processing.py` | 68 | 67.1% | 463 |
| `extract_unmapped_from_offset.py` | 39 | 17.0% | 47 stmts |
| `file_processing.py` | 34 | 10.5% | 38 stmts |
| `variant_parsing.py` | 29 | 53.2% | 62 stmts |

`docker/app/tasks.py` is **63.0% line coverage**, not the 42% carried in the brief.

### 2.3 Mutation survivors in `motif_processing.py`, mapped to functions

The artefact lists 46 genuine survivors as flat line numbers. Mapped onto the module's
seven functions they are not evenly spread:

| Function | Lines | Survivors | Character |
| --- | --- | ---: | --- |
| `preprocessing_insertion` | 60–84 | 6 | column plumbing |
| `preprocessing_deletion` | 87–111 | 6 | column plumbing, **byte-identical twin of the above** |
| `_apply_uniform_filtering_right_motif` | 145–216 | 3 | the branch #186 says stays off |
| `_prioritize_frameshift_and_dedupe` | 217–281 | 3 | dedupe |
| `motif_correction_and_annotation` | 282–463 | **28** | 182-line function fusing config read, plumbing, split and filtering |

---

## 3. The seventeen issues

### 3.1 Decided, no behaviour change — close with a documented specification

Each of these lands as a docs/test-only commit. None can move a genotype; none needs
the golden-cohort gate.

**Scope rule, added after review: upgrade the individual assertions the issue decides,
never a section banner or a module docstring.** The first draft replaced whole-file
headers, which would have ratified behaviour nobody decided — `test_scoring.py`'s
banner also covers the sample-field parsing tests, `test_advntr_frameshift_filter.py`
also characterises parsing columns and import-time config, and
`test_confidence_boundaries.py`'s header covers the entire 54-cell matrix **including
the cell #184 changes**. Each upgraded test carries its own docstring naming the issue
and quoting the decision; everything else in those files stays characterisation.

| # | @hassansaei's decision (quoted) | What lands |
| --- | --- | --- |
| **#181** | *"this was a considered MUC1-specific choice, not accidental asymmetry… only the +1 frame yields the disease-associated MUC1-fs product. The other frame has not been established as pathogenic in patients"* | `tests/unit/test_scoring.py:227` "Characterisation tests (issue #179)" → **specification**, quoting him and citing #181. The biology into `extract_frameshifts`'s docstring and `docs/pipeline/scoring-and-confidence.md`. Rule unchanged. |
| **#182** | *"keep the same 3n+1 / 3n+2 rule for adVNTR as for Kestrel (#181). This is intentional shared convention… We should keep the filtering logic harmonized"* | Same upgrade in `tests/unit/test_advntr_frameshift_filter.py`. **Plus a new cross-module test** asserting Kestrel's `condition_insertion`/`condition_deletion` and adVNTR's `ins_frame`/`del_frame` encode the same convention, so the two cannot drift apart silently — which is the specific risk his "harmonized" answer creates. |
| **#183** | *"Keep the current 2.x last-wins logic — do not restore the absolute region-depth ≤200 cap from 1.3… the intentional behaviour going forward is the 2.x sequential assignment as implemented today"* | `tests/unit/test_confidence_boundaries.py:12` currently reads "Everything here is **characterisation**… do not read it as a claim". That becomes specification citing #183. Verify `docs/pipeline/scoring-and-confidence.md` already agrees (it was corrected on #198). |
| **#186** | *"Do not repopulate `motifs_for_alt_gg` with `["X"]`… Do not change the already-implemented v2.x motif logic on this point. Leave `use_uniform_filtering` off"* | **No config change, no logic change.** One addition: a fail-fast guard that raises if `use_uniform_filtering` is true while `motifs_for_alt_gg` is empty, naming the dupC deletion as the reason. See §3.4 for why this is labelled an inference. |
| **#187** | *"SHARK is sequence-based, not coordinate-based; keep one MUC1 region FASTA; document that assembly does not change SHARK; optionally warn or ignore `reference_assembly`… rather than building a second FASTA"* | Keep the parameter (API compatibility). Docstring and `docs/` state that it does not select a region. One `logger.warning` when a non-hg19 assembly is passed. **No hg38 region FASTA.** `tests/unit/test_shark_filtering.py`'s `TestReferenceAssemblyIsAccceptedAndIgnored` upgrades from "characterisation of a live bug" to specification. |

### 3.2 Decided, behaviour changes — golden-cohort gated

Each lands as its own individually revertible commit, all in the final phase, with a
single golden-cohort run at the tip.

#### #184 — mid-band `Depth_Score` demotion

> *"Any variant with Depth_Score between 0.00469 and 0.00515 (inclusive) must be
> Low_Precision, even when alternate depth is ≥21 (or higher). That mid-range
> Depth_Score demotion from 1.3 is still the desired behaviour… either by making this
> rule take precedence, or by excluding Depth_Score ≤ 0.00515 from the High tiers."*

**His answer and the code disagree only at one point, and the spec records that rather
than picking a side.** In `confidence_assignment.py` the assignments are last-wins and
`cond6` (`low < DS < high` → `Low_Precision`) is the **last** one, at line 154. So the
interior of the band is already demoted at every alt depth. `cond1` covers `DS == low`
and nothing overwrites it. The single input where the code contradicts his sentence is:

```
Depth_Score == 0.00515 exactly, alt_depth >= 21
  -> cond5 (alt in [21,100), DS >= high)  -> High_Precision
  -> cond2 (alt >= 100,       DS >= high)  -> High_Precision*
```

`cond3` therefore looks inert because `cond6` is doing its job, not because the intent
was lost.

**His second suggested implementation is unsafe as written.** Changing `cond2`/`cond5`
from `>= high` to `> high` alone leaves `alt_depth >= 100, DS == 0.00515` matching no
condition at all, so it falls through to `NEGATIVE_LABEL` — turning a reported call
into a non-call. That is a strictly worse outcome than the bug.

**Change.** Move the mid-band `Low_Precision` assignment to last and make it inclusive
at both ends:

```python
cond_midband = df["Depth_Score"].between(low_threshold, high_threshold, inclusive="both")
# ... applied after cond2 and cond5
```

`cond3` is retained as the expression that names the intent — he wrote *"do not remove
this intent"*.

**Blast radius, enumerated rather than estimated.** `Depth_Score` is `alt / region`, a
quotient of two integers, and `0.00515` is exactly `103/20000`. Enumerating every
integer pair with `region <= 200000`:

```
pairs where alt/region == 0.00515 exactly: 10
first: (103, 20000), (206, 40000), (309, 60000), (412, 80000), (515, 100000), ...
minimum alt depth: 103
pairs with alt depth in [21, 100): none
```

Three consequences the first draft of this spec missed:

1. **`alt` must be a multiple of 103, so `alt >= 103 > alt_mid_high (100)`.** `cond5`
   (`alt in [21, 100)`) is therefore **unreachable at the boundary**. Only `cond2` can
   fire there, so the only label that moves is `High_Precision*` → `Low_Precision`.
   `High_Precision` is never produced at `DS == high` at all.
2. `region` must be a multiple of 20000, i.e. an active region of at least 20 000 reads.
3. @hassansaei's "even when alternate depth is ≥21" describes a case that is arithmetically
   unreachable; the realisable form of his concern starts at alt depth 103.

**Two corrections from review.** First, "essentially unreachable" is too strong and is
withdrawn: `(103, 20000)` is an ordinary integer pair, and a MUC1 VNTR active region of
20 000 reads is not exotic. The right statement is *narrow and enumerable*, not
*unreachable*. Second, and more important, **the label is not the only thing that
moves.** `select_single_best_variant` sorts by confidence first
(`kestrel_genotyping.py:756`, priority order documented at `:680`), so demoting one row
from `High_Precision*` to `Low_Precision` can change **which variant is reported** when
a sample carries more than one candidate. The blast radius is therefore "a boundary row
and, transitively, the selection among its siblings" — still narrow, but not cosmetic.

The golden cohort is expected to show zero delta, and this spec says in advance that a
PASS there is weak evidence. The load-bearing evidence is a boundary table test over the
enumerated pairs above, plus `0.00469`, plus their float neighbours, **plus a
multi-candidate selection test** proving the demotion changes which row wins.

**Reporting obligation.** The finding that the band interior was already correct is
posted back to #184 so his answer sits on the record beside what the code does.

#### #185 — missing gate column raises

> *"Prefer fail loud: a missing required gate column should raise (abort the run), not
> be skipped… That is not acceptable for this pipeline… If a softer path is ever needed
> for empty early-return frames, that should be an explicit, documented empty-result
> path — not silent omission of safety gates on a non-empty candidate table."*

**Change.** In `filter_final_dataframe` (`kestrel_genotyping.py:777`), a filter column
absent from a non-empty frame raises `ValueError` after `logger.error`, following the
repo's no-custom-exceptions convention.

**Blast radius, established by reading the call graph rather than asserting it.**
`filter_final_dataframe` has exactly one caller, `kestrel_genotyping.py:594`, and the
six stages before it each end in `if df.empty: return df` (lines 550, 555, 560, 565,
573, 578). A frame reaching line 594 is therefore non-empty and has traversed every
stage, each of which adds its column unconditionally on a non-empty frame. His
empty-frame carve-out is already structurally satisfied: an empty frame never reaches
this function.

**Named risk, must be tested end to end.** `run_pipeline` uses `except Exception` at
stage boundaries by convention (AGENTS.md code style). A raise swallowed there and
converted into a silent no-result would be *worse* than the fail-open it replaces. The
acceptance test asserts the process exits 1 and the message reaches the log — not
merely that the function raises.

#### #192 — adVNTR `Insertion_len` for compound states

> *"use the **sum** of inserted lengths and the **sum** of deleted lengths when
> computing the net length that feeds the pathogenic frameshift filter… Example:
> `I9_2_A_LEN9&I50_2_A_LEN3` to Insertion_len = 9 + 3 = 12 (not first-LEN-only)…
> Do not keep first-LEN-wins as the defined semantics."*

**Correction: the current behaviour is not "first-LEN-wins".** Issue #192's text, the
#198 PR body and the first draft of this spec all say it is. Measured by running
`advntr_genotyping.py:163–169` directly:

| State | `Insertion_len` today |
| --- | ---: |
| `I22_2_G_LEN1` | 1 |
| `I9_2_A_LEN9&I50_2_A_LEN3` | **0** — not 9 |
| `I9_2_A_LEN2&D50_2` | **0** — not 2 |
| `I9_2_A_LEN9&` | **0** — not 9 |
| `I50_2&I9_2_A_LEN3` | 3 — the `LEN` is terminal, so it parses |
| `I50_2&D9_2&I80_2_A_LEN7` | 7 — likewise |

The greedy extract takes everything from the first `LEN` to the end, and the bounded
split then yields the remainder. **The value is zero exactly when any text follows the
first `LEN` token** — not, as the first draft of this spec claimed, for every compound
state. A compound whose only `LEN` is terminal parses correctly.

Where it does collapse to zero, `frame` becomes `0` for a pure-insertion compound, and
`0` is in neither `ins_frame` (3n+1) nor `del_frame` (3n+2), so **those states are
dropped in silence**. That is a materially worse defect than the issue describes — the
issue calls the behaviour "first-LEN-wins", which is true of no input at all — and it is
the real reason the change matters.

**Expected-difference oracle for the sweep, stated in these terms:** a state differs
before/after **iff material follows its first `LEN` token** — a second `LEN`, a further
`&` part, or a trailing `&`. States with no `LEN`, a single terminal `LEN`, or a `LEN`
only in the last part are unchanged. Anything else differing is a regression.

**Change.** Replace the greedy-extract-then-bounded-split
(`advntr_genotyping.py:163–169`, duplicated at 202–208) with a sum over all `LEN\d+`
tokens in the `Variant` string. His documented fallback for malformed strings is
adopted verbatim: if the summed lengths satisfy the frameshift rule, report the
variant as adVNTR emitted it; do not synthesise a simplified single event.

**Direction of the change, stated per path because they differ:**

- **Insertion path** (`Deletion_length == 0`): for the *affected* states — those whose
  value collapses to zero — `frame == 0` today and they are always dropped. Summation can
  only make those **appear**. States that already parse (terminal `LEN`) are unchanged.
  Monotone-additive on this path; no insertion call is lost.
- **Deletion path** (`Deletion_length >= 1`): `frame = abs(sum - Deletion_length)` can
  move either way. `I9_2_A_LEN3&D50_2&D51_2` is kept today (`frame = |0-2| = 2`) and
  **dropped** after (`frame = |3-2| = 1`). A call *can* be lost here, and the plan
  asserts that case explicitly rather than discovering it later.

**His worked example does not demonstrate the change.** `LEN9 + LEN3 = 12`, and
`12 % 3 == 0`, so it is filtered under both semantics. A case that genuinely moves is
`I9_2_A_LEN9&I50_2_A_LEN1` (sum 10, `10 % 3 == 1`, kept after and dropped before).

**The golden cohort cannot exercise it.** The only compound state in the cohort is
`example_dfc3`'s `D17_2&D18_2&D19_2&D20_2&D21_2`, which contains no `LEN` token at all;
`Insertion_len` is 0 under both semantics. A PASS from the gate is therefore *not*
evidence for this change. The required evidence is a **differential sweep** against the
pre-change parser over generated State strings — the same instrument used to clear
`d144505` (2380 probes, 572 non-crashing inputs byte-identical) — with an
expected-difference oracle stated in advance: differences only on states with two or
more `LEN` tokens, or one `LEN` plus a further `&` part.

**Separately: the existing attestation over-claims here.**
`docs/development/golden-cohort-gate.md:121` and `:219` report comparing
`Insertion_len` across runs 2 and 3. That column is not in the adVNTR output schema —
`final_columns` at `advntr_genotyping.py:365` is `VID, Variant,
NumberOfSupportingReads, MeanCoverage, Pvalue, RU, POS, REF, ALT, Flag`, and
`base_columns` at `:432` is shorter still. The row-set comparison remains valid and
does detect survival changes; the `Insertion_len` claim specifically is unsupported and
must be struck from that page rather than repeated. The spec states this explicitly so no attestation claims more than
it measured.

**Open sub-question to record, not to resolve unilaterally.** `Deletion_length` is
`Variant.str.count("D")`, i.e. a count of deletion *events*, used arithmetically as a
length in bp. That is only "the sum of deleted lengths" if every adVNTR `D` event is a
single base. That reading is consistent with the State grammar and with his example,
but it is an inference. It is written into the commit message and posted to #192.

### 3.3 Decided, config edit, provably inert

**#197** — > *"Fall back to the 1.3 Depth_Score-only rule: sort by `Depth_Score`
> descending when marking duplicates… Do not use `Motifs` or `Motif` as a sort key…
> Leave `duplicate_flagging.enabled` as `false` in the shipped config (as now)."*

`duplicate_flagging.sort_by` becomes `[{"column": "Depth_Score", "ascending": false}]`.
`enabled` stays `false`, so the block never executes in production and no genotype can
move. This is the one `kestrel_config.json` edit in this spec and it is authorised by
the quotation above.

**Implementation detail beyond his decision, flagged as such.** A single-key
`sort_values` is not stable, so with `Depth_Score` alone the row that survives a tie is
non-deterministic. `mark_potential_duplicates` will pass `kind="stable"`, making the
outcome reproducible from input order without adding a sort key he ruled out.

**Scope of "inert", corrected after review.** The *config* edit is inert under the
shipped configuration, because `flagging.py:149` reads `sort_by` only inside the
`enabled` block. The *code* change is inert only under that same disabled
configuration: which tied row stays unflagged determines which row is preferred by
`select_single_best_variant` (`kestrel_genotyping.py:749` prefers unflagged), so with
the toggle on it is genotype-affecting. The claim is therefore "inert under the shipped
disabled configuration", never "universally inert", and a tied-`Depth_Score`
enabled-path test lands with it.

### 3.4 Guard added beyond the letter of a decision

**#186's guard is an inference and is labelled one.** His words describe the hazard —
*"This is a high risk issue. It can remove the true positive from results. Especially
dangerous if someone enables uniform filtering while `motifs_for_alt_gg` is still
`[]`"* — but prescribe no mitigation. The tree ships that armed combination with
nothing preventing it.

The guard raises at config load when `use_uniform_filtering` is true **and**
`motifs_for_alt_gg` is empty. It cannot fire on the shipped config (`false`), so it is
genotype-neutral by construction and needs no gate. The uniform branch and its 588-LOC
test file stay, because he wrote *"unless/until that branch is redesigned and
validated"* — which reads as keeping it.

### 3.5 Engineering, no domain answer needed

| # | Disposition |
| --- | --- |
| **#190** | **Already fixed** on the #198 branch at `50d7968`: `cohort_summary.py:657` passes `autoescape=select_autoescape(["html","xml"])`. Closes on merge. No work. |
| **#188** | The issue's three-production-file diagnosis is **stale** — `--cram` already exists (`cli_parser.py:92`), `cli_handlers.py:302` threads it, `pipeline.py:180/313/328` branches on it. But it is **not a one-line fix and not zero-risk**; see §3.5.1. Moves to the gated phase. |
| **#194** | Set `mypy_path = "docker"`, fix the 22 annotation problems in `docker/app`, gate `tests/`. The three preconditions are already recorded as comments in `pyproject.toml`. |
| **#195** | §4.2 below. |
| **#196** | §4.1 below. |

#### 3.5.1 #188 needs two more fixes than the issue implies

Switching `docker/app/tasks.py:143` from `--bam` to `--cram` moves CRAM uploads onto a
different processing branch, and two things downstream are not ready for it. Both were
found by the adversarial review and verified directly:

1. **The sample name would regress to the literal string `"sample"`.**
   `cli_handlers.py:276-281` derives it from `args.bam` or `args.fastq1` and **never**
   from `args.cram`. Today a CRAM reaches the pipeline through `--bam`, so it *does* get
   its file stem; after the fix every CRAM run without an explicit `--sample-name` would
   be named `sample`. **The one-line fix introduces this regression.** Both must land
   together.
2. **The generated index name is wrong for CRAM.** `tasks.py:104` falls back to
   `f"{bam_path}.bai"`, but `samtools index` on a CRAM writes `.crai`. The cleanup path
   then removes a file that does not exist and leaves the one that does.

**The golden cohort contains no CRAM input**, so it cannot attest this — its own page
says so under "What this gate does not cover". The evidence has to be an end-to-end
CRAM run through the worker path, asserting the flag, the derived sample name, the
index filename and the cleanup. If no CRAM fixture can be produced, #188 parks on the
same trigger as #191 rather than shipping on an argument.

### 3.6 Parked, each with a written trigger

| # | Why parked | Trigger to unpark |
| --- | --- | --- |
| **#193** + `install_references.py` quoting | Both files are base-image content-hash inputs; editing either turns the app Docker build red until a new base publishes. `9104f64` reverted the quoting for exactly this reason. | The next `docker-base.yml` publish. Land both in the same commit as the rebuild. |
| **#189** | Per-job auth changes the public HTTP contract and breaks the `vntyper online` subcommand. No API-compatibility decision exists. | A recorded decision on whether per-job access requires a credential distinct from the job id. |
| **#191** | Its own text names the blocker: replacing the process substitution with a plain pipe "could alter which reads end up in `unmapped.bam`", and **the cohort contains no CRAM sample**, so the change cannot be verified. | A CRAM sample added to `tests/data/` and to the golden-cohort matrix. |

---

## 4. Coverage and quality

### 4.1 Enable branch coverage (#196)

**160 more covered units clears the floor of 70.** #196 concludes that splitting
`cohort_summary.py` and `install_references.py` is "likely a prerequisite". **The fresh
measurement says it is not.** Five modules, every one under 400 LOC and none on the
oversized list, hold 345 uncovered units between them:

| Module | Missing units | Branch-incl. |
| --- | ---: | ---: |
| `utils.py` | 147 | 38.8% |
| `cross_match.py` | 96 | 9.4% |
| `extract_unmapped_from_offset.py` | 39 | 17.0% |
| `file_processing.py` | 34 | 10.5% |
| `variant_parsing.py` | 29 | 53.2% |
| `docker/app/tasks.py` | 98 | 61.4% |
| **total** | **443** | |

Capturing under half of that clears 160 with margin. The two 900-LOC splits remain
worth doing (§4.4) but are **not** on the critical path for #196, and the plan says so.

**Two honest qualifications.** First, 160 units reaches exactly `4690/6700 = 70.00%`,
which lets branch coverage be *enabled* against the existing floor of 70 — it does not
raise the floor. `scripts/coverage_gate.py:166` advises a higher integer only when
`int(total) > floor`, so **227** additional units are needed before the floor can go to
71. Both numbers belong in the plan; conflating them would make the ratchet look stuck.
Second, "345 units are available in five small modules" is an upper bound on what those
modules *contain*, not a claim that all of it is cheaply coverable. The plan therefore
re-measures after each module and treats the shortfall, if any, as the trigger to bring
a split forward — rather than asserting the prerequisite is gone before it is measured.

**Sequencing is forced.** The floor is a ratchet, so `branch = true` and the floor
raise land in one commit, *after* the tests, set to the figure
`scripts/coverage_gate.py` prints — never the rounded `TOTAL` column.

### 4.2 `motif_processing.py` mutation gaps (#195)

30.9% (21/68), 46 genuine survivors. Attacked by shape, not by line number:

**a. The twin preprocessors — 12 survivors, one shared function.**
`preprocessing_insertion` (60–84) and `preprocessing_deletion` (87–111) are identical
except the `Variant` string literal. All 12 of their survivors are column-plumbing
mechanics — the `#CHROM`→`Motifs` rename, the exact five-element drop list, the
last-column→`Sample` rename — that no test observes because the existing tests assert
on the merged output rather than on the column contract. Extract one
`_preprocess_vcf_frame(df, muc1_ref, variant_label)`, and test the contract directly:
the rename, the drop list membership, the last-column identification with a
deliberately misordered frame, and the merge key. Twelve mutants, one test module.

**b. `motif_correction_and_annotation` — 28 survivors in 182 lines.** The file is 463
LOC so AGENTS.md rule 2 is satisfied, but rule 3 applies to the *function*: it fuses
config read, column plumbing, the left/right positional split and the motif filtering
decisions. Extract the pure decision layer (`motif_decisions.py`) — left/right split at
`position_threshold`, the right-motif exclusions, the GG handling, the combined
exclusions — leaving the merge/annotate plumbing behind. The pure layer is table-testable;
that is what makes the 28 killable.

**c. `_prioritize_frameshift_and_dedupe` (217–281) — 3 survivors**: lines 263, 265, 268,
all `False`→`True`. Direct tests.

**d. The uniform branch (145–216) — 3 survivors (189, 204, 212) — are deferred, and
they are NOT equivalent.** The first draft proposed classifying them equivalent because
the shipped toggle is off. That is wrong by the artefact's own definition: an equivalent
mutant is one that *cannot change behaviour any test could legitimately observe*, and
`tests/unit/test_motif_filtering_issue_136.py:352,376` calls
`_apply_uniform_filtering_right_motif` **directly**, with `use_uniform_filtering=True`
and a non-empty allowlist. A test can kill these; none currently does.

They stay recorded as **genuine gaps, deferred with a reason** — killing them buys no
production regression protection while #186 keeps the branch off, and #186 leaves open
that the branch may be redesigned. Misclassifying them as equivalent would inflate the
adjusted score with an untrue claim, which is the specific failure mode the artefact's
"if one turns out to be wrong the entry should be deleted, not the score explained away"
paragraph exists to prevent.

**e. Line 332 is a live defect the mutation run points at and no issue covers.**
`working_df["Motifs"].str.count("-").max() != 1` is a column-wide gate: one malformed
motif ID sets `motif_filter_pass = False` for **every** row of the sample. This is W8
from #179's §3.1, never filed. The `or`→`and` mutant at 332 survives. Disposition: file
it as a new issue with a characterisation test; **do not fix it here** — it changes
which rows survive motif filtering and no decision exists.

### 4.3 `cross_match.py` — the largest untested wrong-answer surface

Not on any issue list. 190 LOC, **13.9% line coverage, 9.4% branch-inclusive, and no
test file exists**. It `eval()`s a `match_logic` rule string from config
(`cross_match.py:137`) — AGENTS.md trap 3, the class that silently disabled a flag for
months — and it produces the `"Cross-Match Variant Comparison"` step that
`generate_report.py` and `cohort_summary.py` consume by exact string match.

Its five functions (`determine_variant_type`, `compute_allele_change`,
`cross_match_variants`, `write_results_tsv`, `extract_results_from_pipeline_summary`)
are pure or near-pure. This is the highest coverage-per-hour target in the repository
and it is a wrong-answer surface, not a crash surface.

**Confirmed while writing this spec: the trap-3 fail-open is live here.** The `eval` at
`cross_match.py:137` is wrapped in `except Exception` at `:138-140`, which sets
`match = False` and logs at ERROR. A `match_logic` rule naming a column the records do
not carry therefore reports **"no match"** rather than failing — the same shape as the
`Poylmorhic_Call` typo that disabled a flag for months and the `RU == 7` rule that never
fired. Unlike `flagging.py:84` it at least logs at ERROR rather than WARNING.

Disposition: **characterise, do not fix.** Turning it fail-closed changes what the
cross-match step reports, and no decision exists. Pin it with a test named for what it
does, and file it. Two further behaviours to pin rather than assume, both of which
contradicted this spec's first draft: `determine_variant_type` returns `"Other"` (not
`"Substitution"`) when REF and ALT are the same length (`:48-49`), and
`cross_match_variants` returns `overall_match` as the **string** `"Yes"`/`"No"`, not a
boolean (`:146`).

### 4.4 Oversized files

`cohort_summary.py` is now **911** LOC (AGENTS.md records 856) and
`install_references.py` 901. Both are worth splitting per rule 3, and together they
hold 772 uncovered units — but they are **not** blocking #196. They land after branch
coverage is on, so the split's coverage effect is measured on the sharper instrument.

`docker/app/main.py` at 1151 LOC is the worst offender by size, but at 88.8%
branch-inclusive it is not a coverage risk; it is a maintainability one. Out of scope.

### 4.5 `docker/app/tasks.py`

63.0% line / 61.4% branch-inclusive, 98 missing units. Covered in the same phase as
#188, which touches it.

---

## 5. Things worth more than the list

1. **`cross_match.py`** (§4.3). Highest value per hour in the repo; on no issue.

2. **PROPOSED, NOT IMPLEMENTED — the `.get()` fallback divergence should be deleted, not classified.**

   **Held back after review.** The first draft justified this as an extension of #185's
   "prefer fail loud". The review's objection is correct and is accepted: #185 decided
   *gate columns*, not *calibration keys*, and this spec's own binding constraint says
   no confidence logic changes without an explicit decision. Stretching one decision to
   cover a different site is exactly the overreach this spec is supposed to catch in
   others. It is therefore **written up and posted as a question**, and implemented only
   once answered. The six mutants stay hand-classified in the meantime.

   The case, for the issue:
   `confidence_assignment.py` reads its six calibration constants as
   `thresholds.get("low", 0.2)` etc., where the shipped config supplies `0.00469` — a
   43× divergence. The mutation artefact classifies all six survivors as *equivalent*
   on the grounds that a partial config is unsupported (AGENTS.md trap 2). That is
   defensible, but it means six wrong calibration values still exist in the source
   where a future reader can trust them. Replacing `.get(key, WRONG)` with `config[key]`
   makes them stop existing and **deletes six mutants rather than classifying them**.
   Behaviour is identical under the shipped config, which supplies all six keys; it
   differs only for a partial config, which AGENTS.md trap 2 already documents as
   unsupported input that fails with `KeyError` elsewhere in the pipeline.

3. **AGENTS.md's LOC/coverage table is stale and its central claim no longer holds.**
   It states "every module over 650 lines is under 30% covered". Measured now:
   `generate_report.py` is 574 LOC (was 861, split on #198) at 64.8%;
   `cohort_summary.py` is 911 not 856; `kestrel_genotyping.py` is 843 at 41.9%;
   `motif_processing.py` is 463. A table that is wrong stops being trusted, and this
   one is load-bearing — it is the argument for rule 2. Refresh it from measurement.

   Four further stale statements, all found by the adversarial review and all in the
   same refresh:

   | Site | Says | Reality |
   | --- | --- | --- |
   | `docs/pipeline/scoring-and-confidence.md:68` | tier precedence is "unspecified" | #183 decided it; and this spec's first draft claimed the page already agreed |
   | `docs/pipeline/flagging.md:38` | the three-key duplicate sort | #197 replaces it with `Depth_Score` alone |
   | `tests/unit/test_flagging.py:400` | pins the three-key sort | must move with the config |
   | `AGENTS.md:151` | equivalent mutants "have not been hand-classified" | `mutation-testing.md:21` hand-classifies ten |
   | `AGENTS.md` trap 11 | `vntyper report` is broken | fixed on #198; remove or mark resolved |
   | `docs/development/golden-cohort-gate.md:121,219` | runs 2 and 3 compared `Insertion_len` | not in the adVNTR output schema (§3.2) — strike the claim |

4. **Four further defects, found while writing the implementation plan.** All are in
   modules this spec was already going to test, all are silent-failure surfaces, and
   **none is fixed here** — each changes what is reported and no decision exists. They
   are characterised, filed, and named so the coverage work is not mistaken for a
   correctness claim.

   | Site | Defect |
   | --- | --- |
   | `cross_match.py:138-140` | The `eval` of `match_logic` is wrapped in `except Exception`, so a rule naming a missing column reports **"no match"** instead of failing. Trap 3, live. |
   | `file_processing.py:28-30` | A row is kept only if **exactly one** of REF/ALT has length 1, so `REF="AC", ALT="ACGGG"` — a real 3-base insertion — is discarded before Kestrel post-processing ever sees it. |
   | `file_processing.py:60-63` | The insertion test is `len(ref) == 1 and len(alt) > 1` with an unconditional `else`, so a multi-base-REF insertion is routed to the **deletion** file, where the 3n+2 rule is applied instead of 3n+1. Under #181 those are not interchangeable. Currently masked by the row above. |
   | `extract_unmapped_from_offset.py:38-53` | `f.read(4)` at EOF returns `b""` and `int.from_bytes` gives `0`, so a truncated BAI parses as "no more references" and the unmapped scan starts at offset 0 rather than raising. |

5. **`tests/unit/test_marker_hygiene.py` fails under any filtered run, `-k` included.**
   Recorded so parallel execution does not rediscover it.

---

## 6. Sequencing and gating

| Phase | Contents | Gate | Genotype risk |
| --- | --- | --- | --- |
| **0** | Merge #198. Close #179, #190. | CI on #198 | none |
| **1** | #181, #182, #183, #186, #187, #194, #197. Refresh AGENTS.md and the four stale doc sites. File the line-332 and §5.4 issues. | `make check-all` | none |
| **2** | `cross_match`, `utils`, `file_processing`, `extract_unmapped_from_offset`, `variant_parsing`, `tasks.py` | `make check-all` | none |
| **3** | `branch = true`; floor raised to the printed figure (#196) | `make ci-local` | none |
| **4** | Split `cohort_summary.py`; `install_references.py` stays parked | `make check-all`; floor re-raise | none |
| **5** | The gated set, **one commit and one gate run each**: `motif_processing` refactor (a–c), #184, #185, #192, #188 | **golden cohort per commit**, plus the #192 differential sweep and a CRAM end-to-end for #188 | yes |
| **6** | Final attestation at the PR tip; close the issues | golden cohort at the tip; `make ci-local` | — |

**Three corrections to the first draft's sequencing, all from the adversarial review:**

1. **The gate attests one commit and nothing after it**, so a phase after the gate
   invalidates it. The first draft put the two file splits *after* the gated phase while
   asserting the gated phase was last. The splits move to phase 4 and the gated work is
   genuinely last.
2. **A single cumulative run cannot attribute a failure**, and two changes can produce
   offsetting deltas that a joint run reports as clean. Each genotype-affecting commit
   gets its own run — five runs at ~25 min is ~2 hours, which is affordable for the
   attribution it buys — followed by one final run at the PR tip.
3. **The `motif_processing` refactor is not genotype-neutral** and moves from phase 2 to
   the gated phase. It extracts the code that computes `motif_filter_pass`
   (`motif_processing.py:430`), which `filter_final_dataframe` requires
   (`kestrel_genotyping.py:806`). A unit oracle over one synthetic frame is weaker
   evidence than 58 real samples; "intended to be behaviour-preserving" is not evidence
   at all.

---

## 7. Acceptance criteria

1. Every one of #181–#197 is closed, or open with a written reason and a trigger
   condition (#189, #191, #193 only).
2. Every characterisation test whose issue @hassansaei has decided is upgraded to a
   specification test that **quotes him and cites the issue number**. No test docstring
   claims specification where only characterisation was measured.
3. `branch = true` is enabled and `fail_under` is raised to the exact figure
   `scripts/coverage_gate.py` prints. The floor is never lowered.
4. **One golden-cohort run per genotype-affecting commit**, each naming its own
   candidate SHA, plus a final run at the PR tip. Every run reports an empty diff on
   Kestrel variant set, `Confidence`, `Flag`, adVNTR genotype fields and `Flag`, and
   every run's record states what it did **not** exercise. No run is cited for a commit
   that landed after it.
5. A differential sweep substantiates #192 independently of the gate, against an
   expected-difference oracle written before the sweep runs.
5b. A CRAM end-to-end run substantiates #188, or #188 parks on the #191 trigger.
6. `make check-all` and `make ci-local` are green, with output shown.
7. `vntyper/scripts/kestrel_config.json` differs from its pre-branch state in exactly
   one key, `duplicate_flagging.sort_by`, authorised by the #197 quotation.

## 8. Non-goals

- Fixing `motif_processing.py:332` (W8). Characterised and filed; no decision exists.
- Fixing any of the four defects in §5.4. Every one changes what is reported —
  which variants reach post-processing, which side of the indel split they land on,
  where the unmapped scan starts, or whether a cross-match rule failure reads as
  "no match". All are characterised and filed. Fixing them is the next spec's work,
  after a decision.
- Populating `motifs_for_alt_gg` or enabling `use_uniform_filtering`. #186 forbids it.
- Restoring the region-depth cap. #183 forbids it.
- Building an hg38 SHARK region FASTA. #187 rules it out.
- Per-job authentication (#189), the CRAM process substitution (#191), the base-image
  items (#193). Parked with triggers.
- `docker/app/main.py`'s size.
