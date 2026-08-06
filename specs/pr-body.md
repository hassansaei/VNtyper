# Follow-ups to #179: close #181–#188, #192, #194–#197

Branch `fix/issue-181-197-followups`, on top of `4fd638a` (v2.0.6), gated by
golden-cohort **run 5** (`4fd638a` → `9816f86`).

> **On the commit count and tip SHA:** they are deliberately not stated here. Three times
> in this PR a hardcoded count and tip went stale, because *a number describing a branch,
> written inside a commit on that branch, invalidates itself the moment it is written* —
> the fix commit is the (n+1)th. GitHub renders both live at the top of this page; read
> them there. This is the same defect class the PR is about, and the durable repair is to
> delete the self-referential number rather than keep correcting it.

The final commits answer an adversarial review of this PR — see
**[Adversarial review](#adversarial-review)** at the end, which includes a **merge
blocker that review found and this branch had introduced**.

This finishes the work PR #198 opened. #198 raised unit coverage from 25.68% to 70.57%
and, in doing so, filed seventeen issues describing everything it had found but not
fixed. This branch resolves thirteen of them, fixes eleven defects found while doing so,
commits the golden-cohort gate that had never been committed, and ratchets the coverage
floor from 70 to 80 on a strictly harder measurement.

Every genotype-affecting change is its own revertible commit with its own evidence. The
evidence table below states, for each one, **what does not attest it** — because for two
of them the golden cohort provably cannot.

---

## 1. What landed, by phase

### Phase 0 — prerequisite (already on `main`)

PR #198 merged as `2672dfa`; #179 and #190 closed.

### Phase 1 — one issue per lane (`e36f3ad..65e033c`)

| SHA | Issue | Change |
| --- | --- | --- |
| `99affb2` | #186 | `motif_processing.py` refuses `use_uniform_filtering` with an empty GG allowlist. Additive guard, 17 lines, cannot fire on the shipped config. |
| `88cb88e` | #181 | The +1-frame rule recorded as **specification** in `docs/pipeline/scoring-and-confidence.md` and in one upgraded test docstring. No executable change. |
| `513890e` | #182 | Kestrel and adVNTR asserted to share one frameshift convention, driving the real `advntr_processing_ins`/`_del` against `scoring.extract_frameshifts`. No executable line of `advntr_genotyping.py` changed. |
| `39fe0eb` | #187 | SHARK's `reference_assembly` documented as accepted-and-ignored, plus a warning on a non-hg19 value. |
| `6e7cda2` | #197 | `duplicate_flagging.sort_by` reduced to `Depth_Score` alone; `flagging.py` sorts with `kind="stable"`. One config key changed. |
| `451a4b7` | #183 | Last-wins tier precedence specified, with a 54-cell matrix and a docstring that names the future collision (#184) rather than leaving it implicit. |
| `65e033c` | #194 | `docker/app/` gated under mypy: 23 errors → 0, zero `# type: ignore`, zero relaxed settings. |

### Phase 2 — coverage of the modules that had none

| SHA | Module | Before → after (branch-inclusive) |
| --- | --- | --- |
| `dca508c` | `file_processing.py` | 0% → 100% (the module had never been imported by any unit test) |
| `813eaea` | `cross_match.py` | 13.9% / 9.4% → 100% / 100% |
| `47f0f63` | `extract_unmapped_from_offset.py` | 17.0% → 100% |
| `b497de8` | `variant_parsing.py` | 53.2% → 100%, written to kill three named surviving mutants |
| `ba5a247` + `5bb2463` | `utils.py` | 38.8% → 100% |
| `4819817` + `4ce1973` | `docker/app/tasks.py` | 61% → 98% |
| `a07de32` | — | `[tool.ruff] include` pinned so the bare `ruff format` agrees with the gate (ruff 0.16 widened discovery to Python fenced in Markdown). |

`5bb2463` is the most useful commit in this phase and the one to read if you read only
one: `utils.py` reached **100% branch coverage with a surviving mutant behind it**.
Dropping `cwd=cwd` at `utils.py:334` left the whole unit tier passing, byte-identical.
That is AGENTS.md trap 7 — every path in `config.json` is relative to the process CWD —
and 100% was the number that hid it.

### Phase 3 — the ratchet (`debee23`)

Branch coverage enabled (`branch = true`), floor 70 → 74, taken verbatim from what
`scripts/coverage_gate.py` printed rather than from the rounded `TOTAL` column.

`debee23` also ports a correction into `scripts/mutation_test.py` rather than into
`docs/development/mutation-testing.md`. That page is machine-generated: editing the
Markdown is silently reverted by the next `make mutation-render`, and
`tests/unit/test_mutation_test.py` asserted against the generator, so the test would have
kept passing the whole time. A correction that reverts on the next sweep is worse than no
correction. The round trip is proved: `make mutation-render` leaves the committed file
byte-identical.

### Phase 4 — splitting `cohort_summary.py` (`1c6c9d6`, `2184805`, `c6d967a`)

911 → 455 LOC across five focused modules: `cohort_rules` (162), `cohort_categories`
(134), `cohort_tables` (198), `cohort_inputs` (259), `cohort_exports` (108). Coverage of
the six together 38% → 95%; all five extracted modules at 100%; `cohort_summary.py`
itself at 84%, its remaining misses all matplotlib/Jinja2/filesystem. Floor 74 → 79.

The evidence is an end-to-end oracle whose fingerprint was recorded **before** any
extraction and never edited. Its limits are stated in §5.

### Phase 5 — the gated phase

| SHA | Change |
| --- | --- |
| `b89a4f3` | D4 — `extract_unmapped_from_offset` raises on a truncated BAI instead of scanning from offset 0 |
| `2873ad3` | D6 + D7 — `validate_bam_file` honours its documented `ValueError`; the kestrel help flag is no longer duplicated |
| `a4264f8` | D10 — `report_config.json` adVNTR rule 2 guarded on VID, matching rule 1 |
| `5ee1e4a` | D1 — `cross_match` fails **closed** when a `match_logic` rule cannot be evaluated |
| `cde3aca` | The golden-cohort harness, committed for the first time, with cohort-mode cases |
| `102c46f` | D2 + D3 — indels classified by length difference; Kestrel's allele contract asserted |
| `c489d49` | #185 — a missing required gate column raises instead of failing open |
| `2a267aa` | #192 — every `LEN` token in a compound adVNTR state is summed |
| `5f5222e` | #188 — `--cram` passed for CRAM uploads, plus the two regressions that alone would cause |
| `6f03ddf` | #184 — mid-band `Depth_Score` stays `Low_Precision`; D11 — the calibration constants are required, not defaulted |
| `90f61fa` | D8 + D9 — working columns no longer leak into exports; sample order is deterministic **for fixed directory inputs** (ZIP inputs are not: #205) |
| `11e2300` | #195 — motif decision layer extracted, twin preprocessors merged, D5 column-wide dash gate fixed per row |
| `ec67fff` | Floor 79 → 80 |

---

## 2. The golden-cohort gate

`cde3aca` commits it. This matters on its own: the gate is the instrument backing every
genotype claim on this project, runs 1–3 are documented in prose in
`docs/development/golden-cohort-gate.md`, and **the scripts that produced them lived in
`/tmp` and are gone**. The gate had to be reconstructed from a description each time, by
a different person, with no guarantee it was the same instrument. It also had never
covered cohort mode — which became load-bearing the moment `1c6c9d6` moved 456 lines the
gate could not see.

The committed harness derives its 58-case matrix and 3 probes from what is actually in
`tests/data/` rather than from a hardcoded list, so it cannot silently drift. The
`launcher` is the part that makes the rest worth anything: the `vntyper` console script
resolves through setuptools' editable finder, which is appended to `sys.meta_path` and
points at whichever worktree the editable install was made from irrespective of cwd — so
a baseline checkout can silently execute candidate code. The `probe` subcommand
reproduces that leak on demand and every run refuses to start (exit 97) unless its
resolved `vntyper.__file__` agrees with its side.

**One deliberate deviation from the plan, and its residual risk.** The plan required one
golden-cohort run per genotype-affecting commit (6+ runs, ~25 min each). One run at the
tip was taken instead, bisecting only on a delta: six unconditional runs pay 2.5 h for
attribution needed only in the failure case, and one-run-plus-bisect has the identical
worst case. The argument this does **not** answer is that two changes could produce
offsetting deltas a joint run reads as clean. The comparison is per-sample and per-field
across 58 cases, so cancellation would have to be exact on the same field of the same
sample — but the attestation claims joint attestation, not per-commit attribution, and
should be read that way.

The identity of the 5 non-`--fast-mode` cases is not recoverable from the published page —
only their per-assembly counts. The reconstruction here reproduces every published count
(50/5/3 = 58, and 20/9/8/7/7/7) with zero mismatches, but it is a reconstruction,
documented as such and overridable with `--non-fast-cases`.

---

## 3. Genotype-affecting commits, and what does not attest them

| Commit | What it changes | What attests it | **What that evidence does not cover** |
| --- | --- | --- | --- |
| `102c46f` | `filter_vcf` keeps indels with multi-base alleles on both sides; `filter_indel_vcf` routes by length difference; a new `_assert_kestrel_allele_contract` guard | `javap` on the vendored `kestrel.jar`: `VariantInsertion.getVcfRef()` ends in `Character.toString(char)` on all three return paths, so REF is always exactly one base and Kestrel 1.0.1 emits only 1-vs-1, 1-vs-N, N-vs-1. Both defects are therefore **latent**. The post-fix classification table is measured and every 1-vs-N and N-vs-1 row is identical to the pre-fix table: production blast radius zero. | Nothing attests behaviour under a different jar — the guarantee is a property of *this* jar and is enforced nowhere else in VNtyper (AGENTS.md trap 8). Equal-length multi-base rows (MNPs) move from **silently dropped to hard failure**; deliberate, but it is a silent-drop-to-raise transition. `filter_indel_vcf` no longer reads `config.json`. |
| `c489d49` (#185) | A missing gate column raises instead of being skipped | End-to-end test driving the real `run_pipeline` with a stage made to stop emitting `motif_filter_pass`, asserting `SystemExit`, code 1, and the message at ERROR — not merely reading the call graph. Before the fix the same test failed with "the run did not exit; it returned None". All five gates are parametrised. All five stages were read: every gate assignment is unconditional on each non-empty return path. | A cohort PASS proves only that no shipped stage omits a gate **on the cohort's inputs**. A cohort *difference* here would have been a finding, not a reason to weaken the raise. Web-visible: a job that previously produced a wrong result can now surface as failed. |
| `2a267aa` (#192) | `Insertion_len` is the sum of every `LEN` token, not the remainder after a greedy split | A **differential sweep**: 52,511 probes; 13,563 / 13,563 previously-parsing inputs byte-identical; 38,943 differing, every one predicted by an oracle stated *before* the sweep ran; 0 violations in the three unchanged classes; reported-call delta **+11,467 / −3,145** (restated at `ad515c6`: both sides of the comparison are now evaluated through the corrected, narrower filter). Induced-failure proof: tightening `LEN(\d+)` to `^LEN(\d+)` produced 21 failures and a non-zero exit. | **The golden cohort provably cannot attest this.** The cohort's only compound state, dfc3's `D17_2&D18_2&D19_2&D20_2&D21_2`, carries no `LEN` token, so `Insertion_len` is 0 under both semantics — measured, and confirmed byte-identical on both sides. The sweep covers a swept grammar, not observed adVNTR output. And `Deletion_length` is still `Variant.str.count("D")`, a count of *events* used arithmetically as a length in bp (unresolved, posted to #192, drafted as a new issue). |
| `5f5222e` (#188) | `--cram` for CRAM uploads; CRAM sample-name derivation; `.crai` cleanup | An end-to-end run through the real Celery worker with only Redis faked, on a CRAM built at test time into scratch: the command carries `--cram`, `pipeline_summary.json` records `input_files={"cram": …}`, the VCF sample column is the file stem, the `.crai` is produced and cleaned up, and `kestrel_result.tsv` is byte-identical to the BAM run of the same sample apart from the timestamp (93 / 9883 / 0.009410098148335525 / High_Precision). Two induced reverts: pre-#188 leaves the `.crai` on the volume; fix 1 alone yields a VCF sample column of `sample`. | **The golden cohort contains no CRAM input at all** — that is #191's trigger. The genotype-equality assertion is **hand-run, not in CI**: it needs samtools, reference data and a real pipeline subprocess, all barred from the unit tier. #188 does not park on #191, but it stays unattested in CI until an integration case exists. |
| `6f03ddf` (#184) | `cond6` becomes `cond_midband`, inclusive at both ends, applied last | A 54-cell boundary matrix; 116 enumerated integer depth pairs reaching 0.00515 for alt ≤ 12000, every one with alt a multiple of 103 and therefore > `alt_mid_high`, so only `High_Precision*` moves on integer depths; a selection test measured in both directions (before POS 67 / `High_Precision*`, after POS 68 / `High_Precision`). | **A golden-cohort PASS is weak evidence here and should not be read as strong.** Exact float equality with 0.00515 requires alt depth to be an exact multiple of 103, essentially unreachable on real data, so the cohort very likely contains no boundary row at all. The boundary table and the selection test are the load-bearing evidence. `cond3` is now dead in effect and is retained deliberately, as the expression that names the intent; its operands will surface as equivalent-mutant survivors in a future sweep. |
| `11e2300` (#195) | `motif_processing.py:332` dash gate applied per row; decision layer extracted to `motif_decisions.py`; twin preprocessors merged | An end-to-end oracle hashed against the pre-extraction code and passed unmodified. The oracle was strengthened before any code moved: a `<` → `<=` mutation survived the first draft, so a POS-60 boundary test was added and the hash re-derived. `motif_decisions.py` is registered in `scripts/mutation_test.py` TARGETS — a module absent from TARGETS is not mutated at all. | The expected cohort delta is none **because real motif IDs carry exactly one dash** — an expectation, not evidence, and the attestation must say so. The oracle is an equivalence proof for the *extraction*, not for the dash-gate fix. The fix also uncovered an opposite face nothing had recorded: a **dashless** row rides through on a sibling's dash and is reported as a passing call on a motif not in the reference. That direction is more dangerous than the one #179's audit recorded, and had no instrument before this branch. |
| `a4264f8` | `report_config.json` adVNTR rule 2 guarded on VID | `test_advntr_rule_symmetry.py` enumerates all five reachable combinations; genotype-neutral by construction. | Nothing for the cohort to show: the misreport was proved unreachable twice (the only producers of `VID == "Negative"` set `Flag = "Not applicable"`, which rule 2 excludes, and `add_flags` is never reached on that path). |
| `90f61fa` | `sample_categories` stops mutating the caller's frame; `discover_sample_directories` returns a sorted list | The oracle fingerprint does **not** move, and the reason is in the source rather than in luck: `cohort_tables` builds from the fixed `KESTREL_DISPLAY_COLUMNS`/`ADVNTR_DISPLAY_COLUMNS` tuples, so the leaked columns never reached the HTML the fingerprint measures. Determinism is tested **cross-process** at five `PYTHONHASHSEED` values — a single-process comparison cannot see a per-process hash seed, and five interpreters produced five distinct orders before the fix. | Report-affecting only; no genotype moves. The oracle exercises 3 of the 5 Kestrel verdict paths and 2 of 3 on the adVNTR side (see §5). |

---

## 4. Coverage

| Point | Figure | Floor |
| --- | --- | --- |
| Before #198 | 25.68% statement | — |
| Branch point (`4fd638a`) | 70.57% statement | 70 |
| After phase 3 (`debee23`) | 74.22% **branch-inclusive** | 74 |
| After the cohort split (`c6d967a`) | 79.02% branch-inclusive | 79 |
| `ec67fff` | 80.24% branch-inclusive | 80 |
| After the review fixes | **80.37% branch-inclusive** | **80** |

Statement-only coverage of the same suite is 80.77%, so the two are not interchangeable
and deleting `branch = true` would *raise* the reported total while covering strictly
less. The ratchet cannot catch that regression — the number moves the wrong way — so
`tests/unit/test_coverage_gate.py::test_branch_coverage_is_enabled` exists to.

80.24% meets the 80% target `COVERAGE_TARGET` records. The target stays a warning rather
than becoming the floor: they measure the same number for different purposes, and
collapsing them would leave no headroom above the gate.

Every floor raise is a separate, revertible commit, and every value was taken from what
`scripts/coverage_gate.py` printed, never from the rounded `TOTAL` column. The floor has
never been lowered.

Tests: 1744 at the branch point → **2290**.

---

## 5. Where the evidence is weakest

Read this section before approving anything in §3.

**The cohort oracle's real limit, found by mutation and not by argument.** Changing
`unify_kestrel_result` so `Low_Precision` maps to `Negative` instead of `Positive` — a
one-token change that would alter every cohort donut chart in the field — **passed the
oracle and passed the escaping suite.** The fixture never produces a bare `Low_Precision`
verdict. The oracle exercises 3 of 5 Kestrel verdict paths and 2 of 3 adVNTR paths. Only
`test_cohort_categories.py` caught it, and that is a *dependent* instrument written by the
same agent from the same reading, so a mis-transcription would propagate into test and
module together and both would agree.

What rescues the refactor is two facts held **together**, never the oracle alone:
(1) the reviewer reverted `cohort_summary.py` to its pre-split 911-line state and ran the
new oracle against it — 27 passed, so the fingerprint is what the *original* file
produces; and (2) an AST body comparison found every whole-function move byte-identical,
including the mutated function. A branch the oracle never exercises cannot have been
mis-transcribed if the body is unchanged. Residual risk is confined to the twelve
functions lifted from inline blocks, which carry no byte-identity guarantee; all twelve
were read line-by-line against their originals.

**A covered line is not a tested line, demonstrated on this branch's own work.**
See `5bb2463` above.

**`vntyper report`, `install_references.py` and `fastq_bam_processing.py` are still
largely untested** (26.0% and 60.7%). Nothing here changes that.

---

## 6. The fourteen defects

Eleven found during the coverage work were **fixed** on this branch, each as its own
revertible commit. Three remain open and are drafted as issues.

| # | Defect | Status |
| --- | --- | --- |
| D1 | `cross_match:138-140` — an `eval` failure on a `match_logic` rule reports "no match" instead of raising (AGENTS.md trap 3; same shape as the `Poylmorhic_Call` typo and the `RU == 7` rule) | **Fixed** `5ee1e4a` |
| D2 | `filter_vcf` kept a row only if *exactly one* of REF/ALT was one base, discarding real multi-base-REF insertions before Kestrel post-processing | **Fixed** `102c46f` |
| D3 | `filter_indel_vcf` routed such a row to the *deletion* file | **Fixed** `102c46f` (had to land with D2; fixing D2 alone would have exposed D3) |
| D4 | `extract_unmapped_from_offset` returned 0 on a truncated BAI while its own docstring promised to raise — **a violated documented contract** | **Fixed** `b89a4f3` |
| D5 | `motif_processing:332` used a column-wide `.max()`, so one malformed motif ID set `motif_filter_pass=False` for **every row of the sample** and the run reported Negative while exiting 0 (W8 from the #179 audit, never filed) | **Fixed** `11e2300` |
| D6 | `get_tool_versions` duplicated `-jar "<path>" -h`; argv was `['java','-jar',jar,'-h','-jar',jar,'-h']` | **Fixed** `2873ad3` |
| D7 | `validate_bam_file`'s docstring promised `ValueError`; the code raised `RuntimeError`, because `critical=True` made the `ValueError` branch unreachable by construction — **the second violated documented contract** | **Fixed** `2873ad3` |
| D8 | `__row_result` / `__unified` leaked into every cohort CSV/TSV/JSON export | **Fixed** `90f61fa` |
| D9 | Cohort sample order was non-deterministic (`Path.__hash__` delegates to a string hash, so `PYTHONHASHSEED` randomised it) | **Fixed** `90f61fa` |
| D10 | `report_config.json` adVNTR rule 2 lacked the VID guard rule 1 has | **Fixed** `a4264f8` |
| D11 | Six `.get()` calibration defaults in `confidence_assignment.py` encoded a **second, wrong calibration** — 0.4 where the shipped config says 0.00515, a factor of 78, reachable by dropping one key | **Fixed** `6f03ddf` |
| D12 | `vntyper cohort` over a directory with no `pipeline_summary.json` exits 0 and writes no report | **Open** — issue drafted |
| D13 | `scripts/` is in `RUFF_PATHS` but in **no** mypy target, so ~3200 lines added by this branch are type-checked by nothing | **Open** — issue drafted; see §8 |
| D14 | Every job leaves `<alignment>.quickcheck.log` in its per-job input directory, so the directory is never reclaimed | **Open** — issue drafted |

Two further findings from the gated phase, noted in their commits and drafted as issues:
`Deletion_length` is a count of deletion *events* used arithmetically as a length in bp,
and the `POS_fasta` threshold rebase in `motif_processing.py` updates a column nothing
reads afterwards.

**A pattern worth naming:** D4 and D7 are the same failure — a docstring promising an
exception type the code does not raise. Both were found by writing tests that read the
docstring as a contract. There is no reason to think they are the last two.

---

## 7. Corrections made to the plan, the issues and this branch's own commits

**This is the most useful section for a reviewer**, because it tells you which of your
own assumptions to re-check. Ten places where a document asserted something the source
contradicted were found and corrected *before* anything reached GitHub. Every one was
found by reading the source rather than by reasoning from the document.

1. **#192's "first-LEN-wins" diagnosis is wrong, and describes no input at all.** The
   issue, the #198 PR body and this branch's own spec all say the old behaviour kept the
   first `LEN`. It did not. The greedy extract took everything from the first `LEN` to the
   end and the bounded split returned the remainder, so the value **collapsed to 0
   whenever material followed the first `LEN` token**, while a compound whose only `LEN`
   is terminal parsed correctly. Where it collapsed, `frame` became 0 for a pure-insertion
   compound, and 0 is in neither `ins_frame` (3n+1) nor `del_frame` (3n+2) — so those
   states were **dropped in silence**. Measured against the real transform at `e36f3ad`
   across nine states. That is a materially worse defect than the issue describes.

2. **The disappearing-call example does not disappear.**
   `I9_2_A_LEN3&D50_2&D51_2`, carried by the spec and by this branch until `2a267aa`, is
   dropped by the deletion half but **kept** by the insertion half — frame |3−2| = 1, i.e.
   3n+1 — so the reported output is unchanged and the row merely moves halves. The state
   that genuinely stops being reported is `I9_2_A_LEN2&D50_2&D51_2`, whose net length is 0
   and which is therefore in frame and dropped by both halves. Both are now asserted.

3. **The `javap`-disproven reachability claim (D2/D3), the only error that originated with
   the orchestrator rather than with a document.** The brief asserted D2 was live and that
   the misrouting applied 3n+2 where 3n+1 belonged. Three separate parts were wrong.
   (a) `javap` on the vendored jar shows `VariantInsertion.getVcfRef()` ends in
   `Character.toString(char)` on all three return paths, so `REF="AC"`, `ALT="ACGGG"` is
   **structurally unreachable** and both defects are latent. (b) `scoring.extract_frameshifts`
   derives direction and amount from `ref_len`/`alt_len` directly and never reads the
   `Variant` label or the source file, so a misrouted row is judged **identically** —
   demonstrated empirically. The real harm is a wrong label plus a corrupted
   `Allele_Change`. (c) The worked example is Δ = +3, i.e. **in frame**, so it could not
   have illustrated a frame consequence even had the mechanism been real. The implementer
   refused to write the specification test rather than enshrine a disproven claim, and was
   right to.

4. **A fixture flag string that cannot exist.** `test_cohort_summary_oracle.py`'s fixture
   gave `sample_two` an adVNTR row of exactly `VID="Negative"`, `Flag="Coverage_flagged"` —
   a combination the pipeline **cannot produce** (both sentinel producers set
   `Flag="Not applicable"`, and `add_flags` is never reached on that path). Because of it,
   the genotype-neutral D10 fix still moved the fingerprint and flipped the adVNTR donut.
   An oracle whose fixture models an impossible state has its fingerprint hostage to rules
   about impossible inputs. The fixture was corrected to a reachable row rather than the
   hash re-derived around an impossible one; the guard stays covered by
   `test_advntr_rule_symmetry.py`, which enumerates all five reachable combinations.

5. **A selection test that would have passed before the change.** #184's drafted
   "a demotion changes which variant is reported" test gave the competitor row
   `High_Precision*` as well, so the `Depth_Score` tie-break decided the outcome either
   way and the test proved nothing. Rewritten with a `High_Precision` competitor and now
   measured in both directions.

6. **#184's blast radius was understated by 5×.** The spec said one matrix cell moves.
   **Five do** — the 54-cell matrix constructs *fractional* region depths, so the boundary
   is reachable at every alt depth there, not only at multiples of 103. On integer depths
   the narrow analysis does hold, and is now enumerated rather than assumed.

7. **The adVNTR VID defect (D10) was drafted as a live misreport; it is unreachable.**
   Investigating rather than deferring produced a better finding: VNtyper 1.3's adVNTR
   output had **no `Flag` column at all** — flagging is a 2.x invention, so there is no 1.3
   behaviour to restore — and `VID` was always 25561, hardcoded in the invocation. **It
   identifies the MUC1 VNTR locus; it was never a genotype call.** `VID == "Negative"` is a
   2.x sentinel overloading an identifier column to carry a verdict, and the rule-table
   asymmetry is a symptom of that.

8. **#188 diagnoses a three-production-file change; `--cram` already existed.** It has been
   in `cli_parser.py:92` all along, `cli_handlers.py` threads it and `pipeline.py` branches
   on it — only the task layer was wrong. But it is not "one line" either: it is **three
   lines across two files, and fixes 2 and 3 exist only because of fix 1.** Fixing the flag
   alone would have made every CRAM run without an explicit `--sample-name` be called the
   literal string `sample`, in the report and in the filenames — a regression the one-line
   fix would have introduced, confirmed by induction.

9. **Four issue drafts named mechanisms that cannot occur, and were rewritten before
   filing.** (a) A truncated BAI's partly written `chunk_end` was claimed to "read back
   larger than any real value"; `int.from_bytes(buf, "little")` on a short buffer yields
   the true value **mod 2^(8·len)**, which is always ≤ the true value at every truncation
   depth. The direction matters to the argument: too small restarts the scan early, too
   large silently skips reads. (b) The kestrel duplicate-flag draft blamed a
   "swallowed `except Exception`"; `subprocess.run(..., check=False)` never raises, and the
   real outcome is a **garbage version string** — a usage banner split on `": "` yields
   `kestrel [options]`. (c) A sentence attributed to "the module's own docstring" appeared
   nowhere in the source; it was lifted from an internal planning document. (d) The cohort
   export-leak draft quoted a CSV header with `Sample` first; the real header has it
   **twelfth** — the exact display-vs-export confusion the issue is about.

10. **Two "equivalent" fixes were not equivalent.** The D6 draft offered dropping the
    pre-built `command` as an alternative; doing so breaks `get_tool_version`'s
    `"java" in command and "kestrel" in command` dispatch, because the shipped `java_path`
    is literally `"java"`, making the reported version `Usage: kestrel [options]`.
    Demonstrated, then reverted.

Beyond those, seven drafted tests in the plan could not have failed as written, and were
corrected rather than committed: a tautological parity test that asserted properties of
its own construction; a 3-row determinism test (pandas' quicksort only reorders ties above
~260 rows on this build); a fixture missing the `Variant` column that `keep_cols`
requires, which would have died on an unrelated `KeyError`; a call to a nonexistent
`_shipped_config()` helper; a snippet using `logging.ERROR` without importing `logging`;
a mutation (`>` → `>=`) that broke none of the tests meant to catch it; and a two-argument
call to a one-argument function. One drafted fixture (`["X-5", "BROKEN"]`) did not
reproduce the defect it was written for — with that input **both rows pass**, because
suppression needs a value with two or more dashes.

The lesson this branch keeps re-learning, and the reason the section exists: **the source
is the authority, not the plan.**

---

## 8. Still open

| Issue | Trigger |
| --- | --- |
| **#189** — per-job web routes are unauthenticated; the job id is the only access-control factor | Needs a decision on whether per-job access requires a credential **distinct from the job id**. That is a product decision about how job links are shared, not a code change waiting to be written. |
| **#191** — the CRAM unmapped-read branch does not wait for its process substitution and can merge a partial BAM | Needs a **CRAM sample in the golden-cohort matrix**. #188 shows a CRAM fixture is producible (`samtools view -T … -C`), so this is now buildable; without it, any fix is unattestable. |
| **#193** — declare `bcrypt` directly in `docker/requirements-web.txt` | Unparks at the **next base-image rebuild**. `requirements-web.txt` is a base-image content-hash input (AGENTS.md trap 10), so editing it now turns the app Docker build red until a new base publishes. |

Plus the five defects drafted as new issues in §6.

**Not done, and deliberately so:** `scripts/` was to be added to a mypy target in this
PR (D13). Doing so surfaces exactly one error, in a script this branch does not touch —
`scripts/download_test_data.py:147: Need type annotation for "dir_counts"`. Fixing an
unrelated pre-existing file inside this PR would have been a silent scope widening, so the
Makefile change was reverted and the gap is documented as AGENTS.md trap 16 and filed.

---

## 9. Documentation updated

- `AGENTS.md` — the LOC/coverage table re-measured and **its argument rewritten**: the
  claim that "every module over 650 lines is under 30% covered" is now false in both
  directions (`docker/app/main.py` is 1151 lines at 88.8%, `pipeline.py` 721 at 68.5%).
  What the numbers say is that coupling to I/O predicts coverage and size does not. The
  table is now dated and marked as a snapshot, because it has been wrong twice. Also: the
  five cohort modules and `motif_decisions` added to Layout; `docker/app/` no longer
  described as ungated; the floor recorded as 80 and branch-inclusive; the mutation page
  correctly described as hand-classifying equivalent mutants **and** as machine-generated;
  new trap 16 for the `scripts/` mypy gap.
- `docs/pipeline/scoring-and-confidence.md` — the +1-frame rule (#181) and last-wins tier
  precedence (#183) recorded as specification.
- `docs/pipeline/flagging.md` — duplicate-flagging sort key (#197).
- `docs/pipeline/optional-modules.md` — SHARK `reference_assembly` (#187).
- `docs/development/mutation-testing.md` — regenerated through its generator, with the
  superseded conclusion quoted and refuted rather than deleted.

---

## Adversarial review

After the work above was complete and green, this PR was put through an adversarial
review — Codex `gpt-5.6-sol` at `xhigh`, read-only, five scoped lanes rather than one
whole-branch pass, because 18k added lines against a 272k window is a shallow read.
It returned **47 findings, 6 Critical**. Every finding spot-checked against source was
real; no false positives were found.

The recurring class is the one this branch exists to eliminate: **a document or a test
asserting something the source contradicts.** The branch removed old instances and
produced new ones. That is the honest headline.

### The merge blocker, which this branch introduced

`advntr_processing_del/ins` computed `frame = abs(Insertion_len - Deletion_length)` and
tested it against unsigned series. `abs()` discards the direction of the net change, and a
mixed state satisfies both presence guards, so either series could admit it. The
pathogenic ADTKD-MUC1 frame is signed: Δ = +1 mod 3.

Both arms are therefore correct for **pure** states — which is all the Kestrel parity test
ever covered — and wrong for mixed ones, where the sign is the whole question.

The sign error predates this branch. **`2a267aa` (#192) made it reachable**:
`I9_2_A_LEN3&D50_2` had `Insertion_len` coerced to 0 by the old parser and was dropped;
summing every `LEN` token gives 3, Δ becomes **+2** — the opposite, non-pathogenic frame —
`abs()` makes it 2, and the deletion arm reported it. A false-positive clinical call.

Two tests asserted the defect was correct, one of them named
`test_a_deletion_with_a_compensating_insertion_is_judged_on_the_net_change` while
asserting the behaviour that ignores the net change.

Fixed at **`ad515c6`**. A verdict changes iff the state is mixed and Δ % 3 == 2; nothing
is ever newly reported (new ⊆ old, `newly reported: 0` over 52,511 probes); no pure state
moves (0 of 23,064). Both are now hard failure conditions of the sweep.

**The golden cohort cannot attest this fix.** The cohort contains no mixed adVNTR state,
so a PASS is silence, not evidence. The sweep is the evidence, plus three states in
`advntr_config.json`'s `Polymorphic_Call` list that change verdict.

### What else the review found, and where it went

| Area | Finding | Outcome |
|---|---|---|
| CRAM | `tee` into `>(…)` — bash does not wait, so the merge raced the writer | `175011e` |
| Gate | `expect_exit` declared in 7 places, read in 0; two failed runs compared `IDENTICAL` | `c8baa36` |
| Gate | Nothing enforced that two *different* trees were compared | `c8baa36` |
| Gate | Attestation claimed a commit the harness never recorded | `c8baa36` |
| Cohort | The escaping test could not fail — widening the exemption passed all 2290 tests | `e369dc5` |
| Confidence | Docs said `>= 0.00515` → High; code demotes exactly 0.00515 to Low | `4b2c7d4` |
| Web | Unguarded Redis call in `finally` could skip patient-file deletion | `b91bee0` |
| I/O | `SAMPLE.CRAM` accepted by upload, rejected by the validator | `29d88bd` |
| Docs | 12 repository claims contradicted by the repository | `235b392` |

Filed rather than fixed, each needing a decision or predating the branch:
**#205** cohort ZIP identity, **#206** 20-bit pseudonym collisions, **#207** Flag-tooltip
DOM XSS, **#208** in-place mutation surviving SIGHUP, **#209** the CRAM reference
contract, **#210** the pipeline's second index, **#211** `scripts/` outside coverage.

### CRAM is now exercised for the first time (#188)

`make cram-fixtures` derives a verified CRAM beside every cohort BAM — 50/50, zero lossy,
each proved by digesting the decoded record stream from both sides. `no_ref=1`, because
the cohort's headers carry no `M5` tags and `cram_ref_option` is unconditionally empty, so
a reference-compressed CRAM would be undecodable by this pipeline. `embed_ref=2` was
measured and rejected as not lossless.

That made the CRAM race measurable rather than arguable. On `7a61_hg38_ensembl`
(622,690 unmapped pairs) the old form returned a file **28 bytes short — the BGZF EOF
block — in 3 of 3 trials**, and `samtools merge` accepted it with **exit 0**. Read loss
was intermittent: 105 reads in one trial of three; a synthetic CRAM lost 199,797 of
200,000. The truncation is deterministic, the loss is load-dependent and unbounded, and
nothing reports it. After the fix, BAM and CRAM produce identical genotypes and recover
all 622,690 reads.

These fixtures exercise the container, decoder, `.crai` indexing and unmapped scan. They
deliberately do **not** exercise reference resolution — they need no reference. That
separate failure mode is measured in #209.

### Corrections to this document

Three figures here were wrong and are now fixed: the commit count and tip, the #192
reported-call delta, and the unqualified cohort-determinism claim. They were true when
written and falsified by later commits on this same branch — which is precisely the defect
class above, occurring in the PR body describing the fix for it.
