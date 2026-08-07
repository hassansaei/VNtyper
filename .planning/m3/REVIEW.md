# Milestone 3 — adversarial review of SPEC.md and PLAN.md

**Reviewer:** `codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox read-only`,
302,463 tokens, against `main` @ `b46da80` and the pinned dependency sources.
**Received:** 2026-08-07. 15 findings, 5 HIGH.

Every finding was verified against the code before being accepted or rejected. Verdicts below;
the evidence command is given where the verification is not obvious from the citation.

## Findings I found myself, before the review landed

Recorded separately so the review's independence is visible.

| # | Sev | Finding | Action |
| --- | --- | --- | --- |
| S1 | HIGH | `tests/unit/test_golden_cohort_compare.py:219` asserts the literal source string `return sorted(processed_dirs), temp_dirs`, and `:230` asserts `COHORT_ORDER_WHY` still says the ZIP temp path is in the sort key. Task C3 falsifies both. | Task **C4** added; C's ownership widened to `scripts/golden_cohort/compare.py` and `tests/unit/test_golden_cohort_compare.py`. Codex found the same thing independently (#8). |
| S2 | HIGH | `tests/unit/test_cohort_summary_oracle.py:862-863` pins the literal MD5 pseudonym `anon_65622`. Task C1 changes it. Both A and C reach this file. | Assigned to **C**; A is forbidden from editing it and must stop and report if `EXPECTED_FINGERPRINT` moves. |
| S3 | MEDIUM | `tests/unit/test_cohort_inputs.py:182-186` `_ORDER_PROBE` does `[p.name for p in ...]` — an `AttributeError` against `DiscoveredSample`. | Folded into C4. |
| S4 | MEDIUM | SPEC claimed igv-reports emits a trailing `;`. It does not — `variant_template.html:155-156` has no terminator and `report.py:178-183` substitutes the placeholder including its quotes. | SPEC corrected. Codex found the same (#12). |
| S5 | — | Verified `textContent` is not a display regression: VNtyper's BED is three columns (`kestrel_genotyping.py:682-686`), so `has_name` is false and the one column igv-reports pre-escapes (`bedtable.py:32`) is never rendered. | SPEC corrected with the evidence. |

## Verdicts

### Accepted in full

**1. HIGH — the report fix does not establish "same IGV panel." CONFIRMED, and it changes D's design.**

Verified: `generate_report.py:492` is `if bed_file and os.path.exists(bed_file):`, and
`cli_report.py:104-115` discovers `--bam-file`/`--bed-file` **only** when `--input-dir` is given.
So `vntyper report -o <run>` resolves no BED, skips IGV entirely, and the newly-resolved VCF
changes nothing for that invocation.

My spec explicitly declared widening bam/bed "out of scope". That was wrong: without it, #167's
fix is inert for exactly the invocation `--output-dir`'s own help text describes ("Directory
containing pipeline results"). It is *not* wholly inert — for `--input-dir` users the VCF fix is
real — but shipping a fix that does nothing for the documented common case is not shipping a fix.

**Action:** D resolves BAM, BED and VCF from one effective run directory: explicit flag, then
`--input-dir`, then `--output-dir`. The test asserts `run_igv_report` is actually invoked with
all three, not merely that the spy received `vcf_file`.

**4. HIGH — the cohort identity map is not injective over actual samples. CONFIRMED.**

Two distinct discovered directories can share a basename (`a/sample/` and `b/sample/`), so their
identities are equal *before* hashing. My guard compares `existing != original_sample`, which is
false for that pair, so it passes — and `cohort_categories.sample_categories` then groups the two
patients as one. No digest width fixes identical inputs.

**Action:** Task C2 detects repeated **identities** across distinct `DiscoveredSample.directory`
values, before pseudonymisation, and aborts naming both source directories. This is a
pre-existing defect wider than #206; it is in scope because the injectivity invariant is the one
this milestone claims.

**6. MEDIUM — the B-series verification is broken or non-discriminating. CONFIRMED.**

The B3 assertion `all(str(tmp_path / "input") not in c for c in index_commands)` is simply
broken: a correct `samtools index -o <out> <in>` command must name the input as its operand, so
the test fails *after* the fix. B1 is a characterisation test that passes before the fix (the
plan said so, but did not say it proves no pipeline behaviour). B4 mocks the pipeline subprocess,
so no sidecar is ever created and it cannot distinguish fixed from unfixed.

**Action:** B3 asserts on `shlex.split(command)` — specifically that the argument after `-o` is
under the output directory. B4 gains a subprocess simulator that reproduces the current sidecar,
so it fails on `main`. The pipeline harness gains explicit `log_dir` assertions for BAM and CRAM.

**7. MEDIUM — B misses existing tests and documentation encoding the old behaviour. CONFIRMED, all three.**

- `alignment_processing.py:128` is an unlisted production caller of
  `build_samtools_index_command`. The new parameter is keyword-only and optional, so it is
  source-compatible — but an unlisted caller of a changed signature is exactly the omission this
  review was for.
- `tests/unit/web/test_index_handoff.py:124` (`_index_the_way_the_pipeline_does`) and `:406`
  (`test_the_index_the_pipeline_builds_for_itself_is_removed_too`) *simulate* the pipeline
  writing `f"{in_bam}.bai"`. After B3 the pipeline no longer does that, so the tests keep passing
  while documenting behaviour that no longer exists.
- `tests/docker/conftest.py:234-236` copies the BAM out of the read-only mount with the comment
  *"VNtyper writes log files next to input BAM, but input is read-only."* That workaround is the
  reason the Docker tier cannot currently prove E1 — it masks the very behaviour under test.

**Action:** all three added to B's ownership. The `conftest.py` workaround is removed so the
Docker tier runs against the read-only mount directly, which is also E1's proof mechanism.

**9. MEDIUM — several proposed cohort tests are empty or probabilistic. CONFIRMED.**

Five tests in C2/C3 were given as docstring-only shapes, which pass without asserting anything —
a direct violation of the writing-plans no-placeholder rule I acknowledged and then broke. The
two-ZIP ordering test could also pass on `main` by chance.

**Action:** every test body written out. `tempfile.mkdtemp` is monkeypatched to produce
deliberately reverse-sorting roots, so the ordering test fails deterministically on `main`.

**10. MEDIUM — E2 claims more than its proof establishes. CONFIRMED.**

The cohort HTML carries a report timestamp and Plotly UUIDs (which is why
`test_cohort_summary_oracle.py` has a `_skeleton()` normaliser at all), so "two runs are
byte-comparable" cannot be true of the whole output.

**Action:** E2 restated as **"the machine-readable cohort exports (CSV/TSV/JSON) and the
pseudonymization table are byte-identical across two runs; the HTML is identical after the
oracle's existing normalisation."**

**11. MEDIUM — the A-series threat model names the wrong provenance. CONFIRMED.**

`flagging.py:135-147`: `df_copy["Flag"] = ", ".join(flag_list)` where `flag_list` holds the
**keys** of `flagging_rules` — config-declared identifiers like `False_Positive_4bp_Insertion`.
No DataFrame value is interpolated. My spec's "composes it from DataFrame values derived from the
sample's own reads" is wrong.

Issue #207 itself says the correct thing — *"`Flag` values come from supplied
`pipeline_summary.json` files, so they are not trusted markup"* — and I embellished past it.

**Action:** SPEC restated to the issue's own trust boundary: report input and cohort ZIP contents
are untrusted. Tests use a malicious stored artefact, not a nonexistent rule-interpolation path.

**12. MEDIUM — `ensure_ascii` already escapes the JS line separators. CONFIRMED.**

`json.dumps({"a": "U+2028"})` returns `{"a": "U+2028"}` escaped, because `ensure_ascii` defaults
to True. My two `.replace()` calls are dead code and my test would have passed without them.

**Action:** `ensure_ascii=True` written explicitly with the reason; the two replaces dropped; only
`</` → `<\/` remains. The line-separator test is kept as a property assertion on the output (the
separator must not appear literally), which is still worth pinning.

**13. LOW — the collision probability was wrong by three orders of magnitude. CONFIRMED.**

Computed: `1 - exp(-n(n-1)/2**49)` gives **1.78e-9 at 1,000 samples** and **1.78e-7 at 10,000** —
not the 1.8e-10 I quoted for 10,000.

**Action:** numbers corrected everywhere. **The decision does not change**: 12 hex characters was
chosen against a milestone whose realistic cohorts are hundreds of samples (1.78e-9 at 1,000),
and Task C2 now aborts on a collision rather than merging, so the residual is a failed run rather
than a silent patient merge. Also added: the configured length is validated as a positive integer
within the digest's range, which the plan did not check.

**14. LOW — the ZIP-order mechanism is misidentified. CONFIRMED.**

For ZIP inputs the irreproducibility comes from `mkdtemp`'s random suffix being in the sort key,
not from `PYTHONHASHSEED`. The hash-seed mechanism is real but applies to the *set iteration* the
determinism fix already closed.

**Action:** the ZIP ordering test monkeypatches `mkdtemp` with controlled reverse-sorting roots.
The hash-seed cross-process test stays, testing what it actually tests.

**8. MEDIUM — changing `discover_sample_directories` breaks unlisted contracts.**

Independently found (S1-S3) and already applied before the review landed. Codex adds one item I
had not: many `test_cohort_inputs.py` assertions treat returned items as `Path`
(`dirs == [Path(temp_dirs[0])]` at line 271, `[d.name for d in dirs]` at 292). C owns that file;
the task now names those call sites.

### Accepted in part, with reasoned pushback

**2. HIGH — E1 is impossible as stated. ACCEPTED as a wording defect; the suggested fix is rejected.**

Verified: `docker/app/tasks.py:351` writes `status="completed"` *before* the `finally:` cleanup
block, and cleanup failures are logged and swallowed — pinned deliberately by
`test_tasks.py:715` and `:752`, which require the task not to raise when removal fails. So a
"completed" job can leave a file behind, on `main` and after this plan.

**Rejected:** making verified cleanup a condition of completion, or adding a
`completed_with_cleanup_error` state. That is a change to the web service's error model that no
issue in this milestone asks for, and it would invert two existing tests whose rationale — a job
whose pipeline completed and whose results are written must not be reported as a failure — is
recorded in the code.

**Accepted:** the criterion was overclaimed. E1 is restated as:

> **E1.** On a web job whose cleanup succeeds, no file and no directory remains under the input
> tree — and cleanup no longer fails on the normal path, so the `os.rmdir` that was dead code in
> production now runs. A cleanup that fails for an unrelated filesystem reason is still logged
> and swallowed; that behaviour is deliberate and out of scope.

**3. HIGH — E3 has no proof capable of detecting the DOM exploit. ACCEPTED as an overclaim; the suggested fix is rejected.**

The substance is right: the exploit is a runtime `.text()` → `.html()` transition, and no grep of
rendered HTML observes it. My regex also only matched a single-quoted `.html('literal' + var)`
spelling.

**Rejected:** adding a jsdom/Playwright behavioural test. `AGENTS.md` is explicit that unit tests
must be pure — `tmp_path` + `unittest.mock`, no network, no Docker — and CI runs `pytest -m unit`
only. A browser tier does not exist in this repo, and introducing one is a build-system change
that dwarfs the milestone. That is a real coverage limit and it is stated rather than papered
over.

**Accepted:** the regex is broadened to catch double quotes, template literals,
`insertAdjacentHTML`, `outerHTML` and `innerHTML` assignment from a variable, and E3 is restated
to what the evidence supports:

> **E3.** No template builds markup from a value it did not author: every sample-derived value in
> both templates reaches the DOM through `.attr()`, `.text()` or `textContent`, and every value
> reaching a `<script>` block is the output of `json.dumps`. Proven by a source-text tripwire over
> both templates plus a Python-side test that a malicious `Flag` in a stored
> `pipeline_summary.json` is server-escaped in the rendered HTML. **Not** proven by execution —
> this repo has no browser test tier, and the residual risk is a future sink the tripwire's
> pattern does not match.

**5. HIGH — the "no writes into the input directory" invariant still does not hold. ACCEPTED as an overclaim; the suggested worker rewrite is rejected.**

Verified: `docker/app/tasks.py:290` runs `subprocess.run(["samtools", "index", bam_path])`, which
writes `<upload>.bam.bai` into the job's input directory before the pipeline starts.

**Rejected:** routing the worker's index to a scratch directory. The worker puts it beside the
alignment precisely so the pipeline finds it, #199 already made the worker remove every
deterministically derived index name it created, and rethreading that handoff is milestone-4
work at best.

**Accepted:** the invariant was stated too widely. It is narrowed to what #162 actually reports —
*"the input directory does not regularly allow writing to it to protect the data"*, which is the
CLI-with-a-mounted-dataset case, not the web service's own writable upload volume:

> **Invariant.** No code in the `vntyper` package writes into the input directory. The web
> worker deliberately creates one index beside the upload and deliberately removes it in cleanup
> (#199); the upload volume is writable by design, so that is not a violation. The read-only-input
> guarantee is a guarantee about `vntyper pipeline`.

**15. HIGH — concurrent workstreams on one branch are operationally unsafe. ACCEPTED; the suggested fix is adapted.**

The reasoning is correct and I had missed it entirely: disjoint *files* do not make concurrent
`git add`/`git commit` safe, because `.git/index`, `index.lock`, `HEAD` and the branch ref are
shared. Two subagents committing at the same moment can collide on the lock or stage each other's
work into the wrong commit.

**Rejected:** a separate branch per workstream. The user's constraint is explicit — one branch,
one PR, never stacked — and it is the right constraint here, because CI only fires on PRs into
`main`.

**Accepted, adapted:** **subagents run no git commands at all.** They edit files and run pytest;
the integrator stages and commits each workstream sequentially, in dependency order, after its
tests pass. This removes the race without touching the branching model, and it also means the
per-task `git commit` blocks in the plan become the integrator's checklist rather than the
subagents'.

### Corrections to the reviewer

- Finding 1 cites `vntyper/cli_report.py`; the file is `vntyper/scripts/cli_report.py`. The
  finding itself is correct.
- Finding 15's claim that `tests/support/pipeline_harness.py`, `tests/builders.py` and
  `tests/conftest.py` need no edits agrees with the plan, which already declared the first two
  read-only. Finding 6 then asks for harness assertions on `log_dir` — those go in
  `tests/unit/test_pipeline_cwd.py`, which already reads the harness, rather than in the harness
  itself.

## Result

No HIGH finding remains unaddressed: 1, 4 accepted in full and redesigned; 2, 3, 5, 15 accepted
as overclaims with the criteria restated and the suggested remedies rejected on stated grounds.
All MEDIUM and LOW findings accepted.

One decision the user made in Phase 1 was based on a number that turned out to be wrong (#13).
The corrected figures still support the choice by a wide margin and Task C2 now makes a collision
loud rather than silent, so the decision stands — recorded here rather than re-asked.
