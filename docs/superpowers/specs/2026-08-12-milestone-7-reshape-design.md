# Milestone 7 Reshape: Output Record, Faithful Report, Access Hardening

**Date:** 2026-08-12
**Supersedes:** the milestone shape set by the 2026-08-07 open-issue audit
**Milestone being dissolved:** [#7 — 6. Output model and API (parked, needs decisions)](https://github.com/hassansaei/VNtyper/milestone/7)
**Issues:** [#20](https://github.com/hassansaei/VNtyper/issues/20), [#33](https://github.com/hassansaei/VNtyper/issues/33), [#95](https://github.com/hassansaei/VNtyper/issues/95), [#119](https://github.com/hassansaei/VNtyper/issues/119), [#173](https://github.com/hassansaei/VNtyper/issues/173), [#175](https://github.com/hassansaei/VNtyper/issues/175), [#189](https://github.com/hassansaei/VNtyper/issues/189), [#202](https://github.com/hassansaei/VNtyper/issues/202), [#214](https://github.com/hassansaei/VNtyper/issues/214), [#215](https://github.com/hassansaei/VNtyper/issues/215), [#222](https://github.com/hassansaei/VNtyper/issues/222), [#223](https://github.com/hassansaei/VNtyper/issues/223), [#224](https://github.com/hassansaei/VNtyper/issues/224), [#227](https://github.com/hassansaei/VNtyper/issues/227), [#238](https://github.com/hassansaei/VNtyper/issues/238), [#242](https://github.com/hassansaei/VNtyper/issues/242)

## Decision summary

Milestone 7 held nine issues under a name that described none of the two most recent
arrivals. It is dissolved into **three milestones and one ship-now patch**, on the same
dependency-ordered principle the 2026-08-07 audit used — not on theme.

1. **A ship-now patch** carries the live regression, the last milestone-5 issue, and the two
   contained halves of the access work: #238, #214, #189 part 1, #227 part 1.
2. **M6 — Output record: one honest state per algorithm.** #175, #215, #173, #223, #224,
   #119, #202, #222. One gate: the execution-state taxonomy.
3. **M7 — Report as a faithful record.** #242 plus a new issue for a shipped documentation
   overclaim. Inherits M6's schema version.
4. **M8 — Cohort insight and access hardening.** #33, and the deferred halves of #189 and
   #227.

Two issues leave the milestone structure entirely: **#20** to backlog pending a body rewrite,
**#95** to be decided or closed.

The organising property of M6 is that **only one of its six steps is gated on a decision**.
Three land before anything needs answering.

## Standing constraint: never "clinical"

Owner decision, 2026-08-12: **VNtyper output, issues, milestones, labels and design documents
must not be described in clinical terms.** `README.md:315` and `reference/README.md:196` state
research use only. Naming a defect class "clinical safety" or a report a "clinical record"
asserts an intended use the tool has not been validated for, and a reviewer reads it
literally.

Say what is actually wrong instead: *wrong or misleading result*, *silent negative*,
*result integrity*, *faithful record*, *reported numbers*.

### Renames this requires

| Artefact | Current | Becomes |
| --- | --- | --- |
| GitHub label | `clinical-safety` — "Can produce a wrong or misleading clinical result" | `result-integrity` — "Can produce a wrong or misleading result, or a silent negative" |
| Milestone name | "Report as a clinical record" (proposed) | "Report as a faithful record" |
| `docs/superpowers/specs/2026-08-11-report-template-modernisation-design.md` title | "clinical-safety fixes and presentation modernisation" | "result-integrity fixes and presentation modernisation" |
| `docs/superpowers/plans/2026-08-11-report-template-modernisation-plan.md` | "Phase 1 — Clinical safety" | "Phase 1 — Result integrity" |
| `vntyper/scripts/screening_summary.py:15` | "The clinical text is config-driven" | "The interpretive text is config-driven" |
| `docs/pipeline/reports.md:23` | "summarizes the clinical significance" | "summarizes the combined result" |

The label carries six issues: #242, #223, #202 (open) and #212, #174, #171 (closed).
Renaming through the API preserves all six associations.

### New issue this exposes

`docs/pipeline/scoring-and-confidence.md` **claims clinical fitness that `README.md` explicitly
disclaims**, in the same shipped documentation set:

- `:74` — High_Precision is "suitable for clinical consideration"
- `:200` — a section headed "## Clinical Interpretation"
- `:202` — "The confidence level directly informs clinical decision-making"
- `:204` — "sufficient evidence for clinical reporting … per clinical guidelines"
- `:205` — "Independent validation is essential before clinical action"

This is a self-contradiction between two documents that ship together, and it is the exact
overclaim the standing constraint forbids. It is filed as a new issue in **M7**, because the
report is where a user meets that wording.

## What changed since the 2026-08-07 audit

Nine issues the milestone reasoned around have closed: **#171, #172, #174, #192, #203, #205,
#206, #207, #212**. Two consequences were never re-checked and are load-bearing here:

| Item | Audit said | True on 2026-08-12 |
| --- | --- | --- |
| #173(c) `confidence_tier` | "M, blocked on #172 (→ #171)" | Both closed. The coverage QC verdict exists and #242 already consumes it as `quality_metrics_pass`. Unblocked. |
| #33 Phase A | "Do not schedule ahead of #173, #206 and #205" | #205 and #206 closed. Phase A is buildable today against final calls. |

Two issues arrived after the audit and neither matches the milestone's name: **#242**
(report rendering, with a design and a nine-task plan already written) and **#227**
(pseudonym key management).

## Repo state at the time of this decision

Verified 2026-08-12 against `main` @ `d5fb9fb`.

| | |
| --- | --- |
| Open issues | 16 |
| Open PRs | 2, both dead |
| Milestones | #2, #3, #4, #5 closed; #6 open with 1 remaining; #7 open with 9 |

**Both open PRs are superseded and must be closed, not rebased.**

| PR | Targets | Status |
| --- | --- | --- |
| [#166](https://github.com/hassansaei/VNtyper/pull/166) `fix/detect_naming_convention` | #165 | #165 closed, fixed by PR #230. PR is `CONFLICTING`/`DIRTY`. |
| [#168](https://github.com/hassansaei/VNtyper/pull/168) `fix/report_summary` | #167 | #167 closed by commit `e331cc7`, fixed by PR #228. |

Six issues carried no milestone. Three of them are not backlog and are placed by this design:
#223 (M6), #224 (M6), #222 (M6), #119 (M6), #238 (ship-now patch), #95 (decide or close).

## Ship-now patch

Ahead of all three milestones. Each item is independently revertible.

| Item | Work | Notes |
| --- | --- | --- |
| **#238** | CRAM under Docker fails with ELOOP | **Needs triage before it can be scheduled.** Filed 2026-08-11 against 2.0.12 by a user. In the reported invocation both `-v` mounts point at the same host directory, so the run-local `.input_reference_1/reference.fa` view resolves back through its own tree and `input.cram` itself becomes unreadable. Root cause is not yet established. |
| **#214** | GHCR `latest` and semver tags; decide Docker Hub | Closes milestone #6. |
| **#189 part 1** | Gate `/download/{job_id}/` on cohort membership | See the sizing correction below. |
| **#227 part 1** | Help text and `docs/user-guide/cohort-analysis.md` state that the pseudonym is obfuscation for readability, not a privacy control | Docs only. Removes the overclaim immediately at zero risk. |
| PRs #166, #168 | Close with a pointer to #230 and #228 | Repo hygiene. |

### #189 part 1: sizing correction

The defect, sharpened from the issue's own audit comment: `/download/{job_id}/` does not
consult the cohort at all, so a job submitted into a passphrase-protected cohort is
retrievable by its bare id, forever, including after removal from that cohort. Anyone holding
a cohort passphrase can read the member ids from `/cohort-status/` and each then works
unauthenticated and indefinitely.

The deployed service is **publicly reachable** (owner, 2026-08-12), and job ids are
distributed to submitters by email in plaintext. That combination is why this rides the
patch rather than waiting for M8.

Verified at `main` @ `d5fb9fb`:

- `download_result` (`docker/app/main.py:634-664`) validates the id via `require_job_id`,
  joins it to `DEFAULT_OUTPUT_DIR`, and serves the zip. There is no cohort lookup.
- `authorized_cohort` (`:1105`) and `get_cohort_jobs` (`:1141`) exist and gate the cohort
  routes correctly.
- **`main.py:488` writes `sadd(f"{cohort_key}:jobs", job_id)` — a forward-only set. There is
  no job→cohort reverse index.**

So this is **not** a matter of calling an existing helper. The work is:

1. Write a reverse key at submission (`main.py:488`), and remove it on the enqueue-failure
   path that already calls `srem` at `:534`.
2. Have `download_result` consult it, and require the cohort passphrase when the job belongs
   to a passphrase-protected cohort.
3. **Decide the transition for jobs submitted before the change**, which have no reverse key.
   Ungated-by-default is the compatible reading and the insecure one; a Redis backfill over
   existing `{cohort_key}:jobs` sets is the complete one. This design specifies the backfill,
   because leaving pre-existing cohort jobs permanently ungated preserves exactly the defect
   being fixed.

**Named behaviour change:** a user who submitted into a passphrase-protected cohort will need
that passphrase to download. That is the intent of the fix and it will surprise someone; it
belongs in the release notes, not only in the changelog.

Non-cohort jobs are unaffected. `vntyper online` is unbroken for them. No emailed link is
invalidated except for cohort members, which is the point.

## M6 — Output record: one honest state per algorithm

**Issues:** #175, #215, #173, #223, #224, #119, #202, #222

### The problem, stated once

The pipeline has one word for four different states. A sample where Kestrel ran and found
nothing, a sample where it found only subthreshold candidates, a sample where the stage never
ran, and a sample where it exited 0 while writing an unreadable VCF all produce the same
`Negative`.

Verified mechanics at `main` @ `d5fb9fb`:

- `output_empty_result` (`vntyper/scripts/kestrel_genotyping.py:439`) writes
  `Confidence: "Negative"` with every other field the literal string `"None"`.
- `record_step` (`vntyper/scripts/summary.py:229`) carries `step`, `start`, `end`, `command`,
  `result_file`, `file_type`, `md5sum`, `parsed_result`. It sets `result_file_missing` at
  `:282` **only when the file does not exist**, so a present-but-headerless VCF sails through
  (#223).
- There is **no schema version** anywhere in `summary.py` or `pipeline.py`.
- `kestrel_pre_result.tsv` is written and never read by any non-test code, and never passed to
  `record_step`, so residual signal cannot reach `pipeline_summary.json` (per the #173 audit).

### The grain, which is the thing most likely to get this wrong

`Confidence` is a **per-variant** column, pinned by `tests/unit/test_confidence_boundaries.py`'s
full alt-depth × depth-score product and by the golden cohort. The missing concept is
**per-sample-per-algorithm**. `Confidence` is not touched by this work.

### The record

One flat value per algorithm per sample, plus a reason code:

```
state ∈ { not_run | failed | negative | negative_subthreshold | positive }
state_reason : a controlled vocabulary, not free text — absent when state is
               negative, negative_subthreshold or positive
```

The boundary between the two non-result states is **whether the stage was asked to run**.
`not_run` means the pipeline decided against it (`disabled_fast_mode`, `disabled_by_config`);
`failed` means it was asked, and did not produce a readable result
(`vcf_missing_header`, `result_file_missing`, `nonzero_exit`). The vocabulary is defined in one
module and asserted exhaustively, so a new reason cannot be introduced as an ad-hoc string.

Flat rather than two orthogonal fields (`execution_state` × `finding`), because every consumer
wants the flat value: the report banner, the cohort table, #33's `groupby` (which must exclude
everything but `positive`), and #119's linearization. Two fields would introduce an N/A-valued
column into a TSV whose current failure mode is already the literal string `"None"`.

| Sink | Change |
| --- | --- |
| `pipeline_summary.json` | `state` and `state_reason` per step in `record_step`, behind a new `schema_version` |
| `kestrel_result.tsv` placeholder row | `residual_depth_score`, `residual_alt_depth`, `residual_motif`, so cohort mode sees residual signal through `parsed_result` at no extra plumbing cost |
| Report context | `ScreeningSummary.matched_rule`, currently computed and discarded |

**One schema version serves two milestones.** #242's P5 needs the declared assembly, the
detected assembly and the resolved VNTR region persisted. That is the same bump. It happens
once, here, and M7 inherits it.

### QC stays orthogonal; no synthesised tier

#173(c) proposed a `confidence_tier` combining residual signal with the coverage QC verdict.
**It is not created.** `state`, the residual fields and `quality_metrics_pass` are three facts;
a tier is an interpretation, and synthesising one is the same overclaim the standing
constraint forbids. `report_config.json`'s forty rules already combine them for display, which
is config-driven and is the established convention.

This also answers #242's decision P6 and supplies the fallback its P1 describes.

### Sequence

| Step | Work | Gated |
| --- | --- | --- |
| **6.1** | Remove the import-time config globals: `kestrel_genotyping.py:89`, `advntr_genotyping.py:75`, `vntyper/modules/shark/shark_filtering.py:26`. Add `--threshold-profile`, thread `threads` into `run_advntr`, record the resolved profile in `pipeline_summary.json` (#175, #215) | no |
| **6.2** | #173(a): convert the six order-dependent `df.loc` writes in `confidence_assignment.py` into an explicit ordered rule table. Behaviour-preserving | no |
| **6.3** | #224 region validation and #223 VCF well-formedness check | no |
| **6.4** | **The taxonomy and `schema_version`.** Lands #173(b), #223's second half, #242's P5 and P6 | **yes** |
| **6.5** | #119: linearize the summary TSV against a now-stable record | no |
| **6.6** | #202 and #222, in parallel — see below | see below |

`run_advntr` (`advntr_genotyping.py:78`) takes no `threads` parameter and reads
`advntr_settings.get("threads", 8)` at `:95`, emitted at `:136`. **`advntr_config.json` sets
`1` while the code default says `8`** — the two disagree and no record says which is
intended. 6.1 must resolve that rather than preserve it silently.

### 6.6 — the two domain questions, reframed

Both were previously parked as "needs a decision from @hassansaei". Neither is blocked on him.

- **#202** — "can a single adVNTR `D` token denote more than one deleted base?" The issue
  itself offers the alternative: *"or from adVNTR's State grammar"*. adVNTR is open source.
  This is a **research task**, answerable by reading their source. If the answer is yes, the
  fix is genotype-affecting and needs a golden-cohort gate run; if no, it is a rename plus an
  assumption test. Note the golden cohort **cannot** answer it — its only compound state,
  `D17_2&D18_2&D19_2&D20_2&D21_2`, is five single-base events and is identical under either
  reading.
- **#222** — which interval is the homologous VNTR span across all six region keys, and what
  flank convention applies. That is a MUC1-VNTR domain call and it belongs to the project's
  own domain expertise, not to a routing decision. It moves every hg38 coverage number and can
  move `advntr_result` through `downsample_bam_if_needed`, so it needs a golden-cohort gate
  run and a decision on whether historical numbers are restated.

## M7 — Report as a faithful record

**Issues:** #242, plus the new `scoring-and-confidence.md` overclaim issue.

#242 already carries a design (`docs/superpowers/specs/2026-08-11-report-template-modernisation-design.md`)
and a nine-task plan. This design changes three things about it:

1. **Its P6 is answered** by M6.4's taxonomy, and its **P5 rides M6.4's `schema_version`**.
   M7 therefore depends on M6.4 landing, and on nothing else in M6.
2. **Its P1** — whether a top-level verdict word is created — resolves consistently with the
   standing constraint: **no verdict word.** The banner renders the configured text plus the
   three orthogonal state chips, exactly as the issue's own fallback describes.
3. **Its wording is renamed** per the standing constraint, in the design, the plan and the
   commit messages.

P2 (artifact shape), P3 (sample identity) and P4 (cohort scope) remain open and gate tasks 9,
3 and 7 respectively. Phase 1 tasks 1, 2, 4, 5 and 6 are ungated and can start as soon as
M6.4 lands.

## M8 — Cohort insight and access hardening

**Issues:** #33, #189 (deferred half), #227 (deferred half)

| Item | Work |
| --- | --- |
| **#189 part 2** | Link expiry. Addresses "the capability is mailed in plaintext and never expires" without a full credential model. Needs a TTL store and **does** invalidate outstanding links, which is why it is not in the patch. |
| **#227 part 2** | HMAC-SHA-256 pseudonyms keyed from `VNTYPER_PSEUDONYM_KEY`, absent by default, current unkeyed behaviour as the documented fallback. Keeps the secret out of a config file that gets copied around, and out of `--config-path`'s whole-file replacement trap. |
| **#33** | Rewrite the body to Phase A / Phase B first — the proposal as written is unbuildable, because cohort mode loads only `pipeline_summary.json` and never sees `kestrel_pre_result.tsv`. Phase A is a cohort **call-frequency** table over final calls, named honestly. Phase B falls out of M6.4 at no extra cost. |

**#33 Phase A has a hard requirement:** every cell is sample-derived and #207 was a DOM XSS in
this exact report. Escaping is not optional and must be asserted, not assumed.

## Out of scope

- **#20** returns to backlog until its body is rewritten. It is one sentence; the real content
  is its audit comment, which establishes that roughly 70% of the proposed checkpoint system
  already exists as `pipeline_summary.json` and that what is missing is a *reader*. A second
  manifest must not be added. It also depends on #223 landing first: a `--resume` that trusts
  recorded step completion is only safe once completion is recorded honestly.
- **#95** (`mutation_counter`) is decided or closed. A `needs-decision` label open since
  2025-02-04 is itself a decision.
- Per-job bearer tokens for the web API. The narrow cohort fix plus expiry addresses the
  measured defect; tokens break `vntyper online` and invalidate every outstanding link for a
  122-bit CSPRNG identifier that is not brute-forceable.
- Tuned per-assay threshold values. #175 delivers a reachable and traceable operating point,
  not validated per-kit cut-offs, which need data across kits that does not exist here.

## Verification

| Claim | How it was established |
| --- | --- |
| Both open PRs are superseded | `gh api` issue timelines: #165 cross-referenced by PR #230, #167 closed by commit `e331cc7` via PR #228 |
| `/download/` performs no cohort lookup | Read `docker/app/main.py:634-664` |
| No job→cohort reverse index exists | `main.py:488` is `sadd(f"{cohort_key}:jobs", job_id)`; `cohorts.py` exposes `cohort_job_ids` only |
| No schema version in the summary | `grep -rn "schema_version" vntyper/scripts/summary.py vntyper/scripts/pipeline.py` → no match |
| `result_file_missing` only fires on absence | `summary.py:280-282` |
| Three import-time config globals | `kestrel_genotyping.py:89`, `advntr_genotyping.py:75`, `vntyper/modules/shark/shark_filtering.py:26` |
| adVNTR threads config/code disagreement | `advntr_config.json` `"threads": 1` vs `advntr_genotyping.py:95` `.get("threads", 8)` |
| Documentation overclaim | `docs/pipeline/scoring-and-confidence.md:74,200,202,204,205` vs `README.md:315`, `reference/README.md:196` |

Anything genotype-affecting — #202 if the answer is yes, #222 in every case, and #173(b) if the
residual fields change which rows are written — requires a golden-cohort gate re-run. The gate
is hand-run and not in CI.

## Traps carried forward

- `tests/unit/test_kestrel_filtering.py:72` reads `kestrel_genotyping.py` **as source text** and
  asserts exactly five gate entries. 6.4 changes it; land it once.
- `tests/unit/test_cohort_summary_oracle.py:411` holds a whole-document fingerprint. Any cohort
  output change must deliberately re-baseline it.
- `pyproject.toml:135` declares package data for `vntyper.scripts` only, and `MANIFEST.in:9-10`
  names exactly two templates. New template files would be absent from wheels while every
  editable-checkout test passed.
- `--config-path` replaces the whole config rather than merging it (AGENTS.md trap 2). Any new
  configurable value inherits that.
