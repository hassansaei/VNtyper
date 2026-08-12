# Milestone 7 Reshape: Output Record, Faithful Report, Access Hardening

**Date:** 2026-08-12
**Revision:** 2 — rebuilt after nine adversarial reviews (eight parallel agents plus Codex
gpt-5.6-sol at high effort). Revision 1's taxonomy was **not implementable**; see
[What revision 1 got wrong](#what-revision-1-got-wrong).
**Supersedes:** the milestone shape set by the 2026-08-07 open-issue audit
**Milestone being dissolved:** [#7 — 6. Output model and API (parked, needs decisions)](https://github.com/hassansaei/VNtyper/milestone/7)

## Decision summary

Milestone 7 held nine issues under a name that described neither of its two most recent
arrivals. It is dissolved into **three milestones and one ship-now release**, on the same
dependency-ordered principle the 2026-08-07 audit used.

Two of its issues are already closed by this work: **#215** (adVNTR's `-t` is a proven no-op,
so the requested fix would be a false affordance — replaced by #247) and **#202** (the adVNTR
`D`-token question is answered from source; the code was correct).

1. **A ship-now release, 2.0.14**, carrying the #238 CRAM regression. Users are affected today.
2. **M6 — Output record: one honest state per algorithm.** #173, #175, #223, #224, #119, #222,
   #247. One gate: the execution-state design.
3. **M7 — Report as a faithful record.** #242 plus a new documentation-overclaim issue.
4. **M8 — Cohort insight and access hardening.** #33, #189, #227.

**#20** goes to backlog pending a body rewrite, and its audit comment must be replaced first —
the defect that motivated it was fixed five days after it was written. **#95** stays open by
owner decision (2026-08-12); the evidence for deletion is recorded below but not acted on.

## Standing constraint: never "clinical"

Owner decision, 2026-08-12: **VNtyper output, issues, milestones, labels and design documents
must not be described in clinical terms.** `README.md:315` and `reference/README.md:196` state
research use only. Naming a defect class "clinical safety" or a report a "clinical record"
asserts an intended use the tool has not been validated for.

Say what is actually wrong instead: *wrong or misleading result*, *silent negative*,
*result integrity*, *faithful record*, *reported numbers*.

### Renames this requires

| Artefact | Current | Becomes |
| --- | --- | --- |
| GitHub label | `clinical-safety` — "Can produce a wrong or misleading clinical result" | `result-integrity` — "Can produce a wrong or misleading result, or a silent negative" |
| `docs/superpowers/specs/2026-08-11-report-template-modernisation-design.md` title | "clinical-safety fixes and presentation modernisation" | "result-integrity fixes and presentation modernisation" |
| `docs/superpowers/plans/2026-08-11-report-template-modernisation-plan.md` | "Phase 1 — Clinical safety" | "Phase 1 — Result integrity" |
| `vntyper/scripts/screening_summary.py:15` | "The clinical text is config-driven" | "The interpretive text is config-driven" |
| `docs/pipeline/reports.md:23` | "summarizes the clinical significance" | "summarizes the combined result" |

The label carries six issues: #242, #223 (open) and #212, #202, #174, #171 (closed). Renaming
through the API preserves all six associations.

### New issue this exposes

`docs/pipeline/scoring-and-confidence.md` **claims clinical fitness that `README.md` explicitly
disclaims**, in the same shipped documentation set: `:74` "suitable for clinical consideration";
`:200` a section headed "## Clinical Interpretation"; `:202` "directly informs clinical
decision-making"; `:204` "sufficient evidence for clinical reporting … per clinical guidelines";
`:205` "essential before clinical action". Filed as a new issue in **M7**.

## What revision 1 got wrong

Recorded because the errors are instructive, and because two of them were repeated from the
issues themselves.

| Revision 1 claim | Reality |
| --- | --- |
| Store `state` per step in `record_step` | **Structurally impossible.** `pipeline_kestrel.py:60` runs the runner *before* `record_step` at `:76`; `:86-88` deletes partial records; `pipeline.py:664-668` exits without a summary write. A failed algorithm leaves **no step record at all**. Independently found by two reviewers. |
| Cohort mode will see the state | **No.** `cohort_inputs.py:572-586` reads only `step["parsed_result"]["data"]` and discards the step's own keys. |
| A new flat 5-value enum | **Duplicates and degrades** `ScreeningSummary.kestrel_result`/`.advntr_result`, which is already at this grain, already config-driven, already has `NOT_PERFORMED`, and distinguishes flagged/unflagged and precision tier — which the flat enum collapsed into one `positive`. |
| `test_kestrel_filtering.py:72` asserts five gate entries | `:73`, and it asserts **six**. `flag_filter_pass` was added by #174. Following revision 1 literally would have changed the gate list — a genotype change. |
| `pyproject.toml:135` declares package data for `vntyper.scripts` only | False. A second `"vntyper"` block at `:142-149` already ships both templates. The conclusion (explicit file lists, no glob, so a new partial still would not ship) survives; the premise does not. |
| #242's P5 needs the same schema bump | **Two-thirds already satisfied.** `summary.py:80-89` persists `reference_assembly_requested`, `reference_key_used`, `reference_path`, `reference_source_effective`, and `cli_handlers.py:323-324` falls back to config only when the flag is absent — so it **does** survive `--reference-assembly`. Only the resolved VNTR region and `--custom-regions` are unpersisted. |
| #189: retrievable "forever, including after removal from the cohort" | Both false. There is **no removal mechanism at all** — membership is add-only, with two writers (`main.py:488`, `tasks.py:373`). And archives are reclaimed at `MAX_RESULT_AGE_DAYS`, **7 days** by default (`tasks.py:466-516`, `docker-compose.yml:53,83,107`). |
| Three import-time config globals | **Two.** SHARK already receives `threads` (`pipeline.py:347`) and resolves its only user-relevant setting from the argument; `test_shark_filtering.py:134` asserts *"patching the module global has no effect"*. Its globals are dead code to delete, not plumbing to add. |
| `disabled_fast_mode` as a `not_run` reason | `--fast-mode` skips unmapped-read extraction (`fastq_bam_processing.py:157-158`). It disables neither algorithm. `not_run` is reachable for adVNTR only. |

## Repo state

Verified 2026-08-12 against `main` @ `d5fb9fb`, with the reviews run at `ebcda33`.

Open PRs **#166** and **#168** were closed as superseded (by `6e344b6` and `f9fc27c`/`e331cc7`).
Neither loses anything; both in fact carried a latent defect, recorded in their closing
comments.

## Ship-now: release 2.0.14

**#238 — CRAM under Docker fails with ELOOP. Root-caused; the fix is written but unreleased.**

`reference_binding.py:293` builds the run-local view directory. For a **generated** sidecar the
pre-fix code replaced the destination with a symlink to `/proc/<pid>/fd/<n>` where the
descriptor had been opened on *that same pathname*, so the `os.replace` unlinked the only name
and left the descriptor's recorded path pointing back at itself. Any layer resolving the procfs
link by re-walking its text closes the cycle — ELOOP after 40 hops.

- **Trigger:** the reference FASTA has no `.fai` beside it, so htslib generates one during the
  probe. A shipped index is bound with `generated=False` and is never self-referential — which
  is why candidate 3 succeeded in the reporter's log while candidate 1 failed.
- **The double bind-mount is not the cause.** See the separate defect below.
- **`input.cram` is not actually unreadable** — `--reference` is a format-specific open option,
  so `refs_load_fai` failing makes `cram_open` return NULL with `errno` still `ELOOP`, and three
  layers print it. The reporter's own log shows the CRAM decoding one second later.
- **Affected: v2.0.10, 2.0.11, 2.0.12, 2.0.13** — introduced by `1c861aa`. 2.0.13 is GHCR
  `:latest`, so every Docker user with an unindexed reference is affected now.
- **Workaround, posted on the issue:** `samtools faidx <reference>.fa` before running.
- Not genotype-affecting: pre-fix the run aborts with no output; both link strategies resolve to
  the same inode.

**Also unfiled and needing an issue:** `alignment_target_io.py:490-498` compares paths lexically
and via `resolve()` only, never `os.path.samefile` (which exists at `:29-33` and is used for
protected alignment paths at `:503`). With `-v .:/input -v .:/output` the two are the same host
directory but lexically distinct, so **the "output root must stay outside the patient input
tree" guard silently does not fire**. This is the genuine input-equals-output defect. The #238
plan states it was filed separately; it was not.

## M6 — Output record: one honest state per algorithm

**Issues:** #173, #175, #223, #224, #119, #222, #247

### The problem

The pipeline has one word for four different situations: an algorithm that ran and found
nothing, one that found only subthreshold candidates, one that never ran, and one that exited 0
while writing an unreadable VCF.

Verified mechanics at `main`:

- `output_empty_result` (`kestrel_genotyping.py:439`) writes `Confidence: "Negative"` with every
  other field the **literal string `"None"`** — not `NaN`, which matters for any consumer.
- The two negative paths produce **byte-identical** files (demonstrated: same md5), and
  identical `parsed_result["data"]` unconditionally.
- `record_step` sets `result_file_missing` only when the file is absent (`summary.py:280-282`),
  so a present-but-headerless VCF sails through.
- There is **no schema version** anywhere in `summary.py` or `pipeline.py`.

### Design: extend `ScreeningSummary`, and make execution state pre-emptive

Owner decision, 2026-08-12: **extend the existing config-driven vocabulary rather than
introducing a parallel one.** `ScreeningSummary.kestrel_result`/`.advntr_result`
(`screening_summary.py:78-99`) is already the per-sample-per-algorithm grain.

The combinatorial objection to extending it is real but avoidable. `report_config.json` is
today an exact **5 × 4 × 2 = 40-rule cartesian product with zero gaps** (verified
programmatically). Adding execution states as a third finding axis would cost 84 rules and 44
new sentences.

**It is not a third axis.** A failed Kestrel does not need a different sentence for each of
adVNTR's four findings — it needs one sentence saying the run is incomplete and which stage
produced no result. So:

- **`execution_state` is evaluated first and pre-empts the rule table.** If either algorithm is
  not `completed`, the report renders an incompleteness statement naming the stage and the
  reason. The 40 rules are untouched and continue to describe the both-completed case only.
- **Only `negative_subthreshold` is a genuine new finding value**, so `kestrel_result` goes from
  5 to 6 values: **8 new rules**, not 44. That is the honest cost.

```
execution_state ∈ { completed | not_run | failed }      # pre-emptive
state_reason    : controlled vocabulary, required when execution_state ≠ completed
```

`not_run` means the pipeline decided against the stage (adVNTR only — `--fast-mode` does not
disable either algorithm). `failed` means it was asked and produced no readable result
(`vcf_missing_header`, `result_file_missing`, `nonzero_exit`).

**`state_reason` is also populated for negatives**, sourced from the failing gate name. This is
load-bearing, not decoration: there are **six** gates, and a row killed by `flag_filter_pass`
(#174) is a known artifact whose `Depth_Score` may be excellent. Calling that "subthreshold"
says *weak signal* when the truth is *strong signal, deliberately discarded*. The vocabulary is
free — it is `filter_cols` (`kestrel_genotyping.py:814-821`), already computed and already
pinned.

The gate's own admissibility rules make this decisive: `admissibility.py:67-73` requires
`kestrel_pre_result.tsv` on **every** zero-exit case, so **all ten golden-cohort negatives had
gate-rejected candidates**. The truly-empty negative barely occurs; `negative_subthreshold`
would absorb essentially every real negative, leaving the reason as the only informative field.

### Where it lives

**Not in `steps[]`.** A top-level `algorithms` block in `pipeline_summary.json`, written on both
the success path and the exception path, because that is the only place a `failed` or `not_run`
state can exist at all.

Two consumers must be taught to read it, and neither does today:

- `screening_summary.py:276-290` — `build_screening_summary` must receive the execution states.
- `cohort_inputs.py:572-586` — `parse_pipeline_summary` reads only `parsed_result.data`.

`run_advntr` returns non-zero on failure (`advntr_genotyping.py:110-158`) and its caller
**ignores the return value** (`pipeline.py:545-568`). `nonzero_exit` cannot be derived until
that is consumed, and whether an algorithm failure aborts the pipeline or permits a report is
part of this decision.

### The `negative_subthreshold` classifier must be specified, not implied

Three empty paths reach `output_empty_result` and two of them run **before anything is scored**
(`kestrel_genotyping.py:385`, `:402`) — no residual frame exists there. Only `:413` has one, and
`filter_final_dataframe:801` writes it to disk and returns only the filtered frame.

So the design must state: which gates establish *eligibility* (a non-frameshift is not a
subthreshold candidate), the ranking key and tie-break for selecting the residual row, and
whether artifact and motif failures are excluded. Without that, "subthreshold" can report
evidence from a row the pipeline explicitly says is not a candidate.

### Schema version comes first, not at 6.4

Revision 1 put `schema_version` at 6.4 while 6.1 already added a field to
`pipeline_summary.json`. Introduce the version **before the first structural change**.

Note also that `convert_summary_to_*` flatten `summary["steps"]` only, so **top-level additions
are silently absent from the linearized output** — which is #119's second defect, below.

### Sequence

| Step | Work | Gated |
| --- | --- | --- |
| **6.0** | Introduce `schema_version` and version dispatch with a legacy fallback | no |
| **6.1** | Remove the import-time config globals in `kestrel_genotyping.py:89` and `advntr_genotyping.py:74-75`. Add `--threshold-profile`; record the profile in the summary. Delete SHARK's dead globals and `test_shark_filtering.py:134` | no |
| **6.2** | #173(a): the six order-dependent `df.loc` writes → an explicit **last-wins** rule table | no |
| **6.3** | #223 VCF well-formedness (**in `run_kestrel`**) and #224 region validation | no |
| **6.4** | #119: linearize the summary TSV, and carry the top-level fields into it | no |
| **6.5** | **The execution-state design and `negative_subthreshold`** | **yes** |
| **6.6** | #222 region harmonisation, and #247 | see below |

`--threshold-profile` is not yet an implementable specification. It must decide: name or path,
the profile schema, precedence over the packaged and main configuration, validation behaviour,
and what provenance is recorded (profile identifier plus content digest and effective
thresholds). These choices are genotype-affecting. Note `--config-path` **replaces** the whole
config (AGENTS.md trap 2) and `confidence_assignment.py:104-111` deliberately raises `KeyError`
on a partial config — so a named profile is safe and an arbitrary path is not, unless merged.

**#119 moved from last to 6.4.** It is provably independent of the taxonomy: `flatten_dict` is
generic and the new fields are flat scalars, so they flow through with no code change. Leaving
it behind the gated step made an ungated, zero-blast-radius item inherit schedule risk. It has
no in-repo consumers — the golden gate writes the TSV and never reads it, and the web service
touches neither — so linearizing breaks no code, only an unversioned user spreadsheet. Its
second defect is worse than the one filed: `input_files`, `version` and all four
reference-provenance fields are **absent from every row**, because the converters flatten
`steps` only. A table of record that cannot state which reference the run resolved is a worse
defect than JSON-in-a-cell.

### 6.2 — two provably dead branches

`cond3` is dead: deleting `confidence_assignment.py:195` changed **0 labels across 780,000
integer cells** and ~310,000 fractional ones. `cond3 ⊂ cond_midband`, both write the same label,
and `cond_midband` is written after. `cond1` is dead on every **integer** input, alive only for
a fractional alt depth in `(20, 21)` with region depth ≤ 200 — which
`test_confidence_boundaries.py:322-349` pins precisely. Only 4 of the 6 writes decide any
producible row.

**The table must stay last-wins.** A first-match-wins table in source order restores `cond1`'s
region-depth cap, which #183 explicitly forbids (`test_confidence_boundaries.py:281`, `:352`).

### 6.3 — where the #223 check goes, and why

**In `run_kestrel`, never in `read_vcf_without_comments`.** That function's empty-frame fallback
is frozen in three places: `scripts/ble001_policy.json:929-936` (`fail_open[32]`, disposition
*"preserved-no-authorized-alternative"*), and `test_variant_parsing.py:335` and `:348`, which
assert the current behaviour as **correct**. Fixing it there would break the policy manifest and
collapse 6.3 into 6.5.

**The issue understates itself.** A headerless VCF *carrying real indel records* silently
discards them — demonstrated end-to-end: two well-formed indel rows passed the allele contract,
were correctly routed by `filter_indel_vcf`, and were then dropped, yielding `Negative` with no
exception and **zero ERROR log records**. That is a positive converted to a negative, not merely
an unreadable result. Also: the parse happens on the *derived* `output_insertion.vcf` /
`output_deletion.vcf`, not on Kestrel's VCF directly — a fix aimed at the wrong file would miss.

6.3 raises, consistent with `kestrel_genotyping.py:250` and `:274`. If 6.5 later adopts
record-and-continue, this raise converts to a recorded state and `vcf_missing_header` becomes
its reason code. Cost: 5 of 8 tests in `test_run_kestrel_skip.py` write `"fresh\n"` as a
successful VCF and need a real header line.

#224 is two fixes: the impossible-percentage clamp belongs here, and the **region validation is
a prerequisite for 6.6's #222** — it gives the harmonisation a guardrail. Measured:
`summarise_coverage([0,0,0,0], 2)` returns `percent_uncovered: 200.0`; a reversed region yields
`region_length: -999` with a **PASS** verdict; and `samtools depth -a -r chr1:1-250000000`
emitted 48 million rows and exhausted 30 GB before Python allocated anything — so `-a` made the
*file* unbounded too, and the bound belongs before the samtools call.

### 6.6 — #222, now diagnosed

This was parked as a domain question. It is diagnosable, and the answer points at an oversight
rather than a trade-off.

**The reference VNTR is 840 bp in both builds** — proved from this repo's own
`reference/vntr_db_advntr/{hg19,hg38}_muc1.db`, VNTR id 25561, 14 segments:

| build | window | left flank | VNTR | right flank | total |
| --- | --- | --- | --- | --- | --- |
| hg19 | 155160500-155162000 | 483 | 840 | 178 | 1501 |
| hg38 | 155188000-155192500 | 507 | 840 | **3154** | 4501 |

The left flanks agree to within 24 bp. **The entire ~3 kb asymmetry is the hg38 right flank**,
over-running the reference VNTR into downstream MUC1 sequence. The hg38 window's provenance is
unrecorded (the *start*, 155188507, is documented as found via UCSC BLAT at
`reference/README.md:458-462`; the window is not). Note the hg19 window is itself asymmetric,
so "harmonise to hg19" is not automatically right either.

Corrections to the issue: the coordinates are at `config.json:69` and `:74` (not `:46`/`:51`),
and there are **eight** legacy `vntr_region_*` keys at `:99-106` plus the two coordinate
entries — ten occurrences, not six. **A second definition site the issue misses entirely:**
`reference_registry.py:54-65`, headed *"Biological Truth — Single Source of Truth"*. Any
harmonisation must change both files or they diverge silently. The doctest asserting 1501 is
**not executed** (no `--doctest-modules`); the executed pin is
`test_coverage_stats.py:305-307`.

**Measured impact, from real golden-cohort depth files:** means move **+20% to +29%** under an
hg19-symmetric window and **−21% to −48%** under VNTR-only. `percent_uncovered` moves in **both
directions** across samples, so **no sample-independent correction factor exists** — historical
hg38 numbers cannot be back-corrected, only recomputed. That settles the issue's restatement
question. Downsampling is affected concretely on the Docker/web path, which always sets
`--advntr-max-coverage 300` (`docker/app/tasks.py:224`).

## M7 — Report as a faithful record

**Issues:** #242, plus the new `scoring-and-confidence.md` overclaim issue.

#242's measured claims were re-verified and **every behavioural claim holds**, including the
40-message `<br>` distribution (31 one, 9 two) and the contrast figures. But `1e08875` inserted
a 16-line block on 2026-08-11, so **every `report_template.html` citation past `:161` is stale
by exactly +16**, and `generate_report.py` by +12. The design and plan both need re-anchoring.

Three changes to #242's own decisions:

1. **P6 is answered** by M6.5's execution states — but only once
   `build_screening_summary` and the report actually read them. M7 depends on M6.5.
2. **P5 is two-thirds already satisfied** (see the corrections table). It needs a bump only for
   the resolved VNTR region and `--custom-regions`, which weakens the M6 coupling revision 1
   claimed.
3. **P1 resolves to no verdict word**, per the standing constraint. **The plan still
   contradicts this** — it adds a `verdict` to `ScreeningSummary` (`:7-9`, `:65-67`), titles a
   phase "Clinical safety" (`:57`), and includes a commit `feat(report): lead with the clinical
   verdict` (`:365`). An implementer following it would build the surface this design rejects.
   The plan must be revised before M7 is scheduled.

**Unacknowledged prerequisite:** the report design makes a browser tier mandatory and its plan
gates Tasks 2, 4 and 6 on it, but `pytest.ini:5-10` declares no browser marker and CI runs only
the unit gate. Three supposedly ungated tasks have **no executable acceptance environment**.
Browser-runner selection, dependency pinning, marker, Make target and CI job are an explicit M7
prerequisite.

Defects the issue missed, all worth folding in: offline, the Coverage QC verdict becomes
**unreachable entirely** (it lives inside a `display:none` block whose only reveal is in the
aborted script) — so the failure that discloses more variant rows also *removes* QC evidence;
the deliberate "no IGV data" fallback never renders either; against the even-row stripe the
contrast figures are worse than reported (`Low_Precision` 1.76:1); and `report_config.json`'s
`algorithm_logic` and the template filter **disagree about `"Not applicable"`**, so there is a
state where the narrative says flagged and the row is shown.

Also: the plan's Task 2 Step 2 is mis-sold as a security fix. `escape_frame_cells` already ran
at `generate_report.py:276`, so this is **not** a live XSS; the real difference is the
`if df.empty: return ""` guard.

## M8 — Cohort insight and access hardening

**Issues:** #33, #189, #227

### #189 — the fix is bigger than revision 1 said, and has an unnamed gap

Verified: `download_result` (`main.py:634-664`) touches **no Redis client at all**. Cohort
membership is forward-only (`sadd` at `main.py:488` and `tasks.py:373`); there is no
job→cohort reverse index anywhere in the key space.

Corrections and additions:

- The exposure is bounded at `MAX_RESULT_AGE_DAYS` (7 days), not unlimited. "No backfill, wait
  seven days" is therefore a defensible option, not negligence.
- **A backfill would brick legacy cohorts.** `resolve_cohort` refuses a record with an empty
  `hashed_passphrase` (`cohorts.py:411-415`); backfilling their members makes those downloads
  permanently fail with no recovery path. Skip them, or decide the bricking explicitly.
- **Ordering matters.** The `sadd` at `:488` sits outside the `try` whose `except` calls `srem`.
  Write the reverse key **before** the `sadd` so both failure modes are fail-closed, wrap both
  Redis calls (today the `srem` at `:534` can replace the original exception), and give the key
  a TTL extended by `extend_cohort_retention` — otherwise it outlives `cohort:<id>` and the
  member becomes permanently un-downloadable.
- **A new failure mode on a route that has none.** Adding a Redis read to `download_result`
  means deciding Redis-down behaviour, and the only safe answer — fail closed — stops **all**
  downloads, including non-cohort ones.
- **Unnamed gap:** cohort-*analysis* jobs get a fresh `analysis_job_id` with no cohort
  membership (`main.py:1056-1074`) while their archive is served by the same generic route
  (`tasks.py:541-603`). The derived archive stays retrievable by bare id.
- **Credential transport is undecided**, and it is the first decision, not an implementation
  detail. `online_mode.py:188-204` downloads without the passphrase it already holds; success
  emails carry only a bare URL (`tasks.py:376-396`); the `--resume` path (`:294-299`) has no
  cohort context persisted at all. Prefer a header or form flow over passphrases in URLs.

**A simpler alternative worth weighing:** write cohort members' archives under a cohort-scoped
path so `/download/{job_id}/` *physically cannot* serve them. Fail-closed by construction, no
reverse index, no TTL coupling, no backfill, no Redis on the download path. It costs per-member
download entirely, which is a bigger behaviour change — so it is better only if that capability
is expendable.

### #227 — one verified argument settles the choice

Confirmed unsalted and unkeyed (`cohort_pseudonyms.py:185`, `:204`), with nothing but the sample
name entering the digest.

The decisive fact the issue does not yet make: **stability is pinned by a test literal.**
`test_cohort_pseudonyms.py:45` asserts `pseudonymized_sample_name("anon_", "sample_one") ==
"anon_c788e939395d"`, and the golden gate compares on `Pseudonym`
(`scripts/golden_cohort/compare.py:70`). A random per-run salt breaks both. A per-installation
key breaks them **only if the key is present by default**. So option 2 — `VNTYPER_PSEUDONYM_KEY`,
absent by default, current behaviour as the documented fallback — is the **only** option that
does not force a golden-fixture regeneration and a test-literal rewrite.

The docs fix ships first and independently: `cli_parser.py:268-277` says "to protect sensitive
information", naming a protective purpose with no limit attached, and
`docs/user-guide/cohort-analysis.md:43-72` describes the mechanism without stating what it does
not protect against.

One qualifier for the issue body: the pseudonym is stable across cohorts **only for samples
whose local name is unique in every cohort they appear in** — `qualify_colliding_identities`
(`cohort_inputs.py:279`) rewrites colliding identities to `namespace/name` before hashing.

### #33 — needs a body rewrite, and "free" was overstated

The original mechanism is refuted: cohort mode loads only `pipeline_summary.json`, sees the
post-filter file, and gets exactly one Kestrel row per sample. Phase A — a cohort **call
frequency** table over final calls — is buildable today and was demonstrated end-to-end.

Three traps, all empirically confirmed: negative placeholder cells are the **literal string
`'None'`** (a naive `groupby` emits a `None` variant group at 0.667 frequency); every value
arrives as a **string** via `parse_tsv`; and the denominator cannot be `Sample.nunique()`
because a crashed Kestrel contributes zero rows while still counting as a sample.

**Revision 1's "Phase B falls out at no extra cost" is half true.** The *data path* into
`kestrel_df` genuinely is free. The *grouping semantics* are not: state propagation requires
editing `parse_pipeline_summary`, plus negative exclusion, numeric coercion and the denominator.

A new table must route through `escaped_table_html` with an **empty** `html_columns` — #207 is
closed and its fix is a per-column exemption list, so passing exemptions would reintroduce it.
And it must not end in a column named `Flag`: `cohort_summary_template.html:95-96`, `:137-138`
rewrite any table whose last header is literally `Flag` into a tooltip glyph.

## Issues needing a body rewrite before scheduling

- **#214** — both factual claims are stale and **the defect is inverted**. The workflow never
  ran on tags and no semver rule was ever present; meanwhile `publish-pypi.yml:736` already
  promotes a verified digest to `latest`, `2`, `2.0`, `2.0.13`. Live GHCR confirms all five
  aliases resolve to one digest. The live defect is that four surfaces still *say* the aliases
  do not exist yet (`README.md:133`, `installation.md:69-71`, `docker.md:13-15`,
  `docker/README.md:51-55`). **It is now a docs-only task**, its stated blocker is gone, and it
  is the only open issue in milestone #6 — the cheapest milestone close available.
  `test_release_workflow_contract.py:748-775` hard-asserts the stale prose and changes with it.
- **#20** — its audit comment must be **replaced**: the unsafe skip it cites was removed by
  `0a0af2a` on 2026-08-07. A planner reading it today would chase a fixed bug. The rewrite
  should ask for a *reader*, not a second manifest; `--resume` gated on the recorded `md5sum`
  still matching; refusal to resume across a `version` or `input_files` change; and it should
  resolve `--keep-intermediates`, which is now **an advertised no-op** — accepted at
  `fastq_bam_processing.py:91`, documented at `:109`, and never read. Prerequisite: only 5 of 10
  step names are constants; promote the inline literals first or a resume reader reintroduces
  AGENTS.md trap 5.
- **#33** — see above.

## Out of scope

- **#95** stays open (owner decision, 2026-08-12). For the record: 3 files, 484 lines, no Python,
  zero references anywhere, unchanged since 2025-02-18, invisible to ruff and coverage. Not in
  the wheel, but **present in the Docker image** (`docker/Dockerfile:46`), where `samtools` and
  `awk` exist so `count_variants.sh` would run while `Rscript` is absent. Proven not
  genotype-affecting. If it is ever deleted, salvage the pairwise-Fisher-across-a-cohort idea
  into #33 first — it is recorded nowhere else.
- Per-job bearer tokens. See the M8 alternatives.
- Tuned per-assay threshold values. #175 delivers a reachable and traceable operating point.

## Verification

Every claim in this document was established by reading current code, running it, or both.
Where a claim could not be settled statically it is marked as such: the #242 row-count and
rounding measurements need a browser, and `mediaRules: []` is an unsound measurement (a
cross-origin stylesheet throws `SecurityError`, which enumeration reports as empty) even though
the underlying claim — VNtyper ships no print CSS — is confirmed.

**Golden-cohort gate re-run is required for:** #222 (every hg38 case), #173's residual work, and
**any change that adds columns to `cohort_kestrel.{csv,tsv}`**. That last trigger was missing
from revision 1: `scripts/golden_cohort/artifacts.py:335-342` reads those headers, so a column
addition is a delta on every cohort case even though no row changes. `cohort_exports.py:14-21`
strips nothing by design, and `pd.DataFrame(list_of_dicts)` orders columns by first appearance —
so both presence and position become composition-dependent.

The gate is hand-run and not in CI.

## Traps carried forward

- `test_kestrel_filtering.py:73` reads `kestrel_genotyping.py` as source text and asserts
  **six** gate entries. **M6.5 must not alter the gate list** — residual classification reads
  pre-filter evidence without becoming a seventh gate.
- `test_cohort_summary_oracle.py` holds a whole-document fingerprint, but it is HTML-only and
  the display whitelists drop new columns, so it does **not** move on a column addition. The
  header contract most at risk is `:912`, which survives only because its fixture is hand-built.
- `MANIFEST.in:8-10` names exactly two templates, and `test_packaging_consistency.py:22-58`
  inspects only `.json` files while `make check-all` does not run `make build`. A new template
  or asset fails only at release time.
- `--config-path` replaces the whole config rather than merging it (AGENTS.md trap 2). Every new
  configurable value inherits that, including #227's digest width, which can already change
  silently today.
