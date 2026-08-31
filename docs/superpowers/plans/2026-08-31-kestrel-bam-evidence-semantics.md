# Kestrel BAM Evidence Semantics Implementation Plan

Date: 2026-08-31

Issue: [#295](https://github.com/hassansaei/VNtyper/issues/295)

Specification: [Kestrel BAM Evidence Semantics — Phase 1 Design](../specs/2026-08-31-kestrel-bam-evidence-semantics-design.md)

Scope: Phase 1 semantic containment. Preserve all naming, voting, tie, candidate-selection, threshold, refinement, and tier decisions. The pull request will say `Refs #295`, not close the issue.

## Execution rules

- Work only in the isolated `fix/issue-295-kestrel-bam-evidence-semantics` worktree based on `origin/main`.
- Execute tasks sequentially. Use one fresh implementer per task; implementers may not delegate.
- Before each test edit, state the realistic regression or deliberate production mutation it kills.
- For every production task: write the behavior test, run the named RED command and inspect the expected failure, make the smallest GREEN change, rerun, then refactor without changing behavior.
- Use Python 3.10-compatible syntax, Ruff's 120-column limit, Google docstrings, and existing logging/error conventions.
- New unit-test modules carry `pytestmark = pytest.mark.unit` immediately after imports.
- Each task receives one fresh read-only Opus 5 diff review with separate spec-compliance and code-quality verdicts. Return Critical/Important fixes to the same implementer and request only one scoped re-review when such a finding exists.
- Do not run overlapping implementation agents: Tasks 2–7 share nomenclature/report interfaces.

## Task 1: Publish the approved planning exception without weakening policy

**Files**

- Modify: `docs/superpowers/specs/2026-08-31-kestrel-bam-evidence-semantics-design.md`
- Add: `docs/superpowers/plans/2026-08-31-kestrel-bam-evidence-semantics.md`
- Modify: `mkdocs.yml`
- Modify: `AGENTS.md`
- Modify: `tests/unit/test_coverage_gate.py`
- Modify: `tests/unit/test_issue_233_documentation_contract.py`

**Mutation caught**

An unrestricted `docs/superpowers/` exception could silently republish arbitrary planning artifacts or let MkDocs macro delimiters break the site; an incomplete exception leaves the user-required pages red in the unit gate.

**RED**

1. Extend `test_contributor_docs_match_the_scripts_quality_scope` to accept only the two exact repository-relative paths from the specification, assert their two docs-relative nav paths and titles, and continue rejecting any third planning artifact.
2. Add an assertion that both published pages contain none of the three raw Jinja opening delimiters.
3. Update the #233 documentation-contract test so its prose and assertions recognize this bounded Phase 1 exception.
4. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_coverage_gate.py tests/unit/test_issue_233_documentation_contract.py -q`

   Expected failure before policy/prose changes: the current unconditional ban names `docs/superpowers`, the plan nav entry is absent, and the old #233 module contract says no such pages are published.

**GREEN**

1. Mark the specification `approved for implementation` and record that the final oversized Opus session was unavailable rather than claiming a verdict.
2. Add the plan to Development navigation as `Kestrel BAM Evidence Semantics Plan` immediately after the design page.
3. Replace the MkDocs planning comment with the exact two-path exception and retain the general `.planning/` rule.
4. Update `AGENTS.md` with the exact exception and later focused-module inventory destination without changing contributor policy elsewhere.
5. Run the RED command, then:

   `conda run --no-capture-output -n vntyper python -m mkdocs build --strict`

   Expected: targeted tests pass and strict docs build succeeds.

**Commit**

`docs(nomenclature): approve Kestrel evidence implementation plan`

## Task 2: Introduce source-aware evidence vocabulary and config compatibility

**Files**

- Add: `vntyper/scripts/nomenclature_evidence.py`
- Modify: `vntyper/scripts/nomenclature.py`
- Modify: `vntyper/scripts/nomenclature_config.json`
- Add: `tests/unit/test_nomenclature_evidence.py`
- Modify: `tests/unit/test_nomenclature_reconcile.py`

**Interfaces**

- Constants: `FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT`, `FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT`, `FLAG_LOW_KMER_PATH_SUPPORT`, `FLAG_LOW_READ_SUPPORT`, `FLAG_LOW_EVIDENCE_SUPPORT`.
- `resolve_bam_thin_haplotype_record_support(thresholds: Mapping[str, int]) -> int`.
- Pure source-unit selection used by `reconcile` to select truthful low-support flags while retaining its current all-known gate.
- `nomenclature.py` re-exports every authoritative `FLAG_*` string and retains `CALLER_OF`, source keys, and `MIN_SUPPORT_FOR_TIER_A` behavior.

**Mutations caught**

Tests must fail if Kestrel BAM support is relabeled as reads, if unknown sources are guessed to be reads, if one missing backing value no longer suppresses every low-support token, if a non-backing source adds a flag, if `kestrel_vcf` and `kestrel_bam` become independent callers, or if either threshold changes.

**RED**

1. Add resolver tests for new-only, legacy-only, both-present/new-wins, and neither-present/`KeyError`.
2. Add reconciliation tests for each source unit, mixed low known sources, unknown-source fallback, non-backing exclusion, scalar legacy `support=`, and the paired unknown/fully-known all-known cases.
3. Assert four BAM records are not thin at `< 3` but do receive the Tier-A `< 5` BAM token; keep exact boundary tests at 3 and 5.
4. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature_evidence.py tests/unit/test_nomenclature_reconcile.py -q`

   Expected failures: module/constants absent, config resolver absent, and `reconcile` emits the single misleading `low-read-support` token.

**GREEN**

1. Move only evidence flags, source-unit mapping, and the pure config resolver into the focused module; keep naming/canonicalization in `nomenclature.py`.
2. Change the shipped key to `bam_thin_haplotype_record_support: 3`; support legacy fallback only in the resolver.
3. Make `reconcile` preserve its current effective minimum and all-known calculation, adding truthful tokens only after the unchanged numeric comparison succeeds.
4. Re-export flags from `nomenclature.py`; document `supports` as heterogeneous source-specific evidence and the stable scalar legacy contract.
5. Run the RED command and the full reconciliation family:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature*.py -q`

**Refactor limit**

Do not change allele selection, `_is_corroborated`, `CALLER_OF`, threshold values, render syntax, or unrelated logic in the 1,170-line nomenclature module.

**Commit**

`fix(nomenclature): distinguish evidence support units`

## Task 3: Parse XD as optional observational BAM evidence

**Files**

- Modify: `vntyper/scripts/nomenclature_bam.py`
- Modify: `tests/unit/test_nomenclature_bam.py`
- Modify: `tests/unit/test_nomenclature_perf.py`

**Interfaces**

- `minimum_kmer_depth(record: pysam.AlignedSegment) -> int | None`.
- `BamConsensus.supporting_haplotype_records`, `fetched_haplotype_records`, `distinct_edit_count`, and keyword-only, equality-neutral `supporting_record_minimum_kmer_depths`.
- Read-only compatibility properties `support`, `total`, and `n_distinct`.
- Constant `THIN_HAPLOTYPE_RECORD_SUPPORT`; remove ambiguous internal `THIN_SUPPORT`.

**Mutations caught**

The BAM fixtures model resolved haplotypes, not duplicate raw reads. Parser cases kill tag coercion, missing-tag dropping, negative acceptance, overflow clipping, and subtype rejection. Metamorphic cases kill XD weighting, XD tie-breaking, XD-based filtering, and XD-dependent thinness.

**RED**

1. Replace the raw-read fixture helper with haplotype-record specifications carrying realistic `XD` values 5, 181, 7,416, and 8,704 where appropriate.
2. Add direct parser tests for valid compact integer subtypes, absent, wrong-type/non-integral, zero, negative, signed max, and unsigned-above-max. Use real `pysam.AlignedSegment` objects or indexed temporary BAMs, not mocks of decisions.
3. Add record-vote metamorphic tests: 3-versus-2 across changed XD states, unequal-XD tie still abstains, several weak-XD records beat one max-XD record, and an XD-sum tie does not override the record winner.
4. Assert the XD tuple aligns one-for-one with winning supporting records but equality/hash and `is_thin` ignore its values.
5. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature_bam.py tests/unit/test_nomenclature_perf.py -q`

   Expected failures: truthful fields/helper absent; existing fixtures and assertions use `reads`, `support`, `total`, and no XD evidence.

**GREEN**

1. Parse `XD` with `get_tag("XD", with_value_type=True)`; accept `cCsSiI`, preserve zero, accept through `2_147_483_647`, and return `None` for absent/invalid/negative/above-producer-range values. Log present invalid tags at debug level and let BAM decode/fetch failures follow existing boundaries.
2. Rename loop locals and logs to records/haplotype records.
3. Keep exactly one counter increment per record/edit key and the existing record-count tie rule; collect XD beside, never inside, the counter or comparison.
4. Build the canonical dataclass and compatibility properties; make the observational tuple keyword-only and `compare=False`.
5. Update `from_bam` to emit `thin-haplotype-record-support`, and `refine` to consume that exact token so its disagreement veto is unchanged.
6. Run the RED command, then:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature_bam.py tests/unit/test_nomenclature_perf.py tests/unit/test_nomenclature_reconcile.py -q`

**Commit**

`fix(nomenclature): type Kestrel haplotype evidence`

## Task 4: Propagate haplotype-record support without propagating XD into decisions

**Files**

- Modify: `vntyper/scripts/nomenclature_annotate.py`
- Modify: `tests/unit/test_nomenclature_surfaces.py`
- Modify: `tests/unit/test_nomenclature_reconcile.py`

**Interfaces**

- Rename `_row_read_call` to `_row_haplotype_call` and `_read_calls` to `_haplotype_calls`.
- `_row_haplotype_call` returns the named BAM call plus `supporting_haplotype_records`; it does not return or aggregate XD.
- `annotate_kestrel_frame`, `_open_rescuer`, `_row_verdicts`, and `reconcile_caller_outputs` retain current row order, eligibility, and output columns.

**Mutations caught**

A realistic surface test must fail if XD replaces support, if candidate filtering opens/fetches BAM for ineligible rows, if row order changes a winner, or if renamed thin support no longer blocks the existing `refine` disagreement path.

**RED**

1. Replace surface BAM fixtures with distinct Kestrel haplotype records and realistic XD tags.
2. Add a direct production-path assertion that a four-record BAM passes support `4`, not any XD statistic, to `supports["kestrel_bam"]`.
3. Parameterize identical VCF/adVNTR/BAM records over wildly different XD states and assert byte-equivalent nomenclature columns, flags, tier, name, and candidate/fetch counts.
4. Pin one rescuer per sample, per-row fetch order, no BAM access for non-candidates, and the thin disagreement-veto result.
5. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature_surfaces.py tests/unit/test_nomenclature_reconcile.py -q`

   Expected failures: old private names remain, fixtures lack XD, and BAM support still uses ambiguous field names.

**GREEN**

1. Perform the private-name/local/docstring/log terminology changes.
2. Thread only `supporting_haplotype_records` into the existing minimum support mapping; retain XD solely inside `BamConsensus`.
3. Preserve output column names, `_lesser`, row/candidate order, eligibility predicates, and TSV writing.
4. Run the RED command and:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature_kestrel.py tests/unit/test_nomenclature_surfaces.py tests/unit/test_nomenclature_reconcile.py -q`

**Commit**

`refactor(nomenclature): name BAM records truthfully`

## Task 5: Make report presentation source-aware

**Files**

- Add: `vntyper/scripts/nomenclature_presentation.py`
- Modify: `vntyper/scripts/report_formatting.py`
- Modify: `tests/unit/test_report_formatting.py`

**Interfaces**

- Focused owner for `NOMENCLATURE_TIERS`, `TIER_A_BLOCKERS`, `NOMENCLATURE_FLAG_MEANINGS`, and caller-neutral nomenclature column help.
- `report_formatting.py` re-exports existing presentation names so callers remain compatible.

**Mutations caught**

Direct table tests reject any Kestrel/shared entry that calls support reads, missing meanings for authoritative flags, a Tier-A blocker assigned to the nonblocking thin token, loss of legacy archived-token explanation, or removal of genuine adVNTR read wording.

**RED**

1. Assert all 14 authoritative flags are re-exported, explained, and present in the documentation contract.
2. Assert only reconciliation low-support tokens have Tier-A reasons; `thin-haplotype-record-support` has none.
3. Assert exact approved legacy wording, Kestrel VCF k-mer-path depth wording, Kestrel BAM resolved-haplotype wording, caller-neutral `ALT` help, and retained adVNTR `Supporting Reads` wording.
4. Add lexical guards over presentation tables rejecting `read`/`reads` for Kestrel/shared evidence, including `allele-unrepresentable-in-vcf`.
5. Render legends containing every new token so the dynamic paths are non-vacuous.
6. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_report_formatting.py -q`

   Expected failures: new tokens are absent and current tier/flag/help text universally describes reads.

**GREEN**

1. Extract only the pure nomenclature presentation tables/helpers from the 1,718-line formatting module.
2. Use the exact source-aware legacy explanations from the spec and non-clinical wording.
3. Retain function signatures, HTML escaping, table schemas, and current legend filtering/order.
4. Run the RED command plus:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_report_formatting.py tests/unit/test_report_presentation.py -q`

**Commit**

`fix(report): explain nomenclature evidence units`

## Task 6: Add the concise sample and cohort HTML clarification

**Files**

- Modify: `vntyper/templates/report_template.html`
- Modify: `vntyper/templates/cohort_summary_template.html`
- Modify: `vntyper/scripts/generate_report.py`
- Modify the cohort context owner found by inspection in `vntyper/scripts/cohort_summary.py` or its focused input/presentation module; do not duplicate decision logic.
- Modify: `tests/unit/test_generate_report.py`
- Modify: `tests/unit/test_report_igv_modes.py`
- Modify: `tests/unit/test_cohort_summary_escaping.py`
- Modify: `tests/unit/test_cohort_summary_oracle.py`
- Modify: `tests/unit/test_template_escaping.py`

**Mutation caught**

The rendered-report tests fail if the clarification disappears when IGV is off, if it is hidden behind a missing BAM, if XD is described as decisional, if the cohort injects unsafe markup, if the IGV accessible label says reads, or if the exported cohort schema changes.

**RED**

1. Render the real sample report in embedded, sidecar, and off modes and assert exactly one normalized occurrence of:

   “Kestrel output.bam contains resolved haplotype records, not sequencing reads. Its record counts are haplotype-record support; XD is minimum k-mer depth and does not weight votes or alter names or tiers.”

2. Require static `<code>` markup for `output.bam` and `XD`, ordinary autoescaping, and placement directly below the Reading key independent of IGV/BAM availability.
3. Assert the IGV region's accessible name says resolved haplotype-record alignments.
4. Render cohort reports with and without Kestrel nomenclature tokens. Require `show_kestrel_bam_semantics` and `nomenclature_legend` to produce the concise block after the complete `advntr_missing` conditional and before `additional_stats`; hostile values must be escaped.
5. Pin CSV/TSV/JSON columns unchanged and add the required Move 9 cohort skeleton fingerprint reason.
6. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_generate_report.py tests/unit/test_report_igv_modes.py tests/unit/test_cohort_summary_escaping.py tests/unit/test_cohort_summary_oracle.py tests/unit/test_template_escaping.py -q`

   Expected failures: exact text and context key absent, old IGV label remains, fingerprint differs once markup is added.

**GREEN**

1. Add the exact static sample text with no `safe` filter.
2. Derive the cohort boolean from actual Kestrel BAM evidence tokens and pass only plain strings/booleans to the template.
3. Add the conditional cohort block at the exact anchor and update the reasoned fingerprint history.
4. Do not change tables, charts, filtering, searching, paging, CSV, TSV, or JSON schemas.
5. Run the RED command, then:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_generate_report.py tests/unit/test_report_*.py tests/unit/test_cohort_summary*.py tests/unit/test_template_escaping.py -q`

**Commit**

`fix(report): clarify Kestrel BAM evidence in HTML`

## Task 7: Add repository-wide terminology containment tests

**Files**

- Add: `tests/unit/test_nomenclature_evidence_terminology.py`
- Modify code/docstrings only in already in-scope files named by failures: `vntyper/scripts/nomenclature_bam.py`, `vntyper/scripts/nomenclature_annotate.py`, `vntyper/scripts/nomenclature.py`, `vntyper/scripts/report_formatting.py`, and `vntyper/scripts/generate_report.py`.

**Mutation caught**

A deliberate restoration of `_row_read_call`, `_read_calls`, `reads` locals, “row's reads”, “reads per source”, a read-based tie log, or Kestrel/read wording in authoritative presentation tables must fail while legitimate FASTQ/input-BAM/adVNTR read terminology remains allowed.

**RED**

1. Inspect source text for the full functions `BamRescuer.rescue`, `from_bam`, `refine`, `annotate_kestrel_frame`, `_open_rescuer`, `_row_haplotype_call`, `_row_verdicts`, `_haplotype_calls`, and `reconcile`.
2. Assert old private names and banned Kestrel-BAM phrases are absent, all evidence flags are re-exported, and `nomenclature_evidence` flags are contained in the nomenclature public vocabulary.
3. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature_evidence_terminology.py -q`

   Expected failure: existing private/read terminology remains.

**GREEN**

Make only lexical in-scope corrections revealed by the targeted functions. Do not globally replace `read`, touch true sequencing-read paths, or refactor decisions.

Run the RED command and the complete nomenclature/report unit family.

**Commit**

`test(nomenclature): contain Kestrel BAM terminology`

## Task 8: Replace the vacuous golden approximation with exact display metrics

**Files**

- Modify: `tests/golden/test_nomenclature_golden.py`

**Mutation caught**

Exact per-tier equality catches added wrong Tier-B names, hidden lower-tier errors, tier migration with unchanged totals, a rendered/internal-name mismatch, altered BAM eligibility/fetching, support order drift, or any XD-driven decision.

**RED**

1. Replace `_reconciled` and `_hybrid` with one production-shaped replay implementing both policies from specification Section 2.2.
2. Define a metrics object whose displayed predicate is `call.name is not None and render(call) == call.name`; assert `displayed == exact + wrong` per tier and total.
3. Assert the full exact table: no BAM 138/119/19 with A 53/53/0, B 85/66/19, C 0/0/0; with current record voting 154/136/18 with A 53/53/0, B 101/83/18, C 0/0/0.
4. Assert 200 mutated and 200 normal samples, zero control findings, 83 eligible samples, 68 per-row fetches, all 178 Kestrel rows have usable loci, the exact per-class policy-2 counter, and gains `insCCCC +6`, `insG +4`, `insG_pos54 +7`.
5. Retain explicit non-vacuity checks for Tier A and BAM recovery and remove the obsolete 134-exact narrative.
6. Run with explicit roots:

   `VNTYPER_SIM_ROOT="${VNTYPER_SIM_ROOT:?set VNTYPER_SIM_ROOT to the simulation corpus}" VNTYPER_ADVNTR_ROOT="${VNTYPER_ADVNTR_ROOT:?set VNTYPER_ADVNTR_ROOT to the adVNTR corpus}" conda run --no-capture-output -n vntyper python -m pytest -m golden tests/golden/test_nomenclature_golden.py -q -rs`

   Expected RED while the old helpers remain: new metrics and exact fetch/per-class assertions are absent. After the test rewrite but before any accidental decision change, the exact established values must pass with zero skips.

**GREEN/refactor**

This task changes the oracle only. Share row parsing/replay code, preserve file order, count adVNTR support once per producing row, and avoid copying production selection expressions into expectations. If the exact baseline differs, stop and diagnose; do not update expected values.

**Commit**

`test(nomenclature): pin display-aware golden metrics`

## Task 9: Correct public documentation and changelog terminology

**Files**

- Modify: `docs/pipeline/nomenclature.md`
- Modify: `docs/pipeline/reports.md`
- Modify: `docs/user-guide/output-files.md`
- Modify: `docs/cli/report.md`
- Modify: `docs/getting-started/quickstart.md`
- Modify as consistency requires: `docs/pipeline/kestrel.md`, `docs/pipeline/scoring-and-confidence.md`
- Modify: `vntyper/scripts/README.md`
- Modify: `docs/about/changelog.md`
- Modify: `AGENTS.md`
- Modify: `tests/unit/test_issue_233_documentation_contract.py`
- Modify: `tests/unit/test_nomenclature_evidence_terminology.py`

**Mutation caught**

Documentation guards fail if a Kestrel `output.bam` count is called reads, if XD is called a read count or made decisional, if the published flag vocabulary drifts from all 14 tokens, if the measured 83/68 figures regress to “about a fifth”, or if released changelog history is rewritten.

**RED**

1. Add documentation assertions for resolved haplotype records, haplotype-record support, minimum k-mer depth, explicit non-weighting, exact flag/config names, 83 eligible/68 fetched, and the 14-token table.
2. Assert genuine adVNTR/input sequencing-read language remains.
3. Run:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_issue_233_documentation_contract.py tests/unit/test_nomenclature_evidence_terminology.py -q`

   Expected failures: current public pages say reads/alignments ambiguously and document only 8 flags.

**GREEN**

1. Correct only Kestrel BAM semantic claims; explain the heterogeneous unchanged threshold without calibration claims.
2. Document config fallback, flag compatibility, HTML clarification, and source units exactly.
3. Update only `Unreleased`; state `Refs #295`, Phase 1, and non-goals. Do not rewrite released entries.
4. Update `AGENTS.md` focused-module count from 16 to 18 and list `nomenclature_evidence.py` and `nomenclature_presentation.py`.
5. Run the RED command, strict docs build, and a raw Jinja-delimiter scan over both published superpowers pages.

**Commit**

`docs(nomenclature): correct Kestrel BAM evidence terms`

## Task 10: Integrated behavior-preservation regression pass

**Files**

- Modify only failing in-scope tests or code already listed above; no new interface or behavior is authorized.

**Mutations caught**

This task catches cross-module breakage that isolated RED/GREEN loops miss: compatibility-property consumers, report templates using old tokens, untested modified functions, fixture drift, and source ordering differences.

**Verification and controlled fixes**

1. Run targeted families:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature_bam.py tests/unit/test_nomenclature_evidence.py tests/unit/test_nomenclature_reconcile.py tests/unit/test_nomenclature_surfaces.py tests/unit/test_nomenclature_kestrel.py tests/unit/test_report_formatting.py tests/unit/test_generate_report.py tests/unit/test_report_igv_modes.py tests/unit/test_cohort_summary_escaping.py tests/unit/test_cohort_summary_oracle.py -q`

2. Run complete nomenclature unit family:

   `conda run --no-capture-output -n vntyper python -m pytest tests/unit/test_nomenclature*.py -q`

3. For each failure, apply systematic debugging: reproduce, trace to the earliest incorrect value, and fix only the approved interface. Add a mutation-killing regression before the fix.
4. Run `git diff --check` and inspect `git diff --stat origin/main...HEAD` for unrelated files.
5. Commit only if fixes were required:

   `fix(nomenclature): preserve evidence decision behavior`

## Final verification

Run fresh commands from the worktree root and preserve complete outputs. A skip or missing dependency is recorded as unavailable evidence, not a pass.

1. Targeted nomenclature/report commands from Task 10.
2. Golden test with explicit roots from Task 8; require 200 mutated, 200 normal, zero skips, and the exact metrics.
3. Survey all 200 BAMs and require 1,327,794 integral XD values, none missing/malformed, min 5, median 181, p95 7,416, max 8,704. This is corpus evidence, not calibration.
4. `make format-check`
5. `make lint`
6. `make type-check-all`
7. `make test-unit`
8. `make check-integration-compatibility`
9. `make patch-coverage`
10. `make docs-build`
11. `make test-browser`
12. `make check-all`
13. Scan both published design/plan files for raw Jinja opening delimiters and run `git diff --check`.

Do not run `make ci-local`; workflow files are out of scope and must remain unchanged.

## Final review and branch finish

1. Smoke-test canonical Opus 5 without repository tools, then review bounded core, report/docs, and golden/spec diff partitions supplied on standard input with tools disabled. Require the combined review to try to prove decision drift, XD coupling, stale Kestrel read language, mirrored test logic, or omitted lower-tier wrong names.
2. If and only if it reports validated Critical/Important findings, make one controlled fix wave and rerun affected targeted tests, golden tests when decision surfaces changed, and `make check-all`. At the owner's direction, do not start another review loop after that verified fix wave.
3. Invoke verification-before-completion: compare every claim with fresh command output and report unavailable external evidence separately.
4. Invoke finishing-a-development-branch and prepare, but do not merge, a PR titled `fix(nomenclature): correct Kestrel BAM evidence semantics`.
5. The PR body must include `Refs #295`; Phase 1 scope/non-goals; before/after terms; XD parser contract; exact unchanged golden table; every verification result; any unavailable evidence; Opus findings/dispositions; and follow-up order `#270 -> #267 -> joint #295/#269 calibration`.
