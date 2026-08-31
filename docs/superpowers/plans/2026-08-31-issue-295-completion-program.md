# VNtyper Issue 295 Completion Program Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver typed molecular identity, governed adVNTR evidence, complete explicit decision profiles, and future opt-in dominance/abstention evaluation in four sequential pull requests while leaving every shipped cutoff unchanged.

**Architecture:** Pure typed modules sit between caller artifacts and presentation, while the existing stage files retain I/O and orchestration. Each pull request starts from newly fetched `origin/main`, proves its compatibility projection independently, receives bounded adversarial review, and merges before the next branch is created. The packaged policy remains `legacy-selection-v1`; generated profiles are never installed or selected automatically.

**Tech Stack:** Python 3.10+, pandas, pytest, dataclasses, strict JSON with RFC 8785 canonicalization, SHA-256, Jinja2/MkDocs, Playwright, GitHub CLI, and the existing golden corpora.

---

## Program boundaries

This is one execution plan with four functional PRs, not one combined branch. PR A is the
only branch that exists at the start. After each green PR merges under the owner's standing
authorization, delete its worktree, fetch `origin`, and create the next worktree from the
new `origin/main`. Never stack B on unmerged A, C on unmerged B, or D on unmerged C.

The following values are immutable throughout all four PRs:

```text
Depth_Score low/high                 0.00469 / 0.00515
alternate k-mer-path boundaries     20 / 21 / 100
active-region boundary              200
GG floor                            0.00469
BAM flank                           8
thin BAM haplotype-record support   3
Tier-A Kestrel k-mer-path depth     5
Tier-A adVNTR sequencing reads      5
```

The reporting floor from #266 and the model/version guard from #268 are outside every
custom or generated decision surface. Haplotype records remain records; `XD` remains
optional minimum k-mer depth. PR D is future tooling for #269 and does not alter package
defaults. The available 200 mutated plus 200 paired-normal simulations are development
evidence, so PR D ends with `closure_eligible: false`; #295 remains open and release work
does not start.

## Controller protocol used for every PR

- [ ] Read the relevant task and test mutation aloud in the implementation update before
  running its RED command.
- [ ] Start every new unit file with `pytestmark = pytest.mark.unit`, and add a behavioral
  regression test for every modified production function before changing it.
- [ ] Dispatch one fresh implementer for the bounded task. Implementers may not spawn
  agents, and two implementers never edit concurrently.
- [ ] Have the controller inspect every changed production function and its behavioral
  test before accepting the commit.
- [ ] Run a spec-compliance review first and a code-quality review second, then resolve
  Critical and Important findings in one fix wave.
- [ ] Partition tools-disabled Claude Code Opus 5 review into core code, public surfaces,
  and tests/golden validity. Use the already smoke-tested command:

```bash
claude -p --model claude-opus-5 \
  --tools "" --no-session-persistence --output-format json
```

- [ ] Run one separate GPT-5.6-sol high adversarial review of decision invariants and
  non-vacuous tests. Allow one scoped re-review after the single fix wave; record and
  disposition repeated invalid findings instead of looping.
- [ ] Before opening the PR, run the PR-specific tests followed by the common gate:

```bash
make format-check
make lint
make type-check-all
make test-unit
make test-unit-cov
make check-integration-compatibility
make patch-coverage
make docs-build
make check-all
```

- [ ] Run `make test-browser` when HTML behavior changes. Run `make check-full` only when
  the required integration archive and reference are present. Run advisory mutation tests
  for new decision modules. Run `make ci-local` only if a separately approved workflow
  change exists.
- [ ] Inspect complete command output and report failures, five known unit skips or any
  changed skip count, warnings, coverage, corpus counts, and exact by-tier metrics.
- [ ] Open a focused PR with scope, non-goals, review dispositions, commands/results, and
  unavailable validation. Use only `Refs #N` footers. Wait for GitHub checks; merge without
  squash only when all gates are green.

## PR A — issue 270 molecular identity

**Branch:** `fix/issue-270-molecular-identity`  
**PR title:** `feat(identity): add canonical molecular identity`  
**Issue footer:** `Refs #270`

### PR A file map

Create:

- `vntyper/scripts/molecular_identity.py` — validated value types and stable serialization.
- `vntyper/scripts/molecular_identity_translation.py` — complete motif-pair translation.
- `vntyper/scripts/identity_candidates.py` — pre-projection capture and identity grouping.
- `vntyper/scripts/identity_reconciliation.py` — identity-keyed backing and Tier logic.
- `vntyper/scripts/molecular_identity_presentation.py` — public field contract.
- `vntyper/scripts/nomenclature_bam_evidence.py` — complete replayable BAM record evidence.
- `tests/unit/test_molecular_identity.py`
- `tests/unit/test_molecular_identity_translation.py`
- `tests/unit/test_identity_candidates.py`
- `tests/unit/test_identity_reconciliation.py`
- `tests/unit/test_identity_bam_binding.py`
- `tests/unit/test_molecular_identity_surfaces.py`
- `tests/unit/test_molecular_identity_presentation.py`
- `tests/golden/identity_oracle.py`
- `tests/golden/test_molecular_identity_golden.py`

Modify only as thin adapters:

- `vntyper/scripts/motif_processing.py`
- `vntyper/scripts/kestrel_genotyping.py`
- `vntyper/scripts/nomenclature.py`
- `vntyper/scripts/nomenclature_annotate.py`
- `vntyper/scripts/nomenclature_bam.py`
- `vntyper/modules/advntr/advntr_genotyping.py`
- `vntyper/scripts/advntr_variant_annotations.py`
- `vntyper/scripts/report_formatting.py`
- `vntyper/scripts/cohort_tables.py`
- `vntyper/scripts/summary.py`
- `vntyper/scripts/generate_report.py`
- `vntyper/templates/report_template.html`
- `vntyper/templates/cohort_summary_template.html`
- `tests/builders.py`
- `tests/unit/test_kestrel_filtering.py`
- `tests/unit/test_kestrel_result_publication.py`
- `tests/unit/test_nomenclature_reconcile.py`
- `tests/unit/test_nomenclature_bam.py`
- `tests/unit/test_report_formatting.py`
- `tests/unit/test_cohort_tables.py`
- `tests/unit/test_generate_report.py`
- `docs/pipeline/nomenclature.md`
- `docs/user-guide/output-files.md`

Do not place pure logic in the 1,200-line `nomenclature.py`, 1,109-line
`kestrel_genotyping.py`, 1,532-line `report_formatting.py`, or 1,100-line
`generate_report.py`.

### Task A0: Freeze the independent golden oracle before production changes

**Production mutation caught:** the identity implementation and expected-value calculation
share normalization/grouping logic, a lower tier disappears, or the historical Phase 1
projection is rewritten during migration.

**Files:** create `tests/golden/identity_oracle.py` and
`tests/golden/test_molecular_identity_golden.py`.

- [ ] Implement the independent truth/canonical-sequence oracle and recursively AST-scan
  its complete local import closure. Reject imports of production identity,
  canonicalization, serialization, grouping, reconciliation, and decision code.
- [ ] Write literal historical assertions for 200 mutated/200 controls, `154/136/18`, Tier
  A `53/53/0`, zero control findings, 83 BAM-eligible, 68 fetches, and dupC `96/96` named
  `59dupC`; add per-sample expected molecular identity/name and literal tier keys.
- [ ] Run the explicit-root golden command. Expect RED because current caller rows lack the
  identity quartet, while all historical Phase 1 assertions pass and zero tests skip.
- [ ] Commit only the oracle and failing contract after recording the RED output in the
  execution log: `test(golden): define molecular identity oracle`.

### Task A1: Canonical identity value model

**Production mutation caught:** an insertion uses ordinary interval coordinates, identity
construction accepts lowercase or overlapping edits, or serialization changes canonical
dupC away from `MUC1-X-60-coding-v1|60|59|-|C`.

**Files:** create `vntyper/scripts/molecular_identity.py`; create
`tests/unit/test_molecular_identity.py`.

- [ ] Write the RED tests, including the unit marker:

```python
from collections.abc import Callable

import pytest

from vntyper.scripts.molecular_identity import (
    CodingEdit,
    IdentityTranslation,
    MolecularIdentity,
    make_coding_edit,
    make_molecular_identity,
    parse_molecular_identity,
    serialize_molecular_identity,
)

pytestmark = pytest.mark.unit


def test_canonical_dupc_has_stable_identity() -> None:
    identity = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))
    assert serialize_molecular_identity(identity) == "MUC1-X-60-coding-v1|60|59|-|C"


def test_empty_identity_is_rejected() -> None:
    with pytest.raises(ValueError, match="Molecular identity requires at least one edit"):
        make_molecular_identity(())


def test_direct_invalid_insertion_is_rejected() -> None:
    with pytest.raises(ValueError, match="Insertion coordinates require end == start - 1"):
        CodingEdit(60, 60, "", "C")


@pytest.mark.parametrize(
    "factory",
    [
        lambda: CodingEdit(60, 59, "", "c"),
        lambda: CodingEdit(61, 60, "", "C"),
    ],
)
def test_case_and_range_are_strict(factory: Callable[[], CodingEdit]) -> None:
    with pytest.raises(ValueError):
        factory()


def test_foreign_identity_version_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unsupported molecular identity reference"):
        parse_molecular_identity("MUC1-X-60-coding-v2|60|59|-|C")


def test_translation_state_is_consistent() -> None:
    with pytest.raises(ValueError):
        IdentityTranslation(identity=None, status="resolved", failure=None, context_diverges=False)
```

- [ ] Run `python -m pytest -q tests/unit/test_molecular_identity.py` in the activated
  `vntyper` environment.
  Expect collection failure because the module does not exist.
- [ ] Implement frozen `CodingEdit`, `MolecularIdentity`, `IdentityTranslation`,
  `FrameConsequence`, `KestrelRepresentation`, `AdvntrRepresentation`,
  `IdentityObservation`, `EvidenceDisposition`, and `IdentityDecision`. Restrict alleles to uppercase `[ACGT]*`, positions to 1–60,
  insertion coordinates to `end == start - 1`, and edits to sorted non-overlapping tuples.
  Use this stable serializer:

```python
def serialize_molecular_identity(identity: MolecularIdentity) -> str:
    def allele(value: str) -> str:
        return value or "-"

    encoded = ";".join(
        f"{edit.start}|{edit.end}|{allele(edit.deleted)}|{allele(edit.inserted)}"
        for edit in identity.edits
    )
    return f"{identity.reference}|{encoded}"
```

- [ ] Run the test file again and expect all tests to pass. Run Ruff on both files.
- [ ] Commit:

```bash
git add vntyper/scripts/molecular_identity.py tests/unit/test_molecular_identity.py
git commit -m "feat(identity): add canonical molecular identity values"
```

### Task A2: Complete-context caller translation

**Production mutation caught:** translation validates only edited bases, omits the
neighbouring motif anchor, reverse-complements incorrectly, or makes a wrong anchor equal
to a resolved identity.

**Files:** create `vntyper/scripts/molecular_identity_translation.py`; create
`tests/unit/test_molecular_identity_translation.py`; modify compatibility projection
functions in `vntyper/scripts/nomenclature.py` without moving display behavior.

- [ ] Write table-driven RED tests for `G>GG` and `C>CG` canonical dupC, homopolymer
  anchors, distinct dupA, invalid allele, absent motif context, wrong neighbour anchor,
  pair-boundary edit, non-X repeat unit, and reconstruction mismatch. Pin closed failures:

```python
def test_wrong_neighbour_anchor_is_unresolved(canonical_pair: str) -> None:
    representation = kestrel_representation(
        replace_neighbour(canonical_pair, "A"), motifs="known-pair-key"
    )
    result = translate_kestrel_representation(representation, TEST_MOTIF_MAP)
    assert result.identity is None
    assert result.status == "unresolved"
    assert result.failure == "motif-anchor-mismatch"
```

  Add RED cases for `translate_advntr_representation()` using parsed RU/position State
  changes, insufficient bare State context returning unresolved, and
  `bind_bam_translation()` preserving the complete co-occurring edit set.

- [ ] Run the new file and expect import failure for the translation module.
- [ ] Implement `translate_kestrel_representation()`,
  `translate_advntr_representation()`, `bind_bam_translation()`, and the internal
  `resolve_coding_pair_edit()`. Reconstruct the full 120-bp pair, orient the 60-bp X unit,
  validate the neighbour against the checked-in motif map, 3-prime normalize in coding
  orientation, and reapply the edit to reproduce the complete alternate X unit.
- [ ] Run the new translation tests plus `tests/unit/test_nomenclature.py` and expect green
  display compatibility.
- [ ] Commit:

```bash
git add vntyper/scripts/molecular_identity_translation.py \
  vntyper/scripts/nomenclature.py tests/unit/test_molecular_identity_translation.py
git commit -m "feat(identity): translate complete caller representations"
```

### Task A3: Capture identity before legacy motif projection

**Production mutation caught:** two rows with the same `POS/REF/ALT` and different
`Motifs` collapse before translation, equivalent-row support is summed, or a new sort
changes the selected row.

**Files:** create `vntyper/scripts/identity_candidates.py`; create
`tests/unit/test_identity_candidates.py`; modify `vntyper/scripts/motif_processing.py` and
`vntyper/scripts/kestrel_genotyping.py`; modify `tests/unit/test_kestrel_filtering.py`.

- [ ] Write RED tests using two motif-distinct raw rows plus two equivalent shapes:

```python
def test_same_raw_tuple_under_two_motifs_retains_two_hypotheses() -> None:
    captured = capture_kestrel_observations(SAME_TUPLE_DIFFERENT_MOTIFS, TEST_COMPONENT)
    candidates = overlay_legacy_projection(captured, LEGACY_PASSING_KEYS, LEGACY_SELECTED)
    assert candidates.identity_hypothesis_count == 2
    assert candidates.selected_row_key == LEGACY_SELECTED


def test_equivalent_representations_do_not_sum_support() -> None:
    candidates = build_equivalent_candidates(support=(4, 7), selected_index=0)
    assert candidates.selected_support == 4
    assert candidates.equivalent_representation_count == 2


def test_group_union_cannot_unblock_tier_a() -> None:
    candidates = build_equivalent_candidates(context_diverges=(False, True), selected_index=0)
    assert candidates.selected_group.context_diverges is True
```

  Add an adVNTR RED case proving two positive State rows are translated before
  reconciliation and receive caller-local hypothesis count `2`.

- [ ] Run the new file and the exact legacy selection tests. Expect missing symbols while
  the existing selection tests remain green.
- [ ] Implement `RawRepresentationKey`, `IdentityCandidate`, `IdentityCandidateSet`,
  `capture_kestrel_observations()`, `capture_advntr_observations()`,
  `overlay_legacy_projection()`, and group-level unions
  of blocking gates, flags, and context divergence. Capture after confidence calculation and before
  `add_haplo_count()`. Overlay eligibility only after all six existing gates. Keep
  `select_single_best_variant()` byte-for-byte unchanged.
- [ ] Run `test_identity_candidates.py`, `test_kestrel_filtering.py`, and
  `test_haplo_count_and_selection.py`; expect green with identical legacy ordering and
  ties.
- [ ] Commit:

```bash
git add vntyper/scripts/identity_candidates.py vntyper/scripts/motif_processing.py \
  vntyper/scripts/kestrel_genotyping.py tests/unit/test_identity_candidates.py \
  tests/unit/test_kestrel_filtering.py
git commit -m "feat(kestrel): retain identity hypotheses across legacy projection"
```

### Task A4: Identity-keyed reconciliation and evidence disposition seam

**Production mutation caught:** reconciliation groups by displayed name, unresolved rows
agree by event class, Kestrel VCF and BAM count as separate callers, or Tier-A support is
borrowed from a non-backing source.

**Files:** create `vntyper/scripts/identity_reconciliation.py`; create
`tests/unit/test_identity_reconciliation.py`; modify `vntyper/scripts/nomenclature.py` and
`vntyper/scripts/nomenclature_annotate.py`; modify
`tests/unit/test_nomenclature_reconcile.py`.

- [ ] Write RED tests for display-equal/identity-different, display-different/identity-equal,
  unresolved, tied, divergent-context, Kestrel VCF+BAM, and unchanged PR-A
  `Polymorphic_Call` admissibility cases. Consume the production `EvidenceDisposition`
  seam defined in Task A1:

```python
def test_display_equal_identity_different_cannot_agree() -> None:
    decision = reconcile_identity_observations(DISPLAY_EQUAL_IDENTITY_DIFFERENT, POLICY)
    assert decision.molecular_agreement is False
    assert decision.tier != "A"
```

- [ ] Run the new file and expect missing identity reconciler.
- [ ] Copy the current source/support/tie behavior into the focused module, replacing only
  `str(call.name)` backing keys with resolved `MolecularIdentity`. Keep source-bound
  Kestrel alternate k-mer-path depth and adVNTR sequencing-read support as separately
  named quantities. Reject unresolved backing; retain context-divergent representation
  evidence at lower tiers while applying the existing motif-status Tier-A blocker; count
  Kestrel VCF and BAM as one source.
- [ ] Run both reconciliation files and nomenclature surface tests; expect every legacy
  case green except the intentional identity-key correction.
- [ ] Commit:

```bash
git add vntyper/scripts/identity_reconciliation.py vntyper/scripts/nomenclature.py \
  vntyper/scripts/nomenclature_annotate.py tests/unit/test_identity_reconciliation.py \
  tests/unit/test_nomenclature_reconcile.py
git commit -m "fix(nomenclature): reconcile callers by molecular identity"
```

### Task A5: Bind BAM evidence and retain replay data

**Production mutation caught:** BAM increments equivalent-representation count, becomes a
second caller, weights votes by XD, synthesizes a delins from separate VCF rows, or drops
runner-up evidence needed by future replay.

**Files:** create `vntyper/scripts/nomenclature_bam_evidence.py`; create
`tests/unit/test_identity_bam_binding.py`; modify `vntyper/scripts/nomenclature_bam.py` and
`tests/unit/test_nomenclature_bam.py`.

- [ ] Write RED tests proving one vote per identity per eligible resolved record,
  co-occurring identities can both receive that record, the denominator contains every
  eligible resolved record, missing XD remains `None`, and binding cannot alter caller or
  representation counts:

```python
def test_record_votes_once_for_each_cooccurring_identity() -> None:
    evidence = collect_locus_evidence([RECORD_WITH_DUPLICATE_AND_SECOND_EDIT], LOCUS, COMPONENT)
    assert evidence.counts == {DUPC: 1, SECOND_EDIT: 1}
    assert evidence.eligible_record_count == 1
    assert evidence.xd_by_record == (None,)
```

- [ ] Run the new file and expect the evidence module to be absent.
- [ ] Implement frozen `BamIdentityEvidence` and `BamLocusEvidence` plus
  `collect_locus_evidence()`. Add a compatibility adapter around `BamConsensus`; do not
  edit Phase 1 eligibility, unweighted vote, tie, thin-support `3`, flank `8`, VCF-primary
  refinement, or XD semantics. Persist the complete eligible record identity sets and
  typed XD values for PR D replay.
- [ ] Run new and existing BAM tests; explicitly verify that unequal XD cannot break a
  record-count tie and `pair_4092` receives no fabricated delins identity.
- [ ] Commit:

```bash
git add vntyper/scripts/nomenclature_bam_evidence.py \
  vntyper/scripts/nomenclature_bam.py tests/unit/test_identity_bam_binding.py \
  tests/unit/test_nomenclature_bam.py
git commit -m "feat(identity): bind bam evidence to kestrel representations"
```

### Task A6: Publish identity fields without widening negative output

**Production mutation caught:** a negative Kestrel or adVNTR row gains identity columns,
an unresolved row emits a serialization, a resolved row reports representation count
zero, or a below-floor candidate enters the hypothesis count.

**Files:** create `vntyper/scripts/molecular_identity_presentation.py`; create
`tests/unit/test_molecular_identity_surfaces.py`; modify `tests/builders.py`,
`vntyper/scripts/kestrel_genotyping.py`, `vntyper/modules/advntr/advntr_genotyping.py`,
`vntyper/scripts/nomenclature_annotate.py`, and
`tests/unit/test_kestrel_result_publication.py`.

- [ ] Write RED assertions for the exact quartet and unchanged ten-column negative schema:

```python
IDENTITY_COLUMNS = (
    "Molecular_Identity",
    "Molecular_Identity_Status",
    "Equivalent_Representation_Count",
    "Identity_Hypothesis_Count",
)


def test_unresolved_positive_row_has_consistent_cells() -> None:
    cells = identity_result_cells(UNRESOLVED_SELECTED, RESOLVED_OTHER)
    assert cells == {
        "Molecular_Identity": "",
        "Molecular_Identity_Status": "unresolved",
        "Equivalent_Representation_Count": 0,
        "Identity_Hypothesis_Count": 1,
    }
```

- [ ] Run the focused result-publication tests and expect missing columns/helpers.
- [ ] Implement `identity_result_cells()` in `molecular_identity_presentation.py`. Assert
  `unique`, `legacy-selected-among-multiple`, and `unresolved`, including a resolved unique
  row with `Equivalent_Representation_Count == 1`. Append the quartet only to emitted positive caller tables. Repeat caller-local
  hypothesis count on each positive row. Add only the approved four translation diagnostic
  fields `Molecular_Identity`, `Molecular_Identity_Translation_Status`,
  `Molecular_Identity_Translation_Failure`, and
  `Molecular_Identity_Context_Diverges` to `kestrel_pre_result.tsv`. Preserve negative
  schemas and the single `##` note. Add a below-floor fixture whose identity does not enter
  `Identity_Hypothesis_Count`.
- [ ] Run the focused files and #266 tests; expect green.
- [ ] Commit:

```bash
git add vntyper/scripts/molecular_identity_presentation.py \
  vntyper/scripts/kestrel_genotyping.py vntyper/modules/advntr/advntr_genotyping.py \
  vntyper/scripts/nomenclature_annotate.py \
  tests/builders.py tests/unit/test_molecular_identity_surfaces.py \
  tests/unit/test_kestrel_result_publication.py
git commit -m "feat(outputs): publish molecular identity fields"
```

### Task A7: Report, cohort, and schema-2 compatibility

**Production mutation caught:** a public export omits one identity field, old summaries
render blank or infer identity from legacy cells, or PR A prematurely changes summary
schema version.

**Files:** create `tests/unit/test_molecular_identity_presentation.py`; modify
`vntyper/scripts/molecular_identity_presentation.py`, `report_formatting.py`,
`cohort_tables.py`, `summary.py`, `generate_report.py`, both report templates, related unit
tests, `docs/pipeline/nomenclature.md`, and `docs/user-guide/output-files.md`.

- [ ] Write RED tests that all sample/cohort HTML/TSV/CSV/JSON result tables carry the four
  exact labels, schema 1/2 rows without them render `legacy identity not recorded`, rows
  with the quartet render their values, and `SUMMARY_SCHEMA_VERSION == 2`.
- [ ] Run the listed report/cohort tests and expect the new contract assertions to fail.
- [ ] Add presentation helpers and thin surface wiring. Add
  `decision_policy: legacy-selection-v1` to current summary provenance without changing
  schema version. Do not invent or revise configuration-owned clinical wording.
- [ ] Run unit report tests and `make test-browser`; expect online/offline parity and all
  identity fields visible.
- [ ] Commit:

```bash
git add vntyper/scripts/molecular_identity_presentation.py \
  vntyper/scripts/report_formatting.py vntyper/scripts/cohort_tables.py \
  vntyper/scripts/summary.py vntyper/scripts/generate_report.py vntyper/templates \
  tests/unit/test_molecular_identity_presentation.py tests/unit/test_report_formatting.py \
  tests/unit/test_cohort_tables.py tests/unit/test_generate_report.py \
  docs/pipeline/nomenclature.md docs/user-guide/output-files.md
git commit -m "feat(report): expose molecular identity provenance"
```

### Task A8: Independent golden identity oracle

**Production mutation caught:** expected identity is calculated by importing the production
normalizer, a structural delins is fabricated, a lower tier is omitted, or selected row
order/tie behavior changes unnoticed.

**Files:** modify `tests/golden/identity_oracle.py` and
`tests/golden/test_molecular_identity_golden.py` created in Task A0.

- [ ] Confirm the Task-A0 oracle's recursive import-closure guard still passes and add
  structural-delins rejection plus paired positive controls proving the BAM recovery path
  executed.
- [ ] Write the golden test to assert both roots exist and were loaded, exactly 200 mutated
  and 200 controls entered, and the historical Phase 1 reference fixture is
  `154/136/18`, Tier A `53/53/0`, zero controls,
  83 BAM-eligible, 68 fetches, and canonical dupC `96/96` named `59dupC`. Assert the exact
  Kestrel selected-row projection and report PR-A deltas across every tier. Count the seven
  structural delins carriers and `pair_4092` as misses unless unchanged Phase 1 BAM evidence
  truly represents the edit. Assert current PR-A metrics against the independent oracle;
  they may differ only through identity-keyed reconciliation, while the Kestrel selected-row
  projection remains exact.
- [ ] Before updating any measured aggregate, apply each named mutation to a temporary
  worktree and observe the pre-existing oracle fail on the affected sample/field. Then run
  the unmutated branch with explicit roots, never accepting a skip:

```bash
: "${VNTYPER_SIM_ROOT:?set VNTYPER_SIM_ROOT to the simulation corpus}"
: "${VNTYPER_ADVNTR_ROOT:?set VNTYPER_ADVNTR_ROOT to the adVNTR corpus}"
python -m pytest -q -m golden tests/golden -rs
```

- [ ] Record only metrics derived by the already independent oracle, rerun, and require all
  golden tests pass with zero skips. Record the measured PR-A identity-policy baseline in
  the PR; do not invent it in advance.
- [ ] Commit:

```bash
git add tests/golden/identity_oracle.py tests/golden/test_molecular_identity_golden.py
git commit -m "test(golden): pin molecular identity policy baseline"
```

### Task A9: Verify, review, and merge PR A

- [ ] Run all PR A targeted tests, common gates, explicit-root golden command, browser test,
  and advisory mutation runs for identity reconciliation and candidate grouping.
- [ ] Run the bounded three-part Opus review and GPT invariant review. Apply one controlled
  fix commit such as `fix(identity): address adversarial review findings`, then one scoped
  re-review.
- [ ] Open PR A with the measured full-result identity baseline and exact selected-row
  compatibility evidence. Merge without squash only after required checks are green.
- [ ] Fetch and verify the merge before creating PR B:

```bash
git fetch origin
git merge-base --is-ancestor <pr-a-merge-sha> origin/main
cd "$(git rev-parse --show-toplevel)/../.."
git worktree remove .worktrees/issue-270-molecular-identity
git worktree add .worktrees/issue-267-evidence-policy -b fix/issue-267-evidence-policy origin/main
```

## PR B — issue 267 governed adVNTR artifact evidence

**Branch:** `fix/issue-267-evidence-policy`  
**PR title:** `fix(advntr): govern recurrent-state evidence`  
**Issue footer:** `Refs #267`

### Task B1: Canonical governed evidence artifact

**Production mutation caught:** one of 24 states is duplicated/unreachable, unknown cohort
facts become invented numbers, JSON and digest drift, or noncanonical bytes are accepted.

**Files:** create `vntyper/scripts/canonical_json.py`,
`vntyper/modules/advntr/advntr_artifact_evidence.json`, its `.sha256`,
`vntyper/modules/advntr/artifact_evidence.py`, and
`tests/unit/test_advntr_artifact_evidence.py`; modify
  `tests/unit/test_advntr_polymorphic_calls.py`, `pyproject.toml`, and remove the carried-forward state/provenance
  list from `vntyper/modules/advntr/advntr_calibration.json` so the new artifact is the sole
  production source.

- [ ] Write RED tests for strict duplicate-key rejection, RFC 8785 bytes plus one newline,
  matching SHA-256, exactly 24 distinct active states, one confirmed/23 pending, and JSON
  `null` for count, denominator, frequency, and exact revision. Pin fields
  `cohort_name: "renome"`, `advntr_version_upper_bound_exclusive: "2.0.4"`, and
  `exact_advntr_version: null`. Define reachable as exactly matching a State string through
  the production flagging rule inventory, and assert both sets are equal. Include one
  hard-coded RFC 8785 known-answer byte vector and literal digest rather than recomputing
  expected bytes with the implementation under test.
- [ ] Store an independent literal 24-State tuple in the test module. Run each literal
  through the real flagging entry point and require exactly one hit; run non-listed
  one-edit near-matches and require zero hits. Do not generate this oracle from the
  production evidence file or rule inventory. Before B2, run the golden collision audit and
  record the exact nonzero flagged-row and affected-decision counts in the test fixture.
- [ ] Run both evidence tests and expect missing artifact/loader failures.
- [ ] Add `rfc8785>=0.1.4,<0.2`, strict decoder functions
  `load_strict_json_object()`, `canonical_json_bytes()`, and `canonical_sha256()`. Use
  `object_pairs_hook` to reject duplicate keys and `parse_constant` to reject nonfinite
  values before RFC 8785 processing. Implement
  frozen evidence entry/container types and verified packaged/file loaders. Store only the
  approved renome statement and null unknowns.
- [ ] Run evidence and package-data tests; expect green and a digest matching canonical
  file bytes.
- [ ] Commit `feat(advntr): verify governed artifact evidence`.

### Task B2: Make recurrent-state evidence inadmissible for identity agreement

**Production mutation caught:** only pending rows are filtered, a flagged row promotes
Tier A, Kestrel BAM is counted independently, or the flagged detection disappears.

**Files:** modify `vntyper/modules/advntr/artifact_evidence.py`,
`vntyper/scripts/identity_reconciliation.py`, `nomenclature_annotate.py`, `cross_match.py`,
`tests/unit/test_nomenclature_reconcile.py`, `test_nomenclature_advntr.py`, and
`test_cross_match.py`.

- [ ] Write RED cases for agreeing Kestrel plus flagged adVNTR, flagged-only adVNTR, and an
  otherwise identical unflagged row:

```python
def test_flagged_state_is_visible_but_cannot_back_tier_a() -> None:
    result = reconcile_caller_outputs(KESTREL_AND_FLAGGED_ADVNTR, EVIDENCE)
    assert result.advntr_result == "Positive (Flagged)"
    assert result.molecular_agreement is False
    assert (result.name, result.tier) == kestrel_only_name_and_tier(KESTREL_AND_FLAGGED_ADVNTR)
    assert result.evidence_disposition == "identity-insufficient"
```

- [ ] Run focused tests and observe current molecular agreement/cross-match RED.
- [ ] Return `IDENTITY_INSUFFICIENT` for all listed states regardless of confirmed/pending
  provenance. Thread disposition per observation. Keep row, detection, flag, caller-local
  name, support, and Kestrel finding; block only molecular backing and Tier-A promotion.
  Cross-match may show the comparison but cannot emit molecular agreement. Append
  `Evidence_Disposition` to positive adVNTR rows and decision explanations; do not widen
  negative output.
- [ ] Run focused tests and require every unflagged projection byte-identical.
- [ ] Commit `fix(nomenclature): exclude ambiguous advntr states from agreement`.

### Task B3: Snapshot evidence provenance on every public surface

**Production mutation caught:** a legacy run is labeled with the current package digest,
the digest is truncated, standalone reporting ignores the run snapshot, or one cohort row
loses provenance.

**Files:** modify `vntyper/scripts/summary.py`, `pipeline.py`, `artifact_names.py`,
`generate_report.py`, `cohort_inputs.py`, `cohort_summary.py`, both templates, and their
unit tests.

- [ ] Write RED tests for a canonical
  `provenance/advntr_artifact_evidence.json` snapshot, additive schema-2
  `advntr_evidence_digest`, legacy text `artifact-evidence revision not recorded`, and the
  exact approved assertion that a carried-forward recurrent State is insufficient for
  molecular identity.
- [ ] Run summary/artifact/report/cohort tests and expect missing-key and missing-snapshot
  failures.
- [ ] Resolve evidence before caller stages, atomically snapshot canonical bytes, record
  the full digest, and make report/cohort readers use only recorded run provenance. Preserve
  all existing labels and config-owned wording.
- [ ] Run unit and browser surface tests; expect green.
- [ ] Commit `feat(report): record advntr evidence provenance`.

### Task B4: Governed-policy golden baseline, docs, review, and merge

**Production mutation caught:** all flagged rows are demoted or hidden, an unflagged
decision changes, or #266 comments become findings.

- [ ] Extend the independent golden tests so only agreement/tier fields on
  `Polymorphic_Call` rows may change from A. Assert detection, caller-local names, every
  unflagged row, #266 behavior, complete by-tier metrics, roots, counts, and zero skips.
  Translate independently established pathogenic truth and report every collision with the
  literal 24-State oracle without assigning benignity or inventing cohort frequency.
- [ ] Use the exact counts frozen by B1. After B2, require the same oracle to turn green; do
  not create expectations from the changed production projection.
- [ ] Update `docs/pipeline/flagging.md` and `docs/pipeline/nomenclature.md` with the exact
  limited assertion and explicit non-claims. Commit `docs(advntr): define recurrent-state evidence`.
- [ ] Run all PR B gates, browser tests, bounded reviews, one fix wave/re-review, open the
  PR, and merge after green checks.
- [ ] Fetch the merge and create a fresh worktree:

```bash
git fetch origin
git merge-base --is-ancestor <pr-b-merge-sha> origin/main
git worktree remove .worktrees/issue-267-evidence-policy
git worktree add .worktrees/issue-175-decision-profiles -b feat/issue-175-decision-profiles origin/main
```

## PR C — issue 175 complete decision profiles

**Branch:** `feat/issue-175-decision-profiles`  
**PR title:** `feat(config): add complete decision profiles`  
**Issue footer:** `Refs #175`

### PR C file map

Create `vntyper/profiles/decision_profile.json`, its `.sha256`,
`decision_projection.json`, and focused modules `decision_profile.py`,
`decision_profile_schema.py`, `run_configuration.py`, and `profile_provenance.py`. Create
unit tests for schema, hashing, resolution, import boundaries, provenance, and round trips.
Modify CLI/pipeline consumers, summary/report/cohort surfaces, package data, docs, and
existing behavioral tests. Keep SHARK's decision component `{}` because its current JSON
contains excluded reference paths only.

### Task C1: Exhaustive field inventory and canonical packaged profile

**Production mutation caught:** any decision leaf, unit, comparator, inclusivity, or value
changes; an excluded path enters the profile; object-key order changes the hash; or an
ordered comparator list is treated as unordered.

- [ ] Write RED schema/hash tests. Represent the inventory as an object keyed by JSON
  Pointer:

```json
{
  "/components/kestrel/confidence/depth_score_low": {
    "class": "fixed-safety",
    "value": 0.00469,
    "unit": "depth-score-ratio",
    "comparator": "gte",
    "inclusive": true
  }
}
```

  Assert object-key reordering preserves hash, while reversing a semantic array changes
  it. Assert every current Kestrel, adVNTR, nomenclature, flagging, motif, selection, and
  cross-match decision leaf belongs to exactly one class. Pin independent literals for
  `0.00469`, `0.00515`, `20`, `21`, `100`, `200`, flank `8`, thin record support `3`, both
  separately named/unit-bearing Tier-A values `5`, the low/GG equality link, and every
  comparator/inclusivity value. A profile key for an issue-268 model/version guard must
  fail as unknown.
- [ ] Run new tests and expect missing packaged files/modules.
- [ ] Implement `ValidationClass`, `DecisionField`,
  `validate_complete_inventory()`, component projections, and immutable-field checks. Add
  kinds `packaged`, `explicit-custom`, `generated`; sources `package`, `explicit-cli`.
  Make fixed safety include all listed numeric values, evidence units, comparators,
  rerun-required fields, #266 floor, motif/caller maps, and governed evidence. Keep #268
  guards excluded. Define the closed generated metadata schema (packaged base hash,
  generator name/version, objective, dataset and partition hashes, seed) and neutral
  generated-mutable dominance/whole-locus-abstention fields; validate generated profiles
  now even though only PR D constructs them. Generated profiles must copy every
  explicit-custom leaf equal to packaged.
- [ ] Generate canonical profile bytes, checked-in digest, and exact legacy decision
  projection. Run tests green.
- [ ] Commit `feat(config): add complete decision profile`.

### Task C2: Strict one-profile resolver and CLI entry

**Production mutation caught:** a partial profile overlays defaults, invalid explicit input
falls back, duplicate keys or nonfinite numbers pass, or stage artifacts exist before
resolution fails.

**Files:** create/modify `decision_profile.py`, `run_configuration.py`,
`tests/unit/test_decision_profile_resolution.py`, `cli_parser.py`, `cli.py`,
`cli_handlers.py`, and CLI tests.

- [ ] Write RED cases for absent file, directory, unreadable bytes, duplicate key, NaN,
  infinity, missing/unknown key, unsupported version, invalid enum/range, comparator-rule
  error, rule/evidence mismatch, incomplete component, immutable-field mismatch, and output
  directory non-creation. Add only pipeline syntax:

```text
vntyper pipeline ... --decision-profile /path/to/complete-profile.json
```

- [ ] Run resolution and CLI tests; expect unknown option and missing resolver.
- [ ] Implement `parse_decision_profile()`, `resolve_decision_profile(path)`, frozen
  `RunConfiguration`, and `resolve_run_configuration()`. Resolve once in `cli.main()` before
  log/output creation only when the parsed subcommand is `pipeline`; report/cohort never
  resolve a current package profile. There are no overlays, includes, environment/current-directory discovery, or
  fallback. Argparse usage remains exit 2; validated application failure remains exit 1.
- [ ] Run CLI/resolution tests green and a profile canonical round trip byte-identical.
- [ ] Commit `feat(cli): resolve one explicit decision profile`.

### Task C3: Snapshot profile and introduce summary schema 3

**Production mutation caught:** the hash covers parsed input but snapshot bytes differ,
absolute operator paths leak, PR-B evidence digest is omitted, or schema 3 accepts a
mismatched profile.

- [ ] Write RED tests for exact top-level schema-3 fields from the approved design and the
  relative snapshot `provenance/decision_profile.json`. Schema 1/2 must render `decision
  profile not recorded by legacy run`; schema 3 must fail on a missing field, bad identity
  status/count combination, unknown identity version, or digest mismatch.
- [ ] Run provenance/summary/artifact tests and observe missing schema-3 support.
- [ ] Implement `DecisionProfileProvenance`, `profile_summary_fields()`,
  `snapshot_decision_profile()`, `verify_profile_snapshot()`, and
  `resolve_summary_profile()`. Raise `SUMMARY_SCHEMA_VERSION` once to 3. Persist source as
  `package` or `explicit-cli`, never a path, and require PR-B artifact digest.
- [ ] Run focused tests green.
- [ ] Commit `feat(summary): record decision profile provenance`.

### Task C4: Remove hidden Kestrel/adVNTR/SHARK configuration ownership

**Production mutation caught:** a custom run changes a field but a module-level packaged
global still decides behavior, or a compatibility wrapper silently mixes contexts.

**Files:** modify `kestrel_genotyping.py`, `pipeline_kestrel.py`,
`modules/advntr/advntr_genotyping.py`, `modules/shark/shark_filtering.py`, `pipeline.py`,
and their existing tests; create `tests/unit/test_decision_profile_import_boundaries.py`.

- [ ] Write an AST RED test forbidding decision code from reading decision leaves through
  module-level `kestrel_config`, `advntr_config`, `advntr_settings`, `shark_config`, or
  `shark_settings`; excluded runtime paths/references remain available through the frozen
  runtime component. Add a source-scan assertion that SHARK currently has no independent
  decision constant. Build a
  differential explicit-custom profile changing every mutable leaf and assert each exact
  stage argument changes.
- [ ] Run AST and stage tests; expect current import globals to fail the contract.
- [ ] Pass immutable resolved components through pipeline calls. Compatibility wrappers
  receive `resolved_component: Component | None` plus `custom_context_active: bool`; they
  load packaged components only when both indicate no run context and raise when
  `custom_context_active` is true without the explicit component. Extract pure projection into new modules instead of enlarging
  oversized stage files.
- [ ] Run Kestrel confidence/filtering, adVNTR command/output, SHARK pipeline, and import
  boundary tests green.
- [ ] Commit three focused commits:

```bash
git commit -m "refactor(kestrel): consume resolved decision config"
git commit -m "refactor(advntr): consume resolved decision config"
git commit -m "refactor(shark): remove sidecar config globals"
```

### Task C5: Pass resolved policy through identity and cross-match

**Production mutation caught:** nomenclature reads old globals, one universal support `5`
loses source units, BAM becomes independent, or governed adVNTR disposition is bypassed.

- [ ] Write RED differential tests in nomenclature, BAM, reconciliation, and cross-match
  families. Explicitly assert the two value-5 fields have different units and names.
- [ ] Run focused tests and observe old globals/implicit values.
- [ ] Project and pass resolved nomenclature/cross-match components. Preserve
  `legacy-selection-v1`, identity semantics, evidence disposition, and every fixed value.
- [ ] Run focused tests green and commit
  `refactor(nomenclature): consume resolved decision policy`.

### Task C6: Verify report and cohort provenance

**Production mutation caught:** standalone report substitutes the current package profile,
cohort uses the first sample's profile for all rows, mixed-profile metrics are pooled, or
a tampered schema-3 snapshot renders.

- [ ] Write RED report/cohort tests for package and explicit-cli provenance, tamper
  rejection, legacy text, profile hash exports, visible grouping, and suppression of pooled
  decision-performance aggregates across hashes. Include schema 1, schema 2 without the
  PR-B evidence digest, schema 2 with its recorded digest, and schema 3; no legacy reader may
  infer a current package digest.
- [ ] Run report/cohort/CLI report tests and observe missing profile verification.
- [ ] Make standalone report verify only the recorded snapshot. Pass immutable resolved
  report/cohort components explicitly; extend cohort sample payload
  with profile provenance, group by SHA-256, export `Decision_Profile_ID`,
  `Decision_Profile_Revision`, and `Decision_Profile_SHA256`, and suppress mixed pooled
  metrics.
- [ ] Run unit and browser tests green. Commit:

```bash
git commit -m "feat(report): verify recorded decision profiles"
git commit -m "feat(cohort): group samples by decision profile"
```

### Task C7: Packaged projection, docs, review, and merge

**Production mutation caught:** no-profile output differs from PR B, a generated profile
auto-activates, fixed `3`, either `5`, flank `8`, or any packaged threshold becomes mutable.

- [ ] Extend golden tests to assert the complete packaged manifest/hash and byte/exact-field
  PR-B decision projection. Add one mutation case for every fixed value, comparator, unit,
  #266 floor, plus an unknown-key rejection for every excluded #268 guard. Verify the
  packaged `.sha256` sidecar is 64 lowercase hexadecimal characters plus one newline and
  hashes the exact RFC-8785-plus-newline bytes.
- [ ] Run explicit-root golden tests RED until profile provenance is expected, then require
  the governed PR-B decision projection exactly and zero skips.
- [ ] Update pipeline/report/cohort CLI docs, configuration guide, and nomenclature docs.
  State that explicit profiles are complete files, generated outputs never activate, and
  package defaults are unchanged.
- [ ] Run all gates, browser, `make check-full` if available, bounded reviews, one fix
  wave/re-review, open the PR, and merge after green checks.
- [ ] Fetch the merge and create the fresh PR-D worktree:

```bash
git fetch origin
git merge-base --is-ancestor <pr-c-merge-sha> origin/main
git worktree remove .worktrees/issue-175-decision-profiles
git worktree add .worktrees/issue-269-calibration -b feat/issue-269-calibration origin/main
```

## PR D — future issue 269 dominance/abstention evaluation

**Branch:** `feat/issue-269-calibration`  
**PR title:** `feat(calibration): add opt-in dominance evaluation`  
**Issue footers:** `Refs #269` and `Refs #295`

PR D ships an evaluation engine, not a selected model. `--partitions PARTITIONS` is the
canonical study declaration passed to `extract`; it contains the complete frozen finite candidate grid, maximum free-parameter count,
minimum stratum counts, abstention cap, objective, group-bootstrap seed, and split-role
manifests. The repository supplies no real-study grid and no active generated profile.
Tests use a checked-in synthetic declaration with exact finite values solely to prove
engine behavior. The current simulations may exercise the workflow but are marked
`development` and rejected as closure evidence.

### Task D1: Dominance evidence and closed abstention decisions

**Production mutation caught:** votes count edits instead of records, shares use only
supporting records, missing XD becomes zero, XD chooses a winner, a tie falls through, or
inadmissible adVNTR promotes an identity.

**Files:** reuse `nomenclature_bam_evidence.py`; create
`vntyper/scripts/nomenclature_dominance.py` and
`tests/unit/test_nomenclature_dominance.py`.

- [ ] Write RED tests for tie despite unequal XD, runner-up present during XD veto,
  missing evidence, governed adVNTR disposition, and all abstention reasons. Use closed
  outcomes `selected`, `abstained`, `not-applicable`; use closed reason tokens
  `record-tie`, `insufficient-dominance`, `xd-missingness`, `xd-concentration`,
  `xd-discordance`, and `inadmissible-advntr`.
- [ ] Add RED count/share assertions: one record votes at most once per identity, may vote
  for two co-occurring identities, every eligible record enters the denominator, shares may
  sum above one, count margin is top minus runner-up, and share margin divides by that same
  denominator.
- [ ] Run evidence/dominance tests and expect missing dominance module.
- [ ] Implement `DominanceDecision` and `evaluate_dominance(evidence, component)` in this
  order: applicability; record-count tie; top count/share/margins; governed disposition;
  whole-locus XD veto. Any veto returns abstention immediately. XD never reorders
  identities, breaks a tie, or selects a runner-up.
- [ ] Run tests green and register the module for advisory mutation testing.
- [ ] Commit `feat(calibration): add abstaining dominance policy`.

### Task D2: Strict study declaration, artifacts, and transitive leakage groups

**Production mutation caught:** only adjacent group collisions are checked, a group
namespace is dropped, role/hash mismatches pass, labels enter features, or fit opens
validation/held-out roots.

**Files:** create `vntyper/scripts/calibration_contract.py`,
`vntyper/scripts/calibration_manifest.py`, `vntyper/scripts/calibration_features.py`,
`tests/unit/test_calibration_contract.py`, `tests/unit/test_calibration_manifest.py`, and
`tests/unit/test_calibration_features.py`.

- [ ] Write RED strict-decoder tests for protocol/evidence/policy/metrics/attestation
  objects. The synthetic protocol fixture contains exactly:

```json
{
  "objective": "lexicographic-safety-v1",
  "bootstrap_iterations": 10000,
  "bootstrap_interval": "percentile",
  "multiplicity_method": "holm",
  "seed": 295,
  "maximum_free_parameters": 4,
  "minimum_stratum_count": 2,
  "maximum_abstention_fraction": 0.25,
  "candidate_grid": {
    "minimum_record_count_margin": [1, 2],
    "minimum_record_share": [0.5, 0.75],
    "minimum_record_share_margin": [0.0, 0.25],
    "xd_veto": ["disabled", "missingness"]
  }
}
```

  This is test data, not a shipped recommendation or real-study protocol.
- [ ] Add RED transitive-leakage cases across namespaced individual/family, simulated pair,
  backbone/seed lineage, replicate/rerun, depth-series source, batch, and repeat context.
  Add allowlist failures for truth, path, sample, batch identity, raw tuple, selected name,
  and post-decision tier. Add positive allowlist cases for assay class and the four XD
  summaries; reject batch identity, raw XD weighting, and any unlisted XD aggregation.
  Require distinct `training` and `policy-selection` partitions, reject transforms/caps
  fitted from the latter, reject fit access to validation/held-out roles, and reject a
  `development` provenance artifact relabeled as held-out.
- [ ] Run new tests and expect missing modules.
- [ ] Implement strict immutable decoders, canonical hashes, connected-component leakage
  validation, role access restrictions, separate keyed feature/label artifacts, and the
  approved runtime feature allowlist. Snapshot/hash the complete PARTITIONS declaration
  into EVIDENCE with allowlisted features, separately keyed labels, and a non-feature
  shipped decision projection. `extract` must fail on old runs missing complete PR-A
  replay artifacts and must not reopen BAM or rerun callers.
- [ ] Run tests green and commit `feat(calibration): validate replay evidence and splits`.

### Task D3: Lexicographic objective and statistics

**Production mutation caught:** any objective component is deleted/reversed, F1 becomes a
default, abstentions leave denominators, rows rather than groups are bootstrapped, or a
more sensitive policy with one wrong Tier-A identity wins.

**Files:** create `vntyper/scripts/calibration_objective.py`,
`vntyper/scripts/calibration_statistics.py`, `tests/unit/test_calibration_objective.py`,
`tests/unit/test_calibration_statistics.py`, and add both pure modules to
`scripts/mutation_test.py`.

- [ ] Write two-candidate RED fixtures that independently flip the winner for each tuple
  component:

```python
def lexicographic_safety_key(metrics: CandidateMetrics) -> tuple[object, ...]:
    return (
        metrics.wrong_tier_a_displayed_names,
        metrics.control_findings,
        metrics.wrong_displayed_names_all_tiers,
        -metrics.macro_exact_recovery,
        -metrics.binary_detection_sensitivity,
        metrics.free_parameter_count,
        metrics.profile_sha256,
    )
```

  Define displayed-name errors by comparison with the independent expected rendered name,
  separately from molecular-identity exactness. Make any nonzero
  `wrong_tier_a_displayed_names` or `control_findings` a hard inadmissibility result; add
  cases where all candidates violate each constraint and fitting emits no candidate. Define
  free-parameter count as the number of non-neutral generated rule leaves: zero
  margins, disabled XD veto, and false governed-disposition veto are neutral. Assert all-abstain and Tier-A-unreachable candidates are inadmissible, one wrong Tier A
  always loses, one control finding loses despite fewer wrong names, and selective
  abstention remains in exact-recovery and sensitivity denominators. Add distinct RED cases
  for minimum stratum count, abstention cap, applicability mismatch, a free-parameter-only
  tie, and macro-versus-micro selection on imbalanced assay-by-mutation-class strata. For
  every adjacent tuple boundary, oppose a one-unit degradation in the earlier component
  with arbitrarily large improvements in every later component and require the earlier
  component to decide.
- [ ] Add RED statistics tests distinguishing complete-group bootstrap from row bootstrap,
  10,000 deterministic resamples, one-sided paired lower-bound equality at zero, two-sided
  intervals, Clopper-Pearson zero-event intervals, deterministic ROC/PR, and joint surfaces.
  Assert the bootstrap resampling-unit count equals the transitive connected-component
  count, not the number of raw namespace keys.
- [ ] Run tests and expect missing modules.
- [ ] Implement the tuple literally. Use every mutated truth member as the exact-recovery
  and sensitivity denominator, with controls counted separately; binary detection may count
  a wrong identity as detected while exact recovery cannot. The dominance layer may select
  or abstain but never assign/demote a tier. Represent objective rates as exact
  `fractions.Fraction` values. Pair every candidate against the immutable shipped decision
  projection stored in EVIDENCE, and apply the protocol's Holm correction across the frozen
  candidate family before lexicographic selection.
  Require the objective argument. Make both detection and macro-exact paired lower bounds
  at least the zero-percentage-point non-inferiority margin.
- [ ] Run tests and targeted mutation commands; commit
  `feat(calibration): implement safety-first candidate selection`.

### Task D4: Generate only allowable explicit profiles

**Production mutation caught:** generation changes fixed `3`, either `5`, flank `8`, a
Kestrel cutoff/unit/comparator, #266/#268, explicit-custom fields, package bytes, or active
run context.

**Files:** create `vntyper/scripts/calibration_profiles.py`,
`vntyper/scripts/calibration_workflow.py`, `tests/unit/test_calibration_profiles.py`, and
`tests/unit/test_calibration_workflow.py`;
reuse PR-C's sole profile loader/canonicalizer.

- [ ] Write RED mutation cases for every immutable field and prove the only changing leaves
  are protocol-declared dominance/whole-locus-abstention fields. Assert generated output is
  outside `vntyper/profiles`, reloads through `resolve_decision_profile()`, records base
  hash/generator/objective/dataset/partition/seed, and remains inactive without explicit
  `pipeline --decision-profile PATH`. Hash every packaged profile file before and after
  generation. Reject generated predicates that threshold or otherwise proxy any
  fixed-safety feature, and pin excluded #268 keys as unknown-key errors versus #266/fixed
  changes as immutable-field errors.
- [ ] Run tests and expect missing generation/workflow modules.
- [ ] Implement `build_generated_profile()`, `validate_generated_allowlist()`,
  `extract_evidence()`, `fit_candidate()`, and `validate_candidate()`. Defer
  `evaluate_locked_candidate()` to Task D5 so custody cannot be bypassed. Before fitting a sweep, replay the shipped profile and
  require aggregate and per-tier metrics plus row order/name/confidence/flag/tier/support/
  tie/abstention state to equal PR C.
- [ ] Run tests green; commit `feat(calibration): emit consumable opt-in profiles`.

### Task D5: Custody guards and one-use evaluation

**Production mutation caught:** evidence opens before precommit, validation swaps profiles,
failed candidates remain active, held-out evidence is evaluated twice, or local metadata is
misrepresented as proof of independent custody.

**Files:** create `vntyper/scripts/calibration_custody.py` and
`tests/unit/test_calibration_custody.py`; modify
`vntyper/scripts/calibration_workflow.py` and its test.

- [ ] Write RED tests with an opener that raises unless precommit was durably written;
  assert second-use refusal including after failure, candidate/profile/protocol hash
  binding, validation-set retirement, and held-out-role-only evaluation. Add a two-profile,
  one-evidence-hash case that refuses the second attempt, and patch `open_locked_payload()`
  to raise through `evaluate_locked_candidate()`.
- [ ] Run tests and expect missing custody module.
- [ ] Implement `write_precommit()`, `open_locked_payload()`, `record_consumption()`, and
  `retire_candidate()` plus custody-wired `evaluate_locked_candidate()` with atomic
  append-only local receipts keyed by evidence hash. Held-out extraction occurs inside the
  custodian operation after precommit; development-side extraction cannot mark evidence
  closure-eligible. Label these safeguards as
  local guards; independent custody requires an external named custodian and durable
  authority outside this repository.
- [ ] Run tests green; commit `feat(calibration): guard locked evaluation custody`.

### Task D6: Four-command CLI and atomic artifacts

**Production mutation caught:** objective defaults to F1, argparse/application exit codes
blur, failed commands leave partial directories, or fit can access another role.

**Files:** create `vntyper/scripts/cli_calibrate.py`; modify `cli_parser.py`, `cli.py`, and
CLI tests.

- [ ] Write RED parser/dispatch tests for exactly:

```text
vntyper calibrate extract --truth TRUTH --partitions PARTITIONS --runs RUNS --output EVIDENCE
vntyper calibrate fit --evidence EVIDENCE --objective lexicographic-safety-v1 --output CANDIDATE
vntyper calibrate validate --profile PROFILE --evidence VALIDATION --output VALIDATION_ATTESTATION
vntyper calibrate evaluate --profile PROFILE --evidence LOCKED_HELDOUT --output HELDOUT_ATTESTATION
```

  Assert missing or unknown `--objective` is argparse exit 2, an objective differing from
  the one hashed inside PARTITIONS/EVIDENCE is validated exit 1, and no partial output exists.
- [ ] Run parser/dispatch tests and expect unknown `calibrate` command.
- [ ] Implement `add_calibrate_subparser()`, `handle_calibrate()`, and four private handlers.
  Atomically install output directories only after checksums and all artifacts succeed.
- [ ] Run CLI tests green; commit `feat(cli): add calibration workflow commands`.

### Task D7: Deterministic calibration reports

**Production mutation caught:** lower tiers disappear, abstention reasons are hidden, hashes
or limitations are omitted, curves depend on JavaScript/network, or XD is labeled as reads.

**Files:** create `vntyper/scripts/calibration_report.py`, `vntyper/templates/calibration_report.html`,
`tests/unit/test_calibration_report.py`, `tests/browser/test_calibration_report.py`, and
artifact-writing tests.

- [ ] Write RED tests for atomic output contents: evidence manifest, features, labels,
  groups, baseline, checksums; fit profile/attestation/grid/metrics/HTML; validation or
  held-out metrics, intervals, ROC/PR tables, joint-surface TSVs, abstention table, ledger
  receipt, and static HTML.
- [ ] Assert HTML contains every tier's displayed/exact/wrong counts, abstention reasons,
  profile/protocol/evidence hashes, fitted versus validation versus held-out labels,
  boundary coverage, and small-n limitations. Also assert seeds, software/reference
  versions, sample composition, assay, depth, read length, independently measured
  array-size data, mutation classes, objective, manifest hashes, access attempts, two-sided
  intervals, one-sided paired bounds, and abstention rate/reasons by split. Run it offline
  and with JavaScript disabled. Require the literal local-custody limitation and
  `Reporting an interval is not a clinical safety claim.`
- [ ] Run unit/browser tests and expect missing report/artifacts.
- [ ] Implement static rendering and exact terminology `optional minimum k-mer depth` for
  XD. Keep this separate from the sample report.
- [ ] Run unit and browser tests green; commit
  `feat(report): render calibration evidence and limitations`.

### Task D8: Golden replay, truthful blocker, docs, review, and merge

**Production mutation caught:** the sweep runs before shipped reproduction, one external
root is silently absent, simulations are called held-out, skips count as evidence, or #295
closes without custody-quality data.

- [ ] Create `tests/golden/calibration_oracle.py` with an AST guard forbidding production
  predicates, canonicalization, serialization, grouping, and decision imports, and
  `tests/golden/test_calibration_golden.py`. Assert both roots loaded, 200 mutated/200
  controls, complete row/locus counts, every tier, and exact PR-C governed projection before
  fitting. Reject the current corpus for `evaluate` closure eligibility because it is
  previously examined development evidence.
- [ ] Run explicit-root golden RED, then green with zero skips. Record fitted development
  results separately from validation and held-out fields; do not call them validated.
- [ ] Add `docs/cli/calibrate.md`, `docs/development/calibration-validation.md`, and update
  configuration/changelog/navigation. Every page states explicit opt-in, no shipped cutoff
  or default change, simulations are not external validation, #295 remains blocked, and
  intervals are not clinical safety claims.
- [ ] Commit golden and documentation changes:

```bash
git add tests/golden/calibration_oracle.py tests/golden/test_calibration_golden.py \
  docs/cli/calibrate.md docs/development/calibration-validation.md \
  docs/user-guide/configuration.md docs/about/changelog.md mkdocs.yml
git commit -m "test(calibration): pin development replay and blocker"
```
- [ ] Run all common gates, browser, explicit-root golden, available full integration, and
  advisory mutation for objective/dominance. Run bounded Opus partitions and GPT review,
  resolve one fix wave, and perform at most one scoped re-review.
- [ ] Open PR D with `Refs #269` and `Refs #295`, exact packaged profile hash/projection,
  and the external-evidence blocker. Merge under standing authorization only after all
  green checks and all Critical/Important findings are resolved.

## Post-merge issue state and external validation gate

- [ ] Fetch `origin/main`, verify all four merge SHAs are ancestors, and inspect current
  main CI. Post a concise #295 update comparing Phase 1, PR-A identity-policy, PR-B governed
  policy, PR-C packaged replay, and PR-D development-only results.
- [ ] Keep #295 open and record the exact unavailable evidence: independent real carriers
  with orthogonally established exact alleles, representative independent negatives,
  multiple laboratories/assays/depths, reanalysed renome state counts, and a named
  custodian-held locked cohort.
- [ ] Do not draft a closure claim or release change. If qualifying evidence later becomes
  available, stop before opening it. Obtain approval for a precommitted protocol containing
  the exact finite grid, free-parameter limit, group definitions, minimum per-stratum
  counts, abstention cap, multiplicity method, precision target, custodian, and durable
  precommit/ledger authority.
- [ ] After that approval, evaluate one frozen profile once. Closure requires zero observed
  wrong Tier-A identities, zero observed control findings, no increase in wrong displayed
  identities, and both zero-margin paired non-inferiority bounds. Report all confidence
  intervals and limitations.
- [ ] Only if that evidence passes, prepare a separate release PR from then-current main and
  determine the version then. Stop for explicit authorization before merging the release
  PR, creating/pushing the immutable tag, or sending authenticated
  `repository_dispatch(type=vntyper_release)`. Verify main ancestry, all ten exact required
  checks on the full tagged SHA, and promotion by verified image digest.

## Plan self-check

- [ ] Map every approved design heading to at least one task: identity A1–A8; artifact
  policy B1–B4; profile resolution C1–C7; future dominance/calibration D1–D8; external
  closure gate above.
- [ ] Search this plan for incompleteness markers and vague implementation verbs; replace
  every match with an exact file, symbol, command, test assertion, or documented blocker.
- [ ] Verify type and field consistency: `MolecularIdentity`, `EvidenceDisposition`,
  `IdentityCandidateSet`, `ResolvedDecisionProfile`, `DominanceDecision`, the identity
  quartet, profile kind/source enums, and abstention tokens must have one spelling across
  all tasks.
- [ ] Verify no task changes a shipped cutoff, calls XD reads/molecules, sums equivalent
  support, treats BAM as another caller, auto-selects a profile, counts a skip as evidence,
  closes #295, or begins release work.
