# VNtyper issue 295 completion program design

**Date:** 2026-08-31  
**Scope:** issues 270, 267, 175, 269, and 295  
**Baseline:** `530599ab479b83e45807c4305cb5d5cf3763f7a8`  
**Status:** proposed for owner approval

## Purpose

Phase 1 made Kestrel `output.bam` evidence truthful without changing a decision. This
program completes the remaining architecture needed to distinguish molecular identity
from caller representation, govern artifact evidence, consume explicit future calibration
profiles, evaluate dominance and abstention, and state honestly when external validation
is insufficient.

The program is delivered as four sequential functional pull requests. One reviewed design
and one master implementation plan govern all four; each pull request still has its own
test-first implementation, review, verification, CI, and merge boundary. Release
preparation is a fifth pull request after the functional program and external evidence are
complete.

## Non-negotiable invariants

1. **Configuration and calibration do not change shipped decision behavior.** The
   packaged profile retains:

   - Kestrel `Depth_Score` low `0.00469` and high `0.00515`;
   - alternate k-mer-path-depth boundaries `20`, `21`, and `100`;
   - active-region boundary `200`;
   - linked GG floor `0.00469`;
   - Tier-A minimum evidence value `5`, interpreted separately in each source's unit;
   - thin BAM rescue boundary `3` resolved haplotype records; and
   - BAM flank `8`.

   No implementation, generated artifact, migration, or validation result edits these
   values in the packaged profile. The migration inventory covers every decision field
   from the existing Kestrel, adVNTR, SHARK, nomenclature, flagging, and cross-match
   configurations, including its unit, comparator, inclusivity, and value. A checked-in
   canonical profile hash and golden decision projection prove that the entire
   inventory—not only the values listed above—matches pre-migration behavior. A future
   custom profile can be selected only by an explicit CLI path; generating it never
   installs or activates it. `legacy-selection-v1` remains the packaged default in every
   pull request in this program. Identity-based reconciliation in work package A and the
   governed artifact rule in work package B are the two explicit default semantic
   corrections; neither changes selection, voting, or a numeric cutoff.

2. Resolved Kestrel BAM records remain resolved haplotype records. Their count remains
   haplotype-record support. `XD` remains optional typed minimum k-mer depth. Neither is a
   sequencing-read or molecule count.

3. `XD` never expands votes, selects an edit, breaks a record-count tie, or substitutes
   for haplotype-record support. A future explicit custom policy may use only its declared
   typed summary as a whole-locus abstention veto. A veto emits abstention immediately; it
   never falls through to a runner-up.

4. The six Kestrel row gates continue to mark before the final filter. Evidence rows are
   not dropped early.

5. Issue 266 remains intact: the reporting floor is outside the custom-profile and sweep
   surface and does not move; a subthreshold candidate remains one `##` comment, never a
   result row or finding, and never enters cross-match or nomenclature agreement.

6. Issue 268 remains a coverage/model compatibility boundary. Calibration cannot waive
   its fail-closed adVNTR model/version checks.

7. Kestrel VCF and Kestrel BAM remain one caller. Their agreement is not independent
   corroboration.

8. Existing config-owned clinical-sounding messages are not invented or rewritten by
   this program.

## Selected architecture

Three approaches were considered.

### Selected: typed observations, identity, evidence disposition, and policy

Introduce focused pure models between caller parsing and presentation. Raw representations
remain recoverable, canonical identity becomes a first-class key, and decision policy is
explicit. This separates concepts without rewriting the callers or hiding compatibility
changes.

### Rejected: extend `Nomenclature`

`Nomenclature` already combines a name, event word, tier, flags, support consequences, and
source in a 1,200-line module. Adding identity equality and configuration ownership there
would make display changes capable of altering molecular agreement.

### Rejected: improve the global sort

The issue 270 benchmark provides no support for another ordering heuristic: most wrong
selections are not ranking mistakes that the tested universal sorts repair. Raw `POS`,
largest edit, `Depth_Score`, and support reorderings trade mutation classes against each
other. A new sort would also be an uncalibrated decision rule.

## Concept model

The following concepts are distinct and are never aliases:

| Concept | Meaning | Existing or new representation |
| --- | --- | --- |
| Detection | A caller observed and gate-annotated a row; the gate outcome is separate | Existing caller row and six Kestrel booleans |
| Frame consequence | Net length change and mechanical frame status | Existing scoring values, projected into `FrameConsequence` |
| Representation | The caller's coordinates and alleles in its own reference context | `KestrelRepresentation` or `AdvntrRepresentation` |
| Molecular identity | A normalized coding edit on the canonical MUC1 X repeat | `MolecularIdentity` |
| Display name | Human-facing shorthand such as `59dupC` | Projection from identity, never an equality key |
| Evidence disposition | Whether one observation may support molecular agreement | `EvidenceDisposition` |
| Decision | Selected identity or an explicit abstention reason | `IdentityDecision` |

The pure identity module defines frozen typed values equivalent to:

```python
@dataclass(frozen=True)
class CodingEdit:
    start: int
    end: int
    deleted: str
    inserted: str


@dataclass(frozen=True)
class MolecularIdentity:
    reference: Literal["MUC1-X-60-coding-v1"]
    edits: tuple[CodingEdit, ...]


@dataclass(frozen=True)
class IdentityTranslation:
    identity: MolecularIdentity | None
    status: Literal["resolved", "unresolved"]
    failure: Literal[
        "invalid-allele",
        "missing-motif-context",
        "motif-anchor-mismatch",
        "pair-boundary-edit",
        "non-x-unit",
        "reconstruction-mismatch",
    ] | None
    context_diverges: bool
```

Insertion coordinates use `end == start - 1`. Edits are non-overlapping, normalized, and
sorted by coordinate; an empty edit tuple is invalid. The identity key is the normalized
tuple `(reference, edits)`. Its stable serialization is one edit followed by any
co-occurring edits separated with semicolons:

```text
MUC1-X-60-coding-v1|start|end|deleted-or--|inserted-or--[;start|end|deleted-or--|inserted-or--]
```

Alleles are uppercase `[ACGT]*`; edited bases and insertion starts are within canonical
positions 1–60; an indel that normalizes across a unit junction is unresolved; and direct
construction cannot bypass normalization. Identity version 1 values compare only with
version 1 values. A future reference version requires an explicit migration and never
compares equal by display name. Canonical `59dupC`, displayed in the existing MUC1
X-repeat coordinate frame, therefore has identity
`MUC1-X-60-coding-v1|60|59|-|C`. Display names, event words, ambiguity intervals, repeat
forms, confidence tiers, flags, raw positions, and support values do not participate in
identity equality.

## Work package A: molecular identity

### Translation and equivalence

Every Kestrel row is translated independently using its complete motif-pair context. The
translator resolves the named 120 bp motif pair, validates the raw reference allele,
applies the edit, changes to coding orientation, and identifies the affected 60 bp repeat
unit. Pair-junction edits and edits outside an X-compatible unit are unresolved. The
source unit is aligned to canonical X using the checked-in motif map; the neighbouring
unit is the required anchor and must match that map. The alternate X unit is then compared
with canonical X, common prefix and suffix are trimmed, and an ambiguous indel is shifted
to the coding-strand 3-prime-most representation. Reapplying the resulting edit must
reconstruct the entire 60 bp alternate X unit exactly. This exact 60 bp window, plus the
neighbouring-unit anchor, is the equivalence proof. A wrong motif anchor therefore cannot
become equivalent merely because its `POS/REF/ALT`, event word, or display name resembles
another row.

Equivalent representations include the two observed canonical dupC shapes `G>GG` and
`C>CG`, homopolymer anchors that reconstruct the same alternate sequence, and motif-pair
records that produce the same canonical edit. Genuinely different edits remain distinct
even when depth, frame consequence, or rendered name is equal.

The minimum raw representation key is `(Motifs, POS, REF, ALT)`. It is used only to retain
and count a caller representation; it is never called molecular identity. Translation
occurs before the existing compatibility projection that deduplicates or groups by raw
`(POS, REF, ALT)`. That legacy projection remains decision-compatible in the packaged
policy, but it cannot erase a distinct identity hypothesis from the typed candidate set.
Context-local support is not summed across equivalent motif representations, because
those rows may arise from the same evidence aligned against alternative motif references.

### Selection boundary

Identity-aware selection receives an `IdentityCandidateSet`. It groups equivalent
representations and exposes every distinct identity hypothesis to policy. Within a group,
the existing total ordering chooses the representative and that row alone supplies
reported support and depth; support is never summed. Blocking gates, flags, and context
divergence are the union across the group, so grouping cannot unblock Tier A. The packaged
`legacy-selection-v1` policy nevertheless selects the exact row selected before this
program from the compatibility projection, preserving row cardinality, order, name,
confidence, flag, depth, support, tie behavior, and tier. Identity grouping is additive
evidence until an operator explicitly supplies a custom policy path; no program milestone
activates another default.

This boundary is deliberate. Work package A establishes the correct truth target and
records truly equivalent representations without collapsing the packaged selection or
emitted result. It does not introduce an uncalibrated dominance or blanket-abstention
heuristic. A future explicit custom policy can abstain on the same typed candidate set
without reparsing caller rows.

### Reconciliation and public surfaces

Reconciliation keys caller backing by `MolecularIdentity`, never by a name string. An
unresolved translation cannot agree merely through event class or position. Context
divergence means the source motif differs from canonical X outside the normalized edit;
it remains representation evidence and continues to block Tier A. After work package B,
Tier A requires two independent callers with the same resolved identity, a known-variant
match, source-bound
support at or above separately named `5 alternate k-mer-path depth` and `5 adVNTR
sequencing reads` fields, and none of context divergence, sequence-undetermined, caller
disagreement, or inadmissible-evidence dispositions. Kestrel VCF and BAM together still
count as one caller. In work package A every existing adVNTR row retains its current
admissibility; work package B is the later boundary that assigns `Polymorphic_Call` the
inadmissible disposition.

Positive caller TSVs append these columns while retaining existing columns and meanings:

- `Molecular_Identity`: stable identity serialization, empty when unresolved;
- `Molecular_Identity_Status`: `unique`, `legacy-selected-among-multiple`, or
  `unresolved`;
- `Equivalent_Representation_Count`: positive integer for a resolved identity and `0`
  when unresolved; and
- `Identity_Hypothesis_Count`: number of distinct resolved identities considered.

The four fields are per emitted positive row. Status is `unique` exactly when the selected
row resolves and there is one resolved hypothesis, and
`legacy-selected-among-multiple` exactly when it resolves and there is more than one.
`Molecular_Identity` is empty if and only if status is `unresolved`, and an unresolved row
has representation count `0`; its hypothesis count still reports other resolved
hypotheses, if any. Hypothesis count includes only rows that already pass the issue 266
reporting floor and the unchanged packaged detection policy. The ten-column negative row contract and its
single `##` comment are unchanged and gain no identity text. `kestrel_pre_result.tsv`
retains every raw representation and all six gates, with additive diagnostic identity
columns that never enter result or finding counts. Sample and cohort HTML/TSV/CSV/JSON
surfaces use the same four field names and carry all four values. They never infer identity
from old `POS/REF/ALT` or `Nomenclature` cells.

Work package A adds fields without changing the summary schema version. Work package C
introduces summary schema 3 once, with the complete required identity and decision-profile
field set defined together. Schema 1 and schema 2 remain readable and render the literal
`legacy identity not recorded`; schema 3 rejects missing fields, invalid status/count
combinations, unknown identity versions, or a profile snapshot/hash mismatch. Provenance
records `decision_policy: legacy-selection-v1` under the packaged default.

### Structural limitation

Kestrel VCF cannot represent the known deletion-plus-insertion truth for seven simulated
`delinsAT` carriers or the missing candidate in `pair_4092`. The identity layer never
synthesizes delins from separate VCF rows and never credits these cases as recovered. The
existing BAM path may recover a delins only where the unchanged Phase 1 eligibility,
unweighted record voting, tie, thin-support, and VCF-primary refinement rules already do.
Resolved BAM haplotypes bind to the complete Kestrel candidate representation they refine;
they do not create a second caller identity or increment equivalent-representation counts.
For the seven carriers and `pair_4092`, the golden oracle counts the absent canonical
delins as a miss and rejects any fabricated matching identity.

## Work package B: governed adVNTR artifact evidence

The active `Polymorphic_Call` list remains exactly 24 distinct reachable State strings.
No entry is added, removed, promoted, or supplied with invented cohort counts. The only
retained provenance is the recorded `renome` cohort, an adVNTR version older than 2.0.4
with unknown exact revision, no retained denominator/count/frequency, one
owner-confirmed state, and 23 states pending revalidation.

`Polymorphic_Call` asserts only this:

> The adVNTR State matches a carried-forward recurrent-state list whose bare State string
> is not sufficient to establish molecular identity.

It does not assert benignity, pathogenicity, population frequency, or allele equivalence.
The evidence file stores unknown values explicitly as JSON nulls.

Detection and interpretation remain separate:

- the adVNTR row remains a finding and retains `Positive (Flagged)` presentation;
- the row remains visible in reports and cohorts;
- the row does not independently withdraw a Kestrel finding;
- the row is inadmissible for molecular agreement and cannot promote Tier A;
- a cross-match may display the flagged observation but cannot call it molecular
  agreement; and
- if Kestrel is positive, the packaged decision is exactly the Kestrel-only name and tier;
- if Kestrel is negative, the flagged adVNTR finding remains visible with its existing
  `Positive (Flagged)` result and caller-local nomenclature, but reconciliation emits no
  molecular agreement and no Tier-A name; and
- a future explicit custom policy may abstain but cannot use the flagged state to promote,
  switch, or name an identity.

This policy applies to both confirmed and pending entries because the collision is a
property of the bare State key, not of provenance confidence. The evidence record marks
each entry `confirmed` or `pending`, stores null for every unknown count, denominator,
frequency, and exact tool revision, and has canonical bytes plus a checked-in SHA-256.
The loader verifies that digest, and summary schema 3 records it beside the profile hash.
The legacy `Polymorphic_Call` and `Positive (Flagged)` labels are retained compatibility
terms, not endorsements of polymorphism. Revalidation can change list membership only
through an owner-approved data revision carrying counts, denominators, software/model
provenance, and a new digest.

Between work packages B and C, the existing summary provenance gains the additive
`advntr_evidence_digest` field without a schema-version change. Work package C carries that
same value into schema 3 and makes it mandatory. Legacy summaries without it state that the
artifact-evidence revision was not recorded; they never infer the current digest.

## Work package C: explicit decision profiles

### Complete profile

Add a complete standalone `decision_profile.json` containing the exhaustively inventoried
Kestrel scoring, flagging, motif filtering, selection, adVNTR flagging, SHARK filtering,
nomenclature, and cross-match decision components. Every numeric decision field carries a
declared evidence unit and comparator. Runtime paths, executables, references, coverage
presentation, input routing, report text, the issue 266 reporting floor, and issue 268
model/version compatibility guards remain outside the profile in the existing whole-file
configuration or code. Generated profiles cannot supply clinical-sounding text.

Every configuration field belongs to exactly one validation class:

| Class | Explicit CLI profile | Generated profile | Failure |
| --- | --- | --- | --- |
| Excluded runtime/presentation/#268 field | absent | absent | presence is an unknown-key error |
| Fixed safety field | present and equal to packaged | present and equal to packaged | mismatch is an immutable-field error |
| Explicit-custom decision field | may differ | copied equal to packaged | malformed or missing is an error |
| Generated dominance/abstention field | may differ | may differ within the frozen grid | malformed or missing is an error |

Fixed safety fields include the issue 266 reporting floor, which is exactly the inclusive
`Depth_Score >= 0.00469` low gate; the independent GG minimum field, whose value is
validated equal to that low gate; the `20`, `21`, `100`, `200`, and `0.00515` Kestrel
boundaries; `bam_flank = 8`; `bam_thin_haplotype_record_support = 3`; and two separately
named Tier-A fields that both retain value `5`: Kestrel alternate k-mer-path depth and
adVNTR sequencing-read support. Thus a complete profile declares these values and their
units, but neither an explicit nor generated profile can alter them in this program.

The packaged profile embeds the current component values exactly. A no-argument run loads
that packaged profile. `--decision-profile PATH` selects exactly one complete custom
profile. There are no overlays, includes, inheritance, environment-variable selection,
assay-name discovery, current-directory discovery, or missing-field fallback.

Generated profiles are future explicit opt-ins. Profile kind is exactly `packaged` or
`generated`; persisted source is exactly `package` or `explicit-cli` and is never a path.
Generated profiles record `profile_kind: generated`,
the packaged base-profile hash, generator name/version, objective, dataset manifest hash,
partition manifest hash, and seed. Generation never selects the output.

### Resolution and provenance

The profile is parsed once before stage artifacts are created. Duplicate JSON keys,
nonfinite values, unknown or missing keys, invalid types/ranges/enums, comparator-rule
errors, rule/evidence mismatches, or incomplete components raise with completed-application
exit code 1. An explicit invalid path never warns and falls back.

Canonical profile bytes use RFC 8785 JSON canonicalization followed by one trailing
newline. SHA-256 covers exactly the payload bytes snapshotted under
`provenance/decision_profile.json`; the digest is never embedded in its own payload. The
summary records the profile identifier, revision, kind, source, hash, and relative
snapshot path in schema 3. Absolute operator paths are not persisted.

Kestrel, adVNTR, SHARK, nomenclature annotation, reconciliation, report generation, and
cohort processing receive immutable resolved components. Pipeline paths do not read
module globals after resolution. Compatibility wrappers accept an explicit resolved
component or load the packaged profile only when no run context exists; they raise if a
custom run context would otherwise be mixed with packaged values. A structural import
test and a CLI differential profile in which every explicit-custom and generated-mutable
component changes prove that no import-time or hidden packaged value participates in a
custom run.

Standalone report verifies the recorded snapshot. Legacy summaries state `decision
profile not recorded by legacy run`; they never claim the current packaged profile was
used. Cohorts retain each sample's profile hash and visibly group mixed profiles rather
than recomputing all samples with one profile. Pooled decision-performance aggregates
across different hashes are suppressed; each profile group is reported separately.

## Work package D: future calibration and dominance/abstention

Issue 269 is a future-facing feature, not a default recalibration. Its command operates on
completed runs and emits explicit custom profiles. It does not rerun Kestrel, edit the
packaged profile, select its output, or change any shipped cutoff.

The generated dominance/abstention allowlist contains only fields replayable from the
retained pre-result, identity-candidate, BAM-evidence, adVNTR, and provenance artifacts.
Every excluded, fixed-safety, and explicit-custom field is absent or copied according to
the table above. Caller invocation, reporting floor, record eligibility, upstream filters,
all listed Kestrel boundaries, BAM flank, thin haplotype-record support `3`, Kestrel
Tier-A alternate k-mer-path depth `5`, adVNTR Tier-A sequencing-read support `5`, and any
field whose alternative would require a caller rerun cannot be generated differently.
Thus #269 evaluates future dominance and abstention without recalibrating the current
cutoffs.

Version 1 provides four deterministic operations:

```text
vntyper calibrate extract --truth TRUTH --partitions PARTITIONS --runs RUNS --output EVIDENCE
vntyper calibrate fit --evidence EVIDENCE --objective lexicographic-safety-v1 --output CANDIDATE
vntyper calibrate validate --profile PROFILE --evidence VALIDATION --output VALIDATION_ATTESTATION
vntyper calibrate evaluate --profile PROFILE --evidence LOCKED_HELDOUT --output HELDOUT_ATTESTATION
```

The objective is mandatory. No objective, including F1, is silently selected. The primary
objective is the lexicographic tuple:

1. wrong Tier-A displayed identities;
2. normal-control findings;
3. all wrong displayed identities across every tier;
4. exact-recovery rate multiplied by `-1`, macro-averaged across assay and
   mutation-class strata;
5. binary-detection sensitivity multiplied by `-1`;
6. number of free policy parameters; and
7. canonical profile hash as a deterministic final tie-break.

All rates use the complete predeclared truth population as denominator. Abstention is a
non-recovery and a missed detection for these objectives; it is never removed from a
denominator. Wrong-name counts include every displayed name, and component 3 includes
component 1. A candidate is inadmissible before comparison if Tier A becomes unreachable,
its paired detection or macro exact-recovery difference from the packaged policy has a
one-sided 95% group-bootstrap lower confidence bound below the zero-percentage-point
non-inferiority margin, a required stratum is below its protocol-fixed sample count, its
abstention rate exceeds the protocol-fixed cap, or applicability metadata does not match.
The bootstrap resamples complete leakage groups 10,000 times with the recorded seed. The
candidate grid and maximum free-parameter count are frozen in the protocol before fitting;
fitting cannot add thresholds or split points.

### Admissible features

Allowed runtime-available features are canonical identity, motif/pair context,
context-local representation counts, alternate k-mer-path depth, active-region depth,
GDP/DP-derived values, structural gates, co-occurring identity hypotheses, resolved
haplotype-record counts/shares/margins/ties, typed XD availability and the summaries named
below, governed adVNTR disposition, adVNTR sequencing-read support/p-value/coverage, and
tool/reference/assay-class metadata.

For one identity, haplotype-record count is the number of eligible resolved records whose
canonical co-occurring edit set contains that identity; one record contributes at most one
vote to that identity but may support more than one co-occurring identity. The denominator
for share is all eligible resolved records consulted for the locus, so identity shares are
not probabilities and may sum above one. Count margin is top count minus runner-up count;
share margin divides that value by the same denominator. Missing BAM evidence remains
missing rather than zero. Permitted per-identity XD summaries are availability count,
availability fraction, median, and interquartile range across its supporting records. XD
retains the typed meaning “optional minimum k-mer depth” in every schema and report.

Forbidden features are truth labels, sample/path names, batch identifiers, post-decision
tiers or selected names, simulator-only allele length/position, reference-interval length
presented as a patient allele length, raw XD vote weighting, XD winner/tie breaking,
unresolved flag strings as reliable evidence, and raw `POS/REF/ALT` identity.

The candidate dominance model is the finite declarative rule grid frozen in the protocol.
Record-count ties always abstain. XD can only veto the whole locus for missingness,
predeclared extreme concentration, or predeclared discordance; it cannot promote, switch,
or fall through to another winner. Artifact-flagged adVNTR state and its support statistics
are masked from agreement and promotion; the governed disposition may only trigger a
whole-locus abstention. Abstention reasons are closed tokens, are reported separately from
negative detection, and still count as non-recovery in objective denominators.

### Leakage control and statistics

Training, policy selection, validation, and locked held-out evaluation are distinct roles
in role-bound immutable manifests with separate roots and hashes. Split groups include
biological individual/family, simulated mutated/normal pair, simulator backbone/seed lineage,
technical replicate, duplicated rerun, source sample reused across depth series, assay
batch, and repeat context. Any group crossing a split is a fatal input error. Assay class
may be a feature; batch identity may only be a split key. Transforms, caps, missing-value
rules, and feature selection are fitted on training data only; the frozen candidate grid
is compared on policy-selection data. Labels live in a separately keyed sidecar, and the
feature artifact is enforced against an allowlist. `fit` cannot open validation or
held-out roots, `validate` cannot select another candidate, and `evaluate` accepts only a
held-out-role manifest.

The profile hash and protocol hash are precommitted to an append-only record before the
custodian opens locked evidence. A failed validation retires the candidate; reusing that
set for another candidate makes it selection data. The independent custodian runs held-out
evaluation once. A consumption ledger keyed by evidence hash refuses a second local
evaluation, and any access outside the custodian run invalidates the set for issue closure.
Evaluation records seeds, software and reference versions, sample composition, assay,
depth, read length, array-size information when independently measured, mutation classes,
objectives, profile and protocol hashes, manifest and evidence hashes, access attempts,
ROC/PR curves, relevant joint surfaces, two-sided 95% intervals and the one-sided paired
non-inferiority bounds, boundary coverage, abstention rate and reason distribution by
split, fitted results, validation results, held-out results, and small-sample limitations.

The Phase 1 reference projection is 154 displayed names, 136 exact, 18 wrong; Tier A
53/53/0; zero control findings; 83 BAM-eligible samples; and 68 BAM fetches. Work package A
preserves the Kestrel selected-row projection—row order, caller-local name, confidence,
flag, depth, support, and tie outcome—while measuring every full-result metric change
caused by replacing name-string agreement with molecular-identity agreement. It records
that identity-policy projection as A's baseline. Work package B may then change only
agreement and tier for a row carrying `Polymorphic_Call`; it must leave detection,
displayed caller-local names, and every unflagged sample identical, and records the result
as the governed-policy baseline. Work package C reproduces the governed-policy baseline
exactly. Before reporting a sweep, the calibrator must reproduce that shipped baseline:
all aggregate and per-tier metrics plus selected row order, name, confidence, flag, tier,
support, tie outcome, and abstention state. It reports all per-tier
displayed/exact/wrong metrics, not only Tier A.

## External evidence and closure

The available corpus is loaded and usable: 200 simulated mutated samples, 200 paired
normal controls, complete Kestrel TSV/BAM and adVNTR artifacts, with the Phase 1 golden
suite at 32 passing tests and zero skips. It is development and regression evidence, not
locked external validation. The paired simulations and existing golden oracles have
already informed policy analysis.

The repository currently lacks independent real carriers with orthogonally established
exact alleles, representative independent negatives, multiple laboratories/assays/read
depths, reanalysed renome per-state counts, and a locked custodian-held cohort. Therefore:

- safe unit, profile, extraction, and internal-validation work can complete;
- no fitted profile is called validated from the available simulations alone;
- the packaged cutoffs remain unchanged regardless of exploratory performance;
- issue 295 remains open because locked external evidence is unavailable; and
- release preparation does not claim completion of issue 295 while that blocker remains.

The present program therefore does not close issue 295. If a qualifying external corpus
is later supplied, a new validation protocol must be approved and precommitted before any
member is opened. That protocol fixes the leakage groups, minimum sample count per required
stratum, abstention cap, multiplicity handling, and independently justified precision
target; these values cannot be inferred from or revised after seeing the cohort. Only then
may one frozen policy be evaluated once. Closure additionally requires zero observed
wrong Tier-A identities, zero observed control findings, no increase in wrong displayed
identities, and the zero-margin paired non-inferiority bound for exact recovery. Two-sided
95% Clopper-Pearson intervals for zero-event proportions and group-bootstrap intervals for
paired performance are reported with sample-size limitations. Reporting an interval is
not a clinical safety claim.

## Test and mutation contracts

Before each test, the implementation plan states the realistic production mutation it
catches. Every modified function receives a behavioral test, and every new unit file is
marked `unit`.

Required invariants include:

- the Phase 1 reference remains explicitly pinned at 154 displayed, 136 exact, 18 wrong,
  Tier A 53/53/0, zero control findings, 83 BAM-eligible samples, and 68 BAM fetches;
- work package A preserves the exact Kestrel selected-row projection and independently
  reports every full-result delta caused by identity-based reconciliation;
- work package B changes only the governed agreement/tier fields of rows carrying
  `Polymorphic_Call`, with every unflagged decision projection identical;
- work package C reproduces the governed-policy baseline exactly;
- canonical dupC remains exactly 96/96 and named `59dupC`;
- equivalent anchors and `G>GG`/`C>CG` shapes share one identity;
- `dupA` off-by-one and wrong motif anchors remain different identities;
- same raw tuple under different motifs remains distinct in the typed identity candidate
  set even though the packaged compatibility projection is unchanged;
- identity equality never uses display name, tier, depth, support, or raw tuple alone;
- context-divergent evidence remains unable to reach Tier A;
- structural delinsAT candidates are not fabricated;
- artifact-flagged adVNTR evidence cannot confer agreement or Tier A;
- the 24-state list remains duplicate-free, reachable, provenanced, and synchronized;
- no-profile execution preserves the exhaustive packaged decision manifest and canonical
  hash, including all units, comparators, and values;
- partial/malformed profiles fail without fallback or stage artifacts;
- profile hashing is independent of input key and row order;
- generated profiles never auto-activate;
- a custom CLI differential profile reaches every explicit-custom and generated-mutable
  component and never a packaged import-time value; excluded issue 268 fields are rejected
  as unknown, while changing issue 266 or any fixed-safety field, including the separately
  named `3` and two `5` fields, raises an immutable-field error;
- the complete shipped decision projection is reproduced before any sweep;
- every objective component has a mutation that changes the winner and fails;
- all-abstain and Tier-A-unreachable candidates are inadmissible;
- a more sensitive policy with one extra wrong Tier-A identity loses;
- a policy with fewer wrong names but one control finding loses;
- record ties abstain and unequal XD cannot break them;
- an XD veto with a runner-up abstains rather than selecting the runner-up;
- selective abstention remains in every recovery and sensitivity denominator;
- changing `3`, `5`, BAM flank `8`, any packaged Kestrel cutoff, comparator, unit, active
  region boundary, or linked floor fails; and
- issue 266 comments cannot become rows or findings, while issue 268 guards cannot be
  bypassed by a profile.

Golden verification must prove both external roots were loaded, report corpus counts and
all per-tier displayed/exact/wrong metrics, and count no skip as evidence. Independent
oracles derive expected identity from truth and canonical sequence using an independently
implemented normalizer; they do not import production predicates, canonicalization,
serialization, grouping, or decision code.

## Delivery and merge contract

The functional pull requests are sequential:

1. molecular identity and identity-aware candidate interfaces, referencing issue 270;
2. governed adVNTR evidence, referencing issue 267;
3. explicit complete decision profiles, referencing issue 175; and
4. future calibration plus dominance/abstention evaluation, referencing issues 269 and
   295.

Pull requests 1–4 use only `Refs #N`. In particular, pull request 4 must say `Refs #295`
and cannot contain a closing keyword for issue 295.

Each begins from the newly fetched `origin/main` after its predecessor merges. No stacked
feature branch is used. A non-release pull request may merge under the owner's standing
authorization only after its targeted tests, relevant unit families, explicit-root golden
tests, required repository gates, and GitHub checks are green with Critical and Important
findings resolved. Review is bounded to three tools-disabled Opus partitions and one GPT
decision-invariant partition, with one retry for an invocation failure, one controlled fix
wave, and at most one scoped re-review of that wave. Failed or unavailable evidence is
reported rather than normalized.

Owner approval of this design constitutes standing merge authorization for pull requests
1–4 under those gates; it does not authorize a release action. Release preparation is
separate and does not begin while the external-evidence blocker keeps issue 295 open. After
that blocker is resolved, it begins only when every required functional change is merged,
main CI is green, and issue closure is supported. Merging the release pull request,
creating or pushing a tag, and sending the production repository dispatch each require
fresh explicit authorization.

## Explicit non-goals

- changing any shipped cutoff;
- automatic local or assay-based profile selection;
- raw XD weighting;
- treating haplotype records as reads or molecules;
- using the available simulation as independent held-out evidence;
- rebuilding the adVNTR model or treating issue 268 as calibration;
- expanding subthreshold candidates into calls;
- inventing renome provenance or per-state frequencies;
- claiming recovery of structurally absent delinsAT candidates;
- workflow changes without separate approval; and
- tagging or publishing a release.
