# Calibration validation and the open external-evidence gate

The calibration engine is future-facing and opt-in. It evaluates dominance and
whole-locus abstention from retained artifacts; it does not recalibrate caller cutoffs,
change the packaged profile, alter a default run, or auto-select a generated profile.

## Development replay pinned in this repository

The explicit-root golden tier independently loads both the simulation and paired adVNTR
trees. It currently pins:

- 200 mutated members and 200 paired normal controls;
- 374 public identity rows and 178 selected locus rows;
- shipped displayed-name projection 154 displayed, 136 exact, and 18 wrong;
- the frozen historical Phase 1 caller projection recorded in the selected historical
  summaries, Tier A 53/53/0, Tier B 101/83/18, Tier C 0/0/0, kept as the pre-identity
  reference;
- the measured identity-policy projection of the current code on the same corpus,
  Tier A 0/0/0, Tier B 154/136/18, Tier C 0/0/0, pinned literally by
  `test_pr_a_measured_identity_policy_baseline_is_literal`: identity-keyed reconciliation
  deliberately keeps context-divergent evidence out of Tier A, so the total is unchanged
  while no name reaches Tier A on this corpus;
- zero control findings;
- packaged profile SHA-256
  `0b13d07370491b3ea773e65144891cb30caebcae70b0ef98feb0f2c5ccd2f4a1`; and
- packaged projection SHA-256
  `cfa5ec402a3a20096b76273c4347ff8b5975db942aa6dccf9f2d99474260236d`.

`tests/golden/calibration_oracle.py` is recursively guarded against production decision,
canonicalization, serialization, grouping, and calibration imports. These figures must
reproduce before any candidate result is interpreted.

This is development replay, not validation. The corpus and policies have already been
examined, its existing run trees predate the complete schema-3 replay sidecar set required
for a new calibration extraction, and no independent custodian kept a locked cohort.
Consequently no value from this corpus is reported as validation or held-out performance.

The golden tier therefore does not relabel those historical summaries. It asserts that
the four selected source summaries remain schema 2, then independently materializes the
smallest governed schema-3 fixture from literal facts in both roots:

- `train-pair-3000` supplies training facts;
- `select-pair-3001` supplies policy-selection facts;
- `validate-pair-3003` supplies validation facts; and
- `held-pair-3004` supplies only the value-free locked member/run declaration during
  ordinary extraction.

Simulation `ground_truth.csv` supplies the independently derived canonical identity and
displayed name. The paired adVNTR-root Kestrel result supplies motif context, depth,
support, confidence, flag, and tier. The oracle binds the exact SHA-256 of that truth file
and every selected historical summary/result file. It then invokes the real `vntyper
calibrate extract` and `vntyper calibrate fit` CLI operations, independently verifies the
installed checksum inventories and literal feature/baseline rows, and strictly resolves
the fitted generated profile.

Before evaluating the development grid, the gate proves the packaged neutral replay and
the complete historical 400-member figures above. A one-member grid cannot contest the
objective's ordering or tie-break; those are exercised by the objective's unit tests, not
by this gate. The one-member finite development candidate is
dominance enabled, record-count margin `1`, record share `0.5`, share margin `0.0`, and XD
veto `disabled`; its deterministic generated-profile SHA-256 is
`7fc772512eee789cd0a909fa939fd87305b4eb14a6769bcbf15b340afc7ce6f8`.
That file is loadable only by explicit selection and does not replace the packaged
default.

The fixture persists the classification
`previously-examined-development-simulation`, with
`eligible_for_independent_validation=false` and
`eligible_for_locked_evaluate=false`. It states explicitly that this is neither an
independent external cohort nor custodian-locked heldout evidence.

## What remains blocked

Issue #295 stays open until there are:

- independent real carriers with exact alleles established orthogonally;
- representative independent negatives;
- multiple laboratories, assays, depths, and read lengths;
- independently measured array-size context where available;
- reanalysed renome state counts; and
- a named custodian holding a never-opened locked cohort under durable authority outside
  this repository.

Before anyone opens that evidence, the protocol must be precommitted with the exact finite
grid, free-parameter limit, leakage groups, minimum per-stratum counts, abstention cap,
multiplicity method, precision target, seed, custodian, and one-use ledger authority.
Development-side extraction cannot mark evidence closure-eligible.

One frozen profile may then be evaluated once. Closure requires zero observed wrong Tier-A
identities, zero observed control findings, no increase in wrong displayed identities, and
both zero-margin paired non-inferiority lower bounds. Abstentions stay in every denominator.
All two-sided intervals, one-sided paired bounds, missing boundary cells, and small-n
limitations must be reported. Reporting an interval is not a clinical safety claim.

Until that gate passes, there is no issue-closure claim and no calibration-driven release.

## Running the development golden gate

```bash
: "${VNTYPER_SIM_ROOT:?set the simulation corpus root}"
: "${VNTYPER_ADVNTR_ROOT:?set the paired adVNTR root}"
pytest -m golden tests/golden -q -rs
```

A missing root or a skip is not evidence. The recorded verification must show both roots,
400 members, real `extract` plus `fit`, the fitted digest above, and zero skips. The
historical schema-2 bridge is a development regression proof only; it is not a substitute
for producing schema-3 artifacts prospectively under external custody.
