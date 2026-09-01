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
- Tier A 53/53/0, Tier B 101/83/18, Tier C 0/0/0;
- zero control findings;
- packaged profile SHA-256
  `be6329fb12107a1b6b65e425257be6233c7e2115e299e941c12a63a6a6d59718`; and
- packaged projection SHA-256
  `338fe05d010f623e794dcf93393904fa13bd8713e2d074c8a5b6c72d6efd96fd`.

`tests/golden/calibration_oracle.py` is recursively guarded against production decision,
canonicalization, serialization, grouping, and calibration imports. These figures must
reproduce before any candidate result is interpreted.

This is development replay, not validation. The corpus and policies have already been
examined, its existing run trees predate the complete schema-3 replay sidecar set required
for a new calibration extraction, and no independent custodian kept a locked cohort.
Consequently no value from this corpus is reported as validation or held-out performance.

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
400 members, and zero skips.
