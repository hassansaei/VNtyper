# calibrate

Build and evaluate future, explicit opt-in dominance profiles from completed VNtyper
runs. Calibration never reruns Kestrel or adVNTR, edits the packaged profile, changes a
shipped cutoff, or activates its output automatically.

## Operations

```text
vntyper calibrate extract --truth TRUTH --partitions PARTITIONS --runs RUNS --output EVIDENCE
vntyper calibrate fit --evidence EVIDENCE --objective lexicographic-safety-v1 --output CANDIDATE
vntyper calibrate validate --profile PROFILE --evidence VALIDATION --output VALIDATION_ATTESTATION
vntyper calibrate evaluate --profile PROFILE --evidence LOCKED_HELDOUT --output HELDOUT_ATTESTATION
```

All options shown above are required. `fit` has no default objective: the only version-1
choice is `lexicographic-safety-v1`, and it must equal the objective already hashed into
the evidence protocol. A usage error exits 2; a parsed operation that fails validation
exits 1. Output is built as a sibling staging directory and appears at `--output` only
after every artifact and checksum succeeds.

`extract` accepts a strict study declaration with the finite candidate grid, seed,
bootstrap and multiplicity rules, four role-bound partitions, and all seven leakage-group
namespaces. Truth labels are stored separately from allowlisted runtime features. Runs
must retain schema-3 `pipeline_summary.json`, `kestrel/kestrel_pre_result.tsv`, and
`kestrel/bam_identity_replay.v1.json`; callers and BAMs are never reopened.

`fit` can read only training and policy-selection role payloads. It must reproduce the
shipped aggregate, per-tier, and ordered row projection before evaluating the frozen grid.
The emitted `decision_profile.json` is complete, generated, and hash-bound to the base
profile, protocol, dataset, partitions, seed, objective, and generator version.

`validate` evaluates one fixed profile on validation evidence without selecting another.
`evaluate` accepts only externally attributed, locked-held-out evidence. It writes a
durable precommit before opening the payload and consumes an evidence hash once, including
failed attempts. Local append-only receipts are safeguards against accidental reuse; they
are not proof of independent custody.

## Using a generated profile

A generated profile changes no run by existing on disk. Select it explicitly:

```bash
vntyper pipeline --bam sample.bam -o results/ \
  --decision-profile candidate/decision_profile.json
```

Only the six dominance/whole-locus-abstention leaves may differ. Fixed safety fields,
including `0.00469`, `0.00515`, `20`, `21`, `100`, `200`, BAM flank `8`, thin
haplotype-record support `3`, and both Tier-A support values `5`, remain immutable.
XD is always an optional minimum k-mer depth, never a read or vote weight.

## Evidence status

The available 200-mutated/200-control simulation corpus is previously examined
development evidence. It is useful for shipped-profile reproduction and regression
testing, but it is not independent external validation and must not be passed off as
locked held-out evidence. Issue #295 therefore remains blocked on qualifying external
data and a named external custodian. Reporting an interval is not a clinical safety
claim.
