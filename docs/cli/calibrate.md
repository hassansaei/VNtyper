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

These are the complete four commands; there is no command that creates custodian
authority. All options shown above are required. `fit` has no default objective: the only
version-1 choice is `lexicographic-safety-v1`, and it must equal the objective already
hashed into the evidence protocol. argparse usage errors exit 2. A parsed operation that
cannot complete exits 1 without installing a partial output. A completed failed
`validate` or `evaluate` also exits 1, but only after installing its complete failed
attestation, reports, checksums, and retirement record. Every output is built as a sibling
staging directory and appears at `--output` atomically.

`extract` accepts a strict study declaration with the finite candidate grid, seed,
bootstrap and multiplicity rules, four role-bound partitions, and all seven leakage-group
namespaces. `--truth` is ordinary, nonlocked truth: it contains exactly the canonical
training, policy-selection, and validation label rows, separately from allowlisted runtime
features. Supplying a locked label, an unknown field, a duplicate, or a reordered row is
an error. Runs must retain schema-3 `pipeline_summary.json`,
`kestrel/kestrel_pre_result.tsv`, `kestrel/bam_identity_replay.v1.json`, fixed caller
results, and the recorded profile snapshot; callers, BAMs, and CRAMs are never reopened.

Ordinary extraction does not open the locked run root and does not write locked features,
labels, baseline, plaintext payload, evidence manifest, or custody claim. Its
`roles/locked-heldout/` directory contains only the value-free member declaration, exact
run commitments, and their checksums for a future custodian. A declaration naming an
external-custodian role is not itself external evidence.

`fit` can read only training and policy-selection role payloads. It must reproduce the
shipped aggregate, per-tier, and ordered row projection before evaluating the frozen grid.
The emitted `decision_profile.json` is complete, generated, and hash-bound to the base
profile, protocol, dataset, partitions, seed, objective, and generator version.

Calibration-v1 grids accept only `disabled` and `missingness` XD vetoes. Concentration
and discordance remain valid runtime profile modes, but their per-record decisions cannot
be reproduced losslessly from the retained scalar evidence, so a study that declares
either mode fails before evidence extraction.

`validate` evaluates one fixed generated profile on validation evidence without selecting
another. The profile must bind the exact study, protocol, partition, dataset, objective,
and seed. A completed failure retires that profile/evidence pair.

`evaluate --evidence` is not an ordinary extraction directory. It is a separately
supplied custodian import containing a locked payload, passed validation attestation,
study binding, named external authority attestation, and checksums that bind the exact
study, partition, profile, runs, and payload. A locally minted closure-shaped directory is
rejected. Evaluation writes a durable precommit before opening locked bytes and consumes
the profile/evidence pair once. Success, completed failure, interruption, and exception
all reach a durable terminal state before the operation lock is released. Local file
locks, precommits, receipts, and append-only retirement records prevent accidental local
reuse; they do not prove that the cohort was independently selected or externally held.

## Using a generated profile

A generated profile changes no run by existing on disk. The verified packaged neutral
profile remains the default. Select a generated profile explicitly:

```bash
vntyper pipeline --bam sample.bam -o results/ \
  --decision-profile candidate/decision_profile.json
```

Only the six dominance/whole-locus-abstention leaves may differ. Fixed safety fields,
including `0.00469`, `0.00515`, `20`, `21`, `100`, `200`, BAM flank `8`, thin
haplotype-record support `3`, and both Tier-A support values `5`, remain immutable.
XD is always an optional minimum k-mer depth, never a read or vote weight.

## Evidence status

The available 200-mutated/200-control simulation corpus is classified
`previously-examined-development-simulation`, with
`eligible_for_independent_validation=false` and
`eligible_for_locked_evaluate=false`. It is useful for shipped-profile reproduction and
regression testing, but it is neither an independent external cohort nor
custodian-locked heldout evidence. Issue #295 therefore remains blocked on qualifying
external data and a named external custodian. Reporting an interval is not a clinical
safety claim.
