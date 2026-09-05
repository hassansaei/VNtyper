# calibrate

Build and evaluate explicit opt-in dominance profiles from completed VNtyper runs. Calibration never re-executes Kestrel or adVNTR, never modifies the packaged profile, never alters shipped cutoffs, and never activates outputs automatically.

## Operations

```text
vntyper calibrate extract --truth TRUTH --partitions PARTITIONS --runs RUNS --output EVIDENCE
vntyper calibrate fit --evidence EVIDENCE --objective lexicographic-safety-v1 --output CANDIDATE
vntyper calibrate validate --profile PROFILE --evidence VALIDATION --output VALIDATION_ATTESTATION
vntyper calibrate evaluate --profile PROFILE --evidence LOCKED_HELDOUT --output HELDOUT_ATTESTATION
```

All options listed above are required. `fit` accepts only `lexicographic-safety-v1`, which must match the objective hashed into the evidence protocol. Argument parsing errors exit 2. Unexecutable operations exit 1 without writing partial outputs. Failed `validate` or `evaluate` runs exit 1 after recording attestations, reports, checksums, and retirement logs. Outputs stage in temporary directories and appear at `--output` atomically.

`extract` parses a study protocol containing candidate grids, seeds, bootstrap rules, role-bound partitions, and seven leakage-group namespaces. `--truth` contains canonical training, policy-selection, and validation label records separated from allowlisted runtime features. Duplicate, reordered, or unknown fields trigger rejection. Input runs must retain schema-3 `pipeline_summary.json`, `kestrel/kestrel_pre_result.tsv`, `kestrel/bam_identity_replay.v1.json`, caller outputs, and profile snapshots. Candidate selection evaluates the `dominance.enabled` flag in the verified `provenance/decision_profile.json` snapshot: packaged snapshots must match package profiles byte-for-byte; alternative snapshots re-parse against the packaged base.

Standard extraction does not inspect locked run roots and writes no locked features, labels, or manifests. Its `roles/locked-heldout/` folder contains only member declarations, run commitments, and checksums for custodial review.

`fit` processes training and policy-selection payloads, reproducing aggregate, per-tier, and ordered row projections before searching the grid. The generated `decision_profile.json` binds cryptographically to base profile, protocol, dataset, partitions, seed, objective, and version hashes.

Calibration grids support only `disabled` and `missingness` XD vetoes. Discordance and concentration modes require scalar record inputs not present in retained evidence and trigger rejection if declared.

`validate` scores a generated profile against validation evidence. The profile must match study parameters exactly. Failure retires the profile and evidence pair.

`evaluate --evidence` requires an external custodian package containing locked payloads, validated attestations, study bindings, authority signatures, and cryptographic checksums. Evaluation records a precommit before inspecting locked bytes and evaluates the pair once. Execution states finalize before releasing the advisory file lock. POSIX advisory locks serialize cooperating local processes; they do not guarantee cross-host coordination or non-POSIX lock behavior.

## Using a generated profile

A generated profile remains inactive until explicitly referenced. The neutral packaged profile remains default. To apply a generated profile:

```bash
vntyper pipeline --bam sample.bam -o results/ \
  --decision-profile candidate/decision_profile.json
```

Only six dominance and whole-locus-abstention leaves may vary. Fixed safety parameters (`0.00469`, `0.00515`, `20`, `21`, `100`, `200`, BAM flank `8`, thin haplotype support `3`, and Tier-A support `5`) remain immutable. XD acts solely as an optional minimum k-mer depth: never as a read count, molecule metric, vote weight, or winner selector.

## Evidence status

The 200-mutated and 200-control simulation dataset is classified as `previously-examined-development-simulation`, with `eligible_for_independent_validation=false` and `eligible_for_locked_evaluate=false`. It supports regression testing and profile reproduction, but does not qualify as an independent external cohort or custodian-held validation set. Issue #295 remains blocked pending qualified external data from an independent custodian. Reporting a confidence interval constitutes no clinical safety claim.
