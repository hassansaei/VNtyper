# Caller and VNTR length calibration

The explicit `callers` and `length` targets extend `vntyper calibrate` for issues
[#269](https://github.com/hassansaei/VNtyper/issues/269) and
[#332](https://github.com/hassansaei/VNtyper/issues/332). Omitting `--target` retains
[dominance calibration](calibrate.md). A fit produces a research candidate. Runtime
activation requires a separately approved portable bundle.

## What each target estimates

| Target | Inputs to the scientific model | Output |
| --- | --- | --- |
| `callers` | Complete Kestrel evidence and, when declared, adVNTR capture/native results | One predeclared caller policy, with fixed-denominator sensitivity, specificity, false-positive rate, no-call rate, and identity accuracy |
| `length` | Reference-bound core/invariant or array/flank depth ratios | Total diploid repeat count, including the configured invariant units; an unavailable result when required evidence is missing |

Length estimation does not resolve two allele lengths, assign a mutation to an allele,
or alter the genotype. Its fitted family is affine regression on either ratio, with
an intercept and a training-only mean baseline. Physical A/F hypotheses are reserved
in the protocol but currently remain fit-ineligible and cannot be loaded as runtime
models. No fitted coefficients ship as a default model.

Metadata describes applicability and supports separate confounding checks. Assay,
processing labels, specimen identifiers, and source labels are not predictive features.
Processing groups do not define automatic fitting or validation splits. Specimen,
family, readset, and simulation identities enforce duplicate and leakage exclusions.

## Length measurements and count conventions

For one depth-counting policy:

- **A** is mean core depth divided by mean invariant depth.
- **F** is mean array depth divided by mean external-flank depth, weighted by the
  number of measured flank positions.
- Coverage, depth, and distinct-fragment support determine whether a measurement is
  usable. Missing measurements stay in availability denominators.
- A fit records its training domain. Extrapolation beyond the configured margin is
  unavailable; predictions are not clipped into the training range.

The packaged reviewed GRCh38 annotation is
`vntyper/data/length/grch38-length-annotation-v1.json`. Coordinates are zero-based,
half-open:

| Region | Interval on chr1 |
| --- | --- |
| Complete target | `[155188486, 155192239)` |
| Variable core | `[155188726, 155191939)` |
| Invariant ends | `[155188486, 155188726)` and `[155191939, 155192239)` |
| Array feature interval | `[155188529, 155192010)` |
| External flanks | `[155188296, 155188486)` and `[155192239, 155192429)` |

The reference contains five terminal units on one side and four on the other. The
public method also distinguishes canonical units from these nine terminal units.
See the [method and supplement](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458345/)
and [public repeat dictionary](https://github.com/pristanna/muc1repeats/tree/74a8bab867ef984798991dcd9c848faef4371cc9).

The array feature interval is not identical to the complete target. The core's base
length is not assumed to be an integer multiple of 60. Consequently this annotation
supports the affine features and disables both physical hypotheses.

Truth ingestion requires a declared counting convention. The packaged
`vntyper/data/length/grch38-count-conversions-v2.json` defines a signed **+9 per
allele** conversion for canonical-only counts. Counts already including terminal
units must not receive that addition. A reported value is not assumed to use either
convention merely because it came from a particular assay. There is no fixed upper
allele-length cap.

## Freeze the study before examining model outcomes

A `calibration-study-v2` declaration commits to:

1. The target, finite candidate family, objective, seed, QC rules, and acceptance rules.
2. Training, policy-selection, validation, and external locked-heldout membership.
3. The original caller baseline, or the length training-only baseline recipe.
4. Applicable reference, assay, input scope, preprocessing identity, and tool versions.
5. An external exposure-ledger identity shared across targets.

The objectives are `caller-safety-v1` and `length-total-v1`. The immutable protocol
specifies them; CLI arguments cannot replace the protocol after outcomes are read.

Each `calibration-runs-v2` entry commits to the actual input, policy, producer,
capture policy, process exit status, and named asset bytes. Assets have absolute local
paths, byte sizes, and SHA-256 hashes. Native TSVs, complete captures, and measured
length-feature artifacts are verified when their role is opened. A failed process
is not a negative genotype.

Each role has a `calibration-role-source-v2` document containing its eligible roster,
non-outcome exclusions, audited identity commitments, evidence domain, and a separate
sealed truth asset. QC failures are not truth exclusions. The external custodian
supplies the locked role; ordinary intake does not open locked reads or truth.

The Python initialization APIs create external state before the study is declared:

```python
from pathlib import Path
from vntyper.scripts.calibration_exposure_io import initialize_exposure_ledger
from vntyper.scripts.calibration_target_custody import initialize_target_custody

ledger_id = initialize_exposure_ledger(Path("/private/calibration/exposure.jsonl"))
initialize_target_custody(
    Path("/private/calibration/custody"), exposure_ledger_id=ledger_id
)
```

Use a new absolute location outside repositories and study/output directories. Put
`ledger_id` into the study declaration. The same ledger must accompany subsequent
studies that could reuse specimens. Do not reset it to make examined data appear new.

## Extract and fit

Extraction copies and verifies metadata without opening sealed outcomes:

```bash
vntyper calibrate extract --target length \
  --study study.json --sources sealed-sources/ --runs runs.json \
  --length-annotation annotation.json --output evidence/

vntyper calibrate fit --target length --objective length-total-v1 \
  --evidence evidence/ --exposure-ledger /private/calibration/exposure.jsonl \
  --output candidate/
```

For callers, use `--target callers`, omit `--length-annotation`, and fit with
`--objective caller-safety-v1`. The source directory contains
`roles/<role>/source.json` for all four declared roles. Truth assets remain separate.
When the protocol includes an exact adVNTR policy, fitting also requires
`--advntr-executable /path/to/advntr`. Its installed version and build must match the
frozen capture producer. The controller fits the background from authorized training
captures and negative truth only; known positive training observations contribute
diagnostics. Unknown VNTR lengths remain null diagnostic metadata. Every exact
selection run must use the resulting portable background bytes.

Fitting records exposure before opening training or selection outcomes. Length fitting
uses only training truth to estimate coefficients and the baseline. Selection compares
the predeclared hypotheses on policy-selection evidence. Caller fitting evaluates the
finite declared policies and verifies baseline replay against independently parsed
native results. Candidate changes that cannot be replayed require complete recaptured
runs. adVNTR background and capture adapters use the installed native implementation.

A successful fit writes `candidate.json`, a payload manifest, runtime payload files,
local scientific evidence, and an offline report. If no candidate qualifies, the
completed failed report is retained and the command exits 1. Fit outputs are not
runtime-approved bundles.

## Validate one fixed candidate

A confirmation evidence directory contains `source.json` and `runs.json`; its source
role must be `validation`:

```bash
vntyper calibrate validate --target length \
  --profile candidate/ --evidence validation-inputs/ \
  --exposure-ledger /private/calibration/exposure.jsonl \
  --custody /private/calibration/custody --output validation-result/
```

Validation does not refit, retune, or select a replacement. It records exposure and
claims the candidate before reading outcomes. The output includes `metrics.json`,
`report.html`, `exposure-receipt.json`, `validation-attestation.json`, and checksums.
A scientific failure is recorded as failed. An interrupted or operationally failed
claim cannot be retried as a fresh validation.

## Locked evaluation and export

After passed validation, an external custodian authorizes the exact candidate,
validation attestation, locked source/run commitments, and locked truth payload hash:

```bash
vntyper calibrate evaluate --target length \
  --profile candidate/ --evidence locked-inputs/ \
  --validation validation-result/validation-attestation.json \
  --authority custodian-authority.json \
  --exposure-ledger /private/calibration/exposure.jsonl \
  --custody /private/calibration/custody --output locked-result/

vntyper calibrate export --target length --profile candidate/ \
  --validation validation-result/validation-attestation.json \
  --evaluation locked-result/locked-attestation.json \
  --authority custodian-authority.json \
  --completion locked-result/completion.json --output approved-length/
```

Locked consumption is durably recorded before the truth bytes are read. A successful
locked outcome produces the completion receipt needed for export. Export checks
candidate, validation, locked result, authority, completion, and actual payload bytes.
It exports runtime parameters and aggregate-free approval bindings. Research reports,
identifiers, paths, training rows, and arbitrary background provenance are excluded.

Custody and exposure protect cooperating local workflows against accidental reuse,
concurrent claims, and interrupted retries. They rely on a trusted operator preserving
the external state; hashes do not authenticate a person or prevent an administrator
from replacing all historical files.

## Assess previously examined development evidence

Development assessment is explicitly nonpromotable:

```bash
vntyper calibrate assess --target length --profile candidate/ \
  --intake prepared-development-source/ --runs development-runs.json \
  --exposure-ledger /private/calibration/exposure.jsonl \
  --output development-assessment/
```

`--intake` here is a separately prepared `calibration-development-source-v1` bundle,
not the combined raw intake JSON. It declares previously examined representatives,
their audited identities, and a sealed target-specific truth asset. Assessment records
exposure, evaluates the fixed candidate, and writes a report. It cannot emit promotion
authority or become validation through a filename or role change.

## Reports and acceptance

Caller reports include confusion counts, sensitivity/specificity, false-positive and
no-call rates, exact-identity accuracy, confidence intervals, and cutoff/ROC/PR
summaries for the evaluated finite policies. Operational failures and no-calls retain
explicit status. Multidimensional policy comparisons are not presented as a fabricated
monotonic scalar ROC curve.

Length reports include absolute and relative errors, bias, availability, tolerance
coverage, and paired comparison with the training-only mean baseline. Missing and
out-of-domain predictions remain visible and count against availability.

The frozen caller policy checks wrong Tier-A identities, the upper confidence bound
on false-positive rate, paired sensitivity, and no-call increase. Selection additionally
requires a declared benefit. Length gates check absolute error, improvement over the
baseline, paired error bounds, tolerance coverage, availability, and independent-group
counts. Insufficient evidence is not a pass. Tiny invented test fixtures may declare
reduced thresholds to exercise orchestration; those settings are not evidence of
accuracy on external samples.

## Activate an approved bundle

```bash
vntyper pipeline --bam reads.bam -o results/ \
  --length-model approved-length/ --length-context length-context.json

vntyper pipeline --bam reads.bam -o results/ --extra-modules advntr \
  --calibration-bundle approved-callers/ --calibration-context caller-context.json
```

Caller context must match the bundle's frozen applicability and native capture policy.
Enable the adVNTR extra module explicitly when the approved caller bundle includes it.
Native model bytes and installed adVNTR capabilities are checked independently.
Caller-v2 profiles cannot bypass bundle approval through `--decision-profile`. Run-local
snapshots preserve the verified portable bundle; reports check the recorded bundle
identity and exact profile bytes.

For feature measurement without a fitted model, use
`--measure-vntr-length-features --length-annotation annotation.json --length-context
length-context.json`. Measurement and model contexts record provenance and applicability;
they do not introduce metadata predictors.

All outputs use private atomic directory publication and refuse occupied destinations.
Argument/usage errors exit 2. Completed scientific failures exit 1 with their reports.
Malformed or unexecutable operations exit 1 without installing a partial output tree.
