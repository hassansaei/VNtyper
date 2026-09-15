# Caller and VNTR length calibration

The explicit `callers` and `length` targets extend `vntyper calibrate` for issues
[#269](https://github.com/hassansaei/VNtyper/issues/269) and
[#332](https://github.com/hassansaei/VNtyper/issues/332). Omitting `--target` retains
[dominance calibration](calibrate.md) for the existing extract/fit/validate/evaluate
operations. Strict caller-policy activation requires a separately approved portable
bundle. Standard length estimation and development cohort comparisons are described below.

## Standard pipeline length estimate

The installed configuration enables a packaged research length model on compatible
GRCh38 alignments. It measures the MUC1 locus during the pipeline's existing alignment
lifetime, verifies the reference sequence, and records the estimate in the JSON and
HTML report. It performs one indexed locus traversal; inference requires no additional
variant caller or machine-learning runtime. For BAM input, `--reference-fasta` selects
the measurement FASTA explicitly; otherwise the configured assembly reference is used.
The FASTA must already have its `.fai` index; measurement never creates it implicitly.
CRAM uses the reference proven for decoding. Unsupported assemblies, unavailable reference
sequence, or insufficient denominator evidence produce an unavailable estimate with reasons.

```bash
vntyper pipeline --bam sample.bam --reference-assembly hg38 --output-dir result/
vntyper pipeline --bam sample.bam --no-estimate-vntr-length --output-dir result/
vntyper pipeline --bam sample.bam --reference-assembly hg38 \
  --standard-length-model fitted/length-model.json --output-dir result/
```

`length_estimation.enabled` controls the configuration default. The explicit
`--estimate-vntr-length` and `--no-estimate-vntr-length` options override it.
An explicit approved `--length-model` or measurement-only request selects the strict
path described later; conflicting explicit requests are rejected.

The standard model uses Bayesian ridge regression fitted on 13 depth/read features:
A, F, core zero-depth fraction, core depth variation, variation among eight core bins,
core-half depth balance, invariant-end balance, flank-end balance, invariant depth,
MAPQ-zero fraction, mean MAPQ, soft-clipped-read fraction, and whole-query GC fraction.
Its scaler and regression are collapsed into an intercept and 13 coefficients stored
in closed, checksummed JSON. Extrapolation is flagged; estimates are not clipped to a
training maximum.

The packaged model predicts **assay-reported total diploid repeat counts**. Its
`source-reported` count convention does not assert whether the assay includes invariant
terminal units. The report states that convention explicitly. Reference measurement
geometry includes invariant units, which is distinct from the truth assay's counting
convention. Standard length estimates remain research outputs and do not alter mutation
calls, confidence assignments, or screening conclusions.

## Development calibration with optional truth

Install the training extra when fitting length models:

```bash
pip install 'vntyper[calibration]'
vntyper calibrate cohort --manifest samples.tsv --reference reference.fa \
  --target both --folds 5 --output comparison/
vntyper calibrate cohort --manifest mutation-only.tsv --target callers \
  --caller-runs pipeline-results/ --output mutation-comparison/
```

The local TSV requires `sample_id`, `bam`, and `assembly`. Optional columns are
`genotype` (`positive`, `negative`, `unknown`, or empty), `allele_1`, `allele_2`,
`group_id`, `kestrel_result`, and `advntr_result`. Paths are relative to the TSV.
For example, one exact paired length and one mutation-only observation can coexist:

```text
sample_id	bam	assembly	allele_1	allele_2	genotype
example-a	reads/a.bam	hg38	40	65	positive
example-b	reads/b.bam	hg38			negative
```

Missing length truth leaves mutation eligibility intact. An incomplete allele pair is
audited and excluded from diploid fitting. Missing mutation truth stays unknown.
Duplicate alignment contributions are rejected. `group_id`, when supplied, describes
biological relationships for leakage control; it is not a processing label or predictor.
The default `--target auto` selects targets with available truth. A reference is needed
only for length extraction. `--count-convention source-reported|complete|canonical-only`
labels the supplied counts; this command performs no implicit terminal-unit conversion.

Length output includes held-out MAE, RMSE, R², bias, availability, comparison against
training-fold mean predictions, and a locally fitted research `length-model.json`.
Fitting and scaling use only the training fold. `--folds` controls the number of folds;
the default is five. Per-row extrapolation warnings remain visible without clipping or
excluding estimates. Previously examined data remain exploratory regardless of fold count.

Caller output reports unchanged native baseline calls and, when supplied through
`--caller-policies`, a finite inventory of actual native outputs under alternative
policies. Selection uses training groups only. Reports include held-out confusion counts,
sensitivity, specificity, false positives, no-calls, per-caller results, paired differences,
uncertainty, and descriptive cutoff/ROC operating points. Missing files are no-calls;
absent truth classes have undefined rates. Baseline-only results cannot demonstrate a
benefit from calibration. Native-positive union results are explicitly identified and
are separate from screening conclusions.

This command creates an offline `report.html`, `metrics.json`, checksums, and an optional
research length model. It does not create validation authority or approve a caller bundle.

### Supplying caller policy comparisons

`--caller-runs` reads `<root>/<sample_id>/kestrel/kestrel_result.tsv` for the baseline.
Explicit `kestrel_result` and `advntr_result` manifest columns take precedence. Declaring
adVNTR results or finding a conventional `advntr/output_adVNTR_result.tsv` under any
declared sample includes adVNTR in the required composition for the entire roster.
A missing required caller remains a no-call even when the other caller is positive.

For cutoff comparisons, pass `--caller-policies policies.json`. The closed JSON inventory
contains `schema_version: "cohort-caller-policies-v1"`, a `baseline` policy ID,
`required_callers` (`["kestrel"]` or `["kestrel", "advntr"]`),
`training_scope: "fixed-before-cohort"`, and a nonempty `policies` list. Each policy has:

- `policy_id`: a unique name;
- `cutoff` and `comparison`: a finite cutoff and `<`, `<=`, `>`, or `>=`, or both `null`;
- `samples`: the exact manifest sample roster, each mapped to `kestrel_result` and
  `advntr_result` paths relative to the JSON file, with unavailable paths represented by `null`.

These must be actual native outputs from policies fixed before this comparison. The
command does not reconstruct excluded variants from final TSVs. Use complete captures
and the strict caller workflow below to produce replayable evidence. A background model
trained on this cohort cannot be labelled `fixed-before-cohort`; fit it inside each
training fold or evaluate a model fixed on separate data.

## What each target estimates

| Target | Inputs to the scientific model | Output |
| --- | --- | --- |
| `callers` | Complete Kestrel evidence and, when declared, adVNTR capture/native results | One predeclared caller policy, with fixed-denominator sensitivity, specificity, false-positive rate, no-call rate, and identity accuracy |
| `length` | Reference-bound core/invariant or array/flank depth ratios | Total diploid repeat count, including the configured invariant units; an unavailable result when required evidence is missing |

Length estimation does not resolve two allele lengths, assign a mutation to an allele,
or alter the genotype. Its fitted family is affine regression on either ratio, with
an intercept and a training-only mean baseline. Physical A/F hypotheses are reserved
in the protocol but currently remain fit-ineligible and cannot be loaded as runtime
models. This strict affine family is separate from the packaged 13-feature research model.

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
