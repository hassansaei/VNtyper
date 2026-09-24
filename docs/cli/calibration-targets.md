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

The packaged model predicts **complete total diploid repeat counts**: every repeat unit on
both alleles, including the nine invariant terminal units per allele. Its model file
declares `count_convention: complete`, and its PacBio truth counts follow the same
convention (canonical units plus 18). The report states the convention explicitly.
Standard length estimates remain research outputs and do not alter mutation calls,
confidence assignments, or screening conclusions.

### Sensitivity tiers

Kestrel's depth ratio falls as roughly one over the total array length, so short-read
frameshift detection becomes less sensitive as the total array grows. When a length
estimate is available, the pipeline records a sensitivity tier in
`pipeline_summary.json` for both length paths:

| Tier | Estimate | Summary code | Report |
|---|---|---|---|
| `below` | at or under the caution threshold | none | estimate with the model's typical error |
| `caution` | above `caution_threshold` | `vntr_length_exceeds_sensitivity_cutoff` | badge and help text |
| `high` | above `high_threshold` | the caution code and `vntr_length_high_sensitivity_risk` | badge, help text, and a header notice when the sample has no finding |
| `not-assessed` | no estimate | none | nothing |

The policy lives in the run configuration and is recorded next to the tier as
`length_sensitivity_policy`:

```json
"length_estimation": {
  "sensitivity": {"caution_threshold": 110.0, "high_threshold": 150.0, "typical_error_repeats": 14.3}
}
```

Comparisons are strict (`>`) on the point estimate rounded to the tenth the report
displays, so a report never shows a value that "exceeds" a cutoff it equals. They are made
in the complete-count frame the cutoffs are defined on. A `canonical-only` estimate is moved into that frame by the
packaged conversion (+9 terminal units per allele, +18 in total); a `source-reported`
estimate has no known frame and is recorded as `not-assessed`. Approved-path models are
complete by contract. For the packaged model only, the report shows
`typical_error_repeats` beside the estimate, worded as a typical error. It is never shown
as "± 14", because it is not an interval (see below). A configuration without the block
records no tier. Report wording comes from the `length_sensitivity` block of
`report_config.json`. The cutoffs are not part of the resume identity.

Basis for the defaults, measured on the 76 PacBio-truth exomes with leave-one-out
predictions:

- 82% of samples have a true total above 110 repeats, so 110 is a caution, not a banner.
- Predictions regress toward the mean (slope of prediction on truth 0.66): arrays above
  130 are under-estimated by about 12 repeats on average. The high trigger therefore sits
  below the length it targets: 150 on the estimate catches 7 of the 8 arrays truly above
  160, where simulated Kestrel sensitivity falls to 0.75 and lower.
- `typical_error_repeats` is the leave-one-out RMSE of the packaged model (14.26). It is
  a typical error, not an interval. Only 55 of the 76 leave-one-out residuals fall within
  ±14. The empirical 2.5% and 97.5% residual quantiles are −29.4 and +28.2, so an
  empirical 95% band would be about ±29 repeats. The report therefore words the value as
  "typical error ≈14 repeats" and never attaches it to the estimate as "±". Summaries
  recorded before this rename carry the same value as `uncertainty_repeats`, and the
  report still reads them.

#### Limits of the tier evidence

- **The cutoffs were chosen on the numbers that describe them.** 110 and 150 were picked
  on the same 76 leave-one-out predictions that the figures above report. No held-out
  set confirms them.
- **The high tier is imprecise.** At 150 on the estimate, the tier catches 7 of the 8
  arrays truly above 160 repeats (Wilson 95% CI 53–98%). It flags 15 of the 76 samples,
  and only 7 of those 15 are truly above 160.
- **The pooled length performance is selection-adjusted.** The leave-one-out MAE of 10.80
  repeats (95% bootstrap CI 8.8–13.0) was measured after several candidate models were
  compared on the same predictions. It is therefore optimistic for this model.
- **Transfer between cohorts is asymmetric.** Trained on the German cohort and applied to
  the French one, MAE is 10.1. Trained on the French cohort and applied to the German
  one, MAE is 27.1 with R² −0.69. See
  [#336](https://github.com/hassansaei/VNtyper/issues/336).
- **Out-of-range input still gets an estimate.** Features outside the training range add
  an extrapolation warning, but the estimate and its tier are still reported. The model
  contract binds no assay or aligner, so nothing stops other input from being scored.
  Input other than Twist exome data is out of scope, and its estimates and tiers are not
  supported by this evidence.

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

## Deriving cutoffs from labelled data

The original study selected the depth-score detection threshold `0.00469` by comparing
sensitivity and specificity across cutoffs (Figure 5C). Its separate high-confidence
definition required a depth score above `0.00515` **and** an alternate depth above `20`
(Figure 5B), and that band was chosen as a cautious label rather than swept. The paper
also limits its detection-threshold rationale to dupC. VNtyper's general reporting floor
and its `ALT=GG` gate are distinct decisions, so the workflow below replays the current
production rules instead of reconstructing the historical implementation. See
[Saei et al., iScience 26, 107171 (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10338300/).

### Candidate cutoffs come from the data

A threshold compared with `>=` only changes a decision at a value the data actually take,
so enumerating the observed values yields the complete set of distinct outcomes with no
resolution argument to defend. `calibration_cutoff_axes.derive_axis` collects the measured
statistic from every candidate row that passes all four structural gates
(`is_frameshift`, `is_valid_frameshift`, `motif_filter_pass`, `flag_filter_pass`) and
returns those values as the axis breakpoints, always including the baseline value so the
shipped operating point appears on every curve. Breakpoints the policy decoder refuses are
recorded with the decoder's own message rather than dropped.

Production floors are inclusive: a score equal to the floor passes. A threshold equal to
the largest observed value therefore still passes that value, so the observed values alone
never include the operating point that rejects every row. Each derived axis therefore adds
one **endpoint sentinel**: the next representable value above the largest observation on a
`>=` axis, or below the smallest on a `<=` axis (`+1`/`-1` on integer axes). The report
records it as `endpoint_sentinel` in the axis document. The sentinel is left out when the
baseline value already lies beyond the observed range, and it is recorded as rejected when
the policy decoder refuses it. It is kept when `--max-breakpoints` caps the axis. Without
it, a specificity floor that only "reject everything" can meet would be reported as
unreachable. The sentinel depends on the extreme observed value, so it is derived again
inside every outer fold from that fold's training samples (see
[Choosing an operating point](#choosing-an-operating-point)).

`calibration_cutoff_grid.build_cutoff_grid` remains available for an explicitly declared
grid. Declared values are a supplement to the derived breakpoints, not a replacement.

### What each axis moves

| Axis | Pointers moved | Endpoint it changes |
| --- | --- | --- |
| `depth_floor_linked` | reporting floor, depth-score low, GG depth-score gate | detection |
| `gg_gate_independent` | GG depth-score gate only | detection, link broken on purpose |
| `depth_score_high` | depth-score high | confidence labelling |
| `alt_depth_band` | alternate-depth low and mid-low | confidence labelling |
| `var_active_region` | active-region threshold | confidence labelling |
| `advntr_cutoff` | adVNTR `calibrated_calling/cutoff` | adVNTR detection (`p < cutoff`) |
| `advntr_min_support` | adVNTR `calibrated_calling/minimum_read_support` | adVNTR detection (`read support >= value`) |

The two adVNTR axes are described under [Which caller is searched](#which-caller-is-searched).

#### The Kestrel depth gates move together

Lowering the reporting floor on its own changes nothing. The ordered confidence table
sends a score below the floor to `Negative`, labels a score inside the closed mid-band
interval `Low_Precision`, and requires a score at or above the high boundary for every
remaining rule, so a score between a lowered floor and an unchanged mid-band edge matches
no rule and reaches the `Negative` fallback. The MUC1 dupC candidate is an `ALT=GG` row,
so it is held by the GG gate as well. The linked axis therefore sets the floor and the GG
gate to the breakpoint and the mid-band edge to the smaller of that breakpoint and the
high boundary, which keeps every admitted score inside a labelled band and keeps the
tightening arm of the curve reachable above the high boundary.

Reference, recruitment, frameshift, motif and artifact rules stay fixed. Complete captures
are replayed through every production gate and through selection; a final result TSV
cannot substitute for replay, because it no longer contains the candidates that were
filtered out.

### Which caller is searched

`vntyper calibrate optimize` searches Kestrel axes, adVNTR axes, or both
([#269](https://github.com/hassansaei/VNtyper/issues/269)).

- `--caller kestrel` (the default) searches the Kestrel axes against Kestrel calls. The
  default axis is `depth_floor_linked`.
- `--caller advntr` searches the adVNTR axes and scores adVNTR calls only. The default axis
  is `advntr_cutoff`. Kestrel captures are still required, because they carry the complete
  baseline policy of both callers. Kestrel is replayed at its baseline only, to prove
  Kestrel parity, and is not scored.
- `--caller both` searches the requested axes of both callers, one axis at a time, and
  scores the either-caller union. The default axes are `depth_floor_linked` and
  `advntr_cutoff`. Every candidate moves one caller's axis and holds the other caller at
  its baseline. When only Kestrel axes are requested, the adVNTR arm is replayed at its
  baseline for every candidate, and the report records
  `advntr_policy: "held-at-baseline"`.

`--caller advntr` and `--caller both` require `--advntr-executable`. `--caller kestrel`
refuses an adVNTR axis, and `--caller advntr` refuses a Kestrel axis. The report lists the
searched callers under `search_scope.searched_callers` and states the adVNTR policy
(`searched`, `held-at-baseline` or `not-evaluated`) on the HTML page.

A derived research profile runs adVNTR under fixed capture parameters (see
[Applying a derived profile](#applying-a-derived-profile)), so the adVNTR captures must have
been produced under those same parameters; only the thread count may differ. With
`--caller advntr` or `--caller both`, a run whose captures differ stops before any adVNTR
replay and names the differing fields.

#### The adVNTR axes

adVNTR's legacy frameshift caller scores a candidate only when its read support reaches
`minimum_read_support`, and calls it when the p-value is strictly below `cutoff`. A
sample is called when any assessable locus returns a call. Each axis moves one pointer;
every other adVNTR value, including the `legacy` mode, keeps its baseline value.

| Axis | Comparator | Per-sample statistic |
| --- | --- | --- |
| `advntr_cutoff` | call when `p < cutoff` (strict) | smallest p-value among visits whose read support reaches the baseline support |
| `advntr_min_support` | call when `read support >= value` (integer) | largest read support among visits whose p-value is below the baseline cutoff |

Because the cutoff comparison is strict, the candidate for an observed p-value `q` is the
next float above `q`, which is the smallest cutoff that calls that sample. The sentinel
that calls no sample is the smallest observed p-value itself. adVNTR requires a cutoff
above 0, so a sample with `p = 0` cannot be rejected by any cutoff. The sentinel is then
the smallest positive float, and the cutoff axis document states the number of such
samples as `unrejectable_samples`. The support axis uses the `>=` sentinel rule above
(`+1` beyond the largest observation). Values the decoder refuses are recorded under
`rejected` with the decoder's message.

The capture file cannot supply these statistics: it scores only the visits its own
baseline support admitted. The breakpoints therefore come from **one probe replay per
axis at a permissive projection**: the largest admissible cutoff below 1 at the baseline
support for `advntr_cutoff`, and support 1 at the baseline cutoff for
`advntr_min_support`. The probes and the baseline run in a separate native grid
(`advntr-probe`, at most three executions), and are never reused as candidate
executions. Unassessable samples contribute no statistic.

Every candidate is then replayed natively again, and the run stops unless four checks
hold:

- **Evidence binding.** The probe grid and the candidate grid replayed the same capture
  bytes, the same per-locus records and the same adVNTR tool.
- **Replay consistency.** At every tested candidate, the native calls equal the calls the
  legacy rule predicts from the probe statistics (`replay_consistency` in `report.json`).
  This proves agreement only at the tested values. Completeness is a separate property:
  the observed values are the complete set of distinct outcomes only for an uncapped
  search. Each axis document states `breakpoint_completeness: "complete"`, or
  `"capped-subsample"` when `--max-breakpoints` subsampled the full-data or any fold
  inventory; a capped search may miss distinct operating points.
- **adVNTR baseline parity.** The replay of the baseline candidate reproduces, sample by
  sample, the baseline calls recorded in each capture (`baseline_parity.advntr`).
- **Exact execution count.** The candidate grid executes exactly one native replay per
  distinct adVNTR policy among the candidates: one in all when only Kestrel axes are
  searched, and one plus the number of distinct adVNTR-axis candidates otherwise. The
  report publishes this count (`advntr_distinct_executions`) and the probe count
  (`advntr_probe_executions`) separately. The probe and candidate grid digests, the adVNTR
  tool identity and both wall times are recorded under `provenance.advntr`.

Cost is about one native replay execution per distinct adVNTR candidate. One execution
over about 80 captures took about 20 s on one workstation, and executions run serially,
so an uncapped search of both adVNTR axes can take over half an hour. `--max-breakpoints`
bounds the cost.

Not searched: adVNTR exact mode and its background fitting, the rare-unit coverage
guard, and pinning adVNTR v2.3.0. The replay accepts only a `legacy` baseline for cutoff
axes, and a non-legacy baseline stops the run with an error.

### Choosing an operating point

`calibration_cutoff_selection.SearchSpec` requires an objective. There is no default.

| Objective | Maximises | Notes |
| --- | --- | --- |
| `max-sensitivity-at-specificity` | sensitivity | requires `min_specificity`; refuses construction without it |
| `youden-j` | sensitivity + specificity - 1 | rank-equivalent to balanced accuracy |
| `max-f1` | F1 | |
| `balanced-accuracy`, `sensitivity`, `specificity` | the named rate | retained for compatibility |

For this application the useful objective is usually maximum sensitivity subject to a
specificity constraint, so a balanced score chosen by default would trade sensitivity away
silently. `select_cutoff_policy` reads training observations only. Ties prefer fewer false
positives, then more true positives, then fewer no-calls, then the baseline, then a stable
key. `vntyper calibrate optimize` uses each candidate's policy digest as that key, because
candidate identifiers number the merged inventory, which held-out samples help build;
called without keys, `select_cutoff_policy` compares the identifiers. A missing truth class prevents selection, and unsatisfiable constraints
return an explicit no-selection result rather than a nearest match.

Sensitivity and specificity keep no-calls inside their truth-class denominators, and
samples of unknown truth stay visible in their own denominator. Full-data operating points
describe the searched cohort; pooled held-out predictions assess fold-selected policies
separately. Neither is independent validation once the data or the search design has been
examined.

A data-derived axis is a function of the samples it is derived from: its breakpoints,
its endpoint sentinel and, under `--max-breakpoints`, which breakpoints the cap retains.
If a held-out sample took part in that derivation, its own values could create, remove
or displace a threshold its fold selects, and the held-out estimate would be optimistic.
Every axis of both callers is therefore derived once from all samples and once per outer
fold from that fold's **training** samples alone, each with its own sentinel and cap. The
replayed inventory is the union of these values. Inside a fold, only the baseline and the
candidates of that fold's training-derived inventory are admissible. Each axis document
counts the values that only a fold inventory produced (`fold_only_values`), each fold
record states how many candidates were admissible (`admissible_candidates`), and the
evaluation records `fold_admissibility: "training-derived-inventories"`. The full-data
selection, which produces the exported profile, still searches every replayed candidate.
Because the sentinel and the cap are now derived per fold, Kestrel held-out results can
differ from those of earlier releases on the same cohort.

On the either-caller union of `--caller both`, a positive call from one caller resolves
the other caller's no-call, so the set of no-call samples can change along one axis. No
ROC or PR curve with fixed denominators then exists. The axis curve is reported with
`status: "unavailable"` and the reason, and the page shows the reason instead of a
curve. The operating points stay in the cutoff table.

Outer folds are allocated stratified by truth label (positive, negative, unknown). Each
scored row is the only representative of its group. Groups are shuffled within each label
and dealt across folds in turn, so a class with at least as many groups as folds is present
in every fold. With 27 negatives, no fold can be left without negatives unless more than
27 folds are requested. The labels only steer allocation; each fold's selection still reads
training truth alone. Length evaluation keeps its unstratified allocation.

The HTML report shows the held-out estimate first. It gives the pooled held-out TP, FN,
TN and FP, sensitivity and specificity with exact 95% intervals, and each fold's policy
and admissible-candidate count. The intervals cover the pooled held-out calls and exclude
selection uncertainty. A warning appears when a held-out rate falls below the floor the
objective requested. The floor is enforced on training folds, so held-out performance can
fall short of it. A fold in which no candidate met the floor on its training data uses the
baseline instead, and the warning states how many folds did so. The selection and the table of every tested cutoff follow, labelled as
descriptive searched-cohort points. They are not validated performance.

### Endpoints are reported separately

A binary detection change, a confidence relabelling, an artifact flag and an exact-variant
identity are different outcomes, and a threshold can move one without moving another.
Because confidence participates in candidate ranking, a confidence threshold can change
which variant is selected. The nomenclature tier is assigned downstream of the capture
boundary and is therefore not available to this replay.

### Applying a derived profile

`vntyper pipeline --research-decision-profile PROFILE` runs the pipeline with cutoffs
derived locally. The profile is a complete caller-generated decision profile; no approval
artefact is involved, and the run logs that the values carry no deployment approval. The
flag is exclusive with `--decision-profile` and with `--calibration-bundle`, and neither of
those contracts changes: an explicit profile still may not alter a fixed-safety field, and
an approved portable bundle is still the only approved path.

A research profile that carries adVNTR values runs adVNTR with the derived legacy policy.
The native arguments are rendered from the profile under fixed, CLI-representable capture
parameters, the same ones the calibration captures were produced under: Illumina,
frameshift mode, diploid, maximum error rate 0.05, legacy error rate 0.01, MAPQ 0, base
quality 20, maximum low-quality fraction 0.1, enhanced HMM on, trained HMMs off,
reference alignment on, full-RU-only off, and no minimum read length. Exact mode is
refused for research profiles, because it needs an approved bundle with a fitted
background. `vntyper calibrate optimize` renders the adVNTR arguments of every exported
profile this way before it publishes the run, and with any `--caller`, including
`kestrel`, it refuses to export a policy whose adVNTR mode is not `legacy`. The run then
publishes nothing. A research profile carries no tool pin, so the runtime does not enforce an
adVNTR build. A derived cutoff is tied to the adVNTR build recorded under
`provenance.advntr` in the optimize report, and should be applied with that build.

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
