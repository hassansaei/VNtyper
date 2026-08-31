# Kestrel BAM Evidence Semantics — Phase 1 Design

Date: 2026-08-31

Issue: [#295](https://github.com/hassansaei/VNtyper/issues/295)

Status: proposed for explicit approval after adversarial review

## 1. Purpose

This phase corrects the evidence ontology used for Kestrel's `output.bam` without changing the decisions built on that evidence. The BAM contains Kestrel-resolved haplotype records, not sequencing reads. Its `XD` auxiliary tag is the minimum k-mer depth of the resolved haplotype, not a read count.

The implementation will make those meanings explicit in production code, typed interfaces, configuration, tests, reports, and documentation while preserving the present:

- edit voting and tie abstention;
- VCF-primary/BAM-refinement policy;
- Kestrel caller row/candidate selection and BAM-consultation eligibility;
- numeric thresholds (`3` for thin BAM evidence and `5` for Tier A);
- winning allele and emitted name;
- confidence tier;
- caller-independence model (`kestrel_vcf` and `kestrel_bam` remain one Kestrel caller).

This is semantic containment only. The change will say `Refs #295`; it will not claim to close #295.

## 2. Evidence established during discovery

### 2.1 Pinned Kestrel contract

VNtyper's vendored Kestrel 1.0.1 JAR has SHA-256:

```text
ceca2682578d641323405b3e605c85eebbe181574d056f0fabdf0ea2bf0bd1f1
```

Inspection of that exact artifact establishes that its SAM writer emits one SAM record for each resolved `Haplotype`. The writer emits `XD:i:<value>` from `Haplotype.stats`, whose `RegionStats.min` value is the minimum k-mer depth. The Java field is a signed 32-bit `int`.

Consequently:

- a fetched BAM record is a resolved haplotype record;
- counting fetched or supporting BAM records is haplotype-record counting;
- `XD` is minimum k-mer depth for one resolved haplotype record;
- `XD` is not a sequencing-read count and is not a multiplicity to expand into votes.

### 2.2 Current decision baseline

The available 200-sample simulated cohort and matching adVNTR corpus were loaded explicitly from the local benchmark trees. The golden tests collected 23 tests, skipped none, and passed all 23. An independent display-aware replay applied the emitted-name predicate defined in Section 11.4. That predicate currently agrees with the internal non-`None` name count for every sample, and the oracle asserts the agreement rather than assuming it. The replay measured:

| Policy | Tier | Displayed | Exact | Wrong |
| --- | ---: | ---: | ---: | ---: |
| no BAM rescue | A | 53 | 53 | 0 |
| no BAM rescue | B | 85 | 66 | 19 |
| no BAM rescue | C | 0 | 0 | 0 |
| **no BAM rescue total** |  | **138** | **119** | **19** |
| current record voting | A | 53 | 53 | 0 |
| current record voting | B | 101 | 83 | 18 |
| current record voting | C | 0 | 0 | 0 |
| **current record voting total** |  | **154** | **136** | **18** |

Both policies called zero normal controls. The emitted display agreed with the internal selected name for all measured samples.

Both replays are production-shaped. For each sample they read every non-negative Kestrel row and every kept adVNTR row in file order. `kestrel_vcf` support is the minimum `Estimated_Depth_AlternateVariant` across translated Kestrel rows with the current unknown-propagation rule. Each adVNTR row may translate to multiple events, but its `NumberOfSupportingReads` enters the source minimum once when that row produced at least one translated event, exactly like `_advntr_calls_by_row`; an unparseable/empty row contributes neither a call nor support. Policy 1 stops there and runs `reconcile(*named_vcf, *advntr_calls, supports=supports)`. Policy 2 uses that same preliminary reconciliation to decide BAM-consultation eligibility. If eligible, one `BamRescuer` is opened and each Kestrel row is rescued at that row's own `Motif_fasta` and 1-based `POS_fasta`, preserving an aligned `None` for rows without a consensus. Every named BAM call contributes its unweighted supporting-haplotype-record count to the minimum `kestrel_bam` support. The final call is `reconcile(*named_vcf, *named_bam, *advntr_calls, supports=supports)`, followed by `refine` once for each named BAM call in row order. The displayed-name predicate is defined in Section 11.4.

The current `_hybrid()` golden helper instead reads one locus from `output.bed`; its `test_the_hybrid_total_does_not_regress` docstring describes the resulting 134-exact approximation. That helper, its docstring, and its floor will be replaced by the production-shaped replay above. The 136-exact requirement does not coexist with the obsolete 134 claim.

The production-shaped policy-2 replay marked 83 of 200 samples eligible for BAM consultation. Fifteen carried no Kestrel result row and therefore caused no fetch; the remaining 68 each carried one row, producing exactly 68 per-row BAM fetches. Its exact displayed calls by truth class were: `delGCCCA=3`, `dupA=6`, `dupC=98`, `insCCCC=6`, `insC_pos23=9`, `insG=6`, `insG_pos54=7`, and `insG_pos58=1`, with every other class zero. Relative to policy 1, BAM rescue adds 6 exact `insCCCC`, 4 exact `insG`, and 7 exact `insG_pos54` calls. All 178 non-negative Kestrel records in the corpus carry valid `Motif_fasta` and `POS_fasta`; zero records silently collapse out of the per-row replay. The corpus guard proves the golden row filter selects the same 178 rows as production `_is_negative` before metrics are computed.

The corpus contains 200 Kestrel BAMs and 1,327,794 records with integral `XD`; none was missing or malformed. Observed values ranged from 5 to 8,704, with median 181 and 95th percentile 7,416. These observations establish realistic fixture values, but they do not calibrate a threshold.

These figures are the behavior-preservation oracle for this phase. They were verified from the available data rather than copied from the issue.

## 3. Considered designs

### 3.1 Lexical correction only

Rename comments and report prose while retaining `support`, `total`, `reads`, and an unparsed `XD` internally.

This is the smallest diff, but it leaves the false model in the interfaces that future work will consume. It also cannot prove that `XD` remains observational and cannot express malformed or missing `XD` behavior. Rejected.

### 3.2 Focused semantic containment

Introduce truthful BAM-specific fields and flags, parse `XD` as optional evidence, make low-support wording source-aware, and add decision-invariance tests. Keep source identifiers, numeric comparisons, and algorithms unchanged. Extract only the pure evidence terminology/presentation logic touched inside oversized modules.

This is the selected design. It makes the ontology truthful without introducing a new scoring model.

### 3.3 General evidence framework

Replace the heterogeneous `supports` mapping with a polymorphic or generic evidence hierarchy across Kestrel VCF, Kestrel BAM, and adVNTR.

This would make unit handling structurally stronger, but it is broader than semantic containment and would entangle this phase with calibration and canonicalization work. Rejected for Phase 1.

## 4. Selected data model

### 4.1 `BamConsensus`

`vntyper/scripts/nomenclature_bam.py` will use the following canonical meanings:

```python
from dataclasses import field

@dataclass(frozen=True)
class BamConsensus:
    kind: str
    start: int
    ref_span: int
    inserted: int
    supporting_haplotype_records: int
    fetched_haplotype_records: int
    distinct_edit_count: int
    bases: str = ""
    supporting_record_minimum_kmer_depths: tuple[int | None, ...] = field(
        kw_only=True,
        compare=False,
    )
```

Field meanings:

- `supporting_haplotype_records`: number of fetched resolved haplotype records carrying the winning edit;
- `fetched_haplotype_records`: number of resolved haplotype records fetched for the locus, including records that do not contribute a qualifying edit;
- `distinct_edit_count`: number of distinct qualifying length-changing edit keys seen;
- `supporting_record_minimum_kmer_depths`: one parsed `XD` value per record supporting the winning edit, in encounter order; `None` represents unavailable or invalid `XD` and is retained rather than discarded.

The tuple is observational evidence only. It is not summed, ranked, thresholded, used to choose the winner, or included in dataclass equality/hash comparisons. Its length equals `supporting_haplotype_records` for every consensus produced by `BamRescuer`.

Read-only compatibility properties `support`, `total`, and `n_distinct` will return the corresponding canonical fields. They avoid gratuitous breakage for code reading these internal attributes while making production code and repository tests use the truthful names. Construction using the old keyword names is not preserved: `BamConsensus` is an internal module interface, and retaining false constructor fields would defeat the ontology correction. The new XD tuple is keyword-only, so an old positional construction fails clearly for a missing keyword rather than silently binding `bases` to the XD field.

### 4.2 Parsed XD field

The BAM module will expose a typed pure helper with the semantic name:

```python
def minimum_kmer_depth(record: pysam.AlignedSegment) -> int | None:
    ...
```

It will call `record.get_tag("XD", with_value_type=True)` so validation considers the BAM auxiliary-tag subtype as well as the Python value. The allowed subtypes are the BAM integer subtypes `c`, `C`, `s`, `S`, `i`, and `I`; htslib/pysam may compact an integer written as SAM `i` into any fitting BAM integer subtype. `KeyError` means the tag is absent and returns `None` silently. `TypeError`, `ValueError`, or `OverflowError` while retrieving or validating a present tag returns `None` with a debug diagnostic.

The result is carried separately from the record's vote. No decision function will branch, sort, filter, or compare on XD; functions that receive a `BamConsensus` must ignore its XD evidence when deciding an outcome.

### 4.3 Source-specific support units

The existing `supports: Mapping[str, int | None]` remains in place to contain scope and preserve numeric behavior. Its source-specific meanings will be documented exactly:

| Source key | Unit represented by its integer |
| --- | --- |
| `kestrel_vcf` | Kestrel alternate-allele k-mer-path depth (`Estimated_Depth_AlternateVariant`) |
| `kestrel_bam` | resolved haplotype records carrying the selected BAM edit |
| `advntr` | adVNTR supporting sequencing reads |

The mapping remains heterogeneous and open to future source keys, matching `_is_corroborated`'s current unknown-source behavior. The current minimum comparison remains temporarily in force. That comparison is uncalibrated across unlike units; Phase 1 documents rather than redesigns it. The shared threshold is described as a source-specific evidence-support threshold, never as a read-count threshold.

## 5. XD validation contract

`minimum_kmer_depth` will behave as follows:

| Input state | Typed result | Record vote | Log/error behavior |
| --- | --- | --- | --- |
| integer subtype, value `1..2_147_483_647` | exact `int` | unchanged, one current record contribution | no warning |
| tag absent | `None` | unchanged | absence is allowed; no warning |
| non-integer subtype or non-integral/malformed value | `None` | unchanged | debug diagnostic, no exception from the field parser |
| integer zero | `0` | unchanged | preserved as observed evidence |
| negative integer | `None` | unchanged | debug diagnostic; a minimum count cannot be negative |
| integer `2_147_483_647` | exact `int` | unchanged | no clipping or overflow |
| unsigned integer above `2_147_483_647` | `None` | unchanged | debug diagnostic; impossible from the pinned Java writer |

A corrupt BAM that pysam cannot decode or fetch remains a BAM-I/O failure handled at the existing rescue boundary; the XD helper does not pretend that undecodable bytes are a recoverable field value.

Missing or invalid XD never causes a haplotype record to be dropped. Zero does not suppress a vote. Extreme depth does not create additional votes.

## 6. Decision-preserving data flow

The production flow remains:

```text
Kestrel output.bam
  -> BamRescuer.fetch
  -> CIGAR-derived length-changing edit keys
  -> one existing +1 contribution per record/edit key
  -> highest record count, ties abstain
  -> BamConsensus
  -> from_bam
  -> refine / reconcile
  -> nomenclature columns
  -> TSV, summary JSON, HTML, cohort/report surfaces
```

The only addition beside each record/edit contribution is parsed optional `minimum_kmer_depth`. The vote counter remains an unweighted integer count. The existing stable encounter/order behavior and tie check remain unchanged.

These invariants are binding:

1. Changing only XD values, including adding, removing, corrupting, zeroing, negating, or maximizing them, cannot change the winning edit, fetched-record count, or distinct-edit count.
2. XD cannot break or create a tie.
3. XD cannot change `is_thin`, a low-support flag, Tier A eligibility, Kestrel caller candidate selection, BAM-consultation eligibility (`is_candidate`), `refine`, `reconcile`, or `render`.
4. A high-XD record is still one record contribution; multiple low-XD records still outvote fewer high-XD records.
5. The thin BAM threshold remains `< 3` supporting haplotype records.
6. The Tier A threshold remains `< 5` in each contributing source's current numeric support unit.
7. `kestrel_bam` remains mapped to caller `kestrel`; it never becomes independent corroboration of `kestrel_vcf`.
8. BAM-consultation eligibility and VCF-primary/BAM-refinement rules remain byte-for-byte equivalent in outcome.

## 7. Terminology, flags, and compatibility

### 7.1 Source-aware thin/low-support flags

The closed public flag vocabulary will distinguish evidence units:

- `thin-haplotype-record-support`: Kestrel BAM support is below the unchanged BAM-rescue thinness threshold of 3 records; this is produced by `from_bam` and is not itself a Tier-A blocker;
- `low-haplotype-record-support`: Kestrel BAM support is below the unchanged Tier-A evidence threshold of 5 records; this is produced by `reconcile` and is a Tier-A blocker;
- `low-kmer-path-support`: Kestrel VCF alternate-allele k-mer-path depth is below the unchanged Tier-A evidence threshold of 5;
- `low-read-support`: genuine sequencing-read support, currently adVNTR and the legacy scalar `support=` compatibility input, is below the applicable threshold.
- `low-evidence-support`: a future/unrecognized backing source is below the Tier-A evidence threshold but has no declared unit; this source-neutral fallback preserves the existing flag instead of silently losing the blocker or guessing that its integer counts reads.

The two BAM flags deliberately preserve two existing thresholds that the old `low-read-support` token conflated. `nomenclature_bam.refine` currently uses the `< 3` token as a behavioral guard: a thin BAM consensus cannot settle a caller-disagreement result whose allele is undetermined. That guard will test `thin-haplotype-record-support` after the rename. This consumer is part of decision preservation, not presentation, and receives a direct regression test.

For `supports=`, flag emission retains `reconcile`'s current all-known gate. If any support value relevant to the selected name is `None`, `effective_support` remains `None` and no source-specific low-support flag is added, even when another known source is numerically below 5. Only after all relevant backing sources are known and the unchanged minimum is below 5 are source-specific tokens added for the relevant sources below 5. A source not backing the selected name cannot contribute a low-support flag. If more than one known backing unit is below 5, each truthful token is emitted. The effective numeric support and the Tier-A comparison remain exactly the same current minimum.

For the older scalar `support=` argument, the existing `low-read-support` token remains because the argument's established contract calls it read support and carries no source key. No current production path uses that scalar for Kestrel BAM evidence.

`low-read-support` remains recognized and explained in report formatting so archived summaries and true read-count sources retain their contract. It is not globally renamed. The four new tokens are an acknowledged public vocabulary expansion, and current Kestrel-BAM rows replace the old misleading token with one or both truthful BAM tokens according to the unchanged producer thresholds. This is an intentional terminology/schema correction; it is not a confidence or naming decision change.

The retained legacy token will use these exact source-aware report explanations rather than repeating the old false universal-read claim:

- flag meaning: “Low source support under the emitting version's rule. Current adVNTR or legacy scalar evidence counts sequencing reads; pre-Phase-1 Kestrel BAM rows used this token for resolved haplotype-record support.”
- Tier-A reason: “a retained low-read-support token records support below the emitting version's applicable rule; older Kestrel BAM rows may use it for legacy haplotype-record thinness rather than a sequencing-read count”.

This is intentionally explicit about the historical ambiguity: older Kestrel BAM data used one token for both the `< 3` rescue-thinness producer and the `< 5` Tier-A producer, so archived rows cannot be assigned a single threshold from the token alone.

### 7.2 Config compatibility

The shipped threshold key becomes:

```json
"bam_thin_haplotype_record_support": 3
```

The old `bam_thin_support` key remains accepted as a compatibility fallback for complete custom nomenclature configurations. If both are present, the canonical new key wins. If neither is present, loading raises `KeyError` as the current import-time hard subscript does; no silent default is introduced. The shipped config contains only the canonical key, and comments describe resolved haplotype records and source-specific evidence units. The numeric value does not change.

`nomenclature_evidence.py` will provide the pure resolver `resolve_bam_thin_haplotype_record_support(thresholds: Mapping[str, int]) -> int`. It returns the canonical key when present; otherwise it uses a hard subscript on `bam_thin_support`, which supplies legacy-only behavior and naturally raises `KeyError` when both are absent. It performs no coercion, validation, or defaulting, preserving the current value/type behavior while making all four config cases directly testable without module reloads.

`min_support_for_high_confidence` remains the stable config key to avoid a broader configuration migration. Its comment and code documentation will explicitly state that the value is applied to each source in that source's own evidence unit; it is not a universal read-count threshold.

The module constant becomes `THIN_HAPLOTYPE_RECORD_SUPPORT`; production and repository tests use that canonical name. The old `THIN_SUPPORT` module-global spelling is not a documented/exported interface and is removed rather than retained as another ambiguous production accessor.

### 7.3 Python compatibility

The source key `kestrel_bam`, caller map, nomenclature column names, TSV/JSON shapes, and rendered allele syntax remain unchanged. Read-only `BamConsensus` aliases cover compatibility attribute access, but production call sites in `nomenclature_annotate.py` and the golden replay move from `.support` to `.supporting_haplotype_records`. The new canonical constructor field names are an intentional internal API correction and will be listed in the changelog.

## 8. Focused module boundaries

Two touched production modules already exceed the repository's approximate 650-line limit, so the pure part under change will be extracted rather than growing them:

1. `vntyper/scripts/nomenclature_evidence.py` will own evidence-source unit constants and the pure mapping from a low backing source to its low-support flag, including the `low-evidence-support` fallback for unrecognized sources. `vntyper/scripts/nomenclature.py` will re-export every `FLAG_*` constant declared by the evidence module, with the same value, and will retain `reconcile`; its arithmetic and winner/tier conditions remain unchanged. The closed-vocabulary unit guard will assert this containment relationship, then the existing meanings-coverage guard over `vars(nomenclature)` closes the loop. Constants that remain owned only by `nomenclature.py` do not have to appear in the smaller evidence module.
2. `vntyper/scripts/nomenclature_presentation.py` will own the nomenclature tier wording, flag meanings, Tier-A blocker wording, and the compact Kestrel BAM semantics note. `vntyper/scripts/report_formatting.py` will import/re-export these values so existing imports remain valid.

The extraction is limited to the semantics touched here. It will not move unrelated nomenclature translation, report table formatting, or I/O.

`vntyper/scripts/nomenclature_bam.py` remains below the limit and owns BAM I/O, XD parsing, record/edit voting, and `BamConsensus`. `nomenclature_annotate.py` continues to bind source-specific support to the corresponding source key.

## 9. HTML and public report behavior

The per-sample HTML report will always include one concise static paragraph immediately below its existing `Reading key` heading. This is intentionally outside the `nomenclature_legend` conditional: that conditional suppresses only the tier/flag definition list when no recognized tier or flag is present, not the heading, and it does not prove that BAM rescue ran. The paragraph states the Kestrel artifact contract rather than claiming that the current run used or retained a BAM; it is intentionally present for IGV-off runs and report rerenders where `output.bam` is absent.

The template markup will use literal `<code>output.bam</code>` and `<code>XD</code>` elements with ordinary autoescaped static text; it will not pass generated HTML through `safe`. Normalized visible text must be exactly:

> Kestrel output.bam contains resolved haplotype records, not sequencing reads. Its record counts are haplotype-record support; XD is minimum k-mer depth and does not weight votes or alter names or tiers.

This paragraph is visible, not hover-only, so archived, offline, and printed reports retain it. It appears once rather than being repeated per row.

Report tier and flag explanations will use:

- `haplotype-record support` for `kestrel_bam`;
- `alternate-allele k-mer-path depth` for `kestrel_vcf`;
- `supporting reads` only for adVNTR or other true read-count sources;
- `source-specific evidence support` when describing the common temporary Tier A threshold.

The Kestrel table's `Estimated_Depth_*` help text will agree with the existing scoring documentation and describe k-mer-path/k-mer depth, not reads. The shared `ALT` help becomes the caller-neutral “The alternate allele reported by the caller.” so it is truthful for both Kestrel and adVNTR. The adVNTR `Supporting Reads` heading and explanation remain unchanged because they are legitimate raw-read terminology.

`NOMENCLATURE_FLAG_MEANINGS["allele-unrepresentable-in-vcf"]` will say that the name comes from resolved haplotype records, not “the reads”. This is a public Kestrel-BAM semantic surface even though it does not contain a numeric count.

The embedded IGV panel displays this same Kestrel `output.bam`. Its accessible `aria-label` and any sibling panel prose will say resolved haplotype-record alignments, not read alignments. The visual track behavior does not change.

The cohort HTML tables also carry `Nomenclature_Flags`. `generate_cohort_summary_report` will add context keys named `nomenclature_legend` (the existing helper evaluated over the unformatted Kestrel/adVNTR frames) and `show_kestrel_bam_semantics` (true when either BAM-specific token occurs). `cohort_summary_template.html` will insert an ordinary autoescaped reading-key block immediately after the complete `advntr_missing` conditional block and before the `additional_stats` conditional block. This real anchor is after the existing adVNTR summary plot; no new value uses `safe`. The block defines every present token and conditionally renders the same BAM semantics paragraph. CSV and Excel exports keep the stable `Nomenclature_Flags` column name, ordering, and schema without prose rows; token values change exactly as Section 7.1 describes. The vocabulary is defined in the published nomenclature/output documentation and the cohort HTML artifact.

All TSV, spreadsheet, and summary fields keep their names. Only the truthful source-specific nomenclature flag tokens can change as described in Section 7.

## 10. Error handling

The implementation follows existing conventions:

- missing BAM, missing contig, fetch failure, and unreadable BAM retain their current rescue-abstention behavior;
- missing XD is valid optional evidence and does not log a warning;
- invalid XD is ignored for the observational field, logged at debug level with record context where available, and cannot abort rescue;
- no custom exception type is introduced;
- no default XD value is invented, and `None` is never converted to zero;
- zero is preserved distinctly from missing/invalid evidence.

## 11. Test design

Before each production change, a failing behavior test will be added. Tests will exercise actual pysam BAM records and report rendering, not mocks of the decision logic.

### 11.1 BAM fixtures and XD parser

Synthetic Kestrel BAM fixtures will use haplotype-oriented record names and comments, resolved-haplotype CIGARs, and realistic XD values such as 5, 181, 7,416, and 8,704. Tests will cover:

- valid compact integer subtypes;
- absent XD;
- string/float or otherwise non-integral XD;
- zero;
- negative;
- signed-int maximum;
- unsigned value above the pinned writer's range;
- exact association of optional XD values with records supporting the winning edit;
- truthful `BamConsensus` canonical fields and compatibility properties;
- unchanged fetched-haplotype-record and distinct-edit counts across every XD state.

Config tests will cover canonical-key only, legacy-key only, both keys with the canonical value winning, neither key raising `KeyError`, and the shipped value remaining exactly 3.

### 11.2 Mutation-oriented decision invariants

Each invariant test is chosen to catch a realistic regression:

- two low-XD supporting records beat one extreme-XD opposing record, catching XD weighting;
- equal record counts with unequal XD totals still abstain, catching XD tie-breaking;
- one or two high-XD supporting records remain thin, catching XD substitution for record count;
- the same records with valid, missing, malformed, zero, negative, and extreme XD yield identical winner, support count, fetched-record count, distinct-edit count, flags, tier, and rendered name, catching any XD path into decisions;
- altering XD on a BAM-rescue candidate cannot change BAM-consultation eligibility (`is_candidate`), `refine`, or reconciliation output;
- a record with invalid/missing XD still votes, catching accidental evidence filtering;
- a caller-disagreement result with no allele remains undetermined when the only BAM name has one or two supporting haplotype records, catching any loss of `refine`'s thin-consensus veto during retokenization;
- source-specific low-support tests prove BAM counts cannot emit `low-read-support`, adVNTR reads still can, Kestrel VCF depth uses its k-mer-path token, an unrecognized source uses `low-evidence-support`, dissenting/non-backing sources cannot donate flags, and multiple low backing units are all named truthfully;
- an agreement in which both `kestrel_vcf` and `advntr` back the selected same name but their supports are `{kestrel_vcf: None, advntr: 2}` (and the source-order permutation) emits no low-support token of any kind and remains Tier B, preserving the all-known gate; a paired positive control gives both backing sources known values with one below 5 and requires the corresponding unit-specific token;
- the `< 3` thin-BAM token and `< 5` Tier-A low-BAM token are tested independently, including a four-record BAM that is not thin but is below the Tier-A evidence threshold;
- all currently modified functions receive direct regression coverage.

### 11.3 Surface and HTML tests

Unit tests will render the real report context/template and assert the exact visible explanation from Section 9. They will also assert:

- `XD` is called minimum k-mer depth and explicitly non-decisional;
- Kestrel BAM support is never described as reads;
- Kestrel VCF depth is described as k-mer-path depth;
- adVNTR supporting-read terminology remains present;
- every authoritative `FLAG_*` token is re-exported through `nomenclature` and has a non-empty meaning;
- only low-support tokens emitted by `reconcile` have Tier-A blocker explanations; `thin-haplotype-record-support` has none and cannot produce a false Tier-A reason;
- direct assertions over `NOMENCLATURE_TIERS`, the flag-meaning table, Tier-A blocker table, and column-help entries reject any use of `read`/`reads` for Kestrel or shared source-specific evidence, not merely numeric read-count phrases; this explicitly covers `allele-unrepresentable-in-vcf`, both Kestrel depth entries, and the shared caller-neutral `ALT` help while allowing adVNTR `Supporting Reads`;
- a targeted source-text/identifier contract inspects the full Kestrel-BAM functions in `nomenclature_bam.py`; `annotate_kestrel_frame`, `_open_rescuer`, `_row_haplotype_call`, `_row_verdicts`, and `_haplotype_calls` in `nomenclature_annotate.py`; and support/reconciliation prose in `nomenclature.py`; it rejects the old private names, read-named BAM locals, “row's reads”/“reads per source” prose, and tie-log wording while allowing legitimate file-reading verbs and adVNTR/input-read terminology;
- the real HTML fixtures carry each new token, so the dynamic legend path is non-vacuous;
- the IGV `aria-label` describes resolved haplotype-record alignments, not read alignments;
- the cohort HTML explains every new token it displays, while cohort CSV/Excel schemas remain unchanged;
- the published flag table in `docs/pipeline/nomenclature.md` grows from 8 to all 14 authoritative tokens: the two currently omitted legacy tokens (`known-variant`, `representation-of-caller-call`) plus the four new evidence-unit tokens; a guard requires all 14 and prevents future code/report/docs vocabulary drift;
- legacy `low-read-support` summaries still render intelligibly;
- TSV/JSON/nomenclature column names remain stable.

Because browser-visible report text changes, `make test-browser` is required in addition to unit report tests.

### 11.4 Golden oracle

For this oracle, a **displayed name** means `call.name is not None` and `render(call) == call.name`: the emitted cell contains the selected positional allele name rather than an `allele undetermined` rendering. An exact displayed name equals that sample's truth name; a wrong displayed name is any displayed name that does not. The existing 178 Kestrel VCF calls count is a separate record/caller coverage measure and is not the denominator for these per-sample reconciled display metrics.

`tests/golden/test_nomenclature_golden.py` will replace its single-`output.bed` `_hybrid()` approximation and compute displayed outcomes for both policies using the exact per-row replay defined in Section 2.2:

1. no BAM rescue;
2. current unweighted BAM record voting.

For each policy it will assert:

- external simulation and adVNTR roots were loaded, with 200 mutated and 200 normal samples accounted for;
- total displayed names under the explicit predicate above, exact displayed names, and wrong displayed names;
- the same three counts separately for tiers A, B, and C;
- the sum of tier counts equals the total, preventing omitted lower-tier errors;
- every internal named/displayed relationship follows `render`;
- normal-control findings are zero;
- the exact baseline table in Section 2.2;
- exactly 83 BAM-eligible samples, 68 per-row BAM fetches, all 178 Kestrel records carrying a valid row locus, the exact policy-2 per-class counter in Section 2.2, and the three exact gains relative to policy 1;
- Tier-A wrong is zero without relying on a vacuous empty Tier A.

All four current `_hybrid()` consumers receive an explicit production-shaped disposition:

1. `test_no_tier_a_name_is_wrong_with_bam_rescue_either` remains as a named non-vacuity guard and is also implied by the exact per-tier wrong counts;
2. `test_the_bam_is_consulted_only_for_a_minority_of_samples` is restated as exactly 68 per-row fetches and retains its `< 100` budget, now explicitly counting row fetches rather than one approximate locus per sample;
3. `test_the_hybrid_total_does_not_regress` becomes the exact 154 displayed / 136 exact / 18 wrong policy-2 assertion and removes the obsolete 134 narrative;
4. `test_the_bam_rescue_recovers_alleles_the_vcf_could_not` is retained against the production-shaped per-class counter and strengthened to the exact policy-1 gains listed in Section 2.2.

The older `_reconciled()` policy-1 helper is also replaced by the shared production-shaped replay so its last-record-wins Kestrel support and zero-fallback cannot coexist with the new oracle. Its three consumers remain: wrong Tier-A stays exactly zero, normal-control findings stay exactly zero, and Tier A remains reachable at exactly 53 displayed names.

The oracle will report missing external roots as skips, but verification will treat those skips as unavailable evidence rather than success. The final PR record will distinguish unit-gate completion from external benchmark validation.

In the current verification environment the corpus is reachable through `VNTYPER_SIM_ROOT=/home/bernt-popp/development/vntyper-analyses/results/simulation` and `VNTYPER_ADVNTR_ROOT=/home/bernt-popp/development/vntyper-bench-266-267/advntr`. The final golden command must show 200 mutated samples, 200 normal controls, zero skips, and the exact metric table; a separate XD survey over the same 200 BAMs must show the record count and range stated in Section 2.2.

## 12. Documentation and changelog

The implementation will update all Kestrel-BAM-specific terminology found in:

- code identifiers, local variables, log strings, docstrings, and comments in the BAM, annotation, reconciliation, and presentation flow; specifically `_row_read_call` becomes `_row_haplotype_call`, `_read_calls` becomes `_haplotype_calls`, fetched `reads`/`read` locals become `records`/`record`, and the tie-abstention log reports haplotype records;
- `vntyper/scripts/nomenclature_config.json` comments and the BAM threshold key;
- the HTML nomenclature reading key and report column help;
- `docs/pipeline/nomenclature.md`;
- `docs/pipeline/reports.md`;
- `docs/user-guide/output-files.md`;
- `docs/cli/report.md` for the Kestrel BAM/IGV artifact contract;
- `docs/getting-started/quickstart.md` where `output.bam` is described as reads/alignments;
- `docs/pipeline/kestrel.md` and `docs/pipeline/scoring-and-confidence.md` only where consistency requires it;
- `vntyper/scripts/README.md` and Kestrel-BAM-specific test module/fixture docstrings;
- `tests/unit/test_issue_233_documentation_contract.py`, whose module docstring currently says no published `docs/superpowers` page exists;
- `tests/unit/test_cohort_summary_oracle.py`, whose cohort HTML fingerprint must move because the new legend markup changes the canonical skeleton even when its fixture carries no nomenclature rows; the digest history gains an explicit “Move 9 (#295 — nomenclature evidence reading key)” entry recording that only the intentional autoescaped reading-key block/whitespace changed and that tables, charts, images, filtering, searching, and paging did not;
- the Unreleased section of `docs/about/changelog.md`, including the public flag/config compatibility notes.

Released changelog sections will not be rewritten. The Unreleased entry will explicitly correct the current artifact terminology and, where useful, point back to the older release text it supersedes. Terminology-only consistency edits in `docs/pipeline/scoring-and-confidence.md` may clarify k-mer-path depth, but no equation, threshold narrative, calibration recommendation, or decision is changed.

The user-requested specification and plan pages under `docs/superpowers/` conflict with the repository's ordinary rule that planning artifacts live under `.planning/`. The committed specification already makes `tests/unit/test_coverage_gate.py::test_contributor_docs_match_the_scripts_quality_scope` fail at the branch tip because that test unconditionally forbids `docs/superpowers`. The direct task requirement controls for exactly these two paths:

- `docs/superpowers/specs/2026-08-31-kestrel-bam-evidence-semantics-design.md`, nav title `Kestrel BAM Evidence Semantics Design`;
- `docs/superpowers/plans/2026-08-31-kestrel-bam-evidence-semantics.md`, nav title `Kestrel BAM Evidence Semantics Plan`.

The coverage-gate filesystem allowlist uses the two repo-relative paths above. Its MkDocs navigation assertions use the corresponding docs-relative values, `superpowers/specs/2026-08-31-kestrel-bam-evidence-semantics-design.md` and `superpowers/plans/2026-08-31-kestrel-bam-evidence-semantics.md`, with the stated titles. Every other planning artifact under published `docs/` remains rejected. `AGENTS.md`, `tests/unit/test_issue_233_documentation_contract.py`, and the planning-artifact policy comment at the top of `mkdocs.yml` will be updated in the same task so prose and executable policy agree. `AGENTS.md`'s focused-module inventory will also add `nomenclature_evidence.py` and `nomenclature_presentation.py` and update its count from sixteen to eighteen. This does not reopen `docs/` as a general planning workspace.

Because the MkDocs macros plugin evaluates every published page as Jinja, both required pages must avoid raw Jinja opening delimiters (double-left-brace expressions, percent-brace statements, and hash-brace comments) unless enclosed in an explicit raw block. A dedicated source scan and strict docs build will verify both pages before completion; no plan example may be silently interpolated or parsed as a template comment.

Legitimate read language for input BAM/CRAM/FASTQ processing, adVNTR, SHARK, and sequencing-read recovery is out of scope and remains unchanged.

Root `SPEC.md` was reviewed because it names `output.bam`; its statement only says the artifact is regenerated after a Kestrel rerun and assigns no read/haplotype semantics, so it deliberately remains unchanged.

The benchmark-frequency paragraph in `docs/pipeline/nomenclature.md` will replace its stale “about a fifth of samples” wording with the measured production-shaped result: 83/200 samples are eligible, 15 have no Kestrel row, and 68/200 (34%) produce one BAM row fetch. This is terminology/measurement documentation of unchanged eligibility, not a rescue-policy change.

## 13. Explicit non-goals

This phase does not include:

- molecular identity/canonicalization and selection of which passing Kestrel candidate is reported from #270;
- adVNTR artifact-flag reliability policy from #267;
- threshold calibration or a calibration command from #269;
- XD weighting, ranking, filtering, or thresholding;
- new dominance, abstention, or tie rules;
- threshold value changes;
- caller-agreement, tier, emitted-name, Kestrel caller candidate-selection, BAM-consultation-eligibility, or rescue-policy changes;
- Kestrel/adVNTR workflow changes;
- broad nomenclature or report refactoring;
- version bump, release, tag, publication, or merge;
- a claim to close #295.

The intended follow-up order remains #270, then #267, then joint #295/#269 calibration.

## 14. Acceptance criteria

The phase is acceptable only if:

1. all production and public Kestrel `output.bam` paths call its records resolved haplotype records and its count haplotype-record support;
2. XD is parsed into typed optional minimum-k-mer-depth evidence with every state in Section 5 tested;
3. no XD value can affect voting, ties, thin support, tier, name, Kestrel caller candidate selection, or BAM-consultation eligibility;
4. the BAM-specific thin/low-support tokens are emitted instead of `low-read-support` for BAM counts, `refine` retains its thin-consensus veto, the all-known gate remains intact, and genuine read-count usage remains;
5. numeric thresholds and every baseline decision metric remain exact;
6. synthetic BAM fixtures model Kestrel haplotypes and realistic XD;
7. the HTML report contains the concise visible explanation in Section 9;
8. all affected source, configuration, report, documentation, and changelog terminology is consistent and non-clinical;
9. the two directly required `docs/superpowers/` pages are published and navigated without weakening the rejection of any other planning artifact under `docs/`;
10. targeted, full unit, browser, golden-with-loaded-corpus, formatting, lint, typing, integration-compatibility, patch-coverage, docs, published-page delimiter scans, and `check-all` verification is recorded;
11. fresh Opus 5 adversarial reviews find no unresolved Critical or Important issue in the specification, plan, per-task diffs, or final diff.

## 15. First adversarial specification review disposition

A fresh restricted Claude Code review verified `canonicalModel: claude-opus-5` and made no edits. Its first verdict was `FAIL / BLOCKED`. This revision resolves its findings as follows; a fresh scoped re-review is still required before owner approval.

| ID | Severity | Disposition in this revision |
| --- | --- | --- |
| C1 | Critical | Names `refine` as a consumer of the renamed thin token and requires a direct veto-regression test. |
| C2 | Critical | Records the already-red docs policy test and adds exact-path test, `AGENTS.md`, and MkDocs-comment updates. |
| C3 | Critical | Retains `reconcile`'s all-known gate before emitting any per-source low-support token. |
| I1 | Important | Splits `< 3` thin BAM support from `< 5` Tier-A low BAM support and their report meanings. |
| I2 | Important | Requires the full authoritative flag vocabulary to remain re-exported and guarded across both modules. |
| I3 | Important | Adds the IGV accessible label and sibling panel prose to the report surfaces/tests. |
| I4 | Important | Adds a dynamic cohort HTML legend and explicitly preserves CSV/Excel schemas with published definitions. |
| I5 | Important | Makes neither config key individually mandatory while preserving `KeyError` when both are absent, and specifies all compatibility cases/tests. |
| I6 | Important | Adds the MkDocs macros delimiter constraint and verification scan for both published pages. |
| I7 | Important | Defines displayed/exact/wrong names as explicit predicates over `call.name` and `render`. |
| I8 | Important | Specifies literal template markup, exact normalized text, unconditional sample-report placement, and cohort trigger. |
| I9 | Important | Adds direct vocabulary-table assertions and non-vacuous fixtures in addition to rendered-HTML checks. |
| M1–M4 | Minor | Excludes XD from equality, names parser exceptions, covers unchanged counts, and makes the new tuple keyword-only. |
| M5–M6 | Minor | Disambiguates both candidate-selection concepts and adds the omitted CLI/README/test-docstring surfaces. |

## 16. First scoped re-review disposition

A second fresh restricted review again verified `canonicalModel: claude-opus-5`. It confirmed that C1–C3 and I1–I9 were resolved, but returned `FAIL / BLOCKED` for five remaining Important ambiguities plus six Minor findings. This revision resolves them as follows and requires another fresh scoped re-review:

| ID | Severity | Disposition in this revision |
| --- | --- | --- |
| F1 | Important | Defines both production-shaped replay policies, per-row locus/order/support flow, and replacement of the obsolete single-`output.bed` 134-exact helper. |
| F2 | Important | Names the exact required plan path and MkDocs navigation title. |
| F3 | Important | Adds private identifiers, rescue locals, and tie-log wording to the production terminology scope plus a targeted source contract. |
| F4 | Important | Requires both unknown/low sources to back the selected same name and adds a fully known positive control. |
| F5 | Important | Adds the truthful source-neutral `low-evidence-support` fallback and its behavior test for unrecognized sources. |
| F6 | Minor | Clarifies that export schemas stay stable while token values change under Section 7.1. |
| F7 | Minor | Corrects the first-review I5 disposition to preserve `KeyError` when both keys are absent. |
| F8 | Minor | Defines evidence-module-to-nomenclature flag containment rather than impossible whole-module equality. |
| F9 | Minor | Gives exact source-aware legacy-token meaning and Tier-A-reason text, including historical threshold ambiguity. |
| F10 | Minor | Names cohort context keys, insertion point, and ordinary autoescaping. |
| F11 | Minor | Adds the stale documentation-contract docstring and explicitly records why root `SPEC.md` needs no semantic edit. |

## 17. Second scoped re-review disposition

A third fresh restricted review verified `canonicalModel: claude-opus-5`. It confirmed the preceding blockers resolved and returned two remaining Important test-scope gaps plus six Minor ambiguities. This revision resolves them and requires a final scoped specification re-review:

| ID | Severity | Disposition in this revision |
| --- | --- | --- |
| G1 | Important | Gives all four old `_hybrid()` consumers explicit replacements and pins 68 per-row fetches plus exact per-class gains. |
| G2 | Important | Widens terminology tests to reconciliation/BAM function prose and the nonnumeric `allele-unrepresentable-in-vcf` report meaning. |
| G3 | Minor | States that the display-aware predicate currently agrees with named counts and is asserted rather than assumed. |
| G4 | Minor | Anchors the cohort legend after the adVNTR missing-status block and before the adVNTR plot heading. |
| G5 | Minor | Separates repo-relative filesystem allowlist paths from docs-relative MkDocs nav values. |
| G6 | Minor | Adds a published nomenclature flag-table drift guard. |
| G7 | Minor | Adds both focused modules to the `AGENTS.md` inventory and updates sixteen to eighteen. |
| G8 | Minor | Makes shared `ALT` help caller-neutral and explicitly tests it. |

## 18. Final scoped re-review disposition

A fourth fresh restricted review verified `canonicalModel: claude-opus-5`. It confirmed the decision and XD invariants but returned four remaining Important surface/test omissions plus four Minor details. This revision resolves them and is subject to one final fresh verdict:

| ID | Severity | Disposition in this revision |
| --- | --- | --- |
| H1 | Important | Anchors cohort markup after the real `advntr_missing` block and before `additional_stats`. |
| H2 | Important | Adds `annotate_kestrel_frame` and `_open_rescuer` to the source terminology contract. |
| H3 | Important | Names the cohort fingerprint oracle and requires a reasoned Move 9 digest history entry. |
| H4 | Important | Authorizes/document-tests all 14 flags, including the two legacy omissions and four new tokens. |
| H5 | Minor | Adds `NOMENCLATURE_TIERS` to direct terminology assertions. |
| H6 | Minor | Defines a pure hard-subscript config resolver so all compatibility cases are testable without reloads. |
| H7 | Minor | Replaces `_reconciled()` with the shared production-shaped policy-1 replay and retains all consumers. |
| H8 | Minor | Replaces stale “about a fifth” prose with measured 83 eligible / 68 fetched sample loci. |
