# Kestrel BAM Evidence Semantics — Phase 1 Design

Date: 2026-08-31

Issue: [#295](https://github.com/hassansaei/VNtyper/issues/295)

Status: proposed for explicit approval after adversarial review

## 1. Purpose

This phase corrects the evidence ontology used for Kestrel's `output.bam` without changing the decisions built on that evidence. The BAM contains Kestrel-resolved haplotype records, not sequencing reads. Its `XD` auxiliary tag is the minimum k-mer depth of the resolved haplotype, not a read count.

The implementation will make those meanings explicit in production code, typed interfaces, configuration, tests, reports, and documentation while preserving the present:

- edit voting and tie abstention;
- VCF-primary/BAM-refinement policy;
- candidate selection;
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

The available 200-sample simulated cohort and matching adVNTR corpus were loaded explicitly from the local benchmark trees. The golden tests collected 23 tests, skipped none, and passed all 23. An independent replay of the displayed result, rather than the internal `name` field alone, measured:

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
@dataclass(frozen=True)
class BamConsensus:
    kind: str
    start: int
    ref_span: int
    inserted: int
    supporting_haplotype_records: int
    fetched_haplotype_records: int
    distinct_edit_count: int
    supporting_record_minimum_kmer_depths: tuple[int | None, ...]
    bases: str = ""
```

Field meanings:

- `supporting_haplotype_records`: number of fetched resolved haplotype records carrying the winning edit;
- `fetched_haplotype_records`: number of resolved haplotype records fetched for the locus, including records that do not contribute a qualifying edit;
- `distinct_edit_count`: number of distinct qualifying length-changing edit keys seen;
- `supporting_record_minimum_kmer_depths`: one parsed `XD` value per record supporting the winning edit, in encounter order; `None` represents unavailable or invalid `XD` and is retained rather than discarded.

The tuple is observational evidence only. It is not summed, ranked, thresholded, or used to choose the winner. Its length equals `supporting_haplotype_records`.

Read-only compatibility properties `support`, `total`, and `n_distinct` will return the corresponding canonical fields. They avoid gratuitous breakage for code reading these internal attributes while making production code and repository tests use the truthful names. Construction using the old keyword names is not preserved: `BamConsensus` is an internal module interface, and retaining false constructor fields would defeat the ontology correction.

### 4.2 Parsed XD field

The BAM module will expose a typed pure helper with the semantic name:

```python
def minimum_kmer_depth(record: pysam.AlignedSegment) -> int | None:
    ...
```

It will call `record.get_tag("XD", with_value_type=True)` so validation considers the BAM auxiliary-tag subtype as well as the Python value. The allowed subtypes are the BAM integer subtypes `c`, `C`, `s`, `S`, `i`, and `I`; htslib/pysam may compact an integer written as SAM `i` into any fitting BAM integer subtype.

The result is carried separately from the record's vote. No decision function will branch, sort, filter, or compare on XD; functions that receive a `BamConsensus` must ignore its XD evidence when deciding an outcome.

### 4.3 Source-specific support units

The existing `supports: Mapping[str, int | None]` remains in place to contain scope and preserve numeric behavior. Its source-specific meanings will be documented exactly:

| Source key | Unit represented by its integer |
| --- | --- |
| `kestrel_vcf` | Kestrel alternate-allele k-mer-path depth (`Estimated_Depth_AlternateVariant`) |
| `kestrel_bam` | resolved haplotype records carrying the selected BAM edit |
| `advntr` | adVNTR supporting sequencing reads |

The mapping remains heterogeneous and the current minimum comparison remains temporarily in force. That comparison is uncalibrated across unlike units; Phase 1 documents rather than redesigns it. The shared threshold is described as a source-specific evidence-support threshold, never as a read-count threshold.

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

1. Changing only XD values, including adding, removing, corrupting, zeroing, negating, or maximizing them, cannot change the winning edit.
2. XD cannot break or create a tie.
3. XD cannot change `is_thin`, a low-support flag, Tier A eligibility, the selected candidate, `refine`, `reconcile`, or `render`.
4. A high-XD record is still one record contribution; multiple low-XD records still outvote fewer high-XD records.
5. The thin BAM threshold remains `< 3` supporting haplotype records.
6. The Tier A threshold remains `< 5` in each contributing source's current numeric support unit.
7. `kestrel_bam` remains mapped to caller `kestrel`; it never becomes independent corroboration of `kestrel_vcf`.
8. BAM rescue candidates and VCF-primary/BAM-refinement rules remain byte-for-byte equivalent in outcome.

## 7. Terminology, flags, and compatibility

### 7.1 Source-aware low-support flags

The closed public flag vocabulary will distinguish evidence units:

- `low-haplotype-record-support`: Kestrel BAM support is below the applicable threshold;
- `low-kmer-path-support`: Kestrel VCF alternate-allele k-mer-path depth is below the applicable threshold;
- `low-read-support`: genuine sequencing-read support, currently adVNTR and the legacy scalar `support=` compatibility input, is below the applicable threshold.

For `supports=`, each known backing source below the unchanged threshold contributes its unit-specific flag. A source not backing the selected name cannot contribute a low-support flag, preserving the current relevance rule. If more than one backing source is low, each truthful flag is emitted. The effective numeric support remains the same current minimum, so tier and name decisions do not change.

For the older scalar `support=` argument, the existing `low-read-support` token remains because the argument's established contract calls it read support and carries no source key. No current production path uses that scalar for Kestrel BAM evidence.

`low-read-support` remains recognized and explained in report formatting so archived summaries and true read-count sources retain their contract. It is not globally renamed. The two new tokens are an acknowledged public vocabulary expansion, and current Kestrel-BAM rows may replace the old misleading token with `low-haplotype-record-support`. This is an intentional terminology/schema correction; it is not a confidence or naming decision change.

### 7.2 Config compatibility

The shipped threshold key becomes:

```json
"bam_thin_haplotype_record_support": 3
```

The old `bam_thin_support` key remains accepted as a compatibility fallback for complete custom nomenclature configurations. If both are present, the canonical new key wins. The shipped config contains only the canonical key, and comments describe resolved haplotype records and source-specific evidence units. The numeric value does not change.

`min_support_for_high_confidence` remains the stable config key to avoid a broader configuration migration. Its comment and code documentation will explicitly state that the value is applied to each source in that source's own evidence unit; it is not a universal read-count threshold.

### 7.3 Python compatibility

The source key `kestrel_bam`, caller map, nomenclature column names, TSV/JSON shapes, and rendered allele syntax remain unchanged. Read-only `BamConsensus` aliases cover common attribute access. The new canonical constructor field names are an intentional internal API correction and will be listed in the changelog.

## 8. Focused module boundaries

Two touched production modules already exceed the repository's approximate 650-line limit, so the pure part under change will be extracted rather than growing them:

1. `vntyper/scripts/nomenclature_evidence.py` will own evidence-source unit constants and the pure mapping from a low backing source to its low-support flag. `vntyper/scripts/nomenclature.py` will re-export flag constants needed by existing importers and will retain `reconcile`; its arithmetic and winner/tier conditions remain unchanged.
2. `vntyper/scripts/nomenclature_presentation.py` will own the nomenclature tier wording, flag meanings, Tier-A blocker wording, and the compact Kestrel BAM semantics note. `vntyper/scripts/report_formatting.py` will import/re-export these values so existing imports remain valid.

The extraction is limited to the semantics touched here. It will not move unrelated nomenclature translation, report table formatting, or I/O.

`vntyper/scripts/nomenclature_bam.py` remains below the limit and owns BAM I/O, XD parsing, record/edit voting, and `BamConsensus`. `nomenclature_annotate.py` continues to bind source-specific support to the corresponding source key.

## 9. HTML and public report behavior

The HTML report's visible nomenclature reading key will include this concise explanation whenever nomenclature results are present:

> Kestrel `output.bam` contains resolved haplotype records, not sequencing reads. Its record counts are haplotype-record support; `XD` is minimum k-mer depth and does not weight votes or alter names or tiers.

This sentence is visible text, not hover-only text, so archived, offline, and printed reports retain the explanation. It will appear once near the existing tier/flag reading key rather than being repeated per row. The existing conditional that suppresses the reading key for an entirely unnamed run remains; no empty report gains irrelevant prose.

Report tier and flag explanations will use:

- `haplotype-record support` for `kestrel_bam`;
- `alternate-allele k-mer-path depth` for `kestrel_vcf`;
- `supporting reads` only for adVNTR or other true read-count sources;
- `source-specific evidence support` when describing the common temporary Tier A threshold.

The Kestrel table's `Estimated_Depth_*` help text will agree with the existing scoring documentation and describe k-mer-path/k-mer depth, not reads. The adVNTR `Supporting Reads` heading and explanation remain unchanged because they are legitimate raw-read terminology.

All TSV and summary fields keep their names. Only the truthful source-specific nomenclature flag token can change as described in Section 7.

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
- truthful `BamConsensus` canonical fields and compatibility properties.

### 11.2 Mutation-oriented decision invariants

Each invariant test is chosen to catch a realistic regression:

- two low-XD supporting records beat one extreme-XD opposing record, catching XD weighting;
- equal record counts with unequal XD totals still abstain, catching XD tie-breaking;
- one or two high-XD supporting records remain thin, catching XD substitution for record count;
- the same records with valid, missing, malformed, zero, negative, and extreme XD yield identical winner, support count, flags, tier, and rendered name, catching any XD path into decisions;
- altering XD on a candidate cannot change `is_candidate`, `refine`, or reconciliation output;
- a record with invalid/missing XD still votes, catching accidental evidence filtering;
- source-specific low-support tests prove BAM counts cannot emit `low-read-support`, adVNTR reads still can, Kestrel VCF depth uses its k-mer-path token, dissenting/non-backing sources cannot donate flags, and multiple low backing units are all named truthfully;
- all currently modified functions receive direct regression coverage.

### 11.3 Surface and HTML tests

Unit tests will render the real report context/template and assert the exact visible explanation from Section 9. They will also assert:

- `XD` is called minimum k-mer depth and explicitly non-decisional;
- Kestrel BAM support is never described as reads;
- Kestrel VCF depth is described as k-mer-path depth;
- adVNTR supporting-read terminology remains present;
- all new public flag tokens have non-empty meanings and Tier-A blocker explanations;
- legacy `low-read-support` summaries still render intelligibly;
- TSV/JSON/nomenclature column names remain stable.

Because browser-visible report text changes, `make test-browser` is required in addition to unit report tests.

### 11.4 Golden oracle

`tests/golden/test_nomenclature_golden.py` will compute displayed outcomes for both policies:

1. no BAM rescue;
2. current unweighted BAM record voting.

For each policy it will assert:

- external simulation and adVNTR roots were loaded, with 200 mutated and 200 normal samples accounted for;
- total displayed names, exact displayed names, and wrong displayed names;
- the same three counts separately for tiers A, B, and C;
- the sum of tier counts equals the total, preventing omitted lower-tier errors;
- every internal named/displayed relationship follows `render`;
- normal-control findings are zero;
- the exact baseline table in Section 2.2;
- Tier-A wrong is zero without relying on a vacuous empty Tier A.

The oracle will report missing external roots as skips, but verification will treat those skips as unavailable evidence rather than success. The final PR record will distinguish unit-gate completion from external benchmark validation.

## 12. Documentation and changelog

The implementation will update all Kestrel-BAM-specific terminology found in:

- code docstrings and comments in the BAM, annotation, reconciliation, and presentation flow;
- `vntyper/scripts/nomenclature_config.json` comments and the BAM threshold key;
- the HTML nomenclature reading key and report column help;
- `docs/pipeline/nomenclature.md`;
- `docs/pipeline/reports.md`;
- `docs/user-guide/output-files.md`;
- `docs/getting-started/quickstart.md` where `output.bam` is described as reads/alignments;
- `docs/pipeline/kestrel.md` and `docs/pipeline/scoring-and-confidence.md` only where consistency requires it;
- the Unreleased section of `docs/about/changelog.md`, including the public flag/config compatibility notes.

Historical changelog entries will not be rewritten except where an existing statement falsely describes the current `output.bam` contract and would mislead readers about the artifact. Any such correction will be explicitly marked as a terminology correction rather than silently changing release history.

The user-requested specification and plan pages under `docs/superpowers/` conflict with the repository's ordinary rule that planning artifacts live under `.planning/`. The direct task requirement controls for these two files only. `tests/unit/test_coverage_gate.py` will replace its blanket `docs/superpowers` prohibition with an exact allowlist for this specification and its matching implementation plan, assert that both are present in `mkdocs.yml`, and continue rejecting every other planning artifact under published `docs/`. This does not reopen `docs/` as a general planning workspace.

Legitimate read language for input BAM/CRAM/FASTQ processing, adVNTR, SHARK, and sequencing-read recovery is out of scope and remains unchanged.

## 13. Explicit non-goals

This phase does not include:

- molecular identity or canonicalization from #270;
- adVNTR artifact-flag reliability policy from #267;
- threshold calibration or a calibration command from #269;
- XD weighting, ranking, filtering, or thresholding;
- new dominance, abstention, or tie rules;
- threshold value changes;
- caller-agreement, tier, emitted-name, candidate-selection, or rescue-policy changes;
- Kestrel/adVNTR workflow changes;
- broad nomenclature or report refactoring;
- version bump, release, tag, publication, or merge;
- a claim to close #295.

The intended follow-up order remains #270, then #267, then joint #295/#269 calibration.

## 14. Acceptance criteria

The phase is acceptable only if:

1. all production and public Kestrel `output.bam` paths call its records resolved haplotype records and its count haplotype-record support;
2. XD is parsed into typed optional minimum-k-mer-depth evidence with every state in Section 5 tested;
3. no XD value can affect voting, ties, thin support, tier, name, or candidate selection;
4. the new BAM-specific low-support token is emitted instead of `low-read-support` for BAM counts while genuine read-count usage remains;
5. numeric thresholds and every baseline decision metric remain exact;
6. synthetic BAM fixtures model Kestrel haplotypes and realistic XD;
7. the HTML report contains the concise visible explanation in Section 9;
8. all affected source, configuration, report, documentation, and changelog terminology is consistent and non-clinical;
9. the two directly required `docs/superpowers/` pages are published and navigated without weakening the rejection of any other planning artifact under `docs/`;
10. targeted, full unit, browser, golden-with-loaded-corpus, formatting, lint, typing, integration-compatibility, patch-coverage, docs, and `check-all` verification is recorded;
11. fresh Opus 5 adversarial reviews find no unresolved Critical or Important issue in the specification, plan, per-task diffs, or final diff.
