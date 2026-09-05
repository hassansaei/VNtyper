# MUC1 Nomenclature

!!! warning "Research use only"
    VNtyper is a research tool. Output nomenclature is for investigative use and does not constitute a clinical diagnosis.

Kestrel and adVNTR each report MUC1 VNTR variants in their own internal coordinate frames (for example, `POS 67 G>GG` or `I22_2_G_LEN1`), and neither emits standard literature naming. VNtyper translates both into standard literature-compatible shorthand with an explicit confidence tier.

## The name

Names use bare positional shorthand anchored to a designated reference sequence:

| Example | Meaning |
|---|---|
| `59dupC` | Canonical MUC1 cytosine duplication |
| `60dupA` | Duplication of the `A` at position 60 |
| `58_59insG` | Insertion of `G` between positions 58 and 59 |
| `1_5delGCCCA` | 5 bp deletion |
| `54_56delinsAT` | 3 bases substituted with `AT` |

The reference is the canonical MUC1 60 bp repeat unit in coding orientation, whose poly-C tract of 7 cytosines at positions 53-59 matches Wenzel et al. (2018, PMID:29520014).

## Molecular identity and caller representation

Raw coordinates (`POS`, `REF`, `ALT`, `Variant`, `Nomenclature`) represent internal caller outputs rather than canonical molecular identities. Every positive Kestrel and adVNTR call records a four-field identity record in this order:

1. `Molecular_Identity`: The serialized canonical edit, empty when unresolved.
2. `Molecular_Identity_Status`: Categorized as `unique`, `legacy-selected-among-multiple`, or `unresolved`.
3. `Equivalent_Representation_Count`: Count of caller representations matching the selected identity (0 when unresolved).
4. `Identity_Hypothesis_Count`: Count of distinct resolved candidate identities evaluated for the caller.

Sample and cohort HTML, TSV, CSV, and JSON outputs propagate these four recorded values. Summary rows missing any member display the literal `legacy identity not recorded` in all four cells.

### Why there is no `c.` prefix

A `c.` prefix specifies a coding-DNA reference sequence, and no official transcript places the repeat unit poly-C tract at positions 53-59. Emitting `c.59dupC` would falsely imply reference to an official transcript coordinate system.

### Strand orientation

HGVS rules specify 3-prime positioning relative to the coding sequence. Because MUC1 is transcribed from the negative strand, coding 3-prime corresponds to genomic-left. Running `bcftools norm` directly on genomic coordinates shifts indels in the wrong direction for MUC1. VNtyper reverse-complements variants into coding orientation before applying left-alignment and normalisation.

## Confidence tiers

Tiers govern emission rules for variant reporting:

| Tier | Criteria | Output format |
|---|---|---|
| **A** | Two independent callers agree after normalisation, motif context matches the canonical unit, and each source meets its support threshold | Bare name, e.g. `59dupC` |
| **B** | A variant name is resolved, but one or more Tier A conditions are unmet | Variant name annotated with Tier B and disqualifying flags |
| **C** | Allele sequence cannot be resolved | `frameshift +1, allele undetermined` (no position at all) |

Tier B calls display the candidate name alongside explanatory flags. In benchmark testing of 200 samples, 129 had correct names computed, but only 46 met Tier A criteria. Suppressing accurate names discards actionable investigative evidence; the tier and flags explicitly document confidence.

Support evaluation is strict: missing or non-numeric support values fail the threshold. Evidence from one caller cannot substitute for another.

**No single caller can reach Tier A independently.** Kestrel systematically places the `insG` family one base 3-prime of truth; only an independent caller resolves this ambiguity. Kestrel VCF and BAM records originate from the same algorithm and count as one caller.

adVNTR calls matching any of the 24 active governed recurrent states are marked `identity-insufficient`. While visible in caller outputs, they cannot corroborate Tier A calls or outvote alternative identities.

Tier C reports unresolvable alleles as `frameshift +1, allele undetermined` to prevent publishing speculative names.

## Flags

Stable, kebab-case, and matchable:

The table contains all 14 authoritative tokens. Its support terms are source-specific:
Kestrel VCF support is alternate-allele k-mer-path depth, Kestrel `output.bam` support is
resolved haplotype-record support, and adVNTR support counts sequencing reads. The
unchanged `min_support_for_high_confidence` value (5) is applied in each source's own
unit; it is not a universal read-count threshold.

| Flag | Meaning |
|---|---|
| `position-ambiguous` | The edit can be written at more than one position in the repeat unit. The ambiguity interval gives the range over which the placements are indistinguishable. |
| `spans-unit-junction` | The event crosses the boundary between two repeat units, so it cannot be expressed as one span on a single unit. |
| `motif-context-diverges` | The sequence around the call in the motif the caller assigned differs from the canonical unit where the name lands, so the coordinate projected onto the canonical unit is less certain. |
| `allele-unrepresentable-in-vcf` | The allele cannot be written in Kestrel's VCF shape. The name comes from Kestrel's resolved haplotype records, which preserve the full allele shape. |
| `thin-haplotype-record-support` | Resolved Kestrel BAM haplotype-record support is below the unchanged thinness threshold. |
| `low-haplotype-record-support` | Resolved Kestrel BAM haplotype-record support is below the corroborated tier's source-specific threshold. |
| `low-kmer-path-support` | Kestrel alternate-allele k-mer-path depth is below the corroborated tier's source-specific threshold. |
| `low-read-support` | Low source support under the emitting version's rule. Current adVNTR or legacy scalar sequencing-read support is low; archived pre-Phase-1 Kestrel BAM rows can also carry this token. |
| `low-evidence-support` | Support from a source whose evidence unit is undeclared is below the corroborated tier's source-specific threshold. |
| `caller-disagreement` | Kestrel and adVNTR did not name the same allele. |
| `length-truncated` | The reported length is a lower bound rather than the full extent. |
| `sequence-undetermined` | The inserted or deleted sequence itself could not be determined, so no position is given. |
| `known-variant` | The name matches a MUC1 variant described in the literature. The table is used to check a name, never to produce one. |
| `representation-of-caller-call` | The name represents what the caller reported and matches no described variant. It requires validation. |

The shipped configuration calls the BAM thinness key `bam_thin_haplotype_record_support` (value 3). Complete custom configurations may still supply the former `bam_thin_support` key as a compatibility fallback; the canonical key wins when both are present, and omitting both remains an error. The broader `min_support_for_high_confidence` name stays stable for configuration compatibility even though its value is interpreted in the evidence unit of each source.

## The two companion fields

`Ambiguity_Interval` states the span in which every anchor represents the same biological edit, unifying designations such as `59dupC`, `53dupC`, and legacy `27dupC`.

`Repeat_Form` states what was actually measured: `53C[7]>53C[8]`, indicating repeat tract expansion from 7 copies to 8, instead of implying base-specific placement.

Both fields remain empty when not applicable.

## Where the resolved haplotype records are consulted

Kestrel's VCF can express 1-vs-1, 1-vs-N and N-vs-1 records and nothing else, so a
delins has no representation in it. Where a call cannot be resolved from the VCF,
VNtyper walks the CIGARs in `output.bam` (the resolved haplotype alignment the report's
IGV track already shows) and merges adjacent non-matching blocks: this is how a `1X1I`
block becomes the delins the VCF had no way to write down. Each BAM record represents a
Kestrel-resolved haplotype, not a sequencing read.

The policy is **VCF primary, BAM refines**. The resolved haplotype records may supply an
allele the VCF lacked, and may override a shape the VCF structurally cannot hold, but
they may not veto a name the VCF already has. Low haplotype-record support remains the
same thin-consensus check that was previously described as a read-count check.

In the 200-sample benchmark, 83/200 samples are eligible for BAM consultation. Of those,
15 have no Kestrel result row, while 68/200 (34%) produce one BAM row fetch; the BAM is
not fetched for the other samples. These are measured calls through the unchanged
eligibility policy, not a new threshold.

The BAM's optional `XD` tag is the minimum k-mer depth of that resolved haplotype. It is
recorded separately from haplotype-record support and does not weight votes or alter
names or tiers. Integer values from 1 through 2,147,483,647 are retained exactly, and
zero is retained as zero. Missing or malformed values are unavailable; negative values
and unsigned integers above 2,147,483,647 are unavailable. Every resolved haplotype
record still contributes one unweighted vote in all of these cases.

Two consequences worth stating plainly. A locus where the two callers describe
different events is a conflict, not a gap, so a thin haplotype-record consensus is not
allowed to settle it with a number; `allele undetermined` stands. And because an ordinary
single-caller call is tier B and therefore not a candidate, a run without the optional
adVNTR module rarely opens the BAM at all: a delins that Kestrel's VCF could not
express is unlikely to be recovered unless adVNTR also ran.

## When two sources outvote a third

When two sources spanning two different callers name the same allele and Kestrel's VCF names another, the two agreeing sources outvote the VCF call. Specifically, agreement between adVNTR and Kestrel's resolved haplotype records represents genuine cross-caller corroboration that outvotes the Kestrel VCF call.

The independence requirement is essential: Kestrel's VCF agreeing with Kestrel's own BAM alignment represents one caller, so it never outvotes adVNTR and never promotes a tier.

On the 200-sample benchmark, this rule resolves 6 samples without introducing errors: `insG` increases from 1 to 5 calls, and `insG_pos54` increases from 0 to 2 calls, while `dupA` and `dupC` remain unchanged.

## Both callers' names are kept

When callers disagree, collapsing them into one verdict destroys necessary diagnostic evidence. Both result files therefore record `Nomenclature_Kestrel` and `Nomenclature_adVNTR` alongside the reconciled `Nomenclature` string.

adVNTR is an optional module. When it has not run, `Nomenclature_adVNTR` is empty and the Kestrel result stands as written.

Runs snapshot governing evidence at `provenance/advntr_artifact_evidence.json` and log its full SHA-256 digest in `advntr_evidence_digest`.

Custom decision profiles can be selected via `--decision-profile`. Fixed safety thresholds remain constant across all profiles: BAM flank window 8, thin haplotype record support 3, Kestrel Tier A alternate k-mer-path depth 5, adVNTR Tier A sequencing read support 5.

## A known limitation

Nomenclature names are projected onto the canonical unit. When the motif assigned by the caller diverges from the canonical unit, the projection may misstate the event word in Tier B calls (for example, reading a duplication as an insertion). Such variants carry `motif-context-diverges`, preventing promotion to Tier A. Positions and ambiguity intervals remain unaffected.

## What this does and does not fix

Normalisation fixes description-level discrepancies where callers describe one allele in different coordinate systems. It cannot correct allele-level errors where a caller misidentifies the underlying sequence. Discrepant sequence calls receive Tier C (`allele undetermined`). Assigning variants to specific repeat unit copies within the array requires long-read sequencing.
