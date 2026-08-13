# MUC1 Nomenclature

!!! warning "Research use only"
    VNtyper is a research tool. Nothing on this page is a diagnostic statement.

Kestrel and adVNTR each report MUC1-VNTR variants in their own internal coordinate
frame — `POS 67 G>GG`, `I22_2_G_LEN1` — and neither emits the naming the MUC1
literature uses. VNtyper translates both into a single literature-compatible name,
with an explicit statement of how much confidence that name carries.

## The name

Names are the **bare positional shorthand** on a named reference:

| Example | Meaning |
|---|---|
| `59dupC` | the canonical MUC1 duplication |
| `60dupA` | duplication of the `A` at position 60 |
| `58_59insG` | a `G` inserted between positions 58 and 59 |
| `1_5delGCCCA` | a 5 bp deletion |
| `54_56delinsAT` | 3 bases replaced by `AT` |

The reference is **the canonical MUC1 60 bp repeat unit in coding orientation**, whose
tract of 7 × C at positions 53–59 is the one Wenzel et al. (2018, PMID:29520014)
publish.

### Why there is no `c.` prefix

A `c.` prefix asserts a coding-DNA reference sequence, and **no transcript places this
tract at positions 53–59**. Emitting `c.59dupC` would claim a coordinate system that
does not exist for this locus. A real HGVS prefix would require depositing the 60 bp
unit as an accessioned reference; until then the reference is named in prose rather
than implied by a prefix that would be wrong.

### The strand trap

HGVS's "3′-most" rule means 3′ of the **coding** sequence. MUC1 is transcribed from the
minus strand, so HGVS-3′ is genomic-**left** here. A pipeline that runs `bcftools norm`
and assumes the result is HGVS-ready happens to be right for MUC1 and would be wrong by
the full repeat width on a plus-strand gene. VNtyper reverse-complements into the coding
frame before normalising, so the direction is explicit rather than incidental.

## Confidence tiers

The tier is an **emission rule**, not a label: it decides what may be printed.

| Tier | Condition | What is emitted |
|---|---|---|
| **A** | Two independent callers agree after normalisation, the motif context matches the canonical unit, and read support meets the threshold | the bare name, e.g. `59dupC` |
| **B** | A name was computed, but something above is missing | the event and its ambiguity window — **never a bare number** |
| **C** | No allele could be determined | `frameshift +1, allele undetermined` — no position at all |

**No single caller can reach tier A on its own.** On the benchmark, Kestrel places the
whole `insG` family one position 3′ of truth; those records look perfectly clean in
isolation. Only a second, independent source separates them.

Tier C is the point of the whole design. Where a caller's allele is genuinely
indistinguishable from another, printing a name would be a confident falsehood, and
`allele undetermined` is the honest answer.

## Flags

Stable, kebab-case, and matchable:

| Flag | Meaning |
|---|---|
| `position-ambiguous` | the variant can shift within a tract; see `Ambiguity_Interval` |
| `spans-unit-junction` | the span crosses the boundary between two repeat units |
| `motif-context-diverges` | the assigned motif differs from the canonical unit where the name lands |
| `allele-unrepresentable-in-vcf` | the reads show something the caller's VCF cannot express |
| `sequence-undetermined` | a length is known but the inserted bases are not |
| `low-read-support` | too few reads under the call to promote it |
| `caller-disagreement` | the sources do not agree |
| `length-truncated` | the reported length is shorter than the evidence |

## The two companion fields

`Ambiguity_Interval` states the span in which **every anchor is the same allele**. It is
what makes `59dupC`, `53dupC` and the older `27dupC` recognisable as one event rather
than three.

`Repeat_Form` states what was actually measured — `53C[7]>53C[8]`, the tract went from
seven copies to eight — instead of implying we know *which* base was added. It also
scales: an `insCCCC` reads `53C[11]`, considerably clearer than `56_59dupCCCC`.

Both are empty where they do not apply.

## Where the reads are consulted

Kestrel's VCF can express 1-vs-1, 1-vs-N and N-vs-1 records and nothing else, so a
delins has no representation in it. Where a call cannot be resolved from the VCF,
VNtyper walks the CIGARs in `output.bam` (the alignment the report's IGV track already
shows) and merges adjacent non-matching blocks — which is how a `1X1I` block becomes
the delins the VCF had no way to write down.

The policy is **VCF primary, BAM refines**. The reads may supply an allele the VCF
lacked, and may override a shape the VCF structurally cannot hold, but they may not
veto a name the VCF already has: the alignment splits reads across many reference
records, so a locus often carries only one to three reads, and a thin consensus
overruling a well-supported record loses more than it gains.

The BAM is opened only for calls that need it — on the benchmark, about a fifth of
samples — and never at all for the rest.

## A known limitation

Every name is anchored on the canonical unit, even where the motif the caller
assigned differs from it there. That projection is what makes the canonical
duplication resolve correctly on the motifs whose own C-tract is shorter or sits
elsewhere — anchoring on the assigned motif instead loses roughly a quarter of the
canonical calls.

The cost is confined to tier B. A projected call always carries
`motif-context-diverges`, which blocks promotion, so no *confident* name rests on the
projection. What the projection can get wrong is the tier-B **event word**: an edit
that is a duplication against the assigned motif may read as an insertion against the
canonical one. The position and the ambiguity window are unaffected.

## What this does and does not fix

Normalisation fixes **description-level** disagreement: two callers describing one
allele in two coordinate systems. It cannot fix **allele-level** disagreement, where a
caller's record encodes a sequence that was never present. On the simulated benchmark a
substantial minority of calls are allele-level disagreements — `ins25bp` is reported as
a 1 bp event, and the `insG` family lands one position 3′ of truth. That is a finding
about the callers, not a gap in the naming. Tier C is the response.

Naming *which copy* in the repeat array carries the variant is out of scope: Kestrel is
alignment-free and cannot know it, and resolving it needs long reads.
