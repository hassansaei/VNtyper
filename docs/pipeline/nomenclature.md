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
| **A** | Two independent callers agree after normalisation, the motif context matches the canonical unit, and each source's evidence support meets the threshold | the name, e.g. `59dupC` |
| **B** | A name was computed, but something above is missing | the same name, carrying the tier and the flags that say what is missing |
| **C** | No allele could be determined | `frameshift +1, allele undetermined` — no position at all |

Tier B **does** show its name. Withholding it was the original design and it was measured
to be a bad trade: of 200 benchmark samples, 129 had the correct name computed and only 46
were allowed to display it. Suppressing a name that was right 83 times, to avoid showing
one that was never wrong, discards information a reader can weigh for themselves. The tier
and the flags travel beside the name and say how far it has been checked; they do not
decide whether the reader may see it.

Evidence support is **unknown-hostile**: a source whose support value is blank or
non-numeric makes the whole agreement's support unknown, and unknown never clears the
threshold. One caller's evidence is not evidence about another's.

**No single caller can reach tier A on its own.** On the benchmark, Kestrel places the
whole `insG` family one position 3′ of truth; those records look perfectly clean in
isolation. Only a second, independent source separates them.

"Independent" means *a different caller*, not merely a different file. Kestrel's VCF
and Kestrel's resolved haplotype alignment are two artifacts from one caller, so their
agreement is one opinion rather than two, and it does not promote anything.

Tier C is the point of the whole design. Where a caller's allele is genuinely
indistinguishable from another, printing a name would be a confident falsehood, and
`allele undetermined` is the honest answer.

## Flags

Stable, kebab-case, and matchable:

The table contains all 14 authoritative tokens. Its support terms are source-specific:
Kestrel VCF support is alternate-allele k-mer-path depth, Kestrel `output.bam` support is
resolved haplotype-record support, and adVNTR support counts sequencing reads. The
unchanged `min_support_for_high_confidence` value (5) is applied in each source's own
unit; it is not a universal read-count threshold.

| Flag | Meaning |
|---|---|
| `position-ambiguous` | the variant can shift within a tract; see `Ambiguity_Interval` |
| `spans-unit-junction` | the span crosses the boundary between two repeat units |
| `motif-context-diverges` | the assigned motif differs from the canonical unit where the name lands |
| `allele-unrepresentable-in-vcf` | Kestrel's VCF cannot express the allele shape; the name comes from its resolved haplotype records |
| `thin-haplotype-record-support` | Kestrel BAM haplotype-record support is below the unchanged thinness threshold |
| `low-haplotype-record-support` | Kestrel BAM haplotype-record support is below the corroborated tier's source-specific threshold |
| `low-kmer-path-support` | Kestrel alternate-allele k-mer-path depth is below the corroborated tier's source-specific threshold |
| `low-read-support` | current adVNTR or legacy scalar sequencing-read support is low; archived pre-Phase-1 Kestrel BAM rows can also carry this token |
| `low-evidence-support` | support from a source whose evidence unit is undeclared is below the source-specific threshold |
| `caller-disagreement` | Kestrel and adVNTR do not name the same allele |
| `length-truncated` | the reported length is a lower bound rather than the full extent |
| `sequence-undetermined` | a length is known but the inserted or deleted bases are not |
| `known-variant` | the name matches a MUC1 variant described in the literature; the table checks a name and never produces one |
| `representation-of-caller-call` | the name represents the caller's report but matches no described variant; it requires validation |

The shipped configuration calls the BAM thinness key
`bam_thin_haplotype_record_support` (value 3). Complete custom configurations may still
use the former `bam_thin_support` key as a compatibility fallback; the canonical key wins
when both are present, and omitting both remains an error. The broader
`min_support_for_high_confidence` name stays stable for configuration compatibility even
though its value is interpreted in the evidence unit of each source.

## The two companion fields

`Ambiguity_Interval` states the span in which **every anchor is the same allele**. It is
what makes `59dupC`, `53dupC` and the older `27dupC` recognisable as one event rather
than three.

`Repeat_Form` states what was actually measured — `53C[7]>53C[8]`, the tract went from
seven copies to eight — instead of implying we know *which* base was added. It also
scales: an `insCCCC` reads `53C[11]`, considerably clearer than `56_59dupCCCC`.

Both are empty where they do not apply.

## Where the resolved haplotype records are consulted

Kestrel's VCF can express 1-vs-1, 1-vs-N and N-vs-1 records and nothing else, so a
delins has no representation in it. Where a call cannot be resolved from the VCF,
VNtyper walks the CIGARs in `output.bam` (the resolved haplotype alignment the report's
IGV track already shows) and merges adjacent non-matching blocks — which is how a `1X1I`
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
*different events* is a conflict, not a gap, so a thin haplotype-record consensus is not
allowed to settle it with a number; `allele undetermined` stands. And because an ordinary
single-caller call is tier B and therefore not a candidate, a run **without** the optional
adVNTR module rarely opens the BAM at all — so a delins that Kestrel's VCF could not
express is unlikely to be recovered unless adVNTR also ran.

## When two sources outvote a third

There is one case where the Kestrel VCF does not have the last word: when **two
sources spanning two different callers** name the same allele and the VCF names
another. adVNTR and Kestrel's resolved haplotype records agreeing is genuine
corroboration across two callers, and it outvotes Kestrel's VCF placement.

The independence requirement is doing the work here. Kestrel's VCF agreeing with
Kestrel's own alignment is *not* two sources, so it never outvotes adVNTR and never
promotes a tier. Without that distinction, one caller corroborating itself would look
exactly like the two independent sources tier A asks for.

Measured over the 200-sample benchmark this recovers 6 samples and loses none —
`insG` goes from 1 to 5 and `insG_pos54` from 0 to 2, while `dupA` and `dupC` are
unmoved. The narrowness matters: an earlier attempt that simply preferred
the BAM-derived name whenever it disagreed cost `dupA` 6 correct calls out of 10 and `insCCCC`
6 out of 10. Requiring a second *caller* to agree is what makes the difference.

## Both callers' names are kept

Where the two callers disagree, collapsing them into one verdict destroys the evidence
a reader needs in order to weigh that verdict. Both result files therefore carry
`Nomenclature_Kestrel` and `Nomenclature_adVNTR` alongside the reconciled
`Nomenclature`, so either file on its own says what each caller reported and where
they parted company.

The columns are named for the callers rather than relatively ("the other caller"), so
a row means the same thing in `kestrel_result.tsv`, in `output_adVNTR_result.tsv` and
in a cohort table that merges the two.

adVNTR is an optional module. When it has not run, `Nomenclature_adVNTR` is empty and
the Kestrel result stands exactly as its own stage wrote it.

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
