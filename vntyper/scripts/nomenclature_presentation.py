"""Source-aware presentation text for MUC1 nomenclature evidence.

The nomenclature engine owns decisions and stable tokens. This module owns the
human-readable meanings shown in reports, without filesystem or dataframe I/O.
"""

from __future__ import annotations

from collections.abc import Sequence

from vntyper.scripts import nomenclature

__all__ = [
    "COLUMN_HELP",
    "KESTREL_BAM_SEMANTICS",
    "NOMENCLATURE_FLAG_MEANINGS",
    "NOMENCLATURE_TIERS",
    "TIER_A_BLOCKERS",
    "tier_reason",
]

#: Exact visible explanation shared by sample and cohort reports.
KESTREL_BAM_SEMANTICS = (
    "Kestrel output.bam contains resolved haplotype records, not sequencing reads. Its record counts are "
    "haplotype-record support; XD is minimum k-mer depth and does not weight votes or alter names or tiers."
)

#: What a nomenclature tier means, keyed by the emitted tier letter.
NOMENCLATURE_TIERS: dict[str, tuple[str, str]] = {
    "A": (
        "corroborated name",
        "Two independent callers name the same allele, that name matches a MUC1 variant described in the "
        "literature, support from each backing source meets the source-specific evidence threshold, and no "
        "context or disagreement flag is set. Tier A is the only tier that states a bare position.",
    ),
    "B": (
        "qualified name",
        "A name was computed, but at least one tier-A condition is unmet - a single caller, an allele nobody "
        "has described before, thin or low source-specific evidence support, a diverging motif context, or the "
        "two callers disagreeing. The name is shown so it can be weighed; the flags beside it say which "
        "condition was missed.",
    ),
    "C": (
        "no position",
        "No coordinate could be computed. What is stated is the event and its net length change, without a "
        "position - not a negative, and not a name withheld.",
    ),
}

#: Flag-evidenced reasons a tier-B result did not reach Tier A.
TIER_A_BLOCKERS: dict[str, str] = {
    nomenclature.FLAG_MOTIF_CONTEXT_DIVERGES: (
        "the motif context around the call diverges from the canonical X repeat unit"
    ),
    nomenclature.FLAG_SEQUENCE_UNDETERMINED: "the inserted or deleted sequence could not be determined",
    nomenclature.FLAG_CALLER_DISAGREEMENT: "the two callers did not name the same allele",
    nomenclature.FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT: (
        "resolved Kestrel haplotype-record support is below the corroborated tier's source-specific evidence threshold"
    ),
    nomenclature.FLAG_LOW_KMER_PATH_SUPPORT: (
        "Kestrel alternate-allele k-mer-path depth is below the corroborated tier's source-specific evidence threshold"
    ),
    nomenclature.FLAG_LOW_READ_SUPPORT: (
        "a retained low-read-support token records support below the emitting version's applicable rule; older "
        "Kestrel BAM rows may use it for legacy haplotype-record thinness rather than a sequencing-read count"
    ),
    nomenclature.FLAG_LOW_EVIDENCE_SUPPORT: (
        "support from a source whose evidence unit is not declared is below the corroborated tier's "
        "source-specific evidence threshold"
    ),
    nomenclature.FLAG_REPRESENTATION_ONLY: "no MUC1 variant described in the literature matches this name",
}

#: Closed nomenclature flag vocabulary in report wording.
NOMENCLATURE_FLAG_MEANINGS: dict[str, str] = {
    nomenclature.FLAG_POSITION_AMBIGUOUS: (
        "The edit can be written at more than one position in the repeat unit. The ambiguity interval gives the "
        "range over which the placements are indistinguishable."
    ),
    nomenclature.FLAG_SPANS_UNIT_JUNCTION: (
        "The event crosses the boundary between two repeat units, so it cannot be expressed as one span on a "
        "single unit."
    ),
    nomenclature.FLAG_MOTIF_CONTEXT_DIVERGES: (
        "The sequence around the call in the motif the caller assigned differs from the canonical X unit, so "
        "the coordinate projected onto X is less certain."
    ),
    nomenclature.FLAG_ALLELE_UNREPRESENTABLE: (
        "The allele cannot be written in Kestrel's VCF shape. The name comes from Kestrel's resolved haplotype "
        "records, which preserve the full allele shape."
    ),
    nomenclature.FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT: (
        "Resolved Kestrel haplotype-record support is below the BAM rescue thinness threshold."
    ),
    nomenclature.FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT: (
        "Resolved Kestrel haplotype-record support is below the corroborated tier's source-specific evidence threshold."
    ),
    nomenclature.FLAG_LOW_KMER_PATH_SUPPORT: (
        "Kestrel alternate-allele k-mer-path depth is below the corroborated tier's source-specific evidence threshold."
    ),
    nomenclature.FLAG_LOW_READ_SUPPORT: (
        "Low source support under the emitting version's rule. Current adVNTR or legacy scalar evidence counts "
        "sequencing reads; pre-Phase-1 Kestrel BAM rows used this token for resolved haplotype-record support."
    ),
    nomenclature.FLAG_LOW_EVIDENCE_SUPPORT: (
        "Support from a source whose evidence unit is not declared is below the corroborated tier's "
        "source-specific evidence threshold."
    ),
    nomenclature.FLAG_CALLER_DISAGREEMENT: "Kestrel and adVNTR did not name the same allele.",
    nomenclature.FLAG_LENGTH_TRUNCATED: "The reported length is a lower bound rather than the full extent.",
    nomenclature.FLAG_SEQUENCE_UNDETERMINED: (
        "The inserted or deleted sequence itself could not be determined, so no position is given."
    ),
    nomenclature.FLAG_KNOWN_VARIANT: (
        "The name matches a MUC1 variant described in the literature. The table is used to check a name, never "
        "to produce one - a match says somebody has described this allele before, which is weaker than saying "
        "the call is correct."
    ),
    nomenclature.FLAG_REPRESENTATION_ONLY: (
        "The name represents what the caller reported and matches no described variant. It requires validation."
    ),
}

#: Results-table heading help, shown on hover and in the reading key.
COLUMN_HELP: dict[str, str] = {
    "Supporting Reads": "Reads adVNTR counted in support of this call.",
    "Mean Coverage": "Mean read depth adVNTR measured over the locus. Not the region mean above.",
    "P-value": "adVNTR's significance for the call, to three significant figures.",
    "Repeat Unit": "The adVNTR repeat unit the call was made in.",
    "Motif": "The MUC1 repeat motif this variant was annotated onto.",
    "Motifs": "The raw left-right motif pair Kestrel emitted for the record.",
    "Variant": "The event class the caller reported.",
    "Position": "Position within the 120 bp Kestrel motif pair, not a genomic coordinate.",
    "POS": "Position within the caller's own coordinate frame, not a genomic coordinate.",
    "REF": "The reference allele at the reported position.",
    "ALT": "The alternate allele reported by the caller.",
    "Depth (Variant)": "Kestrel alternate-allele k-mer-path depth.",
    "Depth (Region)": "Kestrel total k-mer depth across the active region.",
    "Depth Score": (
        "Variant depth over region depth. It scales with the inverse of the array length, so it is comparable "
        "within a sample and not between assemblies."
    ),
    "Confidence": "Kestrel's own calibrated label for the call.",
    "Flag": "Whether a configured flagging rule fired on this row, and which one.",
    "Evidence Disposition": (
        "Whether this adVNTR State may contribute to molecular-identity agreement. "
        "Identity-insufficient rows remain visible findings."
    ),
    "Evidence_Disposition": (
        "Whether this adVNTR State may contribute to molecular-identity agreement. "
        "Identity-insufficient rows remain visible findings."
    ),
    "MUC1 Name": "The reconciled MUC1 name for the allele, on the canonical X repeat unit.",
    "Nomenclature": "The reconciled MUC1 name for the allele, on the canonical X repeat unit.",
    "Tier": "How far the name beside it has been checked. See the reading key below the table.",
    "Nomenclature_Tier": "How far the name beside it has been checked. See the reading key below the table.",
    "Flags": "Closed-vocabulary qualifiers on the name. Each one present is spelled out below the table.",
    "Nomenclature_Flags": (
        "Closed-vocabulary qualifiers on the name. Each one present is spelled out below the table."
    ),
    "Kestrel Name": "The name derived from Kestrel's record alone, before reconciliation.",
    "Nomenclature_Kestrel": "The name derived from Kestrel's record alone, before reconciliation.",
    "adVNTR Name": "The name derived from adVNTR's record alone, before reconciliation.",
    "Nomenclature_adVNTR": "The name derived from adVNTR's record alone, before reconciliation.",
    "Ambiguity": "The range of positions over which this edit is indistinguishable.",
    "Ambiguity_Interval": "The range of positions over which this edit is indistinguishable.",
    "Repeat Form": "The affected tract before and after the edit, written as a repeat count.",
    "Repeat_Form": "The affected tract before and after the edit, written as a repeat count.",
    "Naming Note": "What the name is, and how far it has been checked.",
    "Nomenclature_Note": "What the name is, and how far it has been checked.",
    "Motif Sequence": "The 60 bp repeat-unit half named by the Motif column.",
    "VID": "adVNTR's own variant identifier.",
    "NumberOfSupportingReads": "Reads adVNTR counted in support of this call.",
    "MeanCoverage": "Mean read depth adVNTR measured over the locus.",
    "Pvalue": "adVNTR's significance for the call, to three significant figures.",
    "RU": "The adVNTR repeat unit the call was made in.",
}


def tier_reason(tier: str, flags: Sequence[str]) -> str:
    """Say which flag-evidenced Tier-A condition a call missed.

    Args:
        tier: Tier letter carried by the call.
        flags: Individual nomenclature flag tokens.

    Returns:
        str: A reason clause for a blocked Tier-B call, otherwise an empty string.
    """
    if tier != "B":
        return ""
    reasons = [TIER_A_BLOCKERS[flag] for flag in flags if flag in TIER_A_BLOCKERS]
    if not reasons:
        return ""
    if len(reasons) == 1:
        return f"Held below the corroborated tier because {reasons[0]}."
    return "Held below the corroborated tier because " + "; and ".join(reasons) + "."
