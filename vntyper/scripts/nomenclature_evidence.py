"""Source-aware evidence vocabulary for MUC1 nomenclature decisions."""

from __future__ import annotations

from collections.abc import Mapping

FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT = "thin-haplotype-record-support"
FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT = "low-haplotype-record-support"
FLAG_LOW_KMER_PATH_SUPPORT = "low-kmer-path-support"
FLAG_LOW_READ_SUPPORT = "low-read-support"
FLAG_LOW_EVIDENCE_SUPPORT = "low-evidence-support"

_LOW_SUPPORT_FLAG_BY_SOURCE = {
    "kestrel_bam": FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
    "kestrel_vcf": FLAG_LOW_KMER_PATH_SUPPORT,
    "advntr": FLAG_LOW_READ_SUPPORT,
}


def low_support_flag_for_source(source: str) -> str:
    """Return the low-support flag matching a source's evidence unit.

    Args:
        source: Evidence source key from a nomenclature call.

    Returns:
        str: A source-specific flag, or the unit-neutral fallback for an
        unrecognized source.
    """
    return _LOW_SUPPORT_FLAG_BY_SOURCE.get(source, FLAG_LOW_EVIDENCE_SUPPORT)


def resolve_bam_thin_haplotype_record_support(thresholds: Mapping[str, int]) -> int:
    """Resolve the canonical BAM thin-support threshold with legacy fallback.

    Args:
        thresholds: Nomenclature threshold configuration.

    Returns:
        int: The configured haplotype-record threshold.

    Raises:
        KeyError: If neither the canonical nor legacy key is configured.
    """
    if "bam_thin_haplotype_record_support" in thresholds:
        return thresholds["bam_thin_haplotype_record_support"]
    return thresholds["bam_thin_support"]
