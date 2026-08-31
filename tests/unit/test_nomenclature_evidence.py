"""Source-aware evidence vocabulary and threshold compatibility."""

import pytest

from vntyper.scripts.nomenclature import load_nomenclature_config
from vntyper.scripts.nomenclature_evidence import (
    FLAG_LOW_EVIDENCE_SUPPORT,
    FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
    FLAG_LOW_KMER_PATH_SUPPORT,
    FLAG_LOW_READ_SUPPORT,
    low_support_flag_for_source,
    resolve_bam_thin_haplotype_record_support,
)

pytestmark = pytest.mark.unit


def test_canonical_bam_thin_threshold_is_used() -> None:
    thresholds = {"bam_thin_haplotype_record_support": 3}

    assert resolve_bam_thin_haplotype_record_support(thresholds) == 3


def test_legacy_bam_thin_threshold_remains_accepted() -> None:
    thresholds = {"bam_thin_support": 3}

    assert resolve_bam_thin_haplotype_record_support(thresholds) == 3


def test_canonical_bam_thin_threshold_wins_over_legacy_value() -> None:
    thresholds = {"bam_thin_haplotype_record_support": 3, "bam_thin_support": 99}

    assert resolve_bam_thin_haplotype_record_support(thresholds) == 3


def test_missing_bam_thin_threshold_has_no_silent_default() -> None:
    with pytest.raises(KeyError, match="bam_thin_support"):
        resolve_bam_thin_haplotype_record_support({})


def test_shipped_bam_thin_threshold_uses_only_the_canonical_key() -> None:
    thresholds = load_nomenclature_config()["thresholds"]

    assert thresholds["bam_thin_haplotype_record_support"] == 3
    assert "bam_thin_support" not in thresholds


@pytest.mark.parametrize(
    ("source", "expected_flag"),
    [
        ("kestrel_bam", FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT),
        ("kestrel_vcf", FLAG_LOW_KMER_PATH_SUPPORT),
        ("advntr", FLAG_LOW_READ_SUPPORT),
        ("future_caller", FLAG_LOW_EVIDENCE_SUPPORT),
    ],
)
def test_low_support_flag_uses_the_source_evidence_unit(source: str, expected_flag: str) -> None:
    assert low_support_flag_for_source(source) == expected_flag
