"""Pure candidate-policy and reference-coverage tests."""

from __future__ import annotations

from collections.abc import Iterable

import pytest

from vntyper.scripts.reference_resolution import ordered_reference_candidates, uncovered_reference_contigs

pytestmark = pytest.mark.unit


def test_default_candidate_order_maps_all_explicit_reference_sources() -> None:
    """Removing or reordering a default source would change the documented probe policy."""
    config = {
        "reference_data": {
            "cram_reference_hg38": "/refs/cram-hg38.fa",
            "bwa_reference_hg38": "/refs/bwa-hg38.fa",
        }
    }

    candidates = ordered_reference_candidates(config, "hg38", "/refs/cli.fa")

    assert candidates == (
        ("cli", "/refs/cli.fa"),
        ("config_cram_reference", "/refs/cram-hg38.fa"),
        ("config_bwa_reference", "/refs/bwa-hg38.fa"),
    )


@pytest.mark.parametrize(
    ("assembly", "cram_reference", "bwa_reference"),
    [
        ("GRCh37", "/refs/cram37.fa", "/refs/bwa37.fa"),
        ("GRCh38", "/refs/cram38.fa", "/refs/bwa38.fa"),
        ("hg19", "/refs/cram37.fa", "/refs/bwa37.fa"),
        ("hg19_ensembl", "/refs/cram37.fa", "/refs/bwa37.fa"),
        ("hg19_ncbi", "/refs/cram37.fa", "/refs/bwa37.fa"),
        ("hg38", "/refs/cram38.fa", "/refs/bwa38.fa"),
        ("hg38_ensembl", "/refs/cram38.fa", "/refs/bwa38.fa"),
        ("hg38_ncbi", "/refs/cram38.fa", "/refs/bwa38.fa"),
    ],
)
def test_every_supported_assembly_label_resolves_the_configured_coordinate_family_references(
    assembly: str,
    cram_reference: str,
    bwa_reference: str,
) -> None:
    """Accepted aliases cannot silently discard both configured reference candidates."""
    config = {
        "reference_data": {
            "cram_reference_hg19": "/refs/cram37.fa",
            "bwa_reference_hg19": "/refs/bwa37.fa",
            "cram_reference_hg38": "/refs/cram38.fa",
            "bwa_reference_hg38": "/refs/bwa38.fa",
        }
    }

    candidates = ordered_reference_candidates(config, assembly, None)

    assert candidates[1:] == (
        ("config_cram_reference", cram_reference),
        ("config_bwa_reference", bwa_reference),
    )


def test_label_specific_reference_keys_override_the_family_fallback_including_null() -> None:
    """A replacement config can deliberately specialize or disable one accepted label."""
    config = {
        "reference_data": {
            "cram_reference_hg19": "/refs/family-cram.fa",
            "bwa_reference_hg19": "/refs/family-bwa.fa",
            "cram_reference_hg19_ncbi": None,
            "bwa_reference_hg19_ncbi": "/refs/ncbi-bwa.fa",
        }
    }

    candidates = ordered_reference_candidates(config, "hg19_ncbi", None)

    assert candidates[1:] == (
        ("config_cram_reference", None),
        ("config_bwa_reference", "/refs/ncbi-bwa.fa"),
    )


def test_configured_candidate_order_is_preserved_without_silent_completion() -> None:
    """A valid explicit subset remains verbatim rather than being sorted or filled in."""
    config = {
        "cram": {
            "reference_candidate_order": [
                "config_bwa_reference",
                "cli",
                "htslib_resolved",
            ]
        },
        "reference_data": {"bwa_reference_hg19": "/refs/bwa-hg19.fa"},
    }

    candidates = ordered_reference_candidates(config, "hg19", None)

    assert candidates == (("config_bwa_reference", "/refs/bwa-hg19.fa"), ("cli", None))


@pytest.mark.parametrize(
    ("order", "message"),
    [
        (("cli", "htslib_resolved"), "must be a list"),
        (["cli", "cli", "htslib_resolved"], "duplicate explicit"),
        (["none", "htslib_resolved"], "unknown reference candidate source: none"),
        (["header_ur", "htslib_resolved"], "unknown reference candidate source: header_ur"),
        (["htslib", "htslib_resolved"], "unknown reference candidate source: htslib"),
        ([object(), "htslib_resolved"], "unknown reference candidate source"),
        (["htslib_resolved", "cli"], "must end with terminal htslib_resolved"),
        (["cli"], "exactly one terminal htslib_resolved"),
        (["cli", "htslib_resolved", "htslib_resolved"], "exactly one terminal htslib_resolved"),
    ],
)
def test_invalid_candidate_orders_are_rejected_verbatim(order: object, message: str) -> None:
    """Malformed policy cannot be silently reordered, deduplicated, or interpreted as legacy aliases."""
    with pytest.raises(ValueError, match=message):
        ordered_reference_candidates({"cram": {"reference_candidate_order": order}}, "hg19", None)


@pytest.mark.parametrize(
    ("header_contigs", "reference_contigs", "expected"),
    [
        (("chr1", "chr2"), None, None),
        (("chr1", "chr2"), {"chr1", "chr2", "chr3"}, ()),
        (("chr2", "chr1", "chrM"), {"chr1"}, ("chr2", "chrM")),
    ],
)
def test_reference_coverage_distinguishes_unknown_complete_and_uncovered(
    header_contigs: Iterable[str],
    reference_contigs: set[str] | None,
    expected: tuple[str, ...] | None,
) -> None:
    """Missing FAI evidence remains distinct from a known complete reference."""
    assert uncovered_reference_contigs(header_contigs, reference_contigs) == expected


class TestResolveFromMapping:
    """Membership, not truthiness - the rule `configured_reference_candidates` already uses.

    A key present with value None is a deliberate "disabled". Falling through it would
    silently re-enable a reference an operator switched off.
    """

    def test_the_physical_key_wins_when_present(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        mapping = {"bwa_reference_hg38_ensembl": "/refs/ensembl.fa", "bwa_reference_hg38": "/refs/ucsc.fa"}
        resolved = resolve_from_mapping("bwa", "hg38_ensembl", mapping)
        assert resolved.key == "bwa_reference_hg38_ensembl"
        assert resolved.value == "/refs/ensembl.fa"
        assert resolved.is_fallback is False

    def test_the_family_key_is_used_and_flagged_when_the_physical_key_is_absent(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        resolved = resolve_from_mapping("bwa", "hg38_ensembl", {"bwa_reference_hg38": "/refs/ucsc.fa"})
        assert resolved.key == "bwa_reference_hg38"
        assert resolved.value == "/refs/ucsc.fa"
        assert resolved.is_fallback is True

    def test_a_present_null_is_authoritative_and_does_not_fall_through(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        mapping = {"bwa_reference_hg38_ensembl": None, "bwa_reference_hg38": "/refs/ucsc.fa"}
        resolved = resolve_from_mapping("bwa", "hg38_ensembl", mapping)
        assert resolved.key == "bwa_reference_hg38_ensembl"
        assert resolved.value is None
        assert resolved.is_fallback is False

    def test_no_key_present_at_all_returns_none(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        assert resolve_from_mapping("bwa", "hg38_ensembl", {}) is None

    def test_both_ncbi_labels_resolve_the_same_entry(self):
        from vntyper.scripts.reference_resolution import resolve_from_mapping

        mapping = {"bwa_reference_GRCh38": "/refs/ncbi.fna"}
        for label in ("GRCh38", "hg38_ncbi"):
            resolved = resolve_from_mapping("bwa", label, mapping)
            assert resolved.value == "/refs/ncbi.fna", label
            assert resolved.is_fallback is False, label
