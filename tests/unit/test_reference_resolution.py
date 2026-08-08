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
