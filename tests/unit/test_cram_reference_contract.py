"""Unit tests for the pure reference-compressed CRAM fixture contract."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from scripts.cram_reference_contract import (
    HG19_PRIMARY_CONTIGS,
    REGISTERED_B178_INDEXED_REGION_DIGEST,
    REGISTERED_B178_INDEXED_REGION_RECORDS,
    LossyConversionError,
    header_with_hg19_m5,
    normalize_sam_record,
    validate_registered_b178_index_evidence,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
PROVENANCE = REPO_ROOT / "scripts/ucsc_hg19_primary_contigs.tsv"


def test_normalize_sam_record_sorts_only_optional_fields() -> None:
    record = "read-1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tAS:i:4\n"

    assert normalize_sam_record(record) == (b"read-1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:4\tNM:i:0\n")


def test_hg19_header_adds_exact_primary_contig_m5_tags() -> None:
    header = "@HD\tVN:1.3\tSO:coordinate\n@SQ\tSN:chr1\tLN:249250621\n@SQ\tLN:243199373\tSN:chr2\n@RG\tID:sample\n"

    assert header_with_hg19_m5(header) == (
        "@HD\tVN:1.3\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:249250621\tM5:1b22b98cdeb4a9304cb5d48026a85128\n"
        "@SQ\tLN:243199373\tSN:chr2\tM5:a0d9851da00400dec1098a9255ac712e\n"
        "@RG\tID:sample\n"
    )


@pytest.mark.parametrize(
    ("header", "message"),
    [
        ("@HD\tVN:1.6\n", "contains no @SQ"),
        ("@SQ\tSN:chrUn\tLN:1\n", "unsupported hg19 sequence"),
        ("@SQ\tSN:chr1\tLN:1\n", "length mismatch for chr1"),
        (
            "@SQ\tSN:chr1\tLN:249250621\tM5:00000000000000000000000000000000\n",
            "M5 mismatch for chr1",
        ),
    ],
)
def test_hg19_header_rejects_unverifiable_sequence_dictionaries(header: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        header_with_hg19_m5(header)


def test_committed_ucsc_provenance_independently_binds_every_primary_contig() -> None:
    """A pinned source artifact detects any invented or partially checked M5 table."""
    payload = PROVENANCE.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == "9ee5efd61b116a89ad872ceb6a5358851e839848436b5df533ea2059ca14fdd4"
    assert b"# release=February 2009\n" in payload
    assert b"# source_fasta_url=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz\n" in payload
    assert b"# source_md5_manifest_url=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt\n" in payload
    assert b"# source_fasta_archive_md5=806c02398f5ac5da8ffd6da2d1d5d1a9\n" in payload

    observed: dict[str, tuple[int, str]] = {}
    for line in payload.decode().splitlines():
        if line and not line.startswith("#"):
            name, length, md5 = line.split("\t")
            observed[name] = (int(length), md5)

    assert len(observed) == 25
    assert observed == HG19_PRIMARY_CONTIGS


def test_registered_b178_index_evidence_accepts_only_the_pinned_locus() -> None:
    validate_registered_b178_index_evidence(
        REGISTERED_B178_INDEXED_REGION_DIGEST,
        REGISTERED_B178_INDEXED_REGION_RECORDS,
        REGISTERED_B178_INDEXED_REGION_DIGEST,
        REGISTERED_B178_INDEXED_REGION_RECORDS,
    )


def test_registered_b178_index_evidence_rejects_both_empty_results() -> None:
    empty_digest = hashlib.sha256(b"").hexdigest()

    with pytest.raises(LossyConversionError, match="expected 13868 registered source records"):
        validate_registered_b178_index_evidence(empty_digest, 0, empty_digest, 0)


@pytest.mark.parametrize(
    ("source_digest", "source_records"),
    [
        (REGISTERED_B178_INDEXED_REGION_DIGEST, REGISTERED_B178_INDEXED_REGION_RECORDS - 1),
        ("wrong-locus-digest", REGISTERED_B178_INDEXED_REGION_RECORDS),
    ],
)
def test_registered_b178_index_evidence_rejects_a_wrong_source_locus(source_digest: str, source_records: int) -> None:
    with pytest.raises(LossyConversionError, match="registered source locus"):
        validate_registered_b178_index_evidence(source_digest, source_records, source_digest, source_records)


@pytest.mark.parametrize(
    ("decoded_digest", "decoded_records"),
    [
        (REGISTERED_B178_INDEXED_REGION_DIGEST, REGISTERED_B178_INDEXED_REGION_RECORDS - 1),
        ("decoded-subset", REGISTERED_B178_INDEXED_REGION_RECORDS),
    ],
)
def test_registered_b178_index_evidence_rejects_lossy_crai_results(decoded_digest: str, decoded_records: int) -> None:
    with pytest.raises(LossyConversionError, match="indexed region is not lossless"):
        validate_registered_b178_index_evidence(
            REGISTERED_B178_INDEXED_REGION_DIGEST,
            REGISTERED_B178_INDEXED_REGION_RECORDS,
            decoded_digest,
            decoded_records,
        )
