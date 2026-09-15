from __future__ import annotations

import hashlib

import pytest

from vntyper.scripts.length_standard_io import read_standard_length_features

pytestmark = pytest.mark.unit


class _Reference:
    references = ("chr1",)

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def fetch(self, contig, start, end):
        assert (contig, start, end) == ("chr1", 155188296, 155192429)
        return "A" * (end - start)


class _Header:
    def to_dict(self):
        return {"SQ": [{"SN": "chr1", "LN": 248956422, "AS": "GRCh38"}]}


class _Record:
    flag = 0
    query_name = "fragment"
    mapping_quality = 0
    query_sequence = "G" * 4133
    query_qualities = [30] * 4133
    cigartuples = ((4, 1), (0, 4132))

    def has_tag(self, name):
        return False

    def get_aligned_pairs(self, matches_only=False):
        assert matches_only is False
        return tuple((offset, 155188296 + offset) for offset in range(4133))


class _Alignment:
    references = ("chr1",)
    lengths = (248956422,)
    header = _Header()

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def fetch(self, contig, start, end):
        assert (contig, start, end) == ("chr1", 155188296, 155192429)
        return (_Record(),)


@pytest.mark.parametrize("assembly", ["GRCh38", "hg38", "hg38_ensembl", "hg38_ncbi"])
def test_standard_reader_uses_one_indexed_locus_scan_and_explicit_index(tmp_path, monkeypatch, assembly: str) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM\x01")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="ascii")
    index = tmp_path / "sample.bam.bai"
    index.write_bytes(b"index")
    expected = hashlib.sha256(("A" * 4133).encode("ascii")).hexdigest()
    monkeypatch.setattr("vntyper.scripts.length_standard_io.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr("vntyper.scripts.length_standard_features.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.FastaFile", _Reference)
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.AlignmentFile", _Alignment)

    measured = read_standard_length_features(
        alignment.resolve(), reference.resolve(), assembly=assembly, index_path=index.resolve()
    )

    assert measured.values["A"] == 1
    assert measured.values["mapq_zero_fraction"] == 1
    assert measured.values["soft_clipped_read_fraction"] == 1
    assert measured.values["query_sequence_gc_fraction"] == 1


def test_standard_reader_rejects_wrong_assembly_or_missing_reference(tmp_path) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM\x01")
    reference = tmp_path / "reference.fa"
    with pytest.raises(ValueError, match="reference"):
        read_standard_length_features(alignment.resolve(), reference.resolve(), assembly="GRCh38")
    reference.write_text(">chr1\nA\n", encoding="ascii")
    with pytest.raises(ValueError, match="GRCh38"):
        read_standard_length_features(alignment.resolve(), reference.resolve(), assembly="GRCh37")


def test_standard_reader_restores_cram_reference_environment_on_failure(tmp_path, monkeypatch) -> None:
    alignment = tmp_path / "sample.cram"
    alignment.write_bytes(b"CRAMbad")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="ascii")
    restored: list[str | None] = []
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pin_reference_resolution", lambda value: "old")
    monkeypatch.setattr("vntyper.scripts.length_standard_io.restore_reference_resolution", restored.append)
    monkeypatch.setattr(
        "vntyper.scripts.length_standard_io.pysam.FastaFile",
        lambda value: (_ for _ in ()).throw(OSError("bad")),
    )
    with pytest.raises(ValueError):
        read_standard_length_features(alignment.resolve(), reference.resolve(), assembly="GRCh38")
    assert restored == ["old"]
