from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from vntyper.scripts.length_standard_io import read_standard_length_features

pytestmark = pytest.mark.unit


class _Reference:
    references = ("chr1",)
    init_calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def __init__(self, *args, **kwargs):
        self.init_calls.append((args, kwargs))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def fetch(self, contig, start, end):
        assert (contig, start, end) == ("chr1", 155188296, 155192429)
        return "A" * (end - start)


class _Header:
    def __init__(self, assembly: str = "GRCh38") -> None:
        self.assembly = assembly

    def to_dict(self):
        return {"SQ": [{"SN": "chr1", "LN": 248956422, "AS": self.assembly}]}


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


class _CaseRecord:
    def __init__(
        self,
        *,
        flag: int,
        query_name: str,
        mapping_quality: int,
        query_sequence: str | None,
        cigartuples: tuple[tuple[int, int], ...],
        pairs: tuple[tuple[int | None, int | None], ...],
    ) -> None:
        self.flag = flag
        self.query_name = query_name
        self.mapping_quality = mapping_quality
        self.query_sequence = query_sequence
        self.query_qualities = None if query_sequence is None else [30] * len(query_sequence)
        self.cigartuples = cigartuples
        self._pairs = pairs

    def has_tag(self, name: str) -> bool:
        return name == "RG"

    def get_tag(self, name: str) -> str:
        assert name == "RG"
        return "shared-read-group"

    def get_aligned_pairs(self, matches_only: bool = False):
        assert matches_only is False
        return self._pairs


class _CaseAlignment(_Alignment):
    def fetch(self, contig, start, end):
        assert (contig, start, end) == ("chr1", 155188296, 155192429)
        invariant_position = 155188500
        core_position = 155188800
        included = (
            _CaseRecord(
                flag=0x41,
                query_name="overlapping-pair",
                mapping_quality=10,
                query_sequence="AC",
                cigartuples=((0, 2),),
                pairs=((0, invariant_position), (1, core_position), (None, core_position + 1), (0, None)),
            ),
            _CaseRecord(
                flag=0x91,
                query_name="overlapping-pair",
                mapping_quality=20,
                query_sequence="GT",
                cigartuples=((0, 2),),
                pairs=((0, invariant_position), (1, core_position)),
            ),
            _CaseRecord(
                flag=0x800,
                query_name="supplementary",
                mapping_quality=0,
                query_sequence="GG",
                cigartuples=((4, 1), (0, 1)),
                pairs=((0, None), (1, invariant_position)),
            ),
            _CaseRecord(
                flag=0,
                query_name="missing-sequence",
                mapping_quality=60,
                query_sequence=None,
                cigartuples=((0, 1),),
                pairs=((0, invariant_position),),
            ),
        )
        excluded = tuple(
            _CaseRecord(
                flag=flag,
                query_name=f"excluded-{flag}",
                mapping_quality=60,
                query_sequence="CC",
                cigartuples=((0, 1),),
                pairs=((0, invariant_position),),
            )
            for flag in (0x4, 0x100, 0x200, 0x400)
        )
        return (*included, *excluded)


@pytest.mark.parametrize("assembly", ["GRCh38", "hg38", "hg38_ensembl", "hg38_ncbi"])
@pytest.mark.parametrize("compressed", [False, True])
def test_standard_reader_uses_one_indexed_locus_scan_and_explicit_index(
    tmp_path, monkeypatch, assembly: str, compressed: bool
) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM\x01")
    reference = tmp_path / ("reference.fa.gz" if compressed else "reference.fa")
    reference.write_bytes(b"\x1f\x8bcompressed" if compressed else b">chr1\nA\n")
    Path(f"{reference}.fai").write_text("chr1\t248956422\t6\t1\t2\n", encoding="ascii")
    if compressed:
        Path(f"{reference}.gzi").write_bytes(b"compressed-index")
    index = tmp_path / "sample.bam.bai"
    index.write_bytes(b"index")
    expected = hashlib.sha256(("A" * 4133).encode("ascii")).hexdigest()
    monkeypatch.setattr("vntyper.scripts.length_standard_io.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr("vntyper.scripts.length_standard_features.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr(_Reference, "init_calls", [])
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.FastaFile", _Reference)
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.AlignmentFile", _Alignment)
    monkeypatch.setattr(_Alignment, "header", _Header(assembly))

    measured = read_standard_length_features(
        alignment.resolve(), reference.resolve(), assembly=assembly, index_path=index.resolve()
    )

    assert measured.values["A"] == 1
    assert measured.values["mapq_zero_fraction"] == 1
    assert measured.values["soft_clipped_read_fraction"] == 1
    assert measured.values["query_sequence_gc_fraction"] == 1
    assert _Reference.init_calls == [
        (
            (str(reference.resolve()),),
            {
                "filepath_index": f"{reference.resolve()}.fai",
                "filepath_index_compressed": f"{reference.resolve()}.gzi" if compressed else None,
            },
        )
    ]


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
    Path(f"{reference}.fai").write_text("chr1\t248956422\t6\t1\t2\n", encoding="ascii")
    restored: list[str | None] = []
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pin_reference_resolution", lambda value: "old")
    monkeypatch.setattr("vntyper.scripts.length_standard_io.restore_reference_resolution", restored.append)
    monkeypatch.setattr(
        "vntyper.scripts.length_standard_io.pysam.FastaFile",
        lambda value, **kwargs: (_ for _ in ()).throw(OSError("bad")),
    )
    with pytest.raises(ValueError):
        read_standard_length_features(alignment.resolve(), reference.resolve(), assembly="GRCh38")
    assert restored == ["old"]


def test_standard_reader_applies_frozen_record_depth_and_fragment_policy(tmp_path, monkeypatch) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM\x01")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="ascii")
    Path(f"{reference}.fai").write_text("chr1\t248956422\t6\t1\t2\n", encoding="ascii")
    index = tmp_path / "sample.bam.bai"
    index.write_bytes(b"index")
    expected = hashlib.sha256(("A" * 4133).encode("ascii")).hexdigest()
    monkeypatch.setattr("vntyper.scripts.length_standard_io.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr("vntyper.scripts.length_standard_features.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.FastaFile", _Reference)
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.AlignmentFile", _CaseAlignment)

    measured = read_standard_length_features(
        alignment.resolve(), reference.resolve(), assembly="GRCh38", index_path=index.resolve()
    )

    # Forward and reverse mates overlap at both counted bases: depth counts both,
    # while their shared query name/read group contributes one fragment. The
    # supplementary and missing-sequence records bring invariant depth to four
    # with three distinct fragments. Deletion/non-reference pairs and 0x704
    # records contribute neither depth nor read statistics.
    assert measured.values["A"] == pytest.approx((2 / 3213) / (4 / 540))
    assert measured.qc.invariant_mean_depth == pytest.approx(4 / 540)
    assert measured.qc.invariant_supporting_fragments == 3
    assert measured.qc.eligible_read_count == 3
    assert measured.values["mapq_zero_fraction"] == pytest.approx(1 / 3)
    assert measured.values["mean_mapq"] == 10
    assert measured.values["soft_clipped_read_fraction"] == pytest.approx(1 / 3)
    assert measured.values["query_sequence_gc_fraction"] == pytest.approx(2 / 3)


def test_standard_reader_requires_existing_fasta_index_before_opening_reference(tmp_path, monkeypatch) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM\x01")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="ascii")
    implicit_index = Path(f"{reference}.fai")
    opened = False

    def index_creating_fasta(path):
        nonlocal opened
        opened = True
        implicit_index.write_text("implicitly created", encoding="ascii")
        return _Reference(path)

    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.FastaFile", index_creating_fasta)
    with pytest.raises(ValueError, match="reference index"):
        read_standard_length_features(alignment.resolve(), reference.resolve(), assembly="GRCh38")
    assert opened is False
    assert implicit_index.exists() is False


def test_standard_reader_requires_existing_compressed_fasta_index_before_opening_reference(
    tmp_path, monkeypatch
) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM\x01")
    reference = tmp_path / "reference.fa.gz"
    reference.write_bytes(b"\x1f\x8bcompressed")
    Path(f"{reference}.fai").write_text("chr1\t248956422\t6\t1\t2\n", encoding="ascii")
    implicit_index = Path(f"{reference}.gzi")
    opened = False

    def index_creating_fasta(path):
        nonlocal opened
        opened = True
        implicit_index.write_text("implicitly created", encoding="ascii")
        return _Reference(path)

    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.FastaFile", index_creating_fasta)
    with pytest.raises(ValueError, match="compressed reference index"):
        read_standard_length_features(alignment.resolve(), reference.resolve(), assembly="GRCh38")
    assert opened is False
    assert implicit_index.exists() is False


def test_standard_reader_rejects_wrong_header_assembly(tmp_path, monkeypatch) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM\x01")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="ascii")
    Path(f"{reference}.fai").write_text("chr1\t248956422\t6\t1\t2\n", encoding="ascii")
    expected = hashlib.sha256(("A" * 4133).encode("ascii")).hexdigest()
    monkeypatch.setattr("vntyper.scripts.length_standard_io.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr("vntyper.scripts.length_standard_features.STANDARD_REFERENCE_LOCUS_SHA256", expected)
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.FastaFile", _Reference)
    monkeypatch.setattr("vntyper.scripts.length_standard_io.pysam.AlignmentFile", _Alignment)
    monkeypatch.setattr(_Alignment, "header", _Header("GRCh37"))

    with pytest.raises(RuntimeError, match="indexed alignment locus"):
        read_standard_length_features(alignment.resolve(), reference.resolve(), assembly="GRCh38")
