"""The derived single-end BAM keeps every read while clearing pairing flags."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from scripts.single_end_fixture import (
    derive_single_end_bam,
    parse_single_end_fixture,
)
from tests.parametrization import load_test_config

pytestmark = pytest.mark.unit


def _source_bam(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for index, flag in enumerate((75, 155)):
            read = pysam.AlignedSegment()
            read.query_name = f"read-{index}"
            read.query_sequence = "ACGT"
            read.flag = flag
            read.reference_id = 0
            read.reference_start = 100 + index
            read.mapping_quality = 60
            read.cigarstring = "4M"
            read.next_reference_id = 0
            read.next_reference_start = 200 + index
            read.template_length = 104
            read.query_qualities = pysam.qualitystring_to_array("IIII")
            handle.write(read)
    return path


def test_a_portable_single_end_fixture_declaration_is_parsed() -> None:
    spec = parse_single_end_fixture(
        {
            "name": "positive_single_end",
            "kind": "single_end_bam",
            "source_bam": "tests/data/source.bam",
            "output_bam": "tests/data/derived/source_single_end.bam",
        }
    )

    assert spec.name == "positive_single_end"
    assert spec.source_bam == Path("tests/data/source.bam")
    assert spec.output_bam == Path("tests/data/derived/source_single_end.bam")


def test_the_single_end_fixture_and_integration_case_are_registered() -> None:
    config = load_test_config()
    declarations = {entry["name"]: entry for entry in config["derived_fixtures"]}
    entry = declarations["example_b178_hg19_single_end"]

    spec = parse_single_end_fixture(entry)

    assert spec.source_bam == Path("tests/data/example_b178_hg19_subset.bam")
    assert spec.output_bam == Path("tests/data/derived/single_end/example_b178_hg19_single_end.bam")
    cases = config["integration_tests"]["single_end_bam_tests"]
    assert [case["fixture_name"] for case in cases] == [spec.name]


@pytest.mark.parametrize(
    "entry",
    [
        {},
        {"name": "case", "kind": "cram", "source_bam": "source.bam", "output_bam": "output.bam"},
        {"name": "", "kind": "single_end_bam", "source_bam": "source.bam", "output_bam": "output.bam"},
        {"name": "case", "kind": "single_end_bam", "source_bam": "", "output_bam": "output.bam"},
        {"name": "case", "kind": "single_end_bam", "source_bam": "source.bam", "output_bam": ""},
        {
            "name": "case",
            "kind": "single_end_bam",
            "source_bam": "/absolute/source.bam",
            "output_bam": "output.bam",
        },
        {
            "name": "case",
            "kind": "single_end_bam",
            "source_bam": "source.bam",
            "output_bam": "/absolute/output.bam",
        },
        {
            "name": "case",
            "kind": "single_end_bam",
            "source_bam": "same.bam",
            "output_bam": "same.bam",
        },
        {
            "name": "nested/case",
            "kind": "single_end_bam",
            "source_bam": "source.bam",
            "output_bam": "output.bam",
        },
        {
            "name": "case",
            "kind": "single_end_bam",
            "source_bam": "../source.bam",
            "output_bam": "output.bam",
        },
        {
            "name": "case",
            "kind": "single_end_bam",
            "source_bam": "source.bam",
            "output_bam": "output.cram",
        },
    ],
)
def test_malformed_or_nonportable_declarations_fail_closed(entry: dict[str, object]) -> None:
    with pytest.raises(ValueError, match="single-end fixture"):
        parse_single_end_fixture(entry)


def test_derivation_clears_pairing_flags_without_dropping_or_rewriting_reads(tmp_path: Path) -> None:
    source = _source_bam(tmp_path / "source.bam")
    output = tmp_path / "derived" / "single.bam"
    spec = parse_single_end_fixture(
        {
            "name": "single",
            "kind": "single_end_bam",
            "source_bam": str(source.relative_to(tmp_path)),
            "output_bam": str(output.relative_to(tmp_path)),
        },
        root=tmp_path,
    )

    records = derive_single_end_bam(spec)

    assert records == 2
    assert output.with_suffix(".bam.bai").is_file()
    with (
        pysam.AlignmentFile(str(source), "rb") as source_handle,
        pysam.AlignmentFile(str(output), "rb") as output_handle,
    ):
        before = list(source_handle.fetch(until_eof=True))
        after = list(output_handle.fetch(until_eof=True))
    assert [read.query_name for read in after] == [read.query_name for read in before]
    assert [read.query_sequence for read in after] == [read.query_sequence for read in before]
    assert [read.reference_start for read in after] == [read.reference_start for read in before]
    assert [read.flag for read in before] == [75, 155]
    assert [read.flag for read in after] == [0, 16]
    for read in after:
        assert read.is_paired is False
        assert read.is_proper_pair is False
        assert read.mate_is_unmapped is False
        assert read.is_read1 is False
        assert read.is_read2 is False


def test_derivation_refuses_a_missing_source_before_creating_output(tmp_path: Path) -> None:
    spec = parse_single_end_fixture(
        {
            "name": "missing",
            "kind": "single_end_bam",
            "source_bam": "missing.bam",
            "output_bam": "derived/missing.bam",
        },
        root=tmp_path,
    )

    with pytest.raises(FileNotFoundError, match="missing.bam"):
        derive_single_end_bam(spec)

    assert not spec.output_bam.exists()
