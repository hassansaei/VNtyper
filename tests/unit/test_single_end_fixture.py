"""The derived single-end BAM keeps every read while clearing pairing flags."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Any

import pysam
import pytest

from scripts import single_end_fixture
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


def _artifact_state(path: Path) -> tuple[int, str]:
    return path.stat().st_ino, hashlib.sha256(path.read_bytes()).hexdigest()


def _prior_final_pair(output: Path) -> tuple[Path, Path]:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(b"prior BAM")
    index = output.with_suffix(".bam.bai")
    index.write_bytes(b"prior BAI")
    return output, index


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


def test_a_hardlinked_output_is_rejected_without_changing_the_source(tmp_path: Path) -> None:
    source = _source_bam(tmp_path / "source.bam")
    output = tmp_path / "hardlink.bam"
    os.link(source, output)
    before = _artifact_state(source)

    with pytest.raises(ValueError, match="same file"):
        parse_single_end_fixture(
            {
                "name": "hardlink",
                "kind": "single_end_bam",
                "source_bam": source.name,
                "output_bam": output.name,
            },
            root=tmp_path,
        )

    assert _artifact_state(source) == before


def test_an_output_parent_symlink_cannot_escape_the_declared_root(tmp_path: Path) -> None:
    root = tmp_path / "repo"
    outside = tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    _source_bam(root / "source.bam")
    (root / "escape").symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match="contained"):
        parse_single_end_fixture(
            {
                "name": "escape",
                "kind": "single_end_bam",
                "source_bam": "source.bam",
                "output_bam": "escape/output.bam",
            },
            root=root,
        )


def test_an_existing_output_symlink_is_rejected_without_touching_its_target(tmp_path: Path) -> None:
    _source_bam(tmp_path / "source.bam")
    target = tmp_path / "target.bam"
    target.write_bytes(b"outside fixture")
    output = tmp_path / "output.bam"
    output.symlink_to(target)
    before = _artifact_state(target)

    with pytest.raises(ValueError, match="regular file"):
        parse_single_end_fixture(
            {
                "name": "symlink",
                "kind": "single_end_bam",
                "source_bam": "source.bam",
                "output_bam": "output.bam",
            },
            root=tmp_path,
        )

    assert _artifact_state(target) == before


def test_index_failure_keeps_the_prior_final_pair_and_removes_temps(tmp_path: Path, monkeypatch) -> None:
    source = _source_bam(tmp_path / "source.bam")
    output, index = _prior_final_pair(tmp_path / "derived" / "single.bam")
    before = (_artifact_state(output), _artifact_state(index))
    spec = parse_single_end_fixture(
        {
            "name": "single",
            "kind": "single_end_bam",
            "source_bam": source.name,
            "output_bam": str(output.relative_to(tmp_path)),
        },
        root=tmp_path,
    )
    monkeypatch.setattr(
        single_end_fixture.pysam.samtools, "index", lambda *args: (_ for _ in ()).throw(RuntimeError("index failed"))
    )

    with pytest.raises(RuntimeError, match="index failed"):
        derive_single_end_bam(spec)

    assert (_artifact_state(output), _artifact_state(index)) == before
    assert set(output.parent.iterdir()) == {output, index}


def test_midwrite_failure_keeps_the_prior_final_pair_and_removes_temps(tmp_path: Path, monkeypatch) -> None:
    source = _source_bam(tmp_path / "source.bam")
    output, index = _prior_final_pair(tmp_path / "derived" / "single.bam")
    before = (_artifact_state(output), _artifact_state(index))
    spec = parse_single_end_fixture(
        {
            "name": "single",
            "kind": "single_end_bam",
            "source_bam": source.name,
            "output_bam": str(output.relative_to(tmp_path)),
        },
        root=tmp_path,
    )
    real_alignment_file = pysam.AlignmentFile

    class FailingWriter:
        def __init__(self, handle: pysam.AlignmentFile) -> None:
            self.handle = handle
            self.writes = 0

        def __enter__(self) -> FailingWriter:
            return self

        def __exit__(self, *args: object) -> None:
            self.handle.close()

        def write(self, read: pysam.AlignedSegment) -> None:
            self.writes += 1
            if self.writes == 2:
                raise RuntimeError("write failed")
            self.handle.write(read)

    def alignment_file(path: str, mode: Any, **kwargs: Any) -> Any:
        handle = real_alignment_file(path, mode, **kwargs)
        return FailingWriter(handle) if mode == "wb" else handle

    monkeypatch.setattr(single_end_fixture.pysam, "AlignmentFile", alignment_file)

    with pytest.raises(RuntimeError, match="write failed"):
        derive_single_end_bam(spec)

    assert (_artifact_state(output), _artifact_state(index)) == before
    assert set(output.parent.iterdir()) == {output, index}


def test_second_artifact_install_failure_rolls_back_the_prior_pair(tmp_path: Path, monkeypatch) -> None:
    source = _source_bam(tmp_path / "source.bam")
    output, index = _prior_final_pair(tmp_path / "derived" / "single.bam")
    before = (_artifact_state(output), _artifact_state(index))
    spec = parse_single_end_fixture(
        {
            "name": "single",
            "kind": "single_end_bam",
            "source_bam": source.name,
            "output_bam": str(output.relative_to(tmp_path)),
        },
        root=tmp_path,
    )
    real_replace = os.replace
    failed = False

    def fail_first_index_install(source_path: str | Path, destination_path: str | Path) -> None:
        nonlocal failed
        if Path(destination_path) == index and not failed:
            failed = True
            raise OSError("index install failed")
        real_replace(source_path, destination_path)

    monkeypatch.setattr(os, "replace", fail_first_index_install)

    with pytest.raises(OSError, match="index install failed"):
        derive_single_end_bam(spec)

    assert (_artifact_state(output), _artifact_state(index)) == before
    assert set(output.parent.iterdir()) == {output, index}


def test_backup_creation_failure_keeps_the_prior_pair_and_original_error(tmp_path: Path, monkeypatch) -> None:
    source = _source_bam(tmp_path / "source.bam")
    output, index = _prior_final_pair(tmp_path / "derived" / "single.bam")
    before = (_artifact_state(output), _artifact_state(index))
    spec = parse_single_end_fixture(
        {
            "name": "single",
            "kind": "single_end_bam",
            "source_bam": source.name,
            "output_bam": str(output.relative_to(tmp_path)),
        },
        root=tmp_path,
    )
    real_link = os.link

    def fail_index_backup(source_path: str | Path, destination_path: str | Path) -> None:
        if Path(source_path) == index:
            raise OSError("index backup failed")
        real_link(source_path, destination_path)

    monkeypatch.setattr(os, "link", fail_index_backup)

    with pytest.raises(OSError, match="index backup failed"):
        derive_single_end_bam(spec)

    assert (_artifact_state(output), _artifact_state(index)) == before
    assert set(output.parent.iterdir()) == {output, index}


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
