"""`reference_provenance.py`: the per-file install ledger `canonical_reference_keys()`
now requires a record from, rather than mere file existence.
"""

from __future__ import annotations

import json

import pytest

from vntyper.scripts.reference_provenance import (
    PROVENANCE_FILENAME,
    build_record,
    has_record,
    load_provenance,
    merge,
    provenance_path,
    relative_posix,
)

pytestmark = pytest.mark.unit


class TestProvenancePath:
    def test_lives_at_the_tree_root(self, tmp_path):
        assert provenance_path(tmp_path) == tmp_path / PROVENANCE_FILENAME


class TestLoadProvenance:
    def test_missing_ledger_is_an_empty_dict(self, tmp_path):
        assert load_provenance(tmp_path) == {}

    def test_a_written_ledger_round_trips(self, tmp_path):
        merge(tmp_path, {"alignment/chr1.hg19.fa": build_record(sha256="a" * 64, source="bundle")})
        records = load_provenance(tmp_path)
        assert set(records) == {"alignment/chr1.hg19.fa"}
        assert records["alignment/chr1.hg19.fa"]["sha256"] == "a" * 64

    def test_corrupt_json_degrades_to_empty_rather_than_raising(self, tmp_path):
        provenance_path(tmp_path).write_text("{not json", encoding="utf-8")
        assert load_provenance(tmp_path) == {}

    def test_a_json_array_degrades_to_empty_rather_than_raising(self, tmp_path):
        """The ledger's contract is an object keyed by path; anything else is corrupt."""
        provenance_path(tmp_path).write_text("[1, 2, 3]", encoding="utf-8")
        assert load_provenance(tmp_path) == {}


class TestBuildRecord:
    def test_bundle_shape_carries_asset_and_release_tag(self):
        rec = build_record(sha256="f" * 64, source="bundle", asset="refs-hg19.tar.gz", release_tag="refs-v1")
        assert rec["source"] == "bundle"
        assert rec["sha256"] == "f" * 64
        assert rec["asset"] == "refs-hg19.tar.gz"
        assert rec["release_tag"] == "refs-v1"
        assert rec["source_url"] is None
        assert "installed_at" in rec

    def test_from_source_shape_carries_the_source_url(self):
        rec = build_record(sha256="b" * 64, source="from-source", source_url="https://example.com/chr1.fa.gz")
        assert rec["source"] == "from-source"
        assert rec["asset"] is None
        assert rec["release_tag"] is None
        assert rec["source_url"] == "https://example.com/chr1.fa.gz"


class TestMerge:
    def test_empty_new_records_does_not_touch_an_existing_ledger(self, tmp_path):
        ledger = provenance_path(tmp_path)
        merge(tmp_path, {"a.fa": build_record(sha256="1" * 64, source="bundle")})
        before = ledger.read_text(encoding="utf-8")

        merge(tmp_path, {})

        assert ledger.read_text(encoding="utf-8") == before

    def test_a_second_merge_adds_without_dropping_the_first(self, tmp_path):
        merge(tmp_path, {"hg19.fa": build_record(sha256="1" * 64, source="bundle")})
        merge(tmp_path, {"hg38.fa": build_record(sha256="2" * 64, source="bundle")})

        records = load_provenance(tmp_path)
        assert set(records) == {"hg19.fa", "hg38.fa"}

    def test_the_newest_record_for_a_path_wins(self, tmp_path):
        merge(tmp_path, {"hg19.fa": build_record(sha256="1" * 64, source="bundle")})
        merge(tmp_path, {"hg19.fa": build_record(sha256="2" * 64, source="from-source")})

        assert load_provenance(tmp_path)["hg19.fa"]["sha256"] == "2" * 64

    def test_the_write_is_atomic_no_tmp_file_survives(self, tmp_path):
        merge(tmp_path, {"a.fa": build_record(sha256="1" * 64, source="bundle")})
        assert not provenance_path(tmp_path).with_suffix(".json.tmp").exists()
        assert json.loads(provenance_path(tmp_path).read_text(encoding="utf-8"))


class TestHasRecord:
    def test_present_key(self):
        assert has_record({"a.fa": {}}, "a.fa") is True

    def test_absent_key(self):
        assert has_record({"a.fa": {}}, "b.fa") is False

    def test_empty_ledger(self):
        assert has_record({}, "a.fa") is False


class TestRelativePosix:
    def test_a_nested_file_is_expressed_relative_to_root(self, tmp_path):
        target = tmp_path / "alignment" / "chr1.hg19.fa"
        target.parent.mkdir(parents=True)
        target.write_text("x")
        assert relative_posix(target, tmp_path) == "alignment/chr1.hg19.fa"

    def test_a_path_outside_root_raises(self, tmp_path):
        outside = tmp_path.parent / "elsewhere.fa"
        with pytest.raises(ValueError):
            relative_posix(outside, tmp_path)
