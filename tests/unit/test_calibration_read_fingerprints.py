"""Bounded external sorting preserves read multiplicity and cleans private files."""

from importlib import import_module

import pytest

from tests.unit.test_calibration_read_identity import record

pytestmark = pytest.mark.unit


def test_streaming_fingerprints_ignore_order_and_chunk_boundaries(tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")
    records = [record(name=f"synthetic-{i}") for i in range(13)]
    small = m.fingerprint_read_records(iter(records), temporary_parent=tmp_path, chunk_records=2, merge_fan_in=2)
    large = m.fingerprint_read_records(iter(reversed(records)), temporary_parent=tmp_path, chunk_records=100)
    assert small == large
    assert small.primary_record_count == 13
    assert small.sequence_identity_reliable
    assert small.reasons == ()
    assert list(tmp_path.iterdir()) == []


def test_identical_read_occurrences_are_counted_not_collapsed(tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")
    first = m.fingerprint_read_records([record()], temporary_parent=tmp_path)
    repeated = m.fingerprint_read_records([record(), record()], temporary_parent=tmp_path, chunk_records=1)
    assert first.primary_record_count == 1
    assert repeated.primary_record_count == 2
    for field in ("alignment_sha256", "named_sequence_sha256", "unnamed_sequence_sha256"):
        assert getattr(first, field) != getattr(repeated, field)


def test_renamed_reads_are_a_sequence_only_collision_not_an_identical_named_input(tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")
    first = m.fingerprint_read_records([record()], temporary_parent=tmp_path)
    renamed = m.fingerprint_read_records([record(name="renamed")], temporary_parent=tmp_path)
    assert first.alignment_sha256 != renamed.alignment_sha256
    assert first.named_sequence_sha256 != renamed.named_sequence_sha256
    assert first.unnamed_sequence_sha256 == renamed.unnamed_sequence_sha256


def test_nonprimary_records_are_excluded_and_reconstruction_defects_remain_visible(tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")
    result = m.fingerprint_read_records(
        [
            record(flags=321),
            record(flags=2113),
            record(cigar=((5, 2), (0, 4))),
            record(sequence=None, qualities=None),
        ],
        temporary_parent=tmp_path,
        chunk_records=1,
    )
    assert result.primary_record_count == 2
    assert not result.sequence_identity_reliable
    assert result.reasons == ("hard_clipped_sequence", "missing_qualities", "missing_sequence")


def test_stream_failure_removes_all_temporary_chunks_and_never_returns_partial_fingerprint(tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")

    def failing_records():
        yield record()
        yield record(name="second")
        assert len(list(tmp_path.iterdir())) == 1
        directory = next(tmp_path.iterdir())
        assert directory.stat().st_mode & 0o777 == 0o700
        assert list(directory.iterdir())
        raise RuntimeError("synthetic truncated input")

    with pytest.raises(RuntimeError, match="truncated"):
        m.fingerprint_read_records(failing_records(), temporary_parent=tmp_path, chunk_records=1)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("records", [[], [record(flags=321)]])
def test_empty_primary_read_evidence_fails_closed(records, tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")
    with pytest.raises(ValueError, match="primary"):
        m.fingerprint_read_records(records, temporary_parent=tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "kwargs",
    [
        {"chunk_records": True},
        {"chunk_records": 0},
        {"chunk_records": 1000001},
        {"merge_fan_in": 1},
        {"merge_fan_in": True},
        {"merge_fan_in": 65},
    ],
)
def test_resource_bounds_are_strict_and_validated_before_iteration(kwargs, tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")

    def forbidden():
        raise AssertionError("input opened before resource validation")
        yield

    with pytest.raises(ValueError):
        m.fingerprint_read_records(forbidden(), temporary_parent=tmp_path, **kwargs)
    assert list(tmp_path.iterdir()) == []


def test_truncated_sort_chunk_is_detected_before_any_fingerprint_is_returned(tmp_path, monkeypatch):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")
    original = m._write_chunk

    def truncate_chunk(values, directory):
        path = original(values, directory)
        path.write_text(min(values) + "\n", encoding="ascii")
        return path

    monkeypatch.setattr(m, "_write_chunk", truncate_chunk)
    with pytest.raises(ValueError, match="record count"):
        m.fingerprint_read_records([record(), record(name="second")], temporary_parent=tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_many_tiny_chunks_compact_runs_online_with_logarithmic_retained_metadata(tmp_path):
    m = import_module("vntyper.scripts.calibration_read_fingerprints")
    records = [record(name=f"synthetic-{index % 17}") for index in range(128)]
    peak_retained_run_files = 0

    def observed_records():
        nonlocal peak_retained_run_files
        for item in records:
            directories = list(tmp_path.iterdir())
            if directories:
                peak_retained_run_files = max(peak_retained_run_files, len(list(directories[0].iterdir())))
            yield item

    compacted = m.fingerprint_read_records(
        observed_records(), temporary_parent=tmp_path, chunk_records=1, merge_fan_in=2
    )
    baseline = m.fingerprint_read_records(records, temporary_parent=tmp_path, chunk_records=128, merge_fan_in=2)

    assert compacted == baseline
    assert compacted.primary_record_count == 128
    assert peak_retained_run_files <= 24
    assert list(tmp_path.iterdir()) == []
