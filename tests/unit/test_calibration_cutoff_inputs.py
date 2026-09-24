"""Capture association rejects duplicated evidence and preserves declared identity."""

from types import MappingProxyType

import pytest

from vntyper.scripts.calibration_cohort_manifest import CohortSample

pytestmark = pytest.mark.unit


def samples(tmp_path):
    return tuple(CohortSample(key, tmp_path / f"{key}.bam", "hg38", True, None, key) for key in ("a", "b"))


def test_capture_paths_resolve_and_native_overrides_are_optional(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    manifest = tmp_path / "captures.tsv"
    manifest.write_text("sample_id\tkestrel_capture\nk\tbad.json\n")
    with pytest.raises(ValueError, match="roster"):
        read_cutoff_captures(manifest, samples(tmp_path), caller="kestrel")
    manifest.write_text("sample_id\tkestrel_capture\na\ta.json\nb\tb.json\n")
    for key in ("a", "b"):
        (tmp_path / f"{key}.json").write_text(key)
    result = read_cutoff_captures(manifest, samples(tmp_path), caller="kestrel")
    assert result.kestrel == {key: tmp_path / f"{key}.json" for key in ("a", "b")}
    assert result.native_kestrel == {"a": None, "b": None}
    assert result.advntr == {}


@pytest.mark.parametrize("row", ["b\ta.json", "a\tb.json"])
def test_duplicate_paths_and_samples_are_rejected(tmp_path, row):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    for name in ("a", "b"):
        (tmp_path / f"{name}.json").write_text("same capture")
    path = tmp_path / "captures.tsv"
    path.write_text(f"sample_id\tkestrel_capture\na\ta.json\n{row}\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_cutoff_captures(path, samples(tmp_path), caller="kestrel")


def test_distinct_samples_may_have_identical_capture_bytes(tmp_path):
    # `calibration-kestrel-capture-v1` carries no internal sample identity, so any two
    # samples whose candidates are all filtered serialise to byte-identical captures.
    # Identical bytes are therefore legitimate; only a repeated path or a repeated
    # sample_id is a duplicated evidence claim. `test_calibration_cutoff_advntr.py::
    # test_distinct_samples_may_have_identical_complete_capture_bytes` asserts the same
    # for the adVNTR side.
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    for name in ("a", "b"):
        (tmp_path / f"{name}.json").write_text("same capture")
    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tkestrel_capture\na\ta.json\nb\tb.json\n")
    result = read_cutoff_captures(path, samples(tmp_path), caller="kestrel")
    assert result.kestrel == {"a": tmp_path / "a.json", "b": tmp_path / "b.json"}


@pytest.mark.parametrize(
    "text",
    [
        "sample_id\tkestrel_capture\textra\na\ta.json\tx\n",
        "sample_id\tsample_id\na\ta\n",
        "sample_id\tkestrel_capture\na\n",
        "sample_id\tkestrel_capture\n a\ta.json\n",
        "sample_id\tkestrel_capture\na\t\n",
    ],
)
def test_malformed_capture_manifest_fails_before_replay(tmp_path, text):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    path = tmp_path / "captures.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_cutoff_captures(path, samples(tmp_path), caller="kestrel")


def test_ad_only_capture_contract_and_primary_groups(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import primary_samples, read_cutoff_captures

    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tadvntr_capture\na\ta.jsonl\nb\tb.jsonl\n")
    for name in ("a", "b"):
        (tmp_path / f"{name}.jsonl").write_text(name)
    rows = samples(tmp_path)
    assert len(read_cutoff_captures(path, rows, caller="advntr").advntr) == 2
    repeated = CohortSample("c", tmp_path / "c.bam", "hg38", True, 100, "a")
    assert primary_samples((*rows, repeated)) == rows
    conflict = CohortSample("c", tmp_path / "c.bam", "hg38", False, None, "a")
    with pytest.raises(ValueError, match="conflicting"):
        primary_samples((*rows, conflict))
    with pytest.raises(ValueError, match="caller"):
        read_cutoff_captures(path, rows, caller="invalid")


def test_ad_only_capture_carries_no_native_kestrel_roster(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tadvntr_capture\na\ta.jsonl\n")
    (tmp_path / "a.jsonl").write_text("a")
    result = read_cutoff_captures(path, samples(tmp_path), caller="advntr")
    assert (result.kestrel, result.native_kestrel) == ({}, {})
    assert result.advntr == {"a": tmp_path / "a.jsonl"}


def test_both_callers_pair_captures_and_sort_by_declared_identity(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    for name in ("a", "b"):
        (tmp_path / f"{name}.json").write_text(name)
        (tmp_path / f"{name}.jsonl").write_text(name)
    (tmp_path / "a.tsv").write_text("native")
    path = tmp_path / "captures.tsv"
    path.write_text(
        "sample_id\tkestrel_capture\tadvntr_capture\tnative_kestrel\nb\tb.json\tb.jsonl\t\na\ta.json\ta.jsonl\ta.tsv\n"
    )
    result = read_cutoff_captures(path, samples(tmp_path), caller="both")
    assert tuple(result.kestrel) == ("a", "b")
    assert result.kestrel == {key: tmp_path / f"{key}.json" for key in ("a", "b")}
    assert result.advntr == {key: tmp_path / f"{key}.jsonl" for key in ("a", "b")}
    assert result.native_kestrel == {"a": tmp_path / "a.tsv", "b": None}
    assert all(
        isinstance(mapping, MappingProxyType) for mapping in (result.kestrel, result.native_kestrel, result.advntr)
    )


@pytest.mark.parametrize(
    ("text", "caller"),
    [
        ("", "kestrel"),
        ("sample_id\na\n", "kestrel"),
        ("sample_id\tkestrel_capture\n", "kestrel"),
        ("sample_id\tadvntr_capture\tunknown_capture\na\ta.jsonl\tx\n", "advntr"),
        ("sample_id\tkestrel_capture\tkestrel_capture\na\ta.json\tb.json\n", "kestrel"),
        ("sample_id\tkestrel_capture\na\ta.json\tx\n", "kestrel"),
    ],
)
def test_column_contract_and_row_geometry_fail_before_replay(tmp_path, text, caller):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    path = tmp_path / "captures.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_cutoff_captures(path, samples(tmp_path), caller=caller)


def test_one_row_cannot_claim_the_same_path_as_capture_and_native(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    (tmp_path / "a.json").write_text("a")
    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tkestrel_capture\tnative_kestrel\na\ta.json\ta.json\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_cutoff_captures(path, samples(tmp_path), caller="kestrel")


def test_manifest_must_be_a_path(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tkestrel_capture\na\ta.json\n")
    with pytest.raises(ValueError, match="Path"):
        read_cutoff_captures(str(path), samples(tmp_path), caller="kestrel")


@pytest.mark.parametrize("missing", ["a.json", "a.tsv"])
def test_declared_captures_must_be_present_regular_files(tmp_path, missing):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    for name in ("a.json", "a.tsv"):
        if name != missing:
            (tmp_path / name).write_text(name)
    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tkestrel_capture\tnative_kestrel\na\ta.json\ta.tsv\n")
    with pytest.raises(ValueError, match="unreadable"):
        read_cutoff_captures(path, samples(tmp_path), caller="kestrel")


@pytest.mark.parametrize("value", ["missing.jsonl", ""])
def test_kestrel_selection_ignores_a_declared_advntr_column(tmp_path, value):
    # One manifest serves --caller kestrel, advntr and both: the unselected caller's
    # column is never resolved, opened or checked, so it may name an absent file or
    # be empty, and it contributes nothing to the returned mappings.
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    (tmp_path / "a.json").write_text("a")
    path = tmp_path / "captures.tsv"
    path.write_text(f"sample_id\tkestrel_capture\tadvntr_capture\na\ta.json\t{value}\n")
    result = read_cutoff_captures(path, samples(tmp_path), caller="kestrel")
    assert result.kestrel == {"a": tmp_path / "a.json"}
    assert result.native_kestrel == {"a": None}
    assert result.advntr == {}


def test_an_ignored_column_is_not_counted_as_a_duplicate_path(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    for name in ("a", "b"):
        (tmp_path / f"{name}.json").write_text(name)
    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tkestrel_capture\tadvntr_capture\na\ta.json\tb.json\nb\tb.json\tb.json\n")
    result = read_cutoff_captures(path, samples(tmp_path), caller="kestrel")
    assert result.kestrel == {key: tmp_path / f"{key}.json" for key in ("a", "b")}


def test_advntr_selection_ignores_the_kestrel_and_native_kestrel_columns(tmp_path):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    (tmp_path / "a.jsonl").write_text("a")
    path = tmp_path / "captures.tsv"
    path.write_text(
        "sample_id\tkestrel_capture\tnative_kestrel\tadvntr_capture\na\tmissing.json\tmissing.tsv\ta.jsonl\n"
    )
    result = read_cutoff_captures(path, samples(tmp_path), caller="advntr")
    assert (result.kestrel, result.native_kestrel) == ({}, {})
    assert result.advntr == {"a": tmp_path / "a.jsonl"}


@pytest.mark.parametrize("caller", ["kestrel", "advntr"])
def test_an_unknown_column_is_still_refused_beside_the_other_callers_columns(tmp_path, caller):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    (tmp_path / "a.json").write_text("a")
    (tmp_path / "a.jsonl").write_text("a")
    path = tmp_path / "captures.tsv"
    path.write_text("sample_id\tkestrel_capture\tadvntr_capture\textra\na\ta.json\ta.jsonl\tx\n")
    with pytest.raises(ValueError, match="columns"):
        read_cutoff_captures(path, samples(tmp_path), caller=caller)


@pytest.mark.parametrize("header", ["sample_id\tkestrel_capture", "sample_id\tadvntr_capture"])
def test_both_still_requires_both_capture_columns(tmp_path, header):
    from vntyper.scripts.calibration_cutoff_inputs import read_cutoff_captures

    (tmp_path / "a.json").write_text("a")
    path = tmp_path / "captures.tsv"
    path.write_text(f"{header}\na\ta.json\n")
    with pytest.raises(ValueError, match="columns"):
        read_cutoff_captures(path, samples(tmp_path), caller="both")
