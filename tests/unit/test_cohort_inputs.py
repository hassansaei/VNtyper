"""Finding the samples of a cohort on disk, and reading what the pipeline recorded.

A cohort is whatever `vntyper cohort` was pointed at: directories that are themselves a
sample, directories containing many, zip files of either, and anything that turns out to
be none of those. This module locates the `pipeline_summary.json` files in that, and
reads each one into the three pieces the report needs - the Kestrel rows, the adVNTR
rows, and the per-sample statistics.

This was 120 lines of `cohort_summary.py` reachable only by calling the whole cohort
pipeline, and it was the largest single uncovered block in the file. Everything here is
**characterisation**: it records what a cohort run does today, including the two
behaviours named `..._today`.

Step names are matched by exact string comparison against what `pipeline.py` writes
(AGENTS.md trap 5), so this file asserts against `summary_steps.STEP_*` for the same
reason the module does.
"""

from __future__ import annotations

import json
import logging
import zipfile
from pathlib import Path

import pytest

from vntyper.scripts.cohort_inputs import (
    cleanup_temp_dirs,
    discover_sample_directories,
    load_pipeline_summary_for_sample,
    parse_pipeline_summary,
    pseudonymized_sample_name,
)
from vntyper.scripts.summary_steps import (
    STEP_ADVNTR,
    STEP_BAM_HEADER,
    STEP_COVERAGE,
    STEP_KESTREL,
)

pytestmark = pytest.mark.unit


def _write_summary(directory: Path, summary: dict | None = None) -> Path:
    """Create a sample directory holding a `pipeline_summary.json`.

    Args:
        directory: The sample directory to create.
        summary: The summary to write; a minimal one by default.

    Returns:
        Path: The directory that was created.
    """
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "pipeline_summary.json").write_text(json.dumps(summary or {"version": "2.0.6"}), encoding="utf-8")
    return directory


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def test_a_directory_that_is_itself_a_sample_is_taken_as_one(tmp_path) -> None:
    sample = _write_summary(tmp_path / "sample_one")

    dirs, temp_dirs = discover_sample_directories([str(sample)])

    assert dirs == {sample}
    assert temp_dirs == []


def test_a_parent_directory_is_searched_recursively(tmp_path) -> None:
    one = _write_summary(tmp_path / "cohort" / "sample_one")
    two = _write_summary(tmp_path / "cohort" / "nested" / "sample_two")

    dirs, _ = discover_sample_directories([str(tmp_path / "cohort")])

    assert dirs == {one, two}


def test_a_directory_that_is_itself_a_sample_is_not_also_searched(tmp_path) -> None:
    """The recursive search is the `else` branch, so a sample directory containing
    further samples contributes only itself."""
    parent = _write_summary(tmp_path / "sample_one")
    _write_summary(parent / "inner")

    dirs, _ = discover_sample_directories([str(parent)])

    assert dirs == {parent}


def test_the_same_sample_reached_twice_is_only_processed_once(tmp_path) -> None:
    sample = _write_summary(tmp_path / "cohort" / "sample_one")

    dirs, _ = discover_sample_directories([str(sample), str(tmp_path / "cohort")])

    assert dirs == {sample}


def test_a_missing_input_path_is_skipped_with_a_warning(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")

    dirs, _ = discover_sample_directories([str(tmp_path / "absent")])

    assert dirs == set()
    assert "does not exist" in caplog.text


def test_a_directory_holding_no_summary_at_all_is_reported(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    (tmp_path / "empty").mkdir()

    dirs, _ = discover_sample_directories([str(tmp_path / "empty")])

    assert dirs == set()
    assert "No pipeline_summary.json found in directory" in caplog.text


def test_a_file_that_is_neither_a_directory_nor_a_zip_is_reported(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    stray = tmp_path / "notes.txt"
    stray.write_text("not a cohort", encoding="utf-8")

    dirs, _ = discover_sample_directories([str(stray)])

    assert dirs == set()
    assert "Unsupported file type" in caplog.text


def test_discovery_returns_an_unordered_set_today() -> None:
    """Characterisation of a live defect: the cohort's row order is not reproducible.

    Sample directories are collected into a `set`, and `aggregate_cohort` iterates it
    directly. Set iteration order for `Path` objects follows their string hashes, which
    Python randomises per process, so two runs of `vntyper cohort` over the same inputs
    put the rows of the report - and of `cohort_kestrel.csv` and its siblings - in
    different orders. See
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-nondeterministic-sample-order.md`.

    Pinned as the type rather than as an order, because an order this test could assert
    is exactly what the code does not provide.
    """
    dirs, _ = discover_sample_directories([])

    assert isinstance(dirs, set)


# ---------------------------------------------------------------------------
# Zip inputs
# ---------------------------------------------------------------------------


def _zip_of(tmp_path: Path, name: str, members: dict[str, str]) -> Path:
    """Build a zip archive from a name -> contents mapping.

    Args:
        tmp_path: Where to write the archive.
        name: The archive file name.
        members: Archive-relative path -> file contents.

    Returns:
        Path: The archive.
    """
    archive = tmp_path / name
    with zipfile.ZipFile(archive, "w") as handle:
        for member, contents in members.items():
            handle.writestr(member, contents)
    return archive


def test_a_zip_whose_root_is_a_sample_is_extracted_and_used(tmp_path) -> None:
    archive = _zip_of(tmp_path, "sample.zip", {"pipeline_summary.json": '{"version": "2.0.6"}'})

    dirs, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert len(dirs) == 1
        assert len(temp_dirs) == 1
        assert next(iter(dirs)) == Path(temp_dirs[0])
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_zip_of_several_samples_is_searched_recursively(tmp_path) -> None:
    archive = _zip_of(
        tmp_path,
        "cohort.zip",
        {
            "sample_one/pipeline_summary.json": '{"version": "2.0.6"}',
            "sample_two/pipeline_summary.json": '{"version": "2.0.6"}',
        },
    )

    dirs, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert {d.name for d in dirs} == {"sample_one", "sample_two"}
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_zip_holding_no_summary_is_reported_and_its_temp_dir_still_tracked(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    archive = _zip_of(tmp_path, "empty.zip", {"readme.txt": "nothing here"})

    dirs, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert dirs == set()
        assert "No pipeline_summary.json found in extracted zip file" in caplog.text
        assert len(temp_dirs) == 1
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_corrupt_zip_is_reported_and_leaves_no_temporary_directory(tmp_path, caplog) -> None:
    """A damaged upload must not leave the extraction directory behind.

    The archive's central directory is left intact so `is_zipfile` still accepts it -
    that is the interesting case, because a file rejected outright never reaches the
    extraction path at all - and its first local file header is zeroed, so
    `extractall` raises `BadZipFile` partway through.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_inputs")
    archive = _zip_of(tmp_path, "broken.zip", {"pipeline_summary.json": "{}"})
    raw = bytearray(archive.read_bytes())
    raw[0:4] = b"\x00\x00\x00\x00"
    archive.write_bytes(bytes(raw))
    assert zipfile.is_zipfile(archive)

    dirs, temp_dirs = discover_sample_directories([str(archive)])

    assert dirs == set()
    assert temp_dirs == []
    assert "Bad zip file" in caplog.text


def test_a_zip_that_fails_to_extract_for_any_other_reason_is_also_cleaned_up(tmp_path, caplog, monkeypatch) -> None:
    """`BadZipFile` is not the only way extraction fails - a full disk or a refused
    permission arrives as something else, and the temporary directory has to go either
    way."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_inputs")
    archive = _zip_of(tmp_path, "ok.zip", {"pipeline_summary.json": "{}"})

    def _explode(*args, **kwargs):
        raise OSError("No space left on device")

    monkeypatch.setattr(zipfile.ZipFile, "extractall", _explode)

    dirs, temp_dirs = discover_sample_directories([str(archive)])

    assert dirs == set()
    assert temp_dirs == []
    assert "Error extracting zip file" in caplog.text


def test_cleanup_removes_every_temporary_directory(tmp_path) -> None:
    first = tmp_path / "a"
    second = tmp_path / "b"
    first.mkdir()
    second.mkdir()

    cleanup_temp_dirs([str(first), str(second)])

    assert not first.exists()
    assert not second.exists()


def test_cleanup_reports_a_directory_it_cannot_remove_rather_than_raising(tmp_path, caplog) -> None:
    """One undeletable temporary directory must not abort the run after the report has
    already been written."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_inputs")

    cleanup_temp_dirs([str(tmp_path / "never-existed")])

    assert "Failed to remove temporary directory" in caplog.text


# ---------------------------------------------------------------------------
# Reading one sample's summary
# ---------------------------------------------------------------------------


def test_the_four_consumed_steps_are_read_out_of_a_summary() -> None:
    summary = {
        "pipeline_start": "2026-01-01T00:00:00",
        "pipeline_end": "2026-01-01T00:01:30",
        "version": "2.0.6",
        "steps": [
            {"step": STEP_KESTREL, "parsed_result": {"data": [{"Motif": "5"}]}},
            {"step": STEP_ADVNTR, "parsed_result": {"data": [{"VID": "25561"}]}},
            {"step": STEP_BAM_HEADER, "parsed_result": {"assembly_text": "hg38", "alignment_pipeline": "bwa-mem"}},
            {"step": STEP_COVERAGE, "parsed_result": {"data": [{"mean": "31.2"}]}},
        ],
    }

    kestrel, advntr, stats = parse_pipeline_summary(summary)

    assert kestrel == [{"Motif": "5"}]
    assert advntr == [{"VID": "25561"}]
    assert stats == {
        "runtime": "90.00 seconds",
        "version": "2.0.6",
        "assembly": "hg38",
        "pipeline": "bwa-mem",
        "coverage": {"mean": "31.2"},
    }


def test_an_empty_summary_yields_the_documented_defaults() -> None:
    kestrel, advntr, stats = parse_pipeline_summary({})

    assert kestrel == []
    assert advntr == []
    assert stats == {"runtime": "N/A", "version": "N/A", "assembly": "N/A", "pipeline": "N/A", "coverage": {}}


@pytest.mark.parametrize(
    "start,end",
    [(None, None), ("2026-01-01T00:00:00", None), (None, "2026-01-01T00:01:30"), ("", "")],
)
def test_a_run_missing_either_timestamp_reports_no_runtime(start: str | None, end: str | None) -> None:
    _, _, stats = parse_pipeline_summary({"pipeline_start": start, "pipeline_end": end})

    assert stats["runtime"] == "N/A"


def test_the_runtime_is_rendered_to_two_decimal_places() -> None:
    summary = {"pipeline_start": "2026-01-01T00:00:00", "pipeline_end": "2026-01-01T00:00:01.005"}

    _, _, stats = parse_pipeline_summary(summary)

    assert stats["runtime"] == "1.00 seconds"


def test_a_coverage_step_with_no_rows_leaves_coverage_empty() -> None:
    summary = {"steps": [{"step": STEP_COVERAGE, "parsed_result": {"data": []}}]}

    _, _, stats = parse_pipeline_summary(summary)

    assert stats["coverage"] == {}


def test_only_the_first_coverage_row_is_used() -> None:
    summary = {"steps": [{"step": STEP_COVERAGE, "parsed_result": {"data": [{"mean": "1"}, {"mean": "2"}]}}]}

    _, _, stats = parse_pipeline_summary(summary)

    assert stats["coverage"] == {"mean": "1"}


def test_a_step_name_the_cohort_does_not_consume_is_ignored() -> None:
    summary = {"steps": [{"step": "SHARK Filtering", "parsed_result": {"data": [{"x": "1"}]}}]}

    kestrel, advntr, stats = parse_pipeline_summary(summary)

    assert kestrel == []
    assert advntr == []
    assert stats["assembly"] == "N/A"


def test_a_later_step_of_the_same_name_wins() -> None:
    """`pipeline.py` records `BAM Header Parsing` more than once for a re-aligned run;
    the loop keeps overwriting, so the last one recorded is what the report shows."""
    summary = {
        "steps": [
            {"step": STEP_BAM_HEADER, "parsed_result": {"assembly_text": "hg19", "alignment_pipeline": "bwa"}},
            {"step": STEP_BAM_HEADER, "parsed_result": {"assembly_text": "hg38", "alignment_pipeline": "bwa-mem"}},
        ]
    }

    _, _, stats = parse_pipeline_summary(summary)

    assert stats["assembly"] == "hg38"


# ---------------------------------------------------------------------------
# Reading one sample's summary off disk
# ---------------------------------------------------------------------------


def test_a_sample_directory_is_read_from_its_summary_file(tmp_path) -> None:
    sample = _write_summary(
        tmp_path / "sample_one",
        {"version": "2.0.6", "steps": [{"step": STEP_KESTREL, "parsed_result": {"data": [{"Motif": "5"}]}}]},
    )

    kestrel, advntr, stats = load_pipeline_summary_for_sample(sample)

    assert kestrel == [{"Motif": "5"}]
    assert advntr == []
    assert stats["version"] == "2.0.6"


def test_a_sample_directory_with_no_summary_yields_three_empties(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    (tmp_path / "sample_one").mkdir()

    assert load_pipeline_summary_for_sample(tmp_path / "sample_one") == ([], [], {})
    assert "Pipeline summary file not found" in caplog.text


def test_an_unreadable_summary_yields_three_empties_rather_than_aborting_the_cohort(tmp_path, caplog) -> None:
    """One corrupt sample must not take the other 40 down with it."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_inputs")
    sample = tmp_path / "sample_one"
    sample.mkdir()
    (sample / "pipeline_summary.json").write_text("{not json", encoding="utf-8")

    assert load_pipeline_summary_for_sample(sample) == ([], [], {})
    assert "Error loading pipeline summary" in caplog.text


def test_a_malformed_timestamp_yields_three_empties_today(tmp_path, caplog) -> None:
    """Characterisation, not intent: an unparseable timestamp discards the sample.

    `datetime.fromisoformat` raises, and the handler that catches it wraps the whole
    read - so a sample whose `pipeline_start` is malformed loses its Kestrel and adVNTR
    rows too, not just its runtime. It is logged at error level and the cohort carries
    on without it.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_inputs")
    sample = _write_summary(
        tmp_path / "sample_one",
        {
            "pipeline_start": "not-a-timestamp",
            "pipeline_end": "2026-01-01T00:01:30",
            "steps": [{"step": STEP_KESTREL, "parsed_result": {"data": [{"Motif": "5"}]}}],
        },
    )

    assert load_pipeline_summary_for_sample(sample) == ([], [], {})
    assert "Error loading pipeline summary" in caplog.text


def test_the_sample_directory_may_be_given_as_a_string(tmp_path) -> None:
    sample = _write_summary(tmp_path / "sample_one")

    assert load_pipeline_summary_for_sample(str(sample))[2]["version"] == "2.0.6"


# ---------------------------------------------------------------------------
# Pseudonyms
# ---------------------------------------------------------------------------


def test_a_pseudonym_is_the_prefix_and_five_hex_digits() -> None:
    assert pseudonymized_sample_name("anon_", "sample_one") == "anon_65622"


def test_the_same_sample_name_always_gets_the_same_pseudonym() -> None:
    """The mapping has to be stable so a cohort re-run stays comparable to its
    predecessor and to the pseudonymization table written beside it."""
    assert pseudonymized_sample_name("x", "s1") == pseudonymized_sample_name("x", "s1")


def test_different_sample_names_get_different_pseudonyms() -> None:
    assert pseudonymized_sample_name("x", "s1") != pseudonymized_sample_name("x", "s2")


def test_the_prefix_may_be_any_value_the_cli_accepted() -> None:
    """`--pseudonymize` takes a string, and the CLI has passed `True` through in the
    past; the prefix is interpolated rather than concatenated so neither raises."""
    assert pseudonymized_sample_name(True, "s1").startswith("True")
