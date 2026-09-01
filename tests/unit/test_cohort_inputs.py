"""Finding the samples of a cohort on disk, and reading what the pipeline recorded.

A cohort is whatever `vntyper cohort` was pointed at: directories that are themselves a
sample, directories containing many, zip files of either, and anything that turns out to
be none of those. This module locates the `pipeline_summary.json` files in that, and
reads each one into the three pieces the report needs - the Kestrel rows, the adVNTR
rows, and the per-sample statistics.

This was 120 lines of `cohort_summary.py` reachable only by calling the whole cohort
pipeline, and it was the largest single uncovered block in the file. Everything here is
**characterisation** - it records what a cohort run does today - except the discovery-order
and discovery-identity tests, which are specifications and say so in their own docstrings.
No other test in this file has been ratified.

The two `..._today` tests in the zip section, which characterised a zip-rooted sample's
random identity and its irreproducible order as defects rather than guarantees, are gone:
#205 fixed both, and they now pin the fixed behaviour under names that say so. The wider
specification of a sample's identity is next door, in the four modules
`test_cohort_identity.py` became at 1,210 lines: the qualification rule and the injectivity
of the pseudonym map there, a ZIP sample's identity and order in
`test_cohort_zip_identity.py`, the digest and the configuration that chooses it in
`test_cohort_pseudonym_config.py`, and one-physical-sample-is-one-record in
`test_cohort_deduplication.py`. The pseudonym itself left this module with
`pseudonymized_sample_name` and is characterised in `test_cohort_pseudonyms.py`.

Step names are matched by exact string comparison against what `pipeline.py` writes
(AGENTS.md trap 5), so this file asserts against `summary_steps.STEP_*` for the same
reason the module does.
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import pytest

from vntyper.scripts import cohort_inputs
from vntyper.scripts.cohort_inputs import (
    cleanup_temp_dirs,
    discover_sample_directories,
    load_pipeline_summary_for_sample,
    parse_pipeline_summary,
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

    assert [found.directory for found in dirs] == [sample]
    assert [found.identity for found in dirs] == ["sample_one"]
    assert temp_dirs == []


def test_a_parent_directory_is_searched_recursively(tmp_path) -> None:
    """`nested/sample_two` sorts before `sample_one`, because the sort is by path part
    and `"nested" < "sample_one"` - a deeper sample is not pushed to the end."""
    one = _write_summary(tmp_path / "cohort" / "sample_one")
    two = _write_summary(tmp_path / "cohort" / "nested" / "sample_two")

    dirs, _ = discover_sample_directories([str(tmp_path / "cohort")])

    assert [found.directory for found in dirs] == [two, one]


def test_a_directory_that_is_itself_a_sample_is_not_also_searched(tmp_path) -> None:
    """The recursive search is the `else` branch, so a sample directory containing
    further samples contributes only itself."""
    parent = _write_summary(tmp_path / "sample_one")
    _write_summary(parent / "inner")

    dirs, _ = discover_sample_directories([str(parent)])

    assert [found.directory for found in dirs] == [parent]


def test_a_missing_input_path_is_skipped_with_a_warning(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")

    dirs, _ = discover_sample_directories([str(tmp_path / "absent")])

    assert dirs == []
    assert "does not exist" in caplog.text


def test_a_directory_holding_no_summary_at_all_is_reported(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    (tmp_path / "empty").mkdir()

    dirs, _ = discover_sample_directories([str(tmp_path / "empty")])

    assert dirs == []
    assert "No pipeline_summary.json found in directory" in caplog.text


def test_a_file_that_is_neither_a_directory_nor_a_zip_is_reported(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    stray = tmp_path / "notes.txt"
    stray.write_text("not a cohort", encoding="utf-8")

    dirs, _ = discover_sample_directories([str(stray)])

    assert dirs == []
    assert "Unsupported file type" in caplog.text


def test_the_discovered_directories_come_back_sorted(tmp_path) -> None:
    """**Specification**: discovery has an order, and it is the lexicographic one.

    The directories used to come back in a `set`, which `aggregate_cohort` iterated
    directly, so the row order of the report and of `cohort_kestrel.csv` and its
    siblings varied between processes. Sorting happens here rather than at the one call
    site so that every consumer of this function gets the same order, not just the one
    that was noticed. See
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-nondeterministic-sample-order.md`.

    The inputs are named in reverse on purpose: a deliberately adverse ordering is what
    distinguishes "sorted" from "happens to arrive in order".

    **This expectation survived #205 unchanged, and that was a deliberate choice.** The
    sort key stopped being the sample's *extracted* path - whose leading component for a
    zip is `tempfile.mkdtemp`'s random suffix, so zip samples had no reproducible position
    at all - and became `(parts of the input path, path relative to that input's root)`.
    An earlier draft keyed the outer half on the input's *position* in `input_paths`
    instead, which is equally reproducible but would have reordered the rows of every
    existing directory-input cohort report for no reason #205 required.
    """
    first = _write_summary(tmp_path / "cohort" / "sample_a")
    second = _write_summary(tmp_path / "cohort" / "sample_b")
    third = _write_summary(tmp_path / "cohort" / "sample_c")

    dirs, _ = discover_sample_directories([str(third), str(second), str(first)])

    assert [found.directory for found in dirs] == [first, second, third]


def test_the_order_is_lexicographic_by_path_part_rather_than_by_raw_string(tmp_path) -> None:
    """**Specification**: the sort compares path *parts*, not raw path strings.

    The separator therefore never takes part in the comparison. That is the difference
    this pins: `cohort/sample_one` sorts before `cohort-extra/sample_one` because
    `"cohort" < "cohort-extra"`, where the raw strings compare the other way round (`"-"`
    is 0x2d and `"/"` is 0x2f). It is also the order a user reading `ls` would predict,
    which is the point of choosing it.

    It held when the key was `sorted()` over `Path` (`PurePath.__lt__` compares
    `_parts_normcase`) and it holds now that the key is a tuple of parts compared
    element-wise, for the same reason and with the same result - here across two inputs,
    and in `test_cohort_zip_identity.py` within one.
    """
    plain = _write_summary(tmp_path / "cohort" / "sample_one")
    suffixed = _write_summary(tmp_path / "cohort-extra" / "sample_one")
    assert str(suffixed) < str(plain)

    dirs, _ = discover_sample_directories([str(suffixed), str(plain)])

    assert [found.directory for found in dirs] == [plain, suffixed]


def test_an_input_nested_inside_a_later_input_keeps_its_whole_path_position(tmp_path) -> None:
    """**Specification**: the order is the whole-path order, including when inputs nest.

    The nested case is the one a two-part key gets wrong. Given `cohort/a` as a direct
    input and then `cohort` - which holds both `a` and `z` - discovery de-duplicates on
    the sample directory, so `cohort/a` keeps the record the *first* input gave it. Under
    a `(origin parts, relative parts)` key that record's outer half is
    `(..., "cohort", "a")` while `z`'s is `(..., "cohort")`, and tuple comparison decides
    it on the outer half alone: the parent's samples sort before the child's, giving
    `z, a`.

    Sorting on `Path` - what this did before the key was introduced - gives `a, z`, and so
    does `ls`. A single flattened effective path, the input's parts followed by the
    sample's parts below that input, restores it: `a` keys on
    `(..., "cohort", "a")` whichever input claimed it, because the two compose to the same
    tuple.

    Reaching a directory both directly and through its parent is not a contrived input:
    `vntyper cohort run_1 run_1/../` is one shell expansion away, and the web service
    passes a caller-supplied list straight through.
    """
    child = _write_summary(tmp_path / "cohort" / "a")
    sibling = _write_summary(tmp_path / "cohort" / "z")

    dirs, _ = discover_sample_directories([str(child), str(tmp_path / "cohort")])

    assert [found.directory for found in dirs] == [child, sibling]
    # The same set of samples, named the other way round, has to agree.
    reversed_inputs, _ = discover_sample_directories([str(tmp_path / "cohort"), str(child)])
    assert [found.directory for found in reversed_inputs] == [child, sibling]


#: Run in a subprocess: report the identities discovery finds under ``argv[1]``.
#:
#: ``identity`` rather than ``directory.name``: the identity is what the report actually
#: shows, so the cross-process test pins the value that matters rather than a path
#: component that no longer decides it (#205).
_ORDER_PROBE = (
    "import json, sys;"
    "from vntyper.scripts.cohort_inputs import discover_sample_directories;"
    "print(json.dumps([s.identity for s in discover_sample_directories(sys.argv[1:])[0]]))"
)


def _discovery_order_under_hash_seed(root: Path, seed: str) -> list[str]:
    """Discover ``root``'s samples in a fresh interpreter with a chosen hash seed.

    ``PYTHONPATH`` is seeded from this process's ``sys.path`` so the child imports the
    same working tree's ``vntyper`` regardless of its working directory.

    Args:
        root: The cohort directory to discover.
        seed: The ``PYTHONHASHSEED`` value for the child interpreter.

    Returns:
        list[str]: The discovered sample directory names, in discovery order.
    """
    result = subprocess.run(
        [sys.executable, "-c", _ORDER_PROBE, str(root)],
        env={**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": os.pathsep.join(sys.path)},
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


def test_processes_with_different_hash_seeds_discover_the_same_order(tmp_path) -> None:
    """**Specification**, and the defect this closes is *between* processes.

    A single-process comparison cannot see the old behaviour at all: `PYTHONHASHSEED` is
    fixed for the life of an interpreter, so two calls in one process agreed even while
    two `vntyper cohort` runs did not. The directories came back in a `set`, set
    iteration order for `Path` follows `Path.__hash__`, and that is the hash of the path
    string, which Python randomises per process. Five interpreters with five seeds is
    the only way to observe it, so this spawns five - ten samples each, because a
    two-element set can agree by luck and a ten-element one cannot.

    **Scope: directory inputs.** The ten samples are directories under `tmp_path`, whose
    paths are the same in all five children, so the sort is over a fixed set of strings
    and is total. Zip inputs used to be outside that scope entirely - their sort key ended
    in `tempfile.mkdtemp`'s random suffix, so it differed between processes whatever the
    hash seed was - and are now covered by their own cross-process test,
    `test_cohort_zip_identity.py::test_zip_inputs_come_back_in_the_same_order_in_five_processes`.

    The claim this test does **not** support, and used to: that the stable fingerprint in
    `test_cohort_summary_oracle.py` is evidence for the same fix. That fingerprint is
    taken from `generate_cohort_summary_report`, which never calls this function, so it
    was stable across hash seeds before the fix too.
    """
    for index in range(10):
        _write_summary(tmp_path / "cohort" / f"sample_{index:02d}")

    orders = [_discovery_order_under_hash_seed(tmp_path / "cohort", seed) for seed in ("0", "1", "2", "3", "4")]

    assert len({tuple(order) for order in orders}) == 1
    assert orders[0] == [f"sample_{index:02d}" for index in range(10)]


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
        assert len(temp_dirs) == 1
        assert [found.directory for found in dirs] == [Path(temp_dirs[0])]
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
        assert [found.directory.name for found in dirs] == ["sample_one", "sample_two"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_zip_rooted_sample_is_identified_by_the_archive_rather_than_by_the_temp_dir(tmp_path) -> None:
    """**Specification** (#205); this replaces the characterisation of the same layout.

    When `pipeline_summary.json` sits at the root of the archive - which is the layout the
    web worker produces, so the normal path for web cohorts - the discovered sample
    directory *is* the extraction directory, and that is
    `tempfile.mkdtemp(prefix="cohort_zip_")`. `aggregate_cohort` took the sample's
    identity from `Path(sample_dir).name`, so the sample was reported, exported and
    pseudonymised under a random `cohort_zip_XXXXXXXX` string that was different on every
    run and carried nothing of the archive's own name.

    The identity is now carried out of discovery explicitly. It is the stem of the input
    file the run itself recorded, and the archive's own stem when the summary records
    none - which is this fixture, whose summary predates `input_files`. The full identity
    specification, including which recorded input wins, is in `test_cohort_zip_identity.py`.
    """
    archive = _zip_of(tmp_path, "patient_one.zip", {"pipeline_summary.json": '{"version": "2.0.6"}'})

    dirs, temp_dirs = discover_sample_directories([str(archive)])

    try:
        (sample,) = dirs
        assert sample.identity == "patient_one"
        assert sample.directory.name.startswith("cohort_zip_")
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_two_zip_inputs_are_ordered_by_their_archive_paths_not_by_their_temp_dirs(tmp_path, monkeypatch) -> None:
    """**Specification** (#205); this replaces the characterisation of the same defect.

    `discover_sample_directories` sorts, so directory inputs came back in an order two
    processes agreed on - that is what
    `test_processes_with_different_hash_seeds_discover_the_same_order` pins. The sort key
    for a zip-rooted sample used to be the extraction directory's full path, and its last
    component is `tempfile.mkdtemp`'s random suffix, so for zips the sort was total but
    the thing being sorted was different every run. Two archives came back in a different
    order each time, and the report's row order followed. The outer half of the key is now
    the archive path the caller wrote, which is the zip's analogue of the directory path
    the directory-input tests above sort on.

    `mkdtemp` is stubbed to hand out chosen suffixes rather than random ones, because with
    random ones this test would pass about half the time on unfixed code and would
    therefore prove nothing. The archive given **first** is deliberately the one that
    extracts to the suffix sorting **last**, which is exactly the arrangement the old sort
    got wrong.
    """
    first = _zip_of(tmp_path, "aaa_cohort.zip", {"pipeline_summary.json": '{"version": "AAA"}'})
    second = _zip_of(tmp_path, "zzz_cohort.zip", {"pipeline_summary.json": '{"version": "ZZZ"}'})
    suffixes = iter(["zzzzzzzz", "aaaaaaaa"])

    def _mkdtemp(prefix: str = "") -> str:
        directory = tmp_path / "extracted" / f"{prefix}{next(suffixes)}"
        directory.mkdir(parents=True)
        return str(directory)

    monkeypatch.setattr(tempfile, "mkdtemp", _mkdtemp)

    dirs, temp_dirs = discover_sample_directories([str(first), str(second)])

    try:
        assert [found.identity for found in dirs] == ["aaa_cohort", "zzz_cohort"]
        assert [found.directory.name for found in dirs] == ["cohort_zip_zzzzzzzz", "cohort_zip_aaaaaaaa"]
        # And with the order goes the data: the first archive's sample is reported first.
        assert load_pipeline_summary_for_sample(dirs[0].directory)[2]["version"] == "AAA"
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_zip_holding_no_summary_is_reported_and_its_temp_dir_still_tracked(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    archive = _zip_of(tmp_path, "empty.zip", {"readme.txt": "nothing here"})

    dirs, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert dirs == []
        assert "No pipeline_summary.json found in extracted zip file" in caplog.text
        assert len(temp_dirs) == 1
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_bad_archive_is_skipped_and_other_inputs_continue(tmp_path, caplog, monkeypatch) -> None:
    """A malformed archive does not prevent discovery of later valid inputs."""
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    broken = _zip_of(tmp_path, "broken.zip", {"pipeline_summary.json": "{}"})
    valid = _write_summary(tmp_path / "valid")
    original_extractall = zipfile.ZipFile.extractall

    def _blocked_extractall(self, *args, **kwargs):
        if self.filename == str(broken):
            raise zipfile.BadZipFile("blocked")
        return original_extractall(self, *args, **kwargs)

    monkeypatch.setattr(zipfile.ZipFile, "extractall", _blocked_extractall)

    dirs, temp_dirs = discover_sample_directories([str(broken), str(valid)])

    assert [sample.identity for sample in dirs] == ["valid"]
    assert temp_dirs == []
    records = [record for record in caplog.records if record.name == "vntyper.scripts.cohort_inputs"]
    assert len(records) == 1
    assert records[0].levelno == logging.ERROR
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

    assert dirs == []
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


def test_cleanup_attempts_all_directories_after_failure(tmp_path, caplog, monkeypatch) -> None:
    """Cleanup keeps attempting later directories after one filesystem failure."""
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    first = str(tmp_path / "first")
    second = str(tmp_path / "second")
    attempted: list[str] = []

    def _blocked_rmtree(path: str) -> None:
        attempted.append(path)
        if path == first:
            raise OSError("blocked")

    monkeypatch.setattr(cohort_inputs.shutil, "rmtree", _blocked_rmtree)

    cleanup_temp_dirs([first, second])
    assert attempted == [first, second]
    records = [record for record in caplog.records if record.name == "vntyper.scripts.cohort_inputs"]
    assert len(records) == 1
    assert records[0].levelno == logging.ERROR
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

    legacy_identity = {
        "Molecular_Identity": "legacy identity not recorded",
        "Molecular_Identity_Status": "legacy identity not recorded",
        "Equivalent_Representation_Count": "legacy identity not recorded",
        "Identity_Hypothesis_Count": "legacy identity not recorded",
    }
    assert kestrel == [{"Motif": "5", **legacy_identity}]
    assert advntr == [{"VID": "25561", **legacy_identity}]
    assert stats == {
        "runtime": "90.00 seconds",
        "version": "2.0.6",
        "assembly": "hg38",
        "pipeline": "bwa-mem",
        "coverage": {"mean": "31.2"},
    }


@pytest.mark.parametrize("schema_version", [1, 2])
def test_legacy_positive_rows_reach_every_cohort_surface_with_four_explicit_legacy_cells(
    schema_version: int,
) -> None:
    """Catch cohort reconstruction from plausible legacy allele and nomenclature cells."""
    legacy = {"Motif": "5", "POS": 67, "REF": "G", "ALT": "GG", "Nomenclature": "59dupC"}
    summary = {
        "schema_version": schema_version,
        "steps": [{"step": STEP_KESTREL, "parsed_result": {"data": [legacy]}}],
    }

    kestrel, _, _ = parse_pipeline_summary(summary)

    quartet = (
        "Molecular_Identity",
        "Molecular_Identity_Status",
        "Equivalent_Representation_Count",
        "Identity_Hypothesis_Count",
    )
    assert list(kestrel[0])[-4:] == list(quartet)
    assert {key: kestrel[0][key] for key in quartet} == dict.fromkeys(quartet, "legacy identity not recorded")


def test_complete_recorded_identity_values_reach_the_cohort_unchanged() -> None:
    """Catch cohort-side status/count recalculation from the other cells in the row."""
    recorded = {
        "Molecular_Identity": "",
        "Molecular_Identity_Status": "unresolved",
        "Equivalent_Representation_Count": 0,
        "Identity_Hypothesis_Count": 5,
    }
    summary = {
        "schema_version": 2,
        "steps": [
            {
                "step": STEP_ADVNTR,
                "parsed_result": {"data": [{"VID": "25561", "POS": 67, "REF": "G", "ALT": "GG", **recorded}]},
            }
        ],
    }

    _, advntr, _ = parse_pipeline_summary(summary)

    assert {key: advntr[0][key] for key in recorded} == recorded


def test_negative_caller_rows_are_not_widened_during_cohort_loading() -> None:
    """Catch compatibility projection changing either frozen negative result schema."""
    kestrel_negative = {"Motif": "None", "POS": "None", "REF": "None", "ALT": "None", "Confidence": "Negative"}
    advntr_negative = {"VID": "Negative", "Variant": "None", "Flag": "Not applicable"}
    summary = {
        "schema_version": 2,
        "steps": [
            {"step": STEP_KESTREL, "parsed_result": {"data": [kestrel_negative]}},
            {"step": STEP_ADVNTR, "parsed_result": {"data": [advntr_negative]}},
        ],
    }

    kestrel, advntr, _ = parse_pipeline_summary(summary)

    assert kestrel == [kestrel_negative]
    assert advntr == [advntr_negative]


def test_an_empty_summary_yields_the_documented_defaults() -> None:
    kestrel, advntr, stats = parse_pipeline_summary({})

    assert kestrel == []
    assert advntr == []
    assert stats == {"runtime": "N/A", "version": "N/A", "assembly": "N/A", "pipeline": "N/A", "coverage": {}}


def test_cohort_loader_verifies_recorded_advntr_evidence_from_each_sample_snapshot(tmp_path: Path) -> None:
    from vntyper.modules.advntr.artifact_evidence import ASSERTION, load_packaged_artifact_evidence

    evidence = load_packaged_artifact_evidence()
    sample = _write_summary(tmp_path / "sample", {"version": "2.0.6", "advntr_evidence_digest": evidence.digest})
    snapshot = sample / "provenance" / "advntr_artifact_evidence.json"
    snapshot.parent.mkdir()
    snapshot.write_bytes(evidence.canonical_bytes)

    stats = load_pipeline_summary_for_sample(sample)[2]

    assert stats["advntr_evidence_revision"] == evidence.digest
    assert stats["advntr_evidence_assertion"] == ASSERTION


def test_cohort_loader_labels_legacy_advntr_evidence_without_reading_current_package(tmp_path: Path) -> None:
    sample = _write_summary(tmp_path / "sample", {"version": "2.0.6"})

    stats = load_pipeline_summary_for_sample(sample)[2]

    assert stats["advntr_evidence_revision"] == "artifact-evidence revision not recorded"
    assert stats["advntr_evidence_assertion"] == ""


def test_cohort_loader_rejects_a_recorded_digest_without_its_run_snapshot(tmp_path: Path) -> None:
    sample = _write_summary(tmp_path / "sample", {"version": "2.0.6", "advntr_evidence_digest": "0" * 64})

    with pytest.raises(ValueError, match="Invalid decision profile or evidence provenance"):
        load_pipeline_summary_for_sample(sample)


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

    assert kestrel == [
        {
            "Motif": "5",
            "Molecular_Identity": "legacy identity not recorded",
            "Molecular_Identity_Status": "legacy identity not recorded",
            "Equivalent_Representation_Count": "legacy identity not recorded",
            "Identity_Hypothesis_Count": "legacy identity not recorded",
            "Decision_Profile_ID": "decision profile not recorded by legacy run",
            "Decision_Profile_Revision": "decision profile not recorded by legacy run",
            "Decision_Profile_SHA256": "decision profile not recorded by legacy run",
        }
    ]
    assert advntr == []
    assert stats["version"] == "2.0.6"


def test_a_sample_directory_with_no_summary_yields_three_empties(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    (tmp_path / "sample_one").mkdir()

    assert load_pipeline_summary_for_sample(tmp_path / "sample_one") == ([], [], {})
    assert "Pipeline summary file not found" in caplog.text


def test_summary_read_failure_returns_three_empty_results(tmp_path, caplog, monkeypatch) -> None:
    """One unreadable summary leaves a structured empty result for that sample."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_inputs")
    sample = tmp_path / "sample_one"
    sample.mkdir()
    (sample / "pipeline_summary.json").write_text("{}", encoding="utf-8")

    def _blocked_open(*args, **kwargs):
        raise OSError("blocked")

    monkeypatch.setattr(cohort_inputs, "open", _blocked_open, raising=False)

    assert load_pipeline_summary_for_sample(sample) == ([], [], {})
    records = [record for record in caplog.records if record.name == "vntyper.scripts.cohort_inputs"]
    assert len(records) == 1
    assert records[0].levelno == logging.ERROR
    assert "Error loading pipeline summary" in caplog.text


def test_identity_read_failure_uses_directory_fallback(tmp_path, caplog, monkeypatch) -> None:
    """A root archive sample keeps the supplied fallback when its summary is unreadable."""
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_inputs")
    sample = tmp_path / "sample_one"
    sample.mkdir()
    (sample / "pipeline_summary.json").write_text("{}", encoding="utf-8")

    def _blocked_open(*args, **kwargs):
        raise OSError("blocked")

    monkeypatch.setattr(cohort_inputs, "open", _blocked_open, raising=False)

    assert cohort_inputs._identity_from_summary(sample, "archive_fallback") == "archive_fallback"
    records = [record for record in caplog.records if record.name == "vntyper.scripts.cohort_inputs"]
    assert len(records) == 1
    assert records[0].levelno == logging.WARNING
    assert "Could not read the recorded inputs" in caplog.text


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
