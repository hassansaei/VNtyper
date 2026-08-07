"""A ZIP sample's identity and its place in the cohort's order (#205).

**SPECIFICATION.** A ZIP is extracted into ``tempfile.mkdtemp(prefix="cohort_zip_")``. When
its ``pipeline_summary.json`` sits at the archive root - the layout the web worker
produces, so the normal path for web cohorts - that temporary directory *was* the sample
directory, and the reported ``Sample`` was therefore a random ``cohort_zip_*`` string. The
same random component was in the sort key, so two runs over the same archives ordered their
rows differently. Neither depends on the extraction any more: the identity comes from the
input file the run itself recorded, and the sort key is the sample's *effective path* - the
parts of the input path the caller wrote, followed by the sample's path relative to that
input's root, as one flat tuple. (It was briefly a nested pair of those two tuples, which
decided the comparison on the input path alone and therefore reordered a cohort whose
inputs nest; see
``test_cohort_inputs.py::test_an_input_nested_inside_a_later_input_keeps_its_whole_path_position``.)

Extraction itself is unchanged, and so is who owns the directories it produces: discovery
removes what it extracted if it cannot return the handles, and ``aggregate_cohort`` removes
them once it has. Both halves of that are pinned at the end of this module.

Split out of ``test_cohort_identity.py``, which reached 1,210 lines and now holds the
qualification rule alone; the digest is in ``test_cohort_pseudonym_config.py`` and
de-duplication in ``test_cohort_deduplication.py``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from tests.support.cohort_samples import sample_on_disk, summary_json, zip_of
from vntyper.cli import load_config
from vntyper.scripts import cohort_summary
from vntyper.scripts.cohort_inputs import (
    cleanup_temp_dirs,
    discover_sample_directories,
)

pytestmark = pytest.mark.unit


def test_a_root_level_zip_sample_is_identified_by_the_run_s_own_input(tmp_path) -> None:
    """#205: this layout is what the web worker produces, so it is the normal path.

    ``processed_dirs.add(temp_path)`` added the ``tempfile.mkdtemp`` root itself, and
    ``cohort_summary`` took ``Path(sample_dir).name`` as the sample - so the reported
    ``Sample`` was a ``cohort_zip_*`` string that differed on every run, and any pseudonym
    derived from it differed too.
    """
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["patient1"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_root_level_zip_without_input_files_falls_back_to_the_archive_stem(tmp_path) -> None:
    """Older summaries recorded no ``input_files``; the archive's own name is the next best thing."""
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json({})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["job7"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_root_level_zip_whose_summary_is_unreadable_falls_back_to_the_archive_stem(tmp_path) -> None:
    """One bad sample must not abort the cohort, so the identity helper never raises."""
    archive = zip_of(tmp_path / "job8.zip", {"pipeline_summary.json": "{not json"})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["job8"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_cram_run_is_identified_by_its_cram(tmp_path) -> None:
    """``bam`` first, then ``cram``, then ``fastq1`` - the order ``pipeline.py`` writes them."""
    archive = zip_of(tmp_path / "job9.zip", {"pipeline_summary.json": summary_json({"cram": "patient2.cram"})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["patient2"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_fastq_run_is_identified_by_its_first_mate(tmp_path) -> None:
    """``patient3_R1.fastq.gz`` carries two suffixes, and one ``.stem`` strips only one."""
    archive = zip_of(
        tmp_path / "job10.zip",
        {
            "pipeline_summary.json": summary_json(
                {"fastq1": "patient3_R1.fastq.gz", "fastq2": "patient3_R2.fastq.gz"},
            )
        },
    )

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["patient3_R1"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_the_alignment_wins_over_the_reads_when_a_run_recorded_both(tmp_path) -> None:
    """A re-aligned run records the FASTQs it started from and the BAM it produced."""
    archive = zip_of(
        tmp_path / "job11.zip",
        {"pipeline_summary.json": summary_json({"fastq1": "reads_R1.fastq.gz", "bam": "patient4.bam"})},
    )

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["patient4"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_sample_in_a_subdirectory_keeps_its_directory_name(tmp_path) -> None:
    """Only the root-level case changes: elsewhere the directory name is meaningful."""
    archive = zip_of(
        tmp_path / "cohort.zip",
        {
            "sample_one/pipeline_summary.json": summary_json({"bam": "patient1.bam"}),
            "sample_two/pipeline_summary.json": summary_json({"bam": "patient2.bam"}),
        },
    )

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["sample_one", "sample_two"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_the_sort_key_contains_no_temporary_path_component(tmp_path) -> None:
    """The mechanism, stated directly: nothing in the key comes from ``mkdtemp``.

    The key is asserted as a value. For a root-level sample it is exactly the archive path
    the caller wrote - there is nothing below the archive root to append - so the whole key
    is reconstructible from the command line alone, which is what makes it reproducible.
    """
    first = zip_of(tmp_path / "aaa_cohort.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})
    second = zip_of(tmp_path / "zzz_cohort.zip", {"pipeline_summary.json": summary_json({"bam": "patient2.bam"})})

    samples, temp_dirs = discover_sample_directories([str(first), str(second)])

    try:
        assert [sample.sort_key for sample in samples] == [first.parts, second.parts]
        assert [sample.origin for sample in samples] == [str(first), str(second)]
        assert all("cohort_zip_" not in part for sample in samples for part in sample.sort_key)
        # The temporary directory is still where extraction goes; it is only out of the key.
        assert all("cohort_zip_" in str(sample.directory) for sample in samples)
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_two_zips_are_ordered_by_their_archive_path_whatever_mkdtemp_hands_out(tmp_path, monkeypatch) -> None:
    """#205, deterministically: the extraction roots are chosen to sort the wrong way.

    For ZIP inputs the irreproducibility is ``mkdtemp``'s random suffix sitting in the
    sort key, not ``PYTHONHASHSEED`` - so genuinely random suffixes would make this test
    pass by chance about half the time. ``mkdtemp`` is stubbed to hand the archive given
    **first**, and named first, the root that sorts **last**, which is exactly the
    arrangement the old sort got wrong, every time.
    """
    first = zip_of(tmp_path / "aaa_cohort.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})
    second = zip_of(tmp_path / "zzz_cohort.zip", {"pipeline_summary.json": summary_json({"bam": "patient2.bam"})})
    suffixes = iter(["zzzzzzzz", "aaaaaaaa"])

    def _mkdtemp(prefix: str = "") -> str:
        directory = tmp_path / "extracted" / f"{prefix}{next(suffixes)}"
        directory.mkdir(parents=True)
        return str(directory)

    monkeypatch.setattr(tempfile, "mkdtemp", _mkdtemp)

    samples, temp_dirs = discover_sample_directories([str(first), str(second)])

    try:
        assert [sample.identity for sample in samples] == ["patient1", "patient2"]
        assert samples[0].directory == tmp_path / "extracted" / "cohort_zip_zzzzzzzz"
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_two_zips_named_in_reverse_still_come_back_in_archive_path_order(tmp_path) -> None:
    """The outer key is the archive *path*, not the argument's position.

    Naming the archives in reverse is what tells the two apart: a positional key would
    return them in the order they were written on the command line, and this asserts it
    does not. The same rule puts a cohort of directories back in the lexicographic order
    it has always had - see
    ``test_cohort_inputs.py::test_the_discovered_directories_come_back_sorted``.
    """
    first = zip_of(tmp_path / "aaa_cohort.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})
    second = zip_of(tmp_path / "zzz_cohort.zip", {"pipeline_summary.json": summary_json({"bam": "patient2.bam"})})

    samples, temp_dirs = discover_sample_directories([str(second), str(first)])

    try:
        assert [sample.identity for sample in samples] == ["patient1", "patient2"]
    finally:
        cleanup_temp_dirs(temp_dirs)


#: Run in a subprocess: report the identities ``discover_sample_directories`` finds.
_IDENTITY_PROBE = (
    "import json, sys;"
    "from vntyper.scripts.cohort_inputs import cleanup_temp_dirs, discover_sample_directories;"
    "samples, temp_dirs = discover_sample_directories(sys.argv[1:]);"
    "print(json.dumps([s.identity for s in samples]));"
    "cleanup_temp_dirs(temp_dirs)"
)


def _identities_under_hash_seed(archives: list[Path], seed: str) -> list[str]:
    """Discover ``archives`` in a fresh interpreter with a chosen hash seed.

    Args:
        archives: The ZIP inputs, in command-line order.
        seed: The ``PYTHONHASHSEED`` value for the child interpreter.

    Returns:
        list[str]: The discovered identities, in discovery order.
    """
    result = subprocess.run(
        [sys.executable, "-c", _IDENTITY_PROBE, *(str(archive) for archive in archives)],
        env={**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": os.pathsep.join(sys.path)},
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


def test_zip_inputs_come_back_in_the_same_order_in_five_processes(tmp_path) -> None:
    """Three archives, five interpreters, real ``mkdtemp`` roots: one order.

    This is the property the stubbed test above establishes mechanically, measured
    end to end against genuinely random extraction roots. It is the ZIP counterpart of
    ``test_cohort_inputs.py::test_processes_with_different_hash_seeds_discover_the_same_order``,
    which covers directory inputs and could not see this case at all.
    """
    archives = [
        zip_of(tmp_path / f"job_{index}.zip", {"pipeline_summary.json": summary_json({"bam": f"patient{index}.bam"})})
        for index in (2, 0, 1)
    ]

    orders = [_identities_under_hash_seed(archives, seed) for seed in ("0", "1", "2", "3", "4")]

    assert len({tuple(order) for order in orders}) == 1
    assert orders[0] == ["patient0", "patient1", "patient2"]


def test_a_zip_and_a_directory_are_ordered_by_their_input_paths(tmp_path) -> None:
    """The outer key is the input path, so the two input kinds interleave predictably.

    ``cohort/`` sorts before ``job.zip`` because ``"cohort" < "job.zip"``, and the archive
    given first therefore comes back second - the same rule as for two directories, with
    an archive standing in for one of them.
    """
    directory = sample_on_disk(tmp_path / "cohort" / "sample_one")
    archive = zip_of(tmp_path / "job.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})

    samples, temp_dirs = discover_sample_directories([str(archive), str(directory)])

    try:
        assert [sample.identity for sample in samples] == ["sample_one", "patient1"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_root_level_zip_sample_reports_the_same_identity_on_two_runs(tmp_path) -> None:
    """The end of #205 stated at the boundary a user sees: the exports name the patient.

    Two ``aggregate_cohort`` runs over the same archive, into two output directories.
    Before the fix the ``Sample`` column carried the ``cohort_zip_*`` root, which differs
    on every run.

    **What is compared is exactly what is claimed.** The machine-readable exports and
    ``pseudonymization_table.tsv`` are compared byte for byte; the HTML deliberately is
    not, because it carries a report timestamp and Plotly's per-figure UUIDs, which is why
    ``test_cohort_summary_oracle.py`` has a ``_skeleton()`` normaliser at all. The
    cross-*process* half of the claim is
    ``test_zip_inputs_come_back_in_the_same_order_in_five_processes`` above: discovery is
    the only step a per-process hash seed could reach.
    """
    archive = zip_of(
        tmp_path / "job7.zip",
        {
            "pipeline_summary.json": json.dumps(
                {
                    "version": "2.0.6",
                    "input_files": {"bam": "patient1.bam"},
                    "steps": [
                        {
                            "step": "Kestrel Genotyping",
                            "parsed_result": {
                                "data": [{"Motif": "5", "Confidence": "High_Precision", "Flag": "Not flagged"}]
                            },
                        }
                    ],
                }
            )
        },
    )
    compared = ("cohort_kestrel.csv", "cohort_kestrel.tsv", "cohort_kestrel.json", "pseudonymization_table.tsv")
    exports: list[dict[str, bytes]] = []
    for run in ("first", "second"):
        output_dir = tmp_path / run
        output_dir.mkdir()
        cohort_summary.aggregate_cohort(
            input_paths=[str(archive)],
            output_dir=str(output_dir),
            summary_file="cohort_summary.html",
            config=load_config(None),
            additional_formats="csv,tsv,json",
            pseudonymize_samples="anon_",
        )
        exports.append({name: (output_dir / name).read_bytes() for name in compared})

    assert exports[0] == exports[1]
    assert b"anon_" + hashlib.sha256(b"patient1").hexdigest()[:12].encode() in exports[0]["cohort_kestrel.csv"]
    assert b"patient1" in exports[0]["pseudonymization_table.tsv"]
    assert not any(b"cohort_zip_" in payload for payload in exports[0].values())


def test_the_extraction_still_goes_to_a_temporary_directory(tmp_path) -> None:
    """Only the temp path's role in identity and order is gone; extraction is unchanged."""
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert len(temp_dirs) == 1
        assert samples[0].directory == Path(temp_dirs[0])
        assert Path(temp_dirs[0]).name.startswith("cohort_zip_")
    finally:
        cleanup_temp_dirs(temp_dirs)
    assert not Path(temp_dirs[0]).exists()


def test_a_directory_sample_sorts_by_its_path_below_the_input_root(tmp_path) -> None:
    """Within one input the key is the path's *parts*, so the separator never sorts.

    ``cohort/sample_one`` comes before ``cohort-extra/sample_one`` because
    ``"cohort" < "cohort-extra"``, where the raw strings compare the other way round
    (``"-"`` is 0x2d and ``"/"`` is 0x2f). It is the order ``ls`` would show, which is the
    point of choosing it, and it survives both the move from whole paths to
    input-relative ones and the flattening of the key that followed.

    For a directory input the flattened key reconstructs the sample's own path exactly,
    which is the property the nesting case in ``test_cohort_inputs.py`` turns on.
    """
    plain = sample_on_disk(tmp_path / "cohort" / "sample_one")
    suffixed = sample_on_disk(tmp_path / "cohort-extra" / "sample_one")
    assert str(suffixed) < str(plain)

    samples, _ = discover_sample_directories([str(tmp_path)])

    assert [sample.directory for sample in samples] == [plain, suffixed]
    assert [sample.sort_key for sample in samples] == [plain.parts, suffixed.parts]
    assert [sample.sort_key for sample in samples] == [
        (*tmp_path.parts, "cohort", "sample_one"),
        (*tmp_path.parts, "cohort-extra", "sample_one"),
    ]


def test_the_cohort_reports_the_zip_identity_rather_than_the_extraction_directory(tmp_path) -> None:
    """`aggregate_cohort` consumes the identity rather than re-deriving it from the path."""
    archive = zip_of(
        tmp_path / "job7.zip",
        {
            "pipeline_summary.json": json.dumps(
                {
                    "version": "2.0.6",
                    "input_files": {"bam": "patient1.bam"},
                    "steps": [
                        {
                            "step": "Kestrel Genotyping",
                            "parsed_result": {
                                "data": [{"Motif": "5", "Confidence": "High_Precision", "Flag": "Not flagged"}]
                            },
                        }
                    ],
                }
            )
        },
    )
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(archive)],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    expected = "anon_" + hashlib.sha256(b"patient1").hexdigest()[:12]
    assert f"{expected}\tpatient1" in table
    assert "cohort_zip_" not in table


def test_the_discovery_result_carries_a_frozen_record_per_sample(tmp_path) -> None:
    """`DiscoveredSample` is a value, so no consumer can rename a sample in place."""
    sample = sample_on_disk(tmp_path / "cohort" / "sample_one")

    samples, _ = discover_sample_directories([str(sample)])

    (found,) = samples
    assert found.directory == sample
    assert found.identity == "sample_one"
    assert found.origin == str(sample)
    assert found.sort_key == sample.parts
    with pytest.raises(dataclasses.FrozenInstanceError):
        found.identity = "renamed"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# extraction cleanup
# ---------------------------------------------------------------------------


def test_a_failure_between_discovery_and_the_sample_loop_still_removes_the_extractions(tmp_path, monkeypatch) -> None:
    """The cleanup ``try`` has to start at discovery, not at the loop.

    Everything between the two - the digest settings read out of the config, the
    duplicate-identity guard - used to sit outside it, so anything raising there left one
    ``cohort_zip_*`` directory per archive on disk. That is how the null-``cohort``
    ``AttributeError`` leaked - see
    ``test_cohort_pseudonym_config.py::test_a_pseudonym_block_that_is_not_a_mapping_falls_back_to_the_defaults``,
    which is the shape of configuration that raised it. The guard is stated generally rather than against the
    one line that raised: ``duplicate_identity`` is made to fail here because it is the
    last thing between discovery and the loop, and no exception from that region may skip
    the cleanup.
    """
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})
    extraction_root = tmp_path / "extracted"

    def _mkdtemp(prefix: str = "") -> str:
        directory = extraction_root / f"{prefix}fixed"
        directory.mkdir(parents=True)
        return str(directory)

    def _explode(_samples: object) -> None:
        raise RuntimeError("raised between discovery and the sample loop")

    monkeypatch.setattr(tempfile, "mkdtemp", _mkdtemp)
    monkeypatch.setattr(cohort_summary, "duplicate_identity", _explode)

    with pytest.raises(RuntimeError, match="between discovery and the sample loop"):
        cohort_summary.aggregate_cohort(
            input_paths=[str(archive)],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort_summary.html",
            config=load_config(None),
            pseudonymize_samples="anon_",
        )

    assert not (extraction_root / "cohort_zip_fixed").exists(), (
        "the extracted archive was left behind: the cleanup try does not cover this region"
    )
