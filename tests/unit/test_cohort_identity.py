"""A cohort pseudonym must identify exactly one sample, and a sample must keep its name (#206, #205).

SPECIFICATION, in three parts.

**The pseudonym.** It was ``<prefix>`` plus the first **five** hex characters of the MD5
of the sample name - 20 bits. ``sample_mapping`` is keyed on the pseudonym, so a collision
silently overwrote the earlier original name: two patients' rows became indistinguishable
across every export, ``sample_categories()`` counted them as one sample, and
``pseudonymization_table.tsv`` mapped the shared pseudonym to whichever original was
written last. Birthday probability of at least one collision is ``1 - exp(-n(n-1)/2**21)``:
~37.9% at 1,000 samples. The verified first collision in ``sample_0..sample_19999`` was
``sample_42`` and ``sample_919``, both MD5-prefixing to ``168eb``; that exact pair is the
probe below.

**The map.** Widening the digest does not make the map injective, because the *inputs* can
already be equal: two discovered directories ``a/sample`` and ``b/sample`` share a
basename, so they share an identity before anything is hashed. No digest width fixes
identical inputs, so ``aggregate_cohort`` refuses to run rather than merging two patients.

**The identity and the order of a ZIP sample.** A ZIP is extracted into
``tempfile.mkdtemp(prefix="cohort_zip_")``. When its ``pipeline_summary.json`` sits at the
archive root - the layout the web worker produces, so the normal path for web cohorts -
that temporary directory *was* the sample directory, and the reported ``Sample`` was
therefore a random ``cohort_zip_*`` string. The same random component was in the sort key,
so two runs over the same archives ordered their rows differently. Neither depends on the
inputs any more: the identity comes from the input file the run itself recorded, and the
sort key is ``(index of the input on the command line, path relative to that input's
root)``.

At 48 bits the collision guard is a tripwire rather than an operational risk:
``1 - exp(-n(n-1)/2**49)`` is 1.78e-9 at 1,000 samples and 1.78e-7 at 10,000.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import pytest

from vntyper.cli import load_config
from vntyper.scripts import cohort_summary
from vntyper.scripts.cohort_inputs import (
    DEFAULT_PSEUDONYM_ALGORITHM,
    DEFAULT_PSEUDONYM_LENGTH,
    cleanup_temp_dirs,
    discover_sample_directories,
    pseudonymized_sample_name,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# C1 - the digest
# ---------------------------------------------------------------------------


def test_the_verified_md5_collision_no_longer_collides() -> None:
    """`sample_42` and `sample_919` were the first colliding pair under the old scheme."""
    assert pseudonymized_sample_name("anon_", "sample_42") != pseudonymized_sample_name("anon_", "sample_919")


def test_pseudonyms_are_injective_over_twenty_thousand_names() -> None:
    """The exact probe that found the bug, at the configured width."""
    names = [f"sample_{index}" for index in range(20000)]
    pseudonyms = {pseudonymized_sample_name("anon_", name) for name in names}

    assert len(pseudonyms) == len(names)


def test_the_pseudonym_is_stable_across_calls() -> None:
    """The pseudonymization table written beside the report has to stay meaningful."""
    assert pseudonymized_sample_name("anon_", "s1") == pseudonymized_sample_name("anon_", "s1")


def test_the_default_digest_is_twelve_characters_of_sha256() -> None:
    """The literal is recorded rather than recomputed with the code under test.

    ``sha256`` has no per-process salt, so a recorded value is the whole cross-process
    stability guarantee - the same reason the MD5 literal was recorded before it.
    """
    assert DEFAULT_PSEUDONYM_ALGORITHM == "sha256"
    assert DEFAULT_PSEUDONYM_LENGTH == 12
    assert pseudonymized_sample_name("anon_", "s1") == "anon_e8bc163c82ee"
    assert pseudonymized_sample_name("anon_", "sample_one") == "anon_c788e939395d"


def test_the_digest_width_is_honoured() -> None:
    assert pseudonymized_sample_name("", "s1", length=5) == "e8bc1"
    assert len(pseudonymized_sample_name("", "s1", length=5)) == 5
    assert len(pseudonymized_sample_name("", "s1", length=64)) == 64


def test_a_non_string_prefix_is_interpolated_rather_than_concatenated() -> None:
    """Preserved behaviour: ``--pseudonymize-samples`` may carry a non-string."""
    assert pseudonymized_sample_name(True, "s1") == "Truee8bc163c82ee"


def test_an_unavailable_algorithm_is_refused_by_name() -> None:
    """A silent fallback would change every pseudonym in the report without saying so."""
    with pytest.raises(ValueError, match="not-a-hash"):
        pseudonymized_sample_name("anon_", "s1", algorithm="not-a-hash")


@pytest.mark.parametrize("length", [0, -1, 2.5, "12", None])
def test_a_digest_length_that_is_not_a_positive_integer_is_refused(length: object) -> None:
    """``digest_characters`` comes out of a JSON file, so it can be anything at all."""
    with pytest.raises(ValueError, match="positive integer"):
        pseudonymized_sample_name("anon_", "s1", length=length)  # type: ignore[arg-type]


def test_a_digest_length_wider_than_the_digest_is_refused() -> None:
    """`sha256` produces 64 hex characters; asking for 65 would silently give 64."""
    with pytest.raises(ValueError, match="65"):
        pseudonymized_sample_name("anon_", "s1", length=65)


def test_a_boolean_length_is_refused_rather_than_read_as_one() -> None:
    """``True`` is an ``int`` in Python, and a one-character pseudonym is not a pseudonym."""
    with pytest.raises(ValueError, match="positive integer"):
        pseudonymized_sample_name("anon_", "s1", length=True)  # type: ignore[arg-type]


def test_an_algorithm_that_needs_a_digest_length_is_refused_by_name() -> None:
    """`shake_128` is in ``algorithms_available`` but its ``hexdigest()`` takes an argument.

    Left unguarded that arrives as a ``TypeError`` from inside ``hashlib``, which is not
    the repository's error convention and does not name the configuration key at fault.
    """
    assert "shake_128" in hashlib.algorithms_available
    with pytest.raises(ValueError, match="shake_128"):
        pseudonymized_sample_name("anon_", "s1", algorithm="shake_128")


def test_the_shipped_config_declares_the_pseudonym_settings() -> None:
    """Config-driven, never hardcoded: the values live in config.json."""
    config = json.loads((REPO_ROOT / "vntyper" / "config.json").read_text(encoding="utf-8"))

    assert config["cohort"]["pseudonym"] == {"algorithm": "sha256", "digest_characters": 12}


def test_the_cohort_reads_the_digest_settings_out_of_the_configuration(tmp_path) -> None:
    """`--config-path` is how an operator changes the width, so the read has to work.

    The narrowed digest is asserted through the pseudonymization table rather than
    through the function, because the point is that ``aggregate_cohort`` threads the
    configured pair down to it.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    config = load_config(None)
    config["cohort"]["pseudonym"] = {"algorithm": "sha1", "digest_characters": 4}

    cohort_summary.aggregate_cohort(
        input_paths=[str(_sample_on_disk(tmp_path / "cohort" / "sample_one"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=config,
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    expected = "anon_" + hashlib.sha1(b"sample_one").hexdigest()[:4]  # noqa: S324 - not a security use
    assert f"{expected}\tsample_one" in table


def test_a_configuration_without_the_pseudonym_block_falls_back_to_the_defaults(tmp_path) -> None:
    """AGENTS.md trap 2: ``--config-path`` replaces the whole config rather than merging it.

    A config that never heard of ``cohort.pseudonym`` - which is every config written
    before this milestone - must produce the default pseudonym rather than a ``KeyError``.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_sample_on_disk(tmp_path / "cohort" / "sample_one"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config={"paths": {"template_dir": "vntyper/templates"}},
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert "anon_c788e939395d\tsample_one" in table


# ---------------------------------------------------------------------------
# C2 - the map is injective, or the run stops
# ---------------------------------------------------------------------------


def _summary(input_files: dict[str, str] | None = None) -> str:
    """Render a minimal ``pipeline_summary.json``.

    Args:
        input_files: The ``input_files`` mapping the run recorded.

    Returns:
        str: The serialised summary.
    """
    return json.dumps({"version": "2.0.6", "input_files": input_files or {}, "steps": []})


def _sample_on_disk(directory: Path, input_files: dict[str, str] | None = None) -> Path:
    """Write one sample directory holding a ``pipeline_summary.json``.

    Args:
        directory: The sample directory to create.
        input_files: The ``input_files`` mapping the run recorded.

    Returns:
        Path: The directory that was created.
    """
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "pipeline_summary.json").write_text(_summary(input_files), encoding="utf-8")
    return directory


def _zip_of(archive: Path, members: dict[str, str]) -> Path:
    """Build a zip archive from a member -> contents mapping.

    Args:
        archive: Where to write the archive.
        members: Archive-relative path -> file contents.

    Returns:
        Path: The archive.
    """
    archive.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, "w") as handle:
        for member, contents in members.items():
            handle.writestr(member, contents)
    return archive


def test_two_directories_sharing_a_basename_abort_the_cohort(tmp_path) -> None:
    """The collision that no digest width can fix, because it is upstream of the digest.

    ``a/sample`` and ``b/sample`` are two patients whose identities are already equal,
    so a guard comparing the two *originals* sees them as the same sample and passes.
    ``sample_categories()`` would then count one result where there are two.
    """
    first = _sample_on_disk(tmp_path / "a" / "sample")
    second = _sample_on_disk(tmp_path / "b" / "sample")

    with pytest.raises(ValueError) as excinfo:
        cohort_summary.aggregate_cohort(
            input_paths=[str(first), str(second)],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort_summary.html",
            config=load_config(None),
        )

    message = str(excinfo.value)
    assert str(first) in message
    assert str(second) in message
    assert "sample" in message


def test_a_shared_basename_aborts_even_without_pseudonymisation(tmp_path) -> None:
    """The merge happens in ``sample_categories``, which never sees the pseudonym."""
    _sample_on_disk(tmp_path / "a" / "sample")
    _sample_on_disk(tmp_path / "b" / "sample")

    with pytest.raises(ValueError, match="Duplicate cohort sample identity"):
        cohort_summary.aggregate_cohort(
            input_paths=[str(tmp_path / "a"), str(tmp_path / "b")],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort_summary.html",
            config=load_config(None),
            pseudonymize_samples=False,
        )


def test_two_web_jobs_that_uploaded_the_same_filename_abort_and_name_their_archives(tmp_path) -> None:
    """The consequence of #205's identity change that an operator will actually meet.

    A web job's archive is a root-level zip, so its identity is now the stem of the file
    that was uploaded. Two jobs that both uploaded ``sample.bam`` are therefore two
    patients with one identity - where before the fix they got two different random
    ``cohort_zip_*`` names and were merely unreadable rather than merged. Aborting is the
    milestone's decision, so what this test defends is that the message is *actionable*:
    the sample directories are extraction directories and say nothing, so each is named
    together with the archive it came out of.
    """
    first = _zip_of(tmp_path / "job_a.zip", {"pipeline_summary.json": _summary({"bam": "sample.bam"})})
    second = _zip_of(tmp_path / "job_b.zip", {"pipeline_summary.json": _summary({"bam": "sample.bam"})})

    with pytest.raises(ValueError) as excinfo:
        cohort_summary.aggregate_cohort(
            input_paths=[str(first), str(second)],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort_summary.html",
            config=load_config(None),
        )

    message = str(excinfo.value)
    assert str(first) in message
    assert str(second) in message
    assert "'sample'" in message


def test_an_aborted_cohort_removes_the_archives_it_extracted(tmp_path) -> None:
    """The guard fires after extraction, so refusing to run must not leave the zips behind."""
    first = _zip_of(tmp_path / "job_a.zip", {"pipeline_summary.json": _summary({"bam": "sample.bam"})})
    second = _zip_of(tmp_path / "job_b.zip", {"pipeline_summary.json": _summary({"bam": "sample.bam"})})
    before = set(Path(tempfile.gettempdir()).glob("cohort_zip_*"))

    with pytest.raises(ValueError):
        cohort_summary.aggregate_cohort(
            input_paths=[str(first), str(second)],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort_summary.html",
            config=load_config(None),
        )

    assert set(Path(tempfile.gettempdir()).glob("cohort_zip_*")) == before


def test_two_originals_mapping_to_one_pseudonym_abort_the_cohort(tmp_path, monkeypatch) -> None:
    """A cohort that silently merges two patients is worse than one that refuses to run.

    At 48 bits this is unreachable in practice, so it is forced here by stubbing the
    digest. That is the point: no digest width makes a silent merge acceptable, so the
    guard does not depend on the width.
    """
    monkeypatch.setattr(
        cohort_summary,
        "pseudonymized_sample_name",
        lambda prefix, name, **kwargs: f"{prefix}same",
    )
    _sample_on_disk(tmp_path / "cohort" / "sample_a")
    _sample_on_disk(tmp_path / "cohort" / "sample_b")

    with pytest.raises(ValueError) as excinfo:
        cohort_summary.aggregate_cohort(
            input_paths=[str(tmp_path / "cohort")],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort_summary.html",
            config=load_config(None),
            pseudonymize_samples="anon_",
        )

    message = str(excinfo.value)
    assert "sample_a" in message
    assert "sample_b" in message
    assert "anon_same" in message


def test_the_same_sample_reached_twice_is_not_a_collision(tmp_path) -> None:
    """De-duplication is not a collision: one directory named by two input paths is one sample."""
    sample = _sample_on_disk(tmp_path / "cohort" / "sample_one")
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(sample), str(tmp_path / "cohort")],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert table.count("sample_one") == 1


def test_two_distinct_samples_are_reported_under_two_pseudonyms(tmp_path) -> None:
    """The control: the guard must not fire on the ordinary two-sample cohort."""
    _sample_on_disk(tmp_path / "cohort" / "sample_one")
    _sample_on_disk(tmp_path / "cohort" / "sample_two")
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(tmp_path / "cohort")],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert "anon_c788e939395d\tsample_one" in table
    assert "anon_" + hashlib.sha256(b"sample_two").hexdigest()[:12] + "\tsample_two" in table


# ---------------------------------------------------------------------------
# C3 - a ZIP sample's identity and its place in the order
# ---------------------------------------------------------------------------


def test_a_root_level_zip_sample_is_identified_by_the_run_s_own_input(tmp_path) -> None:
    """#205: this layout is what the web worker produces, so it is the normal path.

    ``processed_dirs.add(temp_path)`` added the ``tempfile.mkdtemp`` root itself, and
    ``cohort_summary`` took ``Path(sample_dir).name`` as the sample - so the reported
    ``Sample`` was a ``cohort_zip_*`` string that differed on every run, and any pseudonym
    derived from it differed too.
    """
    archive = _zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": _summary({"bam": "patient1.bam"})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["patient1"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_root_level_zip_without_input_files_falls_back_to_the_archive_stem(tmp_path) -> None:
    """Older summaries recorded no ``input_files``; the archive's own name is the next best thing."""
    archive = _zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": _summary({})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["job7"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_root_level_zip_whose_summary_is_unreadable_falls_back_to_the_archive_stem(tmp_path) -> None:
    """One bad sample must not abort the cohort, so the identity helper never raises."""
    archive = _zip_of(tmp_path / "job8.zip", {"pipeline_summary.json": "{not json"})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["job8"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_cram_run_is_identified_by_its_cram(tmp_path) -> None:
    """``bam`` first, then ``cram``, then ``fastq1`` - the order ``pipeline.py`` writes them."""
    archive = _zip_of(tmp_path / "job9.zip", {"pipeline_summary.json": _summary({"cram": "patient2.cram"})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["patient2"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_fastq_run_is_identified_by_its_first_mate(tmp_path) -> None:
    """``patient3_R1.fastq.gz`` carries two suffixes, and one ``.stem`` strips only one."""
    archive = _zip_of(
        tmp_path / "job10.zip",
        {
            "pipeline_summary.json": _summary(
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
    archive = _zip_of(
        tmp_path / "job11.zip",
        {"pipeline_summary.json": _summary({"fastq1": "reads_R1.fastq.gz", "bam": "patient4.bam"})},
    )

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["patient4"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_sample_in_a_subdirectory_keeps_its_directory_name(tmp_path) -> None:
    """Only the root-level case changes: elsewhere the directory name is meaningful."""
    archive = _zip_of(
        tmp_path / "cohort.zip",
        {
            "sample_one/pipeline_summary.json": _summary({"bam": "patient1.bam"}),
            "sample_two/pipeline_summary.json": _summary({"bam": "patient2.bam"}),
        },
    )

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert [sample.identity for sample in samples] == ["sample_one", "sample_two"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_the_sort_key_contains_no_temporary_path_component(tmp_path) -> None:
    """The mechanism, stated directly: nothing in the key comes from ``mkdtemp``.

    Both halves are asserted as values. The outer half is the archive path the caller
    wrote and the inner half is empty for a root-level sample, so the whole key is
    reconstructible from the command line alone - which is what makes it reproducible.
    """
    first = _zip_of(tmp_path / "aaa_cohort.zip", {"pipeline_summary.json": _summary({"bam": "patient1.bam"})})
    second = _zip_of(tmp_path / "zzz_cohort.zip", {"pipeline_summary.json": _summary({"bam": "patient2.bam"})})

    samples, temp_dirs = discover_sample_directories([str(first), str(second)])

    try:
        assert [sample.sort_key for sample in samples] == [(first.parts, ()), (second.parts, ())]
        assert [sample.origin for sample in samples] == [str(first), str(second)]
        assert all("cohort_zip_" not in part for sample in samples for half in sample.sort_key for part in half)
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
    first = _zip_of(tmp_path / "aaa_cohort.zip", {"pipeline_summary.json": _summary({"bam": "patient1.bam"})})
    second = _zip_of(tmp_path / "zzz_cohort.zip", {"pipeline_summary.json": _summary({"bam": "patient2.bam"})})
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
    first = _zip_of(tmp_path / "aaa_cohort.zip", {"pipeline_summary.json": _summary({"bam": "patient1.bam"})})
    second = _zip_of(tmp_path / "zzz_cohort.zip", {"pipeline_summary.json": _summary({"bam": "patient2.bam"})})

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
        _zip_of(tmp_path / f"job_{index}.zip", {"pipeline_summary.json": _summary({"bam": f"patient{index}.bam"})})
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
    directory = _sample_on_disk(tmp_path / "cohort" / "sample_one")
    archive = _zip_of(tmp_path / "job.zip", {"pipeline_summary.json": _summary({"bam": "patient1.bam"})})

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
    archive = _zip_of(
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
    archive = _zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": _summary({"bam": "patient1.bam"})})

    samples, temp_dirs = discover_sample_directories([str(archive)])

    try:
        assert len(temp_dirs) == 1
        assert samples[0].directory == Path(temp_dirs[0])
        assert Path(temp_dirs[0]).name.startswith("cohort_zip_")
    finally:
        cleanup_temp_dirs(temp_dirs)
    assert not Path(temp_dirs[0]).exists()


def test_a_directory_sample_sorts_by_its_path_below_the_input_root(tmp_path) -> None:
    """Within one input the key is the relative path's *parts*, so the separator never sorts.

    ``cohort/sample_one`` comes before ``cohort-extra/sample_one`` because
    ``"cohort" < "cohort-extra"``, where the raw strings compare the other way round
    (``"-"`` is 0x2d and ``"/"`` is 0x2f). It is the order ``ls`` would show, which is the
    point of choosing it, and it survives the move from whole paths to relative ones.
    """
    plain = _sample_on_disk(tmp_path / "cohort" / "sample_one")
    suffixed = _sample_on_disk(tmp_path / "cohort-extra" / "sample_one")
    assert str(suffixed) < str(plain)

    samples, _ = discover_sample_directories([str(tmp_path)])

    assert [sample.directory for sample in samples] == [plain, suffixed]
    assert [sample.sort_key for sample in samples] == [
        (tmp_path.parts, ("cohort", "sample_one")),
        (tmp_path.parts, ("cohort-extra", "sample_one")),
    ]


def test_the_cohort_reports_the_zip_identity_rather_than_the_extraction_directory(tmp_path) -> None:
    """`aggregate_cohort` consumes the identity rather than re-deriving it from the path."""
    archive = _zip_of(
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
    sample = _sample_on_disk(tmp_path / "cohort" / "sample_one")

    samples, _ = discover_sample_directories([str(sample)])

    (found,) = samples
    assert found.directory == sample
    assert found.identity == "sample_one"
    assert found.origin == str(sample)
    assert found.sort_key == (sample.parts, ())
    with pytest.raises(dataclasses.FrozenInstanceError):
        found.identity = "renamed"  # type: ignore[misc]
