"""One physical sample is one record, however many ways the caller spelled the path to it.

**SPECIFICATION.** ``discover_sample_directories`` de-duplicates, and what it de-duplicates
on has to be the sample's *physical* location rather than the ``Path`` it was reached by.
Keyed on the path as constructed, one sample named through an absolute parent and again
through a relative direct input produced two records that are unequal as values and
identical on disk.

Before identity qualification that was accidentally harmless - both records were named
``sample``, so ``sample_categories()`` grouped them back into one row and the count came
out right for the wrong reason. Qualification made it visible and made it worse: the two
qualify to ``cohort/sample`` and ``sample/sample``, ``duplicate_identity`` sees two
distinct identities and passes them both, and one patient is reported as two samples in
every table, plot and export.

A repeated ZIP is the same defect at the other end. Its sample directory is its own fresh
``mkdtemp`` root, so two extractions of one archive can never compare equal however the key
is canonicalised - ``-i job.zip -i job.zip`` wrote both copies to disk and then aborted the
whole cohort, because the two samples shared a namespace *and* a local value. The canonical
form of each input is therefore remembered before anything is extracted.

``origin`` still reports the path **as the caller wrote it**: canonicalisation is for
identity, and a diagnostic has to quote what was typed.

Split out of ``test_cohort_identity.py`` and ``test_cohort_inputs.py``; the qualification
rule is in the former, ZIP identity and ordering in ``test_cohort_zip_identity.py``.
"""

from __future__ import annotations

import os
import tempfile

import pytest

from tests.support.cohort_samples import sample_on_disk, summary_json, zip_of
from vntyper.cli import load_config
from vntyper.scripts import cohort_inputs, cohort_summary
from vntyper.scripts.cohort_inputs import (
    cleanup_temp_dirs,
    discover_sample_directories,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# directory inputs
# ---------------------------------------------------------------------------


def test_the_same_sample_reached_twice_is_only_processed_once(tmp_path) -> None:
    """The case that always worked: a parent and its child name one sample, spelled alike.

    It moved here from ``test_cohort_inputs.py`` because it is the base case of this
    module's rule, and because reading it beside the spellings that *did not* work is what
    shows the difference is canonicalisation and nothing else.
    """
    sample = sample_on_disk(tmp_path / "cohort" / "sample_one")

    dirs, _ = discover_sample_directories([str(sample), str(tmp_path / "cohort")])

    assert [found.directory for found in dirs] == [sample]


def test_a_sample_reached_by_an_absolute_parent_and_a_relative_child_is_one_sample(tmp_path, monkeypatch) -> None:
    """**Specification**: de-duplication keys the *physical* directory, not its spelling.

    De-duplication was keyed on the `Path` exactly as it was constructed, so one sample
    reached through an absolute parent and again through a relative direct input produced
    two `Path` keys that are unequal as values and identical on disk. Two records came out.

    Before identity qualification that was accidentally harmless: both records were named
    `sample`, so `sample_categories()` grouped them back into one row and the cohort's
    count came out right for the wrong reason. Qualification made it visible and made it
    worse - the two now qualify to `cohort/sample` and `sample/sample`, `duplicate_identity`
    sees two distinct identities and lets both through, and one patient is reported as two
    samples in every table, plot and export.

    `vntyper cohort -i "$PWD/cohort" -i cohort/sample` is one shell expansion away, and the
    web service passes a caller-supplied list straight through.
    """
    sample_on_disk(tmp_path / "cohort" / "sample")
    monkeypatch.chdir(tmp_path)

    dirs, _ = discover_sample_directories([str(tmp_path / "cohort"), "cohort/sample"])

    assert [found.identity for found in dirs] == ["sample"]
    assert [found.directory for found in dirs] == [tmp_path / "cohort" / "sample"]


def test_a_directory_named_twice_with_and_without_a_trailing_slash_is_one_sample(tmp_path) -> None:
    """`Path` drops a trailing separator, so this held already; it is pinned as a guard.

    Canonicalising the key must not be the *only* thing that makes it hold - a regression
    that stopped normalising the trailing slash would otherwise be invisible.
    """
    sample = sample_on_disk(tmp_path / "cohort" / "sample_one")

    dirs, _ = discover_sample_directories([str(tmp_path / "cohort"), str(tmp_path / "cohort") + os.sep])

    assert [found.directory for found in dirs] == [sample]


def test_a_directory_reached_through_a_symlink_is_not_a_second_sample(tmp_path) -> None:
    """A link and its target are one physical directory, so they are one sample.

    `Path.resolve()` follows the link, which is exactly why the key is resolved rather
    than compared as written: without it the cohort reports `cohort/sample_one` and
    `link/sample_one` - one patient twice.
    """
    sample = sample_on_disk(tmp_path / "cohort" / "sample_one")
    link = tmp_path / "link"
    link.symlink_to(tmp_path / "cohort", target_is_directory=True)

    dirs, _ = discover_sample_directories([str(tmp_path / "cohort"), str(link)])

    assert [found.identity for found in dirs] == ["sample_one"]
    assert [found.directory for found in dirs] == [sample]


def test_the_same_sample_reached_twice_is_not_a_collision(tmp_path) -> None:
    """De-duplication is not a collision: one directory named by two input paths is one sample."""
    sample = sample_on_disk(tmp_path / "cohort" / "sample_one")
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


# ---------------------------------------------------------------------------
# archive inputs, and the extractions they own
# ---------------------------------------------------------------------------


def test_the_same_archive_named_twice_is_extracted_once(tmp_path) -> None:
    """**Specification**: a repeated input is de-duplicated *before* it is extracted.

    De-duplication happened on the sample directory alone, and a zip's sample directory is
    its own fresh `tempfile.mkdtemp` root - so naming one archive twice extracted it twice
    into two directories that can never compare equal. The two samples then carried the
    same namespace *and* the same local value, which is precisely the shape
    `duplicate_identity` refuses: `vntyper cohort -i job.zip -i job.zip` aborted the whole
    cohort over one archive named twice, having already written both copies to disk.
    """
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json()})

    dirs, temp_dirs = discover_sample_directories([str(archive), str(archive)])

    try:
        assert [found.identity for found in dirs] == ["job7"]
        assert len(temp_dirs) == 1
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_the_same_archive_reached_by_two_spellings_is_extracted_once(tmp_path, monkeypatch) -> None:
    """The canonical form is what is compared, so an absolute and a relative spelling agree."""
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json()})
    monkeypatch.chdir(tmp_path)

    dirs, temp_dirs = discover_sample_directories([str(archive), "job7.zip"])

    try:
        assert [found.identity for found in dirs] == ["job7"]
        assert [found.origin for found in dirs] == [str(archive)]
        assert len(temp_dirs) == 1
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_failure_after_extraction_removes_what_discovery_had_already_extracted(tmp_path, monkeypatch) -> None:
    """Discovery owns its extractions until it returns them, so it must clean up its own.

    The caller's cleanup `try` opens on the temporary directories `discover_sample_directories`
    *returns*, so it cannot cover anything that raises before the return - and the work
    after the first extraction is not trivial: further inputs are extracted and searched,
    the list is sorted, and identities are qualified over the whole cohort. Anything
    raising in there left one `cohort_zip_*` directory per archive extracted so far, with
    no handle on them anywhere.

    `qualify_colliding_identities` is made to fail because it is the last thing discovery
    does before returning; the guard is about the region, not about that one call.
    """
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json()})
    extraction_root = tmp_path / "extracted"

    def _mkdtemp(prefix: str = "") -> str:
        directory = extraction_root / f"{prefix}fixed"
        directory.mkdir(parents=True)
        return str(directory)

    def _explode(_samples: object) -> None:
        raise RuntimeError("raised after the archive was extracted")

    monkeypatch.setattr(tempfile, "mkdtemp", _mkdtemp)
    monkeypatch.setattr(cohort_inputs, "qualify_colliding_identities", _explode)

    with pytest.raises(RuntimeError, match="after the archive was extracted"):
        discover_sample_directories([str(archive)])

    assert list(extraction_root.glob("cohort_zip_*")) == [], (
        "the extracted archive was left behind: discovery does not clean up its own extractions"
    )
