"""A cohort pseudonym must identify exactly one sample, and a sample must keep its name (#206, #205).

SPECIFICATION, in two parts.

**The map.** Widening the digest does not make the map injective, because the *inputs* can
already be equal: two discovered directories called ``sample`` share a name, so they share
a local identity before anything is hashed. No digest width fixes identical inputs.

**The namespace, and what is left to refuse.** The first answer to that was to abort, and
it was too broad: two web jobs that each uploaded ``sample.bam`` are a commonplace input
on the normal web path, and refusing a legitimate cohort is a regression on the behaviour
before #205, where both jobs ran. A shared *value* is not a shared *identifier*. HL7
FHIR's ``Identifier`` makes uniqueness a property of ``(system, value)`` together - the
value is only "unique within the context of the system" - so a name means nothing until
you say which namespace issued it. Each sample's namespace is its input's own name (the
archive stem, or the directory name; the *name*, so that ``-i /data/run1`` and ``-i ./run1``
agree and pseudonyms survive relocation), and ``qualify_colliding_identities`` rewrites a
shared local value to ``namespace/value`` **for the colliding samples only** - a unique
name, and the pseudonym computed from it, never moves.

What ``aggregate_cohort`` still refuses is what no namespace separates: two samples equal
in namespace *and* local value (two archives with one stem, two directory inputs with one
name), and a digest collision between two identities that are already distinct. The
asymmetry between those two is deliberate - for a digest collision the names are perfectly
good and only the width is wrong, so qualifying would be gratuitous; for a name collision
the names are the ambiguity, so qualifying is the only way to proceed at all.

At 48 bits the collision guard is a tripwire rather than an operational risk:
``1 - exp(-n(n-1)/2**49)`` is 1.78e-9 at 1,000 samples and 1.78e-7 at 10,000.

This module was 1,210 lines and is now the qualification rule alone. The digest and the
configuration it is read from moved to ``test_cohort_pseudonym_config.py``, a ZIP sample's
identity and its place in the order to ``test_cohort_zip_identity.py``, and the rule that
one physical sample is one record to ``test_cohort_deduplication.py``.
"""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

import pytest

from tests.support.cohort_samples import pseudonym_of, sample_on_disk, summary_json, zip_of
from vntyper.cli import load_config
from vntyper.scripts import cohort_summary
from vntyper.scripts.cohort_inputs import (
    DiscoveredSample,
    cleanup_temp_dirs,
    discover_sample_directories,
    identity_namespace,
    qualify_colliding_identities,
)

pytestmark = pytest.mark.unit


def test_two_inputs_that_share_a_namespace_and_a_local_name_abort_the_cohort(tmp_path) -> None:
    """The one collision qualification cannot reduce, because the namespaces collide too.

    ``a/sample`` and ``b/sample`` named **as the two inputs** are two patients whose
    identities are already equal, so a guard comparing the two *originals* sees them as the
    same sample and passes; ``sample_categories()`` would then count one result where there
    are two. Qualification does not help here and is not applied: each sample's namespace
    is *its input's own name*, and both inputs are named ``sample``, so both would qualify
    to ``sample/sample``. The namespace is the input's name rather than its whole path on
    purpose - that is what makes a pseudonym survive relocation - and this irreducible case
    is the price of it. Aborting names both input paths so the operator can tell them apart.
    """
    first = sample_on_disk(tmp_path / "a" / "sample")
    second = sample_on_disk(tmp_path / "b" / "sample")

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


def test_two_directory_inputs_holding_a_shared_basename_are_qualified_rather_than_refused(tmp_path) -> None:
    """**This test used to claim the opposite, and the claim was wrong.**

    It was ``test_a_shared_basename_aborts_even_without_pseudonymisation``: two directory
    inputs ``a`` and ``b`` each holding a ``sample`` directory aborted the whole cohort,
    on the reasoning that the two identities were equal and ``sample_categories()`` would
    merge two patients into one row.

    What that reasoning discarded is the *system* the value came from. HL7 FHIR's
    ``Identifier`` makes uniqueness a property of ``(system, value)`` together: the value is
    "unique within the context of the system", never on its own. ``sample`` is a value, and
    ``a`` and ``b`` are two systems, so these are two identifiers and not a collision at
    all. The identity is qualified to ``a/sample`` and ``b/sample`` and both samples are
    reported - which is what the pre-#205 code did by accident, under names that changed
    every run.
    """
    sample_on_disk(tmp_path / "a" / "sample")
    sample_on_disk(tmp_path / "b" / "sample")
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(tmp_path / "a"), str(tmp_path / "b")],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert pseudonym_of("a/sample") + "\ta/sample" in table
    assert pseudonym_of("b/sample") + "\tb/sample" in table


def test_two_web_jobs_that_uploaded_the_same_filename_are_reported_under_both_archives(tmp_path) -> None:
    """**This test used to claim the opposite, and the claim was wrong.**

    It was ``test_two_web_jobs_that_uploaded_the_same_filename_abort_and_name_their_archives``
    and it asserted that two web jobs both uploading ``sample.bam`` abort the entire cohort,
    with an actionable message naming each archive. The message was the right shape; the
    abort was not. ``sample.bam`` is a commonplace upload name and a root-level ZIP is
    exactly what the web worker produces, so that abort is on the normal path, and before
    #205 both jobs ran - under meaningless ``cohort_zip_*`` names that differed every run.
    Refusing a legitimate cohort is a regression on both.

    It inverted because ``sample`` is a **value**, not an identifier. HL7 FHIR's
    ``Identifier`` fixes uniqueness on ``(system, value)`` together - the value is only
    "unique within the context of the system" - and these two values sit in two different
    systems, ``job_a`` and ``job_b``. So they are two identifiers, they are qualified to
    ``job_a/sample`` and ``job_b/sample``, and the cohort completes with both patients in it.
    """
    first = zip_of(tmp_path / "job_a.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
    second = zip_of(tmp_path / "job_b.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(first), str(second)],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert pseudonym_of("job_a/sample") + "\tjob_a/sample" in table
    assert pseudonym_of("job_b/sample") + "\tjob_b/sample" in table
    assert "cohort_zip_" not in table


def test_three_archives_sharing_one_upload_name_are_all_qualified(tmp_path) -> None:
    """Qualification is over the whole discovered list, so a three-way collision resolves too.

    A pairwise fix - rename the second of each colliding pair - would leave the third
    sample sharing a name with one of the other two.
    """
    archives = [
        zip_of(tmp_path / f"job_{letter}.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
        for letter in ("a", "b", "c")
    ]

    samples, temp_dirs = discover_sample_directories([str(archive) for archive in archives])

    try:
        assert [sample.identity for sample in samples] == ["job_a/sample", "job_b/sample", "job_c/sample"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_only_the_samples_that_collide_are_qualified(tmp_path) -> None:
    """A sample whose local name is already unique keeps it, so its pseudonym never moves.

    This is the property the whole design turns on: qualifying every sample would rename
    every sample in every cohort that happens to contain one collision, and every pseudonym
    in ``pseudonymization_table.tsv`` with it.
    """
    first = zip_of(tmp_path / "job_a.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
    second = zip_of(tmp_path / "job_b.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
    third = zip_of(tmp_path / "job_c.zip", {"pipeline_summary.json": summary_json({"bam": "patient9.bam"})})

    samples, temp_dirs = discover_sample_directories([str(first), str(second), str(third)])

    try:
        assert [sample.identity for sample in samples] == ["job_a/sample", "job_b/sample", "patient9"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_two_archives_with_the_same_stem_still_abort_the_cohort(tmp_path) -> None:
    """Same namespace, same local value: qualification produces one identity twice.

    Two archives called ``job.zip`` in two directories share a namespace, because the
    namespace is the archive's *name*. Both samples qualify to ``job/sample``, so the
    ambiguity is irreducible and the cohort stops rather than merging two patients. The
    message has to name the archives, since the sample directories are extraction
    directories and say nothing.
    """
    first = zip_of(tmp_path / "run1" / "job.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
    second = zip_of(tmp_path / "run2" / "job.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})

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
    assert "'job/sample'" in message


def test_an_aborted_cohort_removes_the_archives_it_extracted(tmp_path) -> None:
    """The guard fires after extraction, so refusing to run must not leave the zips behind.

    The fixture is two archives sharing a stem rather than two archives recording one
    upload name: the latter is qualified and completes now, so it no longer reaches the
    abort this is about.
    """
    first = zip_of(tmp_path / "run1" / "job.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
    second = zip_of(tmp_path / "run2" / "job.zip", {"pipeline_summary.json": summary_json({"bam": "sample.bam"})})
    before = set(Path(tempfile.gettempdir()).glob("cohort_zip_*"))

    with pytest.raises(ValueError):
        cohort_summary.aggregate_cohort(
            input_paths=[str(first), str(second)],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort_summary.html",
            config=load_config(None),
        )

    assert set(Path(tempfile.gettempdir()).glob("cohort_zip_*")) == before


def test_a_cohort_with_no_collision_keeps_the_identities_and_pseudonyms_it_had(tmp_path) -> None:
    """The regression guard for qualifying too eagerly, stated against recorded literals.

    Every identity here is already unique, so nothing may be qualified: the identities are
    the bare local names and the pseudonyms are the same twelve SHA-256 characters they
    were before qualification existed. ``anon_c788e939395d`` is the literal
    ``test_cohort_summary_oracle.py`` pins as well - if it moves, both files fail and the
    implementation has qualified a sample it had no business touching.
    """
    sample_on_disk(tmp_path / "cohort" / "sample_one")
    sample_on_disk(tmp_path / "cohort" / "sample_two")
    archive = zip_of(tmp_path / "job7.zip", {"pipeline_summary.json": summary_json({"bam": "patient1.bam"})})
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    samples, temp_dirs = discover_sample_directories([str(tmp_path / "cohort"), str(archive)])
    try:
        assert [sample.identity for sample in samples] == ["sample_one", "sample_two", "patient1"]
    finally:
        cleanup_temp_dirs(temp_dirs)

    cohort_summary.aggregate_cohort(
        input_paths=[str(tmp_path / "cohort"), str(archive)],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert "anon_c788e939395d\tsample_one" in table
    assert "anon_" + hashlib.sha256(b"sample_two").hexdigest()[:12] + "\tsample_two" in table
    assert "anon_" + hashlib.sha256(b"patient1").hexdigest()[:12] + "\tpatient1" in table


def _discovered(identity: str, origin: str, directory: str) -> DiscoveredSample:
    """Build one ``DiscoveredSample`` without touching the filesystem.

    Args:
        identity: The local identity discovery gave the sample.
        origin: The input path as the caller wrote it.
        directory: Where the sample's ``pipeline_summary.json`` is.

    Returns:
        DiscoveredSample: The record, with a sort key that mirrors discovery's and the
        local identity discovery would have given it - the two are equal until
        ``qualify_colliding_identities`` rewrites the identity.
    """
    return DiscoveredSample(
        directory=Path(directory),
        identity=identity,
        local_identity=identity,
        origin=origin,
        sort_key=Path(origin).parts,
    )


@pytest.mark.parametrize(
    ("origin", "expected"),
    [
        ("/data/run1", "run1"),
        ("./run1", "run1"),
        ("run1", "run1"),
        ("/data/run1/", "run1"),
        ("/data/job_a.zip", "job_a"),
        ("./job_a.zip", "job_a"),
        ("/data/job_a.ZIP", "job_a"),
        ("/data/release.v2", "release.v2"),
        (".", "."),
        ("/", "/"),
    ],
    ids=[
        "absolute-dir",
        "relative-dir",
        "bare-dir",
        "trailing-slash",
        "zip",
        "relative-zip",
        "upper-zip",
        "dotted-dir",
        "cwd",
        "root",
    ],
)
def test_the_namespace_is_the_input_s_own_name_rather_than_its_path(origin: str, expected: str) -> None:
    """**Specification**: the namespace is the input's *name*, and a ZIP loses its suffix.

    Using the name rather than the whole path is what lets ``-i /data/run1`` and
    ``-i ./run1`` name the same system, so a qualified identity - and therefore the
    pseudonym derived from it - survives the cohort being relocated or re-run from another
    working directory. A directory keeps its whole name, dots included: only a ``.zip``
    suffix is stripped, because that suffix is the container and not part of the name the
    run was known by.

    ``vntyper cohort -i .`` names no component of its own -- ``Path(".").name`` is the
    empty string -- and an empty namespace would qualify to a leading ``/``. The input
    string stands in for it instead, so a qualified identity is never nameless. It buys
    nothing when both colliding samples came through that same input, which is the usual
    shape of it: they share the namespace too and the abort at the top of this module is
    the answer.
    """
    assert identity_namespace(origin) == expected


def test_qualification_leaves_a_unique_identity_untouched() -> None:
    """A sample whose local value is already unique is returned unchanged, object and all."""
    samples = [
        _discovered("patient1", "/data/job_a.zip", "/tmp/cohort_zip_aaa"),
        _discovered("patient2", "/data/job_b.zip", "/tmp/cohort_zip_bbb"),
    ]

    assert qualify_colliding_identities(samples) == samples


def test_qualification_replaces_only_the_colliding_members_and_keeps_the_order() -> None:
    """The pure rule, stated over a hand-built list: order in, order out; locals preserved.

    ``origin`` and ``directory`` must survive untouched as well - they are what the abort
    message and the summary reader use, and a qualified sample is still the same sample.
    """
    samples = [
        _discovered("sample", "/data/job_a.zip", "/tmp/cohort_zip_aaa"),
        _discovered("keeper", "/data/job_b.zip", "/tmp/cohort_zip_bbb"),
        _discovered("sample", "/elsewhere/job_c.zip", "/tmp/cohort_zip_ccc"),
    ]

    qualified = qualify_colliding_identities(samples)

    assert [sample.identity for sample in qualified] == ["job_a/sample", "keeper", "job_c/sample"]
    assert [sample.origin for sample in qualified] == [sample.origin for sample in samples]
    assert [sample.directory for sample in qualified] == [sample.directory for sample in samples]
    assert [sample.sort_key for sample in qualified] == [sample.sort_key for sample in samples]


def test_qualification_does_not_mutate_the_list_it_was_given() -> None:
    """``DiscoveredSample`` is frozen, so the rule has to build new records rather than edit."""
    samples = [
        _discovered("sample", "/data/job_a.zip", "/tmp/cohort_zip_aaa"),
        _discovered("sample", "/data/job_b.zip", "/tmp/cohort_zip_bbb"),
    ]

    qualify_colliding_identities(samples)

    assert [sample.identity for sample in samples] == ["sample", "sample"]


def test_qualification_leaves_an_irreducible_collision_for_the_guard_to_find() -> None:
    """Two samples with one namespace *and* one local value come back still equal.

    That is deliberate: the rule resolves what it can and refuses to invent a distinction
    it does not have. ``duplicate_identity`` is what turns the leftover into an abort.
    """
    samples = [
        _discovered("sample", "/run1/job.zip", "/tmp/cohort_zip_aaa"),
        _discovered("sample", "/run2/job.zip", "/tmp/cohort_zip_bbb"),
    ]

    qualified = qualify_colliding_identities(samples)

    assert [sample.identity for sample in qualified] == ["job/sample", "job/sample"]


def test_qualification_is_idempotent_over_a_resolvable_collision() -> None:
    """Applying the rule twice must not qualify an already-qualified identity again.

    It is applied once today, at the end of discovery, so this is latent rather than live -
    and it is the kind of latency that is cheap to close and expensive to find later, since
    a second application is what any future caller would naturally reach for.
    """
    samples = [
        _discovered("sample", "/data/job_a.zip", "/tmp/cohort_zip_aaa"),
        _discovered("keeper", "/data/job_b.zip", "/tmp/cohort_zip_bbb"),
        _discovered("sample", "/elsewhere/job_c.zip", "/tmp/cohort_zip_ccc"),
    ]

    once = qualify_colliding_identities(samples)
    twice = qualify_colliding_identities(once)

    assert [sample.identity for sample in once] == ["job_a/sample", "keeper", "job_c/sample"]
    assert twice == once


def test_qualification_is_idempotent_over_an_irreducible_collision() -> None:
    """The case that actually broke: two `job/sample` results became `job/job/sample`.

    An irreducible pair comes back still equal by design - that is what leaves the abort to
    `duplicate_identity` - so a second application sees a shared identity again and, keying
    on the identity it had just written, qualified the qualification. The rule reads the
    sample's **local** value, which no application changes, so the second pass is a no-op.
    """
    samples = [
        _discovered("sample", "/run1/job.zip", "/tmp/cohort_zip_aaa"),
        _discovered("sample", "/run2/job.zip", "/tmp/cohort_zip_bbb"),
    ]

    once = qualify_colliding_identities(samples)
    twice = qualify_colliding_identities(once)

    assert [sample.identity for sample in once] == ["job/sample", "job/sample"]
    assert [sample.identity for sample in twice] == ["job/sample", "job/sample"]
    assert twice == once


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
    sample_on_disk(tmp_path / "cohort" / "sample_a")
    sample_on_disk(tmp_path / "cohort" / "sample_b")

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


def test_two_distinct_samples_are_reported_under_two_pseudonyms(tmp_path) -> None:
    """The control: the guard must not fire on the ordinary two-sample cohort."""
    sample_on_disk(tmp_path / "cohort" / "sample_one")
    sample_on_disk(tmp_path / "cohort" / "sample_two")
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
