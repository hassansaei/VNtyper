"""
vntyper/scripts/cohort_inputs.py

Module Purpose:
---------------
Locate the samples of a cohort and read what the pipeline recorded for each of them.

``vntyper cohort`` is pointed at whatever the caller has: a directory that is itself one
sample's output, a directory holding many, a zip file of either, or something that turns
out to be none of those. :func:`discover_sample_directories` reduces all of that to the
sorted list of directories holding a ``pipeline_summary.json``, extracting any zip files
on the way; :func:`load_pipeline_summary_for_sample` then reads one such directory into
the three pieces the report needs.

One failure mode is deliberate and is characterised rather than changed: a sample that
cannot be read - a missing summary, malformed JSON, an unparseable timestamp - is logged
and dropped so one bad sample cannot abort a 40-sample cohort.

**A sample's identity and its position in the output depend only on the inputs, never on
where they were unpacked** (#205). Both used to depend on the extraction directory, and a
zip is extracted into ``tempfile.mkdtemp(prefix="cohort_zip_")``:

* Discovery de-duplicates through a mapping and then **sorts**, because iterating that
  mapping in insertion order would follow the argument order and iterating a ``set`` -
  what it was before - followed ``Path.__hash__``, which is the hash of the path string
  and is randomised per process. The sort key is now
  :attr:`DiscoveredSample.sort_key`: one flat tuple, the parts of the input path the
  caller wrote followed by the sample's path relative to *that input's* root. The
  sample's extracted path was the old key, so two runs over the same archives ordered
  their rows differently.
* A zip whose ``pipeline_summary.json`` sits at the archive root - the layout the web
  worker produces - had the ``mkdtemp`` root itself as its sample directory, so the
  reported sample was a random ``cohort_zip_*`` string. :attr:`DiscoveredSample.identity`
  carries the name instead, read from the input file the run itself recorded.

Extraction still goes to a temporary directory; only its role in identity and order is
gone.

**An identity is a (namespace, value) pair, and only the value is usually written.** HL7
FHIR's ``Identifier`` puts uniqueness on ``system`` and ``value`` together - the value is
"unique within the context of the system" and says nothing on its own - and that is the
rule this module follows. A sample's namespace is its input's own name (the archive stem,
or the directory name); its local value is what the run recorded or what its directory is
called. Two web jobs that each uploaded ``sample.bam`` are one value in two systems, which
is two subjects. :func:`qualify_colliding_identities` therefore rewrites a *shared* local
value to ``namespace/value`` and leaves every unique one exactly as it was, and
:func:`duplicate_identity` aborts only on what remains equal after that.

The four step names this module matches are compared by exact string against what
``pipeline.py`` writes. A typo does not fail - it silently drops a section of the report
(AGENTS.md trap 5) - so they are imported from ``summary_steps``, never spelled out.

Extracted from ``cohort_summary.py`` in Task 22 of the #181-#197 follow-ups.
"""

import json
import logging
import os
import shutil
import tempfile
import zipfile
from collections import Counter
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Any

from vntyper.scripts.summary_steps import (
    STEP_ADVNTR,
    STEP_BAM_HEADER,
    STEP_COVERAGE,
    STEP_KESTREL,
)

logger = logging.getLogger(__name__)

#: The file every sample directory is identified by.
PIPELINE_SUMMARY_FILENAME = "pipeline_summary.json"

#: The ``input_files`` keys a zip-rooted sample's identity is read from, most specific
#: first. ``pipeline.py`` records the alignment it was handed under ``bam`` or ``cram`` and
#: the reads under ``fastq1``/``fastq2``; a run that aligned its own reads records both, and
#: the alignment is the file the genotype was called from.
IDENTITY_INPUT_KEYS = ("bam", "cram", "fastq1")

#: Stripped before ``Path.stem`` takes the last suffix, so ``patient_R1.fastq.gz``
#: identifies ``patient_R1`` rather than ``patient_R1.fastq``.
COMPRESSION_SUFFIXES = (".gz", ".bgz")

#: The container suffix an archive input carries. Stripped from the namespace because it
#: describes how the run was *packed*, not what the run was called: ``job_a.zip`` and a
#: directory ``job_a`` are the same system under two transports.
ARCHIVE_SUFFIX = ".zip"

#: Separates the namespace from the local value in a qualified identity. Neither half can
#: contain it - a namespace is one path component and a local value is a directory name or
#: a filename stem - so the pair is recoverable from the identity by splitting once.
NAMESPACE_SEPARATOR = "/"


@dataclass(frozen=True)
class DiscoveredSample:
    """One sample the cohort found, with an identity that does not move between runs."""

    #: Where the sample's ``pipeline_summary.json`` is right now. For a zip input this is
    #: under ``tempfile.mkdtemp``'s output and so differs on every run, which is precisely
    #: why it is neither the identity nor the sort key.
    directory: Path

    #: The name the sample is reported under, before any pseudonymisation. Discovery sets
    #: it to the sample's **local** value - what the run recorded, or the directory's name -
    #: and :func:`qualify_colliding_identities` rewrites it to ``namespace/value`` if and
    #: only if another discovered sample carries the same local value.
    identity: str

    #: The sample's **local** value, kept as discovery found it and never rewritten. It is
    #: what :func:`qualify_colliding_identities` compares and what it qualifies, which is
    #: what makes that rule idempotent: reading :attr:`identity` instead meant a second
    #: application saw ``job/sample`` twice, took it for a shared local value and qualified
    #: it again into ``job/job/sample``. The rule is applied once today, so keeping the two
    #: apart costs one field and removes the whole class of mistake.
    local_identity: str

    #: The input path this sample was reached through, exactly as the caller wrote it - a
    #: directory or a zip archive. Kept for diagnostics only: it is the one thing that can
    #: name a zip sample's origin in a message, because its ``directory`` is an extraction
    #: directory. It is **not** part of the ordering - see :attr:`sort_key`.
    origin: str

    #: What the ordering is total on: the sample's **effective path**, a single flat tuple
    #: of the parts of :attr:`origin` followed by the parts of the sample's path relative
    #: to that input's root. Every component is user-supplied or archive-internal, so none
    #: can be a ``tempfile.mkdtemp`` name and two runs agree.
    #:
    #: Flat rather than the pair ``(origin parts, relative parts)`` it briefly was, and the
    #: difference shows when one input is a prefix of another. Given ``cohort/a`` as a
    #: direct input and ``cohort`` after it, the pair compares on its outer half alone -
    #: ``("cohort",) < ("cohort", "a")`` - so the parent's samples all sort before the
    #: child's and ``cohort/z`` came out ahead of ``cohort/a``. Concatenating instead makes
    #: the key the sample's own location whichever input reached it, which is the
    #: whole-path order this had before the key existed and the order ``ls`` shows.
    sort_key: tuple[str, ...]


def _stem_of_recorded_input(recorded: str) -> str:
    """Reduce a recorded input filename to the name a clinician would recognise.

    Args:
        recorded: A value out of a summary's ``input_files``, e.g. ``patient1.bam``.

    Returns:
        str: The stem, with one compression suffix stripped first so
        ``patient_R1.fastq.gz`` yields ``patient_R1``.
    """
    name = Path(recorded).name
    for suffix in COMPRESSION_SUFFIXES:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return Path(name).stem


def _identity_from_summary(sample_dir: Path, fallback: str) -> str:
    """Name a sample from the input file its own run recorded.

    The directory name is the right identity almost everywhere, but not for a zip whose
    ``pipeline_summary.json`` sits at the root: that directory is the randomised
    ``tempfile.mkdtemp`` root, so the reported sample differed on every run (#205). That
    layout is what the web worker produces, so it is the normal path for web cohorts rather
    than an edge case.

    Args:
        sample_dir: The directory holding ``pipeline_summary.json``.
        fallback: The identity to use when the summary records no input files - the
            archive's filename stem.

    Returns:
        str: The sample's identity. Anything that goes wrong yields ``fallback``; this
        helper never raises, matching the module's "one bad sample must not abort the
        cohort" contract.
    """
    summary_path = sample_dir / PIPELINE_SUMMARY_FILENAME
    if not summary_path.exists():
        return fallback
    try:
        with open(summary_path) as handle:
            summary = json.load(handle)
        input_files = summary.get("input_files") or {}
        for key in IDENTITY_INPUT_KEYS:
            recorded = input_files.get(key)
            if recorded:
                return _stem_of_recorded_input(str(recorded))
    except Exception as e:
        logger.warning(f"Could not read the recorded inputs of {sample_dir}, naming it {fallback}: {e}")
    return fallback


def _samples_under_root(root: Path, origin: str, root_identity: str) -> list[DiscoveredSample]:
    """Find the samples one already-resolved input root contains.

    A directory holding a ``pipeline_summary.json`` is taken as one sample and is not
    searched further; any other directory is searched recursively.

    Args:
        root: The input directory, or the directory a zip was extracted into.
        origin: The input path as the caller wrote it. Its parts lead the sort key, so a
            zip's samples sort under the archive's own path rather than under the
            extraction directory.
        root_identity: The identity to give a sample sitting at ``root`` itself. For a
            directory that is the directory's name; for a zip it is what the run recorded.

    Returns:
        list[DiscoveredSample]: The samples found, unsorted. Each carries an effective
        path as its sort key: ``origin``'s parts followed by the sample's parts below
        ``root``. For a directory input the two concatenate back to the sample's own path,
        which is what makes the order independent of whether the caller named the sample
        or one of its ancestors; for a zip they compose the archive's path with the
        member's path inside it, and the extraction directory appears in neither.
    """
    origin_parts = Path(origin).parts
    if (root / PIPELINE_SUMMARY_FILENAME).exists():
        logger.info(f"Found {PIPELINE_SUMMARY_FILENAME} in {root}")
        return [
            DiscoveredSample(
                directory=root,
                identity=root_identity,
                local_identity=root_identity,
                origin=origin,
                sort_key=origin_parts,
            )
        ]

    samples: list[DiscoveredSample] = []
    for summary_file_path in root.rglob(PIPELINE_SUMMARY_FILENAME):
        sample_dir = summary_file_path.parent
        logger.info(f"Found {PIPELINE_SUMMARY_FILENAME} in {sample_dir}")
        samples.append(
            DiscoveredSample(
                directory=sample_dir,
                identity=sample_dir.name,
                local_identity=sample_dir.name,
                origin=origin,
                sort_key=origin_parts + sample_dir.relative_to(root).parts,
            )
        )
    return samples


def identity_namespace(origin: str) -> str:
    """Name the namespace an input establishes for the identities discovered under it.

    HL7 FHIR's ``Identifier`` splits an identifier into a ``system`` - "an absolute URL
    that describes a set [of] values that are unique" - and a ``value``, "the portion of
    the identifier ... unique **within the context of the system**". A cohort input is the
    system: everything one archive or one directory contains came out of one run or one
    upload, so a name is unique inside it and means nothing outside it.

    The **name** of the input is used, not its path. That is what makes ``-i /data/run1``
    and ``-i ./run1`` the same system, so a qualified identity - and the pseudonym derived
    from it - survives the cohort being moved or re-run from a different working directory.
    The cost is that two inputs whose names are equal are one namespace; see
    :func:`duplicate_identity` for what happens then.

    Args:
        origin: The input path exactly as the caller wrote it, i.e.
            :attr:`DiscoveredSample.origin`.

    Returns:
        str: The archive stem for a ``.zip`` input and the directory name otherwise. An
        input naming no component of its own (``.``, ``/``) has no name to give, so the
        path itself is returned rather than an empty namespace.
    """
    name = Path(origin).name
    if name.lower().endswith(ARCHIVE_SUFFIX):
        name = name[: -len(ARCHIVE_SUFFIX)]
    return name or origin


def qualify_colliding_identities(samples: list[DiscoveredSample]) -> list[DiscoveredSample]:
    """Qualify by namespace exactly those identities that more than one sample shares.

    A repeated local value is **not** a repeated identifier. ``sample.bam`` is a value, and
    two web jobs that each uploaded one are that value in two different systems - two
    subjects, not one. Comparing values while discarding the system either merges the two
    (what happened before their identity came out of discovery) or refuses to run (what
    happened after); both are wrong for the same reason, and ``sample.bam`` is common
    enough that a root-level ZIP cohort meets it on the normal path.

    Only the samples that actually collide are qualified. A sample whose local value is
    already unique keeps it untouched, so its identity, its row and its pseudonym are
    exactly what they would have been in a cohort with no collision in it at all - stable
    whichever other samples it is aggregated with.

    This is a pure function of the discovered list and is applied at the end of
    :func:`discover_sample_directories`, so every consumer sees one set of identities.

    **It is idempotent**, because it reads and writes different fields:
    :attr:`DiscoveredSample.local_identity` is what it compares and qualifies, and only
    :attr:`DiscoveredSample.identity` is what it writes. Reading the identity it had just
    written made a second application see an irreducible pair's two ``job/sample`` results
    as a shared local value and qualify them again into ``job/job/sample``. Only one
    application happens today, so that was latent; it is closed here rather than left as a
    trap for the next caller.

    Args:
        samples: The samples discovery found, in the order it found them.

    Returns:
        list[DiscoveredSample]: The same samples in the same order. A sample whose local
        value is shared carries ``f"{namespace}/{local}"`` as its identity instead;
        everything else about it, and every other sample, is returned unchanged.
        ``DiscoveredSample`` is frozen, so the replacements are new records and the
        caller's list is not touched.

        Two samples sharing a namespace *and* a local value come back still equal: the rule
        resolves what the inputs distinguish and invents nothing. :func:`duplicate_identity`
        is what turns that leftover into an abort.
    """
    occurrences = Counter(sample.local_identity for sample in samples)
    return [
        replace(sample, identity=f"{identity_namespace(sample.origin)}{NAMESPACE_SEPARATOR}{sample.local_identity}")
        if occurrences[sample.local_identity] > 1
        else sample
        for sample in samples
    ]


def discover_sample_directories(input_paths: list[str]) -> tuple[list[DiscoveredSample], list[str]]:
    """Resolve the cohort's input paths to the samples they contain.

    A directory holding a ``pipeline_summary.json`` is taken as one sample and is not
    searched further; any other directory is searched recursively. A zip file is
    extracted to a temporary directory and then treated the same way. Anything else -
    a path that does not exist, a file that is not a zip - is logged and skipped, so a
    cohort of forty runs is not lost to one bad argument.

    Args:
        input_paths: Directories or zip files, as ``vntyper cohort`` received them.

    Returns:
        tuple[list[DiscoveredSample], list[str]]: The samples in **sorted** order, and
        the temporary directories the caller must pass to :func:`cleanup_temp_dirs`
        once the report has been written.

        Samples are accumulated in a mapping keyed on the sample directory **as
        ``Path.resolve()`` reports it**, and repeated inputs are dropped by the same
        canonical comparison before anything is extracted; see the two paragraphs after
        the sort. They are sorted on
        :attr:`DiscoveredSample.sort_key` on the way out - the sample's **effective path**,
        the parts of the input path the caller wrote followed by the sample's path
        relative to that input's own root. No component of it can move between runs. The
        key used to be the sample's full *extracted* path, whose leading component for a
        zip is ``tempfile.mkdtemp``'s random suffix - so two runs over the same archives
        sorted their samples into different positions (#205), and every CSV/TSV/JSON row
        order followed. Composing from the input path rather than from its position in
        ``input_paths`` is deliberate: a cohort of directories keeps the lexicographic
        order it has always had, and only the temporary component goes away.

        The two halves are **concatenated rather than nested**. Nested, the comparison is
        decided by the input path alone whenever two inputs differ, which reorders a
        cohort whose inputs nest: ``cohort/a`` named directly and ``cohort`` named after it
        put ``cohort/z`` first, because ``("cohort",) < ("cohort", "a")`` and the relative
        half was never reached. Concatenated, a sample keys on its own location whichever
        input found it. :attr:`DiscoveredSample.origin` is kept for diagnostics and for the
        duplicate-identity message, and takes no part in the order.

        **De-duplication is on the physical location, not on the spelling.** The key used
        to be the ``Path`` exactly as it was constructed, so one sample reached through an
        absolute parent and again through a relative direct input - ``vntyper cohort -i
        "$PWD/cohort" -i cohort/sample`` - produced two keys that are unequal as values and
        identical on disk, and the cohort reported one patient as two samples. A symlinked
        input did the same. ``Path.resolve()`` settles all three, and the record still
        carries :attr:`DiscoveredSample.origin` exactly as the caller wrote it, because
        canonicalisation is for identity and a diagnostic must quote what was typed.

        **A repeated input is dropped before it is extracted.** A zip's sample directory is
        its own fresh ``tempfile.mkdtemp`` root, so two extractions of one archive can never
        compare equal however the key is canonicalised: ``-i job.zip -i job.zip`` wrote both
        copies to disk and then aborted the cohort, because the two samples shared a
        namespace *and* a local value and that is exactly what :func:`duplicate_identity`
        refuses. The canonical form of each input path is therefore remembered and a repeat
        is skipped outright, which also spares the second extraction.

        Identities are then passed through :func:`qualify_colliding_identities`, so what
        comes out is the cohort's final naming and every consumer agrees on it. A local
        value shared by two or more samples becomes ``namespace/value`` for those samples
        only; a unique one is untouched.

        Tuples of path *parts* are compared element by element, so the separator never
        participates and ``cohort/sample`` sorts before ``cohort-extra/sample`` - the
        order ``ls`` would show. Sorting here rather than at the one call site means every
        consumer gets the same order. Pinned by
        ``tests/unit/test_cohort_inputs.py::test_the_discovered_directories_come_back_sorted``,
        ``::test_an_input_nested_inside_a_later_input_keeps_its_whole_path_position``
        and by the cross-process tests beside them.

    Raises:
        Exception: Whatever the work after an extraction raises, re-raised unchanged after
            the archives extracted so far have been removed. The caller's cleanup ``try``
            can only open on the ``temp_dirs`` this function *returns*, so nothing it does
            covers the region between the first extraction and that return - and this
            function extracts every remaining input, sorts, and qualifies in there. Owning
            the cleanup until the handles are handed over is the only place it can be done.
    """
    temp_dirs: list[str] = []
    # Keyed on the resolved directory, so a sample reached by two spellings of one path -
    # absolute and relative, through a symlink, or through a parent and again directly -
    # appears once, and the first input that named it decides its place in the order.
    processed_dirs: dict[Path, DiscoveredSample] = {}
    # The canonical form of every input already handled. A zip's samples land under a fresh
    # mkdtemp root, so a repeated archive is only de-duplicable *before* it is extracted.
    seen_origins: set[Path] = set()

    try:
        for path_str in input_paths:
            path = Path(path_str)
            found: list[DiscoveredSample] = []
            if not path.exists():
                logger.warning(f"Input path does not exist and will be skipped: {path}")
                continue
            canonical_origin = path.resolve()
            if canonical_origin in seen_origins:
                logger.info(f"Input path already processed under another spelling, skipping: {path}")
                continue
            seen_origins.add(canonical_origin)
            if path.is_dir():
                found = _samples_under_root(path, path_str, path.name)
                if not found:
                    logger.warning(f"No {PIPELINE_SUMMARY_FILENAME} found in directory {path}")
            elif zipfile.is_zipfile(path):
                logger.info(f"Extracting zip file: {path}")
                temp_dir = tempfile.mkdtemp(prefix="cohort_zip_")
                try:
                    with zipfile.ZipFile(path, "r") as zip_ref:
                        zip_ref.extractall(temp_dir)
                    temp_path = Path(temp_dir)
                    # The archive's own stem is the fallback identity: it is the only thing
                    # about a zip-rooted sample that is not the extraction directory.
                    found = _samples_under_root(temp_path, path_str, _identity_from_summary(temp_path, path.stem))
                    if not found:
                        logger.warning(f"No {PIPELINE_SUMMARY_FILENAME} found in extracted zip file: {path}")
                    temp_dirs.append(temp_dir)
                except zipfile.BadZipFile as e:
                    logger.error(f"Bad zip file {path}: {e}")
                    shutil.rmtree(temp_dir)
                    continue
                except Exception as e:
                    logger.error(f"Error extracting zip file {path}: {e}")
                    shutil.rmtree(temp_dir)
                    continue
            else:
                logger.warning(f"Unsupported file type (not a directory or zip): {path}")
                continue

            for sample in found:
                processed_dirs.setdefault(sample.directory.resolve(), sample)

        # Sorted on the origin key, not on the extracted path: see the Returns section above.
        ordered = sorted(processed_dirs.values(), key=lambda sample: sample.sort_key)
        # Qualification is last, and over the whole list, because "shared" is a property of the
        # cohort rather than of any one sample: nothing can be qualified until every sample is
        # known. It does not touch sort_key, so the order above survives it.
        return qualify_colliding_identities(ordered), temp_dirs
    except Exception:
        # Everything extracted so far is unreachable once this frame unwinds: the caller has
        # no handle on it, because the handles are the return value it never receives. The
        # two `continue`s above are the extraction failures that clean up after themselves;
        # this covers every other way out of the region.
        cleanup_temp_dirs(temp_dirs)
        raise


def duplicate_identity(samples: list[DiscoveredSample]) -> tuple[DiscoveredSample, DiscoveredSample] | None:
    """Find the first two samples reporting the same identity after qualification.

    Widening the pseudonym digest cannot make the reported sample injective, because the
    *inputs* to it can already be equal. ``cohort_categories.sample_categories`` groups on
    the reported sample, so an undetected pair is counted as one sample and one of the two
    genotypes is lost - with or without pseudonymisation (#206).

    What reaches this guard is **narrower than it used to be**, because
    :func:`qualify_colliding_identities` has already run: a shared local value is no longer
    a collision, only a shared ``(namespace, value)`` pair is. That leaves the single case
    the inputs genuinely cannot separate - two archives with the same stem, or two
    directory inputs with the same name - where qualification produces one identity twice
    and there is nothing left to distinguish the two samples by. Refusing to run is the
    only honest answer; the message names both inputs so the operator can supply one.

    Discovery de-duplicates on the sample directory's canonical path, so two entries here
    always come from two physically distinct directories and no further comparison is
    needed. That canonicalisation is load-bearing for this guard in both directions: one
    directory named twice used to arrive as two records - which an operator would meet as
    an abort naming one archive against itself - and one *sample* reached by two spellings
    used to arrive as two records with two qualified identities, which this guard cannot
    see at all.

    Args:
        samples: The samples :func:`discover_sample_directories` returned.

    Returns:
        tuple[DiscoveredSample, DiscoveredSample] | None: The first colliding pair in
        discovery order, or None when every identity is distinct.
    """
    seen: dict[str, DiscoveredSample] = {}
    for sample in samples:
        previous = seen.get(sample.identity)
        if previous is not None:
            return previous, sample
        seen[sample.identity] = sample
    return None


def cleanup_temp_dirs(temp_dirs: list[str]) -> None:
    """Remove the directories the zip inputs were extracted into.

    Called after the report has been written, so a directory that will not delete is
    logged rather than raised - the cohort's output is already on disk by then.

    Args:
        temp_dirs: The temporary directories from :func:`discover_sample_directories`.
    """
    for temp_dir in temp_dirs:
        try:
            shutil.rmtree(temp_dir)
            logger.debug(f"Cleaned up temporary directory: {temp_dir}")
        except Exception as e:
            logger.error(f"Failed to remove temporary directory {temp_dir}: {e}")


def parse_pipeline_summary(summary: dict[str, Any]) -> tuple[list[dict], list[dict], dict[str, Any]]:
    """Extract the cohort's three inputs from one parsed ``pipeline_summary.json``.

    Steps this cohort does not consume are ignored, and a step recorded more than once
    leaves the last occurrence in place.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.

    Returns:
        tuple[list[dict], list[dict], dict[str, Any]]: The Kestrel rows, the adVNTR
        rows, and the per-sample statistics (``runtime``, ``version``, ``assembly``,
        ``pipeline``, ``coverage``).

    Raises:
        ValueError: If ``pipeline_start`` or ``pipeline_end`` is present but is not an
            ISO-8601 timestamp. :func:`load_pipeline_summary_for_sample` catches this
            and drops the sample.
    """
    kestrel_data: list[dict] = []
    advntr_data: list[dict] = []
    additional_stats: dict[str, Any] = {}

    # Compute runtime from top-level timestamps if available
    pipeline_start = summary.get("pipeline_start")
    pipeline_end = summary.get("pipeline_end")
    if pipeline_start and pipeline_end:
        start_dt = datetime.fromisoformat(pipeline_start)
        end_dt = datetime.fromisoformat(pipeline_end)
        runtime_sec = (end_dt - start_dt).total_seconds()
        additional_stats["runtime"] = f"{runtime_sec:.2f} seconds"
    else:
        additional_stats["runtime"] = "N/A"

    # Pipeline version from top-level field
    additional_stats["version"] = summary.get("version", "N/A")

    # Initialize defaults for assembly, pipeline and coverage
    additional_stats["assembly"] = "N/A"
    additional_stats["pipeline"] = "N/A"
    additional_stats["coverage"] = {}

    for step in summary.get("steps", []):
        if step.get("step") == STEP_KESTREL:
            kestrel_data = step.get("parsed_result", {}).get("data", [])
        elif step.get("step") == STEP_ADVNTR:
            advntr_data = step.get("parsed_result", {}).get("data", [])
        elif step.get("step") == STEP_BAM_HEADER:
            parsed = step.get("parsed_result", {})
            additional_stats["assembly"] = parsed.get("assembly_text", "N/A")
            additional_stats["pipeline"] = parsed.get("alignment_pipeline", "N/A")
        elif step.get("step") == STEP_COVERAGE:
            parsed = step.get("parsed_result", {})
            data_list = parsed.get("data", [])
            if data_list:
                additional_stats["coverage"] = data_list[0]
    return kestrel_data, advntr_data, additional_stats


def load_pipeline_summary_for_sample(
    sample_dir: str | os.PathLike[str],
) -> tuple[list[dict], list[dict], dict[str, Any]]:
    """
    Load the pipeline_summary.json from a sample directory and extract Kestrel,
    adVNTR data and additional statistics (runtime, coverage, version, assembly, pipeline).

    For the adVNTR step, the algorithm result will later be computed based on the logic
    defined in report_config.json.

    Anything that goes wrong - the file is absent, the JSON does not parse, a timestamp
    does not - is logged and yields three empties, so one unreadable sample does not
    abort the cohort.

    Parameters
    ----------
    sample_dir : str or Path
        Directory containing the pipeline_summary.json file.

    Returns
    -------
    tuple
        Three elements: (kestrel_data, advntr_data, additional_stats).
        additional_stats is a dict with keys:
          - runtime: pipeline run duration (in seconds)
          - version: pipeline version
          - assembly: assembly text from BAM Header Parsing
          - pipeline: alignment pipeline from BAM Header Parsing
          - coverage: coverage metrics dict (mean, median, stdev, min, max)
    """
    sample_dir = Path(sample_dir)
    summary_path = sample_dir / PIPELINE_SUMMARY_FILENAME
    if not summary_path.exists():
        logger.warning(f"Pipeline summary file not found in {sample_dir}")
        return [], [], {}
    try:
        with open(summary_path) as f:
            summary = json.load(f)
        return parse_pipeline_summary(summary)
    except Exception as e:
        logger.error(f"Error loading pipeline summary from {sample_dir}: {e}")
        return [], [], {}
