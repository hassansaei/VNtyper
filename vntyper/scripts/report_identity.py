"""
report_identity.py

Module Purpose:
---------------
Who and what the per-sample HTML report is about: the sample name in its
``<title>`` and ``<h1>``, the input files it was run on, and the provenance rows
that say which assembly and which span the run actually used.

Every report a lab produced used to be titled ``Summary Report`` -- the same
constant string for every run -- and the template branched on two of the four
input shapes :func:`~vntyper.scripts.pipeline_inputs.resolve_pipeline_input`
emits, so a CRAM run and a single-end FASTQ run rendered an empty ``Input
Files:`` line (#242).

Like :mod:`vntyper.scripts.report_formatting` and
:mod:`vntyper.scripts.screening_summary`, this is the pure half: values in,
strings out, no filesystem and no pandas. It is a module of its own rather than
more of ``report_formatting`` (793 lines) because identity and provenance are a
different concern from table cell formatting, and AGENTS.md rule 3 asks for the
seam rather than the growth.

Two rules govern everything here and they point the same way -- **say what was
recorded, or say that nothing was**:

* The sample name is derived only from shapes the rule recognises. Anything else
  -- an unrecognised extension, a mate marker that is not trailing, two mates
  whose stems disagree -- is printed **verbatim**. A wrong sample name on a
  report is worse than an ugly one.
* A provenance value the run did not record renders :data:`NOT_RECORDED`, never
  a configured default. Reading ``config["default_values"]["reference_assembly"]``
  would mislabel every ``--reference-assembly`` override and cannot reconstruct
  ``--custom-regions`` at all.

Functions:
    derive_sample_name: One input basename to the sample it names
    resolve_sample_name: An explicit name, or the fallback over the run's inputs
    input_file_names: The file names an ``input_files`` mapping records
    format_input_files: The ``input_files`` mapping as one displayed line
    format_region: A resolved region span, with thousands separators
    format_run_timestamp: A recorded ISO timestamp as a labelled UTC time
    recorded_or_not: A recorded value, or the words that say it is absent
"""

from __future__ import annotations

import logging
import os
import re
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

#: What every provenance row says when the run recorded nothing for it. This is a
#: fallback, never a guess: the alternative -- printing the configured default --
#: makes a report that ran with ``--reference-assembly hg38`` claim ``hg19``.
NOT_RECORDED = "not recorded by this run"

#: The sample name for a run whose summary names no input files at all. A blank
#: title names nobody, which is the defect rather than a smaller version of it.
UNNAMED_SAMPLE = "unnamed sample"

#: What this report is a report of. Descriptive, not a claim about its use:
#: VNtyper is research use only and the interpretive text is config-driven.
ASSAY_NAME = "MUC1 coding VNTR genotyping"

#: The fixed half of the report's title; the sample name follows it.
REPORT_TITLE_PREFIX = "MUC1 VNTR report"

#: Compound extensions the sample-name rule recognises, longest first so that
#: ``.fastq.gz`` is matched before ``.fastq`` could be. A basename ending in
#: anything else is not something this pipeline reads, so the rule declines to
#: strip it rather than guessing where the name ends.
INPUT_EXTENSIONS: tuple[str, ...] = (
    ".fastq.gz",
    ".fastq.bz2",
    ".fq.gz",
    ".fq.bz2",
    ".fastq",
    ".fq",
    ".cram",
    ".bam",
    ".sam",
)

#: Mate markers stripped from the end of a stem, and only from the end. An
#: implementation that looked for ``_R1`` anywhere would rename
#: ``S1_S1_L001_R1_001`` to ``S1_S1_L001_001``, which is not what any other
#: artefact of the run calls that sample.
MATE_SUFFIXES: tuple[str, ...] = ("_R1", "_R2")

#: A samtools-style region. Only a span is reformatted; a custom region name is
#: printed exactly as the run recorded it.
REGION_PATTERN = re.compile(r"^(?P<contig>[^\s:]+):(?P<start>\d+)-(?P<end>\d+)$")

#: How a recorded run timestamp is displayed. ``start_summary`` writes a naive
#: ISO string that *is* UTC, so saying UTC is a statement of fact rather than an
#: assumption -- and it is what makes the time comparable with anything else.
RUN_TIMESTAMP_FORMAT = "%Y-%m-%d %H:%M:%S UTC"


def derive_sample_name(basename: str) -> str:
    """Derive the sample a single input file names.

    Strips one recognised compound extension (:data:`INPUT_EXTENSIONS`) and then
    a single **trailing** ``_R1`` or ``_R2``. Anything the rule does not
    recognise comes back unchanged.

    Args:
        basename: The input file's basename. A path is accepted and its
            directory components are dropped.

    Returns:
        str: The derived sample name, the basename verbatim when the rule does
        not recognise its shape, or :data:`UNNAMED_SAMPLE` for an empty name.
    """
    name = os.path.basename(str(basename)).strip()
    if not name:
        return UNNAMED_SAMPLE

    lowered = name.lower()
    stem: str | None = None
    for extension in INPUT_EXTENSIONS:
        if lowered.endswith(extension) and len(name) > len(extension):
            stem = name[: -len(extension)]
            break
    if stem is None:
        logger.debug("%r is not a recognised input file name; using it as the sample name.", name)
        return name

    for suffix in MATE_SUFFIXES:
        if stem.endswith(suffix) and len(stem) > len(suffix):
            return stem[: -len(suffix)]
    return stem


def resolve_sample_name(explicit: str | None, *sources: str | Path) -> str:
    """Choose the name the report calls its sample.

    An explicit ``--sample-name`` always wins. Otherwise every source is put
    through :func:`derive_sample_name` and the result is used only when they all
    agree; mates whose stems disagree fall back to the first basename verbatim,
    because naming one of the two would present a guess as a decision.

    Args:
        explicit: The ``--sample-name`` value, or None. A blank string is not a
            name and does not win.
        *sources: The run's input files, in the order the summary records them.

    Returns:
        str: The sample name, or :data:`UNNAMED_SAMPLE` when there is nothing to
        derive one from.
    """
    if explicit is not None and str(explicit).strip():
        return str(explicit).strip()

    basenames = [os.path.basename(str(source)).strip() for source in sources]
    basenames = [name for name in basenames if name]
    if not basenames:
        logger.debug("No input files to derive a sample name from; the report is titled %r.", UNNAMED_SAMPLE)
        return UNNAMED_SAMPLE

    derived = {derive_sample_name(name) for name in basenames}
    if len(derived) == 1:
        return derived.pop()

    logger.info(
        "Input file names %s do not agree on a sample name; using %r verbatim rather than choosing between them.",
        basenames,
        basenames[0],
    )
    return basenames[0]


def input_file_names(input_files: object) -> list[str]:
    """The file names a summary's ``input_files`` records, in recorded order.

    The mapping is **iterated**, never branched on. The template used to test for
    ``fastq1 and fastq2`` or ``bam``, which is why the two shapes nobody wrote a
    branch for -- CRAM and single-end FASTQ -- rendered an empty line; a fifth
    shape added to ``resolve_pipeline_input`` would have done the same.

    Args:
        input_files: The summary's ``input_files``. ``vntyper report`` and
            ``vntyper cohort`` both consume a *supplied* ``pipeline_summary.json``
            (#207), so this is untrusted: anything that is not a mapping yields no
            names rather than raising. An empty report section is a usable
            diagnostic; a traceback in place of the report is not.

    Returns:
        list[str]: The non-empty file names, in the order the summary records them.
    """
    if not isinstance(input_files, Mapping):
        logger.error(
            "Summary `input_files` is %s, not a mapping; the report can name no input files.",
            type(input_files).__name__,
        )
        return []
    return [str(value).strip() for value in input_files.values() if str(value).strip()]


def format_input_files(input_files: object) -> str:
    """List the files a run was given, whatever shape its ``input_files`` has.

    Args:
        input_files: The summary's ``input_files``, keyed by input kind.

    Returns:
        str: The file names in recorded order, comma-separated, or
        :data:`NOT_RECORDED` when the summary names none.
    """
    names = input_file_names(input_files)
    if not names:
        return NOT_RECORDED
    return ", ".join(names)


def format_region(region: str | None) -> str:
    """Render the region the run actually analysed.

    Args:
        region: The recorded ``region_resolved``, or None.

    Returns:
        str: A span with thousands separators, the recorded text unchanged when
        it is not a span (a ``--custom-regions`` name, say), or
        :data:`NOT_RECORDED` when nothing was recorded.
    """
    if region is None or not str(region).strip():
        return NOT_RECORDED
    text = str(region).strip()
    match = REGION_PATTERN.match(text)
    if match is None:
        logger.debug("Recorded region %r is not a span; printing it as recorded.", text)
        return text
    return f"{match['contig']}:{int(match['start']):,}-{int(match['end']):,}"


def format_run_timestamp(value: str | None) -> str:
    """Render the time the run started, as distinct from the time this file was made.

    ``report_date`` is ``datetime.now()`` at render, so re-running ``vntyper
    report`` over an archived run restamped it and the artefact silently claimed
    to be newer than the analysis. The report now prints both, and this is the
    half that belongs to the run.

    Args:
        value: The summary's ``pipeline_start``, or None.

    Returns:
        str: The timestamp in UTC, the recorded text unchanged when it cannot be
        parsed, or :data:`NOT_RECORDED` when nothing was recorded.
    """
    if value is None or not str(value).strip():
        return NOT_RECORDED
    text = str(value).strip()
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        logger.debug("Recorded pipeline start %r is not an ISO timestamp; printing it as recorded.", text)
        return text
    if parsed.tzinfo is not None:
        parsed = parsed.astimezone(timezone.utc).replace(tzinfo=None)
    return parsed.strftime(RUN_TIMESTAMP_FORMAT)


def recorded_or_not(value: object) -> str:
    """Render a recorded provenance value, or say that it was not recorded.

    Args:
        value: The value read out of the summary.

    Returns:
        str: The value as text, or :data:`NOT_RECORDED` when it is absent or blank.
    """
    if value is None:
        return NOT_RECORDED
    text = str(value).strip()
    return text or NOT_RECORDED
