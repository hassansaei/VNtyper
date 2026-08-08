"""
vntyper/scripts/alignment_index.py

Module Purpose:
---------------
Where an alignment's index already is, decided without running samtools.

``process_bam_to_fastq`` reuses an existing BAM index when there is one and builds a new
one into the run's *output* directory when there is not - the input directory holds
patient data and is routinely mounted read-only (#162, #210). Only the second half of
that needs a subprocess. The first half is a pure question about filenames: which of the
two names a **BAI** can carry, ``<file>.bai`` or ``<stem>.bai``, exists on disk. CSI is
deliberately out of scope - see :func:`resolve_bam_index` for why the offset extractor
downstream makes that the correct answer rather than a gap.

That question lived in ``fastq_bam_processing.py``, where every other function shells out
through ``run_command``, so it could only be exercised by driving a stage that runs
samtools. It is here so it can be called on its own, and because #210's fix pushed that
file past the ~650-line guideline AGENTS.md rule 2 sets - the rule's answer to which is
to pull the pure logic out and leave the I/O behind.

Nothing about the resolution changed in the move. ``tests/unit/test_alignment_index.py``
pins it, and ``tests/unit/test_input_tree_is_never_written.py`` keeps the invariant it
serves: no index the ``vntyper`` package builds is written beside the input.
"""

from __future__ import annotations

import logging
from pathlib import Path

from vntyper.scripts.alignment_contract import index_candidate_names

logger = logging.getLogger(__name__)


def resolve_any_index(in_path: str | Path, file_format: str) -> str | None:
    """Find the first existing index accepted for an alignment format.

    Args:
        in_path (str | Path): The alignment whose index is wanted.
        file_format (str): The alignment format, either ``"bam"`` or ``"cram"``.

    Returns:
        str | None: The first existing candidate in the contract-defined order, or
        None when no candidate exists.
    """
    for candidate in index_candidate_names(str(in_path), file_format):
        if Path(candidate).exists():
            return candidate
    return None


def resolve_bam_index(in_bam: str | Path) -> str | None:
    """Find an existing **BAI** index, under either of the two names it can carry.

    Both spellings are tried, ``<file>.bai`` first and then ``<stem>.bai``. The
    pipeline used to reconstruct only the first, so given the ``sample.bai`` that
    the upload endpoint and the worker both deliberately accept it found nothing
    and built a second index beside the input that nothing else knew about (#210).

    **This is deliberately not htslib's resolution order, and must not be widened
    to match it.** htslib tries CSI before BAI - ``<file>.csi`` and ``<stem>.csi``
    - and a CSI is ignored here in both spellings. The only consumer of the
    returned path is ``extract_unmapped_from_offset.get_last_chunk_end``, which
    walks the BAI container itself and rejects any file whose first four bytes are
    not ``BAI\\x01``. Handing it a CSI would replace a working run with a
    ``ValueError`` mid-stage, whereas returning None makes the caller build the
    BAI it can actually read - into the run's *output* directory, never beside the
    input. Resolving CSI is a change to the offset extractor first and to this
    function second. Pinned by
    ``tests/unit/test_alignment_index.py::test_a_csi_index_is_ignored_and_a_bai_is_built_instead``.

    Args:
        in_bam (str | Path): The alignment whose index is wanted.

    Returns:
        str | None: The existing BAI, or None when there is none this function
        resolves - which includes an alignment indexed only by a CSI.
    """
    bam = Path(in_bam)
    for candidate in (Path(f"{bam}.bai"), bam.with_suffix(".bai")):
        if candidate.exists():
            return str(candidate)
    return None
