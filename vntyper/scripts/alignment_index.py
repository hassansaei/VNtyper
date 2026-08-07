"""
vntyper/scripts/alignment_index.py

Module Purpose:
---------------
Where an alignment's index already is, decided without running samtools.

``process_bam_to_fastq`` reuses an existing BAM index when there is one and builds a new
one into the run's *output* directory when there is not - the input directory holds
patient data and is routinely mounted read-only (#162, #210). Only the second half of
that needs a subprocess. The first half is a pure question about filenames: which of the
two names htslib itself resolves exists on disk.

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

logger = logging.getLogger(__name__)


def resolve_bam_index(in_bam: str | Path) -> str | None:
    """Find an existing BAM index the way htslib itself does.

    htslib tries ``<file>.bai`` and then ``<stem>.bai``. The pipeline used to
    reconstruct only the first, so given the ``sample.bai`` that the upload
    endpoint and the worker both deliberately accept it found nothing and built a
    second index beside the input that nothing else knew about (#210).

    Args:
        in_bam (str | Path): The alignment whose index is wanted.

    Returns:
        str | None: The existing index, or None when there is none.
    """
    bam = Path(in_bam)
    for candidate in (Path(f"{bam}.bai"), bam.with_suffix(".bai")):
        if candidate.exists():
            return str(candidate)
    return None
