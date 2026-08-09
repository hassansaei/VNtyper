"""
vntyper/scripts/alignment_index.py

Module Purpose:
---------------
Where an alignment's index already is, decided without running samtools.

Alignment preflight enumerates existing indexes as protected patient inputs, but it does
not trust them for retrieval: BAI/CRAI/CSI carry no binding to the alignment bytes, and
round-3 measurement proved a valid wrong-sample BAI can return an empty target with exit
status zero. The trusted index is rebuilt beside the run-local view. The remaining pure
question here is which candidate name permitted by the alignment format exists, in the
order from :func:`alignment_contract.index_candidate_names`. BAM CSI is supported by
:func:`resolve_any_index`, as is CRAM CSI after pinned samtools/htslib acceptance was
measured directly. CSI is deliberately out of scope only for
:func:`resolve_bam_index`, which remains the legacy BAI-only preflight/protection
contract. The optional legacy offset parser remains independently tested but has no
production consumer.

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

    General resolution includes BAM and CRAM CSI. The legacy BAI-only BAM
    preflight/protection path uses :func:`resolve_bam_index` instead.

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
    - and a CSI is ignored here in both spellings. This function remains the legacy
    BAI-only preflight/protection contract; production indexed recovery uses the
    common htslib literal-``'*'`` fetch against the fresh run-local view. The
    optional legacy offset parser is independently tested but has no production
    consumer. Widening this resolver is a separate preflight-contract change. Pinned by
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
