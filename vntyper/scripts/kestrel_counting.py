"""Attempt-owned layout and cleanup for the standalone KAnalyze counting step.

``kestrel_genotyping.py`` is well over the ~650-line guideline, so the lifecycle policy
the split introduces lives here rather than growing it further (AGENTS.md rule 3). This
module is pure path arithmetic plus one removal, so it is directly testable; the loop
and the subprocess calls stay in ``run_kestrel``.

Two decisions are encoded here and neither is arbitrary.

**An attempt owns a directory, not a filename.** ``kestrel_kmers_k20.ikc`` beside the
other outputs still collides between two concurrent or retried attempts at the same k,
and KAnalyze's offloaded segment files land next to whatever it is told to write. A
directory per attempt makes the collision structurally impossible and makes cleanup one
recursive removal rather than a glob nobody re-checks.

**Cleanup never raises.** It runs in a ``finally`` on every exit path -- count failure
leaving a partial IKC, call failure, and success -- so an error here would replace the
count or call failure that caused it with a filesystem error about a temporary file.
That is the same precedence ``pipeline_cleanup.close_alignment_plan`` keeps, and it is
the opposite of ``_discard_attempt_artifacts``, which *must* raise: carrying a stale VCF
or SAM into the next k-mer size would report it against the wrong one, while a leftover
IKC is only wasted disk.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)

#: The IKC file name inside an attempt directory. Kestrel adopts a supplied count file
#: when its format is ``ikc`` or when the format is ``auto`` and the name matches
#: ``.*\.ikc``; VNtyper passes ``-f ikc`` explicitly, but the suffix keeps both routes
#: open and makes the file self-describing on disk.
IKC_FILENAME = "kestrel_kmers.ikc"

#: Where ``kanalyze.jar`` ships. The single source for this default, so the runner and
#: the provenance recorder cannot disagree about whether a run will split.
DEFAULT_KANALYZE_PATH = "vntyper/dependencies/kestrel/kanalyze.jar"


def attempt_directory(output_dir: str | Path, kmer_size: int) -> Path:
    """Return the directory one k-mer attempt owns.

    Args:
        output_dir: The Kestrel stage's output directory.
        kmer_size: The k-mer size this attempt is for.

    Returns:
        Path: ``<output_dir>/kmer_<size>``.
    """
    return Path(output_dir) / f"kmer_{kmer_size}"


def ikc_path(output_dir: str | Path, kmer_size: int) -> Path:
    """Return the IKC path one k-mer attempt writes.

    Args:
        output_dir: The Kestrel stage's output directory.
        kmer_size: The k-mer size this attempt is for.

    Returns:
        Path: The attempt's IKC file.
    """
    return attempt_directory(output_dir, kmer_size) / IKC_FILENAME


def count_log_path(output_dir: str | Path, kmer_size: int) -> Path:
    """Return the log destination for one attempt's counting step.

    Separate from the Kestrel call's log on purpose: ``run_command`` opens its log with
    ``"w"``, so reusing one file for both steps would erase the counting diagnostics
    with the call's output -- exactly when counting is what failed.

    Args:
        output_dir: The Kestrel stage's output directory.
        kmer_size: The k-mer size this attempt is for.

    Returns:
        Path: ``<output_dir>/kanalyze_count_kmer_<size>.log``.
    """
    return Path(output_dir) / f"kanalyze_count_kmer_{kmer_size}.log"


def remove_count_artifacts(attempt_dir: str | Path) -> None:
    """Remove one attempt's counting artefacts, never masking a primary error.

    Kestrel deliberately does not delete an IKC it was handed --
    ``IkcCountMap.preModuleRun`` sets ``rmLastTemp = false`` for an adopted file -- so
    nothing else will. Removing the whole attempt directory also takes KAnalyze's
    offloaded segment files with it, which a filename glob would have to keep in step
    with KAnalyze's own naming.

    This runs in a ``finally`` on every exit path, so it must not raise: a failure to
    remove a temporary file must not replace the count or call failure that caused it.

    Args:
        attempt_dir: The directory returned by :func:`attempt_directory`.
    """
    try:
        shutil.rmtree(attempt_dir)
    except FileNotFoundError:
        # Never created, or already removed. Absent is the desired state.
        return
    except OSError as exc:
        logger.warning(f"Could not remove the k-mer counting directory {attempt_dir} ({exc}).")
