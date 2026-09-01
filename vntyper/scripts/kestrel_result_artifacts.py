"""Focused publication of Kestrel no-call result and identity artifacts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from vntyper.scripts.identity_candidate_persistence import IDENTITY_CAPTURE_COLUMNS
from vntyper.scripts.molecular_identity_presentation import IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS
from vntyper.scripts.nomenclature_bam_replay import BamReplayArtifact, write_bam_replay_artifact


def write_empty_kestrel_artifacts(
    output_dir: str | Path,
    header: list[str],
    *,
    note: str | None = None,
    preserve_pre_result: bool = False,
) -> Path:
    """Publish a Negative result plus complete canonical identity sidecars.

    Args:
        output_dir: Existing Kestrel output directory.
        header: Comment lines written before the result table.
        note: Optional additional ``##`` comment.
        preserve_pre_result: Retain the current run's scored pre-result. Early
            no-call branches replace any stale file with a canonical header-only file.

    Returns:
        Path to the published final Kestrel result.

    Raises:
        OSError: If any completed-run artifact cannot be written.
        ValueError: If the output directory or replay artifact is invalid.
    """
    root = Path(output_dir)
    result = root / "kestrel_result.tsv"
    empty_result_data = {
        "Motif": ["None"],
        "Variant": ["None"],
        "POS": ["None"],
        "REF": ["None"],
        "ALT": ["None"],
        "Motif_sequence": ["None"],
        "Estimated_Depth_AlternateVariant": ["None"],
        "Estimated_Depth_Variant_ActiveRegion": ["None"],
        "Depth_Score": ["None"],
        "Confidence": ["Negative"],
    }
    banner = [*header, *([f"## {note}"] if note else [])]
    with result.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(banner) + "\n")
        pd.DataFrame(empty_result_data).to_csv(handle, sep="\t", index=False)
    pre_result = root / "kestrel_pre_result.tsv"
    if not preserve_pre_result:
        pre_result.write_text(
            "\t".join((*IDENTITY_CAPTURE_COLUMNS, *IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS)) + "\n",
            encoding="utf-8",
        )
    write_bam_replay_artifact(root, BamReplayArtifact(()))
    return result
