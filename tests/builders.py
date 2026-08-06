"""Shared builders for VNtyper domain objects used in tests.

Purpose
-------
Before these existed, every test file hand-rolled its own DataFrame literals
with a slightly different column set, so writing a failing test for a new stage
meant reverse-engineering the column contract first. These builders make that
contract explicit and reusable.

All builders are pure: no filesystem, no network, no mocks.
"""

import copy
import json
from pathlib import Path
from typing import Any

import pandas as pd

_SCRIPTS = Path(__file__).resolve().parents[1] / "vntyper" / "scripts"

# chr1 lengths, the decisive signal for distinguishing the two builds.
_CHR1_LENGTH = {"GRCh37": 249250621, "GRCh38": 248956422}
_CHR1_NAME = {"ucsc": "chr1", "ensembl": "1", "ncbi": {"GRCh37": "NC_000001.10", "GRCh38": "NC_000001.11"}}

STAGE_COLUMNS: dict[str, tuple[str, ...]] = {
    "raw": ("Motifs", "Variant", "POS", "REF", "ALT", "Sample", "Motif_sequence"),
    "scored": (
        "Motifs",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Sample",
        "Motif_sequence",
        "Del",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "ref_len",
        "alt_len",
        "Frame_Score",
        "is_frameshift",
        "direction",
        "frameshift_amount",
        "is_valid_frameshift",
    ),
    "confidence": (
        "Motifs",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Sample",
        "Motif_sequence",
        "Del",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "ref_len",
        "alt_len",
        "Frame_Score",
        "is_frameshift",
        "direction",
        "frameshift_amount",
        "is_valid_frameshift",
        "Depth_Score",
        "Confidence",
        "depth_confidence_pass",
    ),
    "flagged": (
        "Motifs",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Sample",
        "Motif_sequence",
        "Del",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "ref_len",
        "alt_len",
        "Frame_Score",
        "is_frameshift",
        "direction",
        "frameshift_amount",
        "is_valid_frameshift",
        "Depth_Score",
        "Confidence",
        "depth_confidence_pass",
        "haplo_count",
        "alt_filter_pass",
        "motif_filter_pass",
        "Motif",
        "Flag",
    ),
    "final": (
        "Motifs",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Sample",
        "Motif_sequence",
        "Del",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "ref_len",
        "alt_len",
        "Frame_Score",
        "is_frameshift",
        "direction",
        "frameshift_amount",
        "is_valid_frameshift",
        "Depth_Score",
        "Confidence",
        "depth_confidence_pass",
        "haplo_count",
        "alt_filter_pass",
        "motif_filter_pass",
        "Motif",
        "Flag",
        "Motif_fasta",
        "POS_fasta",
    ),
}


def kestrel_row(
    *,
    pos: int = 67,
    ref: str = "G",
    alt: str = "GG",
    motif: str = "X",
    motifs: str = "X-5",
    depth_alt: int = 120,
    depth_region: int = 12000,
    variant: str = "Insertion",
    motif_sequence: str = "SEQ1",
    **extra: Any,
) -> dict[str, Any]:
    """Build one Kestrel variant row with coherent defaults.

    Defaults describe a canonical single-base insertion at position 67 with
    ample depth -- i.e. a clear positive call.

    Args:
        pos: 1-based VNTR position.
        ref: Reference allele.
        alt: Alternate allele.
        motif: Resolved motif label.
        motifs: Raw ``left-right`` motif pair as Kestrel emits it.
        depth_alt: Alternate-allele depth.
        depth_region: Variant-active-region depth.
        variant: Variant-type label.
        motif_sequence: Motif sequence string.
        **extra: Additional columns merged into the row verbatim.

    Returns:
        dict[str, Any]: The row.
    """
    row: dict[str, Any] = {
        "Motifs": motifs,
        "Motif": motif,
        "Variant": variant,
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "Sample": f"Del:{depth_alt}:{depth_region}",
        "Motif_sequence": motif_sequence,
        "Estimated_Depth_AlternateVariant": depth_alt,
        "Estimated_Depth_Variant_ActiveRegion": depth_region,
    }
    row.update(extra)
    return row


def kestrel_df(*rows: dict[str, Any]) -> pd.DataFrame:
    """Build a DataFrame from :func:`kestrel_row` dicts, preserving order.

    Args:
        *rows: Row mappings. Defaults to a single default row when empty.

    Returns:
        pd.DataFrame: Frame with ``POS`` and the depth columns as integers.
    """
    records = list(rows) or [kestrel_row()]
    frame = pd.DataFrame(records)
    for column in ("POS", "Estimated_Depth_AlternateVariant", "Estimated_Depth_Variant_ActiveRegion"):
        if column in frame.columns:
            frame[column] = frame[column].astype(int)
    return frame


def kestrel_stage_frame(stage: str, rows: int = 1, **overrides: Any) -> pd.DataFrame:
    """Build a frame carrying exactly the columns present at a pipeline stage.

    Rows are made distinguishable by walking ``POS`` forward one per row: a
    builder that returns N identical rows silently passes a dedup test that
    should fail.

    Args:
        stage: One of ``raw``, ``scored``, ``confidence``, ``flagged``, ``final``.
        rows: How many rows to produce.
        **overrides: Passed through to :func:`kestrel_row`. A ``pos`` override
            sets the position of the first row; later rows walk forward from it.

    Returns:
        pd.DataFrame: Frame whose columns are exactly ``STAGE_COLUMNS[stage]``.

    Raises:
        ValueError: If ``stage`` is not a known stage name.
    """
    if stage not in STAGE_COLUMNS:
        raise ValueError(f"Unknown stage {stage!r}; expected one of {sorted(STAGE_COLUMNS)}")

    records = []
    for index in range(rows):
        base = kestrel_row(**overrides)
        base["POS"] = base["POS"] + index
        ref_len, alt_len = len(base["REF"]), len(base["ALT"])
        delta = alt_len - ref_len
        depth_alt = base["Estimated_Depth_AlternateVariant"]
        depth_region = base["Estimated_Depth_Variant_ActiveRegion"]
        depth_score = depth_alt / depth_region if depth_region else float("nan")
        enriched = {
            **base,
            "Del": "Del",
            "ref_len": ref_len,
            "alt_len": alt_len,
            "Frame_Score": delta / 3,
            "is_frameshift": delta % 3 != 0,
            "direction": (delta > 0) - (delta < 0),
            "frameshift_amount": abs(delta) % 3,
            "is_valid_frameshift": (delta > 0 and abs(delta) % 3 == 1) or (delta < 0 and abs(delta) % 3 == 2),
            "Depth_Score": depth_score,
            "Confidence": "High_Precision*",
            "depth_confidence_pass": True,
            "haplo_count": 1,
            "alt_filter_pass": True,
            "motif_filter_pass": True,
            "Flag": "Not flagged",
            "Motif_fasta": base["Motif"],
            "POS_fasta": base["POS"],
        }
        records.append(enriched)

    frame = pd.DataFrame(records)
    return frame[list(STAGE_COLUMNS[stage])]


def kestrel_config(**dotted_overrides: Any) -> dict[str, Any]:
    """Return a deep copy of the shipped Kestrel config with dotted overrides.

    Args:
        **dotted_overrides: Keys like
            ``**{"confidence_assignment.depth_score_thresholds.low": 0.5}``.

    Returns:
        dict[str, Any]: The modified copy. The on-disk config is never touched.
    """
    config = json.loads((_SCRIPTS / "kestrel_config.json").read_text(encoding="utf-8"))
    for dotted, value in dotted_overrides.items():
        node = config
        parts = dotted.split(".")
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = value
    return copy.deepcopy(config)


def bam_header(convention: str = "ucsc", assembly: str = "GRCh37", n_contigs: int = 25) -> str:
    """Build a SAM ``@SQ`` header in the requested naming convention.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``; sets the chr1 length.
        n_contigs: How many contigs to emit.

    Returns:
        str: The header text.

    Raises:
        ValueError: If the convention or assembly is unknown.
    """
    if assembly not in _CHR1_LENGTH:
        raise ValueError(f"Unknown assembly {assembly!r}")
    if convention not in _CHR1_NAME:
        raise ValueError(f"Unknown convention {convention!r}")

    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    for number in range(1, n_contigs + 1):
        if convention == "ucsc":
            name = f"chr{number}"
        elif convention == "ensembl":
            name = str(number)
        else:
            name = f"NC_{number:06d}.{'10' if assembly == 'GRCh37' else '11'}"
        length = _CHR1_LENGTH[assembly] if number == 1 else 100_000_000 + number
        lines.append(f"@SQ\tSN:{name}\tLN:{length}")
    return "\n".join(lines) + "\n"
