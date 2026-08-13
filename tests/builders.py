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
_CONVENTIONS = ("ucsc", "ensembl", "ncbi")

# `_construct_ncbi_accession` documents its second argument as "hg19 or hg38".
_UCSC_ALIAS = {"GRCh37": "hg19", "GRCh38": "hg38"}

# A chr1 length that belongs to neither build. 249_250_621 is GRCh37 and
# 248_956_422 is GRCh38; this sits between them and matches nothing.
UNKNOWN_CHR1_LENGTH = 249_000_000

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
        "flag_filter_pass",
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
        "flag_filter_pass",
    ),
    # Nomenclature is computed after filtering, so it forms a terminal stage rather
    # than widening ``final``. Every earlier stage's columns stay a subset of every
    # later stage's, which is the invariant ``test_builders`` asserts.
    "named": (
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
        "flag_filter_pass",
        "Nomenclature",
        "Nomenclature_Tier",
        "Nomenclature_Flags",
        "Ambiguity_Interval",
        "Repeat_Form",
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
            # #174's artifact gate. True by default because the default row is
            # "Not flagged", and the gate is False only for a *declared* artifact.
            "flag_filter_pass": True,
            "Motif_fasta": base["Motif"],
            "POS_fasta": base["POS"],
        }
        records.append(enriched)

    frame = pd.DataFrame(records)
    if stage == "named":
        # Use the real annotator rather than literal values: the nomenclature module
        # is pure, so this stays a pure builder, and a builder that invented its own
        # names would drift from what the pipeline actually writes.
        from vntyper.scripts.nomenclature_annotate import annotate_kestrel_frame

        frame = annotate_kestrel_frame(frame)
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


def _validate_naming(convention: str, assembly: str) -> None:
    """Reject an unknown convention or assembly.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``.

    Raises:
        ValueError: If either is unknown.
    """
    if assembly not in _CHR1_LENGTH:
        raise ValueError(f"Unknown assembly {assembly!r}")
    if convention not in _CONVENTIONS:
        raise ValueError(f"Unknown convention {convention!r}")


def _contig_name(number: int, convention: str, assembly: str) -> str:
    """Build one contig name the way production builds it.

    Chromosome numbering follows ``chromosome_utils._build_chromosome_name``:
    1-22 are themselves, 23 is X, 24 is Y and 25 is the mitochondrion.

    NCBI accessions are delegated to production's own version table rather than
    reconstructed here. They are *chromosome-specific* -- GRCh37 chr1 is
    ``NC_000001.10`` but GRCh37 chr5 is ``NC_000005.9`` -- and an earlier version
    of this builder stamped one version on every chromosome, producing accessions
    that do not exist. That kind of wrong builder makes every consumer green for
    the wrong reason, so the table is imported, not copied.

    Args:
        number: Chromosome number, 1-25.
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``.

    Returns:
        str: The contig name.
    """
    if convention == "ncbi":
        # `_construct_ncbi_accession` takes the UCSC-style alias, not the build
        # name. See the note in tests/unit/test_builders.py.
        from vntyper.scripts.chromosome_utils import _construct_ncbi_accession

        return _construct_ncbi_accession(number, _UCSC_ALIAS[assembly])

    if number == 23:
        suffix = "X"
    elif number == 24:
        suffix = "Y"
    elif number == 25:
        suffix = "M" if convention == "ucsc" else "MT"
    else:
        suffix = str(number)

    return f"chr{suffix}" if convention == "ucsc" else suffix


def _contig_length(number: int, assembly: str) -> int:
    """Return the length for a contig.

    Only chr1 carries a real length: it is the single decisive signal for telling
    the two builds apart (``chromosome_utils.detect_assembly_from_chr1_length``).
    Every other contig gets a distinct filler length so that contigs stay
    distinguishable without implying an accuracy this builder does not have.

    Args:
        number: Chromosome number, 1-25.
        assembly: ``GRCh37`` or ``GRCh38``.

    Returns:
        int: The contig length in bases.
    """
    return _CHR1_LENGTH[assembly] if number == 1 else 100_000_000 + number


def bam_contigs(convention: str = "ucsc", assembly: str = "GRCh37", n_contigs: int = 25) -> list[dict[str, Any]]:
    """Build the parsed contig list a BAM header yields.

    This is the shape ``fastq_bam_processing.parse_contigs_from_header`` returns
    and the shape every assembly-reconciliation consumer takes: a list of
    ``{"name": str, "length": int}``. :func:`bam_header` renders exactly this
    list, so the two can never disagree.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``; sets the chr1 length.
        n_contigs: How many contigs to emit, starting at chr1.

    Returns:
        list[dict[str, Any]]: The contigs, in chromosome order.

    Raises:
        ValueError: If the convention or assembly is unknown.
    """
    _validate_naming(convention, assembly)
    return [
        {"name": _contig_name(number, convention, assembly), "length": _contig_length(number, assembly)}
        for number in range(1, n_contigs + 1)
    ]


def bam_header(
    convention: str = "ucsc",
    assembly: str = "GRCh37",
    n_contigs: int = 25,
    contigs: list[dict[str, Any]] | None = None,
) -> str:
    """Build a SAM ``@SQ`` header in the requested naming convention.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``; sets the chr1 length.
        n_contigs: How many contigs to emit.
        contigs: Render these contigs verbatim instead of building them. An entry
            whose ``length`` is None emits an ``@SQ`` line with no ``LN:`` field;
            a string ``length`` is written out as-is, which is how the malformed
            fixtures below are produced.

    Returns:
        str: The header text.

    Raises:
        ValueError: If the convention or assembly is unknown.
    """
    if contigs is None:
        contigs = bam_contigs(convention=convention, assembly=assembly, n_contigs=n_contigs)

    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    for contig in contigs:
        length = contig.get("length")
        if length is None:
            lines.append(f"@SQ\tSN:{contig['name']}")
        else:
            lines.append(f"@SQ\tSN:{contig['name']}\tLN:{length}")
    return "\n".join(lines) + "\n"


def bam_contigs_unknown_chr1_length(convention: str = "ucsc", n_contigs: int = 25) -> list[dict[str, Any]]:
    """Build contigs whose chr1 length matches neither build.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        n_contigs: How many contigs to emit.

    Returns:
        list[dict[str, Any]]: The contigs; chr1's length is
        :data:`UNKNOWN_CHR1_LENGTH`.
    """
    contigs = bam_contigs(convention=convention, assembly="GRCh37", n_contigs=n_contigs)
    if contigs:
        contigs[0]["length"] = UNKNOWN_CHR1_LENGTH
    return contigs


def bam_contigs_without_chr1(
    convention: str = "ucsc", assembly: str = "GRCh37", n_contigs: int = 25
) -> list[dict[str, Any]]:
    """Build contigs with chr1 absent entirely.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``.
        n_contigs: How many contigs to build before chr1 is dropped.

    Returns:
        list[dict[str, Any]]: The contigs, chr1 removed.

    Raises:
        ValueError: If ``n_contigs`` leaves nothing behind.
    """
    if n_contigs < 2:
        raise ValueError(f"n_contigs={n_contigs} leaves no contigs once chr1 is dropped")
    return bam_contigs(convention=convention, assembly=assembly, n_contigs=n_contigs)[1:]


def bam_contigs_mixed_conventions(assembly: str = "GRCh37") -> list[dict[str, Any]]:
    """Build contigs whose names do not agree on a single naming convention.

    chr1 stays UCSC-named and keeps its decisive length, so the *build* is still
    recoverable from this fixture even though the *convention* is not: those two
    questions are separate and must not be conflated.

    Args:
        assembly: ``GRCh37`` or ``GRCh38``.

    Returns:
        list[dict[str, Any]]: A UCSC chr1 followed by ENSEMBL- and NCBI-named
        contigs, none of which reaches a majority.
    """
    _validate_naming("ucsc", assembly)
    ucsc = bam_contigs(convention="ucsc", assembly=assembly, n_contigs=1)
    ensembl = bam_contigs(convention="ensembl", assembly=assembly, n_contigs=3)[1:]
    ncbi = bam_contigs(convention="ncbi", assembly=assembly, n_contigs=5)[3:]
    return ucsc + ensembl + ncbi


def bam_header_missing_chr1_length(convention: str = "ucsc", assembly: str = "GRCh37", n_contigs: int = 25) -> str:
    """Build a header whose chr1 ``@SQ`` line carries no ``LN:`` field.

    ``parse_contigs_from_header`` keeps only contigs with an integer length, so
    this parses to a contig list with no chr1 at all rather than to a chr1 with a
    missing length. A consumer that only ever sees the parsed list cannot tell
    the two apart.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``.
        n_contigs: How many contigs to emit.

    Returns:
        str: The header text.
    """
    contigs = bam_contigs(convention=convention, assembly=assembly, n_contigs=n_contigs)
    if contigs:
        contigs[0]["length"] = None
    return bam_header(contigs=contigs)


def bam_header_malformed_chr1_length(
    convention: str = "ucsc",
    assembly: str = "GRCh37",
    n_contigs: int = 25,
    length: str = "not-a-number",
) -> str:
    """Build a header whose chr1 ``LN:`` field is not an integer.

    ``parse_contigs_from_header`` catches the ``ValueError``, sets the length to
    None and then drops the contig -- see
    :func:`bam_header_missing_chr1_length` for why that matters.

    Args:
        convention: One of ``ucsc``, ``ensembl``, ``ncbi``.
        assembly: ``GRCh37`` or ``GRCh38``.
        n_contigs: How many contigs to emit.
        length: The malformed value to write into chr1's ``LN:`` field.

    Returns:
        str: The header text.
    """
    contigs = bam_contigs(convention=convention, assembly=assembly, n_contigs=n_contigs)
    if contigs:
        contigs[0]["length"] = length
    return bam_header(contigs=contigs)
