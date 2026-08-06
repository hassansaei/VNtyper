# vntyper/scripts/motif_decisions.py

"""
motif_decisions.py

Module Purpose:
---------------
The pure decision layer of MUC1 motif filtering. Every function here takes a
DataFrame plus the relevant slice of ``kestrel_config.json``'s
``motif_filtering`` block and returns a new DataFrame; none of them reads a
file, reads the config dict, logs, or mutates its argument.

They were extracted verbatim from ``motif_processing.motif_correction_and
_annotation``, which had grown to 182 lines fusing the config read, the column
plumbing, the positional split and the filtering rules into one function that
could not be called without a full Kestrel frame. AGENTS.md rule 3: extract the
pure part, leave the I/O. ``motif_processing`` keeps the merge, the annotation
and the pass/fail bookkeeping and calls into here for each decision.

The rules, in the order ``motif_correction_and_annotation`` applies them:

1. ``split_left_right`` - split each motif ID on its dash and divide the rows
   into the two halves of the repeat unit by position.
2. ``has_gg_alternate`` - gate the legacy right-motif branch on the presence of
   the GG alternate anywhere in the right half.
3. ``apply_right_motif_exclusions`` - drop conserved motifs, where a call is
   more likely an artifact than a variant.
4. ``apply_gg_alt_rule`` - narrow to the motifs on the GG allowlist, but only
   if at least one row is on it.
5. ``apply_combined_exclusions`` - the final blocklists over ALT and motif,
   applied to both halves together.

Nothing here changes a motif decision. The extraction is behaviour-preserving
and is pinned by the end-to-end oracle in
``tests/unit/test_motif_decisions.py``.

References:
-----------
- Saei et al., iScience 26, 107171 (2023).
- GitHub issues #136 (uniform filtering), #179 (audit), #186 (GG allowlist).
"""

import logging
from collections.abc import Sequence

import pandas as pd

logger = logging.getLogger(__name__)


def split_left_right(df: pd.DataFrame, position_threshold: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split motif IDs on their dash and divide rows into left and right motifs.

    A MUC1 motif ID is two half-motif names joined by a dash, e.g. ``'X-3'``.
    Which half a variant belongs to depends on where in the repeat unit it sits:
    a variant at or above ``position_threshold`` is in the *right* half of the
    unit, so its motif is the name **before** the dash; below the threshold the
    motif is the name **after** it. ``motif_correction_and_annotation`` performs
    that rename, which is why this function returns both halves under both
    names.

    ``POS`` is coerced to int here because Kestrel VCF positions arrive as
    strings. A value that will not parse becomes ``-1``, which sorts below any
    real threshold and so lands in the left frame.

    Args:
        df (pd.DataFrame): Rows carrying at least 'Motifs' and 'POS'. Every
            'Motifs' value must contain exactly one dash; callers are expected
            to have dropped the rest (``motif_processing`` does this per row).
        position_threshold (int): Boundary between the halves, exclusive below
            and inclusive at and above.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: ``(motif_left, motif_right)``, each a
        copy carrying the added 'Motif_left' and 'Motif_right' columns and an
        integer 'POS'. Together they partition the input rows exactly once.
    """
    working_df = df.copy(deep=True)
    working_df[["Motif_left", "Motif_right"]] = working_df["Motifs"].str.split("-", expand=True)
    working_df["POS"] = pd.to_numeric(working_df["POS"], errors="coerce").fillna(-1).astype(int)

    motif_left = working_df[working_df["POS"] < position_threshold].copy()
    motif_right = working_df[working_df["POS"] >= position_threshold].copy()
    return motif_left, motif_right


def has_gg_alternate(motif_right: pd.DataFrame, alt_for_motif_right_gg: str) -> bool:
    """
    Report whether any right-motif row carries the GG alternate.

    This gates the whole legacy right-motif branch: the exclusions and the
    allowlist narrowing below only run when a GG alternate is present somewhere
    in the frame. The match is word-bounded, so a longer allele that merely
    contains 'GG' (``'CGGCG'``) does not open the branch.

    Args:
        motif_right (pd.DataFrame): Right-motif rows, carrying 'ALT'.
        alt_for_motif_right_gg (str): The alternate to look for, from
            ``motif_filtering.alt_for_motif_right_gg``. Typically ``'GG'``.

    Returns:
        bool: True if at least one row's ALT matches the alternate as a whole
        word.
    """
    return bool(motif_right["ALT"].str.contains(r"\b" + alt_for_motif_right_gg + r"\b").any())


def apply_right_motif_exclusions(motif_right: pd.DataFrame, exclude_motifs_right: Sequence[str]) -> pd.DataFrame:
    """
    Drop right-motif rows whose motif is on the exclusion list.

    Biological rationale: the listed motifs are conserved, so a variant called
    against one of them is more likely a mapping or k-mer artifact than a real
    event. An empty list excludes nothing.

    Args:
        motif_right (pd.DataFrame): Right-motif rows, carrying 'Motif'.
        exclude_motifs_right (Sequence[str]): Motif names to remove, from
            ``motif_filtering.exclude_motifs_right``.

    Returns:
        pd.DataFrame: The surviving rows. Column set and order are unchanged;
        the input is not modified.
    """
    return motif_right[~motif_right["Motif"].isin(exclude_motifs_right)]


def apply_gg_alt_rule(motif_right: pd.DataFrame, motifs_for_alt_gg: Sequence[str]) -> pd.DataFrame:
    """
    Narrow right-motif rows to the GG allowlist, but only if it matches anything.

    The ``.any()`` guard is the whole rule. Applied unconditionally the
    allowlist would be a filter, and an empty one would delete every row -
    including the canonical dupC call at POS 67, which is the failure mode #186
    is about. Guarded, an allowlist that matches nothing is a no-op, which is
    what makes the shipped config (``motifs_for_alt_gg: []``) inert.

    Callers gate this on :func:`has_gg_alternate` against the frame as it stood
    *before* the exclusions. That ordering is deliberate and is not re-checked
    here: re-testing for a GG alternate at this point would change the result
    whenever the exclusions had just removed the only GG row.

    Args:
        motif_right (pd.DataFrame): Right-motif rows, carrying 'Motif'.
        motifs_for_alt_gg (Sequence[str]): Motifs in which a GG alternate is
            allowed, from ``motif_filtering.motifs_for_alt_gg``.

    Returns:
        pd.DataFrame: The rows on the allowlist if any row is on it, otherwise
        the frame unchanged. The input is not modified.
    """
    if motif_right["Motif"].isin(motifs_for_alt_gg).any():
        return motif_right[motif_right["Motif"].isin(motifs_for_alt_gg)]
    return motif_right


def apply_combined_exclusions(
    df: pd.DataFrame,
    exclude_alts_combined: Sequence[str],
    exclude_motifs_combined: Sequence[str],
) -> pd.DataFrame:
    """
    Apply the final ALT and motif blocklists to the combined left/right frame.

    Both are blocklists, not keeplists: a row survives unless it is named. Two
    empty lists therefore remove nothing, which is the property that makes the
    gate safe to widen.

    Matching is exact, not substring - ``'CCGCCA'`` is a different allele from
    ``'CCGCC'`` and is not caught by listing the latter.

    Args:
        df (pd.DataFrame): Left and right motif rows concatenated, carrying
            'ALT' and 'Motif'.
        exclude_alts_combined (Sequence[str]): Alternates to remove, from
            ``motif_filtering.exclude_alts_combined``.
        exclude_motifs_combined (Sequence[str]): Motifs to remove, from
            ``motif_filtering.exclude_motifs_combined``.

    Returns:
        pd.DataFrame: The surviving rows; the input is not modified.
    """
    combined_df = df[~df["ALT"].isin(exclude_alts_combined)]
    combined_df = combined_df[~combined_df["Motif"].isin(exclude_motifs_combined)]
    return combined_df
