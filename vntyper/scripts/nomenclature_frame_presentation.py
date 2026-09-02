"""Present caller result rows with nomenclature and positive molecular identity.

This focused adapter is the pandas boundary for adVNTR result frames. The naming
and molecular-identity decisions remain in their pure owning modules.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from vntyper.scripts.molecular_identity_presentation import IDENTITY_COLUMNS, advntr_identity_result_rows
from vntyper.scripts.nomenclature import Nomenclature, confidence_note, from_advntr, reconcile, render

if TYPE_CHECKING:
    from vntyper.modules.advntr.artifact_evidence import ArtifactEvidence
    from vntyper.scripts.identity_candidates import IdentityTranslator
    from vntyper.scripts.nomenclature_decision_config import NomenclatureDecisionConfig

#: The columns, in this order, on every surface that carries them.
#:
#: ``Nomenclature`` is the reconciled verdict. The two per-caller columns beside it
#: are what each caller reported on its own, kept so a disagreement stays legible:
#: collapsing both files to a single verdict destroyed the evidence a reader needs
#: in order to weigh that verdict at all. They are named absolutely rather than
#: relatively ("the other caller"), so a row means the same thing in either file and
#: in a cohort table that merges them.
NOMENCLATURE_COLUMNS: tuple[str, ...] = (
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
)

#: A negative Kestrel run writes a different, 10-column schema whose first column is
#: singular ``Motif``. It carries no variant, so there is nothing to name.
_NEGATIVE_MARKERS = ("None", "Negative", "Not applicable")
_CANONICAL_KESTREL_NEGATIVE = {
    "Motif": "None",
    "Variant": "None",
    "POS": "None",
    "REF": "None",
    "ALT": "None",
    "Motif_sequence": "None",
    "Estimated_Depth_AlternateVariant": "None",
    "Estimated_Depth_Variant_ActiveRegion": "None",
    "Depth_Score": "None",
    "Confidence": "Negative",
}


def is_negative_result_row(row: pd.Series) -> bool:
    """Return whether a result row is a caller's negative placeholder.

    Args:
        row: One result row.

    Returns:
        True for the negative-run placeholder of either caller.
    """
    if str(row.get("Confidence", "")) == "Negative":
        return True
    if str(row.get("VID", "")) == "Negative":
        return True
    return str(row.get("Motif", "")) == "None" and "Motifs" not in row.index


def is_canonical_kestrel_negative_frame(frame: pd.DataFrame) -> bool:
    """Recognize the exact one-row negative Kestrel production schema.

    Args:
        frame: Parsed Kestrel result frame.

    Returns:
        Whether columns, values, and order equal the canonical negative artifact.
    """
    return (
        len(frame) == 1
        and tuple(frame.columns) == tuple(_CANONICAL_KESTREL_NEGATIVE)
        and all(
            str(frame.iloc[0][column]) == expected or (expected == "None" and pd.isna(frame.iloc[0][column]))
            for column, expected in _CANONICAL_KESTREL_NEGATIVE.items()
        )
    )


def nomenclature_result_cells(
    call: Nomenclature,
    kestrel: str = "",
    advntr: str = "",
    decision_config: NomenclatureDecisionConfig | None = None,
) -> dict[str, str]:
    """Project a call onto the nomenclature columns.

    The displayed name comes from :func:`render`, never from ``call.name``: the tier
    decides what may be shown, and reading the field directly would leak a bare
    number from a tier that is not entitled to state one.

    Args:
        call: The reconciled call.
        kestrel: What Kestrel reported, or ``""`` when it reported nothing.
        advntr: What adVNTR reported, or ``""`` when it did not run or reported
            nothing. adVNTR is an optional module, so empty is the ordinary case.
        decision_config: Explicit resolved nomenclature values for a run.

    Returns:
        The nomenclature column values.
    """
    interval = ""
    if call.ambiguity is not None:
        low, high = call.ambiguity
        interval = f"{low}_{high}"

    return {
        "Nomenclature": render(call),
        "Nomenclature_Tier": call.tier,
        "Nomenclature_Flags": ";".join(call.flags),
        "Ambiguity_Interval": interval,
        "Repeat_Form": call.repeat_form or "",
        "Nomenclature_Note": confidence_note(call, decision_config),
        "Nomenclature_Kestrel": kestrel,
        "Nomenclature_adVNTR": advntr,
    }


def annotate_advntr_frame(
    frame: pd.DataFrame,
    *,
    identity_component: IdentityTranslator | None = None,
    artifact_evidence: ArtifactEvidence | None = None,
    decision_config: NomenclatureDecisionConfig | None = None,
) -> pd.DataFrame:
    """Add nomenclature and optional identity columns to an adVNTR result frame.

    Unlike Kestrel, adVNTR may report several events for one sample. Each row is
    named on its own here; cross-row reconciliation belongs to whatever combines the
    two callers, which is the only place that can see both.

    Args:
        frame: The processed adVNTR frame, as written to ``output_adVNTR_result.tsv``.
        identity_component: Current-run complete-context translation component. When
            supplied, public identity cells are appended to positive rows.
        artifact_evidence: Verified governed State evidence. The packaged artifact is
            used by the standalone compatibility boundary when omitted.
        decision_config: Explicit resolved nomenclature values for a run.

    Returns:
        A copy with the nomenclature columns appended, empty on the negative
        placeholder row. The identity quartet is appended only when a component is
        supplied for a positive result frame.
    """
    from vntyper.modules.advntr.artifact_evidence import (
        ASSERTION,
        EVIDENCE_DISPOSITION_COLUMN,
        evidence_disposition_for_state,
        load_packaged_artifact_evidence,
    )

    if frame.empty or "Variant" not in frame.columns:
        return frame

    annotated = frame.copy()
    resolved_evidence = artifact_evidence or load_packaged_artifact_evidence()
    cells: list[dict[str, str]] = []
    dispositions: list[str] = []
    has_positive_row = False
    for _, row in annotated.iterrows():
        if is_negative_result_row(row):
            cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
            dispositions.append("")
            continue
        has_positive_row = True
        support: int | None
        try:
            support = int(float(row.get("NumberOfSupportingReads", "")))
        except (TypeError, ValueError):
            support = None
        calls = from_advntr(str(row["Variant"]), decision_config)
        if not calls:
            cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
            dispositions.append("admissible")
            continue
        disposition = evidence_disposition_for_state(str(row["Variant"]), resolved_evidence)
        merged = reconcile(*calls, support=support, decision_config=decision_config)
        # Kestrel's column is filled by the cross-caller stage, which is the only
        # place that can see it.
        row_cells = nomenclature_result_cells(merged, advntr=render(merged), decision_config=decision_config)
        if disposition.value == "identity-insufficient":
            row_cells["Nomenclature_Note"] = append_decision_explanation(
                row_cells["Nomenclature_Note"],
                ASSERTION,
            )
        cells.append(row_cells)
        dispositions.append(disposition.value)

    for column in NOMENCLATURE_COLUMNS:
        annotated[column] = [cell[column] for cell in cells]
    if has_positive_row:
        annotated[EVIDENCE_DISPOSITION_COLUMN] = dispositions
    if identity_component is not None:
        result_cells = advntr_identity_result_rows(annotated.to_dict("records"), identity_component)
        for column in IDENTITY_COLUMNS:
            annotated[column] = [cell[column] for cell in result_cells]
    return annotated


def append_decision_explanation(note: str, explanation: str) -> str:
    """Append one governed sentence without duplicating existing text."""
    if explanation in note:
        return note
    return f"{note}; {explanation}" if note else explanation
