"""Attach MUC1 nomenclature columns to caller result frames.

Kept separate from :mod:`vntyper.scripts.nomenclature` so the translator stays pure
and pandas-free: this is the only place the two worlds meet.

The five columns are added in one pass for both callers, so the result TSVs, the
summary JSON built from those frames, and the HTML report all inherit them from a
single edit.

Nullable by design: ``Ambiguity_Interval`` and ``Repeat_Form`` are written only when
they mean something -- an interval wider than 1 bp, a variant inside a detectable
tract -- and are the empty string otherwise, never a placeholder like ``N/A``.

Research use only.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from contextlib import suppress
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from vntyper.scripts.molecular_identity import IdentityTranslation
from vntyper.scripts.molecular_identity_presentation import (
    IDENTITY_COLUMNS,
    persisted_identity_result_rows,
)
from vntyper.scripts.nomenclature import (
    Nomenclature,
    from_advntr,
    from_kestrel,
    reconcile,
    render,
)
from vntyper.scripts.nomenclature_bam import is_candidate, refine
from vntyper.scripts.nomenclature_bam_adapter import (
    group_identity_rows as _group_identity_rows,
)
from vntyper.scripts.nomenclature_bam_adapter import (
    haplotype_calls as _haplotype_calls,
)
from vntyper.scripts.nomenclature_bam_adapter import (
    observe_identity_groups as _observe_identity_groups,
)
from vntyper.scripts.nomenclature_bam_adapter import (
    open_rescuer as _open_rescuer,
)
from vntyper.scripts.nomenclature_bam_adapter import (
    row_haplotype_call as _row_haplotype_call,
)
from vntyper.scripts.nomenclature_bam_replay import (
    BamReplayArtifact,
    BamReplayLocus,
    invalidate_bam_replay_artifact,
    read_bam_replay_artifact,
    write_bam_replay_artifact,
)
from vntyper.scripts.nomenclature_decision_config import NomenclatureDecisionConfig, decision_config_from_component
from vntyper.scripts.nomenclature_dominance_runtime import (
    DOMINANCE_ABSTENTION_COLUMN,
    dominance_abstention_note,
    reconcile_retained_dominance,
    retain_bam_replay,
)
from vntyper.scripts.nomenclature_dominance_runtime import (
    production_identity_observations as _production_identity_observations,
)
from vntyper.scripts.nomenclature_dominance_runtime import (
    rows_carry_identity_metadata as _rows_carry_identity_metadata,
)
from vntyper.scripts.nomenclature_frame_presentation import (
    NOMENCLATURE_COLUMNS,
    append_decision_explanation,
)
from vntyper.scripts.nomenclature_frame_presentation import (
    annotate_advntr_frame as annotate_advntr_frame,
)
from vntyper.scripts.nomenclature_frame_presentation import (
    is_negative_result_row as _is_negative,
)
from vntyper.scripts.nomenclature_frame_presentation import (
    nomenclature_result_cells as _cells,
)
from vntyper.scripts.nomenclature_identity_projection import (
    IdentityAwareNomenclatureResult,
    reconcile_identity_aware_nomenclature,
)
from vntyper.scripts.reconciled_identity_persistence import (
    RECONCILED_IDENTITY_COLUMN,
    encode_reconciled_identity,
)
from vntyper.scripts.run_configuration import resolve_compatibility_component

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from vntyper.modules.advntr.artifact_evidence import ArtifactEvidence


def _summarise(calls: list[Nomenclature]) -> str:
    """One caller's verdict as display text.

    Distinct renderings joined, so a caller reporting several events says so rather
    than being reduced to its first one. Empty means the caller contributed nothing,
    which is different from contributing a call it could not name -- that renders as
    a statement of what is known.

    Args:
        calls: Every call from one caller.

    Returns:
        str: The joined renderings, or ``""`` when there are none.
    """
    seen: list[str] = []
    for call in calls:
        text = render(call)
        if text and text not in seen:
            seen.append(text)
    return ";".join(seen)


def annotate_kestrel_frame(
    frame: pd.DataFrame,
    output_dir: str | Path | None = None,
    *,
    identity_component: Any = None,
    resolved_component: Mapping[str, object] | None = None,
    custom_context_active: bool = False,
) -> pd.DataFrame:
    """Add the nomenclature columns to a Kestrel result frame.

    Args:
        frame: The final Kestrel frame, as written to ``kestrel_result.tsv``.
        output_dir: The Kestrel output directory, when resolved haplotype records
            are available. Given one, rows the VCF cannot resolve are refined
            against ``output.bam`` -- the rescue path that recovers a delins.
            Omitted, the VCF result stands and no BAM is opened.
        identity_component: Current-run complete-context translation component.
        resolved_component: Immutable nomenclature component for this run.
        custom_context_active: Whether an explicit custom profile owns this run.

    Returns:
        pd.DataFrame: A copy with the five columns appended. The negative-run frame
        is returned unchanged: it has no variant, and padding it with five empty
        columns would widen a schema that deliberately differs.
    """
    decision_component = resolve_compatibility_component(
        "nomenclature",
        resolved_component,
        custom_context_active=custom_context_active,
    )
    decision_config = decision_config_from_component(decision_component)
    if output_dir is not None:
        invalidate_bam_replay_artifact(output_dir)
    if frame.empty or "Motifs" not in frame.columns:
        return frame

    annotated = frame.copy()
    rows = [row for _, row in annotated.iterrows()]
    result_identity_cells: list[dict[str, str | int]] | None = None
    if _rows_carry_identity_metadata(rows):
        result_identity_cells = persisted_identity_result_rows(rows)
    identity_aware = output_dir is not None and _rows_carry_identity_metadata(rows)
    groups = None
    if identity_aware:
        assert output_dir is not None
        groups = _group_identity_rows(rows)
        if identity_component is None:
            from vntyper.scripts.identity_candidates import translation_component_from_config

            identity_component = translation_component_from_config(decision_component)

    merged_calls: list[Nomenclature | None] = []
    eligible_indices: set[int] = set()
    for row_index, row in enumerate(rows):
        if _is_negative(row):
            merged_calls.append(None)
            continue
        call = _kestrel_row_call(row, decision_config)
        merged = reconcile(call, decision_config=decision_config) if call is not None else None
        merged_calls.append(merged)
        if merged is not None and is_candidate(merged):
            eligible_indices.add(row_index)

    rescuer = _open_rescuer(output_dir, decision_config)
    replay_artifact: BamReplayArtifact | None = None

    try:
        if groups is not None:
            projection = _observe_identity_groups(
                groups,
                len(rows),
                frozenset(eligible_indices),
                rescuer,
                identity_component,
                decision_config,
            )
            bam_calls = list(projection.calls)
            replay_artifact = projection.artifact
        else:
            bam_calls = [None] * len(rows)
            if rescuer is not None:
                for row_index in eligible_indices:
                    bam_calls[row_index], _ = _row_haplotype_call(rows[row_index], rescuer, decision_config)

        cells: list[dict[str, Any]] = []
        for row_index, merged in enumerate(merged_calls):
            if merged is None:
                cells.append(dict.fromkeys(NOMENCLATURE_COLUMNS, ""))
                continue
            bam_call = bam_calls[row_index]
            merged = refine(merged, bam_call)
            # adVNTR is optional and has not run at this point, so its column is
            # empty; the cross-caller stage fills it in when it does run.
            cells.append(_cells(merged, kestrel=render(merged), decision_config=decision_config))
    finally:
        if rescuer is not None:
            rescuer.close()

    for column in NOMENCLATURE_COLUMNS:
        annotated[column] = [cell[column] for cell in cells]
    if result_identity_cells is not None:
        for column in IDENTITY_COLUMNS:
            annotated[column] = [cell[column] for cell in result_identity_cells]
    if replay_artifact is not None:
        write_bam_replay_artifact(output_dir, replay_artifact)  # type: ignore[arg-type]
    return annotated


def _kestrel_row_call(
    row: pd.Series,
    decision_config: NomenclatureDecisionConfig | None = None,
) -> Nomenclature | None:
    """Translate one Kestrel VCF row, or ``None`` when it cannot be translated.

    Args:
        row: One Kestrel result row.
        decision_config: Explicit resolved nomenclature values for a run.

    Returns:
        Nomenclature | None: The named record.
    """
    # A missing cell in a frame read with `dtype=str` is NaN, and `str(NaN)` is the
    # text "nan" -- which every emptiness check accepts as an allele. Reject the row
    # here, where the frame is still available to ask.
    if any(pd.isna(row.get(field)) for field in ("Motifs", "POS", "REF", "ALT")):
        logger.debug("Kestrel row missing a required field; not translatable: %s", dict(row))
        return None
    try:
        return from_kestrel(
            str(row["Motifs"]),
            int(row["POS"]),
            str(row["REF"]),
            str(row["ALT"]),
            decision_config,
        )
    except (TypeError, ValueError):
        logger.debug("Kestrel row not translatable to nomenclature: %s", dict(row))
        return None


def reconcile_caller_outputs(
    kestrel_tsv: str | Path,
    advntr_tsv: str | Path | None,
    kestrel_dir: str | Path | None = None,
    *,
    artifact_evidence: ArtifactEvidence | None = None,
    resolved_component: Mapping[str, object] | None = None,
    dominance_component: Mapping[str, object] | None = None,
    custom_context_active: bool = False,
) -> bool:
    """Reconcile the two callers' written results and rewrite both in place.

    Each caller names its own rows as it writes them, which is what puts the columns
    on every surface. But this is the only place that can see *both*, so it is the
    only place that can do three things no caller stage can:

    1. Promote to tier A, which requires two independent sources agreeing.
    2. Let two independent sources outvote a third -- the reason ``insG`` families
       are recovered rather than lost to a Kestrel placement one base 3' of truth.
    3. Record what each caller said, so a disagreement stays legible instead of
       being collapsed into one verdict.

    adVNTR is optional. Without it, the packaged neutral path leaves Kestrel bytes
    untouched; an explicitly enabled dominance component still evaluates this final
    whole-locus seam from retained Kestrel evidence.

    Args:
        kestrel_tsv: Path to ``kestrel_result.tsv``.
        advntr_tsv: Path to ``output_adVNTR_result.tsv``, or None when adVNTR
            was not selected and an enabled dominance profile owns the final seam.
        kestrel_dir: Directory holding ``output.bam``. Defaults to the directory the
            Kestrel TSV sits in, which is where the pipeline writes it.
        artifact_evidence: Verified governed State evidence. The packaged artifact is
            used by this standalone compatibility boundary when omitted.
        resolved_component: Immutable nomenclature component for this run.
        dominance_component: Immutable whole-locus dominance component for this run.
        custom_context_active: Whether an explicit custom profile owns this run.

    Returns:
        bool: True when both files were read and rewritten.
    """
    decision_component = resolve_compatibility_component(
        "nomenclature",
        resolved_component,
        custom_context_active=custom_context_active,
    )
    decision_config = decision_config_from_component(decision_component)
    resolved_dominance = resolve_compatibility_component(
        "dominance",
        dominance_component,
        custom_context_active=custom_context_active,
    )

    from vntyper.modules.advntr.artifact_evidence import (
        ASSERTION,
        EVIDENCE_DISPOSITION_COLUMN,
        evidence_disposition_for_state,
        load_packaged_artifact_evidence,
    )

    kestrel_path = Path(kestrel_tsv)
    advntr_path = Path(advntr_tsv) if advntr_tsv is not None else None
    if not kestrel_path.is_file() or (advntr_path is not None and not advntr_path.is_file()):
        logger.debug("Cross-caller reconciliation skipped; one of the result files is absent.")
        return False
    if advntr_path is None and resolved_dominance.get("enabled") is False:
        return False

    try:
        header = [line for line in kestrel_path.read_text().splitlines() if line.startswith("##")]
        kestrel = pd.read_csv(kestrel_path, sep="\t", comment="#", dtype=str)
        advntr = pd.read_csv(advntr_path, sep="\t", dtype=str) if advntr_path is not None else pd.DataFrame()
    except (OSError, ValueError, pd.errors.ParserError) as error:
        logger.warning("Cross-caller reconciliation skipped; could not read a result file: %s", error)
        return False

    if kestrel.empty or (advntr_path is not None and advntr.empty) or "Motifs" not in kestrel.columns:
        return False

    resolved_artifact_evidence = artifact_evidence or load_packaged_artifact_evidence()

    supports: dict[str, int | None] = {}
    kestrel_keep = [not _is_negative(row) for _, row in kestrel.iterrows()]
    kestrel_rows = [row for (_, row), keep in zip(kestrel.iterrows(), kestrel_keep, strict=True) if keep]
    identity_aware = _rows_carry_identity_metadata(kestrel_rows)
    existing_replay: BamReplayArtifact | None = None
    if identity_aware:
        with suppress(FileNotFoundError):
            existing_replay = read_bam_replay_artifact(kestrel_dir or kestrel_path.parent)

    # Aligned to `kestrel_rows`, `None` where a row could not be translated. Position
    # is what ties a call back to the record it came from, so nothing is compacted.
    vcf_calls: list[Nomenclature | None] = []
    for row in kestrel_rows:
        call = _kestrel_row_call(row, decision_config)
        vcf_calls.append(call)
        if call is None:
            continue
        # The thinnest row sets the depth, as everywhere else. Recording only the
        # first row's depth made the tier depend on the order rows happened to be
        # written in.
        depth = _as_int(row.get("Estimated_Depth_AlternateVariant"))
        existing = supports.get("kestrel_vcf")
        supports["kestrel_vcf"] = depth if "kestrel_vcf" not in supports else _lesser(existing, depth)

    advntr_keep = [not _is_negative(row) for _, row in advntr.iterrows()]
    positive_advntr_rows = [row for (_, row), keep in zip(advntr.iterrows(), advntr_keep, strict=True) if keep]
    advntr_dispositions = [
        evidence_disposition_for_state(str(row.get("Variant", "")), resolved_artifact_evidence)
        for row in positive_advntr_rows
    ]
    advntr_calls_by_row = _advntr_calls_by_row(advntr, advntr_keep, supports, decision_config)
    advntr_calls = [call for calls in advntr_calls_by_row for call in calls]

    named_vcf = [call for call in vcf_calls if call is not None]
    if not named_vcf and not advntr_calls:
        return False

    identity_component: Any = None
    if identity_aware:
        from vntyper.scripts.identity_candidates import translation_component_from_config

        identity_component = translation_component_from_config(decision_component)
    bam_translations: list[IdentityTranslation | None] | None = [] if identity_component is not None else None
    bam_replay_loci: list[BamReplayLocus] | None = [] if identity_component is not None else None
    bam_calls = _haplotype_calls(
        kestrel_rows,
        kestrel_dir or kestrel_path.parent,
        named_vcf,
        advntr_calls,
        supports,
        bam_translations=bam_translations,
        identity_component=identity_component,
        bam_replay_loci=bam_replay_loci,
        decision_config=decision_config,
    )
    named_bam = [call for call in bam_calls if call is not None]
    retained_replay: BamReplayArtifact | None = None
    if bam_replay_loci is not None:
        from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells

        expected_ordinals = tuple(
            sorted(parse_selected_candidate_cells(row).selected_observation_ordinal for row in kestrel_rows)
        )
        retained_replay = retain_bam_replay(existing_replay, tuple(bam_replay_loci), expected_ordinals)

    # Order matters: the Kestrel VCF is offered first so that whenever nothing
    # outvotes it, it is the call that stands. adVNTR is optional; Kestrel is not.
    ordered_calls = (*named_vcf, *named_bam, *advntr_calls)
    identity_inputs = _production_identity_observations(
        kestrel_rows,
        vcf_calls,
        bam_calls,
        positive_advntr_rows,
        advntr_calls_by_row,
        bam_translations,
        decision_config=decision_config,
        translation_component=identity_component,
        artifact_evidence=resolved_artifact_evidence,
    )
    merged: Nomenclature | None
    if identity_inputs is None:
        if resolved_dominance.get("enabled") is True:
            raise ValueError("enabled dominance requires current-run canonical identity metadata")
        merged = reconcile(*ordered_calls, supports=supports, decision_config=decision_config)
        reconciled_identity = None
        identity_aware_result = None
    else:
        identity_observations, policy = identity_inputs
        identity_aware_result = reconcile_identity_aware_nomenclature(
            *ordered_calls,
            supports=supports,
            identity_observations=identity_observations,
            identity_policy=policy,
            decision_config=decision_config,
        )
        merged = identity_aware_result.call
        reconciled_identity = identity_aware_result.selected_identity
    for bam_call in named_bam:
        # Still applied after the vote, because it carries one rule the vote cannot:
        # A delins is unrepresentable in Kestrel's VCF, so one seen in the
        # haplotype records is better evidence than the closest VCF shape.
        merged = refine(merged, bam_call)

    # Each row keeps its *own* caller's name. Broadcasting one joined string to every
    # row destroyed row identity: a file reporting two variants said both names on
    # both rows, so neither row described itself any more.
    kestrel_row_calls = _row_verdicts(vcf_calls, bam_calls, decision_config)
    kestrel_row_names = [render(call) if call is not None else "" for call in kestrel_row_calls]
    advntr_row_names = [_summarise(calls) for calls in advntr_calls_by_row]

    kestrel_summary = _summarise([call for call in kestrel_row_calls if call is not None])
    advntr_summary = _summarise(advntr_calls)

    dominance_result = None
    if identity_aware_result is not None and retained_replay is not None:
        from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells

        assert merged is not None
        row_identities = tuple(parse_selected_candidate_cells(row).translation.identity for row in kestrel_rows)
        dominance_result = reconcile_retained_dominance(
            IdentityAwareNomenclatureResult(merged, reconciled_identity),
            row_identities,
            tuple(kestrel_row_calls),
            identity_observations,
            ordered_calls,
            retained_replay,
            tuple(advntr_dispositions),
            resolved_dominance,
        )
        if dominance_result.decision.outcome != "not-applicable":
            merged = dominance_result.call
            reconciled_identity = dominance_result.selected_identity
        elif advntr_path is None:
            return False

    if merged is None:
        cells = dict.fromkeys(NOMENCLATURE_COLUMNS, "")
        cells["Nomenclature_Kestrel"] = kestrel_summary
        cells["Nomenclature_adVNTR"] = advntr_summary
        assert dominance_result is not None and dominance_result.abstention_reason is not None
        cells["Nomenclature_Note"] = dominance_abstention_note(dominance_result.abstention_reason)
    else:
        cells = _cells(
            merged,
            kestrel=kestrel_summary,
            advntr=advntr_summary,
            decision_config=decision_config,
        )
    if identity_inputs is not None:
        cells[RECONCILED_IDENTITY_COLUMN] = encode_reconciled_identity(reconciled_identity)
    if dominance_result is not None and dominance_result.decision.outcome != "not-applicable":
        cells[DOMINANCE_ABSTENTION_COLUMN] = dominance_result.abstention_reason or ""
    if any(disposition.value == "identity-insufficient" for disposition in advntr_dispositions):
        cells["Nomenclature_Note"] = append_decision_explanation(cells["Nomenclature_Note"], ASSERTION)

    surfaces: tuple[tuple[pd.DataFrame, Path, list[bool], list[str], str, list[str]], ...] = (
        (kestrel, kestrel_path, kestrel_keep, header, "Nomenclature_Kestrel", kestrel_row_names),
    )
    if advntr_path is not None:
        surfaces += ((advntr, advntr_path, advntr_keep, [], "Nomenclature_adVNTR", advntr_row_names),)
    for frame, path, keep, file_header, own_column, own_names in surfaces:
        updated = frame.copy()
        mask = pd.Series(keep, index=updated.index)
        for column, value in cells.items():
            if column not in updated.columns:
                # The per-caller columns are this step's own output. A file written
                # before they existed must gain them rather than silently drop the
                # information they carry.
                updated[column] = ""
            updated.loc[mask, column] = value
        updated.loc[mask, own_column] = own_names
        if own_column == "Nomenclature_adVNTR":
            if EVIDENCE_DISPOSITION_COLUMN not in updated.columns:
                updated[EVIDENCE_DISPOSITION_COLUMN] = ""
            updated.loc[mask, EVIDENCE_DISPOSITION_COLUMN] = [disposition.value for disposition in advntr_dispositions]
        _write_tsv(updated, path, file_header)

    if retained_replay is not None:
        write_bam_replay_artifact(kestrel_dir or kestrel_path.parent, retained_replay)

    logger.info(
        "Cross-caller nomenclature reconciled: tier %s (kestrel=%r advntr=%r).",
        merged.tier if merged is not None else "abstained",
        kestrel_summary,
        advntr_summary,
    )
    return True


def _lesser(left: int | None, right: int | None) -> int | None:
    """The smaller of two depths, where ``None`` means unknown and wins.

    Unknown is not sufficient, so an unknown depth anywhere in an agreement makes the
    agreement's depth unknown rather than letting a known one stand in for it.
    """
    if left is None or right is None:
        return None
    return min(left, right)


def _row_verdicts(
    vcf_calls: list[Nomenclature | None],
    bam_calls: list[Nomenclature | None],
    decision_config: NomenclatureDecisionConfig | None = None,
) -> list[Nomenclature | None]:
    """What Kestrel concludes per row, refined by that row's haplotype records.

    Both lists are aligned to the Kestrel rows, so the pairing is by position. An
    earlier version compacted the haplotype consensuses, which pushed row 2's
    records onto row 1 whenever row 1 had none.

    Args:
        vcf_calls: Per Kestrel row, ``None`` where untranslatable.
        bam_calls: Per Kestrel row, ``None`` where the haplotype records said nothing.
        decision_config: Explicit resolved nomenclature values for a run.

    Returns:
        list: Per Kestrel row, ``None`` where untranslatable.
    """
    verdicts: list[Nomenclature | None] = []
    for index, call in enumerate(vcf_calls):
        if call is None:
            verdicts.append(None)
            continue
        haplotype_call = bam_calls[index] if index < len(bam_calls) else None
        verdicts.append(refine(reconcile(call, decision_config=decision_config), haplotype_call))
    return verdicts


def _advntr_calls_by_row(
    advntr: pd.DataFrame,
    keep: list[bool],
    supports: dict[str, int | None],
    decision_config: NomenclatureDecisionConfig | None = None,
) -> list[list[Nomenclature]]:
    """Every adVNTR event, grouped by the row that reported it.

    Every event, not just the first: a sample reporting several simultaneous events is
    not describing one simple allele, and hiding the rest would let a wrong name be
    promoted. Grouped by row rather than flattened, so each row can still state its
    own call.

    Args:
        advntr: The adVNTR result frame.
        keep: Per row, whether it is a real call rather than a placeholder.
        supports: Per-source depths, updated in place.
        decision_config: Explicit resolved nomenclature values for a run.

    Returns:
        list[list[Nomenclature]]: One inner list per kept row, in row order.
    """
    grouped: list[list[Nomenclature]] = []
    for (_, row), wanted in zip(advntr.iterrows(), keep, strict=True):
        if not wanted:
            continue
        parsed = list(from_advntr(str(row.get("Variant", "")), decision_config))
        grouped.append(parsed)
        if not parsed:
            continue
        depth = _as_int(row.get("NumberOfSupportingReads"))
        existing = supports.get("advntr")
        supports["advntr"] = depth if "advntr" not in supports else _lesser(existing, depth)
    return grouped


def _as_int(value: Any) -> int | None:
    """Read a legacy depth cell using its established float-compatible coercion."""
    try:
        return int(float(value))
    except (TypeError, ValueError, OverflowError):
        return None


def _write_tsv(frame: pd.DataFrame, path: Path, header: list[str]) -> None:
    """Rewrite a result TSV, preserving its ``##`` header lines.

    Args:
        frame: The frame to write.
        path: Destination.
        header: ``##`` lines to re-emit above the table.
    """
    with path.open("w") as handle:
        if header:
            handle.write("\n".join(header) + "\n")
        frame.to_csv(handle, sep="\t", index=False)
