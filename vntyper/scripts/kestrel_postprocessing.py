"""Pure production Kestrel postprocessing shared by live calls and replay."""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass

import pandas as pd

from vntyper.scripts.flagging import (
    KESTREL_FLAG_COLUMNS,
    CompiledFlagRules,
    add_artifact_gate,
    add_flags,
    compile_flag_rules,
    validate_duplicate_flagging_config,
)
from vntyper.scripts.identity_candidate_persistence import (
    IDENTITY_CAPTURE_COLUMNS,
    IDENTITY_SELECTION_COLUMNS,
    candidate_capture_cells,
    complete_candidate_projection_cells,
    selected_candidate_cells,
)
from vntyper.scripts.identity_candidates import (
    IdentityTranslationComponent,
    capture_kestrel_observations,
    overlay_legacy_projection,
    with_candidate_evidence,
)
from vntyper.scripts.kestrel_decision_config import KestrelSelection
from vntyper.scripts.molecular_identity_presentation import (
    IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS,
    identity_translation_diagnostic_cells,
)
from vntyper.scripts.motif_processing import motif_correction_and_annotation
from vntyper.scripts.scoring import extract_frameshifts, split_depth_and_calculate_frame_score, split_frame_score
from vntyper.scripts.variant_parsing import filter_by_alt_values_and_finalize

logger = logging.getLogger(__name__)

FILTER_COLUMNS: tuple[str, ...] = (
    "is_frameshift",
    "is_valid_frameshift",
    "depth_confidence_pass",
    "alt_filter_pass",
    "motif_filter_pass",
    "flag_filter_pass",
)

AddHaplotypeCount = Callable[[pd.DataFrame], pd.DataFrame]
SelectSingleVariant = Callable[[pd.DataFrame, KestrelSelection], pd.DataFrame]
AnnotateMotifs = Callable[[pd.DataFrame, pd.DataFrame, dict[str, object]], pd.DataFrame]
PrefilterObserver = Callable[[pd.DataFrame], None]


@dataclass(frozen=True)
class KestrelEvaluation:
    """Complete production prefilter evidence and the selected result."""

    prefilter: pd.DataFrame
    selected: pd.DataFrame
    reached_final_filter: bool


def _early(frame: pd.DataFrame) -> KestrelEvaluation:
    return KestrelEvaluation(frame, frame, False)


def filter_and_select_kestrel_candidates(
    frame: pd.DataFrame,
    *,
    selection: KestrelSelection,
    select_single_best_variant_fn: SelectSingleVariant,
) -> pd.DataFrame:
    """Apply all production gates and the configured deterministic selection.

    Args:
        frame: Complete annotated candidate population.
        selection: Validated production selection policy.
        select_single_best_variant_fn: Production selection entry point.

    Returns:
        Passing candidates reduced to at most one selected row.

    Raises:
        ValueError: If a nonempty frame is missing any mandatory gate.
    """
    if frame.empty:
        return frame
    final_mask = pd.Series(True, index=frame.index)
    for column in selection.final_filter_columns:
        if column not in frame.columns:
            message = (
                f"Required filter column '{column}' is missing from a non-empty Kestrel result frame. "
                "An upstream stage stopped emitting it, so its safety gate would silently become a "
                "permit. Aborting rather than reporting unfiltered variants. See issue #185."
            )
            logger.error(message)
            raise ValueError(message)
        final_mask &= frame[column]
    passing = frame[final_mask].copy()
    if len(passing) > 1:
        return select_single_best_variant_fn(passing, selection)
    return passing


def evaluate_kestrel_candidates(
    combined_df: pd.DataFrame,
    merged_motifs: pd.DataFrame,
    kestrel_config: dict[str, object],
    *,
    selection: KestrelSelection,
    add_haplo_count_fn: AddHaplotypeCount,
    select_single_best_variant_fn: SelectSingleVariant,
    compiled_flag_rules: CompiledFlagRules | None = None,
    identity_component: IdentityTranslationComponent | None = None,
    retain_complete_identity_candidates: bool = False,
    motif_annotation_fn: AnnotateMotifs = motif_correction_and_annotation,
    prefilter_observer: PrefilterObserver | None = None,
) -> KestrelEvaluation:
    """Run the production decision chain without filesystem effects.

    Args:
        combined_df: Ordered complete post-VCF candidate population.
        merged_motifs: Parsed motif annotation table.
        kestrel_config: Complete Kestrel decision component.
        selection: Validated production selection rules.
        add_haplo_count_fn: Production haplotype multiplicity helper.
        select_single_best_variant_fn: Production selection entry point.
        compiled_flag_rules: Optionally precompiled production flag rules.
        identity_component: Optional frozen identity translation component.
        retain_complete_identity_candidates: Whether to retain every eligible identity projection.
        motif_annotation_fn: Production motif annotation entry point.
        prefilter_observer: Optional stage-boundary evidence writer invoked before gating.

    Returns:
        Complete annotated prefilter evidence and selected result.
    """
    if compiled_flag_rules is None:
        flagging = kestrel_config.get("flagging_rules", {})
        if not isinstance(flagging, dict):
            raise ValueError("Kestrel flagging_rules must be a mapping")
        compiled_flag_rules = compile_flag_rules(flagging, KESTREL_FLAG_COLUMNS)
    duplicates_config = kestrel_config.get("duplicate_flagging", {})
    if not isinstance(duplicates_config, dict):
        raise ValueError("Kestrel duplicate_flagging must be a mapping")
    validate_duplicate_flagging_config(duplicates_config, compiled_flag_rules)

    if combined_df.empty:
        return _early(combined_df)
    frame = split_depth_and_calculate_frame_score(combined_df, modulus=selection.modulus)
    if frame.empty:
        return _early(frame)
    frame = split_frame_score(frame, modulus=selection.modulus)
    if frame.empty:
        return _early(frame)
    frame = extract_frameshifts(
        frame,
        frameshift={
            "insertion_remainder": selection.insertion_remainder,
            "deletion_remainder": selection.deletion_remainder,
        },
    )
    if frame.empty:
        return _early(frame)

    from vntyper.scripts.confidence_assignment import calculate_depth_score_and_assign_confidence

    frame = calculate_depth_score_and_assign_confidence(frame, kestrel_config)
    if frame.empty:
        return _early(frame)

    identity_candidates = None
    if identity_component is not None:
        capture_records = frame.to_dict("records")
        for record in capture_records:
            record["Motif_sequence"] = str(record["Motif_sequence"])
        identity_candidates = capture_kestrel_observations(capture_records, identity_component)
        capture_rows = [candidate_capture_cells(candidate) for candidate in identity_candidates.candidates]
        frame = frame.copy()
        for column in IDENTITY_CAPTURE_COLUMNS:
            frame[column] = [cells[column] for cells in capture_rows]
        diagnostics = [
            identity_translation_diagnostic_cells(candidate.observation.translation)
            for candidate in identity_candidates.candidates
        ]
        for column in IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS:
            frame[column] = [cells[column] for cells in diagnostics]

    frame = add_haplo_count_fn(frame)
    frame = filter_by_alt_values_and_finalize(frame, kestrel_config)
    if frame.empty:
        return _early(frame)
    frame = motif_annotation_fn(frame, merged_motifs, kestrel_config)
    if frame.empty:
        return _early(frame)
    if compiled_flag_rules.rules or duplicates_config.get("enabled", False):
        frame = add_flags(frame, compiled_flag_rules, duplicates_config=duplicates_config)
    artifacts = kestrel_config.get("artifact_flags", [])
    if not isinstance(artifacts, list):
        raise ValueError("Kestrel artifact_flags must be a list")
    frame = add_artifact_gate(frame, artifacts)

    evidenced_candidates = None
    passing_identity_ordinals: tuple[int, ...] = ()
    if identity_candidates is not None:
        evidenced_candidates = with_candidate_evidence(identity_candidates, frame.to_dict("records"))
        passing_mask = frame[list(selection.final_filter_columns)].all(axis=1)
        passing_identity_ordinals = tuple(
            int(serialized) for serialized in frame.loc[passing_mask, IDENTITY_CAPTURE_COLUMNS[5]]
        )
        if retain_complete_identity_candidates:
            projections = complete_candidate_projection_cells(evidenced_candidates, passing_identity_ordinals)
            for column in IDENTITY_SELECTION_COLUMNS:
                frame[column] = ""
            for row_index in frame.index[passing_mask]:
                ordinal = int(frame.loc[row_index, IDENTITY_CAPTURE_COLUMNS[5]])
                for column, value in projections[ordinal].items():
                    frame.loc[row_index, column] = value

    prefilter = frame.copy(deep=True)
    if prefilter_observer is not None:
        prefilter_observer(prefilter)
    selected = filter_and_select_kestrel_candidates(
        frame,
        selection=selection,
        select_single_best_variant_fn=select_single_best_variant_fn,
    )
    if evidenced_candidates is not None and not selected.empty:
        selected = selected.drop(columns=list(IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS))
        selected_ordinal = int(selected.iloc[0][IDENTITY_CAPTURE_COLUMNS[5]])
        selected_candidates = overlay_legacy_projection(
            evidenced_candidates,
            passing_identity_ordinals,
            selected_ordinal,
        )
        for column, value in selected_candidate_cells(selected_candidates).items():
            selected[column] = value
    return KestrelEvaluation(prefilter, selected, True)
