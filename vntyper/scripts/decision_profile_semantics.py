"""Semantic validation for fully projected decision-profile components."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from vntyper.scripts.comparator_rules import validate_rule
from vntyper.scripts.cross_match import CROSS_MATCH_COLUMNS
from vntyper.scripts.flagging import (
    ADVNTR_FLAG_COLUMNS,
    KESTREL_FLAG_COLUMNS,
    compile_flag_rules,
    validate_duplicate_flagging_config,
)
from vntyper.scripts.kestrel_decision_config import project_kestrel_selection
from vntyper.scripts.nomenclature_decision_config import decision_config_from_component

_CONFIDENCE_LABELS = {"Negative", "Low_Precision", "High_Precision", "High_Precision*"}
_SELECTION_COLUMNS = {"confidence_priority", "is_unflagged", "Depth_Score", "haplo_count", "POS"}


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be an object")
    return value


def _string_sequence(value: object, label: str) -> tuple[str, ...]:
    if not isinstance(value, Sequence) or isinstance(value, str):
        raise ValueError(f"{label} must be an array of strings")
    if any(not isinstance(item, str) or not item for item in value):
        raise ValueError(f"{label} must be an array of non-empty strings")
    return tuple(value)


def _validate_kestrel(component: Mapping[str, object]) -> None:
    rules = compile_flag_rules(component.get("flagging_rules"), KESTREL_FLAG_COLUMNS)
    duplicate = _mapping(component.get("duplicate_flagging"), "kestrel duplicate_flagging")
    validate_duplicate_flagging_config(duplicate, rules)
    if not isinstance(duplicate.get("enabled"), bool):
        raise ValueError("kestrel duplicate_flagging enabled must be boolean")
    group_by = _string_sequence(duplicate.get("group_by"), "kestrel duplicate_flagging group_by")
    if not group_by or any(column not in KESTREL_FLAG_COLUMNS for column in group_by):
        raise ValueError("kestrel duplicate_flagging group_by must contain known Kestrel columns")
    sort_by = duplicate.get("sort_by")
    if not isinstance(sort_by, Sequence) or isinstance(sort_by, str) or not sort_by:
        raise ValueError("kestrel duplicate_flagging sort_by must be a non-empty array")
    for index, raw_field in enumerate(sort_by):
        sort_field = _mapping(raw_field, f"kestrel duplicate_flagging sort_by[{index}]")
        if set(sort_field) != {"ascending", "column"}:
            raise ValueError("kestrel duplicate_flagging sort_by entries require ascending and column")
        if sort_field.get("column") not in KESTREL_FLAG_COLUMNS or not isinstance(sort_field.get("ascending"), bool):
            raise ValueError("kestrel duplicate_flagging sort_by entries require a known column and boolean ascending")

    artifact_flags = _string_sequence(component.get("artifact_flags"), "kestrel artifact_flags")
    rule_names = {name for name, _rule in rules.rules}
    if not artifact_flags or any(flag not in rule_names for flag in artifact_flags):
        raise ValueError("kestrel artifact_flags must name configured flagging rules")

    alt = _mapping(component.get("alt_filtering"), "kestrel alt_filtering")
    if not isinstance(alt.get("gg_alt_value"), str) or not alt["gg_alt_value"]:
        raise ValueError("kestrel alt_filtering gg_alt_value must be a non-empty string")
    _string_sequence(alt.get("exclude_alts"), "kestrel alt_filtering exclude_alts")

    motif = _mapping(component.get("motif_filtering"), "kestrel motif_filtering")
    threshold = motif.get("position_threshold")
    if isinstance(threshold, bool) or not isinstance(threshold, int) or threshold <= 0:
        raise ValueError("kestrel motif_filtering position_threshold must be a positive integer")
    if not isinstance(motif.get("use_uniform_filtering"), bool):
        raise ValueError("kestrel motif_filtering use_uniform_filtering must be boolean")
    for motif_field in (
        "exclude_motifs_right",
        "motifs_for_alt_gg",
        "exclude_alts_combined",
        "exclude_motifs_combined",
    ):
        _string_sequence(motif.get(motif_field), f"kestrel motif_filtering {motif_field}")
    if not isinstance(motif.get("alt_for_motif_right_gg"), str) or not motif["alt_for_motif_right_gg"]:
        raise ValueError("kestrel motif_filtering alt_for_motif_right_gg must be a non-empty string")

    selection_raw = _mapping(component.get("selection"), "kestrel selection")
    selection = project_kestrel_selection(selection_raw)
    if not selection.unflagged_value:
        raise ValueError("kestrel selection unflagged_value must be a non-empty string")
    if set(selection.confidence_priority) != _CONFIDENCE_LABELS:
        raise ValueError("kestrel selection confidence_priority labels differ from the closed result vocabulary")
    ranks = tuple(selection.confidence_priority.values())
    if any(isinstance(rank, bool) or not isinstance(rank, int) or rank < 0 for rank in ranks):
        raise ValueError("kestrel selection confidence_priority values must be non-negative integers")
    if len(set(ranks)) != len(ranks):
        raise ValueError("kestrel selection confidence_priority values must be distinct")
    sort_columns = tuple(field.column for field in selection.sort_order)
    if set(sort_columns) != _SELECTION_COLUMNS or len(sort_columns) != len(_SELECTION_COLUMNS):
        raise ValueError("kestrel selection sort_order must contain each closed selection column exactly once")


def _validate_advntr(component: Mapping[str, object]) -> None:
    settings = _mapping(component.get("settings"), "adVNTR settings")
    for field in ("max_frameshift", "frameshift_multiplier", "vid"):
        value = settings.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"adVNTR settings {field} must be a positive integer")
    if settings.get("output_format") not in {"vcf", "tsv"}:
        raise ValueError("adVNTR settings output_format must be 'vcf' or 'tsv'")

    flagging = _mapping(component.get("flagging_rules"), "adVNTR flagging_rules")
    compile_flag_rules(flagging, ADVNTR_FLAG_COLUMNS)
    evidence = _mapping(component.get("artifact_evidence"), "adVNTR artifact_evidence")
    active_states = _string_sequence(evidence.get("active_states"), "adVNTR artifact_evidence active_states")
    active_statuses = _string_sequence(evidence.get("active_statuses"), "adVNTR artifact_evidence active_statuses")
    if len(active_states) != len(active_statuses):
        raise ValueError("adVNTR artifact evidence active state and status arrays must have equal length")
    polymorphic = _mapping(flagging.get("Polymorphic_Call"), "adVNTR Polymorphic_Call")
    all_nodes = polymorphic.get("all")
    if not isinstance(all_nodes, Sequence) or isinstance(all_nodes, str) or len(all_nodes) != 1:
        raise ValueError("adVNTR Polymorphic_Call must contain exactly one governed predicate")
    predicate = _mapping(all_nodes[0], "adVNTR Polymorphic_Call predicate")
    right = _mapping(predicate.get("right"), "adVNTR Polymorphic_Call right operand")
    governed_states = _string_sequence(right.get("literal"), "adVNTR Polymorphic_Call right literal")
    if governed_states != active_states:
        raise ValueError("adVNTR Polymorphic_Call states must equal governed artifact evidence active_states")


def _validate_cross_match(component: Mapping[str, object]) -> None:
    if component.get("required_advntr_evidence_disposition") != "admissible":
        raise ValueError("cross_match required_advntr_evidence_disposition must remain admissible")
    validate_rule(
        component.get("match_rule"),
        allowed_columns=CROSS_MATCH_COLUMNS,
        context="cross_match.match_rule",
    )


def validate_component_semantics(components: Mapping[str, object]) -> None:
    """Reject profiles whose complete values are unsafe for their consumers.

    Args:
        components: Reconstructed complete component projection.

    Raises:
        ValueError: If a consumer would reject, truncate, or reinterpret a value.
    """
    kestrel = _mapping(components.get("kestrel"), "kestrel component")
    advntr = _mapping(components.get("advntr"), "adVNTR component")
    nomenclature = _mapping(components.get("nomenclature"), "nomenclature component")
    cross_match = _mapping(components.get("cross_match"), "cross_match component")
    shark = _mapping(components.get("shark"), "SHARK component")
    if shark:
        raise ValueError("SHARK decision component must remain empty")
    _validate_kestrel(kestrel)
    _validate_advntr(advntr)
    decision_config_from_component(nomenclature)
    _validate_cross_match(cross_match)
