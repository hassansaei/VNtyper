"""Mutation-resistant boundary contracts for decision-profile semantics."""

from __future__ import annotations

from pathlib import Path

import pytest

from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.decision_profile_semantics import validate_component_semantics

pytestmark = pytest.mark.unit

PROJECTION_PATH = Path(__file__).parents[2] / "vntyper" / "profiles" / "decision_projection.json"


def _components() -> dict[str, object]:
    return load_strict_json_object(PROJECTION_PATH.read_bytes())


def _replace(components: dict[str, object], path: tuple[str | int, ...], value: object) -> None:
    current: object = components
    for token in path[:-1]:
        if isinstance(token, int):
            assert isinstance(current, list)
            current = current[token]
        else:
            assert isinstance(current, dict)
            current = current[token]
    final = path[-1]
    if isinstance(final, int):
        assert isinstance(current, list)
        current[final] = value
    else:
        assert isinstance(current, dict)
        current[final] = value


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (("kestrel", "alt_filtering", "exclude_alts"), {}, "array of strings"),
        (("kestrel", "alt_filtering", "exclude_alts"), [1], "non-empty strings"),
        (("kestrel", "duplicate_flagging", "group_by"), [], "group_by"),
        (("kestrel", "duplicate_flagging", "sort_by"), {}, "sort_by"),
        (("kestrel", "duplicate_flagging", "sort_by"), 7, "must be a non-empty array"),
        (("kestrel", "duplicate_flagging", "sort_by"), "P", "must be a non-empty array"),
        (("kestrel", "duplicate_flagging", "sort_by"), [], "sort_by"),
        (("kestrel", "alt_filtering", "gg_alt_value"), 7, "gg_alt_value"),
        (("kestrel", "alt_filtering", "gg_alt_value"), "", "gg_alt_value"),
        (("kestrel", "motif_filtering", "position_threshold"), True, "position_threshold"),
        (("kestrel", "motif_filtering", "position_threshold"), 1.0, "position_threshold"),
        (("kestrel", "motif_filtering", "position_threshold"), 0, "position_threshold"),
        (("kestrel", "motif_filtering", "alt_for_motif_right_gg"), 7, "alt_for_motif_right_gg"),
        (("kestrel", "motif_filtering", "alt_for_motif_right_gg"), "", "alt_for_motif_right_gg"),
        (("advntr", "settings", "max_frameshift"), True, "max_frameshift"),
        (("advntr", "settings", "max_frameshift"), 1.0, "max_frameshift"),
        (("advntr", "settings", "max_frameshift"), 0, "max_frameshift"),
        (("advntr", "flagging_rules", "Polymorphic_Call", "all"), {}, "Polymorphic_Call"),
        (
            ("advntr", "flagging_rules", "Polymorphic_Call", "all"),
            7,
            "must be a non-empty list",
        ),
        (
            ("advntr", "flagging_rules", "Polymorphic_Call", "all"),
            "x",
            "must be a non-empty list",
        ),
        (("advntr", "flagging_rules", "Polymorphic_Call", "all"), [], "Polymorphic_Call"),
    ],
)
def test_each_compound_semantic_guard_rejects_its_individual_boundaries(
    path: tuple[str | int, ...],
    value: object,
    message: str,
) -> None:
    components = _components()
    _replace(components, path, value)

    with pytest.raises(ValueError, match=message):
        validate_component_semantics(components)


@pytest.mark.parametrize(
    "path",
    [
        ("kestrel", "motif_filtering", "position_threshold"),
        ("advntr", "settings", "max_frameshift"),
        ("advntr", "settings", "frameshift_multiplier"),
        ("advntr", "settings", "vid"),
    ],
)
def test_positive_integer_semantic_boundaries_include_one(path: tuple[str | int, ...]) -> None:
    components = _components()
    _replace(components, path, 1)

    validate_component_semantics(components)
