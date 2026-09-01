"""Pure projection of the Kestrel selection decision component."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class KestrelSortField:
    """One logical Kestrel selection sort field."""

    column: str
    ascending: bool


@dataclass(frozen=True)
class KestrelSelection:
    """Typed values consumed by Kestrel frame filtering and row selection."""

    confidence_priority: Mapping[str, int]
    final_filter_columns: tuple[str, ...]
    modulus: int
    insertion_remainder: int
    deletion_remainder: int
    sort_order: tuple[KestrelSortField, ...]
    unflagged_value: str


def project_kestrel_selection(selection: Mapping[str, object]) -> KestrelSelection:
    """Project a validated profile selection mapping into exact stage arguments.

    Args:
        selection: Complete ``components.kestrel.selection`` mapping.

    Returns:
        Typed immutable Kestrel selection values.

    Raises:
        ValueError: If a direct compatibility caller supplies an incomplete mapping.
    """
    try:
        priority_raw = selection["confidence_priority"]
        columns_raw = selection["final_filter_columns"]
        frameshift_raw = selection["frameshift"]
        sort_raw = selection["sort_order"]
        unflagged_value = selection["unflagged_value"]
    except KeyError as exc:
        raise ValueError(f"Kestrel selection is missing {exc.args[0]!r}") from exc
    if not isinstance(priority_raw, Mapping) or not isinstance(frameshift_raw, Mapping):
        raise ValueError("Kestrel selection priority and frameshift values must be mappings")
    if not isinstance(columns_raw, Sequence) or isinstance(columns_raw, str):
        raise ValueError("Kestrel final_filter_columns must be an ordered sequence")
    if not isinstance(sort_raw, Sequence) or isinstance(sort_raw, str):
        raise ValueError("Kestrel sort_order must be an ordered sequence")
    if not isinstance(unflagged_value, str):
        raise ValueError("Kestrel unflagged_value must be a string")

    priority = {str(label): int(value) for label, value in priority_raw.items()}
    sort_order: list[KestrelSortField] = []
    for item in sort_raw:
        if not isinstance(item, Mapping):
            raise ValueError("Kestrel sort_order entries must be mappings")
        column = item.get("column")
        ascending = item.get("ascending")
        if not isinstance(column, str) or not isinstance(ascending, bool):
            raise ValueError("Kestrel sort_order entries require string column and boolean ascending")
        sort_order.append(KestrelSortField(column=column, ascending=ascending))

    return KestrelSelection(
        confidence_priority=priority,
        final_filter_columns=tuple(str(column) for column in columns_raw),
        modulus=int(frameshift_raw["modulus"]),
        insertion_remainder=int(frameshift_raw["insertion_remainder"]),
        deletion_remainder=int(frameshift_raw["deletion_remainder"]),
        sort_order=tuple(sort_order),
        unflagged_value=unflagged_value,
    )
