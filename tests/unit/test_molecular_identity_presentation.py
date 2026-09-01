"""Compatibility presentation of the public molecular-identity quartet."""

from __future__ import annotations

import pytest

from vntyper.scripts import molecular_identity_presentation as presentation

pytestmark = pytest.mark.unit

IDENTITY_VALUES = {
    "Molecular_Identity": "muc1-mi-v1|INS|c.59_60|C",
    "Molecular_Identity_Status": "legacy-selected-among-multiple",
    "Equivalent_Representation_Count": 2,
    "Identity_Hypothesis_Count": 3,
}


def test_a_complete_quartet_is_copied_without_recalculation() -> None:
    """Catch replacing recorded counts or status with values inferred from legacy alleles."""
    row = {"POS": 67, "REF": "G", "ALT": "GG", "Nomenclature": "misleading", **IDENTITY_VALUES}

    projected = presentation.identity_compatibility_cells(row, schema_version=2)

    assert projected == IDENTITY_VALUES
    assert list(projected) == list(presentation.IDENTITY_COLUMNS)
    assert projected["Equivalent_Representation_Count"] == 2
    assert projected["Identity_Hypothesis_Count"] == 3


def test_an_unresolved_complete_quartet_preserves_empty_identity_zero_and_nonzero_hypotheses() -> None:
    """Catch truthiness defaults that replace the three meaningful falsy/unequal values."""
    recorded = {
        "Molecular_Identity": "",
        "Molecular_Identity_Status": "unresolved",
        "Equivalent_Representation_Count": 0,
        "Identity_Hypothesis_Count": 4,
    }

    assert presentation.identity_compatibility_cells(recorded, schema_version=1) == recorded


@pytest.mark.parametrize("schema_version", [1, 2])
def test_legacy_schema_rows_do_not_infer_identity_from_misleading_cells(schema_version: int) -> None:
    """Catch reconstruction from POS/REF/ALT/Nomenclature in either readable legacy schema."""
    row = {
        "POS": 67,
        "REF": "G",
        "ALT": "GG",
        "Variant": "Insertion",
        "Nomenclature": "59dupC",
    }

    projected = presentation.identity_compatibility_cells(row, schema_version=schema_version)

    assert projected == dict.fromkeys(presentation.IDENTITY_COLUMNS, "legacy identity not recorded")


@pytest.mark.parametrize("missing", presentation.IDENTITY_COLUMNS)
def test_removing_any_one_recorded_field_makes_the_entire_quartet_legacy(missing: str) -> None:
    """Catch partial completion or inference from the three identity cells that remain."""
    partial = {key: value for key, value in IDENTITY_VALUES.items() if key != missing}
    partial.update({"POS": 67, "REF": "G", "ALT": "GG", "Nomenclature": "59dupC"})

    projected = presentation.identity_compatibility_cells(partial, schema_version=2)

    assert projected == dict.fromkeys(presentation.IDENTITY_COLUMNS, "legacy identity not recorded")


def test_downstream_projection_appends_the_exact_quartet_and_leaves_source_unchanged() -> None:
    """Catch an omitted label, reordered label, or in-place mutation of summary data."""
    row = {"Sample": "s1", "POS": 67, **IDENTITY_VALUES}

    projected = presentation.identity_compatible_result_row(row, schema_version=2, positive=True)

    assert list(projected)[-4:] == list(presentation.IDENTITY_COLUMNS)
    assert projected == row
    assert projected is not row


def test_downstream_projection_does_not_widen_a_negative_caller_row() -> None:
    """Catch adding the quartet to the frozen caller-negative schema downstream."""
    negative = {"VID": "Negative", "Variant": "None", "Flag": "Not applicable"}

    projected = presentation.identity_compatible_result_row(negative, schema_version=2, positive=False)

    assert projected == negative
    assert not set(presentation.IDENTITY_COLUMNS).intersection(projected)
