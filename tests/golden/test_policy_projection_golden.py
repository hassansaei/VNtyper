"""Focused golden contracts for policy-stable public projections."""

import pytest

from tests.golden.policy_projection_oracle import public_projection_fingerprint, record_policy_projections

pytestmark = pytest.mark.golden


def test_policy_projection_excludes_only_the_internal_reconciled_identity_handoff() -> None:
    """Private identity persistence is not public, while an arbitrary new field remains locked."""

    def projections(
        row: dict[str, str],
    ) -> tuple[
        dict[str, tuple[tuple[tuple[str, str], ...], ...]],
        dict[str, tuple[tuple[tuple[str, str], ...], ...]],
    ]:
        stable: dict[str, tuple[tuple[tuple[str, str], ...], ...]] = {}
        unaffected: dict[str, tuple[tuple[tuple[str, str], ...], ...]] = {}
        explanations: dict[str, tuple[tuple[str, str, str, str], ...]] = {}
        record_policy_projections(
            stable,
            unaffected,
            explanations,
            "experiment/sample",
            "mutated",
            {"kestrel": (row,), "advntr": ()},
        )
        return stable, unaffected

    baseline = projections({"Nomenclature": "59dupC"})
    with_internal_handoff = projections(
        {
            "Nomenclature": "59dupC",
            "__Reconciled_Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|C",
        }
    )
    with_new_public_field = projections({"Nomenclature": "59dupC", "New_Public_Field": "changed"})

    assert with_internal_handoff == baseline
    assert public_projection_fingerprint(with_new_public_field[0]) != public_projection_fingerprint(baseline[0])
    assert public_projection_fingerprint(with_new_public_field[1]) != public_projection_fingerprint(baseline[1])
