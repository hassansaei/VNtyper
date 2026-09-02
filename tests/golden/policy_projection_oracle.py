"""Independent policy-stable public projections for the molecular golden corpus."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping

# Independent literal transcription of the governed #267 recurrent-State evidence.
# Never import the production artifact here: this oracle must detect its drift.
RECURRENT_STATE_EVIDENCE = frozenset(
    {
        "I10_2_A_LEN1",
        "D8_2&D9_2&I9_2_A_LEN9",
        "D2_2&I2_2_C_LEN5",
        "I39_2_A_LEN4",
        "I52_2_A_LEN7",
        "I45_2_A_LEN4",
        "D45_2&I45_2_A_LEN2",
        "D14_2&I14_2_G_LEN14",
        "D58_2&D59_2",
        "I60_2_A_LEN10",
        "I14_2_G_LEN16",
        "I18_2_T_LEN1",
        "I21_2_G_LEN4",
        "D29_2&I29_2_A_LEN2",
        "D8_2&I8_2_A_LEN20",
        "D20_2&D21_2",
        "D21_2&D22_2",
        "I14_2_A_LEN1",
        "I11_2_G_LEN1",
        "I26_7_A_LEN25",
        "D17_2&D18_2&D19_2&D20_2&D21_2",
        "I14_2_C_LEN4",
        "I23_6_G_LEN1",
        "I21_2_T_LEN1",
    }
)
GOVERNED_ASSERTION = "A carried-forward recurrent adVNTR State is insufficient for molecular identity."
_POLICY_EXPLANATION_COLUMNS = frozenset({"Nomenclature_Flags", "Nomenclature_Note", "Evidence_Disposition"})
_NONPUBLIC_PERSISTENCE_COLUMNS = frozenset({"__Reconciled_Molecular_Identity"})


def selected_projection_fingerprint(projection: Mapping[str, tuple[tuple[str, ...], ...]]) -> str:
    """Return a deterministic fingerprint of every literal selected-row cell."""
    payload = json.dumps(projection, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def public_projection_fingerprint(projection: Mapping[str, object]) -> str:
    """Return a deterministic fingerprint for a nested public-row projection."""
    payload = json.dumps(projection, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def sample_sets_fingerprint(sample_sets: dict[str, frozenset[str]]) -> str:
    """Return a deterministic fingerprint of named literal sample-key sets."""
    payload = json.dumps(
        {name: sorted(keys) for name, keys in sample_sets.items()},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def record_policy_projections(
    stable: dict[str, tuple[tuple[tuple[str, str], ...], ...]],
    unaffected: dict[str, tuple[tuple[tuple[str, str], ...], ...]],
    explanations: dict[str, tuple[tuple[str, str, str, str], ...]],
    key: str,
    condition: str,
    public_rows: dict[str, tuple[dict[str, str], ...]],
) -> None:
    """Record B1-stable cells and independently scoped B2 explanation cells."""
    decision_key = f"{key}/{condition}"
    collision = any((row.get("Variant") or "").strip() in RECURRENT_STATE_EVIDENCE for row in public_rows["advntr"])
    for caller, rows in public_rows.items():
        projection_key = f"{decision_key}/{caller}"
        stable[projection_key] = tuple(
            tuple(
                sorted(
                    (column, value)
                    for column, value in row.items()
                    if column not in _POLICY_EXPLANATION_COLUMNS | _NONPUBLIC_PERSISTENCE_COLUMNS
                )
            )
            for row in rows
        )
        unaffected[projection_key] = tuple(
            tuple(
                sorted(
                    (column, value)
                    for column, value in row.items()
                    if column != "Evidence_Disposition" and column not in _NONPUBLIC_PERSISTENCE_COLUMNS
                )
            )
            for row in rows
            if not (
                (caller == "advntr" and (row.get("Variant") or "").strip() in RECURRENT_STATE_EVIDENCE)
                or (caller == "kestrel" and collision)
            )
        )
        explanations[projection_key] = tuple(
            (
                (row.get("Variant") or "").strip(),
                (row.get("Nomenclature_Flags") or "").strip(),
                (row.get("Nomenclature_Note") or "").strip(),
                row.get("Evidence_Disposition", "<missing>").strip(),
            )
            for row in rows
        )
