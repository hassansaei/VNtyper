"""Pure shipped-verdict projection from final fixed caller rows."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells
from vntyper.scripts.molecular_identity import parse_molecular_identity, serialize_molecular_identity
from vntyper.scripts.reconciled_identity_persistence import (
    RECONCILED_IDENTITY_COLUMN,
    decode_reconciled_identity,
)

_KESTREL_NEGATIVE = {
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


def is_kestrel_negative_placeholder(row: Mapping[str, object]) -> bool:
    """Recognize and strictly validate the production Kestrel Negative row.

    Args:
        row: One parsed final Kestrel row.

    Returns:
        ``True`` for the exact production placeholder, otherwise ``False``.

    Raises:
        ValueError: If a row claims Negative confidence but differs from the
            closed placeholder schema or values.
    """
    if row.get("Confidence") != "Negative":
        return False
    if dict(row) != _KESTREL_NEGATIVE:
        raise ValueError("calibration Kestrel Negative row differs from the exact production placeholder")
    return True


def retained_kestrel_ordinals(rows: Sequence[Mapping[str, object]]) -> tuple[int, ...]:
    """Return the exact selected ordinals consulted by final Kestrel candidates.

    Args:
        rows: Independently parsed final Kestrel rows.

    Returns:
        Unique increasing retained observation ordinals, or an empty tuple for
        the exact Negative placeholder.

    Raises:
        ValueError: If final rows mix placeholders and findings or repeat an ordinal.
    """
    if len(rows) == 1 and is_kestrel_negative_placeholder(rows[0]):
        return ()
    if any(is_kestrel_negative_placeholder(row) for row in rows):
        raise ValueError("calibration Kestrel result cannot mix Negative and finding rows")
    ordinals = tuple(parse_selected_candidate_cells(row).selected_observation_ordinal for row in rows)
    if len(ordinals) != len(set(ordinals)):
        raise ValueError("calibration final Kestrel selected ordinals must be unique")
    return tuple(sorted(ordinals))


def build_shipped_projection(
    manifest_key: str,
    kestrel_rows: Sequence[Mapping[str, object]],
    advntr_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Build one whole-locus baseline row from both final caller projections.

    Args:
        manifest_key: Opaque study-member key.
        kestrel_rows: Final parsed Kestrel rows from one representation.
        advntr_rows: Final parsed adVNTR rows from the same representation.

    Returns:
        Closed whole-locus verdict plus a non-feature identity rendering map.

    Raises:
        ValueError: If identities, rendered verdicts, placeholders, or selected
            whole-locus projection are incomplete or ambiguous.
    """
    entries = [*_kestrel_entries(kestrel_rows), *_advntr_entries(advntr_rows)]
    if not entries:
        return _empty_projection(manifest_key)
    rendered = {(entry[2], entry[3]) for entry in entries}
    if len(rendered) != 1:
        raise ValueError("calibration final caller rows disagree on the shipped whole-locus name and tier")
    selected_name, selected_tier = next(iter(rendered))
    sources = {entry[4] for entry in entries}
    if sources == {"kestrel", "advntr"}:
        persisted = {decode_reconciled_identity(entry[5].get(RECONCILED_IDENTITY_COLUMN)) for entry in entries}
        if len(persisted) != 1 or None in persisted:
            raise ValueError("calibration cross-caller result requires one persisted reconciled canonical identity")
        selected_identity = serialize_molecular_identity(next(iter(persisted)))  # type: ignore[arg-type]
    else:
        identities = {entry[0] for entry in entries}
        if len(identities) != 1:
            raise ValueError("calibration single-caller result requires one canonical identity")
        selected_identity = next(iter(identities))
    if not any(entry[0] == selected_identity and entry[1] == selected_name for entry in entries):
        raise ValueError("calibration selected canonical identity disagrees with the rendered whole-locus verdict")
    projection: dict[str, dict[str, str | None]] = {}
    for identity, own_name, _global_name, tier, _source, _row in entries:
        value = {"name": own_name, "tier": tier}
        previous = projection.setdefault(identity, value)
        if previous != value:
            raise ValueError("calibration canonical identity has inconsistent caller rendering")
    selected_entries = [entry for entry in entries if entry[0] == selected_identity and entry[1] == selected_name]
    selected_entries.sort(key=lambda entry: 0 if entry[4] == "kestrel" else 1)
    source = selected_entries[0][4]
    row = selected_entries[0][5]
    support_column = "Estimated_Depth_AlternateVariant" if source == "kestrel" else "NumberOfSupportingReads"
    return {
        "manifest_key": manifest_key,
        "order": 0,
        "canonical_identity": selected_identity,
        "name": selected_name,
        "confidence": _optional(row.get("Confidence")),
        "flag": _optional(row.get("Flag")),
        "tier": selected_tier,
        "support": _support(row.get(support_column), source),
        "tie": False,
        "abstention": None,
        "identity_projection": {identity: projection[identity] for identity in sorted(projection)},
    }


def with_complete_kestrel_candidate_projections(
    shipped: Mapping[str, object],
    candidate_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Close a shipped rendering map over persisted complete Kestrel candidates.

    Args:
        shipped: Fixed whole-locus shipped projection.
        candidate_rows: Passing generated-run pre-result rows with fixed names.

    Returns:
        A copy whose rendering map includes every typed candidate at the shipped tier.

    Raises:
        ValueError: If the shipped map or a candidate's fixed presentation is malformed
            or conflicts with an existing canonical projection.
    """
    raw_projection = shipped.get("identity_projection")
    if not isinstance(raw_projection, Mapping):
        raise ValueError("calibration shipped identity projection must be an object")
    if any(not isinstance(value, Mapping) for value in raw_projection.values()):
        raise ValueError("calibration shipped candidate projections must be objects")
    projection = {str(identity): dict(value) for identity, value in raw_projection.items()}  # type: ignore[arg-type]
    tier = shipped.get("tier")
    if tier is not None and (not isinstance(tier, str) or not tier):
        raise ValueError("calibration shipped candidate tier must be non-empty or null")
    for row in candidate_rows:
        persisted = parse_selected_candidate_cells(row)
        identity = persisted.translation.identity
        if identity is None:
            continue
        serialized = serialize_molecular_identity(identity)
        name = row.get("Nomenclature")
        if not isinstance(name, str) or not name:
            raise ValueError("calibration complete Kestrel candidate requires one fixed presentation")
        value = {"name": name, "tier": tier}
        previous = projection.setdefault(serialized, value)
        if previous != value:
            raise ValueError("calibration canonical candidate has conflicting fixed presentations")
    result = dict(shipped)
    result["identity_projection"] = {identity: projection[identity] for identity in sorted(projection)}
    return result


def _kestrel_entries(
    rows: Sequence[Mapping[str, object]],
) -> list[tuple[str, str, str, str | None, str, Mapping[str, object]]]:
    if len(rows) == 1 and is_kestrel_negative_placeholder(rows[0]):
        return []
    entries = []
    for row in rows:
        if is_kestrel_negative_placeholder(row):
            raise ValueError("calibration Kestrel result cannot mix Negative and finding rows")
        persisted = parse_selected_candidate_cells(row)
        if persisted.translation.identity is None:
            raise ValueError("calibration selectable Kestrel row requires a canonical identity")
        entries.append(
            (
                serialize_molecular_identity(persisted.translation.identity),
                _required(row.get("Nomenclature_Kestrel"), "Kestrel caller name"),
                _required(row.get("Nomenclature"), "Kestrel whole-locus name"),
                _optional(row.get("Nomenclature_Tier")),
                "kestrel",
                row,
            )
        )
    return entries


def _advntr_entries(
    rows: Sequence[Mapping[str, object]],
) -> list[tuple[str, str, str, str | None, str, Mapping[str, object]]]:
    if len(rows) == 1 and rows[0].get("VID") == "Negative":
        return []
    entries = []
    for row in rows:
        if row.get("VID") == "Negative":
            raise ValueError("calibration adVNTR result cannot mix Negative and finding rows")
        identity_text = _required(row.get("Molecular_Identity"), "adVNTR canonical identity")
        identity = parse_molecular_identity(identity_text)
        if serialize_molecular_identity(identity) != identity_text:
            raise ValueError("calibration adVNTR canonical identity must use canonical serialization")
        entries.append(
            (
                identity_text,
                _required(row.get("Nomenclature_adVNTR"), "adVNTR caller name"),
                _required(row.get("Nomenclature"), "adVNTR whole-locus name"),
                _optional(row.get("Nomenclature_Tier")),
                "advntr",
                row,
            )
        )
    return entries


def _empty_projection(manifest_key: str) -> dict[str, object]:
    return {
        "manifest_key": manifest_key,
        "order": 0,
        "canonical_identity": None,
        "name": None,
        "confidence": None,
        "flag": None,
        "tier": None,
        "support": None,
        "tie": False,
        "abstention": None,
        "identity_projection": {},
    }


def _required(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"calibration {label} must be a non-empty string")
    return value


def _optional(value: object) -> str | None:
    return value if isinstance(value, str) and value else None


def _support(value: object, source: str) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise ValueError(f"calibration {source} baseline support must be numeric")
    if isinstance(value, str) and ("_" in value or value.strip() != value or not value):
        raise ValueError(f"calibration {source} baseline support must be numeric")
    try:
        if isinstance(value, float) or (isinstance(value, str) and any(char in value for char in ".eE")):
            number: int | float = float(value)
        else:
            number = int(value)
    except ValueError as error:
        raise ValueError(f"calibration {source} baseline support must be numeric") from error
    if not math.isfinite(number):
        raise ValueError(f"calibration {source} baseline support must be finite")
    if number < 0:
        raise ValueError(f"calibration {source} baseline support must be non-negative")
    return number
