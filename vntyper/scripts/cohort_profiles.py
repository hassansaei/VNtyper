"""Pure decision-profile provenance projection and grouping for cohorts."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from vntyper.scripts.profile_provenance import (
    LEGACY_DECISION_PROFILE_REVISION,
    DecisionProfileProvenance,
)

PROFILE_ID_COLUMN = "Decision_Profile_ID"
PROFILE_REVISION_COLUMN = "Decision_Profile_Revision"
PROFILE_SHA256_COLUMN = "Decision_Profile_SHA256"
PROFILE_EXPORT_COLUMNS = (PROFILE_ID_COLUMN, PROFILE_REVISION_COLUMN, PROFILE_SHA256_COLUMN)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class CohortDecisionProfileGroup:
    """Samples whose runs recorded one decision-profile SHA-256."""

    profile_id: str
    revision: str
    sha256: str
    samples: tuple[str, ...]


def profile_export_fields(provenance: DecisionProfileProvenance) -> dict[str, str]:
    """Project verified run provenance onto the three cohort export columns.

    Args:
        provenance: Verified schema-3 provenance or an explicit legacy value.

    Returns:
        The profile ID, revision, and SHA-256 export cells.

    Raises:
        ValueError: If the value did not come from summary-profile resolution.
    """
    if not isinstance(provenance, DecisionProfileProvenance):
        raise ValueError("cohort decision-profile fields require verified provenance")
    legacy = LEGACY_DECISION_PROFILE_REVISION
    return {
        PROFILE_ID_COLUMN: provenance.profile_id or legacy,
        PROFILE_REVISION_COLUMN: provenance.revision,
        PROFILE_SHA256_COLUMN: provenance.sha256 or legacy,
    }


def annotate_profile_rows(rows: Sequence[Mapping[str, object]], fields: Mapping[str, str]) -> list[dict[str, object]]:
    """Append one sample's verified profile fields to each result row.

    Args:
        rows: Parsed result rows for one sample.
        fields: Exact fields returned by :func:`profile_export_fields`.

    Returns:
        Copied rows carrying the three profile export fields.

    Raises:
        ValueError: If the fields are incomplete or malformed.
    """
    if set(fields) != set(PROFILE_EXPORT_COLUMNS) or any(
        not isinstance(value, str) or not value for value in fields.values()
    ):
        raise ValueError("cohort decision-profile export fields are incomplete")
    return [{**row, **fields} for row in rows]


def _required_entry_string(entry: Mapping[str, object], field: str) -> str:
    value = entry.get(field)
    if not isinstance(value, str) or not value:
        raise ValueError(f"cohort decision-profile group requires non-empty {field}")
    return value


def group_decision_profiles(entries: Sequence[Mapping[str, object]]) -> tuple[CohortDecisionProfileGroup, ...]:
    """Group sample provenance by exact SHA-256 without pooling metadata.

    Args:
        entries: Per-sample mappings carrying ``Sample`` and the three export fields.

    Returns:
        Groups in first-seen hash order, with samples in input order.

    Raises:
        ValueError: If fields are missing or one hash carries conflicting metadata.
    """
    grouped: dict[str, tuple[str, str, list[str]]] = {}
    for entry in entries:
        sample = _required_entry_string(entry, "Sample")
        profile_id = _required_entry_string(entry, PROFILE_ID_COLUMN)
        revision = _required_entry_string(entry, PROFILE_REVISION_COLUMN)
        sha256 = _required_entry_string(entry, PROFILE_SHA256_COLUMN)
        if sha256 != LEGACY_DECISION_PROFILE_REVISION and _SHA256.fullmatch(sha256) is None:
            raise ValueError("cohort decision-profile SHA-256 must be verified lowercase hexadecimal or legacy text")
        existing = grouped.get(sha256)
        if existing is None:
            grouped[sha256] = (profile_id, revision, [sample])
            continue
        if existing[:2] != (profile_id, revision):
            raise ValueError(f"decision-profile hash {sha256} has conflicting metadata across cohort samples")
        existing[2].append(sample)
    return tuple(
        CohortDecisionProfileGroup(profile_id, revision, sha256, tuple(samples))
        for sha256, (profile_id, revision, samples) in grouped.items()
    )
