"""Run-local decision-profile provenance and summary-schema validation."""

from __future__ import annotations

import hashlib
import logging
import os
import re
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from vntyper.scripts.artifact_names import DECISION_PROFILE_SNAPSHOT_RELATIVE
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.decision_profile import (
    ProfileKind,
    ProfileSource,
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)
from vntyper.scripts.decision_profile_schema import component_projection, validate_complete_inventory
from vntyper.scripts.molecular_identity import parse_molecular_identity, serialize_molecular_identity
from vntyper.scripts.molecular_identity_presentation import IDENTITY_COLUMNS
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_KESTREL

logger = logging.getLogger(__name__)

LEGACY_DECISION_PROFILE_REVISION = "decision profile not recorded by legacy run"
HISTORICAL_REVISION_1_PACKAGED_SHA256 = "be6329fb12107a1b6b65e425257be6233c7e2115e299e941c12a63a6a6d59718"

_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_PROFILE_FIELDS = frozenset(
    {
        "decision_profile_id",
        "decision_profile_revision",
        "decision_profile_kind",
        "decision_profile_source",
        "decision_profile_sha256",
        "decision_profile_snapshot",
    }
)
_PROFILE_KINDS = frozenset({"packaged", "explicit-custom", "generated"})
_PROFILE_SOURCES = frozenset({"package", "explicit-cli"})
_IDENTITY_STATUSES = frozenset({"unique", "legacy-selected-among-multiple", "unresolved"})
_EVIDENCE_DISPOSITIONS = frozenset({"admissible", "identity-insufficient"})


@dataclass(frozen=True)
class DecisionProfileProvenance:
    """Profile identity recorded by a run after its snapshot is verified."""

    profile_id: str | None
    profile_revision: str | None
    profile_kind: ProfileKind | None
    source: ProfileSource | None
    sha256: str | None
    snapshot_path: str | None

    @property
    def revision(self) -> str:
        """Return the recorded revision or explicit legacy compatibility text."""
        return self.profile_revision or LEGACY_DECISION_PROFILE_REVISION


def profile_summary_fields(profile: ResolvedDecisionProfile) -> dict[str, str]:
    """Project one resolved profile onto the complete schema-3 provenance fields.

    Args:
        profile: Profile resolved for the current run.

    Returns:
        Exact top-level summary fields, with a fixed relative snapshot path.

    Raises:
        ValueError: If ``profile`` is not a resolved decision profile.
    """
    if not isinstance(profile, ResolvedDecisionProfile):
        raise ValueError("Decision profile summary provenance requires a resolved profile")
    return {
        "decision_profile_id": profile.profile_id,
        "decision_profile_revision": profile.profile_revision,
        "decision_profile_kind": profile.profile_kind,
        "decision_profile_source": profile.source,
        "decision_profile_sha256": profile.digest,
        "decision_profile_snapshot": DECISION_PROFILE_SNAPSHOT_RELATIVE.as_posix(),
    }


def snapshot_decision_profile(profile: ResolvedDecisionProfile, destination: str | Path) -> None:
    """Atomically snapshot the exact canonical bytes resolved for a run.

    Args:
        profile: Profile resolved before stage artifacts are created.
        destination: Run-local final snapshot path.

    Raises:
        ValueError: If ``profile`` is not a resolved profile.
        RuntimeError: If the canonical candidate cannot be written or installed.
    """
    if not isinstance(profile, ResolvedDecisionProfile):
        raise ValueError("Decision profile snapshot requires a resolved profile")
    destination_path = Path(destination)
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    installed = False
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination_path.parent,
            prefix=f".{destination_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as candidate:
            temporary_path = Path(candidate.name)
            os.fchmod(candidate.fileno(), 0o644)
            candidate.write(profile.canonical_bytes)
            candidate.flush()
            os.fsync(candidate.fileno())
        os.replace(temporary_path, destination_path)
        installed = True
    except Exception as error:
        message = f"Failed to snapshot decision profile at {destination_path}: {error}"
        logger.error(message)
        raise RuntimeError(message) from error
    finally:
        if temporary_path is not None and not installed:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError as cleanup_error:
                logger.error(f"Failed to remove incomplete decision profile snapshot {temporary_path}: {cleanup_error}")


def _nonempty_string(value: object, *, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"summary schema 3 {field} must be a non-empty string")
    return value


def _parse_recorded_provenance(summary: Mapping[str, object]) -> DecisionProfileProvenance:
    present = {key for key in summary if isinstance(key, str) and key.startswith("decision_profile_")}
    missing = _PROFILE_FIELDS - present
    unknown = present - _PROFILE_FIELDS
    if missing or unknown:
        raise ValueError(
            "summary schema 3 required decision-profile fields differ: "
            f"missing {sorted(missing)}, unknown {sorted(unknown)}"
        )

    profile_id = _nonempty_string(summary["decision_profile_id"], field="decision_profile_id")
    revision = _nonempty_string(summary["decision_profile_revision"], field="decision_profile_revision")
    kind_value = _nonempty_string(summary["decision_profile_kind"], field="decision_profile_kind")
    source_value = _nonempty_string(summary["decision_profile_source"], field="decision_profile_source")
    digest = _nonempty_string(summary["decision_profile_sha256"], field="decision_profile_sha256")
    snapshot = _nonempty_string(summary["decision_profile_snapshot"], field="decision_profile_snapshot")
    if kind_value not in _PROFILE_KINDS:
        raise ValueError(f"summary schema 3 decision profile kind is unsupported: {kind_value}")
    if source_value not in _PROFILE_SOURCES:
        raise ValueError(f"summary schema 3 decision profile source is unsupported: {source_value}")
    if source_value == "package" and kind_value != "packaged":
        raise ValueError("summary schema 3 package source requires packaged decision profile kind")
    if source_value == "explicit-cli" and kind_value not in {"explicit-custom", "generated"}:
        raise ValueError("summary schema 3 explicit-cli source requires explicit-custom or generated profile kind")
    if _SHA256_RE.fullmatch(digest) is None:
        raise ValueError("summary schema 3 decision profile SHA-256 must be 64 lowercase hexadecimal characters")
    expected_snapshot = DECISION_PROFILE_SNAPSHOT_RELATIVE.as_posix()
    if Path(snapshot).is_absolute() or snapshot != expected_snapshot:
        raise ValueError(f"summary schema 3 relative snapshot must be exactly {expected_snapshot}")
    return DecisionProfileProvenance(
        profile_id=profile_id,
        profile_revision=revision,
        profile_kind=cast(ProfileKind, kind_value),
        source=cast(ProfileSource, source_value),
        sha256=digest,
        snapshot_path=snapshot,
    )


def _parse_count(value: object, *, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"schema 3 identity {field} must be a non-negative integer")
    if isinstance(value, int):
        count = value
    elif isinstance(value, str) and value.isdigit() and (value == "0" or not value.startswith("0")):
        count = int(value)
    else:
        raise ValueError(f"schema 3 identity {field} must be a non-negative integer")
    if count < 0:
        raise ValueError(f"schema 3 identity {field} must be a non-negative integer")
    return count


def _validate_identity_row(row: Mapping[str, object], *, positive: bool) -> None:
    present = set(IDENTITY_COLUMNS) & set(row)
    if not positive:
        if present:
            raise ValueError("summary schema 3 negative caller rows must retain their legacy schema")
        return
    if present != set(IDENTITY_COLUMNS):
        raise ValueError("summary schema 3 positive caller row requires all four molecular identity fields")

    identity_value = row["Molecular_Identity"]
    status = row["Molecular_Identity_Status"]
    if not isinstance(identity_value, str):
        raise ValueError("summary schema 3 Molecular_Identity must be a string")
    if not isinstance(status, str) or status not in _IDENTITY_STATUSES:
        raise ValueError(f"summary schema 3 molecular identity status is unsupported: {status}")
    equivalent_count = _parse_count(row["Equivalent_Representation_Count"], field="Equivalent_Representation_Count")
    hypothesis_count = _parse_count(row["Identity_Hypothesis_Count"], field="Identity_Hypothesis_Count")

    if status == "unresolved":
        if identity_value or equivalent_count != 0:
            raise ValueError("summary schema 3 unresolved identity requires an empty identity and zero representations")
        return
    if not identity_value or equivalent_count == 0 or hypothesis_count == 0:
        raise ValueError("summary schema 3 resolved identity requires identity and positive counts")
    identity = parse_molecular_identity(identity_value)
    if serialize_molecular_identity(identity) != identity_value:
        raise ValueError("summary schema 3 molecular identity serialization is not canonical")
    if status == "unique" and hypothesis_count != 1:
        raise ValueError("summary schema 3 unique identity requires exactly one identity hypothesis")
    if status == "legacy-selected-among-multiple" and hypothesis_count <= 1:
        raise ValueError("summary schema 3 legacy-selected-among-multiple identity requires multiple hypotheses")


def _validate_schema_three_identity_rows(
    summary: Mapping[str, object],
    advntr_component: Mapping[str, object],
    evidence_digest: str | None,
) -> None:
    artifact_evidence = advntr_component.get("artifact_evidence")
    if not isinstance(artifact_evidence, Mapping):
        raise ValueError("recorded decision profile adVNTR artifact evidence must be an object")
    active_states_raw = artifact_evidence.get("active_states")
    if not isinstance(active_states_raw, (list, tuple)) or any(
        not isinstance(state, str) or not state for state in active_states_raw
    ):
        raise ValueError("recorded decision profile adVNTR active states must be non-empty strings")
    active_states = frozenset(active_states_raw)
    matching_disposition = artifact_evidence.get("matching_disposition")
    nonmatching_disposition = artifact_evidence.get("nonmatching_disposition")
    if matching_disposition not in _EVIDENCE_DISPOSITIONS or nonmatching_disposition not in _EVIDENCE_DISPOSITIONS:
        raise ValueError("recorded decision profile adVNTR dispositions are unsupported")
    raw_steps = summary.get("steps", [])
    if not isinstance(raw_steps, list):
        raise ValueError("summary schema 3 steps must be an array")
    from vntyper.scripts.report_formatting import is_empty_result_row

    for raw_step in raw_steps:
        if not isinstance(raw_step, Mapping):
            raise ValueError("summary schema 3 step must be an object")
        step_name = raw_step.get("step")
        if step_name not in {STEP_KESTREL, STEP_ADVNTR}:
            continue
        if step_name == STEP_ADVNTR and evidence_digest is None:
            raise ValueError("summary schema 3 adVNTR caller step requires a non-null evidence digest")
        parsed_result = raw_step.get("parsed_result")
        if not isinstance(parsed_result, Mapping):
            raise ValueError("summary schema 3 caller parsed_result must be an object")
        rows = parsed_result.get("data", [])
        if not isinstance(rows, list):
            raise ValueError("summary schema 3 caller step data must be an array")
        for raw_row in rows:
            if not isinstance(raw_row, Mapping):
                raise ValueError("summary schema 3 caller result row must be an object")
            caller_row = {key: value for key, value in raw_row.items() if key not in IDENTITY_COLUMNS}
            positive = (
                not is_empty_result_row(caller_row) if step_name == STEP_KESTREL else raw_row.get("VID") != "Negative"
            )
            if step_name == STEP_ADVNTR and positive:
                disposition = raw_row.get("Evidence_Disposition")
                if disposition not in _EVIDENCE_DISPOSITIONS:
                    raise ValueError("summary schema 3 positive adVNTR row requires a supported Evidence_Disposition")
                state = raw_row.get("Variant")
                if not isinstance(state, str) or not state:
                    raise ValueError("summary schema 3 positive adVNTR row requires a non-empty Variant state")
                expected_disposition = matching_disposition if state in active_states else nonmatching_disposition
                if disposition != expected_disposition:
                    raise ValueError(
                        "summary schema 3 adVNTR Evidence_Disposition differs from the governed Variant state"
                    )
            _validate_identity_row(raw_row, positive=positive)


def verify_profile_snapshot(
    provenance: DecisionProfileProvenance,
    snapshot_path: str | Path,
    *,
    advntr_evidence_digest: str | None = None,
    schema_three_summary: Mapping[str, object] | None = None,
) -> DecisionProfileProvenance:
    """Verify one recorded snapshot without substituting a current profile.

    Args:
        provenance: Complete schema-3 values recorded by the run.
        snapshot_path: Run-local snapshot resolved from its fixed relative path.
        advntr_evidence_digest: Recorded PR-B evidence digest, or ``None`` when adVNTR did not run.
        schema_three_summary: Complete summary whose caller rows must be checked against the verified profile.

    Returns:
        The same provenance after strict bytes, hash, schema, and metadata checks.

    Raises:
        OSError: If the recorded snapshot cannot be read.
        ValueError: If bytes or recorded metadata are inconsistent.
    """
    if not isinstance(provenance, DecisionProfileProvenance) or provenance.sha256 is None:
        raise ValueError("Decision profile snapshot verification requires recorded schema-3 provenance")
    raw = Path(snapshot_path).read_bytes()
    document = load_strict_json_object(raw)
    canonical = canonical_json_bytes(document)
    if raw != canonical:
        raise ValueError("recorded decision profile snapshot must use canonical RFC 8785 bytes plus one newline")
    actual_digest = hashlib.sha256(raw).hexdigest()
    if actual_digest != provenance.sha256:
        raise ValueError(f"decision profile digest mismatch: expected {provenance.sha256}, got {actual_digest}")

    if provenance.source == "package":
        validate_complete_inventory(document)
        packaged = load_packaged_decision_profile()
        if document.get("profile_revision") == "1":
            if actual_digest != HISTORICAL_REVISION_1_PACKAGED_SHA256:
                raise ValueError("recorded package-source snapshot does not match the checked-in packaged profile")
        elif raw != packaged.canonical_bytes:
            raise ValueError("recorded package-source snapshot does not match the checked-in packaged profile")
        advntr_component = component_projection(document, "advntr")

        profile_id = document["profile_id"]
        revision = document["profile_revision"]
        kind = document["profile_kind"]
    else:
        packaged = load_packaged_decision_profile()
        resolved = parse_decision_profile(raw, packaged_document=packaged.document)
        advntr_component = resolved.components["advntr"]
        profile_id = resolved.profile_id
        revision = resolved.profile_revision
        kind = resolved.profile_kind
    if (profile_id, revision, kind) != (
        provenance.profile_id,
        provenance.profile_revision,
        provenance.profile_kind,
    ):
        raise ValueError("recorded decision profile metadata differs from its verified snapshot")
    if not isinstance(advntr_component, Mapping):
        raise ValueError("recorded decision profile adVNTR component must be an object")
    if advntr_evidence_digest is not None:
        artifact_evidence = advntr_component.get("artifact_evidence")
        if not isinstance(artifact_evidence, Mapping):
            raise ValueError("recorded decision profile adVNTR artifact evidence must be an object")
        if artifact_evidence.get("digest") != advntr_evidence_digest:
            raise ValueError("recorded profile adVNTR evidence digest differs from the summary evidence digest")
    if schema_three_summary is not None:
        _validate_schema_three_identity_rows(schema_three_summary, advntr_component, advntr_evidence_digest)
    return provenance


def resolve_summary_profile(
    summary: Mapping[str, object],
    run_root: str | Path,
) -> DecisionProfileProvenance:
    """Resolve only the profile provenance recorded by a pipeline summary.

    Args:
        summary: Parsed ``pipeline_summary.json`` mapping.
        run_root: Directory containing the summary and its provenance directory.

    Returns:
        Verified schema-3 provenance, or explicit legacy provenance for schema 1/2.

    Raises:
        OSError: If a schema-3 snapshot cannot be read.
        ValueError: If the summary or snapshot violates the closed contract.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("pipeline summary must be an object")
    schema_version = summary.get("schema_version")
    if isinstance(schema_version, bool):
        raise ValueError(f"unsupported pipeline summary schema version: {schema_version}")
    if schema_version is None or (isinstance(schema_version, int) and schema_version in {1, 2}):
        forbidden = sorted(_PROFILE_FIELDS & set(summary))
        if forbidden:
            raise ValueError(f"legacy summary cannot contain schema 3 decision profile fields: {forbidden}")
        return DecisionProfileProvenance(None, None, None, None, None, None)
    if not isinstance(schema_version, int) or schema_version != 3:
        raise ValueError(f"unsupported pipeline summary schema version: {schema_version}")
    if "advntr_evidence_digest" not in summary:
        raise ValueError("summary schema 3 requires advntr_evidence_digest")
    evidence_digest = summary["advntr_evidence_digest"]
    if evidence_digest is not None and (
        not isinstance(evidence_digest, str) or _SHA256_RE.fullmatch(evidence_digest) is None
    ):
        raise ValueError("summary schema 3 adVNTR evidence digest must be 64 lowercase hexadecimal characters or null")
    if "decision_policy" not in summary:
        raise ValueError("summary schema 3 requires a decision_policy string")
    _nonempty_string(summary["decision_policy"], field="decision_policy")

    provenance = _parse_recorded_provenance(summary)
    assert provenance.snapshot_path is not None
    return verify_profile_snapshot(
        provenance,
        Path(run_root) / provenance.snapshot_path,
        advntr_evidence_digest=evidence_digest,
        schema_three_summary=summary,
    )
