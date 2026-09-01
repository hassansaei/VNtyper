"""Closed JSON codec and atomic I/O for Kestrel BAM replay evidence.

Current identity-aware runs write ``bam_identity_replay.v1.json`` with schema
``bam-identity-replay-v1``. Each deterministic candidate group records exactly one
``not-consulted``, ``unavailable``, or complete ``observed`` locus.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

from vntyper.scripts.molecular_identity import parse_molecular_identity, serialize_molecular_identity
from vntyper.scripts.nomenclature_bam_evidence import BamIdentityEvidence, BamLocusEvidence

BAM_REPLAY_FILENAME = "bam_identity_replay.v1.json"
BAM_REPLAY_SCHEMA_VERSION = "bam-identity-replay-v1"

BamReplayState = Literal["not-consulted", "unavailable", "observed"]
_REPLAY_STATES = frozenset({"not-consulted", "unavailable", "observed"})


@dataclass(frozen=True)
class BamReplayLocus:
    """Replay state for one exact persisted Kestrel candidate group.

    Attributes:
        candidate_observation_ordinals: Stable persisted A3 observation ordinals
            sharing the exact BAM locus and complete pair context.
        state: Whether BAM was skipped, unavailable, or observed.
        evidence: Complete evidence for ``observed`` only, including zero records.
    """

    candidate_observation_ordinals: tuple[int, ...]
    state: BamReplayState
    evidence: BamLocusEvidence | None

    def __post_init__(self) -> None:
        """Validate the closed state/evidence pairing."""
        ordinals = self.candidate_observation_ordinals
        if not isinstance(ordinals, tuple) or not ordinals:
            raise ValueError("BAM replay candidate observation ordinals must be a nonempty tuple")
        for ordinal in ordinals:
            _nonnegative_integer(ordinal, "BAM replay candidate observation ordinal")
        if ordinals != tuple(sorted(ordinals)) or len(ordinals) != len(set(ordinals)):
            raise ValueError("BAM replay candidate observation ordinals must be unique and increasing")
        if self.state not in _REPLAY_STATES:
            raise ValueError(f"Unsupported BAM replay state: {self.state}")
        if self.state == "observed":
            if not isinstance(self.evidence, BamLocusEvidence):
                raise ValueError("Observed BAM replay loci require BamLocusEvidence")
            bound_ordinals = {
                ordinal
                for record in self.evidence.records
                for bindings in record.candidate_observation_ordinals
                for ordinal in bindings
            }
            if not bound_ordinals <= set(ordinals):
                raise ValueError("BAM replay evidence bindings must belong to its candidate group")
        elif self.evidence is not None:
            raise ValueError("Unobserved BAM replay loci require null evidence")


@dataclass(frozen=True)
class BamReplayArtifact:
    """Versioned replay loci in unique deterministic ordinal order.

    Attributes:
        loci: One locus per exact persisted Kestrel candidate group.
    """

    loci: tuple[BamReplayLocus, ...]

    def __post_init__(self) -> None:
        """Validate immutable, unique, deterministic candidate groups."""
        if not isinstance(self.loci, tuple) or any(not isinstance(locus, BamReplayLocus) for locus in self.loci):
            raise ValueError("BAM replay loci must be a tuple of BamReplayLocus values")
        groups = tuple(locus.candidate_observation_ordinals for locus in self.loci)
        if groups != tuple(sorted(groups)):
            raise ValueError("BAM replay candidate groups must be deterministically increasing")
        ordinals = tuple(ordinal for group in groups for ordinal in group)
        if len(ordinals) != len(set(ordinals)):
            raise ValueError("BAM replay candidate groups must not overlap")


def validate_bam_replay_artifact_ordinals(
    artifact: BamReplayArtifact,
    expected_candidate_observation_ordinals: tuple[int, ...],
) -> None:
    """Require every persisted candidate ordinal to occur in exactly one group.

    Args:
        artifact: Validated grouped replay artifact.
        expected_candidate_observation_ordinals: Complete increasing persisted set.

    Raises:
        ValueError: If the expected ordinals are malformed or coverage differs.
    """
    if not isinstance(artifact, BamReplayArtifact):
        raise ValueError("BAM replay artifact must be a BamReplayArtifact")
    expected = expected_candidate_observation_ordinals
    if not isinstance(expected, tuple) or any(
        isinstance(ordinal, bool) or not isinstance(ordinal, int) or ordinal < 0 for ordinal in expected
    ):
        raise ValueError("Expected BAM replay candidate ordinals must be a tuple of non-negative integers")
    if expected != tuple(sorted(expected)) or len(expected) != len(set(expected)):
        raise ValueError("Expected BAM replay candidate ordinals must be unique and increasing")
    retained = tuple(sorted(ordinal for locus in artifact.loci for ordinal in locus.candidate_observation_ordinals))
    if retained != expected:
        raise ValueError("BAM replay groups must cover persisted candidate ordinals exactly")


def merge_bam_replay_artifacts(
    existing: BamReplayArtifact,
    current: BamReplayArtifact,
) -> BamReplayArtifact:
    """Merge a later lifecycle observation without downgrading prior evidence.

    Args:
        existing: Earlier Kestrel-stage replay artifact.
        current: Later reconciliation-stage artifact for the exact same groups.

    Returns:
        Merged artifact using ``observed > unavailable > not-consulted``. A new
        observed result replaces prior evidence for the same validated group.

    Raises:
        ValueError: If artifacts do not carry the exact same candidate partition.
    """
    if not isinstance(existing, BamReplayArtifact) or not isinstance(current, BamReplayArtifact):
        raise ValueError("BAM replay merge inputs must be BamReplayArtifact values")
    existing_groups = tuple(locus.candidate_observation_ordinals for locus in existing.loci)
    current_groups = tuple(locus.candidate_observation_ordinals for locus in current.loci)
    if existing_groups != current_groups:
        raise ValueError("BAM replay merge requires identical candidate groups")
    merged: list[BamReplayLocus] = []
    for previous, later in zip(existing.loci, current.loci, strict=True):
        if later.state == "observed" or previous.state == "not-consulted":
            merged.append(later)
        elif previous.state == "observed" or later.state == "not-consulted":
            merged.append(previous)
        else:
            merged.append(later)
    return BamReplayArtifact(tuple(merged))


def encode_bam_replay_artifact(artifact: BamReplayArtifact) -> dict[str, object]:
    """Encode replay evidence into canonical JSON-compatible primitives.

    Args:
        artifact: Validated replay artifact.

    Returns:
        Closed versioned primitive object without proxy or dataclass values.

    Raises:
        ValueError: If ``artifact`` is not a validated replay artifact.
    """
    if not isinstance(artifact, BamReplayArtifact):
        raise ValueError("BAM replay artifact must be a BamReplayArtifact")
    return {
        "schema_version": BAM_REPLAY_SCHEMA_VERSION,
        "loci": [_encode_locus(locus) for locus in artifact.loci],
    }


def decode_bam_replay_artifact(value: object) -> BamReplayArtifact:
    """Decode and strictly validate canonical primitive replay evidence.

    Args:
        value: Parsed JSON-compatible value.

    Returns:
        Validated replay artifact.

    Raises:
        ValueError: If keys, version, values, ordering, or derived counts differ
            from the closed canonical schema.
    """
    root = _exact_object(value, {"schema_version", "loci"}, "BAM replay artifact")
    if root["schema_version"] != BAM_REPLAY_SCHEMA_VERSION:
        raise ValueError(f"Unsupported BAM replay schema version: {root['schema_version']}")
    raw_loci = root["loci"]
    if not isinstance(raw_loci, list):
        raise ValueError("BAM replay loci must be a list")
    artifact = BamReplayArtifact(tuple(_decode_locus(raw_locus) for raw_locus in raw_loci))
    if encode_bam_replay_artifact(artifact) != value:
        raise ValueError("BAM replay artifact must use canonical ordering and values")
    return artifact


def write_bam_replay_artifact(directory: str | Path, artifact: BamReplayArtifact) -> Path:
    """Atomically write the replay sidecar through a sibling temporary file.

    Args:
        directory: Existing Kestrel run directory.
        artifact: Validated artifact to persist.

    Returns:
        Final sidecar path.

    Raises:
        ValueError: If the directory or artifact is invalid.
        OSError: If writing or atomic installation fails.
    """
    root = Path(directory)
    if not root.is_dir():
        raise ValueError("BAM replay destination must be an existing directory")
    payload = encode_bam_replay_artifact(artifact)
    destination = root / BAM_REPLAY_FILENAME
    temporary = root / f"{BAM_REPLAY_FILENAME}.tmp"
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, destination)
    except OSError:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def read_bam_replay_artifact(directory: str | Path) -> BamReplayArtifact:
    """Read a replay sidecar and fail closed on malformed or duplicate JSON.

    Args:
        directory: Kestrel run directory containing the exact sidecar filename.

    Returns:
        Strictly decoded replay artifact.

    Raises:
        FileNotFoundError: If the artifact is absent.
        OSError: If it cannot be read.
        ValueError: If JSON or the closed schema is malformed.
    """
    path = Path(directory) / BAM_REPLAY_FILENAME
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_object)
    except json.JSONDecodeError as error:
        raise ValueError("BAM replay artifact must contain valid JSON") from error
    return decode_bam_replay_artifact(value)


def invalidate_bam_replay_artifact(directory: str | Path) -> None:
    """Remove only the exact replay sidecar before an identity-aware attempt.

    Args:
        directory: Kestrel run directory that may contain a stale sidecar.
    """
    (Path(directory) / BAM_REPLAY_FILENAME).unlink(missing_ok=True)


def _encode_locus(locus: BamReplayLocus) -> dict[str, object]:
    return {
        "candidate_observation_ordinals": list(locus.candidate_observation_ordinals),
        "state": locus.state,
        "evidence": _encode_evidence(locus.evidence) if locus.evidence is not None else None,
    }


def _encode_evidence(evidence: BamLocusEvidence) -> dict[str, object]:
    return {
        "eligible_record_count": evidence.eligible_record_count,
        "records": [
            {
                "identities": [serialize_molecular_identity(identity) for identity in record.identities],
                "candidate_observation_ordinals": [
                    list(ordinals) for ordinals in record.candidate_observation_ordinals
                ],
                "minimum_kmer_depth": record.minimum_kmer_depth,
            }
            for record in evidence.records
        ],
        "counts": [
            {"identity": serialize_molecular_identity(identity), "record_count": count}
            for identity, count in evidence.counts.items()
        ],
    }


def _decode_locus(value: object) -> BamReplayLocus:
    raw = _exact_object(value, {"candidate_observation_ordinals", "state", "evidence"}, "BAM replay locus")
    raw_ordinals = raw["candidate_observation_ordinals"]
    if not isinstance(raw_ordinals, list):
        raise ValueError("BAM replay candidate observation ordinals must be a list")
    ordinals = tuple(raw_ordinals)
    state = raw["state"]
    if not isinstance(state, str):
        raise ValueError("BAM replay state must be a string")
    evidence = None if raw["evidence"] is None else _decode_evidence(raw["evidence"])
    return BamReplayLocus(cast(tuple[int, ...], ordinals), cast(BamReplayState, state), evidence)


def _decode_evidence(value: object) -> BamLocusEvidence:
    raw = _exact_object(value, {"eligible_record_count", "records", "counts"}, "BAM replay evidence")
    raw_records = raw["records"]
    raw_counts = raw["counts"]
    if not isinstance(raw_records, list):
        raise ValueError("BAM replay records must be a list")
    if not isinstance(raw_counts, list):
        raise ValueError("BAM replay counts must be a list")
    records = tuple(_decode_record(record) for record in raw_records)
    counts = {}
    for raw_count in raw_counts:
        count_entry = _exact_object(raw_count, {"identity", "record_count"}, "BAM replay count")
        identity = _decode_identity(count_entry["identity"])
        if identity in counts:
            raise ValueError("BAM replay counts must contain unique identities")
        counts[identity] = count_entry["record_count"]
    evidence = BamLocusEvidence(records, raw["eligible_record_count"], counts)  # type: ignore[arg-type]
    if _encode_evidence(evidence) != value:
        raise ValueError("BAM replay evidence must use canonical ordering and values")
    return evidence


def _decode_record(value: object) -> BamIdentityEvidence:
    raw = _exact_object(
        value,
        {"identities", "candidate_observation_ordinals", "minimum_kmer_depth"},
        "BAM replay record",
    )
    raw_identities = raw["identities"]
    raw_bindings = raw["candidate_observation_ordinals"]
    if not isinstance(raw_identities, list):
        raise ValueError("BAM replay record identities must be a list")
    if not isinstance(raw_bindings, list) or any(not isinstance(binding, list) for binding in raw_bindings):
        raise ValueError("BAM replay candidate bindings must be lists")
    identities = tuple(_decode_identity(identity) for identity in raw_identities)
    bindings = tuple(tuple(binding) for binding in raw_bindings)
    return BamIdentityEvidence(identities, bindings, raw["minimum_kmer_depth"])  # type: ignore[arg-type]


def _decode_identity(value: object):
    if not isinstance(value, str):
        raise ValueError("BAM replay molecular identities must be strings")
    identity = parse_molecular_identity(value)
    if serialize_molecular_identity(identity) != value:
        raise ValueError("BAM replay molecular identities must use canonical serialization")
    return identity


def _exact_object(value: object, keys: set[str], name: str) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != keys:
        raise ValueError(f"{name} must contain exactly {sorted(keys)}")
    return value


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"BAM replay JSON contains duplicate key: {key}")
        result[key] = value
    return result


def _nonnegative_integer(value: object, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
