"""Strict read-only extraction from immutable completed pipeline runs."""

from __future__ import annotations

import csv
import hashlib
import io
import math
import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

from vntyper.modules.advntr.artifact_evidence import load_artifact_evidence
from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.identity_candidate_persistence import (
    IDENTITY_CAPTURE_COLUMNS,
    PersistedIdentityCandidate,
    PersistedIdentityCapture,
    parse_candidate_capture_cells,
    parse_selected_candidate_cells,
)
from vntyper.scripts.identity_candidates import LEGACY_GATE_COLUMNS
from vntyper.scripts.molecular_identity import serialize_molecular_identity
from vntyper.scripts.molecular_identity_presentation import IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS
from vntyper.scripts.nomenclature_bam_replay import (
    BamReplayArtifact,
    decode_bam_replay_artifact,
    validate_bam_replay_artifact_ordinals,
)
from vntyper.scripts.profile_provenance import resolve_summary_profile
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_KESTREL

CORE_RUN_ARTIFACTS = (
    "pipeline_summary.json",
    "provenance/decision_profile.json",
    "kestrel/kestrel_pre_result.tsv",
    "kestrel/bam_identity_replay.v1.json",
    "kestrel/kestrel_result.tsv",
)
ADVNTR_RUN_ARTIFACTS = ("provenance/advntr_artifact_evidence.json", "advntr/output_adVNTR_result.tsv")
_SHA256_HEX = frozenset("0123456789abcdef")


@dataclass(frozen=True)
class RunArtifactDeclaration:
    """One precommitted run root and its exact fixed-relative artifact hashes."""

    root: Path
    artifact_sha256: Mapping[str, str]


@dataclass(frozen=True)
class RunExtraction:
    """Derived runtime features and two independently parsed shipped projections."""

    manifest_key: str
    features: Mapping[str, object]
    expected_row: Mapping[str, object]
    observed_row: Mapping[str, object]
    artifact_sha256: Mapping[str, str]


def decode_run_artifact_declaration(value: object) -> RunArtifactDeclaration:
    """Decode one path-free-in-output precommit for a completed run.

    Args:
        value: Closed ``root`` plus structured per-artifact SHA-256 object.
    Returns:
        Validated immutable declaration.
    Raises:
        ValueError: If fields, paths, artifact membership, or hashes are invalid.
    """
    raw = _exact_mapping(value, {"root", "artifacts"}, "calibration run declaration")
    root = raw["root"]
    if not isinstance(root, str) or not root:
        raise ValueError("calibration run root must be a non-empty string")
    artifact_values = raw["artifacts"]
    if not isinstance(artifact_values, Mapping):
        raise ValueError("calibration run artifacts must be a structured SHA-256 object")
    artifact_sha256 = validate_run_artifact_hashes(artifact_values)
    return RunArtifactDeclaration(Path(root), MappingProxyType(artifact_sha256))


def validate_run_artifact_hashes(value: Mapping[object, object]) -> dict[str, str]:
    """Validate the exact fixed-relative artifact hash map for one current run."""
    artifact_sha256: dict[str, str] = {}
    for relative, digest in value.items():
        if not isinstance(relative, str) or not relative:
            raise ValueError("calibration run artifact names must be non-empty strings")
        if Path(relative).is_absolute() or relative != Path(relative).as_posix() or ".." in Path(relative).parts:
            raise ValueError("calibration run artifacts must use fixed safe relative paths")
        if not _is_sha256(digest):
            raise ValueError(f"calibration run artifact SHA-256 is invalid: {relative}")
        assert isinstance(digest, str)
        artifact_sha256[relative] = digest
    allowed = {frozenset(CORE_RUN_ARTIFACTS), frozenset((*CORE_RUN_ARTIFACTS, *ADVNTR_RUN_ARTIFACTS))}
    if frozenset(artifact_sha256) not in allowed:
        raise ValueError("calibration run artifact membership differs from the closed current-run contract")
    return dict(sorted(artifact_sha256.items()))


def decode_run_hashes(value: Mapping[str, object]) -> dict[str, dict[str, str]]:
    """Decode member keys mapped to exact structured current-run artifact hashes."""
    hashes: dict[str, dict[str, str]] = {}
    for key, raw_artifacts in value.items():
        if not isinstance(key, str) or not key or not isinstance(raw_artifacts, Mapping):
            raise ValueError("calibration run hashes must map member keys to structured artifact hashes")
        hashes[key] = validate_run_artifact_hashes(raw_artifacts)
    return hashes


def extract_completed_run(
    manifest_key: str,
    assay_class: str,
    declaration: RunArtifactDeclaration,
) -> RunExtraction:
    """Verify and derive one member without opening alignments or rerunning callers.

    Args:
        manifest_key: Opaque partition identity bound to this declaration.
        assay_class: Canonical study-member assay dimension.
        declaration: Precommitted fixed artifact map and run root.
    Returns:
        Runtime-only feature row, independent baseline projections, and exact hashes.
    Raises:
        ValueError: If any artifact is absent, stale, malformed, or cross-inconsistent.
    """
    if not isinstance(manifest_key, str) or not manifest_key:
        raise ValueError("calibration run manifest key must be a non-empty string")
    if not isinstance(assay_class, str) or not assay_class:
        raise ValueError("calibration run assay class must be a non-empty string")
    if not isinstance(declaration, RunArtifactDeclaration):
        raise ValueError("calibration run extraction requires a RunArtifactDeclaration")
    artifacts = _read_precommitted_artifacts(declaration)
    summary = load_strict_json_object(artifacts["pipeline_summary.json"])
    profile = resolve_summary_profile(summary, declaration.root)

    steps = _caller_steps(summary)
    has_advntr = STEP_ADVNTR in steps
    expected_artifacts = (
        frozenset((*CORE_RUN_ARTIFACTS, *ADVNTR_RUN_ARTIFACTS)) if has_advntr else frozenset(CORE_RUN_ARTIFACTS)
    )
    if frozenset(artifacts) != expected_artifacts:
        raise ValueError("calibration run artifact membership differs from its verified caller stages")
    evidence_digest = summary.get("advntr_evidence_digest")
    if has_advntr:
        if not isinstance(evidence_digest, str):
            raise ValueError("calibration adVNTR stage requires a governed evidence digest")
        load_artifact_evidence(
            declaration.root / ADVNTR_RUN_ARTIFACTS[0],
            expected_digest=evidence_digest,
        )
    elif evidence_digest is not None:
        raise ValueError("calibration run without an adVNTR step cannot record governed adVNTR evidence")

    expected_kestrel = _summary_tsv(steps[STEP_KESTREL], artifacts["kestrel/kestrel_result.tsv"], STEP_KESTREL)
    observed_kestrel = _parse_tsv(artifacts["kestrel/kestrel_result.tsv"], "Kestrel final result")
    if expected_kestrel != observed_kestrel:
        raise ValueError("calibration summary parsed result differs from the fixed final Kestrel result")
    expected_rows = _rows(expected_kestrel, "summary Kestrel result", allow_empty=True)
    observed_rows = _rows(observed_kestrel, "fixed Kestrel result", allow_empty=True)
    if len(expected_rows) > 1 or len(observed_rows) > 1:
        raise ValueError("calibration shipped Kestrel projection permits at most one final result row")

    pre_result = _parse_tsv(artifacts["kestrel/kestrel_pre_result.tsv"], "Kestrel pre-result")
    pre_rows = _rows(pre_result, "Kestrel pre-result", allow_empty=True)
    pre_by_ordinal = _validate_pre_result(pre_rows)
    expected_persisted: PersistedIdentityCandidate | None = None
    observed_persisted: PersistedIdentityCandidate | None = None
    if expected_rows:
        expected_persisted = parse_selected_candidate_cells(expected_rows[0])
        observed_persisted = parse_selected_candidate_cells(observed_rows[0])
        if expected_persisted != observed_persisted:
            raise ValueError("calibration summary identity projection differs from the final Kestrel result")
        selected_ordinal = observed_persisted.selected_observation_ordinal
        if selected_ordinal not in pre_by_ordinal:
            raise ValueError("calibration final Kestrel selection is absent from the complete pre-result")
        _crosscheck_selected_pre_result(observed_rows[0], pre_by_ordinal[selected_ordinal].row)

    replay_value = load_strict_json_object(artifacts["kestrel/bam_identity_replay.v1.json"])
    replay = decode_bam_replay_artifact(replay_value)
    validate_bam_replay_artifact_ordinals(replay, tuple(pre_by_ordinal))
    _crosscheck_bam_identity_bindings(replay, pre_by_ordinal)
    bam_features = _bam_features(replay)

    expected_advntr: dict[str, object] | None = None
    observed_advntr: dict[str, object] | None = None
    if has_advntr:
        advntr_bytes = artifacts["advntr/output_adVNTR_result.tsv"]
        expected_advntr = _summary_tsv(steps[STEP_ADVNTR], advntr_bytes, STEP_ADVNTR)
        observed_advntr = _parse_tsv(advntr_bytes, "adVNTR final result")
        if expected_advntr != observed_advntr:
            raise ValueError("calibration summary parsed result differs from the fixed final adVNTR result")

    expected_row = _baseline_row(manifest_key, expected_rows[0] if expected_rows else None)
    observed_row = _baseline_row(manifest_key, observed_rows[0] if observed_rows else None)
    features = _runtime_features(
        assay_class,
        summary,
        profile.sha256,
        pre_by_ordinal,
        bam_features,
        observed_advntr,
    )
    return RunExtraction(
        manifest_key,
        MappingProxyType(features),
        MappingProxyType(expected_row),
        MappingProxyType(observed_row),
        declaration.artifact_sha256,
    )


def _read_precommitted_artifacts(declaration: RunArtifactDeclaration) -> dict[str, bytes]:
    artifacts: dict[str, bytes] = {}
    for relative, expected in declaration.artifact_sha256.items():
        path = declaration.root / relative
        try:
            raw = path.read_bytes()
        except OSError as error:
            raise ValueError(f"calibration replay artifact is missing or unreadable: {path}") from error
        actual = hashlib.sha256(raw).hexdigest()
        if actual != expected:
            raise ValueError(
                f"calibration run artifact SHA-256 differs for {relative}: expected {expected}, got {actual}"
            )
        artifacts[relative] = raw
    return artifacts


def _caller_steps(summary: Mapping[str, object]) -> dict[str, Mapping[str, object]]:
    raw_steps = summary.get("steps")
    if not isinstance(raw_steps, list):
        raise ValueError("calibration schema-3 summary steps must be an array")
    steps: dict[str, Mapping[str, object]] = {}
    for raw_step in raw_steps:
        if not isinstance(raw_step, Mapping):
            raise ValueError("calibration schema-3 summary step must be an object")
        name = raw_step.get("step")
        if name not in {STEP_KESTREL, STEP_ADVNTR}:
            continue
        assert isinstance(name, str)
        if name in steps:
            raise ValueError(f"calibration schema-3 summary contains duplicate {name} steps")
        steps[name] = raw_step
    if STEP_KESTREL not in steps:
        raise ValueError("calibration schema-3 summary requires the Kestrel caller stage")
    return steps


def _summary_tsv(step: Mapping[str, object], final_bytes: bytes, name: str) -> dict[str, object]:
    if step.get("file_type") != "tsv":
        raise ValueError(f"calibration {name} summary step must record a TSV result")
    expected_md5 = step.get("md5sum")
    if not isinstance(expected_md5, str) or expected_md5 != hashlib.md5(final_bytes).hexdigest():
        raise ValueError(f"calibration {name} summary checksum differs from its fixed final result")
    parsed = step.get("parsed_result")
    if not isinstance(parsed, Mapping) or set(parsed) != {"comments", "data"}:
        raise ValueError(f"calibration {name} summary parsed result violates the closed TSV projection")
    comments = parsed["comments"]
    data = parsed["data"]
    if not isinstance(comments, list) or any(not isinstance(item, str) for item in comments):
        raise ValueError(f"calibration {name} summary comments must be strings")
    if not isinstance(data, list) or any(not isinstance(item, Mapping) for item in data):
        raise ValueError(f"calibration {name} summary data must contain row objects")
    return {"comments": list(comments), "data": [dict(item) for item in data]}


def _parse_tsv(raw: bytes, label: str) -> dict[str, object]:
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"calibration {label} must be UTF-8") from error
    reader = csv.reader(io.StringIO(text, newline=""), delimiter="\t")
    comments: list[str] = []
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    for values in reader:
        if not values or not any(value.strip() for value in values):
            continue
        if values[0].startswith("#"):
            comments.append("\t".join(values).lstrip("#").strip())
            continue
        if header is None:
            if any(not value for value in values) or len(values) != len(set(values)):
                raise ValueError(f"calibration {label} header must be non-empty and unique")
            header = values
            continue
        if len(values) != len(header):
            raise ValueError(f"calibration {label} contains a malformed TSV row")
        rows.append(dict(zip(header, values, strict=True)))
    if header is None:
        raise ValueError(f"calibration {label} is empty")
    return {"comments": comments, "data": rows}


def _rows(document: Mapping[str, object], label: str, *, allow_empty: bool = False) -> list[Mapping[str, object]]:
    rows = document.get("data")
    if (
        not isinstance(rows, list)
        or (not rows and not allow_empty)
        or any(not isinstance(row, Mapping) for row in rows)
    ):
        raise ValueError(f"calibration {label} requires non-empty row data")
    return rows


@dataclass(frozen=True)
class _PreResultRow:
    capture: PersistedIdentityCapture
    row: Mapping[str, object]


def _validate_pre_result(rows: Sequence[Mapping[str, object]]) -> dict[int, _PreResultRow]:
    required = {*IDENTITY_CAPTURE_COLUMNS, *IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS}
    ordinals: dict[int, _PreResultRow] = {}
    for row in rows:
        if not required <= set(row):
            raise ValueError("calibration Kestrel pre-result lacks complete identity capture artifacts")
        capture = parse_candidate_capture_cells(row)
        if capture.row_key.source != "kestrel":
            raise ValueError("calibration Kestrel pre-result raw key must record the Kestrel source")
        ordinal = capture.observation_ordinal
        if ordinal in ordinals:
            raise ValueError("calibration Kestrel pre-result observation ordinals must be unique")
        translation = capture.translation
        diagnostic = {
            "Molecular_Identity": (
                serialize_molecular_identity(translation.identity) if translation.identity is not None else ""
            ),
            "Molecular_Identity_Translation_Status": translation.status,
            "Molecular_Identity_Translation_Failure": translation.failure or "",
            "Molecular_Identity_Context_Diverges": "True" if translation.context_diverges else "False",
        }
        if any(row[column] != diagnostic[column] for column in IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS):
            raise ValueError("calibration Kestrel pre-result identity diagnostics differ from captured identity")
        ordinals[ordinal] = _PreResultRow(capture, row)
    if tuple(sorted(ordinals)) != tuple(range(len(ordinals))):
        raise ValueError("calibration Kestrel pre-result observation ordinals must be contiguous from zero")
    return ordinals


def _crosscheck_selected_pre_result(final: Mapping[str, object], pre: Mapping[str, object]) -> None:
    for column in IDENTITY_CAPTURE_COLUMNS:
        if final.get(column) != pre.get(column):
            raise ValueError("calibration final Kestrel selection differs from its pre-result identity capture")


def _crosscheck_bam_identity_bindings(replay: BamReplayArtifact, pre_by_ordinal: Mapping[int, _PreResultRow]) -> None:
    for locus in replay.loci:
        if locus.evidence is None:
            continue
        for record in locus.evidence.records:
            for identity, ordinals in zip(record.identities, record.candidate_observation_ordinals, strict=True):
                if any(pre_by_ordinal[ordinal].capture.translation.identity != identity for ordinal in ordinals):
                    raise ValueError("calibration BAM identity binding differs from its captured pre-result identity")


def _bam_features(replay: BamReplayArtifact) -> dict[str, object]:
    counts: dict[str, int] = {}
    eligible = 0
    depths: dict[str, list[int]] = {}
    for locus in replay.loci:
        if locus.state != "observed" or locus.evidence is None:
            continue
        eligible += locus.evidence.eligible_record_count
        for identity, count in locus.evidence.counts.items():
            serialized = serialize_molecular_identity(identity)
            counts[serialized] = counts.get(serialized, 0) + count
        for record in locus.evidence.records:
            if record.minimum_kmer_depth is None:
                continue
            for identity in record.identities:
                depths.setdefault(serialize_molecular_identity(identity), []).append(record.minimum_kmer_depth)
    if not counts or eligible == 0:
        return {}
    top_count = max(counts.values())
    winners = sorted(identity for identity, count in counts.items() if count == top_count)
    selected = winners[0]
    runner_up = max((count for identity, count in counts.items() if identity != selected), default=0)
    selected_depths = sorted(depths.get(selected, []))
    features: dict[str, object] = {
        "canonical_identity": selected,
        "haplotype_record_count": top_count,
        "haplotype_record_count_margin": top_count - runner_up,
        "haplotype_record_share": top_count / eligible,
        "haplotype_record_share_margin": (top_count - runner_up) / eligible,
        "haplotype_record_tie": len(winners) > 1,
        "xd_availability_count": len(selected_depths),
        "xd_availability_fraction": len(selected_depths) / top_count,
    }
    if selected_depths:
        features["xd_median"] = float(statistics.median(selected_depths))
        features["xd_interquartile_range"] = float(_interquartile_range(selected_depths))
    return features


def _interquartile_range(values: Sequence[int]) -> float:
    if len(values) < 2:
        return 0.0
    midpoint = len(values) // 2
    lower = values[:midpoint]
    upper = values[-midpoint:]
    return float(statistics.median(upper) - statistics.median(lower))


def _runtime_features(
    assay_class: str,
    summary: Mapping[str, object],
    profile_sha256: str | None,
    pre_by_ordinal: Mapping[int, _PreResultRow],
    bam_features: Mapping[str, object],
    advntr: Mapping[str, object] | None,
) -> dict[str, object]:
    features: dict[str, object] = {
        "assay_class": assay_class,
        "advntr_evidence_disposition": _advntr_disposition(advntr),
    }
    features.update(bam_features)
    resolved = {
        serialize_molecular_identity(candidate.capture.translation.identity)
        for candidate in pre_by_ordinal.values()
        if candidate.capture.translation.identity is not None
    }
    selected_identity = features.get("canonical_identity")
    if selected_identity is None and len(resolved) == 1:
        selected_identity = next(iter(resolved))
        features["canonical_identity"] = selected_identity
    selected_rows = [
        candidate.row
        for candidate in pre_by_ordinal.values()
        if candidate.capture.translation.identity is not None
        and serialize_molecular_identity(candidate.capture.translation.identity) == selected_identity
    ]
    if selected_rows:
        row = max(
            selected_rows,
            key=lambda item: _nonnegative_number(
                item.get("Estimated_Depth_AlternateVariant"), "Kestrel alternate k-mer-path depth"
            ),
        )
        features.update(_kestrel_features(row, len(selected_rows), len(resolved)))
    features.update(_advntr_support_features(advntr, selected_identity))
    version = summary.get("version")
    if isinstance(version, str) and version:
        features["tool_version"] = version
    reference = summary.get("reference_assembly_requested")
    if isinstance(reference, str) and reference:
        features["reference_version"] = reference
    if isinstance(profile_sha256, str):
        features["decision_profile_sha256"] = profile_sha256
    policy = summary.get("decision_policy")
    if isinstance(policy, str) and policy:
        features["decision_policy"] = policy
    return features


def _kestrel_features(row: Mapping[str, object], representation_count: int, identity_count: int) -> dict[str, object]:
    gate_values = tuple(row.get(column) for column in LEGACY_GATE_COLUMNS)
    if any(value not in {"True", "False"} for value in gate_values):
        raise ValueError("calibration Kestrel pre-result structural gates must use canonical booleans")
    return {
        "motif_context": _nonempty(row.get("Motif"), "Kestrel motif context"),
        "pair_context": _nonempty(row.get("Motifs"), "Kestrel pair context"),
        "context_local_representation_count": representation_count,
        "alternate_kmer_path_depth": _nonnegative_number(
            row.get("Estimated_Depth_AlternateVariant"), "Kestrel alternate k-mer-path depth"
        ),
        "active_region_depth": _nonnegative_number(
            row.get("Estimated_Depth_Variant_ActiveRegion"), "Kestrel active-region depth"
        ),
        "depth_score": _nonnegative_number(row.get("Depth_Score"), "Kestrel depth score"),
        "structural_gate": all(value == "True" for value in gate_values),
        "cooccurring_identity_count": identity_count,
    }


def _advntr_disposition(document: Mapping[str, object] | None) -> str:
    if document is None:
        return "admissible"
    rows = document.get("data")
    if not isinstance(rows, list):
        raise ValueError("calibration adVNTR result data must be an array")
    dispositions = {
        row.get("Evidence_Disposition") for row in rows if isinstance(row, Mapping) and row.get("VID") != "Negative"
    }
    if not dispositions:
        return "admissible"
    if not dispositions <= {"admissible", "identity-insufficient"}:
        raise ValueError("calibration adVNTR result contains an unsupported governed disposition")
    return "identity-insufficient" if "identity-insufficient" in dispositions else "admissible"


def _advntr_support_features(document: Mapping[str, object] | None, selected_identity: object) -> dict[str, object]:
    if document is None or not isinstance(selected_identity, str):
        return {}
    rows = document.get("data")
    if not isinstance(rows, list):
        raise ValueError("calibration adVNTR result data must be an array")
    matching = [
        row
        for row in rows
        if isinstance(row, Mapping)
        and row.get("VID") != "Negative"
        and row.get("Molecular_Identity") == selected_identity
    ]
    if len(matching) > 1:
        raise ValueError("calibration adVNTR result has multiple support rows for the selected identity")
    if not matching or matching[0].get("Evidence_Disposition") != "admissible":
        return {}
    row = matching[0]
    return {
        "advntr_sequencing_read_support": _nonnegative_number(
            row.get("NumberOfSupportingReads"), "adVNTR sequencing-read support"
        ),
        "advntr_p_value": _nonnegative_number(row.get("Pvalue"), "adVNTR p-value"),
        "advntr_coverage": _nonnegative_number(row.get("Coverage"), "adVNTR coverage"),
    }


def _baseline_row(manifest_key: str, row: Mapping[str, object] | None) -> dict[str, object]:
    if row is None:
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
        }
    persisted = parse_selected_candidate_cells(row)
    identity = persisted.translation.identity
    return {
        "manifest_key": manifest_key,
        "order": 0,
        "canonical_identity": serialize_molecular_identity(identity) if identity is not None else None,
        "name": _optional_nonempty(row.get("Nomenclature")),
        "confidence": _optional_nonempty(row.get("Confidence")),
        "flag": _optional_nonempty(row.get("Flag")),
        "tier": _optional_nonempty(row.get("Nomenclature_Tier")),
        "support": _nonnegative_number(row.get("Estimated_Depth_AlternateVariant"), "Kestrel baseline support"),
        "tie": False,
        "abstention": None,
    }


def _nonnegative_number(value: object, label: str) -> int | float:
    if isinstance(value, bool):
        raise ValueError(f"calibration {label} must be numeric")
    if isinstance(value, (int, float)):
        parsed = value
    elif isinstance(value, str):
        try:
            parsed = float(value) if any(character in value for character in ".eE") else int(value)
        except ValueError as error:
            raise ValueError(f"calibration {label} must be numeric") from error
    else:
        raise ValueError(f"calibration {label} must be numeric")
    if not math.isfinite(parsed) or parsed < 0:
        raise ValueError(f"calibration {label} must be finite and non-negative")
    return parsed


def _nonempty(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"calibration {label} must be a non-empty string")
    return value


def _optional_nonempty(value: object) -> str | None:
    return value if isinstance(value, str) and value else None


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and set(value) <= _SHA256_HEX


def _exact_mapping(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    return value


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"calibration {label} must be an object")
    return value
