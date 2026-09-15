"""Pure immutable contracts for normalized calibration intake."""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Literal, Protocol, TypeVar, cast

from vntyper.scripts.calibration_contract import CalibrationRole, EvidenceProvenance
from vntyper.scripts.calibration_manifest import GROUP_NAMESPACES
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256

IdentityStatus = Literal["confirmed", "unresolved"]
ArtifactFormat = Literal["BAM", "CRAM", "FASTQ_PAIR"]
InputScope = Literal["full", "regional"]
AliasEvidence = Literal["explicit-user", "source-crosswalk"]
Genotype = Literal["positive", "negative", "unknown"]
TruthStatus = Literal["confirmed", "disputed", "missing"]
LengthUnit = Literal["repeat-count", "bp"]
LengthMeasurement = Literal["exact", "interval", "censored"]

_ROOT_FIELDS = {"schema_version", "specimens", "artifacts", "aliases", "truth", "assignments"}
_SPECIMEN_FIELDS = {"key", "individual_key", "family_key", "identity_status", "previously_examined"}
_ARTIFACT_FIELDS = {
    "key",
    "specimen_key",
    "path",
    "format",
    "mate_path",
    "assembly",
    "assay_class",
    "input_scope",
    "preprocessing_id",
    "replicate_group",
    "expected_sha256",
}
_ALIAS_FIELDS = {"alias", "specimen_key", "evidence"}
_TRUTH_FIELDS = {
    "specimen_key",
    "genotype",
    "variants",
    "length",
    "method",
    "source_digest",
    "source_row",
    "status",
}
_LENGTH_FIELDS = {
    "allele_1",
    "allele_2",
    "unit",
    "repeat_unit_bp",
    "boundary_definition",
    "measurement",
    "lower_bound",
    "upper_bound",
    "conversion_id",
}
_ASSIGNMENT_FIELDS = {"specimen_key", "role", "provenance", "groups"}
_IDENTITY_STATUSES = frozenset({"confirmed", "unresolved"})
_ARTIFACT_FORMATS = frozenset({"BAM", "CRAM", "FASTQ_PAIR"})
_INPUT_SCOPES = frozenset({"full", "regional"})
_ALIAS_EVIDENCE = frozenset({"explicit-user", "source-crosswalk"})
_GENOTYPES = frozenset({"positive", "negative", "unknown"})
_TRUTH_STATUSES = frozenset({"confirmed", "disputed", "missing"})
_LENGTH_UNITS = frozenset({"repeat-count", "bp"})
_MEASUREMENTS = frozenset({"exact", "interval", "censored"})
_ROLES = frozenset({"training", "policy-selection", "validation", "locked-heldout"})
_PROVENANCE = frozenset({"development", "external-custodian"})
_CONFIRMATORY_ROLES = frozenset({"validation", "locked-heldout"})
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")

_DecodedRow = TypeVar("_DecodedRow")


class _SpecimenJoined(Protocol):
    @property
    def specimen_key(self) -> str: ...


@dataclass(frozen=True)
class Specimen:
    """One local opaque biological specimen declaration."""

    key: str
    individual_key: str | None
    family_key: str | None
    identity_status: IdentityStatus
    previously_examined: bool


@dataclass(frozen=True)
class InputArtifact:
    """One declared local input artifact without filesystem resolution."""

    key: str
    specimen_key: str
    path: str
    format: ArtifactFormat
    mate_path: str | None
    assembly: str
    assay_class: str
    input_scope: InputScope
    preprocessing_id: str
    replicate_group: str
    expected_sha256: str | None


@dataclass(frozen=True)
class SpecimenAlias:
    """One explicit source alias joined directly to a specimen key."""

    alias: str
    specimen_key: str
    evidence: AliasEvidence


@dataclass(frozen=True)
class LengthTruth:
    """One preserved paired length observation and its measurement semantics."""

    allele_1: Fraction | None
    allele_2: Fraction | None
    unit: LengthUnit
    repeat_unit_bp: int
    boundary_definition: str
    measurement: LengthMeasurement
    lower_bound: Fraction | None
    upper_bound: Fraction | None
    conversion_id: str | None


@dataclass(frozen=True)
class TruthRecord:
    """Optional genotype and length truth kept independently of eligibility."""

    specimen_key: str
    genotype: Genotype
    variants: tuple[str, ...]
    length: LengthTruth | None
    method: str
    source_digest: str
    source_row: int
    status: TruthStatus


@dataclass(frozen=True)
class Assignment:
    """One role assignment and its predeclared leakage groups."""

    specimen_key: str
    role: CalibrationRole
    provenance: EvidenceProvenance
    groups: Mapping[str, tuple[str, ...]]


@dataclass(frozen=True)
class IntakeDeclaration:
    """Normalized immutable intake with a digest of its canonical encoding."""

    specimens: tuple[Specimen, ...]
    artifacts: tuple[InputArtifact, ...]
    aliases: tuple[SpecimenAlias, ...]
    truth: tuple[TruthRecord, ...]
    assignments: tuple[Assignment, ...]
    sha256: str


def decode_intake(value: object) -> IntakeDeclaration:
    """Decode a strict normalized calibration intake value.

    Duplicate JSON member names must be rejected while parsing with
    :func:`vntyper.scripts.canonical_json.load_strict_json_object` before this
    already-decoded value reaches the typed contract.

    Args:
        value: Parsed calibration intake JSON value.

    Returns:
        The normalized immutable declaration.

    Raises:
        ValueError: If fields, values, joins, or role protections are invalid.
    """
    root = _exact_object(value, _ROOT_FIELDS, "calibration intake")
    if root["schema_version"] != "calibration-intake-v1":
        raise ValueError("calibration intake schema version must be calibration-intake-v1")
    specimens = tuple(sorted(_decode_rows(root["specimens"], _decode_specimen, "specimens"), key=lambda row: row.key))
    artifacts = tuple(sorted(_decode_rows(root["artifacts"], _decode_artifact, "artifacts"), key=lambda row: row.key))
    aliases = tuple(sorted(_decode_rows(root["aliases"], _decode_alias, "aliases"), key=lambda row: row.alias))
    truth = tuple(
        sorted(
            _decode_rows(root["truth"], _decode_truth, "truth"),
            key=lambda row: (row.specimen_key, row.source_digest, row.source_row),
        )
    )
    assignments = tuple(
        sorted(_decode_rows(root["assignments"], _decode_assignment, "assignments"), key=lambda row: row.specimen_key)
    )
    _require_unique((row.key for row in specimens), "calibration specimen keys")
    _require_unique((row.key for row in artifacts), "calibration artifact keys")
    _require_unique((row.alias for row in aliases), "calibration aliases")
    _require_unique((row.specimen_key for row in assignments), "calibration assignment specimen keys")
    _require_unique(
        ((row.specimen_key, row.source_digest, row.source_row) for row in truth),
        "calibration truth source rows",
    )
    specimen_keys = {row.key for row in specimens}
    if not specimen_keys:
        raise ValueError("calibration intake specimens must be a non-empty list")
    if {row.specimen_key for row in assignments} != specimen_keys:
        raise ValueError("every calibration specimen must have exactly one assignment")
    _require_known_specimens(artifacts, specimen_keys, "artifact")
    _require_known_specimens(aliases, specimen_keys, "alias")
    _require_known_specimens(truth, specimen_keys, "truth")
    alias_values = {row.alias for row in aliases}
    if alias_values & specimen_keys:
        raise ValueError("calibration alias cannot shadow a specimen key")
    _validate_role_boundaries(specimens, artifacts, truth, assignments)
    declaration = IntakeDeclaration(specimens, artifacts, aliases, truth, assignments, "")
    return IntakeDeclaration(
        specimens, artifacts, aliases, truth, assignments, canonical_sha256(encode_intake(declaration))
    )


def encode_intake(declaration: IntakeDeclaration) -> dict[str, object]:
    """Project an immutable intake into its normalized JSON document.

    Args:
        declaration: Validated intake declaration.

    Returns:
        A fresh JSON-compatible normalized document.

    Raises:
        ValueError: If ``declaration`` is not an intake declaration.
    """
    if not isinstance(declaration, IntakeDeclaration):
        raise ValueError("calibration intake must be an IntakeDeclaration")
    return {
        "schema_version": "calibration-intake-v1",
        "specimens": [_encode_specimen(row) for row in declaration.specimens],
        "artifacts": [_encode_artifact(row) for row in declaration.artifacts],
        "aliases": [_encode_alias(row) for row in declaration.aliases],
        "truth": [_encode_truth(row) for row in declaration.truth],
        "assignments": [_encode_assignment(row) for row in declaration.assignments],
    }


def canonical_intake_bytes(declaration: IntakeDeclaration) -> bytes:
    """Return the canonical normalized JSON bytes for an intake declaration.

    Args:
        declaration: Validated intake declaration.

    Returns:
        RFC 8785 JSON bytes terminated by one newline.
    """
    return canonical_json_bytes(encode_intake(declaration))


def _decode_specimen(value: object) -> Specimen:
    raw = _exact_object(value, _SPECIMEN_FIELDS, "calibration specimen")
    status = raw["identity_status"]
    if status not in _IDENTITY_STATUSES:
        raise ValueError(f"unsupported calibration specimen identity status: {status!r}")
    previously_examined = raw["previously_examined"]
    if not isinstance(previously_examined, bool):
        raise ValueError("calibration specimen previously examined must be a boolean")
    return Specimen(
        _nonempty_string(raw["key"], "calibration specimen key"),
        _optional_string(raw["individual_key"], "calibration specimen individual key"),
        _optional_string(raw["family_key"], "calibration specimen family key"),
        cast(IdentityStatus, status),
        previously_examined,
    )


def _decode_artifact(value: object) -> InputArtifact:
    raw = _exact_object(value, _ARTIFACT_FIELDS, "calibration artifact")
    artifact_format = raw["format"]
    if artifact_format not in _ARTIFACT_FORMATS:
        raise ValueError(f"unsupported calibration artifact format: {artifact_format!r}")
    scope = raw["input_scope"]
    if scope not in _INPUT_SCOPES:
        raise ValueError(f"unsupported calibration artifact input scope: {scope!r}")
    mate_path = _optional_string(raw["mate_path"], "calibration artifact mate path")
    if artifact_format == "FASTQ_PAIR" and mate_path is None:
        raise ValueError("calibration FASTQ pair artifact requires a mate path")
    if artifact_format != "FASTQ_PAIR" and mate_path is not None:
        raise ValueError("calibration BAM/CRAM artifact forbids a mate path")
    digest = raw["expected_sha256"]
    if digest is not None and (not isinstance(digest, str) or _SHA256.fullmatch(digest) is None):
        raise ValueError("calibration artifact expected sha256 must be a lowercase SHA-256 or null")
    return InputArtifact(
        _nonempty_string(raw["key"], "calibration artifact key"),
        _nonempty_string(raw["specimen_key"], "calibration artifact specimen key"),
        _nonempty_string(raw["path"], "calibration artifact path"),
        cast(ArtifactFormat, artifact_format),
        mate_path,
        _nonempty_string(raw["assembly"], "calibration artifact assembly"),
        _nonempty_string(raw["assay_class"], "calibration artifact assay class"),
        cast(InputScope, scope),
        _nonempty_string(raw["preprocessing_id"], "calibration artifact preprocessing id"),
        _nonempty_string(raw["replicate_group"], "calibration artifact replicate group"),
        cast(str | None, digest),
    )


def _decode_alias(value: object) -> SpecimenAlias:
    raw = _exact_object(value, _ALIAS_FIELDS, "calibration alias")
    evidence = raw["evidence"]
    if evidence not in _ALIAS_EVIDENCE:
        raise ValueError(f"unsupported calibration alias evidence: {evidence!r}")
    return SpecimenAlias(
        _nonempty_string(raw["alias"], "calibration alias"),
        _nonempty_string(raw["specimen_key"], "calibration alias specimen key"),
        cast(AliasEvidence, evidence),
    )


def _decode_truth(value: object) -> TruthRecord:
    raw = _exact_object(value, _TRUTH_FIELDS, "calibration truth")
    genotype = raw["genotype"]
    if genotype not in _GENOTYPES:
        raise ValueError(f"unsupported calibration truth genotype: {genotype!r}")
    status = raw["status"]
    if status not in _TRUTH_STATUSES:
        raise ValueError(f"unsupported calibration truth status: {status!r}")
    digest = raw["source_digest"]
    if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
        raise ValueError("calibration truth source digest must be a lowercase SHA-256")
    length = None if raw["length"] is None else _decode_length(raw["length"])
    if status == "missing" and length is not None:
        raise ValueError("missing calibration truth requires null length")
    return TruthRecord(
        _nonempty_string(raw["specimen_key"], "calibration truth specimen key"),
        cast(Genotype, genotype),
        _string_list(raw["variants"], "calibration truth variants", allow_empty=True),
        length,
        _nonempty_string(raw["method"], "calibration truth method"),
        digest,
        _positive_integer(raw["source_row"], "calibration truth source row"),
        cast(TruthStatus, status),
    )


def _decode_length(value: object) -> LengthTruth:
    raw = _exact_object(value, _LENGTH_FIELDS, "calibration length truth")
    unit = raw["unit"]
    if unit not in _LENGTH_UNITS:
        raise ValueError(f"unsupported calibration length unit: {unit!r}")
    measurement = raw["measurement"]
    if measurement not in _MEASUREMENTS:
        raise ValueError(f"unsupported calibration length measurement: {measurement!r}")
    allele_1 = _optional_positive_number(raw["allele_1"], "calibration length allele 1")
    allele_2 = _optional_positive_number(raw["allele_2"], "calibration length allele 2")
    lower = _optional_positive_number(raw["lower_bound"], "calibration length lower bound")
    upper = _optional_positive_number(raw["upper_bound"], "calibration length upper bound")
    if lower is not None and upper is not None and lower > upper:
        raise ValueError("calibration length lower bound must not exceed upper bound")
    if measurement == "exact" and (lower is not None or upper is not None):
        raise ValueError("exact calibration length truth forbids interval bounds")
    if measurement != "exact" and (allele_1 is not None or allele_2 is not None):
        raise ValueError("interval or censored calibration length truth forbids exact allele values")
    if measurement != "exact" and lower is None and upper is None:
        raise ValueError("interval or censored calibration length truth requires a directional bound")
    if unit == "repeat-count":
        for allele in (allele_1, allele_2):
            if allele is not None and allele.denominator != 1:
                raise ValueError("exact calibration repeat-count truth must be integral")
    conversion_id = _optional_string(raw["conversion_id"], "calibration length conversion id")
    if unit == "bp" and conversion_id is None:
        raise ValueError("calibration bp length truth requires a conversion id")
    return LengthTruth(
        allele_1,
        allele_2,
        cast(LengthUnit, unit),
        _positive_integer(raw["repeat_unit_bp"], "calibration length repeat unit bp"),
        _nonempty_string(raw["boundary_definition"], "calibration length boundary definition"),
        cast(LengthMeasurement, measurement),
        lower,
        upper,
        conversion_id,
    )


def _decode_assignment(value: object) -> Assignment:
    raw = _exact_object(value, _ASSIGNMENT_FIELDS, "calibration assignment")
    role = raw["role"]
    if role not in _ROLES:
        raise ValueError(f"unsupported calibration assignment role: {role!r}")
    provenance = raw["provenance"]
    if provenance not in _PROVENANCE:
        raise ValueError(f"unsupported calibration assignment provenance: {provenance!r}")
    if role == "locked-heldout" and provenance != "external-custodian":
        raise ValueError("locked calibration assignment requires external custodian provenance")
    raw_groups = _exact_object(raw["groups"], set(GROUP_NAMESPACES), "calibration assignment groups")
    groups = {
        namespace: _string_list(raw_groups[namespace], f"calibration group {namespace}", allow_empty=False)
        for namespace in GROUP_NAMESPACES
    }
    return Assignment(
        _nonempty_string(raw["specimen_key"], "calibration assignment specimen key"),
        cast(CalibrationRole, role),
        cast(EvidenceProvenance, provenance),
        MappingProxyType(groups),
    )


def _validate_role_boundaries(
    specimens: tuple[Specimen, ...],
    artifacts: tuple[InputArtifact, ...],
    truth: tuple[TruthRecord, ...],
    assignments: tuple[Assignment, ...],
) -> None:
    specimens_by_key = {row.key: row for row in specimens}
    artifacts_by_specimen = {row.specimen_key for row in artifacts}
    truth_by_specimen = {row.specimen_key for row in truth}
    for assignment in assignments:
        specimen = specimens_by_key[assignment.specimen_key]
        if assignment.role in _CONFIRMATORY_ROLES and specimen.identity_status == "unresolved":
            raise ValueError(f"unresolved identity cannot enter {assignment.role}")
        if assignment.role in _CONFIRMATORY_ROLES and specimen.previously_examined:
            raise ValueError(f"previously examined specimen cannot enter {assignment.role}")
        if assignment.role == "locked-heldout":
            if assignment.specimen_key in artifacts_by_specimen:
                raise ValueError("locked held-out membership forbids artifact paths")
            if assignment.specimen_key in truth_by_specimen:
                raise ValueError("locked held-out membership forbids truth")
        elif assignment.specimen_key not in artifacts_by_specimen:
            raise ValueError("non-locked calibration specimen requires an input artifact")


def _encode_specimen(row: Specimen) -> dict[str, object]:
    return {
        "key": row.key,
        "individual_key": row.individual_key,
        "family_key": row.family_key,
        "identity_status": row.identity_status,
        "previously_examined": row.previously_examined,
    }


def _encode_artifact(row: InputArtifact) -> dict[str, object]:
    return {
        "key": row.key,
        "specimen_key": row.specimen_key,
        "path": row.path,
        "format": row.format,
        "mate_path": row.mate_path,
        "assembly": row.assembly,
        "assay_class": row.assay_class,
        "input_scope": row.input_scope,
        "preprocessing_id": row.preprocessing_id,
        "replicate_group": row.replicate_group,
        "expected_sha256": row.expected_sha256,
    }


def _encode_alias(row: SpecimenAlias) -> dict[str, object]:
    return {"alias": row.alias, "specimen_key": row.specimen_key, "evidence": row.evidence}


def _encode_truth(row: TruthRecord) -> dict[str, object]:
    length = row.length
    encoded_length = None
    if length is not None:
        encoded_length = {
            "allele_1": _json_number(length.allele_1),
            "allele_2": _json_number(length.allele_2),
            "unit": length.unit,
            "repeat_unit_bp": length.repeat_unit_bp,
            "boundary_definition": length.boundary_definition,
            "measurement": length.measurement,
            "lower_bound": _json_number(length.lower_bound),
            "upper_bound": _json_number(length.upper_bound),
            "conversion_id": length.conversion_id,
        }
    return {
        "specimen_key": row.specimen_key,
        "genotype": row.genotype,
        "variants": list(row.variants),
        "length": encoded_length,
        "method": row.method,
        "source_digest": row.source_digest,
        "source_row": row.source_row,
        "status": row.status,
    }


def _encode_assignment(row: Assignment) -> dict[str, object]:
    return {
        "specimen_key": row.specimen_key,
        "role": row.role,
        "provenance": row.provenance,
        "groups": {namespace: list(row.groups[namespace]) for namespace in GROUP_NAMESPACES},
    }


def _decode_rows(
    value: object,
    decoder: Callable[[object], _DecodedRow],
    label: str,
) -> tuple[_DecodedRow, ...]:
    if not isinstance(value, list):
        raise ValueError(f"calibration intake {label} must be a list")
    return tuple(decoder(row) for row in value)


def _exact_object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    if any(not isinstance(key, str) for key in value):
        raise ValueError(f"{label} keys must be strings")
    return value


def _nonempty_string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a non-empty string")
    return value


def _optional_string(value: object, label: str) -> str | None:
    if value is None:
        return None
    return _nonempty_string(value, label)


def _string_list(value: object, label: str, *, allow_empty: bool) -> tuple[str, ...]:
    if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
        raise ValueError(f"{label} must be a string list")
    values = tuple(value)
    if not allow_empty and not values:
        raise ValueError(f"{label} must be non-empty")
    if len(values) != len(set(values)):
        raise ValueError(f"{label} values must be unique")
    return tuple(sorted(values))


def _positive_integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _optional_positive_number(value: object, label: str) -> Fraction | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be a finite positive number or null")
    return Fraction(str(value))


def _require_unique(values: Iterable[object], label: str) -> None:
    sequence = tuple(values)
    if len(sequence) != len(set(sequence)):
        raise ValueError(f"{label} must be unique")


def _require_known_specimens(rows: Iterable[_SpecimenJoined], specimen_keys: set[str], label: str) -> None:
    unknown = sorted({row.specimen_key for row in rows} - specimen_keys)
    if unknown:
        raise ValueError(f"calibration {label} rows reference unknown specimens: {unknown}")


def _json_number(value: Fraction | None) -> int | float | None:
    if value is None:
        return None
    if value.denominator == 1:
        return value.numerator
    return float(value)
