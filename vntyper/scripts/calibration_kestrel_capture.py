"""Strict immutable complete-prefilter Kestrel replay captures."""

from __future__ import annotations

import hashlib
import logging
import numbers
import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import NoReturn, cast

import pandas as pd
from Bio.Seq import Seq

from vntyper.scripts.calibration_caller_policy import (
    KESTREL_CALLER_POLICY_POINTERS,
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.identity_candidates import IdentityTranslationComponent
from vntyper.scripts.kestrel_decision_config import KestrelSelection, KestrelSortField

logger = logging.getLogger(__name__)

KESTREL_RAW_COLUMNS: tuple[str, ...] = (
    "Motifs",
    "POS",
    "REF",
    "ALT",
    "Sample",
    "Motif_sequence",
    "Variant",
)
KESTREL_MOTIF_COLUMNS: tuple[str, ...] = ("Motif", "Motif_sequence")

_ROOT_FIELDS = {
    "schema_version",
    "capture_stage",
    "population",
    "comparator_semantics",
    "source_row_count",
    "raw_columns",
    "rows",
    "motif_columns",
    "motif_rows",
    "kestrel_config",
    "capture_policy",
    "baseline_policy",
    "selection",
    "identity_translation",
    "provenance",
}
_RAW_FIELDS = {"source_row_ordinal", *KESTREL_RAW_COLUMNS}
_MOTIF_FIELDS = set(KESTREL_MOTIF_COLUMNS)
_IDENTITY_FIELDS = {
    "kestrel_motifs",
    "advntr_repeat_unit_motifs",
    "advntr_rotation_offset",
    "permit_boundary_insertions",
}
_PROVENANCE_FIELDS = {
    "decision_profile_sha256",
    "kestrel_config_sha256",
    "reference_file_sha256",
    "motif_reference_file_sha256",
    "parsed_motif_table_sha256",
    "identity_table_sha256",
    "kestrel_jar_sha256",
    "capture_policy_sha256",
    "selection_sha256",
}
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_COMPARATORS = {
    "/components/kestrel/alt_filtering/gg_depth_score_threshold": "non-GG or Depth_Score >= threshold",
    "/components/kestrel/confidence_assignment/reporting_floor": "outer gate Depth_Score >= threshold",
    "/components/kestrel/confidence_assignment/var_active_region_threshold": (
        "conditional active_depth <= threshold after high score and open alt-depth gap"
    ),
    "/components/kestrel/confidence_assignment/depth_score_thresholds/low": "closed midband starts at threshold",
    "/components/kestrel/confidence_assignment/depth_score_thresholds/high": (
        "closed midband ends at threshold; promoted tiers require greater score"
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": (
        "high-score low tier includes alt_depth <= threshold"
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": (
        "high tier includes alt_depth >= threshold"
    ),
    "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": (
        "star tier includes alt_depth >= threshold; high tier excludes equality"
    ),
}


@dataclass(frozen=True)
class KestrelRawCandidate:
    """One ordered raw post-VCF candidate, before any scalar decision gate."""

    source_row_ordinal: int
    motifs: str
    variant: str
    position: int
    reference_allele: str
    alternate_allele: str
    sample: str
    motif_sequence: str


@dataclass(frozen=True)
class KestrelMotifRecord:
    """One ordered parsed motif annotation record."""

    motif: str
    motif_sequence: str


@dataclass(frozen=True)
class KestrelCaptureProvenance:
    """Distinct byte-asset commitments and derived parsed-table identities."""

    decision_profile_sha256: str
    kestrel_config_sha256: str
    reference_file_sha256: str
    motif_reference_file_sha256: str
    parsed_motif_table_sha256: str
    identity_table_sha256: str
    kestrel_jar_sha256: str
    capture_policy_sha256: str
    selection_sha256: str


@dataclass(frozen=True)
class KestrelCapture:
    """Immutable complete Kestrel replay input and its canonical identity."""

    rows: tuple[KestrelRawCandidate, ...]
    motifs: tuple[KestrelMotifRecord, ...]
    kestrel_config_json: bytes
    capture_policy_json: bytes
    baseline_policy: CallerPolicyValues
    selection: KestrelSelection
    identity_component: IdentityTranslationComponent
    provenance: KestrelCaptureProvenance
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"Kestrel capture {label} fields differ from the closed contract")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        _fail(f"Kestrel capture {label} must be a lowercase SHA256 digest")
    return value


def _json_object(value: object, label: str) -> tuple[dict[str, object], bytes]:
    if not isinstance(value, Mapping):
        _fail(f"Kestrel capture {label} must be an object")
    try:
        encoded = canonical_json_bytes(dict(value))
        decoded = load_strict_json_object(encoded)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Kestrel capture {label} must be finite JSON data") from error
    return cast(dict[str, object], decoded), encoded


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        _fail(f"Kestrel capture {label} must be a non-empty string")
    return value


def _position(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value < 1:
        _fail("Kestrel capture POS must be a positive integer")
    return int(value)


def _builder_position(value: object) -> int:
    if isinstance(value, str):
        if not value.isascii() or not value.isdecimal() or (len(value) > 1 and value.startswith("0")):
            _fail("Kestrel capture POS must be a canonical positive integer")
        value = int(value)
    return _position(value)


def _sample(value: object) -> str:
    sample = _text(value, "Sample")
    pieces = sample.split(":")
    if (
        len(pieces) != 3
        or not pieces[0]
        or any(
            not item.isascii() or not item.isdecimal() or (len(item) > 1 and item.startswith("0"))
            for item in pieces[1:]
        )
    ):
        _fail("Kestrel capture Sample must contain label:alternate_depth:active_depth with canonical integers")
    return sample


def _raw_row(value: object, expected_ordinal: int) -> KestrelRawCandidate:
    row = _object(value, _RAW_FIELDS, "raw row")
    ordinal = row["source_row_ordinal"]
    if isinstance(ordinal, bool) or not isinstance(ordinal, int) or ordinal != expected_ordinal:
        _fail("Kestrel capture source_row_ordinal must be contiguous and ordered from zero")
    variant = _text(row["Variant"], "Variant")
    if variant not in {"Insertion", "Deletion"}:
        _fail("Kestrel capture Variant must be Insertion or Deletion")
    return KestrelRawCandidate(
        ordinal,
        _text(row["Motifs"], "Motifs"),
        variant,
        _position(row["POS"]),
        _text(row["REF"], "REF"),
        _text(row["ALT"], "ALT"),
        _sample(row["Sample"]),
        _text(row["Motif_sequence"], "Motif_sequence"),
    )


def _motif_row(value: object) -> KestrelMotifRecord:
    row = _object(value, _MOTIF_FIELDS, "motif row")
    return KestrelMotifRecord(_text(row["Motif"], "Motif"), _text(row["Motif_sequence"], "Motif_sequence"))


def _identity_document(component: IdentityTranslationComponent) -> dict[str, object]:
    return {
        "kestrel_motifs": dict(sorted(component.kestrel_motifs.items())),
        "advntr_repeat_unit_motifs": dict(sorted(component.advntr_repeat_unit_motifs.items())),
        "advntr_rotation_offset": component.advntr_rotation_offset,
        "permit_boundary_insertions": component.permit_boundary_insertions,
    }


def _identity(value: object) -> IdentityTranslationComponent:
    raw = _object(value, _IDENTITY_FIELDS, "identity translation")
    first, second = raw["kestrel_motifs"], raw["advntr_repeat_unit_motifs"]
    if not isinstance(first, Mapping) or not isinstance(second, Mapping):
        _fail("Kestrel capture identity tables must be objects")
    if type(raw["permit_boundary_insertions"]) is not bool:
        _fail("Kestrel capture identity permit_boundary_insertions must be boolean")
    return IdentityTranslationComponent(
        cast(Mapping[str, str], first),
        cast(Mapping[str, str], second),
        cast(int, raw["advntr_rotation_offset"]),
        cast(bool, raw["permit_boundary_insertions"]),
    )


def _raw_document(row: KestrelRawCandidate) -> dict[str, object]:
    return {
        "source_row_ordinal": row.source_row_ordinal,
        "Motifs": row.motifs,
        "Variant": row.variant,
        "POS": row.position,
        "REF": row.reference_allele,
        "ALT": row.alternate_allele,
        "Sample": row.sample,
        "Motif_sequence": row.motif_sequence,
    }


def _motif_document(row: KestrelMotifRecord) -> dict[str, object]:
    return {"Motif": row.motif, "Motif_sequence": row.motif_sequence}


def _provenance_document(value: KestrelCaptureProvenance) -> dict[str, object]:
    return {field: getattr(value, field) for field in sorted(_PROVENANCE_FIELDS)}


def _document(capture: KestrelCapture) -> dict[str, object]:
    return {
        "schema_version": "calibration-kestrel-capture-v1",
        "capture_stage": "post-vcf-pre-scoring-v1",
        "population": "complete-prefilter-v1",
        "comparator_semantics": dict(_COMPARATORS),
        "source_row_count": len(capture.rows),
        "raw_columns": list(KESTREL_RAW_COLUMNS),
        "rows": [_raw_document(row) for row in capture.rows],
        "motif_columns": list(KESTREL_MOTIF_COLUMNS),
        "motif_rows": [_motif_document(row) for row in capture.motifs],
        "kestrel_config": load_strict_json_object(capture.kestrel_config_json),
        "capture_policy": load_strict_json_object(capture.capture_policy_json),
        "baseline_policy": caller_policy_values_document(capture.baseline_policy),
        "selection": _selection_document(capture.selection),
        "identity_translation": _identity_document(capture.identity_component),
        "provenance": _provenance_document(capture.provenance),
    }


def _pointer_value(config: Mapping[str, object], pointer: str) -> object:
    value: object = config
    for part in pointer.removeprefix("/").split("/")[2:]:
        if not isinstance(value, Mapping) or part not in value:
            _fail(f"Kestrel capture config is missing baseline policy pointer {pointer}")
        value = value[part]
    return value


def _selection_document(selection: KestrelSelection) -> dict[str, object]:
    return {
        "confidence_priority": dict(selection.confidence_priority),
        "final_filter_columns": list(selection.final_filter_columns),
        "frameshift": {
            "modulus": selection.modulus,
            "insertion_remainder": selection.insertion_remainder,
            "deletion_remainder": selection.deletion_remainder,
        },
        "sort_order": [{"column": field.column, "ascending": field.ascending} for field in selection.sort_order],
        "unflagged_value": selection.unflagged_value,
        "strategy": selection.strategy,
    }


def _selection(value: object) -> KestrelSelection:
    raw = _object(
        value,
        {"confidence_priority", "final_filter_columns", "frameshift", "sort_order", "unflagged_value", "strategy"},
        "selection",
    )
    priorities = raw["confidence_priority"]
    columns = raw["final_filter_columns"]
    frameshift = _object(
        raw["frameshift"], {"modulus", "insertion_remainder", "deletion_remainder"}, "selection frameshift"
    )
    sorting = raw["sort_order"]
    if not isinstance(priorities, Mapping) or not priorities:
        _fail("Kestrel capture selection confidence_priority must be a nonempty object")
    priority: dict[str, int] = {}
    for key, item in priorities.items():
        if not isinstance(key, str) or not key or isinstance(item, bool) or not isinstance(item, int):
            _fail("Kestrel capture selection priorities require text keys and integer values")
        priority[key] = item
    expected_gates = (
        "is_frameshift",
        "is_valid_frameshift",
        "depth_confidence_pass",
        "alt_filter_pass",
        "motif_filter_pass",
        "flag_filter_pass",
    )
    if columns != list(expected_gates):
        _fail("Kestrel capture selection must retain the exact six production gates")
    integers = []
    for field in ("modulus", "insertion_remainder", "deletion_remainder"):
        item = frameshift[field]
        if isinstance(item, bool) or not isinstance(item, int):
            _fail("Kestrel capture selection frameshift values must be integers")
        integers.append(item)
    if integers[0] < 2 or any(not 0 <= item < integers[0] for item in integers[1:]):
        _fail("Kestrel capture selection frameshift remainders must fit the modulus")
    if not isinstance(sorting, list) or not sorting:
        _fail("Kestrel capture selection sort_order must be a nonempty list")
    sort_order = []
    for item in sorting:
        row = _object(item, {"column", "ascending"}, "selection sort row")
        if not isinstance(row["column"], str) or not row["column"] or type(row["ascending"]) is not bool:
            _fail("Kestrel capture selection sort rows require column text and boolean direction")
        sort_order.append(KestrelSortField(cast(str, row["column"]), cast(bool, row["ascending"])))
    unflagged, strategy = raw["unflagged_value"], raw["strategy"]
    if not isinstance(unflagged, str) or not unflagged or strategy not in {"legacy", "identity_dominance"}:
        _fail("Kestrel capture selection unflagged value or strategy is invalid")
    return KestrelSelection(
        MappingProxyType(priority),
        expected_gates,
        integers[0],
        integers[1],
        integers[2],
        tuple(sort_order),
        unflagged,
        cast(str, strategy),
    )


def _validate_baseline(config: Mapping[str, object], policy: CallerPolicyValues) -> None:
    discrete = {
        "/components/kestrel/confidence_assignment/var_active_region_threshold",
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/low",
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low",
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high",
    }
    for pointer in KESTREL_CALLER_POLICY_POINTERS:
        configured = _pointer_value(config, pointer)
        declared = policy.values[pointer]
        wrong_type = (
            isinstance(configured, bool)
            or (pointer in discrete and not isinstance(configured, int))
            or (pointer not in discrete and not isinstance(configured, (int, float)))
        )
        if wrong_type or configured != declared:
            _fail("Kestrel capture config differs from its frozen baseline policy")


def decode_kestrel_capture(value: object) -> KestrelCapture:
    """Decode a strict complete Kestrel replay capture.

    Args:
        value: Closed JSON-compatible capture document.

    Returns:
        Immutable validated replay evidence.

    Raises:
        ValueError: If shape, cells, commitments or derived digests differ.
    """
    raw = _object(value, _ROOT_FIELDS, "root")
    fixed = {
        "schema_version": "calibration-kestrel-capture-v1",
        "capture_stage": "post-vcf-pre-scoring-v1",
        "population": "complete-prefilter-v1",
        "comparator_semantics": _COMPARATORS,
        "raw_columns": list(KESTREL_RAW_COLUMNS),
        "motif_columns": list(KESTREL_MOTIF_COLUMNS),
    }
    if any(type(raw[field]) is not type(expected) or raw[field] != expected for field, expected in fixed.items()):
        _fail("Kestrel capture schema, stage, population, columns or comparator semantics differ")
    count = raw["source_row_count"]
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        _fail("Kestrel capture source_row_count must be a nonnegative integer")
    raw_rows = raw["rows"]
    motif_rows = raw["motif_rows"]
    if not isinstance(raw_rows, list) or not isinstance(motif_rows, list) or len(raw_rows) != count:
        _fail("Kestrel capture rows must agree with source_row_count")
    rows = tuple(_raw_row(row, ordinal) for ordinal, row in enumerate(raw_rows))
    motifs = tuple(_motif_row(row) for row in motif_rows)
    config, config_json = _json_object(raw["kestrel_config"], "kestrel config")
    capture_policy, capture_policy_json = _json_object(raw["capture_policy"], "capture policy")
    baseline = decode_caller_policy_values(raw["baseline_policy"])
    _validate_baseline(config, baseline)
    selection = _selection(raw["selection"])
    identity = _identity(raw["identity_translation"])
    provenance_raw = _object(raw["provenance"], _PROVENANCE_FIELDS, "provenance")
    provenance = KestrelCaptureProvenance(
        **{field: _digest(provenance_raw[field], field) for field in _PROVENANCE_FIELDS}
    )
    motif_document = {"columns": list(KESTREL_MOTIF_COLUMNS), "rows": [_motif_document(row) for row in motifs]}
    if provenance.kestrel_config_sha256 != canonical_sha256(config):
        _fail("Kestrel capture config digest differs from its canonical config")
    if provenance.parsed_motif_table_sha256 != canonical_sha256(motif_document):
        _fail("Kestrel capture parsed motif table digest differs from its rows")
    if provenance.identity_table_sha256 != canonical_sha256(_identity_document(identity)):
        _fail("Kestrel capture identity table digest differs from its rows")
    if provenance.capture_policy_sha256 != canonical_sha256(capture_policy):
        _fail("Kestrel capture policy digest differs from its canonical policy")
    if provenance.selection_sha256 != canonical_sha256(_selection_document(selection)):
        _fail("Kestrel capture selection digest differs from its frozen decisions")
    capture = KestrelCapture(
        rows, motifs, config_json, capture_policy_json, baseline, selection, identity, provenance, ""
    )
    return KestrelCapture(
        rows,
        motifs,
        config_json,
        capture_policy_json,
        baseline,
        selection,
        identity,
        provenance,
        canonical_sha256(_document(capture)),
    )


def build_kestrel_capture(
    raw_frame: pd.DataFrame,
    motif_frame: pd.DataFrame,
    *,
    kestrel_config: Mapping[str, object],
    baseline_policy: CallerPolicyValues,
    selection: KestrelSelection,
    identity_component: IdentityTranslationComponent,
    decision_profile_sha256: str,
    reference_file_bytes: bytes,
    motif_reference_file_bytes: bytes,
    kestrel_jar_bytes: bytes,
    capture_policy: Mapping[str, object],
) -> KestrelCapture:
    """Capture complete raw production inputs and exact asset commitments.

    Args:
        raw_frame: Exact ordered post-VCF, pre-scoring candidate frame.
        motif_frame: Exact ordered parsed motif annotation table.
        kestrel_config: Complete frozen production Kestrel component.
        baseline_policy: Complete caller policy represented by the config.
        selection: Frozen production filtering, ordering and dominance decisions.
        identity_component: Frozen production identity translation tables.
        decision_profile_sha256: External full-profile commitment.
        reference_file_bytes: Exact original Kestrel reference FASTA bytes.
        motif_reference_file_bytes: Exact original motif annotation FASTA bytes.
        kestrel_jar_bytes: Exact invoked Kestrel JAR bytes.
        capture_policy: Recruitment/assembly policy requiring recapture when changed.

    Returns:
        Strict immutable capture suitable for production replay.

    Raises:
        ValueError: If frames, policy, config or commitments are invalid.
    """
    if not isinstance(raw_frame, pd.DataFrame) or set(raw_frame.columns) != set(KESTREL_RAW_COLUMNS):
        _fail("Kestrel capture raw columns differ from the exact pre-scoring contract")
    if not isinstance(motif_frame, pd.DataFrame) or tuple(motif_frame.columns) != KESTREL_MOTIF_COLUMNS:
        _fail("Kestrel capture motif columns differ from the exact parsed-table contract")
    if not isinstance(selection, KestrelSelection):
        _fail("Kestrel capture selection must be a validated KestrelSelection")
    if not isinstance(identity_component, IdentityTranslationComponent):
        _fail("Kestrel capture identity component must be a validated translation component")
    selection_document = _selection_document(_selection(_selection_document(selection)))
    _digest(decision_profile_sha256, "decision_profile_sha256")
    for content, label in (
        (reference_file_bytes, "reference file"),
        (motif_reference_file_bytes, "motif reference file"),
        (kestrel_jar_bytes, "Kestrel JAR"),
    ):
        if not isinstance(content, bytes):
            _fail(f"Kestrel capture {label} content must be bytes")
    config, _ = _json_object(kestrel_config, "kestrel config")
    policy_document, _ = _json_object(capture_policy, "capture policy")
    rows = []
    for ordinal, record in enumerate(raw_frame[list(KESTREL_RAW_COLUMNS)].to_dict("records")):
        record["source_row_ordinal"] = ordinal
        record["POS"] = _builder_position(record["POS"])
        if isinstance(record["Motif_sequence"], Seq):
            record["Motif_sequence"] = str(record["Motif_sequence"])
        _raw_row(record, ordinal)
        rows.append(record)
    motifs = motif_frame.to_dict("records")
    for record in motifs:
        _motif_row(record)
    identity_document = _identity_document(identity_component)
    motif_document = {"columns": list(KESTREL_MOTIF_COLUMNS), "rows": motifs}
    document = {
        "schema_version": "calibration-kestrel-capture-v1",
        "capture_stage": "post-vcf-pre-scoring-v1",
        "population": "complete-prefilter-v1",
        "comparator_semantics": dict(_COMPARATORS),
        "source_row_count": len(rows),
        "raw_columns": list(KESTREL_RAW_COLUMNS),
        "rows": rows,
        "motif_columns": list(KESTREL_MOTIF_COLUMNS),
        "motif_rows": motifs,
        "kestrel_config": config,
        "capture_policy": policy_document,
        "baseline_policy": caller_policy_values_document(baseline_policy),
        "selection": selection_document,
        "identity_translation": identity_document,
        "provenance": {
            "decision_profile_sha256": decision_profile_sha256,
            "kestrel_config_sha256": canonical_sha256(config),
            "reference_file_sha256": hashlib.sha256(reference_file_bytes).hexdigest(),
            "motif_reference_file_sha256": hashlib.sha256(motif_reference_file_bytes).hexdigest(),
            "parsed_motif_table_sha256": canonical_sha256(motif_document),
            "identity_table_sha256": canonical_sha256(identity_document),
            "kestrel_jar_sha256": hashlib.sha256(kestrel_jar_bytes).hexdigest(),
            "capture_policy_sha256": canonical_sha256(policy_document),
            "selection_sha256": canonical_sha256(selection_document),
        },
    }
    return decode_kestrel_capture(document)


def kestrel_capture_document(capture: KestrelCapture) -> dict[str, object]:
    """Project a validated capture as fresh canonical JSON-compatible content.

    Args:
        capture: Previously decoded immutable capture.

    Returns:
        Fresh closed capture document.

    Raises:
        ValueError: If typed content or its canonical digest was forged.
    """
    if not isinstance(capture, KestrelCapture):
        _fail("Kestrel capture projection requires KestrelCapture")
    decoded = decode_kestrel_capture(_document(capture))
    if decoded != capture:
        _fail("Kestrel capture differs from its canonical content or digest")
    return _document(decoded)
