"""Strict caller truth and native result projections for target calibration."""

from __future__ import annotations

import hashlib
import logging
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, NoReturn, cast

from vntyper.modules.advntr.advntr_calibration_policy import advntr_canonical_sha256
from vntyper.modules.advntr.advntr_capture import _CAPTURE_FIELDS
from vntyper.modules.advntr.advntr_variant_annotations import derive_ru_and_pos
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_caller_roster import EligibleCallerMember
from vntyper.scripts.calibration_callers import CallerEvidenceRow, EvidenceDisposition
from vntyper.scripts.calibration_kestrel_capture import KestrelCapture
from vntyper.scripts.calibration_kestrel_replay import kestrel_replay_selected_frame, replay_kestrel_capture
from vntyper.scripts.calibration_run_extraction import _parse_tsv, _rows
from vntyper.scripts.calibration_run_projection import build_shipped_projection
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import load_packaged_decision_profile
from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells
from vntyper.scripts.identity_reconciliation import (
    IdentityReconciliationObservation,
    IdentityReconciliationPolicy,
    reconcile_identity_observations,
)
from vntyper.scripts.molecular_identity import AdvntrRepresentation, IdentityTranslation, serialize_molecular_identity
from vntyper.scripts.molecular_identity import EvidenceDisposition as IdentityDisposition
from vntyper.scripts.nomenclature import from_advntr, from_kestrel, nomenclature_config
from vntyper.scripts.nomenclature_decision_config import decision_config_from_component

logger = logging.getLogger(__name__)

TruthGenotype = Literal["positive", "negative", "unknown"]
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_ROOT_FIELDS = {"schema_version", "rows"}
_ROW_FIELDS = {"key", "genotype", "variants"}


@dataclass(frozen=True)
class CallerTruthRow:
    """One independently established genotype and normalized identity set."""

    key: str
    genotype: TruthGenotype
    variants: tuple[str, ...] | None


@dataclass(frozen=True)
class CallerTruth:
    """Exact sealed role truth in its predeclared specimen order."""

    rows: tuple[CallerTruthRow, ...]
    by_key: Mapping[str, CallerTruthRow]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        _fail(f"caller {label} must be nonempty trimmed text")
    return value


def _variants(value: object, genotype: TruthGenotype) -> tuple[str, ...] | None:
    if value is None:
        if genotype == "negative":
            _fail("caller negative truth requires an empty variant list")
        return None
    if not isinstance(value, list):
        _fail("caller truth variants must be null or a sorted unique list")
    values = tuple(_text(item, "truth variant") for item in value)
    if values != tuple(sorted(set(values))):
        _fail("caller truth variants must be sorted and unique")
    if genotype == "positive" and not values:
        _fail("caller positive truth requires identities or explicit missing identity")
    if genotype == "unknown":
        _fail("caller unknown truth cannot provide variant identities")
    if genotype == "negative" and values:
        _fail("caller negative truth requires an empty variant list")
    return values


def decode_caller_truth(value: object, expected_keys: tuple[str, ...]) -> CallerTruth:
    """Decode one role's sealed caller truth against its exact eligible roster."""
    if not isinstance(value, Mapping) or set(value) != _ROOT_FIELDS:
        _fail("caller truth fields differ from the closed contract")
    if value["schema_version"] != "calibration-caller-truth-v1":
        _fail("caller truth schema is unsupported")
    if not isinstance(expected_keys, tuple) or not expected_keys or expected_keys != tuple(sorted(set(expected_keys))):
        _fail("caller truth expected keys must be a sorted unique tuple")
    raw_rows = value["rows"]
    if not isinstance(raw_rows, list):
        _fail("caller truth rows must be a list")
    rows: list[CallerTruthRow] = []
    for value_row in raw_rows:
        if not isinstance(value_row, Mapping) or set(value_row) != _ROW_FIELDS:
            _fail("caller truth row fields differ from the closed contract")
        genotype = value_row["genotype"]
        if not isinstance(genotype, str) or genotype not in {"positive", "negative", "unknown"}:
            _fail("caller truth genotype must be positive, negative, or unknown")
        rows.append(
            CallerTruthRow(
                _text(value_row["key"], "truth key"),
                cast(TruthGenotype, genotype),
                _variants(value_row["variants"], cast(TruthGenotype, genotype)),
            )
        )
    keys = tuple(row.key for row in rows)
    if len(keys) != len(set(keys)):
        _fail("caller truth contains a duplicate key")
    if keys != expected_keys:
        _fail("caller truth does not match the exact eligible roster")
    document = {"schema_version": "calibration-caller-truth-v1", "rows": [_truth_row(row) for row in rows]}
    return CallerTruth(tuple(rows), MappingProxyType({row.key: row for row in rows}), canonical_sha256(document))


def _truth_row(row: CallerTruthRow) -> dict[str, object]:
    return {
        "key": row.key,
        "genotype": row.genotype,
        "variants": None if row.variants is None else list(row.variants),
    }


def caller_truth_document(truth: CallerTruth) -> dict[str, object]:
    """Project caller truth after revalidating typed content and canonical digest."""
    if not isinstance(truth, CallerTruth) or not isinstance(truth.rows, tuple):
        _fail("caller truth projection requires immutable CallerTruth")
    document: dict[str, object] = {
        "schema_version": "calibration-caller-truth-v1",
        "rows": [_truth_row(row) for row in truth.rows],
    }
    if decode_caller_truth(document, tuple(row.key for row in truth.rows)) != truth:
        _fail("caller truth differs from its canonical content or digest")
    return document


def _truth_values(row: CallerTruthRow) -> tuple[bool | None, tuple[str, ...] | None]:
    if row.genotype == "positive":
        return True, row.variants
    if row.genotype == "negative":
        return False, ()
    return None, None


def native_caller_observation(
    member: EligibleCallerMember,
    truth: CallerTruth,
    *,
    kestrel_tsv: bytes,
    advntr_tsv: bytes | None,
    disposition: EvidenceDisposition,
    source_evidence_sha256: str,
) -> CallerEvidenceRow:
    """Parse raw native TSV bytes and bind one final production outcome to truth."""
    caller_truth_document(truth)
    if not isinstance(member, EligibleCallerMember) or member.key not in truth.by_key:
        _fail("native caller result member is absent from sealed truth")
    if not isinstance(source_evidence_sha256, str) or _SHA256.fullmatch(source_evidence_sha256) is None:
        _fail("native caller source evidence requires a lowercase SHA256 digest")
    if disposition not in {"called", "negative", "no-call", "zero-candidate", "unsupported"}:
        _fail("native caller evidence disposition is unsupported")
    kestrel = _rows(_parse_tsv(kestrel_tsv, "native Kestrel result"), "native Kestrel result", allow_empty=True)
    advntr = (
        []
        if advntr_tsv is None
        else _rows(_parse_tsv(advntr_tsv, "native adVNTR result"), "native adVNTR result", allow_empty=True)
    )
    try:
        projection = build_shipped_projection(member.key, kestrel, advntr)
    except (KeyError, TypeError):
        _fail("native caller result cannot be projected to a validated identity")
    identity = projection["canonical_identity"]
    tier = projection["tier"]
    if identity is not None and (not isinstance(identity, str) or not identity):
        _fail("native caller projection identity is invalid")
    called = identity is not None
    if called != (disposition == "called"):
        _fail("native caller projection differs from its capture disposition")
    called_positive: bool | None
    if disposition == "called":
        called_positive = True
    elif disposition in {"negative", "zero-candidate"}:
        called_positive = False
    else:
        called_positive = None
    variants = (cast(str, identity),) if called else ()
    truth_positive, truth_variants = _truth_values(truth.by_key[member.key])
    observation = CallerObservation(
        member.key,
        member.group_key,
        truth_positive,
        truth_variants,
        called_positive,
        variants,
        variants if tier == "A" else (),
    )
    return CallerEvidenceRow(observation, disposition, source_evidence_sha256)


def _known_names() -> frozenset[str]:
    values = nomenclature_config.get("known_variants")
    if not isinstance(values, Mapping):
        _fail("caller replay nomenclature known-variant authority is invalid")
    return frozenset(str(key) for key in values)


def _reconciliation_policy() -> IdentityReconciliationPolicy:
    value = load_packaged_decision_profile().components.get("nomenclature")
    if not isinstance(value, Mapping):
        _fail("caller replay identity reconciliation policy is invalid")
    return decision_config_from_component(value).identity_reconciliation


def _flag_set(value: object) -> frozenset[str]:
    if value is None or value in {"", "Not flagged", "Not applicable"}:
        return frozenset()
    if not isinstance(value, str):
        _fail("caller replay Kestrel flags must be text")
    return frozenset(part.strip() for part in value.split(",") if part.strip())


def _kestrel_reconciliation(
    capture: KestrelCapture, selected: Mapping[str, object]
) -> IdentityReconciliationObservation:
    persisted = parse_selected_candidate_cells(selected)
    motifs, position, reference, alternate = (
        selected.get("Motifs"),
        selected.get("POS"),
        selected.get("REF"),
        selected.get("ALT"),
    )
    if (
        not isinstance(motifs, str)
        or isinstance(position, bool)
        or not isinstance(position, int)
        or not isinstance(reference, str)
        or not isinstance(alternate, str)
    ):
        _fail("caller replay Kestrel selected representation is invalid")
    call = from_kestrel(motifs, position, reference, alternate)
    support = selected.get("Estimated_Depth_AlternateVariant")
    if isinstance(support, bool) or not isinstance(support, int) or support < 0:
        _fail("caller replay Kestrel support must be a nonnegative integer")
    return IdentityReconciliationObservation(
        persisted.translation,
        call.name,
        "kestrel_vcf",
        call.event,
        call.net_length,
        _flag_set(selected.get("Flag")) | persisted.flags,
        IdentityDisposition("admissible"),
        call.name in _known_names() if call.name is not None else False,
        kestrel_alternate_kmer_path_depth=support,
        blocking_gates=persisted.blocking_gates,
        presentation_call_index=0,
    )


def _number(value: object, label: str) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(f"caller replay adVNTR {label} must be numeric")
    if isinstance(value, float) and not math.isfinite(value):
        _fail(f"caller replay adVNTR {label} must be finite")
    return value


def _advntr_call(value: object) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != {"state", "read_support", "mean_coverage", "pvalue"}:
        _fail("caller replay adVNTR call fields differ from the closed contract")
    state = _text(value["state"], "replay adVNTR state")
    support = _number(value["read_support"], "read support")
    coverage = _number(value["mean_coverage"], "mean coverage")
    pvalue = _number(value["pvalue"], "pvalue")
    if not isinstance(support, int) or support < 0 or coverage < 0 or not 0 <= pvalue <= 1:
        _fail("caller replay adVNTR numeric call values are outside their domains")
    return {"state": state, "read_support": support, "mean_coverage": coverage, "pvalue": pvalue}


def _capture_records(raw: bytes, expected_vntr_ids: tuple[int, ...]) -> tuple[dict[str, object], ...]:
    if not isinstance(raw, bytes) or not raw or not raw.endswith(b"\n"):
        _fail("caller adVNTR capture must be complete newline-terminated JSONL")
    records: list[dict[str, object]] = []
    for line in raw.splitlines():
        value = load_strict_json_object(line)
        if set(value) != _CAPTURE_FIELDS or value.get("schema_version") != "advntr-frameshift-capture-v2":
            _fail("caller adVNTR capture fields or schema differ")
        if hashlib.sha256(line).hexdigest() != advntr_canonical_sha256(value):
            _fail("caller adVNTR capture records must use canonical upstream bytes")
        records.append(value)
    loci = tuple(record.get("locus") for record in records)
    if any(not isinstance(locus, Mapping) for locus in loci):
        _fail("caller adVNTR capture locus must be an object")
    ids = tuple(cast(Mapping[str, object], locus).get("vntr_id") for locus in loci)
    if ids != expected_vntr_ids:
        _fail("caller adVNTR capture differs from the exact target roster")
    return tuple(records)


def decode_advntr_baseline_calls(
    capture_raw: bytes, expected_vntr_ids: tuple[int, ...]
) -> tuple[tuple[dict[str, object], ...], bool]:
    """Project independently recorded native baseline calls from completed captures."""
    calls: list[dict[str, object]] = []
    assessable = True
    for record in _capture_records(capture_raw, expected_vntr_ids):
        visits, warnings = record.get("decision_visits"), record.get("warnings")
        if not isinstance(visits, list) or not isinstance(warnings, list):
            _fail("caller adVNTR capture decisions and warnings must be lists")
        for warning in warnings:
            if not isinstance(warning, Mapping):
                _fail("caller adVNTR capture warning must be an object")
            if warning.get("origin") == "calibration-audit":
                assessable = False
        for visit in visits:
            if not isinstance(visit, Mapping):
                _fail("caller adVNTR capture decision visit must be an object")
            statistic, plan = visit.get("statistic"), visit.get("plan")
            if statistic is None:
                continue
            if (
                not isinstance(statistic, Mapping)
                or not isinstance(plan, Mapping)
                or type(statistic.get("called")) is not bool
            ):
                _fail("caller adVNTR capture decision receipt is invalid")
            if statistic["called"]:
                calls.append(
                    _advntr_call(
                        {
                            "state": plan.get("state"),
                            "read_support": plan.get("read_support"),
                            "mean_coverage": visit.get("mean_coverage"),
                            "pvalue": statistic.get("pvalue"),
                        }
                    )
                )
    return tuple(calls), assessable


def decode_advntr_replay_calls(
    raw: bytes,
    *,
    manifest_raw: bytes,
    policy_raw: bytes,
    capture_raw: bytes,
    expected_key: str,
    expected_vntr_ids: tuple[int, ...],
    expected_background_sha256: str | None,
) -> tuple[tuple[dict[str, object], ...], bool]:
    """Decode a stored upstream replay while checking every file/content binding."""
    records = _capture_records(capture_raw, expected_vntr_ids)
    manifest = load_strict_json_object(manifest_raw)
    policy = load_strict_json_object(policy_raw)
    output = load_strict_json_object(raw)
    output_fields = {
        "schema_version",
        "manifest_file_sha256",
        "manifest_sha256",
        "policy_file_sha256",
        "policy_sha256",
        "background_file_sha256",
        "replay_producer",
        "results",
    }
    if set(output) != output_fields or output.get("schema_version") != "advntr-frameshift-replay-output-v1":
        _fail("caller adVNTR replay output fields differ")
    bindings = (
        output.get("manifest_file_sha256") == hashlib.sha256(manifest_raw).hexdigest(),
        output.get("manifest_sha256") == advntr_canonical_sha256(manifest),
        output.get("policy_file_sha256") == hashlib.sha256(policy_raw).hexdigest(),
        output.get("policy_sha256") == advntr_canonical_sha256(policy),
        output.get("background_file_sha256") == expected_background_sha256,
    )
    if not all(bindings):
        _fail("caller adVNTR replay file or semantic binding differs")
    captures = manifest.get("captures")
    if not isinstance(captures, list) or len(captures) != 1 or not isinstance(captures[0], Mapping):
        _fail("caller adVNTR replay manifest must contain one exact member")
    capture = captures[0]
    capture_sha = hashlib.sha256(capture_raw).hexdigest()
    if (
        capture.get("key") != expected_key
        or capture.get("sha256") != capture_sha
        or capture.get("vntr_ids") != list(expected_vntr_ids)
    ):
        _fail("caller adVNTR replay manifest differs from the member capture")
    results = output.get("results")
    if not isinstance(results, list) or len(results) != 1 or not isinstance(results[0], Mapping):
        _fail("caller adVNTR replay must contain one exact member result")
    result = results[0]
    if (
        set(result) != {"key", "capture_sha256", "vntrs"}
        or result.get("key") != expected_key
        or result.get("capture_sha256") != capture_sha
    ):
        _fail("caller adVNTR replay result capture binding differs")
    vntrs = result.get("vntrs")
    if not isinstance(vntrs, list):
        _fail("caller adVNTR replay VNTR results must be a list")
    calls: list[dict[str, object]] = []
    assessable = True
    seen: list[int] = []
    locus_fields = {
        "schema_version",
        "vntr_id",
        "capture_record_sha256",
        "policy_sha256",
        "capture_producer",
        "capture_assets",
        "loaded_background_sha256",
        "baseline_parity",
        "decision_visits",
        "calls",
        "warnings",
        "capture_audit",
    }
    for index, row in enumerate(vntrs):
        if not isinstance(row, Mapping) or set(row) != {"vntr_id", "result"} or not isinstance(row["result"], Mapping):
            _fail("caller adVNTR replay VNTR row fields differ")
        locus = row["result"]
        vntr_id = row["vntr_id"]
        seen.append(cast(int, vntr_id))
        if (
            set(locus) != locus_fields
            or locus.get("schema_version") != "advntr-frameshift-replay-result-v1"
            or locus.get("vntr_id") != vntr_id
            or locus.get("capture_record_sha256") != advntr_canonical_sha256(records[index])
            or locus.get("policy_sha256") != advntr_canonical_sha256(policy)
            or locus.get("capture_producer") != records[index].get("producer")
            or locus.get("baseline_parity") is not True
        ):
            _fail("caller adVNTR replay locus identity or baseline binding differs")
        audit = locus.get("capture_audit")
        if not isinstance(audit, Mapping) or set(audit) != {
            "attribution_outside_trials",
            "calibrated_policy_domain_errors",
        }:
            _fail("caller adVNTR replay audit fields differ")
        for entries in audit.values():
            if not isinstance(entries, list) or any(not isinstance(item, str) for item in entries):
                _fail("caller adVNTR replay audit entries must be text lists")
            assessable = assessable and not entries
        raw_calls = locus.get("calls")
        if not isinstance(raw_calls, list):
            _fail("caller adVNTR replay calls must remain a list")
        calls.extend(_advntr_call(item) for item in raw_calls)
    if tuple(seen) != expected_vntr_ids:
        _fail("caller adVNTR replay result differs from the exact target roster")
    return tuple(calls), assessable


def _advntr_reconciliation(
    capture: KestrelCapture, value: object, presentation_index: int
) -> IdentityReconciliationObservation:
    if not isinstance(value, Mapping) or set(value) != {"state", "read_support", "mean_coverage", "pvalue"}:
        _fail("caller replay adVNTR call fields differ from the closed contract")
    state = _text(value["state"], "replay adVNTR state")
    support = _number(value["read_support"], "read support")
    coverage = _number(value["mean_coverage"], "mean coverage")
    pvalue = _number(value["pvalue"], "pvalue")
    if not isinstance(support, int) or support < 0 or coverage < 0 or not 0 <= pvalue <= 1:
        _fail("caller replay adVNTR numeric call values are outside their domains")
    raw_ru, raw_pos = derive_ru_and_pos((state,))
    repeat_units: tuple[str, ...] | None = tuple(raw_ru[0].split(","))
    positions: tuple[int, ...] | None
    try:
        positions = tuple(int(item) for item in raw_pos[0].split(","))
    except ValueError:
        repeat_units = None
        positions = None
    representation = AdvntrRepresentation(state, repeat_units, positions)
    translation: IdentityTranslation = capture.identity_component.translate_advntr(representation)
    calls = from_advntr(state)
    net = sum(call.net_length for call in calls)
    event = calls[0].event if calls and len({call.event for call in calls}) == 1 else "delins"
    flags = frozenset(flag for call in calls for flag in call.flags)
    names = tuple(call.name for call in calls if call.name is not None)
    return IdentityReconciliationObservation(
        translation,
        names[0] if len(names) == 1 else None,
        "advntr",
        event,
        net,
        flags,
        IdentityDisposition("admissible"),
        len(names) == 1 and names[0] in _known_names(),
        advntr_sequencing_read_support=support,
        advntr_mean_coverage=coverage,
        presentation_call_index=presentation_index if len(names) == 1 else None,
    )


def replayed_caller_observation(
    member: EligibleCallerMember,
    truth: CallerTruth,
    *,
    kestrel_capture: KestrelCapture,
    policy: CallerPolicyValues,
    capture_policy_sha256: str,
    advntr_calls: tuple[object, ...],
    advntr_assessable: bool,
    source_evidence_sha256: str,
) -> CallerEvidenceRow:
    """Replay production caller decisions and bind the result to one truth row."""
    caller_truth_document(truth)
    if not isinstance(member, EligibleCallerMember) or member.key not in truth.by_key:
        _fail("caller replay member is absent from sealed truth")
    if _SHA256.fullmatch(source_evidence_sha256) is None:
        _fail("caller replay source evidence requires a lowercase SHA256 digest")
    if type(advntr_assessable) is not bool or not isinstance(advntr_calls, tuple):
        _fail("caller replay adVNTR evidence must be a typed tuple and audit decision")
    replay = replay_kestrel_capture(
        kestrel_capture,
        policy,
        capture_policy_sha256=capture_policy_sha256,
    )
    if not advntr_assessable:
        disposition: EvidenceDisposition = "unsupported"
        called_positive: bool | None = None
        variants: tuple[str, ...] = ()
        tier_a: tuple[str, ...] = ()
    else:
        observations: list[IdentityReconciliationObservation] = []
        selected = kestrel_replay_selected_frame(replay)
        if not selected.empty:
            observations.append(_kestrel_reconciliation(kestrel_capture, selected.iloc[0].to_dict()))
        observations.extend(
            _advntr_reconciliation(kestrel_capture, value, len(observations) + index)
            for index, value in enumerate(advntr_calls)
        )
        if observations:
            result = reconcile_identity_observations(tuple(observations), _reconciliation_policy())
            identity = result.decision.identity
            variants = () if identity is None else (serialize_molecular_identity(identity),)
            called_positive = True
            disposition = "called"
            tier_a = variants if result.decision.tier == "A" else ()
        elif replay.source_row_count == 0:
            disposition, called_positive, variants, tier_a = "zero-candidate", False, (), ()
        else:
            disposition, called_positive, variants, tier_a = "no-call", None, (), ()
    truth_positive, truth_variants = _truth_values(truth.by_key[member.key])
    return CallerEvidenceRow(
        CallerObservation(
            member.key,
            member.group_key,
            truth_positive,
            truth_variants,
            called_positive,
            variants,
            tier_a,
        ),
        disposition,
        source_evidence_sha256,
    )
