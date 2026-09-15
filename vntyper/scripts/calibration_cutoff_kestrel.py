"""Efficient evidence-bound Kestrel scalar-policy grid replay."""

from __future__ import annotations

import hashlib
import logging
import math
import multiprocessing
import re
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_caller_policy import (
    KESTREL_CALLER_POLICY_POINTERS,
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_kestrel_capture import decode_kestrel_capture
from vntyper.scripts.calibration_kestrel_replay import (
    KestrelReplayResult,
    decode_kestrel_replay_result,
    kestrel_replay_document,
    kestrel_replay_selected_frame,
    replay_kestrel_capture,
)
from vntyper.scripts.calibration_run_extraction import _parse_tsv, _rows
from vntyper.scripts.calibration_run_projection import is_kestrel_negative_placeholder
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object
from vntyper.scripts.identity_candidate_persistence import parse_selected_candidate_cells
from vntyper.scripts.molecular_identity import parse_molecular_identity, serialize_molecular_identity

logger = logging.getLogger(__name__)

GridDisposition = Literal["called", "no-call", "unassessable-no-candidates"]
BaselineParity = Literal["native-exact", "capture-replay-authoritative"]

_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
# Replay-only bookkeeping: the capture row ordinal exists solely inside calibration.
_PRIVATE_COLUMN_PREFIX = "__Calibration_"
# The native comparison is worthless unless it covers the actual selection decision:
# the variant coordinates, the score the floor is applied to and the assigned band.
_NATIVE_PARITY_REQUIRED = ("POS", "REF", "ALT", "Depth_Score", "Confidence")
_OBSERVATION_FIELDS = {
    "key",
    "policy_id",
    "policy_sha256",
    "kestrel_policy_sha256",
    "capture_sha256",
    "replay_sha256",
    "disposition",
    "called_positive",
    "confidence",
    "flag",
    "canonical_identity",
    "translation_status",
    "translation_failure",
    "context_diverges",
}


@dataclass(frozen=True)
class KestrelGridObservation:
    """One selected Kestrel endpoint under one complete caller policy."""

    key: str
    policy_id: str
    policy_sha256: str
    kestrel_policy_sha256: str
    capture_sha256: str
    replay_sha256: str
    disposition: GridDisposition
    called_positive: bool | None
    confidence: str | None
    flag: str | None
    canonical_identity: str | None
    translation_status: str | None
    translation_failure: str | None
    context_diverges: bool | None


@dataclass(frozen=True)
class KestrelGridReplay:
    """Immutable policy-major replay evidence and exact local input commitments."""

    policy_ids: tuple[str, ...]
    sample_keys: tuple[str, ...]
    observations: Mapping[str, Mapping[str, KestrelGridObservation]]
    policy_sha256: Mapping[str, str]
    kestrel_policy_sha256: Mapping[str, str]
    capture_file_sha256: Mapping[str, str]
    native_file_sha256: Mapping[str, str | None]
    baseline_parity: Mapping[str, BaselineParity]
    baseline_policy_sha256: str
    unique_kestrel_parameter_count: int
    input_sha256: str
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        _fail(f"Kestrel cutoff {label} must be a lowercase SHA256 digest")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip() or "\x00" in value:
        _fail(f"Kestrel cutoff {label} must be nonempty trimmed text")
    return value


def _read(path: object, label: str) -> bytes:
    if not isinstance(path, Path):
        _fail(f"Kestrel cutoff {label} must be a Path")
    try:
        return read_regular_path(path)
    except (OSError, ValueError) as error:
        raise ValueError(f"Kestrel cutoff {label} is missing, unsafe, or unreadable") from error


def _capture_documents(capture_paths: Mapping[str, Path]) -> tuple[tuple[str, ...], dict[str, object], dict[str, str]]:
    if not isinstance(capture_paths, Mapping) or not capture_paths:
        _fail("Kestrel cutoff capture paths require a nonempty mapping")
    validated_keys = tuple(_text(key, "sample key") for key in capture_paths)
    keys = tuple(sorted(validated_keys))
    if len(keys) != len(capture_paths):
        _fail("Kestrel cutoff sample keys must be unique")
    documents: dict[str, object] = {}
    digests: dict[str, str] = {}
    for key in keys:
        raw = _read(capture_paths[key], f"capture for {key}")
        digests[key] = hashlib.sha256(raw).hexdigest()
        try:
            documents[key] = load_strict_json_object(raw)
        except ValueError as error:
            raise ValueError(f"Kestrel cutoff capture for {key} is not strict JSON") from error
    return keys, documents, digests


def _candidate_documents(
    policies: Mapping[str, CallerPolicyValues],
) -> tuple[tuple[str, ...], dict[str, dict[str, object]], dict[str, str]]:
    if not isinstance(policies, Mapping):
        _fail("Kestrel cutoff policies require a mapping")
    documents: dict[str, dict[str, object]] = {}
    digests: dict[str, str] = {}
    policy_ids = tuple(sorted(_text(policy_id, "policy ID") for policy_id in policies))
    for policy_id in policy_ids:
        name = _text(policy_id, "policy ID")
        if name == "baseline":
            _fail("Kestrel cutoff policy ID baseline is reserved")
        policy = policies[policy_id]
        if not isinstance(policy, CallerPolicyValues):
            _fail("Kestrel cutoff policies require immutable CallerPolicyValues")
        documents[name] = caller_policy_values_document(policy)
        digests[name] = policy.sha256
    return ("baseline", *tuple(documents)), documents, digests


def _kestrel_parameter_sha(policy: CallerPolicyValues) -> str:
    return canonical_sha256(
        {
            "schema_version": "calibration-kestrel-parameters-v1",
            "values": {pointer: policy.values[pointer] for pointer in KESTREL_CALLER_POLICY_POINTERS},
        }
    )


def _effective_policy(baseline: CallerPolicyValues, candidate: CallerPolicyValues) -> CallerPolicyValues:
    if candidate.required_callers != baseline.required_callers:
        _fail("Kestrel cutoff candidate caller set differs from the captured baseline")
    values = dict(baseline.values)
    for pointer in KESTREL_CALLER_POLICY_POINTERS:
        values[pointer] = candidate.values[pointer]
    return decode_caller_policy_values(
        {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": list(baseline.required_callers),
            "values": values,
        }
    )


def _selected_fields(result: KestrelReplayResult) -> dict[str, object] | None:
    """Reduce one replay to the selected row's decision-bearing fields.

    Only these fields reach the parent. The prefilter population stays inside the
    worker, because carrying it back would make peak memory scale with samples times
    policies times candidate rows.

    Args:
        result: One complete replay of one capture under one policy.

    Returns:
        The selected row's fields, or None when the policy selected nothing.

    Raises:
        ValueError: If a selected row carries empty confidence or flag text.
    """
    selected = kestrel_replay_selected_frame(result)
    if selected.empty:
        return None
    row = selected.iloc[0].to_dict()
    persisted = parse_selected_candidate_cells(row)
    identity = persisted.translation.identity
    confidence, flag = row.get("Confidence"), row.get("Flag")
    if not isinstance(confidence, str) or not confidence or not isinstance(flag, str) or not flag:
        _fail("Kestrel cutoff selected confidence and Flag must be nonempty text")
    return {
        "confidence": confidence,
        "flag": flag,
        "canonical_identity": None if identity is None else serialize_molecular_identity(identity),
        "translation_status": persisted.translation.status,
        "translation_failure": persisted.translation.failure,
        "context_diverges": persisted.translation.context_diverges,
    }


def _replay_case(
    item: tuple[str, object, Mapping[str, dict[str, object]], bool],
) -> tuple[str, str, str, dict[str, str], dict[str, dict[str, object]], dict[str, object] | None]:
    key, capture_document, candidate_documents, needs_native_baseline = item
    capture = decode_kestrel_capture(capture_document)
    candidates = {name: decode_caller_policy_values(value) for name, value in candidate_documents.items()}
    effective: dict[str, CallerPolicyValues] = {"baseline": capture.baseline_policy}
    effective.update({name: _effective_policy(capture.baseline_policy, policy) for name, policy in candidates.items()})
    parameter_by_policy = {name: _kestrel_parameter_sha(policy) for name, policy in effective.items()}
    unique: dict[str, CallerPolicyValues] = {}
    for name in ("baseline", *tuple(sorted(candidates))):
        unique.setdefault(parameter_by_policy[name], effective[name])
    baseline_parameter = parameter_by_policy["baseline"]
    replayed: dict[str, dict[str, object]] = {}
    baseline_document: dict[str, object] | None = None
    for parameter_sha, policy in unique.items():
        result = replay_kestrel_capture(
            capture,
            policy,
            capture_policy_sha256=capture.provenance.capture_policy_sha256,
        )
        replayed[parameter_sha] = {
            "policy_sha256": result.policy_sha256,
            "capture_sha256": result.capture_sha256,
            "replay_sha256": result.sha256,
            "disposition": result.disposition,
            "selected": _selected_fields(result),
        }
        if needs_native_baseline and parameter_sha == baseline_parameter:
            baseline_document = kestrel_replay_document(result)
    return (
        key,
        capture.sha256,
        capture.baseline_policy.sha256,
        parameter_by_policy,
        replayed,
        baseline_document,
    )


def _native_inputs(
    keys: tuple[str, ...], native_paths: Mapping[str, Path | None] | None
) -> tuple[dict[str, bytes | None], dict[str, str | None]]:
    if native_paths is None:
        return (dict.fromkeys(keys), dict.fromkeys(keys))
    if not isinstance(native_paths, Mapping) or set(native_paths) != set(keys):
        _fail("Kestrel cutoff native paths must cover the exact capture roster")
    raw: dict[str, bytes | None] = {}
    digests: dict[str, str | None] = {}
    for key in keys:
        path = native_paths[key]
        if path is None:
            raw[key], digests[key] = None, None
        else:
            value = _read(path, f"native Kestrel result for {key}")
            raw[key], digests[key] = value, hashlib.sha256(value).hexdigest()
    return raw, digests


def _cell_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value)


def validate_native_kestrel_baseline_parity(raw: bytes, baseline: KestrelReplayResult) -> Literal["called", "negative"]:
    """Require a native final Kestrel TSV to equal its baseline capture replay.

    The replay's selected frame is compared over the columns it shares with the
    native header, excluding calibration-private ``__Calibration_*`` columns that
    only the replay defines and production never publishes. The shared set must
    still contain every decision-bearing column, so the comparison can never
    degrade into agreeing about nothing.

    Args:
        raw: Exact independently retained native result bytes.
        baseline: Baseline replay from the same complete capture.

    Returns:
        Whether the exact native endpoint is called or negative.

    Raises:
        ValueError: If bytes are malformed, the shared columns omit a required
            decision column, or a shared native field differs from the replay.
    """
    if not isinstance(raw, bytes):
        _fail("Kestrel cutoff native baseline must be exact bytes")
    kestrel_replay_document(baseline)
    native_rows = _rows(_parse_tsv(raw, "native Kestrel baseline"), "native Kestrel baseline", allow_empty=True)
    selected = kestrel_replay_selected_frame(baseline)
    if selected.empty:
        if len(native_rows) != 1 or not is_kestrel_negative_placeholder(native_rows[0]):
            _fail("Kestrel cutoff native Kestrel baseline differs from capture replay")
        return "negative"
    if len(native_rows) != 1 or is_kestrel_negative_placeholder(native_rows[0]):
        _fail("Kestrel cutoff native Kestrel baseline differs from capture replay")
    observed = native_rows[0]
    comparable = {
        column: value
        for column, value in selected.iloc[0].to_dict().items()
        if not column.startswith(_PRIVATE_COLUMN_PREFIX) and column in observed
    }
    missing = tuple(column for column in _NATIVE_PARITY_REQUIRED if column not in comparable)
    if missing:
        _fail(f"Kestrel cutoff native Kestrel baseline is missing required comparable columns: {', '.join(missing)}")
    if any(observed[column] != _cell_text(value) for column, value in comparable.items()):
        _fail("Kestrel cutoff native Kestrel baseline differs from capture replay selected fields")
    return "called"


def _observation(
    key: str,
    policy_id: str,
    policy_sha256: str,
    record: Mapping[str, object],
    *,
    empty_native_negative: bool,
) -> KestrelGridObservation:
    """Compose one endpoint from a worker's reduced replay record.

    Args:
        key: Sample key the capture was associated with.
        policy_id: Requested policy identifier.
        policy_sha256: Digest of the requested complete policy.
        record: Reduced replay fields returned by :func:`_replay_case`.
        empty_native_negative: Whether an exact native Negative placeholder proved the
            run completed and found nothing.

    Returns:
        The immutable per-sample, per-policy observation.

    Raises:
        ValueError: If the reduced record is malformed.
    """
    if set(record) != {"policy_sha256", "capture_sha256", "replay_sha256", "disposition", "selected"}:
        _fail("Kestrel cutoff replay record fields differ from the reduced contract")
    disposition = cast(GridDisposition, record["disposition"])
    replay_policy_sha256 = cast(str, record["policy_sha256"])
    capture_sha256 = cast(str, record["capture_sha256"])
    replay_sha256 = cast(str, record["replay_sha256"])
    selected = record["selected"]
    if selected is None:
        called = False if disposition == "no-call" or empty_native_negative else None
        return KestrelGridObservation(
            key,
            policy_id,
            policy_sha256,
            replay_policy_sha256,
            capture_sha256,
            replay_sha256,
            disposition,
            called,
            None,
            None,
            None,
            None,
            None,
            None,
        )
    fields = cast(Mapping[str, object], selected)
    return KestrelGridObservation(
        key,
        policy_id,
        policy_sha256,
        replay_policy_sha256,
        capture_sha256,
        replay_sha256,
        disposition,
        True,
        cast(str, fields["confidence"]),
        cast(str, fields["flag"]),
        cast("str | None", fields["canonical_identity"]),
        cast("str | None", fields["translation_status"]),
        cast("str | None", fields["translation_failure"]),
        cast("bool | None", fields["context_diverges"]),
    )


def _observation_document(value: KestrelGridObservation) -> dict[str, object]:
    if not isinstance(value, KestrelGridObservation):
        _fail("Kestrel cutoff observation requires KestrelGridObservation")
    for field in ("key", "policy_id"):
        _text(getattr(value, field), f"observation {field}")
    for field in ("policy_sha256", "kestrel_policy_sha256", "capture_sha256", "replay_sha256"):
        _digest(getattr(value, field), f"observation {field}")
    if value.disposition not in {"called", "no-call", "unassessable-no-candidates"}:
        _fail("Kestrel cutoff observation disposition is invalid")
    if value.disposition == "called":
        if value.called_positive is not True:
            _fail("Kestrel cutoff called replay must be a positive endpoint")
        for field in ("confidence", "flag", "translation_status"):
            _text(getattr(value, field), f"observation {field}")
        if value.context_diverges is None or type(value.context_diverges) is not bool:
            _fail("Kestrel cutoff called replay requires boolean identity context")
        if value.canonical_identity is not None:
            identity = parse_molecular_identity(value.canonical_identity)
            if serialize_molecular_identity(identity) != value.canonical_identity:
                _fail("Kestrel cutoff observation canonical identity is noncanonical")
        if value.translation_failure is not None:
            _text(value.translation_failure, "observation translation failure")
    elif any(
        item is not None
        for item in (
            value.confidence,
            value.flag,
            value.canonical_identity,
            value.translation_status,
            value.translation_failure,
            value.context_diverges,
        )
    ):
        _fail("Kestrel cutoff unselected replay cannot carry selected metadata")
    if value.disposition == "no-call" and value.called_positive is not False:
        _fail("Kestrel cutoff complete filtered no-call must be negative")
    if value.disposition == "unassessable-no-candidates" and (
        value.called_positive is not None and value.called_positive is not False
    ):
        _fail("Kestrel cutoff empty capture endpoint is invalid")
    return {field: getattr(value, field) for field in _OBSERVATION_FIELDS}


def _input_document(result: KestrelGridReplay) -> dict[str, object]:
    return {
        "schema_version": "calibration-kestrel-grid-input-v1",
        "baseline_policy_sha256": result.baseline_policy_sha256,
        "policies": [
            {
                "policy_id": policy_id,
                "policy_sha256": result.policy_sha256[policy_id],
                "kestrel_policy_sha256": result.kestrel_policy_sha256[policy_id],
            }
            for policy_id in result.policy_ids
        ],
        "samples": [
            {
                "key": key,
                "capture_file_sha256": result.capture_file_sha256[key],
                "native_file_sha256": result.native_file_sha256[key],
            }
            for key in result.sample_keys
        ],
    }


def _document(result: KestrelGridReplay) -> dict[str, object]:
    return {
        "schema_version": "calibration-kestrel-grid-replay-v1",
        "input_sha256": result.input_sha256,
        "baseline_policy_sha256": result.baseline_policy_sha256,
        "unique_kestrel_parameter_count": result.unique_kestrel_parameter_count,
        "policy_ids": list(result.policy_ids),
        "sample_keys": list(result.sample_keys),
        "policy_sha256": dict(result.policy_sha256),
        "kestrel_policy_sha256": dict(result.kestrel_policy_sha256),
        "capture_file_sha256": dict(result.capture_file_sha256),
        "native_file_sha256": dict(result.native_file_sha256),
        "baseline_parity": dict(result.baseline_parity),
        "observations": [
            {
                "policy_id": policy_id,
                "samples": [_observation_document(result.observations[policy_id][key]) for key in result.sample_keys],
            }
            for policy_id in result.policy_ids
        ],
    }


def _require_result(result: KestrelGridReplay) -> None:
    if not isinstance(result, KestrelGridReplay):
        _fail("Kestrel cutoff projection requires KestrelGridReplay")
    if (
        not isinstance(result.policy_ids, tuple)
        or not isinstance(result.sample_keys, tuple)
        or type(result.observations) is not _MAPPING_PROXY_TYPE
        or any(type(value) is not _MAPPING_PROXY_TYPE for value in result.observations.values())
        or any(
            type(value) is not _MAPPING_PROXY_TYPE
            for value in (
                result.policy_sha256,
                result.kestrel_policy_sha256,
                result.capture_file_sha256,
                result.native_file_sha256,
                result.baseline_parity,
            )
        )
    ):
        _fail("Kestrel cutoff projection requires immutable decoded content")
    expected_policies = ("baseline", *tuple(sorted(set(result.policy_ids) - {"baseline"})))
    if result.policy_ids != expected_policies or result.sample_keys != tuple(sorted(set(result.sample_keys))):
        _fail("Kestrel cutoff projection policy or sample ordering differs")
    if any(set(value) != set(result.sample_keys) for value in result.observations.values()) or set(
        result.observations
    ) != set(result.policy_ids):
        _fail("Kestrel cutoff projection observation roster differs")
    if any(set(value) != set(result.policy_ids) for value in (result.policy_sha256, result.kestrel_policy_sha256)):
        _fail("Kestrel cutoff projection policy digest roster differs")
    if any(
        set(value) != set(result.sample_keys)
        for value in (result.capture_file_sha256, result.native_file_sha256, result.baseline_parity)
    ):
        _fail("Kestrel cutoff projection input roster differs")
    for policy_id in result.policy_ids:
        _digest(result.policy_sha256[policy_id], "policy digest")
        _digest(result.kestrel_policy_sha256[policy_id], "effective Kestrel policy digest")
        for key in result.sample_keys:
            row = result.observations[policy_id][key]
            _observation_document(row)
            if (
                row.policy_id != policy_id
                or row.key != key
                or row.policy_sha256 != result.policy_sha256[policy_id]
                or row.kestrel_policy_sha256 != result.kestrel_policy_sha256[policy_id]
            ):
                _fail("Kestrel cutoff observation differs from its policy or sample binding")
    if result.baseline_policy_sha256 != result.policy_sha256["baseline"]:
        _fail("Kestrel cutoff baseline digest differs from baseline policy")
    if (
        isinstance(result.unique_kestrel_parameter_count, bool)
        or not isinstance(result.unique_kestrel_parameter_count, int)
        or result.unique_kestrel_parameter_count != len(set(result.kestrel_policy_sha256.values()))
    ):
        _fail("Kestrel cutoff unique parameter count differs from the effective policies")
    for key in result.sample_keys:
        _digest(result.capture_file_sha256[key], "capture file digest")
        native_digest = result.native_file_sha256[key]
        if native_digest is not None:
            _digest(native_digest, "native file digest")
        if (native_digest is None) != (result.baseline_parity[key] == "capture-replay-authoritative"):
            _fail("Kestrel cutoff native digest and baseline parity differ")
        capture_shas = {result.observations[policy_id][key].capture_sha256 for policy_id in result.policy_ids}
        if len(capture_shas) != 1:
            _fail("Kestrel cutoff policy observations differ in capture identity")
    if (
        canonical_sha256(_input_document(result)) != result.input_sha256
        or canonical_sha256(_document(result)) != result.sha256
    ):
        _fail("Kestrel cutoff replay differs from its canonical content or digest")


def kestrel_grid_replay_document(result: KestrelGridReplay) -> dict[str, object]:
    """Project a validated immutable grid replay as fresh JSON-compatible content."""
    _require_result(result)
    return _document(result)


def replay_kestrel_grid(
    capture_paths: Mapping[str, Path],
    policies: Mapping[str, CallerPolicyValues],
    *,
    workers: int = 1,
    native_paths: Mapping[str, Path | None] | None = None,
) -> KestrelGridReplay:
    """Read each complete capture once and replay a finite Kestrel policy grid.

    The external mapping key is the sample/capture association: capture-v1 has no
    internal sample identity. A complete nonempty capture whose candidates are all
    filtered is a scientific negative. An empty capture is unavailable unless an
    independently stored exact native Negative result proves completed execution.

    Args:
        capture_paths: Opaque sample keys mapped to complete capture JSON files.
        policies: Named complete caller policies. ``baseline`` is reserved and
            derived from the captures.
        workers: Positive process count; one runs directly in the caller process.
        native_paths: Optional exact roster mapping to native Kestrel TSVs or null.

    Returns:
        Immutable policy-major endpoint evidence and exact input digests.

    Raises:
        ValueError: For unsafe/malformed inputs, incompatible policies, divergent
            baselines, replay failures, or native baseline parity failures.
    """
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        _fail("Kestrel cutoff workers must be a positive integer")
    keys, capture_documents, capture_hashes = _capture_documents(capture_paths)
    policy_ids, candidate_documents, candidate_hashes = _candidate_documents(policies)
    native_raw, native_hashes = _native_inputs(keys, native_paths)
    tasks = [(key, capture_documents[key], candidate_documents, native_raw[key] is not None) for key in keys]
    if workers == 1:
        case_rows = [_replay_case(task) for task in tasks]
    else:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(tasks)), mp_context=multiprocessing.get_context("spawn")
        ) as executor:
            case_rows = list(executor.map(_replay_case, tasks))
    baselines = {row[2] for row in case_rows}
    if len(baselines) != 1:
        _fail("Kestrel cutoff capture baseline policies differ")
    baseline_sha = next(iter(baselines))
    policy_hashes = {"baseline": baseline_sha, **candidate_hashes}
    parameter_by_case = {key: values for key, _capture, _baseline, values, _replays, _document in case_rows}
    first_parameters = parameter_by_case[keys[0]]
    if any(values != first_parameters for values in parameter_by_case.values()):
        _fail("Kestrel cutoff effective policy projections differ across captures")
    result_by_case: dict[str, dict[str, dict[str, object]]] = {}
    capture_semantic: dict[str, str] = {}
    parity: dict[str, BaselineParity] = {}
    empty_native_negative: dict[str, bool] = {}
    for key, capture_sha, _baseline, parameter_map, replay_records, baseline_document in case_rows:
        capture_semantic[key] = capture_sha
        result_by_case[key] = replay_records
        native = native_raw[key]
        if native is None:
            parity[key] = "capture-replay-authoritative"
            empty_native_negative[key] = False
        else:
            if baseline_document is None:
                _fail(f"Kestrel cutoff baseline replay is missing for native parity on {key}")
            baseline = decode_kestrel_replay_result(baseline_document)
            if baseline.policy_sha256 != replay_records[parameter_map["baseline"]]["policy_sha256"]:
                _fail(f"Kestrel cutoff baseline replay identity differs from its reduced record on {key}")
            empty_native_negative[key] = validate_native_kestrel_baseline_parity(native, baseline) == "negative"
            parity[key] = "native-exact"
    observations = {
        policy_id: MappingProxyType(
            {
                key: _observation(
                    key,
                    policy_id,
                    policy_hashes[policy_id],
                    result_by_case[key][parameter_by_case[key][policy_id]],
                    empty_native_negative=empty_native_negative[key],
                )
                for key in keys
            }
        )
        for policy_id in policy_ids
    }
    draft = KestrelGridReplay(
        policy_ids,
        keys,
        MappingProxyType(observations),
        MappingProxyType(policy_hashes),
        MappingProxyType(
            {
                policy_id: cast(str, result_by_case[keys[0]][first_parameters[policy_id]]["policy_sha256"])
                for policy_id in policy_ids
            }
        ),
        MappingProxyType(capture_hashes),
        MappingProxyType(native_hashes),
        MappingProxyType(parity),
        baseline_sha,
        len(set(first_parameters.values())),
        "",
        "",
    )
    bound = replace(draft, input_sha256=canonical_sha256(_input_document(draft)))
    return replace(bound, sha256=canonical_sha256(_document(bound)))
