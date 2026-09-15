"""Production Kestrel replay over complete immutable prefilter captures."""

from __future__ import annotations

import logging
import math
import re
from copy import deepcopy
from dataclasses import dataclass
from typing import Literal, NoReturn, cast

import pandas as pd

from vntyper.scripts.calibration_caller_policy import (
    KESTREL_CALLER_POLICY_POINTERS,
    CallerPolicyValues,
    caller_policy_values_document,
)
from vntyper.scripts.calibration_kestrel_capture import (
    KESTREL_MOTIF_COLUMNS,
    KESTREL_RAW_COLUMNS,
    KestrelCapture,
    kestrel_capture_document,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.kestrel_genotyping import add_haplo_count, select_single_best_variant
from vntyper.scripts.kestrel_postprocessing import FILTER_COLUMNS, evaluate_kestrel_candidates

logger = logging.getLogger(__name__)

_SOURCE_ORDINAL = "__Calibration_Source_Row_Ordinal"
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_ROOT_FIELDS = {
    "schema_version",
    "capture_sha256",
    "policy_sha256",
    "disposition",
    "source_row_count",
    "prefilter",
    "selected",
    "source_dispositions",
}
_TABLE_FIELDS = {"columns", "rows"}
_DISPOSITION_FIELDS = {"source_row_ordinal", "blocking_gates", "selected"}
_POINTER_PATHS = {
    pointer: tuple(pointer.removeprefix("/").split("/")[2:]) for pointer in KESTREL_CALLER_POLICY_POINTERS
}

ReplayDisposition = Literal["called", "no-call", "unassessable-no-candidates"]


@dataclass(frozen=True)
class KestrelSourceDisposition:
    """One preserved source row and the exact gates that excluded it."""

    source_row_ordinal: int
    blocking_gates: tuple[str, ...]
    selected: bool


@dataclass(frozen=True)
class KestrelReplayResult:
    """Immutable replay tables, row dispositions and canonical identity."""

    capture_sha256: str
    policy_sha256: str
    disposition: ReplayDisposition
    source_row_count: int
    prefilter_json: bytes
    selected_json: bytes
    source_dispositions: tuple[KestrelSourceDisposition, ...]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        _fail(f"Kestrel replay {label} must be a lowercase SHA256 digest")
    return value


def _cell(value: object) -> str | bool | int | float | None:
    if value is None or value is pd.NA:
        return None
    item = getattr(value, "item", None)
    if callable(item):
        value = item()
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return value
    if isinstance(value, str):
        return value
    _fail("Kestrel replay tables contain a non-scalar cell")


def _table_document(frame: pd.DataFrame) -> dict[str, object]:
    columns = [str(column) for column in frame.columns]
    return {
        "columns": columns,
        "rows": [[_cell(value) for value in row] for row in frame.itertuples(index=False, name=None)],
    }


def _decode_table(value: object, label: str) -> bytes:
    if not isinstance(value, dict) or set(value) != _TABLE_FIELDS:
        _fail(f"Kestrel replay {label} table fields differ")
    columns, rows = value["columns"], value["rows"]
    if (
        not isinstance(columns, list)
        or any(not isinstance(column, str) or not column for column in columns)
        or len(columns) != len(set(columns))
        or not isinstance(rows, list)
    ):
        _fail(f"Kestrel replay {label} table columns or rows are invalid")
    for row in rows:
        if not isinstance(row, list) or len(row) != len(columns):
            _fail(f"Kestrel replay {label} row width differs from its columns")
        for value_cell in row:
            if isinstance(value_cell, float) and not math.isfinite(value_cell):
                _fail(f"Kestrel replay {label} contains a non-finite cell")
            if value_cell is not None and not isinstance(value_cell, (str, bool, int, float)):
                _fail(f"Kestrel replay {label} contains a non-scalar cell")
    return canonical_json_bytes(value)


def _frame(encoded: bytes) -> pd.DataFrame:
    table = load_strict_json_object(encoded)
    return pd.DataFrame(table["rows"], columns=table["columns"])


def _document(result: KestrelReplayResult) -> dict[str, object]:
    return {
        "schema_version": "calibration-kestrel-replay-v1",
        "capture_sha256": result.capture_sha256,
        "policy_sha256": result.policy_sha256,
        "disposition": result.disposition,
        "source_row_count": result.source_row_count,
        "prefilter": load_strict_json_object(result.prefilter_json),
        "selected": load_strict_json_object(result.selected_json),
        "source_dispositions": [
            {
                "source_row_ordinal": row.source_row_ordinal,
                "blocking_gates": list(row.blocking_gates),
                "selected": row.selected,
            }
            for row in result.source_dispositions
        ],
    }


def decode_kestrel_replay_result(value: object) -> KestrelReplayResult:
    """Decode one strict immutable Kestrel replay result document.

    Args:
        value: Closed JSON-compatible replay document.

    Returns:
        Immutable replay output with a derived canonical digest.

    Raises:
        ValueError: If fields, tables, gates or dispositions are invalid.
    """
    if not isinstance(value, dict) or set(value) != _ROOT_FIELDS:
        _fail("Kestrel replay root fields differ from the closed contract")
    if value["schema_version"] != "calibration-kestrel-replay-v1":
        _fail("Kestrel replay schema_version is unsupported")
    capture_sha = _digest(value["capture_sha256"], "capture_sha256")
    policy_sha = _digest(value["policy_sha256"], "policy_sha256")
    disposition = value["disposition"]
    if disposition not in {"called", "no-call", "unassessable-no-candidates"}:
        _fail("Kestrel replay disposition is invalid")
    count = value["source_row_count"]
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        _fail("Kestrel replay source_row_count must be a nonnegative integer")
    prefilter = _decode_table(value["prefilter"], "prefilter")
    selected = _decode_table(value["selected"], "selected")
    prefilter_document = load_strict_json_object(prefilter)
    selected_document = load_strict_json_object(selected)
    prefilter_columns = cast(list[str], prefilter_document["columns"])
    selected_columns = cast(list[str], selected_document["columns"])
    prefilter_rows = cast(list[list[object]], prefilter_document["rows"])
    selected_rows = cast(list[list[object]], selected_document["rows"])
    if len(prefilter_rows) != count:
        _fail("Kestrel replay prefilter rows differ from source_row_count")
    if len(selected_rows) > 1:
        _fail("Kestrel replay selected table must contain at most one row")
    if count and (_SOURCE_ORDINAL not in prefilter_columns or _SOURCE_ORDINAL not in selected_columns):
        _fail("Kestrel replay prefilter and selected tables must retain source ordinals")
    prefilter_ordinals: list[int] = []
    if count:
        ordinal_column = prefilter_columns.index(_SOURCE_ORDINAL)
        for row in prefilter_rows:
            ordinal = row[ordinal_column]
            if isinstance(ordinal, bool) or not isinstance(ordinal, int):
                _fail("Kestrel replay prefilter source ordinals must be integers")
            prefilter_ordinals.append(ordinal)
        if prefilter_ordinals != list(range(count)):
            _fail("Kestrel replay prefilter source ordinals must be complete and ordered")
    selected_ordinal: int | None = None
    if selected_rows:
        selected_value = selected_rows[0][selected_columns.index(_SOURCE_ORDINAL)]
        if (
            isinstance(selected_value, bool)
            or not isinstance(selected_value, int)
            or selected_value not in range(count)
        ):
            _fail("Kestrel replay selected source ordinal is invalid")
        selected_ordinal = selected_value
    if count and any(gate not in prefilter_columns for gate in FILTER_COLUMNS):
        _fail("Kestrel replay prefilter table is missing a production gate")
    rows = value["source_dispositions"]
    if not isinstance(rows, list) or len(rows) != count:
        _fail("Kestrel replay source dispositions differ from source_row_count")
    dispositions = []
    for expected, item in enumerate(rows):
        if not isinstance(item, dict) or set(item) != _DISPOSITION_FIELDS:
            _fail("Kestrel replay source disposition fields differ")
        ordinal, gates, selected_value = item["source_row_ordinal"], item["blocking_gates"], item["selected"]
        if isinstance(ordinal, bool) or not isinstance(ordinal, int) or ordinal != expected:
            _fail("Kestrel replay source dispositions must be contiguous and ordered")
        if (
            not isinstance(gates, list)
            or any(gate not in FILTER_COLUMNS for gate in gates)
            or gates != [gate for gate in FILTER_COLUMNS if gate in gates]
            or len(gates) != len(set(gates))
            or type(selected_value) is not bool
        ):
            _fail("Kestrel replay blocking gates or selected marker are invalid")
        if count:
            source_row = prefilter_rows[ordinal]
            actual_gates: list[str] = []
            for gate in FILTER_COLUMNS:
                gate_value = source_row[prefilter_columns.index(gate)]
                if type(gate_value) is not bool:
                    _fail("Kestrel replay prefilter production gates must be boolean")
                if not gate_value:
                    actual_gates.append(gate)
            if gates != actual_gates:
                _fail("Kestrel replay source disposition differs from prefilter gates")
        dispositions.append(KestrelSourceDisposition(ordinal, tuple(gates), selected_value))
    selected_markers = [row.source_row_ordinal for row in dispositions if row.selected]
    if selected_markers != ([] if selected_ordinal is None else [selected_ordinal]):
        _fail("Kestrel replay selected table differs from source disposition markers")
    expected_disposition: ReplayDisposition
    if count == 0:
        expected_disposition = "unassessable-no-candidates"
    elif selected_ordinal is None:
        expected_disposition = "no-call"
    else:
        expected_disposition = "called"
    if disposition != expected_disposition:
        _fail("Kestrel replay disposition differs from selected evidence")
    result = KestrelReplayResult(
        capture_sha,
        policy_sha,
        cast(ReplayDisposition, disposition),
        count,
        prefilter,
        selected,
        tuple(dispositions),
        "",
    )
    return KestrelReplayResult(
        capture_sha,
        policy_sha,
        cast(ReplayDisposition, disposition),
        count,
        prefilter,
        selected,
        tuple(dispositions),
        canonical_sha256(_document(result)),
    )


def _config(capture: KestrelCapture, policy: CallerPolicyValues) -> dict[str, object]:
    config = deepcopy(load_strict_json_object(capture.kestrel_config_json))
    for pointer, path in _POINTER_PATHS.items():
        node: dict[str, object] = config
        for part in path[:-1]:
            child = node.get(part)
            if not isinstance(child, dict):
                _fail(f"Kestrel replay config is missing policy pointer {pointer}")
            node = child
        node[path[-1]] = policy.values[pointer]
    return cast(dict[str, object], config)


def _raw_frame(capture: KestrelCapture) -> pd.DataFrame:
    rows = [
        {
            "Motifs": row.motifs,
            "Variant": row.variant,
            "POS": row.position,
            "REF": row.reference_allele,
            "ALT": row.alternate_allele,
            "Sample": row.sample,
            "Motif_sequence": row.motif_sequence,
            _SOURCE_ORDINAL: row.source_row_ordinal,
        }
        for row in capture.rows
    ]
    return pd.DataFrame(rows, columns=[*KESTREL_RAW_COLUMNS, _SOURCE_ORDINAL])


def _motif_frame(capture: KestrelCapture) -> pd.DataFrame:
    rows = [{"Motif": row.motif, "Motif_sequence": row.motif_sequence} for row in capture.motifs]
    return pd.DataFrame(rows, columns=KESTREL_MOTIF_COLUMNS)


def _source_dispositions(
    capture: KestrelCapture,
    prefilter: pd.DataFrame,
    selected: pd.DataFrame,
) -> tuple[KestrelSourceDisposition, ...]:
    if not capture.rows:
        return ()
    if _SOURCE_ORDINAL not in prefilter or len(prefilter) != len(capture.rows):
        _fail("Kestrel replay did not retain the complete source population; recapture is required")
    by_ordinal = {int(row[_SOURCE_ORDINAL]): row for _, row in prefilter.iterrows()}
    if set(by_ordinal) != set(range(len(capture.rows))):
        _fail("Kestrel replay source ordinals differ from the complete capture")
    selected_ordinals = set()
    if not selected.empty:
        selected_ordinals = {int(value) for value in selected[_SOURCE_ORDINAL]}
    result = []
    for ordinal in range(len(capture.rows)):
        row = by_ordinal[ordinal]
        gates = []
        for gate in capture.selection.final_filter_columns:
            value = _cell(row[gate])
            if type(value) is not bool:
                _fail("Kestrel replay production gate is not boolean")
            if not value:
                gates.append(gate)
        result.append(KestrelSourceDisposition(ordinal, tuple(gates), ordinal in selected_ordinals))
    return tuple(result)


def replay_kestrel_capture(
    capture: KestrelCapture,
    policy: CallerPolicyValues,
    *,
    capture_policy_sha256: str,
) -> KestrelReplayResult:
    """Replay exact R3 policy values through the shared production evaluator.

    Args:
        capture: Complete validated post-VCF evidence.
        policy: Complete validated candidate policy values.
        capture_policy_sha256: Requested recruitment/assembly commitment.

    Returns:
        Immutable complete prefilter evidence, selected row and dispositions.

    Raises:
        ValueError: If capture/policy is forged, incompatible, or requires recapture.
    """
    kestrel_capture_document(capture)
    caller_policy_values_document(policy)
    if policy.required_callers != capture.baseline_policy.required_callers:
        _fail("Kestrel replay candidate caller set differs from the captured baseline")
    if _digest(capture_policy_sha256, "capture_policy_sha256") != capture.provenance.capture_policy_sha256:
        _fail("Kestrel recruitment or assembly policy changed; recapture is required")
    config = _config(capture, policy)
    evaluation = evaluate_kestrel_candidates(
        _raw_frame(capture),
        _motif_frame(capture),
        config,
        selection=capture.selection,
        add_haplo_count_fn=add_haplo_count,
        select_single_best_variant_fn=lambda frame, selection: select_single_best_variant(frame, selection=selection),
        identity_component=capture.identity_component,
    )
    if capture.rows and not evaluation.reached_final_filter:
        _fail("Kestrel replay could not evaluate the complete population; recapture is required")
    prefilter = evaluation.prefilter
    selected = evaluation.selected
    dispositions = _source_dispositions(capture, prefilter, selected)
    disposition: ReplayDisposition
    if not capture.rows:
        disposition = "unassessable-no-candidates"
    elif selected.empty:
        disposition = "no-call"
    else:
        disposition = "called"
    document = {
        "schema_version": "calibration-kestrel-replay-v1",
        "capture_sha256": capture.sha256,
        "policy_sha256": policy.sha256,
        "disposition": disposition,
        "source_row_count": len(capture.rows),
        "prefilter": _table_document(prefilter),
        "selected": _table_document(selected),
        "source_dispositions": [
            {
                "source_row_ordinal": row.source_row_ordinal,
                "blocking_gates": list(row.blocking_gates),
                "selected": row.selected,
            }
            for row in dispositions
        ],
    }
    return decode_kestrel_replay_result(document)


def _require_result(result: KestrelReplayResult) -> KestrelReplayResult:
    if not isinstance(result, KestrelReplayResult):
        _fail("Kestrel replay boundary requires KestrelReplayResult")
    decoded = decode_kestrel_replay_result(_document(result))
    if decoded != result:
        _fail("Kestrel replay result differs from its canonical content or digest")
    return decoded


def kestrel_replay_document(result: KestrelReplayResult) -> dict[str, object]:
    """Project a validated replay result as fresh canonical JSON content.

    Args:
        result: Previously decoded immutable replay output.

    Returns:
        Fresh closed replay document.

    Raises:
        ValueError: If typed content or its canonical digest was forged.
    """
    return _document(_require_result(result))


def kestrel_replay_prefilter_frame(result: KestrelReplayResult) -> pd.DataFrame:
    """Return a fresh complete prefilter frame from a validated replay result.

    Args:
        result: Previously decoded immutable replay output.

    Returns:
        Complete annotated candidates before the final six gates.

    Raises:
        ValueError: If typed content or its canonical digest was forged.
    """
    return _frame(_require_result(result).prefilter_json)


def kestrel_replay_selected_frame(result: KestrelReplayResult) -> pd.DataFrame:
    """Return a fresh selected frame from a validated replay result.

    Args:
        result: Previously decoded immutable replay output.

    Returns:
        Zero or one production-selected candidate row.

    Raises:
        ValueError: If typed content or its canonical digest was forged.
    """
    return _frame(_require_result(result).selected_json)
