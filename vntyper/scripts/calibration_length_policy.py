"""Pure length calibration protocol, measurement policy and source truth validation."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from typing import NoReturn

from vntyper.scripts.calibration_length_metrics import LengthEligibleRoster
from vntyper.scripts.calibration_length_protocol import (
    LengthProtocol,
    decode_length_protocol,
    length_protocol_document,
)
from vntyper.scripts.calibration_target_contract import (
    LengthBaselinePlan,
    TargetStudy,
    target_study_document,
)
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION

logger = logging.getLogger(__name__)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def length_measurement_policy_sha256(study: TargetStudy) -> str:
    """Derive the complete geometry/counting policy identity for length runs.

    Args:
        study: Frozen length study with a declared measurement baseline.

    Returns:
        Canonical hash of the annotation and counting policy, independent of fit.

    Raises:
        ValueError: If the study does not declare a length baseline.
    """
    target_study_document(study)
    if not isinstance(study.baseline, LengthBaselinePlan):
        _fail("length extraction requires a length study")
    return canonical_sha256(
        {
            "schema_version": "length-measurement-policy-v1",
            "annotation_sha256": study.baseline.annotation_sha256,
            "counting_policy_sha256": study.baseline.counting_policy_sha256,
        }
    )


def decode_length_source_truth(value: object, keys: tuple[str, ...]) -> dict[str, float]:
    """Decode and validate exact source truth repeat counts against eligible keys.

    Args:
        value: Raw JSON document containing calibration-length-truth-v1 rows.
        keys: Declared roster of eligible specimen artifact keys.

    Returns:
        Dictionary mapping specimen keys to their positive integral float counts.

    Raises:
        ValueError: On malformed fields, non-integral values, or mismatched keys.
    """
    fields = {"schema_version", "boundary_definition", "rows"}
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail("length source truth fields differ")
    if (
        value["schema_version"] != "calibration-length-truth-v1"
        or value["boundary_definition"] != TARGET_BOUNDARY_DEFINITION
    ):
        _fail("length source truth requires the declared exact repeat-count boundary")
    raw_rows = value["rows"]
    if not isinstance(raw_rows, list):
        _fail("length source truth rows must be a list")
    result: dict[str, float] = {}
    for row in raw_rows:
        if not isinstance(row, Mapping) or set(row) != {"key", "total_repeat_count"}:
            _fail("length source truth row fields differ")
        key, number = row["key"], row["total_repeat_count"]
        if not isinstance(key, str) or key in result:
            _fail("length source truth contains an invalid or duplicate key")
        if isinstance(number, bool) or not isinstance(number, (int, float)):
            _fail("length source truth requires positive integral counts")
        try:
            numeric = float(number)
        except OverflowError:
            _fail("length source truth count exceeds numeric range")
        if not math.isfinite(numeric) or numeric <= 0 or not numeric.is_integer():
            _fail("length source truth requires positive integral counts")
        result[key] = numeric
    if tuple(result) != keys:
        _fail("length source truth does not match the exact eligible artifact roster")
    return result


def length_protocol_for_roster(protocol: LengthProtocol, roster: LengthEligibleRoster) -> LengthProtocol:
    """Derive the frozen protocol for one separately declared eligible roster.

    Args:
        protocol: Study protocol whose scientific policy remains unchanged.
        roster: Role-specific eligible population frozen before outcomes.

    Returns:
        A strict protocol differing only in its eligible-roster commitment.
    """
    raw = length_protocol_document(protocol)
    raw["eligible_roster_sha256"] = roster.sha256
    return decode_length_protocol(raw)
