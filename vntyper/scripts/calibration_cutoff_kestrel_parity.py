"""Native Kestrel baseline result parity validation against replayed captures."""

from __future__ import annotations

import logging
import math
from typing import Literal, NoReturn

from vntyper.scripts.calibration_kestrel_replay import (
    KestrelReplayResult,
    kestrel_replay_document,
    kestrel_replay_selected_frame,
)
from vntyper.scripts.calibration_run_extraction import _parse_tsv, _rows
from vntyper.scripts.calibration_run_projection import is_kestrel_negative_placeholder

logger = logging.getLogger(__name__)

# Replay-only bookkeeping: the capture row ordinal exists solely inside calibration.
PRIVATE_COLUMN_PREFIX = "__Calibration_"

# The native comparison is worthless unless it covers the actual selection decision:
# the variant coordinates, the score the floor is applied to and the assigned band.
NATIVE_PARITY_REQUIRED = ("POS", "REF", "ALT", "Depth_Score", "Confidence")


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _cell_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value)


def validate_native_kestrel_baseline_parity(
    raw: bytes,
    baseline: KestrelReplayResult,
) -> Literal["called", "negative"]:
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
        if not column.startswith(PRIVATE_COLUMN_PREFIX) and column in observed
    }
    missing = tuple(column for column in NATIVE_PARITY_REQUIRED if column not in comparable)
    if missing:
        _fail(f"Kestrel cutoff native Kestrel baseline is missing required comparable columns: {', '.join(missing)}")
    if any(observed[column] != _cell_text(value) for column, value in comparable.items()):
        _fail("Kestrel cutoff native Kestrel baseline differs from capture replay selected fields")
    return "called"
