"""Private report writers for a derived cutoff decision: JSON, TSV and one offline page.

The JSON is the record; the TSVs exist so the same rows can be joined against anything
else without parsing nested objects; the HTML exists so the decision can be read without
either. All three are written into a directory that holds cohort results, so every file
is created ``0600`` and the caller is expected to have staged the directory ``0700``.

The page itself is rendered by :mod:`calibration_cutoff_report_html`, which keeps it
self-contained and offline; :func:`render_cutoff_report_html` is re-exported here.

Curves and joint operating points are kept apart here exactly as they are upstream: a
curve is a sweep of one axis with everything else held fixed, while a joint point is a
labelled multi-parameter policy with no axis to plot against. They never share a file.
"""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final, NoReturn

from vntyper.scripts.calibration_cutoff_report_html import (
    cutoff_table,
    render_cutoff_report_html,
)
from vntyper.scripts.calibration_cutoff_report_html import (
    document_mapping as _mapping,
)
from vntyper.scripts.calibration_cutoff_report_html import (
    document_sequence as _sequence,
)
from vntyper.scripts.canonical_json import canonical_json_bytes

logger = logging.getLogger(__name__)

SCHEMA_VERSION: Final[str] = "calibration-cutoff-report-v1"
REPORT_JSON: Final[str] = "report.json"
REPORT_HTML: Final[str] = "report.html"
CHECKSUMS_JSON: Final[str] = "checksums.json"

_CURVE_COLUMNS: Final[tuple[str, ...]] = (
    "candidate_id",
    "threshold",
    "true_positives",
    "false_positives",
    "true_negatives",
    "false_negatives",
    "positive_no_calls",
    "negative_no_calls",
    "unknown_truth_count",
    "false_positive_rate",
    "sensitivity",
    "precision",
)

_FOLD_COLUMNS: Final[tuple[str, ...]] = (
    "fold",
    "selected_policy",
    "selection_reason",
    "fallback_reason",
    "training_count",
    "held_out_count",
    "admissible_candidates",
)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _write_private(path: Path, raw: bytes) -> None:
    """Write one artifact byte-exactly with owner-only permissions."""
    path.write_bytes(raw)
    path.chmod(0o600)


def _cell(value: object) -> str:
    """Render one scalar for a TSV cell without ever emitting a delimiter."""
    if value is None:
        return ""
    if value is True or value is False:
        return "true" if value else "false"
    text = str(value)
    if any(character in text for character in "\t\r\n"):
        _fail("cutoff report TSV values must not contain tabs or newlines")
    return text


def _tsv(path: Path, columns: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    """Write one header-plus-rows TSV privately, even when there are no rows."""
    lines = ["\t".join(columns)]
    lines.extend("\t".join(_cell(value) for value in row) for row in rows)
    _write_private(path, ("\n".join(lines) + "\n").encode("utf-8"))


def _curve_rows(document: Mapping[str, Any]) -> list[list[object]]:
    """One row per tested threshold, tagged with the axis that produced it."""
    rows: list[list[object]] = []
    for entry in _sequence(document.get("curves"), "curves"):
        curve = _mapping(entry, "curve")
        if curve.get("status") == "unavailable":
            continue
        for item in _sequence(curve.get("points"), "curve points"):
            point = _mapping(item, "curve point")
            rows.append([curve.get("axis"), curve.get("comparison"), *(point.get(name) for name in _CURVE_COLUMNS)])
    return rows


def _fold_rows(document: Mapping[str, Any]) -> list[list[object]]:
    """The per-fold training decisions, kept apart from the full-data fit."""
    evaluation = _mapping(document.get("evaluation"), "evaluation")
    rows: list[list[object]] = []
    for entry in _sequence(evaluation.get("folds"), "evaluation folds"):
        fold = _mapping(entry, "fold")
        selection = _mapping(fold.get("selection"), "fold selection")
        rows.append(
            [
                fold.get("fold"),
                fold.get("used_policy"),
                selection.get("reason"),
                fold.get("fallback_reason"),
                len(_sequence(fold.get("training_keys"), "fold training keys")),
                len(_sequence(fold.get("held_out_keys"), "fold held-out keys")),
                fold.get("admissible_candidates"),
            ]
        )
    return rows


def _joint_rows(joint: Mapping[str, Any]) -> tuple[tuple[str, ...], list[list[object]]]:
    """The labelled multi-parameter table; it carries no threshold column by design."""
    points = [_mapping(item, "joint point") for item in _sequence(joint.get("points"), "joint points")]
    metric_names: list[str] = []
    for point in points:
        for name in _mapping(point.get("metrics"), "joint metrics"):
            if name not in metric_names:
                metric_names.append(name)
    rows = [
        [
            point.get("label"),
            ";".join(
                f"{pointer}={_mapping(point.get('values'), 'joint values')[pointer]}"
                for pointer in sorted(_mapping(point.get("values"), "joint values"))
            ),
            *(_mapping(point.get("metrics"), "joint metrics").get(name) for name in metric_names),
        ]
        for point in points
    ]
    return ("label", "values", *metric_names), rows


def _checksums(output: Path) -> None:
    """Bind every direct artifact by digest, privately, as the last step."""
    files = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != CHECKSUMS_JSON
    }
    _write_private(
        output / CHECKSUMS_JSON,
        canonical_json_bytes({"schema_version": "calibration-checksums-v1", "files": files}),
    )


def write_cutoff_reports(output: Path, document: Mapping[str, Any]) -> None:
    """Write the canonical report, its joinable TSVs, the offline page and the checksums.

    Args:
        output: Already staged private directory; every file is created ``0600``.
        document: A complete ``calibration-cutoff-report-v1`` object.

    Raises:
        ValueError: If the document is not a cutoff report, a section it must publish is
            missing or malformed, or a value would break the TSV delimiter contract.
    """
    if not isinstance(output, Path):
        _fail("cutoff report output must be a Path")
    body = _mapping(document, "document")
    if body.get("schema_version") != SCHEMA_VERSION:
        _fail(f"cutoff report schema_version must be {SCHEMA_VERSION}")
    columns, rows = cutoff_table(body)
    curve_rows = _curve_rows(body)
    fold_rows = _fold_rows(body)
    rejected = [
        [row.get("axis"), row.get("value"), row.get("reason")]
        for row in (
            _mapping(item, "rejected row") for item in _sequence(body.get("rejected_candidates"), "rejected candidates")
        )
    ]
    pointers = [
        [row.get("pointer"), row.get("baseline_value"), row.get("derived_value"), row.get("changed")]
        for row in (
            _mapping(item, "old-versus-derived row")
            for item in _sequence(body.get("old_versus_derived"), "old versus derived")
        )
    ]
    page = render_cutoff_report_html(body)
    _write_private(output / REPORT_JSON, canonical_json_bytes(dict(body)))
    _tsv(output / "cutoffs.tsv", columns, rows)
    _tsv(output / "roc-pr-curves.tsv", ("axis", "comparison", *_CURVE_COLUMNS), curve_rows)
    _tsv(output / "folds.tsv", _FOLD_COLUMNS, fold_rows)
    _tsv(output / "old-versus-derived.tsv", ("pointer", "baseline_value", "derived_value", "changed"), pointers)
    _tsv(output / "rejected-candidates.tsv", ("axis", "value", "reason"), rejected)
    joint = body.get("joint_points")
    if joint is not None:
        joint_columns, joint_rows = _joint_rows(_mapping(joint, "joint points"))
        _tsv(output / "joint-points.tsv", joint_columns, joint_rows)
    _write_private(output / REPORT_HTML, page.encode("utf-8"))
    _checksums(output)


__all__ = [
    "CHECKSUMS_JSON",
    "REPORT_HTML",
    "REPORT_JSON",
    "SCHEMA_VERSION",
    "render_cutoff_report_html",
    "write_cutoff_reports",
]
