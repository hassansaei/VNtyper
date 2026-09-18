"""Private report writers for a derived cutoff decision: JSON, TSV and one offline page.

The JSON is the record; the TSVs exist so the same rows can be joined against anything
else without parsing nested objects; the HTML exists so the decision can be read without
either. All three are written into a directory that holds cohort results, so every file
is created ``0600`` and the caller is expected to have staged the directory ``0700``.

The page is deliberately self-contained. It embeds its own stylesheet and draws its curves
as inline SVG, and :func:`render_cutoff_report_html` never emits a ``<script>``, a
``<link>`` or an ``@import``. A report read from a private directory on an offline machine
must render completely or it is not a report, and a page that silently degrades to a blank
frame is worse than a TSV.

Curves and joint operating points are kept apart here exactly as they are upstream: a
curve is a sweep of one axis with everything else held fixed, while a joint point is a
labelled multi-parameter policy with no axis to plot against. They never share a file.
"""

from __future__ import annotations

import hashlib
import html
import logging
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final, NoReturn

from vntyper.scripts.canonical_json import canonical_json_bytes

logger = logging.getLogger(__name__)

SCHEMA_VERSION: Final[str] = "calibration-cutoff-report-v1"
REPORT_JSON: Final[str] = "report.json"
REPORT_HTML: Final[str] = "report.html"
CHECKSUMS_JSON: Final[str] = "checksums.json"

#: Count columns published for every tested cutoff, in the order the TSV lists them.
_COUNT_COLUMNS: Final[tuple[str, ...]] = (
    "eligible_count",
    "positive_count",
    "negative_count",
    "unknown_truth_count",
    "true_positives",
    "false_positives",
    "true_negatives",
    "false_negatives",
    "positive_no_calls",
    "negative_no_calls",
    "no_calls",
    "sensitivity",
    "specificity",
    "false_positive_rate",
    "precision",
    "f1",
    "balanced_accuracy",
    "youden_j",
    "no_call_rate",
)

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
)

_SVG_WIDTH: Final[int] = 320
_SVG_HEIGHT: Final[int] = 240
_SVG_PAD: Final[int] = 34

_STYLE: Final[str] = """
:root { color-scheme: light; }
body { font-family: ui-sans-serif, system-ui, sans-serif; margin: 0 auto; max-width: 62rem;
       padding: 1.5rem 1rem 4rem; color: #16191d; background: #fbfbfc; line-height: 1.45; }
h1 { font-size: 1.5rem; margin: 0 0 .25rem; }
h2 { font-size: 1.1rem; margin: 2rem 0 .5rem; border-bottom: 1px solid #d8dce1; padding-bottom: .25rem; }
p.lede { color: #4a525c; margin: 0 0 1.25rem; }
table { border-collapse: collapse; width: 100%; font-size: .82rem; margin: .5rem 0 0; }
th, td { border: 1px solid #d8dce1; padding: .3rem .45rem; text-align: right; }
th { background: #eef1f4; text-align: left; font-weight: 600; }
td.k, th.k { text-align: left; font-family: ui-monospace, monospace; }
dl.facts { display: grid; grid-template-columns: max-content 1fr; gap: .2rem .9rem; margin: .5rem 0 0; font-size: .88rem; }
dl.facts dt { font-weight: 600; color: #4a525c; }
dl.facts dd { margin: 0; font-family: ui-monospace, monospace; }
.warn { background: #fff4e5; border-left: 4px solid #d98324; padding: .55rem .8rem; margin: .6rem 0; font-size: .88rem; }
.fail { background: #fdecec; border-left: 4px solid #c0392b; padding: .55rem .8rem; margin: .6rem 0; font-size: .88rem; }
.note { color: #4a525c; font-size: .82rem; margin: .4rem 0 0; }
.curves { display: flex; flex-wrap: wrap; gap: 1rem; }
figure { margin: 0; }
figcaption { font-size: .8rem; color: #4a525c; margin-top: .2rem; }
svg { background: #fff; border: 1px solid #d8dce1; }
code { font-family: ui-monospace, monospace; background: #eef1f4; padding: .05rem .25rem; border-radius: 3px; }
"""


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


def _mapping(value: object, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _fail(f"cutoff report {label} must be an object")
    return value


def _sequence(value: object, label: str) -> Sequence[Any]:
    if isinstance(value, str) or not isinstance(value, Sequence):
        _fail(f"cutoff report {label} must be a sequence")
    return value


def _cutoff_rows(document: Mapping[str, Any]) -> tuple[tuple[str, ...], list[list[object]]]:
    """Flatten every tested cutoff into one joinable row, counts and bound included."""
    columns = ("policy_id", "axis", "value", "policy_sha256", *_COUNT_COLUMNS, "fpr_one_sided_upper", "parameters")
    rows: list[list[object]] = []
    for entry in _sequence(document.get("cutoffs"), "cutoffs"):
        row = _mapping(entry, "cutoff row")
        counts = _mapping(row.get("counts"), "cutoff counts")
        metrics = _mapping(row.get("metrics"), "cutoff metrics")
        parameters = _mapping(row.get("parameters"), "cutoff parameters")
        rows.append(
            [
                row.get("policy_id"),
                row.get("axis"),
                row.get("value"),
                row.get("policy_sha256"),
                *(counts.get(name) for name in _COUNT_COLUMNS),
                metrics.get("fpr_one_sided_upper"),
                ";".join(f"{pointer}={parameters[pointer]}" for pointer in sorted(parameters)),
            ]
        )
    return columns, rows


def _curve_rows(document: Mapping[str, Any]) -> list[list[object]]:
    """One row per tested threshold, tagged with the axis that produced it."""
    rows: list[list[object]] = []
    for entry in _sequence(document.get("curves"), "curves"):
        curve = _mapping(entry, "curve")
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


def _escape(value: object) -> str:
    """Every value reaches the page as text; nothing reaches it as markup."""
    return html.escape("" if value is None else str(value), quote=True)


def _facts(pairs: Sequence[tuple[str, object]]) -> str:
    items = "".join(f"<dt>{_escape(name)}</dt><dd>{_escape(value)}</dd>" for name, value in pairs)
    return f"<dl class='facts'>{items}</dl>"


def _table(columns: Sequence[str], rows: Sequence[Sequence[object]], *, key_columns: int = 1) -> str:
    head = "".join(
        f"<th class='k'>{_escape(name)}</th>" if index < key_columns else f"<th>{_escape(name)}</th>"
        for index, name in enumerate(columns)
    )
    body = "".join(
        "<tr>"
        + "".join(
            f"<td class='k'>{_escape(value)}</td>" if index < key_columns else f"<td>{_escape(value)}</td>"
            for index, value in enumerate(row)
        )
        + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _polyline(points: Sequence[tuple[float, float]]) -> str:
    """Map unit-square coordinates onto the padded plot box."""
    span_x = _SVG_WIDTH - 2 * _SVG_PAD
    span_y = _SVG_HEIGHT - 2 * _SVG_PAD
    return " ".join(f"{_SVG_PAD + x * span_x:.2f},{_SVG_HEIGHT - _SVG_PAD - y * span_y:.2f}" for x, y in points)


def _plot(title: str, x_label: str, y_label: str, points: Sequence[tuple[float, float]]) -> str:
    """One inline SVG panel: axes, a polyline and a marker per measured point."""
    body = "".join(
        f"<circle cx='{_SVG_PAD + x * (_SVG_WIDTH - 2 * _SVG_PAD):.2f}' "
        f"cy='{_SVG_HEIGHT - _SVG_PAD - y * (_SVG_HEIGHT - 2 * _SVG_PAD):.2f}' r='3' fill='#1f6feb' />"
        for x, y in points
    )
    return (
        f"<figure><svg width='{_SVG_WIDTH}' height='{_SVG_HEIGHT}' role='img' "
        f"aria-label='{_escape(title)}'>"
        f"<rect x='{_SVG_PAD}' y='{_SVG_PAD}' width='{_SVG_WIDTH - 2 * _SVG_PAD}' "
        f"height='{_SVG_HEIGHT - 2 * _SVG_PAD}' fill='none' stroke='#b6bec7' />"
        f"<polyline points='{_polyline(points)}' fill='none' stroke='#1f6feb' stroke-width='1.6' />"
        f"{body}"
        f"<text x='{_SVG_WIDTH / 2:.0f}' y='{_SVG_HEIGHT - 8}' font-size='11' text-anchor='middle' "
        f"fill='#4a525c'>{_escape(x_label)}</text>"
        f"<text x='12' y='{_SVG_HEIGHT / 2:.0f}' font-size='11' text-anchor='middle' fill='#4a525c' "
        f"transform='rotate(-90 12 {_SVG_HEIGHT / 2:.0f})'>{_escape(y_label)}</text>"
        f"</svg><figcaption>{_escape(title)}</figcaption></figure>"
    )


def _curve_figures(curve: Mapping[str, Any]) -> str:
    """The ROC and PR panels for one axis, drawn from its own published points."""
    points = [_mapping(item, "curve point") for item in _sequence(curve.get("points"), "curve points")]
    roc = [(float(point.get("false_positive_rate") or 0.0), float(point.get("sensitivity") or 0.0)) for point in points]
    pr = [
        (float(point.get("sensitivity") or 0.0), float(point.get("precision") or 0.0))
        for point in points
        if point.get("precision") is not None
    ]
    axis = _escape(curve.get("axis"))
    panels = _plot(f"ROC - {curve.get('axis')}", "false positive rate", "sensitivity", roc)
    if pr:
        panels += _plot(f"PR - {curve.get('axis')}", "sensitivity (recall)", "precision", pr)
    return (
        f"<h2>Axis {axis}</h2>"
        f"<p class='note'>statistic <code>{_escape(curve.get('statistic'))}</code>, production comparator "
        f"<code>{_escape(curve.get('comparison'))}</code>, baseline threshold "
        f"<code>{_escape(curve.get('baseline_threshold'))}</code>. Every other decision value is held fixed "
        f"under digest <code>{_escape(curve.get('fixed_policy_sha256'))}</code>.</p>"
        f"<div class='curves'>{panels}</div>"
    )


def _selection_section(document: Mapping[str, Any]) -> str:
    """What was decided, the interval it is decided over, and why it may have failed."""
    selection = _mapping(document.get("selection"), "selection")
    infeasible = selection.get("infeasible")
    if infeasible is not None:
        failed = _mapping(infeasible, "infeasible selection")
        best = _mapping(failed.get("best_achievable"), "best achievable")
        body = (
            f"<div class='fail'><strong>No cutoff satisfied the declared constraints.</strong> "
            f"Unsatisfiable: <code>{_escape(', '.join(str(name) for name in failed.get('unsatisfiable_constraints') or ()))}</code>. "
            f"{_escape(failed.get('note'))}</div>"
            + _table(
                ("constraint", "best achievable"),
                [[name, best[name]] for name in sorted(best)],
            )
        )
        return f"<h2>Selection</h2>{body}"
    plateau = selection.get("plateau")
    body = _facts(
        [
            ("selected policy", selection.get("policy_id")),
            ("axis", selection.get("axis")),
            ("selected value", selection.get("value")),
            ("reason", selection.get("reason")),
            ("policy sha256", selection.get("policy_sha256")),
        ]
    )
    if isinstance(plateau, Mapping):
        body += (
            f"<p class='note'>A threshold sweep is a step function, so the selected value is one member of an "
            f"interval, not a measurement. Values <code>{_escape(plateau.get('equivalent_values'))}</code> all "
            f"reproduce the selected outcome; the outcome changes below "
            f"<code>{_escape(plateau.get('open_below'))}</code> and above "
            f"<code>{_escape(plateau.get('open_above'))}</code>. {_escape(plateau.get('note'))}</p>"
        )
    return f"<h2>Selection</h2>{body}"


def _profile_section(document: Mapping[str, Any]) -> str:
    profile = _mapping(document.get("profile"), "profile")
    if profile.get("status") != "available":
        return (
            f"<h2>Derived research profile</h2><div class='warn'>No profile was exported: "
            f"{_escape(profile.get('reason'))}</div>"
        )
    verified = profile.get("round_trip_matches_selected_policy")
    banner = "" if verified is True else "<div class='fail'>The written profile did not round-trip.</div>"
    return (
        f"<h2>Derived research profile</h2>{banner}"
        + _facts(
            [
                ("file", profile.get("path")),
                ("profile id", profile.get("profile_id")),
                ("sha256", profile.get("sha256")),
                ("round-trips to the selected policy", profile.get("round_trip_matches_selected_policy")),
            ]
        )
        + f"<p class='note'>Apply it with <code>{_escape(document.get('usage_hint'))}</code>. "
        f"Research use only: these cutoffs were derived from a local cohort and carry no approval.</p>"
    )


def render_cutoff_report_html(document: Mapping[str, Any]) -> str:
    """Render the complete decision as one self-contained offline HTML page.

    Args:
        document: A ``calibration-cutoff-report-v1`` object.

    Returns:
        A full HTML document with an inline stylesheet and inline SVG curves, and with
        no script, stylesheet, font or image reference of any kind.

    Raises:
        ValueError: If a section the page renders is missing or malformed.
    """
    objective = _mapping(document.get("objective"), "objective")
    truth = _mapping(document.get("truth_set"), "truth set")
    provenance = _mapping(document.get("provenance"), "provenance")
    parity = _mapping(document.get("baseline_parity"), "baseline parity")
    boundary = _mapping(document.get("boundary_support"), "boundary support")
    columns, rows = _cutoff_rows(document)
    warnings = "".join(
        f"<div class='warn'>{_escape(text)}</div>" for text in _sequence(boundary.get("warnings"), "warnings")
    )
    parity_banner = (
        ""
        if parity.get("proven") is True
        else "<div class='fail'>Baseline parity was not proven for every sample.</div>"
    )
    sections = [
        "<h1>Derived caller cutoffs</h1>",
        "<p class='lede'>Research development artifact. Cutoffs below were derived from a labelled local "
        "cohort by replaying captured evidence; they carry no independent validation and no deployment "
        "approval.</p>",
        parity_banner,
        "<h2>Objective</h2>",
        _facts(
            [
                ("objective", objective.get("objective")),
                ("minimum sensitivity", objective.get("min_sensitivity")),
                ("minimum specificity", objective.get("min_specificity")),
                ("caller", document.get("caller")),
                ("folds requested", document.get("folds_requested")),
                ("seed", document.get("seed")),
            ]
        ),
        _selection_section(document),
        "<h2>Old versus derived</h2>",
        _table(
            ("pointer", "baseline value", "derived value", "changed"),
            [
                [row.get("pointer"), row.get("baseline_value"), row.get("derived_value"), row.get("changed")]
                for row in (
                    _mapping(item, "old-versus-derived row")
                    for item in _sequence(document.get("old_versus_derived"), "old versus derived")
                )
            ],
        ),
        "<h2>Boundary support</h2>",
        warnings,
        _facts(
            [
                ("positives inside the tested band", boundary.get("positives_within_band")),
                ("negatives inside the tested band", boundary.get("negatives_within_band")),
                ("band low", boundary.get("band_low")),
                ("band high", boundary.get("band_high")),
            ]
        ),
        "".join(_curve_figures(_mapping(curve, "curve")) for curve in _sequence(document.get("curves"), "curves")),
        "<h2>Every tested cutoff</h2>",
        _table(columns, rows, key_columns=2),
        "<h2>Rejected candidate values</h2>",
        _table(
            ("axis", "value", "reason"),
            [
                [row.get("axis"), row.get("value"), row.get("reason")]
                for row in (
                    _mapping(item, "rejected row")
                    for item in _sequence(document.get("rejected_candidates"), "rejected candidates")
                )
            ],
            key_columns=2,
        ),
        "<h2>Truth-set composition</h2>",
        _facts(
            [
                ("declared samples", truth.get("sample_count")),
                ("scored representatives", truth.get("primary_count")),
                ("dropped biological duplicates", truth.get("dropped_count")),
                ("by genotype", truth.get("by_genotype")),
                ("by assembly", truth.get("by_assembly")),
                ("unknown truth", truth.get("unknown_truth_count")),
            ]
        ),
        _profile_section(document),
        "<h2>Provenance</h2>",
        _facts(
            [
                ("cohort manifest sha256", provenance.get("cohort_manifest_sha256")),
                ("capture manifest sha256", provenance.get("capture_manifest_sha256")),
                ("baseline policy sha256", provenance.get("baseline_policy_sha256")),
                ("grid replay sha256", provenance.get("grid_replay_sha256")),
                ("generator version", provenance.get("generator_version")),
                ("native-exact baselines", parity.get("native_exact_count")),
            ]
        ),
        f"<p class='note'>{_escape(document.get('limitations'))}</p>",
    ]
    return (
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        "<title>Derived caller cutoffs</title>"
        f"<style>{_STYLE}</style></head><body>{''.join(sections)}</body></html>"
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
    columns, rows = _cutoff_rows(body)
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
