"""The offline HTML page of a derived cutoff decision.

:mod:`calibration_cutoff_report` owns the files that reach disk; this module owns what the
page says and how it draws it. The page embeds its own stylesheet and draws its curves as
inline SVG, and :func:`render_cutoff_report_html` never emits a ``<script>``, a ``<link>``
or an ``@import``: a report read from a private directory on an offline machine must render
completely or it is not a report.

A curve published as unavailable (its no-call set changes across the axis on the
either-caller union) is rendered as its reason, never as an empty plot.
"""

from __future__ import annotations

import html
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Final, NoReturn

logger = logging.getLogger(__name__)

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

_CAPPED_WARNING: Final[str] = (
    "At least one searched inventory was subsampled by --max-breakpoints, so distinct operating points may be missing."
)
_UNREJECTABLE_NOTE: Final[str] = (
    "{n} sample(s) have a legacy p-value of exactly 0; no admissible cutoff can reject them."
)
#: What the held-out floor warning may say about training: a fold that fell back to the
#: baseline never met the floor, so the floor is said to be met only where a fold selected.
FLOOR_MET_EVERY_FOLD: Final[str] = "The floor was met on training folds only."
FLOOR_MET_SOME_FOLDS: Final[str] = "The floor was met on the training data of {met} of {total} folds only."
FOLDS_WITHOUT_A_CONSTRAINED_CANDIDATE: Final[str] = (
    "In {n} of {total} folds no candidate met the constraint on training data; those folds used the baseline."
)
FOLDS_WITHOUT_A_TRUTH_CLASS: Final[str] = (
    "In {n} of {total} folds a truth class was missing from the training data; those folds used the baseline."
)
_NOT_REPLAYED: Final[str] = "not replayed"
_NOT_CHECKED: Final[str] = "not checked"
_NO_CURVE: Final[str] = "No ROC/PR curve: "
_NO_BOUNDARY_SUPPORT: Final[str] = "No boundary support: "

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


def document_mapping(value: object, label: str) -> Mapping[str, Any]:
    """Require one object-valued report field; ``label`` names it in the error."""
    if not isinstance(value, Mapping):
        _fail(f"cutoff report {label} must be an object")
    return value


def document_sequence(value: object, label: str) -> Sequence[Any]:
    """Require one array-valued report field; ``label`` names it in the error."""
    if isinstance(value, str) or not isinstance(value, Sequence):
        _fail(f"cutoff report {label} must be a sequence")
    return value


def cutoff_table(document: Mapping[str, Any]) -> tuple[tuple[str, ...], list[list[object]]]:
    """Flatten every tested cutoff into one joinable row, counts and bound included.

    The page table and ``cutoffs.tsv`` share this one flattening, so they cannot disagree.

    Args:
        document: A ``calibration-cutoff-report-v1`` object.

    Returns:
        The column names and one row per tested cutoff, in published order.

    Raises:
        ValueError: If the cutoff rows are missing or malformed.
    """
    columns = ("policy_id", "axis", "value", "policy_sha256", *_COUNT_COLUMNS, "fpr_one_sided_upper", "parameters")
    rows: list[list[object]] = []
    for entry in document_sequence(document.get("cutoffs"), "cutoffs"):
        row = document_mapping(entry, "cutoff row")
        counts = document_mapping(row.get("counts"), "cutoff counts")
        metrics = document_mapping(row.get("metrics"), "cutoff metrics")
        parameters = document_mapping(row.get("parameters"), "cutoff parameters")
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
    """The ROC and PR panels for one axis, drawn from its own published points, or why there are none."""
    if curve.get("status") == "unavailable":
        return (
            f"<h2>Axis {_escape(curve.get('axis'))}</h2>"
            f"<div class='warn'>{_escape(_NO_CURVE)}{_escape(curve.get('reason'))}.</div>"
        )
    points = [document_mapping(item, "curve point") for item in document_sequence(curve.get("points"), "curve points")]
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
    selection = document_mapping(document.get("selection"), "selection")
    infeasible = selection.get("infeasible")
    if infeasible is not None:
        failed = document_mapping(infeasible, "infeasible selection")
        best = document_mapping(failed.get("best_achievable"), "best achievable")
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
    body += (
        "<p class='note'>This policy was selected on the full searched cohort and is the one exported. "
        "Its full-data counts are descriptive of that search, not validated performance; the held-out "
        "estimate above is the performance estimate.</p>"
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


def _rate_cell(rate: object) -> str:
    """One held-out rate as ``events/total = estimate (95% CI lower-upper)``."""
    if not isinstance(rate, Mapping) or rate.get("estimate") is None:
        return "undefined (no eligible samples)"
    return (
        f"{rate.get('events')}/{rate.get('total')} = {float(rate['estimate']):.3f} "
        f"(95% CI {float(rate.get('lower') or 0.0):.3f}-{float(rate.get('upper') or 0.0):.3f})"
    )


def _training_floor_statement(folds: Sequence[Mapping[str, Any]]) -> str:
    """What the folds' training data say about a requested floor, fold fallbacks included."""
    total = len(folds)
    reasons = [fold.get("fallback_reason") for fold in folds]
    met = reasons.count(None)
    unmet = reasons.count("no-candidate-satisfies-constraints")
    sentences = []
    if met == total:
        sentences.append(FLOOR_MET_EVERY_FOLD)
    elif met:
        sentences.append(FLOOR_MET_SOME_FOLDS.format(met=met, total=total))
    if unmet:
        sentences.append(FOLDS_WITHOUT_A_CONSTRAINED_CANDIDATE.format(n=unmet, total=total))
    if total - met - unmet:
        sentences.append(FOLDS_WITHOUT_A_TRUTH_CLASS.format(n=total - met - unmet, total=total))
    return " ".join(sentences)


def _held_out_section(document: Mapping[str, Any]) -> str:
    """Pooled held-out performance of the fold-selected policies, before any full-data number."""
    evaluation = document_mapping(document.get("evaluation"), "evaluation")
    held = evaluation.get("held_out")
    heading = "<h2>Held-out performance (cross-validated)</h2>"
    if not isinstance(held, Mapping):
        return (
            f"{heading}<div class='warn'>No held-out estimate is available: "
            f"{_escape(evaluation.get('status_reason'))}. Only descriptive full-data points follow.</div>"
        )
    counts = document_mapping(held.get("counts"), "held-out counts")
    exact = document_mapping(held.get("exact"), "held-out exact metrics")
    objective = document_mapping(document.get("objective"), "objective")
    folds = [document_mapping(item, "fold") for item in document_sequence(evaluation.get("folds"), "evaluation folds")]
    training = _escape(_training_floor_statement(folds))
    warnings = ""
    for label, floor_name, rate_name in (
        ("specificity", "min_specificity", "specificity"),
        ("sensitivity", "min_sensitivity", "sensitivity"),
    ):
        floor, reached = objective.get(floor_name), counts.get(rate_name)
        if floor is not None and reached is not None and float(reached) < float(floor):
            warnings += (
                f"<div class='warn'>Held-out {label} {float(reached):.3f} is below the requested floor "
                f"{_escape(floor)}. {training}</div>"
            )
    fold_table = _table(
        ("fold", "used policy", "admissible candidates", "training", "held out", "fallback"),
        [
            [
                fold.get("fold"),
                fold.get("used_policy"),
                fold.get("admissible_candidates"),
                len(document_sequence(fold.get("training_keys"), "fold training keys")),
                len(document_sequence(fold.get("held_out_keys"), "fold held-out keys")),
                fold.get("fallback_reason"),
            ]
            for fold in folds
        ],
    )
    return (
        heading + "<p class='note'>Each outer fold selected a policy on its training samples only and applied it to "
        "its held-out samples; the pooled held-out calls below are the performance estimate of this "
        "procedure. The 95% intervals are exact binomial intervals for the pooled calls and exclude "
        "selection uncertainty.</p>"
        + warnings
        + _facts(
            [
                ("true positives", counts.get("true_positives")),
                ("false negatives", counts.get("false_negatives")),
                ("true negatives", counts.get("true_negatives")),
                ("false positives", counts.get("false_positives")),
                ("no-calls", counts.get("no_calls")),
                ("sensitivity", _rate_cell(exact.get("sensitivity"))),
                ("specificity", _rate_cell(exact.get("specificity"))),
                ("fold admissibility", evaluation.get("fold_admissibility")),
            ]
        )
        + fold_table
    )


def _profile_section(document: Mapping[str, Any]) -> str:
    profile = document_mapping(document.get("profile"), "profile")
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


def _boundary_section(boundary: Mapping[str, Any]) -> str:
    """The selected axis's boundary support, or the reason the axis has none."""
    heading = "<h2>Boundary support</h2>"
    if boundary.get("status") == "unavailable":
        return f"{heading}<div class='warn'>{_escape(_NO_BOUNDARY_SUPPORT)}{_escape(boundary.get('reason'))}.</div>"
    warnings = "".join(
        f"<div class='warn'>{_escape(text)}</div>" for text in document_sequence(boundary.get("warnings"), "warnings")
    )
    return (
        heading
        + warnings
        + _facts(
            [
                ("positives inside the tested band", boundary.get("positives_within_band")),
                ("negatives inside the tested band", boundary.get("negatives_within_band")),
                ("band low", boundary.get("band_low")),
                ("band high", boundary.get("band_high")),
            ]
        )
    )


def _search_warnings(document: Mapping[str, Any]) -> str:
    """What the search could not see: a capped inventory, and samples no cutoff can reject."""
    axes = [document_mapping(item, "axis") for item in document_sequence(document.get("axes"), "axes")]
    warnings = ""
    if any(axis.get("breakpoint_completeness") == "capped-subsample" for axis in axes):
        warnings += f"<div class='warn'>{_escape(_CAPPED_WARNING)}</div>"
    for axis in axes:
        unrejectable = axis.get("unrejectable_samples")
        if unrejectable is None:
            continue
        if isinstance(unrejectable, bool) or not isinstance(unrejectable, int) or unrejectable < 0:
            _fail(f"cutoff report axis {axis.get('axis')} unrejectable_samples must be a non-negative integer")
        if unrejectable > 0:
            warnings += f"<div class='warn'>{_escape(_UNREJECTABLE_NOTE.format(n=unrejectable))}</div>"
    return warnings


def _advntr_facts(document: Mapping[str, Any], parity: Mapping[str, Any]) -> list[tuple[str, object]]:
    """The adVNTR parity and replay-consistency lines, plus the tool and wall times when it ran."""
    advntr_parity = parity.get("advntr")
    consistency = document.get("replay_consistency")
    facts: list[tuple[str, object]] = [
        (
            "adVNTR baseline parity",
            document_mapping(advntr_parity, "adVNTR baseline parity").get("proven") if advntr_parity else _NOT_REPLAYED,
        ),
        (
            "adVNTR replay consistency (candidates checked)",
            document_mapping(consistency, "replay consistency").get("checked_candidates")
            if consistency
            else _NOT_CHECKED,
        ),
    ]
    provenance = document_mapping(document.get("provenance"), "provenance")
    advntr = provenance.get("advntr")
    if advntr is None:
        return facts
    record = document_mapping(advntr, "adVNTR provenance")
    tool = document_mapping(record.get("tool_identity"), "adVNTR tool identity")
    timings = document_mapping(document.get("timings"), "timings")
    return [
        *facts,
        ("adVNTR package version", tool.get("package_version")),
        ("adVNTR build id", tool.get("build_id")),
        ("adVNTR probe grid seconds", timings.get("advntr_probe_seconds")),
        ("adVNTR candidate grid seconds", timings.get("advntr_main_seconds")),
    ]


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
    objective = document_mapping(document.get("objective"), "objective")
    scope = document_mapping(document.get("search_scope"), "search scope")
    truth = document_mapping(document.get("truth_set"), "truth set")
    provenance = document_mapping(document.get("provenance"), "provenance")
    parity = document_mapping(document.get("baseline_parity"), "baseline parity")
    boundary = document_mapping(document.get("boundary_support"), "boundary support")
    columns, rows = cutoff_table(document)
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
                (
                    "searched callers",
                    ", ".join(
                        str(name) for name in document_sequence(scope.get("searched_callers"), "searched callers")
                    ),
                ),
                (
                    "searched axes",
                    ", ".join(str(name) for name in document_sequence(scope.get("searched_axes"), "axes")),
                ),
                ("adVNTR policy", scope.get("advntr_policy")),
                ("folds requested", document.get("folds_requested")),
                ("seed", document.get("seed")),
            ]
        ),
        f"<p class='note'>{_escape(scope.get('note'))}</p>",
        _search_warnings(document),
        _held_out_section(document),
        _selection_section(document),
        "<h2>Old versus derived</h2>",
        _table(
            ("pointer", "baseline value", "derived value", "changed"),
            [
                [row.get("pointer"), row.get("baseline_value"), row.get("derived_value"), row.get("changed")]
                for row in (
                    document_mapping(item, "old-versus-derived row")
                    for item in document_sequence(document.get("old_versus_derived"), "old versus derived")
                )
            ],
        ),
        _boundary_section(boundary),
        "".join(
            _curve_figures(document_mapping(curve, "curve"))
            for curve in document_sequence(document.get("curves"), "curves")
        ),
        "<h2>Every tested cutoff (descriptive, full searched cohort)</h2>",
        "<p class='note'>Descriptive searched-cohort points: every row is scored on the same samples the "
        "search ran on, so these counts and intervals are not validated performance. Use the held-out "
        "section for that.</p>",
        _table(columns, rows, key_columns=2),
        "<h2>Rejected candidate values</h2>",
        _table(
            ("axis", "value", "reason"),
            [
                [row.get("axis"), row.get("value"), row.get("reason")]
                for row in (
                    document_mapping(item, "rejected row")
                    for item in document_sequence(document.get("rejected_candidates"), "rejected candidates")
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
                *_advntr_facts(document, parity),
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


__all__ = [
    "cutoff_table",
    "document_mapping",
    "document_sequence",
    "render_cutoff_report_html",
]
