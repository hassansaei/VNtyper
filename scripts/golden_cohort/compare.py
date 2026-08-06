"""Diff two sides, artefact by artefact, and report it two ways.

The JSON result is the machine-readable one: every artefact of every case, with its status
and - when it differs - what changed, down to the row key and the column. The text result
is the gate page's own shape, a "Compared / Cases with a delta" table, so a run can be
written up without re-deriving anything.

An artefact present on one side and absent on the other is a delta, not a skip. A gate
whose comparison quietly narrows when a file stops being written reports a pass it did not
measure.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from golden_cohort import HARNESS_VERSION, artifacts

logger = logging.getLogger(__name__)

#: Artefact name -> (comparison kind, key columns). The kinds are ``scalar`` (compared with
#: ``==``), ``table`` (row set keyed on the named columns), ``sequence`` (ordered list of
#: strings) and ``opaque`` (any JSON value, compared structurally).
PIPELINE_ARTIFACTS: dict[str, tuple[str, tuple[str, ...]]] = {
    "exit_code": ("scalar", ()),
    "kestrel_result": ("table", artifacts.KESTREL_KEY),
    "kestrel_pre_result": ("table", artifacts.KESTREL_KEY),
    "advntr_result": ("table", artifacts.ADVNTR_KEY),
    "coverage_summary": ("table", ()),
    "screening_summary": ("opaque", ()),
    "cross_match_summary": ("opaque", ()),
    "report_tables": ("opaque", ()),
    "pipeline_steps": ("sequence", ()),
    "pipeline_step_records": ("opaque", ()),
    "executed_commands": ("sequence", ()),
}

#: The cohort artefacts, which no run of this gate before now has compared at all.
COHORT_ARTIFACTS: dict[str, tuple[str, tuple[str, ...]]] = {
    "exit_code": ("scalar", ()),
    "cohort_tables": ("opaque", ()),
    "cohort_category_counts": ("sequence", ()),
    "cohort_category_totals": ("sequence", ()),
    "cohort_kestrel_csv": ("table", ("Sample",)),
    "cohort_kestrel_tsv": ("table", ("Sample",)),
    "cohort_kestrel_json": ("opaque", ()),
    "cohort_advntr_csv": ("table", ("Sample",)),
    "cohort_advntr_tsv": ("table", ("Sample",)),
    "cohort_advntr_json": ("opaque", ()),
    "pseudonymization_table": ("table", ("Pseudonym",)),
    "cohort_output_files": ("sequence", ()),
}

#: Recorded on both sides and deliberately **not** compared, with why. ``vntyper cohort``
#: discovers samples into a set, so this order is not reproducible even between two runs of
#: the same code; comparing it would manufacture a delta on every cohort case.
UNCOMPARED: dict[str, str] = {
    "cohort_sample_order_raw": (
        "cohort sample order comes out of a set and is not reproducible between processes "
        "(tests/unit/test_cohort_inputs.py::test_discovery_returns_an_unordered_set_today); "
        "each side's order is recorded here so the defect stays visible"
    ),
    "launch_line": "the resolution proof, which names its own tree and therefore differs by construction",
}


def _row_key(row: dict[str, str], key_columns: tuple[str, ...]) -> str:
    """Build the comparison key for one table row.

    Args:
        row: The row.
        key_columns: The configured key. Columns the row does not have are skipped, and an
            empty result falls back to the whole row, so a table whose schema changed is
            still compared rather than silently keyed on nothing.

    Returns:
        str: A stable key.
    """
    present = [column for column in key_columns if column in row]
    if not present:
        return json.dumps(row, sort_keys=True)
    return json.dumps({column: row[column] for column in present}, sort_keys=True)


def diff_table(before: Any, after: Any, key_columns: tuple[str, ...]) -> dict[str, Any]:
    """Compare two parsed tables as keyed row sets.

    Args:
        before: The baseline table, or None.
        after: The candidate table, or None.
        key_columns: The columns rows are keyed on.

    Returns:
        dict[str, Any]: ``status`` plus, when it differs, ``rows_added``, ``rows_removed``,
        ``cells_changed``, ``columns_added`` and ``columns_removed``.
    """
    presence = _presence(before, after)
    if presence is not None:
        return presence

    before_rows = {_row_key(row, key_columns): row for row in before["rows"]}
    after_rows = {_row_key(row, key_columns): row for row in after["rows"]}

    added = sorted(set(after_rows) - set(before_rows))
    removed = sorted(set(before_rows) - set(after_rows))
    cells: list[dict[str, Any]] = []
    for key in sorted(set(before_rows) & set(after_rows)):
        left, right = before_rows[key], after_rows[key]
        cells.extend(
            {"key": key, "column": column, "before": left.get(column), "after": right.get(column)}
            for column in sorted(set(left) | set(right))
            if left.get(column) != right.get(column)
        )

    columns_added = sorted(set(after["columns"]) - set(before["columns"]))
    columns_removed = sorted(set(before["columns"]) - set(after["columns"]))
    provenance_changed = before.get("provenance") != after.get("provenance")

    same = not (added or removed or cells or columns_added or columns_removed)
    detail: dict[str, Any] = {
        "status": "same" if same else "differ",
        "row_count_before": len(before_rows),
        "row_count_after": len(after_rows),
    }
    if not same:
        detail.update(
            {
                "rows_added": [after_rows[key] for key in added],
                "rows_removed": [before_rows[key] for key in removed],
                "cells_changed": cells,
                "columns_added": columns_added,
                "columns_removed": columns_removed,
            }
        )
    if provenance_changed:
        detail["provenance_before"] = before.get("provenance")
        detail["provenance_after"] = after.get("provenance")
    return detail


def _presence(before: Any, after: Any) -> dict[str, Any] | None:
    """Classify the case where one or both sides did not produce the artefact.

    Args:
        before: The baseline value.
        after: The candidate value.

    Returns:
        dict[str, Any] | None: The status mapping, or None when both sides have a value.
    """
    if before is None and after is None:
        return {"status": "absent_both"}
    if before is None:
        return {"status": "added_after"}
    if after is None:
        return {"status": "removed_after"}
    return None


def diff_sequence(before: Any, after: Any) -> dict[str, Any]:
    """Compare two ordered lists of strings, reporting order and membership separately.

    A command stream that gains a step and one that merely reorders are different findings,
    and run 3's account of the assembly guard's extra ``samtools view -H`` reads turns on
    exactly that distinction.

    Args:
        before: The baseline list, or None.
        after: The candidate list, or None.

    Returns:
        dict[str, Any]: ``status`` plus counts, ``only_before``, ``only_after`` and whether
        the difference is order alone.
    """
    presence = _presence(before, after)
    if presence is not None:
        return presence

    only_before = sorted(_multiset_difference(before, after))
    only_after = sorted(_multiset_difference(after, before))
    order_only = not only_before and not only_after and list(before) != list(after)
    same = list(before) == list(after)
    detail: dict[str, Any] = {
        "status": "same" if same else "differ",
        "count_before": len(before),
        "count_after": len(after),
    }
    if not same:
        detail.update({"only_before": only_before, "only_after": only_after, "order_only": order_only})
    return detail


def _multiset_difference(left: list[str], right: list[str]) -> list[str]:
    """Items in ``left`` beyond the number of times they occur in ``right``.

    Args:
        left: The list to take items from.
        right: The list to subtract.

    Returns:
        list[str]: The surplus items, one entry per surplus occurrence.
    """
    remaining = list(right)
    surplus = []
    for item in left:
        if item in remaining:
            remaining.remove(item)
        else:
            surplus.append(item)
    return surplus


def diff_opaque(before: Any, after: Any) -> dict[str, Any]:
    """Compare two JSON-shaped values structurally.

    Args:
        before: The baseline value, or None.
        after: The candidate value, or None.

    Returns:
        dict[str, Any]: ``status`` and, when it differs, both values.
    """
    presence = _presence(before, after)
    if presence is not None:
        return presence
    if before == after:
        return {"status": "same"}
    return {"status": "differ", "before": before, "after": after}


def diff_scalar(before: Any, after: Any) -> dict[str, Any]:
    """Compare two scalars.

    Args:
        before: The baseline value.
        after: The candidate value.

    Returns:
        dict[str, Any]: ``status`` and both values.
    """
    if before == after:
        return {"status": "same", "value": before}
    return {"status": "differ", "before": before, "after": after}


def diff_case(
    before: dict[str, Any], after: dict[str, Any], spec: dict[str, tuple[str, tuple[str, ...]]]
) -> dict[str, Any]:
    """Compare every artefact of one case.

    Args:
        before: The baseline side's artefacts.
        after: The candidate side's artefacts.
        spec: The artefact specification, e.g. :data:`PIPELINE_ARTIFACTS`.

    Returns:
        dict[str, Any]: ``artefacts`` (name -> detail), ``deltas`` (the names that differ)
        and ``uncompared`` (the values recorded but deliberately not compared).
    """
    results: dict[str, Any] = {}
    for name, (kind, key_columns) in spec.items():
        left, right = before.get(name), after.get(name)
        if kind == "scalar":
            results[name] = diff_scalar(left, right)
        elif kind == "table":
            results[name] = diff_table(left, right, key_columns)
        elif kind == "sequence":
            results[name] = diff_sequence(left, right)
        else:
            results[name] = diff_opaque(left, right)

    deltas = [name for name, detail in results.items() if detail["status"] not in ("same", "absent_both")]
    uncompared = {
        name: {"before": before.get(name), "after": after.get(name), "why": why}
        for name, why in UNCOMPARED.items()
        if name in before or name in after
    }
    return {"artefacts": results, "deltas": deltas, "uncompared": uncompared}


def compare_sides(
    before_root: Path,
    after_root: Path,
    before_side: dict[str, Any],
    after_side: dict[str, Any],
    normalisation: list[dict[str, str]],
    before_rules: list[Any],
    after_rules: list[Any],
) -> dict[str, Any]:
    """Compare two completed sides case by case.

    Args:
        before_root: The baseline side's run root.
        after_root: The candidate side's run root.
        before_side: The baseline ``side.json``.
        after_side: The candidate ``side.json``.
        normalisation: The manifest to embed in the result.
        before_rules: The baseline side's normalisation rules.
        after_rules: The candidate side's normalisation rules.

    Returns:
        dict[str, Any]: The full comparison document.
    """
    cases: dict[str, Any] = {}
    matrix = artifacts.read_json(after_root / "matrix.json") or {}
    groups = {case["case_id"]: case["group"] for case in [*matrix.get("cases", []), *matrix.get("probes", [])]}

    all_pipeline_ids = sorted(set(before_side["pipeline_results"]) | set(after_side["pipeline_results"]))
    for case_id in all_pipeline_ids:
        before = artifacts.read_pipeline_case(
            before_root / "cases" / case_id, before_root / "logs" / case_id, before_rules
        )
        after = artifacts.read_pipeline_case(after_root / "cases" / case_id, after_root / "logs" / case_id, after_rules)
        cases[case_id] = {
            "kind": "pipeline",
            "group": groups.get(case_id, "unknown"),
            **diff_case(before, after, PIPELINE_ARTIFACTS),
        }

    all_cohort_ids = sorted(set(before_side.get("cohort_results", {})) | set(after_side.get("cohort_results", {})))
    for case_id in all_cohort_ids:
        before_record = before_side.get("cohort_results", {}).get(case_id, {})
        after_record = after_side.get("cohort_results", {}).get(case_id, {})
        if before_record.get("blocked") or after_record.get("blocked"):
            cases[case_id] = {
                "kind": "cohort",
                "group": "cohort",
                "artefacts": {},
                "deltas": ["BLOCKED"],
                "uncompared": {},
                "blocked": True,
                "missing_inputs": {
                    "before": before_record.get("missing_inputs", []),
                    "after": after_record.get("missing_inputs", []),
                },
            }
            continue
        before = artifacts.read_cohort_case(
            before_root / "cohorts" / case_id, before_root / "logs" / case_id, before_rules
        )
        after = artifacts.read_cohort_case(after_root / "cohorts" / case_id, after_root / "logs" / case_id, after_rules)
        cases[case_id] = {"kind": "cohort", "group": "cohort", **diff_case(before, after, COHORT_ARTIFACTS)}

    summary = _summarise(cases)
    blocked = sorted(case_id for case_id, case in cases.items() if case.get("blocked"))
    launch_ok = bool(before_side.get("launch_verified")) and bool(after_side.get("launch_verified"))

    return {
        "harness_version": HARNESS_VERSION,
        "before": {key: before_side[key] for key in ("side", "tree", "marker", "expect_marker", "launch_verified")},
        "after": {key: after_side[key] for key in ("side", "tree", "marker", "expect_marker", "launch_verified")},
        "matrix_check": matrix.get("check"),
        "normalisation": normalisation,
        "uncompared": UNCOMPARED,
        "summary": summary,
        "blocked_cases": blocked,
        "launch_verified_both_sides": launch_ok,
        "verdict": _verdict(summary, blocked, launch_ok),
        "cases": cases,
    }


def _summarise(cases: dict[str, Any]) -> dict[str, Any]:
    """Roll the per-case results up into the page's "cases with a delta" shape.

    Args:
        cases: The per-case comparison results.

    Returns:
        dict[str, Any]: Per-artefact counts and the ids that carry each delta.
    """
    per_artifact: dict[str, dict[str, Any]] = {}
    for case_id, case in cases.items():
        for name, detail in case.get("artefacts", {}).items():
            entry = per_artifact.setdefault(name, {"compared": 0, "delta": 0, "case_ids": [], "statuses": {}})
            if detail["status"] == "absent_both":
                continue
            entry["compared"] += 1
            entry["statuses"][detail["status"]] = entry["statuses"].get(detail["status"], 0) + 1
            if detail["status"] != "same":
                entry["delta"] += 1
                entry["case_ids"].append(case_id)
    for entry in per_artifact.values():
        entry["case_ids"].sort()
    return {
        "per_artifact": dict(sorted(per_artifact.items())),
        "cases_total": len(cases),
        "cases_with_any_delta": sorted(case_id for case_id, case in cases.items() if case.get("deltas")),
    }


def _verdict(summary: dict[str, Any], blocked: list[str], launch_ok: bool) -> str:
    """Turn the summary into one word.

    Args:
        summary: The rolled-up summary.
        blocked: Cohort cases that were refused.
        launch_ok: Whether every run on both sides verified its package resolution.

    Returns:
        str: ``BLOCKED``, ``UNVERIFIED``, ``IDENTICAL`` or ``DELTAS``.
    """
    if blocked:
        return "BLOCKED"
    if not launch_ok:
        return "UNVERIFIED"
    return "IDENTICAL" if not summary["cases_with_any_delta"] else "DELTAS"


def render_text(result: dict[str, Any]) -> str:
    """Render the human summary.

    Args:
        result: The comparison document.

    Returns:
        str: A Markdown report shaped like the gate page's own result tables.
    """
    lines: list[str] = []
    lines.append("# Golden-cohort gate result")
    lines.append("")
    lines.append(f"- Harness version: `{result['harness_version']}`")
    lines.append(
        f"- Before: `{result['before']['tree']}` (marker `{result['before']['marker']}` expected "
        f"{result['before']['expect_marker']})"
    )
    lines.append(
        f"- After: `{result['after']['tree']}` (marker `{result['after']['marker']}` expected "
        f"{result['after']['expect_marker']})"
    )
    lines.append(
        f"- Package resolution verified on every run, both sides: "
        f"**{'yes' if result['launch_verified_both_sides'] else 'NO'}**"
    )
    lines.append(f"- Verdict: **{result['verdict']}**")
    lines.append("")

    check = result.get("matrix_check") or {}
    if check:
        counts = check.get("counts", {})
        lines.append(
            f"Matrix derived from `tests/data`: {counts.get('total')} cases "
            f"({counts.get('base')} base, {counts.get('nonfast')} non-fast, {counts.get('advntr')} adVNTR) "
            f"plus {counts.get('probes')} probes."
        )
        if check.get("skipped"):
            lines.append("A case filter was in force, so the check against the gate page's counts is advisory.")
        lines.extend(f"- **Differs from the gate page**: {mismatch}" for mismatch in check.get("mismatches", []))
        lines.append("")

    lines.append("## Cases with a delta, per compared artefact")
    lines.append("")
    lines.append("| Compared | Cases with a delta | Cases compared |")
    lines.append("| --- | --- | --- |")
    for name, entry in result["summary"]["per_artifact"].items():
        lines.append(f"| `{name}` | {entry['delta']} | {entry['compared']} |")
    lines.append("")

    if result["blocked_cases"]:
        lines.append("## Blocked")
        lines.append("")
        for case_id in result["blocked_cases"]:
            missing = result["cases"][case_id]["missing_inputs"]
            lines.append(
                f"- `{case_id}`: inputs missing (before: {len(missing['before'])}, after: {len(missing['after'])})"
            )
        lines.append("")

    lines.append("## Per-case deltas")
    lines.append("")
    any_delta = False
    for case_id, case in sorted(result["cases"].items()):
        if not case.get("deltas"):
            continue
        any_delta = True
        lines.append(f"- `{case_id}` ({case['group']}): {', '.join(case['deltas'])}")
    if not any_delta:
        lines.append("None. Every compared artefact is identical on both sides for every case.")
    lines.append("")

    lines.append("## What was normalised")
    lines.append("")
    lines.append("A normalisation is a claim that a difference does not matter, so every one is listed.")
    lines.append("")
    lines.append("| Rule | Kind | Why |")
    lines.append("| --- | --- | --- |")
    lines.extend(f"| `{entry['name']}` | {entry['kind']} | {entry['why']} |" for entry in result["normalisation"])
    lines.append("")

    lines.append("## Recorded but not compared")
    lines.append("")
    for name, why in result["uncompared"].items():
        lines.append(f"- `{name}`: {why}")
    lines.append("")
    return "\n".join(lines)
