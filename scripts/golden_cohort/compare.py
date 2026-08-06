"""Diff two sides, artefact by artefact, and report it two ways.

The JSON result is the machine-readable one: every artefact of every case, with its status
and - when it differs - what changed, down to the row key and the column. The text result
is the gate page's own shape, a "Compared / Cases with a delta" table, so a run can be
written up without re-deriving anything.

An artefact present on one side and absent on the other is a delta, not a skip. A gate
whose comparison quietly narrows when a file stops being written reports a pass it did not
measure.

Absent on **both** sides is the one case this module cannot judge on its own, and it is
where the gate's worst failure lived. ``_presence`` classifies it ``absent_both``,
``diff_case`` excludes it from the delta list and ``_summarise`` excludes it from the
compared count - all correct, because a non-adVNTR case legitimately has no adVNTR output
and manufacturing 56 deltas out of that would drown the signal. But it means two runs that
both died before writing anything compared *equal*. The judgement of whether an absence is
legitimate cannot be made here; it is made per case, from the matrix's ``expect_exit`` and
``required_artifacts``, by :mod:`golden_cohort.admissibility`, and folded in as the
``EXPECTATION`` delta below.
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

#: Why cohort sample order is recorded but not compared.
#:
#: The note this replaces said cohort discovery "returns a set, so order is not
#: reproducible", and cited
#: ``tests/unit/test_cohort_inputs.py::test_discovery_returns_an_unordered_set_today``.
#: Both halves were wrong by the time they shipped. That test is not in the repository -
#: it was replaced by ``test_the_discovered_directories_come_back_sorted`` when the
#: determinism fix landed - and ``discover_sample_directories`` now ends
#: ``return sorted(processed_dirs), temp_dirs``, so for fixed input directories the order
#: *is* reproducible.
#:
#: Two narrower reasons survive, and the sort is justified by the first of them alone for
#: any comparison whose baseline predates the fix:
#:
#: 1. **Cross-version.** A baseline older than the determinism fix iterates the ``set``
#:    directly - at ``4fd638a``, run 4's baseline, ``cohort_summary.py`` reads
#:    ``for sample_dir in processed_dirs:`` over an unsorted set. ``Path.__hash__`` is the
#:    hash of the path string and Python randomises string hashing per process, so that
#:    side's order differs between two runs of *itself*. Comparing order across such a
#:    pair measures the interpreter's hash seed.
#: 2. **ZIP inputs, on any version.** Each ZIP is extracted to
#:    ``tempfile.mkdtemp(prefix="cohort_zip_")``, whose random suffix is part of the path
#:    and therefore part of the sort key. Two runs over the same ZIP sort its samples into
#:    different absolute positions.
#:
#: What this costs is stated rather than hidden: because the order is normalised away, a
#: change that fixes or breaks cohort ordering is invisible to this gate. It is attested by
#: the unit tests named above and by
#: ``test_processes_with_different_hash_seeds_discover_the_same_order``, not by any run of
#: this harness.
COHORT_ORDER_WHY = (
    "cohort sample order is not comparable across a version boundary: a baseline predating the "
    "determinism fix iterates the discovery set directly and so differs between two runs of itself, and a "
    "ZIP input extracts to a randomly-named tempdir on any version, which is part of the sort key. "
    "vntyper.scripts.cohort_inputs.discover_sample_directories does return sorted() today, pinned by "
    "tests/unit/test_cohort_inputs.py::test_the_discovered_directories_come_back_sorted and "
    "::test_processes_with_different_hash_seeds_discover_the_same_order - so the ordering fix itself is "
    "attested by those tests and NOT by this gate. Each side's raw order is recorded here uncompared."
)

#: Recorded on both sides and deliberately **not** compared, with why.
UNCOMPARED: dict[str, str] = {
    "cohort_sample_order_raw": COHORT_ORDER_WHY,
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

    The ``##`` provenance banner counts. It is compared *after* normalisation, which has
    already replaced the analysis timestamp and the VNtyper version - the two things that
    differ by construction - with placeholders, so anything still differing in it is a real
    change to what the file says about itself. This used to be computed and then discarded:
    ``provenance_changed`` was attached to the detail while ``status`` stayed ``"same"``,
    so it never became a delta and never reached the verdict. Folding it in costs nothing
    on the recorded runs - run 4 has zero provenance changes across all 118 Kestrel tables,
    3 adVNTR tables and 59 coverage tables - and it is the only way a banner that quietly
    stops naming the pipeline version would be seen.

    Args:
        before: The baseline table, or None.
        after: The candidate table, or None.
        key_columns: The columns rows are keyed on.

    Returns:
        dict[str, Any]: ``status`` plus, when it differs, ``rows_added``, ``rows_removed``,
        ``cells_changed``, ``columns_added``, ``columns_removed`` and, when the banner
        moved, ``provenance_before`` / ``provenance_after``.
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

    same = not (added or removed or cells or columns_added or columns_removed or provenance_changed)
    detail: dict[str, Any] = {
        "status": "same" if same else "differ",
        "row_count_before": len(before_rows),
        "row_count_after": len(after_rows),
        "provenance_changed": provenance_changed,
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
            **_expectation_detail(
                before_side["pipeline_results"].get(case_id), after_side["pipeline_results"].get(case_id)
            ),
        }
        _fold_expectation_into_deltas(cases[case_id])

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
        cases[case_id] = {
            "kind": "cohort",
            "group": "cohort",
            **diff_case(before, after, COHORT_ARTIFACTS),
            **_expectation_detail(before_record, after_record),
        }
        _fold_expectation_into_deltas(cases[case_id])

    summary = _summarise(cases)
    blocked = sorted(case_id for case_id, case in cases.items() if case.get("blocked"))
    launch_ok = bool(before_side.get("launch_verified")) and bool(after_side.get("launch_verified"))
    unmet = sorted(case_id for case_id, case in cases.items() if case.get("expectation_problems"))
    check = matrix.get("check") or {}
    # A matrix that predates `attestation_grade` says nothing either way, so it is treated
    # as attestation-grade: this must not retroactively downgrade runs 1-4.
    attestation_grade = bool(check.get("attestation_grade", True))

    return {
        "harness_version": HARNESS_VERSION,
        "before": _side_summary(before_side),
        "after": _side_summary(after_side),
        "matrix_check": matrix.get("check"),
        "attestation_grade": attestation_grade,
        "normalisation": normalisation,
        "uncompared": UNCOMPARED,
        "summary": summary,
        "blocked_cases": blocked,
        "expectations_unmet": unmet,
        "launch_verified_both_sides": launch_ok,
        "verdict": _verdict(summary, blocked, launch_ok, unmet, attestation_grade),
        "cases": cases,
    }


#: The fields of a side record the result document quotes back. ``revision`` is what turns
#: "this run attests commit X" from an assertion into a record; a side written by an older
#: harness has no such key and reports None rather than a plausible-looking blank.
SIDE_SUMMARY_KEYS: tuple[str, ...] = (
    "side",
    "tree",
    # The version of the harness that *produced* the side, which is not necessarily the one
    # comparing it. A 1.0.0 side enforced no case expectation and recorded no revision, so
    # a result document mixing the two must say which measured what.
    "harness_version",
    "marker",
    "expect_marker",
    "launch_verified",
    "revision",
    "cases_launched",
    "expectations_met",
)


def _side_summary(side: dict[str, Any]) -> dict[str, Any]:
    """Quote the fields of one side record the result document carries.

    Args:
        side: The side's ``side.json``.

    Returns:
        dict[str, Any]: The quoted fields, with None for anything an older harness did not
        record.
    """
    return {key: side.get(key) for key in SIDE_SUMMARY_KEYS}


def _expectation_detail(before_record: dict[str, Any] | None, after_record: dict[str, Any] | None) -> dict[str, Any]:
    """Collect what each side's runner said about whether this case did what it declared.

    Args:
        before_record: The baseline side's runner record for this case, or None.
        after_record: The candidate side's runner record for this case, or None.

    Returns:
        dict[str, Any]: ``expectation_problems`` keyed by side, empty when both sides met
        their expectation or when neither side recorded one (an older run root).
    """
    problems: dict[str, list[str]] = {}
    for side, record in (("before", before_record), ("after", after_record)):
        found = (record or {}).get("expectation_problems") or []
        if found:
            problems[side] = list(found)
    return {"expectation_problems": problems}


def _fold_expectation_into_deltas(case: dict[str, Any]) -> None:
    """Make an unmet expectation a delta, in place.

    This is the fix for the gate's worst failure. Two sides that both exit 1 without
    writing a genotype artefact produce nothing but ``absent_both`` statuses, which
    ``diff_case`` correctly excludes from the delta list - so the case, and then the run,
    came back ``IDENTICAL``. An unmet expectation is now a delta in its own right,
    independent of whether the two sides failed identically.

    Args:
        case: The per-case comparison result, mutated in place.
    """
    if case.get("expectation_problems") and "EXPECTATION" not in case.get("deltas", []):
        case.setdefault("deltas", []).insert(0, "EXPECTATION")


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


def _verdict(
    summary: dict[str, Any],
    blocked: list[str],
    launch_ok: bool,
    unmet: list[str],
    attestation_grade: bool,
) -> str:
    """Turn the summary into one word.

    The order is worst-first and is not arbitrary. ``BLOCKED`` means a case never ran;
    ``UNVERIFIED`` means it ran but the harness cannot say which package it ran;
    ``EXPECTATIONS_UNMET`` means it ran the right package and then did not do what the
    matrix said it would, which makes every "same" beneath it meaningless; only after those
    three is a difference in output the interesting finding.

    ``REDUCED`` is the fourth word and is the point of C4: a clean run over a matrix that
    is not the documented one must not be reported with the same word as a clean run over
    the documented one. It is exactly as free of deltas as ``IDENTICAL`` and attests
    strictly less, and the two used to be indistinguishable.

    Args:
        summary: The rolled-up summary.
        blocked: Cohort cases that were refused.
        launch_ok: Whether every run on both sides verified its package resolution.
        unmet: Cases where either side failed its declared exit code or artefact set.
        attestation_grade: Whether the matrix was the documented one, unfiltered.

    Returns:
        str: ``BLOCKED``, ``UNVERIFIED``, ``EXPECTATIONS_UNMET``, ``DELTAS``, ``REDUCED``
        or ``IDENTICAL``.
    """
    if blocked:
        return "BLOCKED"
    if not launch_ok:
        return "UNVERIFIED"
    if unmet:
        return "EXPECTATIONS_UNMET"
    if summary["cases_with_any_delta"]:
        return "DELTAS"
    return "IDENTICAL" if attestation_grade else "REDUCED"


def _render_revision(side: dict[str, Any]) -> str:
    """Render one side's recorded revision for the human summary.

    Args:
        side: The quoted side summary.

    Returns:
        str: The commit and its cleanliness, or an explicit statement that the harness did
        not record one - never a blank that reads as "clean".
    """
    revision = side.get("revision") or {}
    head = revision.get("head")
    if not head:
        return "**revision not recorded**"
    state = "dirty" if revision.get("dirty_relevant") else "clean"
    if revision.get("dirty_relevant") is None:
        state = "cleanliness unknown"
    return f"`{head[:12]}` ({revision.get('branch') or 'detached'}, {state})"


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
    for label in ("before", "after"):
        side = result[label]
        lines.append(
            f"- {label.capitalize()}: {_render_revision(side)} at `{side['tree']}` "
            f"(marker `{side['marker']}` expected {side['expect_marker']})"
        )
    lines.append(
        f"- Package resolution verified on every run, both sides: "
        f"**{'yes' if result['launch_verified_both_sides'] else 'NO'}**"
    )
    lines.append(f"- Verdict: **{result['verdict']}**")
    if not result.get("attestation_grade", True):
        lines.append(
            "- **Not attestation-grade**: the matrix was filtered or deviates from the documented contract, "
            "so a clean result is reported as `REDUCED` rather than `IDENTICAL`."
        )
    lines.append("")

    check = result.get("matrix_check") or {}
    if check:
        counts = check.get("counts", {})
        lines.append(
            f"Matrix: {counts.get('total')} cases ({counts.get('base')} base **derived from `tests/data`**, "
            f"{counts.get('nonfast')} non-fast and {counts.get('advntr')} adVNTR **selected by declared policy**) "
            f"plus {counts.get('probes')} probes, also policy."
        )
        if check.get("skipped"):
            lines.append("A case filter was in force, so this run is not attestation-grade.")
        lines.extend(f"- **Differs from the gate page**: {mismatch}" for mismatch in check.get("mismatches", []))
        lines.append("")

    if result.get("expectations_unmet"):
        lines.append("## Cases that did not do what the matrix declared")
        lines.append("")
        lines.append(
            "A case here did not exit as declared, or exited as declared and wrote none of the artefacts it "
            "must produce. Every comparison beneath it is comparing two absences."
        )
        lines.append("")
        for case_id in result["expectations_unmet"]:
            for side, problems in sorted(result["cases"][case_id]["expectation_problems"].items()):
                lines.extend(f"- `{case_id}` ({side}): {problem}" for problem in problems)
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
