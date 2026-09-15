"""Exploratory comparison of actual fixed native caller policy arms."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import cast

from vntyper.scripts.calibration_caller_metrics import CallerObservation, calculate_caller_metrics
from vntyper.scripts.calibration_cohort_manifest import CohortSample
from vntyper.scripts.calibration_cohort_metrics import caller_metrics_document, group_folds, paired_caller_differences
from vntyper.scripts.calibration_run_extraction import _parse_tsv, _rows
from vntyper.scripts.calibration_run_projection import build_shipped_projection, is_kestrel_negative_placeholder
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object


def _native_rows(path: Path | None) -> list[Mapping[str, object]] | None:
    if path is None or not path.is_file():
        return None
    return _rows(_parse_tsv(read_regular_path(path), "cohort native result"), "cohort native result", allow_empty=True)


def read_native_observation(
    sample: CohortSample,
    kestrel_path: Path | None,
    advntr_path: Path | None,
    required: tuple[str, ...],
) -> CallerObservation:
    """Read native final outputs; missing, empty or nonassessable results remain no-calls.

    Args:
        sample: Declared truth and biological identity.
        kestrel_path: Actual final Kestrel TSV, if available.
        advntr_path: Actual processed adVNTR TSV, if available.
        required: Prespecified caller composition.

    Returns:
        Observed binary call with molecular identity only when native projection supports it.

    Raises:
        ValueError: If supplied TSVs are malformed or contradict native placeholders.
    """
    if required not in {("kestrel",), ("advntr",), ("kestrel", "advntr")}:
        raise ValueError("cohort caller composition is unsupported")
    kestrel = _native_rows(kestrel_path) if "kestrel" in required else []
    advntr = _native_rows(advntr_path) if "advntr" in required else []
    calls: list[bool | None] = []
    if "kestrel" in required:
        if not kestrel:
            calls.append(None)
        else:
            negative = len(kestrel) == 1 and is_kestrel_negative_placeholder(kestrel[0])
            if any(row.get("Confidence") == "Negative" for row in kestrel) and not negative:
                raise ValueError("cohort native result mixes negative and positive rows")
            positive = all(
                row.get("Confidence") in {"Low_Precision", "High_Precision", "High_Precision*"} for row in kestrel
            )
            calls.append(False if negative else True if positive else None)
    if "advntr" in required:
        if not advntr:
            calls.append(None)
        elif len(advntr) == 1 and advntr[0].get("VID") == "Negative":
            calls.append(False)
        else:
            if any(row.get("VID") == "Negative" for row in advntr):
                raise ValueError("cohort native adVNTR result mixes negative and finding rows")
            calls.append(
                True
                if all(row.get("VID") not in {None, ""} and row.get("State") not in {None, ""} for row in advntr)
                else None
            )
    call = None if None in calls else any(calls)
    identity: tuple[str, ...] = ()
    tier: tuple[str, ...] = ()
    if call:
        try:
            projection = build_shipped_projection(sample.sample_id, kestrel or [], advntr or [])
        except (KeyError, TypeError, ValueError):
            projection = {}
        value = projection.get("canonical_identity")
        if isinstance(value, str) and value:
            identity = (value,)
            tier = identity if projection.get("tier") == "A" else ()
    return CallerObservation(
        sample.sample_id,
        sample.group_id,
        sample.genotype,
        () if sample.genotype is False else None,
        call,
        identity,
        tier,
    )


def choose_policy(
    arms: Mapping[str, tuple[CallerObservation, ...]], training_keys: Sequence[str], baseline: str
) -> str:
    """Select a fixed native arm using training groups only and conservative count safeguards.

    Args:
        arms: Full native observations at each predeclared policy.
        training_keys: The only sample keys whose truth may drive selection.
        baseline: Existing unchanged policy identity.

    Returns:
        A policy with no increased false positives/no-calls and improved training sensitivity,
        otherwise baseline. This exploratory rule conveys no approval or validated safety claim.
    """
    keys = set(training_keys)
    base = calculate_caller_metrics(tuple(row for row in arms[baseline] if row.key in keys))
    if not base.sensitivity.total or not base.false_positive_rate.total:
        return baseline
    eligible = []
    for name, rows in arms.items():
        metric = calculate_caller_metrics(tuple(row for row in rows if row.key in keys))
        if (
            metric.false_positives <= base.false_positives
            and metric.no_calls <= base.no_calls
            and metric.true_positives >= base.true_positives
        ):
            eligible.append((metric.false_positives, -metric.true_positives, metric.no_calls, name != baseline, name))
    return min(eligible)[-1] if eligible else baseline


def read_policy_arms(
    samples: tuple[CohortSample, ...],
    policies_path: Path | None,
    caller_runs: Path | None,
    *,
    required_override: tuple[str, ...] | None = None,
) -> tuple[str, dict[str, tuple[CallerObservation, ...]], dict[str, dict[str, object]]]:
    """Read explicit fixed policy inventories or existing standard baseline run outputs.

    Args:
        samples: Complete declared sample roster.
        policies_path: Optional closed policy JSON with actual native output paths.
        caller_runs: Optional baseline output root containing sample_id run directories.

    Returns:
        Baseline identity, complete policy observations, and declared scalar metadata.

    Raises:
        ValueError: For incomplete/duplicate inventories, invalid cutoffs or training claims.
    """
    if policies_path is None:
        has_advntr = any(
            row.advntr_result is not None
            or (caller_runs is not None and (caller_runs / row.sample_id / "advntr/output_adVNTR_result.tsv").is_file())
            for row in samples
        )
        composition: tuple[str, ...] = ("kestrel", "advntr") if has_advntr else ("kestrel",)
        composition = required_override or composition
        rows = []
        for sample in samples:
            root = caller_runs / sample.sample_id if caller_runs else None
            kestrel = sample.kestrel_result or (root / "kestrel/kestrel_result.tsv" if root else None)
            advntr = sample.advntr_result or (root / "advntr/output_adVNTR_result.tsv" if root else None)
            rows.append(read_native_observation(sample, kestrel, advntr, composition))
        return (
            "baseline",
            {"baseline": tuple(rows)},
            {"baseline": {"cutoff": None, "comparison": None, "required_callers": list(composition)}},
        )
    document = load_strict_json_object(read_regular_path(policies_path))
    if set(document) != {"schema_version", "baseline", "required_callers", "training_scope", "policies"} or (
        document["schema_version"] != "cohort-caller-policies-v1" or document["training_scope"] != "fixed-before-cohort"
    ):
        raise ValueError("cohort policies require a closed fixed-before-cohort native policy inventory")
    required = document["required_callers"]
    if required not in [["kestrel"], ["kestrel", "advntr"]]:
        raise ValueError("cohort policies have unsupported caller composition")
    required = list(required_override) if required_override else required
    entries = document["policies"]
    if not isinstance(entries, list) or not entries:
        raise ValueError("cohort policies must be a nonempty list")
    arms: dict[str, tuple[CallerObservation, ...]] = {}
    metadata: dict[str, dict[str, object]] = {}
    for entry in entries:
        if not isinstance(entry, dict) or set(entry) != {"policy_id", "cutoff", "comparison", "samples"}:
            raise ValueError("cohort policy fields differ")
        name, cutoff, comparison = entry["policy_id"], entry["cutoff"], entry["comparison"]
        if not isinstance(name, str) or not name or name in arms:
            raise ValueError("cohort policy identities must be unique nonempty strings")
        if cutoff is not None and (type(cutoff) not in (int, float) or not math.isfinite(cutoff)):
            raise ValueError("cohort policy cutoff must be finite or null")
        if comparison not in {None, "<", "<=", ">", ">="} or ((cutoff is None) != (comparison is None)):
            raise ValueError("cohort policy cutoff and comparison must be supplied together")
        inventory = entry["samples"]
        if not isinstance(inventory, dict) or set(inventory) != {sample.sample_id for sample in samples}:
            raise ValueError("cohort policy must account for the exact sample roster")
        observations = []
        for sample in samples:
            paths = inventory[sample.sample_id]
            if not isinstance(paths, dict) or set(paths) != {"kestrel_result", "advntr_result"}:
                raise ValueError("cohort policy result path fields differ")
            resolved = []
            for key in ("kestrel_result", "advntr_result"):
                value = paths[key]
                if value is not None and (not isinstance(value, str) or not value):
                    raise ValueError("cohort policy result paths must be nonempty strings or null")
                resolved.append(None if value is None else (policies_path.parent / value).resolve())
            observations.append(read_native_observation(sample, resolved[0], resolved[1], tuple(required)))
        arms[name] = tuple(observations)
        metadata[name] = {"cutoff": cutoff, "comparison": comparison, "required_callers": required}
    baseline = document["baseline"]
    if not isinstance(baseline, str) or baseline not in arms:
        raise ValueError("cohort policy baseline must name an available arm")
    return baseline, arms, metadata


def compare_caller_arms(
    samples: tuple[CohortSample, ...],
    policies_path: Path | None,
    caller_runs: Path | None,
    *,
    folds: int,
    seed: int,
) -> dict[str, object]:
    """Compare complete fixed-policy outcomes with group-held-out policy selection.

    Args:
        samples: Declared complete cohort.
        policies_path: Explicit native policy inventory, if available.
        caller_runs: Optional baseline result root.
        folds: Requested outer folds.
        seed: Prespecified fold seed.

    Returns:
        Exploratory counts, native operating points and fold-selected held-out metrics.
    """
    # One outcome-independent primary per declared biological group for binomial metrics.
    primary: dict[str, CohortSample] = {}
    for sample in sorted(samples, key=lambda row: row.sample_id):
        primary.setdefault(sample.group_id, sample)
    representatives = tuple(primary.values())
    baseline, all_arms, metadata = read_policy_arms(samples, policies_path, caller_runs)
    keys = {row.sample_id for row in representatives}
    arms = {name: tuple(row for row in observations if row.key in keys) for name, observations in all_arms.items()}
    assignments = group_folds({row.sample_id: row.group_id for row in representatives}, folds=folds, seed=seed)
    chosen: dict[str, str] = {}
    selected: list[CallerObservation] = []
    if assignments:
        for fold in sorted(set(assignments.values())):
            training = tuple(key for key, value in assignments.items() if value != fold)
            name = choose_policy(arms, training, baseline)
            chosen[str(fold)] = name
            selected.extend(row for row in arms[name] if assignments[row.key] == fold)
    baseline_metrics = calculate_caller_metrics(arms[baseline])
    operating = [
        {"policy_id": name, **metadata[name], "metrics": caller_metrics_document(calculate_caller_metrics(rows))}
        for name, rows in arms.items()
    ]
    per_caller = {}
    for caller in cast(list[str], metadata[baseline]["required_callers"]):
        _, single_arms, _ = read_policy_arms(samples, policies_path, caller_runs, required_override=(caller,))
        single_baseline = tuple(row for row in single_arms[baseline] if row.key in keys)
        single_selected = tuple(
            row
            for fold, name in chosen.items()
            for row in single_arms[name]
            if row.key in keys and assignments[row.key] == int(fold)
        )
        per_caller[caller] = {
            "baseline": caller_metrics_document(calculate_caller_metrics(single_baseline)),
            "selected_out_of_fold": caller_metrics_document(calculate_caller_metrics(single_selected))
            if single_selected
            else None,
        }
    by_key = {row.key: row for row in selected}
    return {
        "status": "available"
        if any(row.truth_positive is not None and row.called_positive is not None for row in arms[baseline])
        else "unavailable",
        "evidence_status": "exploratory; fixed native arms, no fitted-background retraining",
        "composition": "native-positive union (OR); requires all declared callers assessable"
        if len(per_caller) > 1
        else "native " + next(iter(per_caller)),
        "per_caller": per_caller,
        "paired_differences": paired_caller_differences(arms[baseline], selected, seed=seed) if selected else None,
        "predictions": [
            {
                "sample_id": row.key,
                "group_id": row.group_key,
                "truth": row.truth_positive,
                "baseline": row.called_positive,
                "selected": by_key[row.key].called_positive if row.key in by_key else None,
                "fold": assignments.get(row.key),
            }
            for row in arms[baseline]
        ],
        "baseline_policy": baseline,
        "baseline": caller_metrics_document(baseline_metrics),
        "selected_out_of_fold": caller_metrics_document(calculate_caller_metrics(selected)) if selected else None,
        "selected_by_fold": chosen,
        "folds": assignments,
        "operating_points": operating,
        "operating_points_scope": "descriptive full cohort; never used to select the outer held-out policy",
        "roc_scope": "observed native operating points only; no interpolated curve or fabricated endpoints",
        "omitted_related_samples": [row.sample_id for row in samples if row.sample_id not in keys],
        "selection_available": len(arms) > 1 and bool(assignments),
        "identity_truth_status": "not supplied; binary mutation truth only",
    }
