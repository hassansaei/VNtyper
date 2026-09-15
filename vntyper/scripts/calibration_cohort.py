"""Accessible exploratory cohort fitting without fabricated confirmation authority."""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

from jinja2 import Environment, PackageLoader, select_autoescape

from vntyper.scripts.calibration_artifact_io import write_checksums, write_json
from vntyper.scripts.calibration_cohort_callers import compare_caller_arms
from vntyper.scripts.calibration_cohort_manifest import CohortSample, audit_alignment_duplicates, read_cohort_manifest
from vntyper.scripts.calibration_cohort_metrics import group_folds, length_metrics
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)


def _length_services() -> tuple[Callable[..., Any], ...]:
    # Training dependencies are absent from mutation-only and pipeline-inference paths.
    from vntyper.scripts.calibration_standard_length import fit_standard_length_model, standard_length_training_reasons
    from vntyper.scripts.length_standard_io import read_standard_length_features
    from vntyper.scripts.length_standard_model import encode_standard_length_model, predict_standard_length

    return (
        read_standard_length_features,
        fit_standard_length_model,
        predict_standard_length,
        encode_standard_length_model,
        standard_length_training_reasons,
    )


def _length_comparison(
    samples: tuple[CohortSample, ...],
    reference: Path,
    *,
    folds: int,
    seed: int,
    convention: str,
) -> tuple[dict[str, object], object | None]:
    read, fit, predict, encode, training_reasons = _length_services()
    eligible = tuple(sample for sample in samples if sample.length_total is not None)
    assignments = group_folds({sample.sample_id: sample.group_id for sample in eligible}, folds=folds, seed=seed)
    measured: dict[str, Any] = {}
    reasons: dict[str, list[str]] = {}
    for sample in eligible:
        try:
            measurement = read(sample.bam, reference, assembly=sample.assembly)
            qc_reasons = training_reasons(measurement)
            if qc_reasons:
                reasons[sample.sample_id] = list(qc_reasons)
            else:
                measured[sample.sample_id] = measurement
        except (OSError, ValueError, RuntimeError) as error:
            logger.warning("Length measurement unavailable for %s: %s", sample.sample_id, error)
            reasons[sample.sample_id] = ["measurement_failed"]
    predictions: dict[str, float | None] = {sample.sample_id: None for sample in eligible}
    baseline: dict[str, float | None] = dict(predictions)
    fold_records: list[dict[str, object]] = []
    for fold in sorted(set(assignments.values())):
        training = tuple(sample for sample in eligible if assignments[sample.sample_id] != fold)
        fitting = tuple(sample for sample in training if sample.sample_id in measured)
        heldout = tuple(sample for sample in eligible if assignments[sample.sample_id] == fold)
        targets = [sample.length_total for sample in fitting]
        # The baseline uses training truth regardless of feature availability.
        training_targets = [float(sample.length_total) for sample in training if sample.length_total is not None]
        for sample in heldout:
            baseline[sample.sample_id] = sum(training_targets) / len(training_targets) if training_targets else None
        record: dict[str, object] = {
            "fold": fold,
            "training": [sample.sample_id for sample in training],
            "fitting": [sample.sample_id for sample in fitting],
            "held_out": [sample.sample_id for sample in heldout],
        }
        try:
            model = fit(tuple(measured[sample.sample_id] for sample in fitting), targets, count_convention=convention)
        except (ValueError, RuntimeError) as error:
            logger.warning("Length training fold unavailable: %s", error)
            record["status"] = "unavailable"
            record["reason"] = "fold_fit_failed"
        else:
            record["status"] = "fitted"
            for sample in heldout:
                if sample.sample_id in measured:
                    result = predict(measured[sample.sample_id], model)
                    predictions[sample.sample_id] = result.estimated_repeat_count
                    if result.reasons:
                        reasons[sample.sample_id] = list(result.reasons)
        fold_records.append(record)
    final_model = None
    fitting = tuple(sample for sample in eligible if sample.sample_id in measured)
    if assignments:
        try:
            model = fit(
                tuple(measured[sample.sample_id] for sample in fitting),
                [sample.length_total for sample in fitting],
                count_convention=convention,
            )
            final_model = encode(model)
        except (ValueError, RuntimeError) as error:
            logger.warning("Final development length model unavailable: %s", error)
    metric = length_metrics(
        [float(sample.length_total) for sample in eligible if sample.length_total is not None],
        [predictions[sample.sample_id] for sample in eligible],
        [baseline[sample.sample_id] for sample in eligible],
        [sample.group_id for sample in eligible],
        seed=seed,
    )
    return {
        "status": "available" if any(value is not None for value in predictions.values()) else "unavailable",
        "count_convention": convention,
        "target": "total diploid repeat count in the declared convention",
        "metrics": metric,
        "folds": fold_records,
        "measurement_reasons": reasons,
        "rows": [
            {
                "sample_id": sample.sample_id,
                "group_id": sample.group_id,
                "truth": sample.length_total,
                "prediction": predictions[sample.sample_id],
                "training_mean_baseline": baseline[sample.sample_id],
                "fold": assignments.get(sample.sample_id),
            }
            for sample in eligible
        ],
        "missing_length_truth": len(samples) - len(eligible),
        "final_model_status": "research-development",
        "approval": None,
    }, final_model


def _validate_arguments(args: object, output: Path) -> tuple[Path, Path | None, str, int, int, str]:
    manifest = getattr(args, "manifest", None)
    reference = getattr(args, "reference", None)
    target = getattr(args, "target", "auto")
    folds = getattr(args, "folds", 5)
    seed = getattr(args, "seed", 20260915)
    convention = getattr(args, "count_convention", "source-reported")
    if not isinstance(manifest, Path) or (reference is not None and not isinstance(reference, Path)):
        raise ValueError("cohort manifest and optional reference must be Paths")
    if target not in {"auto", "length", "callers", "both"}:
        raise ValueError("unsupported cohort target")
    if convention not in {"source-reported", "complete", "canonical-only"}:
        raise ValueError("unsupported cohort counting convention")
    group_folds({}, folds=folds, seed=seed)
    for field in ("caller_runs", "caller_policies"):
        if getattr(args, field, None) is not None and not isinstance(getattr(args, field), Path):
            raise ValueError("cohort caller asset arguments must be Paths")
    if not output.is_dir() or output.is_symlink() or any(output.iterdir()):
        raise ValueError("cohort output must be an empty staged directory")
    return manifest, reference, target, folds, seed, convention


def run_cohort_calibration(args: object, output: Path) -> bool:
    """Evaluate optional length and mutation truth in an atomic CLI staging directory.

    Args:
        args: manifest/reference Paths; target, folds, seed, count_convention;
            optional caller_runs and caller_policies Paths.
        output: Empty private staging directory owned by the existing CLI atomic adapter.

    Returns:
        Whether any requested scientific target produced assessable measurements.
        Unavailable outcomes still write a complete report and return False.

    Raises:
        ValueError: For malformed inputs or unsupported options before dependent I/O.
    """
    manifest, reference, target, folds, seed, convention = _validate_arguments(args, output)
    samples = read_cohort_manifest(manifest)
    targets = [
        name
        for name in ("length", "callers")
        if target in {name, "both"}
        or (
            target == "auto"
            and any(
                sample.length_total is not None if name == "length" else sample.genotype is not None
                for sample in samples
            )
        )
    ]
    if "length" in targets and any(sample.length_total is not None for sample in samples) and reference is None:
        raise ValueError("cohort length extraction requires an explicit reference")
    alignment_hashes = audit_alignment_duplicates(samples)
    document: dict[str, Any] = {
        "schema_version": "development-cohort-comparison-v1",
        "evidence_status": "exploratory-development",
        "targets": targets,
        "folds_requested": folds,
        "seed": seed,
        "count_convention": convention,
        "count_conversion": None,
        "sample_count": len(samples),
        "approval": None,
        "alignment_content_sha256": alignment_hashes,
        "identity_audit": "declared biological groups; identical file bytes refused; distinct bytes do not prove independence",
        "length_truth_unavailable": {
            sample.sample_id: sample.length_reason or "missing_allele_pair"
            for sample in samples
            if sample.length_total is None
        },
        "manifest_sha256": hashlib.sha256(read_regular_path(manifest)).hexdigest(),
        "limitations": "Previously examined data are exploratory. No independent validation or deployment approval.",
    }
    if "length" in targets:
        if not any(sample.length_total is not None for sample in samples):
            document["length"] = {"status": "unavailable", "reason": "no_length_truth"}
        else:
            assert reference is not None
            result, model = _length_comparison(samples, reference, folds=folds, seed=seed, convention=convention)
            document["length"] = result
            if model is not None:
                write_json(output / "length-model.json", model)
    if "callers" in targets:
        document["callers"] = compare_caller_arms(
            samples, getattr(args, "caller_policies", None), getattr(args, "caller_runs", None), folds=folds, seed=seed
        )
    document["comparison_sha256"] = canonical_sha256(document)
    write_json(output / "metrics.json", document)
    environment = Environment(loader=PackageLoader("vntyper", "templates"), autoescape=select_autoescape(["html"]))
    (output / "report.html").write_text(
        environment.get_template("calibration_cohort_report.html").render(result=document), encoding="utf-8"
    )
    write_checksums(output)
    return any(document[name]["status"] == "available" for name in targets)
