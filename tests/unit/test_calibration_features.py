"""Separated calibration feature and truth-label artifacts."""

import pytest

from vntyper.scripts.calibration_features import (
    decode_feature_artifact,
    decode_label_artifact,
    validate_artifact_alignment,
    validate_training_recipe,
)
from vntyper.scripts.calibration_manifest import decode_partition_manifest

pytestmark = pytest.mark.unit


def _features(**changes: object) -> dict[str, object]:
    values: dict[str, object] = {
        "assay_class": "capture-short-read",
        "xd_availability_count": 2,
        "xd_availability_fraction": 0.5,
        "xd_median": 5.0,
        "xd_interquartile_range": 1.0,
    }
    values.update(changes)
    return {
        "schema_version": "calibration-features-v1",
        "rows": [{"feature_key": "f-1", "manifest_key": "sample-1", "features": values}],
    }


def _labels() -> dict[str, object]:
    return {
        "schema_version": "calibration-labels-v1",
        "rows": [
            {
                "label_key": "l-1",
                "manifest_key": "sample-1",
                "truth_status": "mutated",
                "expected_identity": "identity-json",
                "expected_display_name": "59dupC",
                "mutation_class": "duplication",
            }
        ],
    }


def test_assay_class_and_four_named_xd_summaries_are_allowlisted() -> None:
    artifact = decode_feature_artifact(_features())

    assert artifact.rows[0].features["assay_class"] == "capture-short-read"
    assert set(artifact.rows[0].features) == {
        "assay_class",
        "xd_availability_count",
        "xd_availability_fraction",
        "xd_median",
        "xd_interquartile_range",
    }


@pytest.mark.parametrize(
    "forbidden",
    [
        "truth",
        "path",
        "sample",
        "batch_identity",
        "raw_tuple",
        "selected_name",
        "post_decision_tier",
        "raw_xd_weight",
        "xd_sum",
        "simulator_allele_length",
        "reference_interval_length",
    ],
)
def test_labels_identifiers_raw_xd_and_unlisted_aggregations_cannot_be_features(forbidden: str) -> None:
    with pytest.raises(ValueError, match="feature"):
        decode_feature_artifact(_features(**{forbidden: "leak"}))


def test_labels_are_separate_strictly_keyed_and_never_accepted_as_features() -> None:
    labels = decode_label_artifact(_labels())

    assert labels.rows[0].expected_display_name == "59dupC"
    with pytest.raises(ValueError):
        decode_feature_artifact(_labels())
    duplicate = _labels()
    rows = duplicate["rows"]
    assert isinstance(rows, list)
    rows.append(dict(rows[0]))
    with pytest.raises(ValueError, match="key"):
        decode_label_artifact(duplicate)


def test_transforms_fit_only_on_training_and_candidates_compare_only_on_policy_selection() -> None:
    recipe = validate_training_recipe(
        {
            "transform_fit_role": "training",
            "candidate_comparison_role": "policy-selection",
            "missing_value_fit_role": "training",
            "cap_fit_role": "training",
            "feature_selection_fit_role": "training",
        }
    )
    assert recipe.candidate_comparison_role == "policy-selection"

    for field in (
        "transform_fit_role",
        "missing_value_fit_role",
        "cap_fit_role",
        "feature_selection_fit_role",
    ):
        raw = {
            "transform_fit_role": "training",
            "candidate_comparison_role": "policy-selection",
            "missing_value_fit_role": "training",
            "cap_fit_role": "training",
            "feature_selection_fit_role": "training",
        }
        raw[field] = "policy-selection"
        with pytest.raises(ValueError, match="training"):
            validate_training_recipe(raw)


def test_feature_and_label_rows_align_by_manifest_but_use_separate_key_namespaces() -> None:
    groups = {
        namespace: [f"{namespace}:sample-1"]
        for namespace in (
            "individual-family",
            "simulated-pair",
            "backbone-seed-lineage",
            "replicate-rerun",
            "depth-series-source",
            "batch",
            "repeat-context",
        )
    }
    manifest = decode_partition_manifest(
        {
            "schema_version": "calibration-partitions-v1",
            "members": [
                {
                    "key": "sample-1",
                    "role": "training",
                    "provenance": "development",
                    "assay_class": "capture-short-read",
                    "groups": groups,
                }
            ],
        }
    )
    features = decode_feature_artifact(_features())
    labels = decode_label_artifact(_labels())

    validate_artifact_alignment(features, labels, manifest)

    shared_key = _labels()
    rows = shared_key["rows"]
    assert isinstance(rows, list) and isinstance(rows[0], dict)
    rows[0]["label_key"] = "f-1"
    with pytest.raises(ValueError, match="separate"):
        validate_artifact_alignment(features, decode_label_artifact(shared_key), manifest)
