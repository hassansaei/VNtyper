from __future__ import annotations

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.scripts.length_features import DepthPosition
from vntyper.scripts.length_standard_features import (
    STANDARD_ANNOTATION_SHA256,
    STANDARD_FEATURE_DEFINITION_SHA256,
    STANDARD_FEATURE_ORDER,
    STANDARD_REFERENCE_LOCUS_SHA256,
    StandardReadSummary,
    extract_standard_length_features,
)
from vntyper.scripts.length_standard_model import (
    decode_standard_length_model,
    encode_standard_length_model,
    load_standard_length_model,
    predict_standard_length,
)

pytestmark = pytest.mark.unit


def _depth():
    return tuple(
        DepthPosition(
            "chr1",
            position,
            20 if 155188726 <= position < 155191939 else 10,
            (f"fragment-{position}",),
        )
        for position in range(155188296, 155192429)
    )


def _reads() -> StandardReadSummary:
    return StandardReadSummary(4, 1, 120, 2, 150, 300)


def _model_document() -> dict[str, object]:
    return {
        "schema_version": "standard-length-linear-model-v1",
        "model_version": "synthetic-v1",
        "target": {
            "name": "source-reported-diploid-repeat-count",
            "count_convention": "source-reported",
            "unit": "repeat_units",
        },
        "model_kind": "linear-standard13",
        "feature_order": list(STANDARD_FEATURE_ORDER),
        "intercept": 10.0,
        "coefficients": [50.0] + [0.0] * 12,
        "feature_bounds": {name: {"minimum": -100.0, "maximum": 100.0} for name in STANDARD_FEATURE_ORDER},
        "qc": {
            "minimum_denominator_mean_depth": 1.0,
            "minimum_denominator_covered_fraction": 0.5,
            "minimum_denominator_supporting_fragments": 1,
            "minimum_eligible_reads": 1,
        },
        "assembly": "GRCh38",
        "accepted_contigs": ["chr1", "1"],
        "annotation_sha256": STANDARD_ANNOTATION_SHA256,
        "reference_locus_sha256": STANDARD_REFERENCE_LOCUS_SHA256,
        "feature_definition_sha256": STANDARD_FEATURE_DEFINITION_SHA256,
        "model_source": "local-research",
    }


def _measurement():
    return extract_standard_length_features(
        _depth(),
        _reads(),
        assembly="GRCh38",
        contig="chr1",
        reference_locus_sha256=STANDARD_REFERENCE_LOCUS_SHA256,
    )


def test_standard_model_predicts_raw_linear_value_without_rounding() -> None:
    model = decode_standard_length_model(_model_document())
    result = predict_standard_length(_measurement(), model)
    assert result.status == "estimated"
    assert result.estimated_repeat_count == 110.0
    assert result.reasons == ()
    assert result.warnings == ()
    assert decode_standard_length_model(encode_standard_length_model(model)) == model


def test_feature_bounds_are_warnings_not_unavailability() -> None:
    document = _model_document()
    document["feature_bounds"]["A"] = {"minimum": 0, "maximum": 1}  # type: ignore[index]
    result = predict_standard_length(_measurement(), decode_standard_length_model(document))
    assert result.status == "estimated"
    assert result.estimated_repeat_count == 110
    assert result.warnings == ("feature_A_outside_training_range",)


def test_standard_prediction_checks_each_denominator_side_and_read_qc() -> None:
    document = _model_document()
    document["qc"] = {
        "minimum_denominator_mean_depth": 11,
        "minimum_denominator_covered_fraction": 0.9,
        "minimum_denominator_supporting_fragments": 200,
        "minimum_eligible_reads": 5,
    }
    result = predict_standard_length(_measurement(), decode_standard_length_model(document))
    assert result.status == "unavailable"
    assert result.estimated_repeat_count is None
    assert "low_left_flank_mean_depth" in result.reasons
    assert "low_left_flank_support" in result.reasons
    assert "low_right_flank_support" in result.reasons
    assert "low_eligible_read_count" in result.reasons


def test_standard_model_is_closed_and_rejects_nonfinite_or_wrong_order() -> None:
    for mutate in (
        lambda row: row.update(extra=True),
        lambda row: row.update(intercept=float("nan")),
        lambda row: row.update(intercept=True),
        lambda row: row.update(intercept=10**400),
        lambda row: row.update(feature_order=list(reversed(STANDARD_FEATURE_ORDER))),
        lambda row: row["target"].update(count_convention="unconfirmed"),  # type: ignore[union-attr]
    ):
        document = _model_document()
        mutate(document)
        with pytest.raises(ValueError):
            decode_standard_length_model(document)


def test_load_standard_model_uses_strict_json_and_packaged_default(tmp_path) -> None:
    model = load_standard_length_model()
    assert model.model_source == "packaged-research"
    path = tmp_path / "model.json"
    document = _model_document()
    document["model_source"] = "local-research"
    path.write_bytes(canonical_json_bytes(document))
    assert load_standard_length_model(path) == decode_standard_length_model(document)
    path.write_text('{"schema_version":"standard-length-linear-model-v1","schema_version":"x"}', encoding="utf-8")
    with pytest.raises(ValueError):
        load_standard_length_model(path)


def test_packaged_standard_model_companion_matches_raw_model() -> None:
    model = load_standard_length_model()
    assert model.model_version == "standard13-bayesian-v1"
    assert model.count_convention == "source-reported"
