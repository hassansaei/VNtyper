from __future__ import annotations

import pytest

from vntyper.scripts.calibration_standard_length import (
    collapse_scaled_linear_model,
    fit_standard_length_model,
    standard_length_training_reasons,
)
from vntyper.scripts.length_features import DepthPosition
from vntyper.scripts.length_standard_features import (
    STANDARD_FEATURE_ORDER,
    STANDARD_REFERENCE_LOCUS_SHA256,
    StandardReadSummary,
    extract_standard_length_features,
)

pytestmark = pytest.mark.unit


def _measurement(
    core_depth: int,
    *,
    denominator_depth: int = 10,
    fragments_available: bool = True,
    reads: StandardReadSummary | None = None,
):
    depth = tuple(
        DepthPosition(
            "chr1",
            position,
            core_depth if 155188726 <= position < 155191939 else denominator_depth,
            (f"fragment-{position}",) if fragments_available and denominator_depth > 0 else None,
        )
        for position in range(155188296, 155192429)
    )
    return extract_standard_length_features(
        depth,
        reads or StandardReadSummary(4, 1, 120, 2, 150, 300),
        assembly="GRCh38",
        contig="chr1",
        reference_locus_sha256=STANDARD_REFERENCE_LOCUS_SHA256,
    )


def test_collapse_scaled_linear_model_preserves_predictions() -> None:
    intercept, coefficients = collapse_scaled_linear_model(
        scaled_intercept=5,
        scaled_coefficients=(2, -3),
        means=(10, 20),
        scales=(2, 5),
    )
    assert coefficients == (1, -0.6)
    assert intercept == 7
    assert intercept + coefficients[0] * 12 + coefficients[1] * 25 == pytest.approx(4)


def test_fit_standard_model_uses_training_only_bayesian_result(monkeypatch) -> None:
    rows = (_measurement(10), _measurement(20), _measurement(30))
    targets = (60.0, 110.0, 160.0)

    def fake_fit(matrix, targets):
        assert len(matrix) == 3
        assert targets == [60.0, 110.0, 160.0]
        return 10.0, (50.0,) + (0.0,) * 12

    monkeypatch.setattr("vntyper.scripts.calibration_standard_length._fit_raw_bayesian", fake_fit)
    model = fit_standard_length_model(rows, targets, count_convention="source-reported")
    assert model.intercept == 10
    assert model.coefficients[0] == 50
    assert model.count_convention == "source-reported"


def test_fit_allows_identical_vectors_but_rejects_nonintegral_truth(monkeypatch) -> None:
    monkeypatch.setattr(
        "vntyper.scripts.calibration_standard_length._fit_raw_bayesian",
        lambda matrix, targets: (1.0, (1.0,) * len(STANDARD_FEATURE_ORDER)),
    )
    row = _measurement(20)
    assert (
        fit_standard_length_model((row, row), (10, 20), count_convention="source-reported").model_kind
        == "linear-standard13"
    )
    with pytest.raises(ValueError, match="integral"):
        fit_standard_length_model((row, _measurement(30)), (10.5, 20), count_convention="source-reported")
    with pytest.raises(ValueError, match="integral"):
        fit_standard_length_model((row, _measurement(30)), (10**400, 20), count_convention="source-reported")


def test_public_training_reasons_applies_frozen_qc() -> None:
    row = _measurement(20)
    assert standard_length_training_reasons(row) == ()

    low_evidence = _measurement(
        20,
        denominator_depth=0,
        fragments_available=False,
        reads=StandardReadSummary(0, 0, 0, 0, 0, 0),
    )
    reasons = standard_length_training_reasons(low_evidence)
    for region in ("invariant", "left_flank", "right_flank"):
        assert f"low_{region}_mean_depth" in reasons
        assert f"low_{region}_covered_fraction" in reasons
        assert f"missing_{region}_support" in reasons
    assert "low_eligible_read_count" in reasons


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"scaled_intercept": 0, "scaled_coefficients": (), "means": (), "scales": ()}, "same non-zero"),
        (
            {"scaled_intercept": 0, "scaled_coefficients": (1,), "means": (0,), "scales": (True,)},
            "numeric",
        ),
        (
            {"scaled_intercept": 0, "scaled_coefficients": (1,), "means": (0,), "scales": (0,)},
            "positive scales",
        ),
        (
            {"scaled_intercept": 1e308, "scaled_coefficients": (1e308,), "means": (-1e308,), "scales": (1,)},
            "collapsed",
        ),
    ],
)
def test_collapse_scaled_linear_model_rejects_invalid_parameters(arguments, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        collapse_scaled_linear_model(**arguments)


def test_fit_standard_model_rejects_malformed_inputs_and_fitted_shape(monkeypatch) -> None:
    row = _measurement(20)
    with pytest.raises(ValueError, match="convention"):
        fit_standard_length_model((row, row), (10, 20), count_convention="unknown")
    with pytest.raises(ValueError, match="at least two"):
        fit_standard_length_model((row,), (10,), count_convention="source-reported")
    with pytest.raises(ValueError, match="measurement count"):
        fit_standard_length_model((row, row), (10,), count_convention="source-reported")
    with pytest.raises(ValueError, match="typed measurements"):
        fit_standard_length_model((row, object()), (10, 20), count_convention="source-reported")  # type: ignore[arg-type]

    low_evidence = _measurement(20, denominator_depth=0, fragments_available=False)
    with pytest.raises(ValueError, match="frozen QC"):
        fit_standard_length_model((row, low_evidence), (10, 20), count_convention="source-reported")

    monkeypatch.setattr(
        "vntyper.scripts.calibration_standard_length._fit_raw_bayesian",
        lambda matrix, targets: (1.0, (1.0,)),
    )
    with pytest.raises(ValueError, match="coefficient count"):
        fit_standard_length_model((row, _measurement(30)), (10, 20), count_convention="source-reported")
