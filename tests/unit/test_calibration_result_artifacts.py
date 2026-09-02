"""Pure calibration metrics, evaluation, and attestation encoders."""

from fractions import Fraction

import pytest

from vntyper.scripts.calibration_contract import CandidateMetrics, decode_attestation, decode_metrics
from vntyper.scripts.calibration_objective import CandidateEvaluation
from vntyper.scripts.calibration_result_artifacts import (
    encode_attestation_hashes,
    encode_evaluation,
    encode_metrics,
)

pytestmark = pytest.mark.unit


def _metrics() -> CandidateMetrics:
    return CandidateMetrics(
        "a" * 64,
        0,
        0,
        0,
        Fraction(1),
        Fraction(1),
        1,
        Fraction(0),
        True,
        True,
        "b" * 64,
    )


def test_result_encoders_preserve_exact_metrics_bounds_and_failed_status() -> None:
    metrics = _metrics()
    metrics_document = encode_metrics(metrics)
    evaluation_document = encode_evaluation(
        CandidateEvaluation(metrics, Fraction(3, 4), Fraction(2, 3), (4, 5), Fraction(1, 20))
    )
    attestation_document = encode_attestation_hashes("validation", "a" * 64, "b" * 64, "c" * 64, "d" * 64, False)

    assert decode_metrics(metrics_document).macro_exact_recovery == Fraction(1)
    assert evaluation_document["detection_lower_bound"] == "3/4"
    assert evaluation_document["macro_exact_lower_bound"] == "2/3"
    assert evaluation_document["stratum_counts"] == [4, 5]
    assert evaluation_document["holm_adjusted_p_value"] == "1/20"
    assert decode_attestation(attestation_document).status == "failed"
