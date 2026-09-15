"""Training-only affine length fitting and frozen baseline contracts."""

from dataclasses import FrozenInstanceError, replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_length_protocol import protocol_document
from tests.unit.test_length_estimation import features
from tests.unit.test_length_model import model_document

pytestmark = pytest.mark.unit


def roster_for(rows):
    metrics = import_module("vntyper.scripts.calibration_length_metrics")
    return metrics.decode_length_eligible_roster(
        [
            {"key": row.key, "group_key": row.group_key, "strata": ["training"]}
            for row in sorted(rows, key=lambda item: item.group_key)
        ]
    )


def protocol_for(*kinds, qc=None):
    protocol = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    raw["candidate_grid"] = sorted(
        ({"candidate_id": kind.lower(), "model_kind": kind} for kind in kinds), key=lambda item: item["candidate_id"]
    )
    raw["maximum_candidate_count"] = 4
    if qc is not None:
        raw["qc"].update(qc)
    return protocol.decode_length_protocol(raw)


def training_row(key, x, truth, *, feature="A", role="training", domain="synthetic"):
    training = import_module("vntyper.scripts.calibration_length")
    if feature == "A":
        measured = features(depths=(100, 100 * x, 100 * x, 100, 100, 100), manifest_key=key)
    else:
        measured = features(depths=(100, 100 * x, 100 * x, 100 * x, 100 * x, 100), manifest_key=key)
    return training.LengthTrainingRow(
        key=key,
        group_key=f"group-{key}",
        role=role,
        features=measured,
        truth_boundary_definition="complete-core-plus-invariant-units-v1",
        total_truth_repeat_count=truth,
        evidence_domain=domain,
    )


def metadata_document(rows, roster, *, maximum_condition_number=1e12, **changes):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    measured = rows[0].features
    template = model_document()
    applicability = template["applicability"]
    context = measured.provenance.measurement_context
    applicability.update(
        assemblies=[measured.assembly],
        assay_classes=[measured.assay_class],
        input_scopes=[measured.input_scope],
        preprocessing_ids=[context.preprocessing_id],
        aligner_name=context.aligner.name,
        aligner_version=context.aligner.version,
        aligner_arguments_sha256=context.aligner.arguments_sha256,
        primary_secondary_marking=context.aligner.primary_secondary_marking,
        counting_policy_sha256=measured.counting_policy_sha256,
    )
    producer = candidates.decode_candidate_producer(template["producer"])
    raw = {
        "schema_version": "length-training-metadata-v1",
        "study_sha256": "d" * 64,
        "training_evidence_sha256": "e" * 64,
        "training_roster_sha256": roster.sha256,
        "annotation_sha256": measured.annotation_sha256,
        "counting_policy_sha256": measured.counting_policy_sha256,
        "applicability": applicability,
        "qc": {
            "minimum_denominator_mean_depth": 10,
            "minimum_denominator_covered_fraction": 0.9,
            "minimum_denominator_supporting_fragments": 100,
        },
        "producer": candidates.candidate_producer_document(producer),
        "maximum_condition_number": maximum_condition_number,
    }
    raw.update(changes)
    return raw


def metadata_for(rows, roster, **changes):
    training = import_module("vntyper.scripts.calibration_length")
    return training.decode_length_training_metadata(metadata_document(rows, roster, **changes), roster)


def fit(rows, *kinds, metadata_changes=None, qc=None):
    training = import_module("vntyper.scripts.calibration_length")
    roster = roster_for(rows)
    metadata = metadata_for(rows, roster, **(metadata_changes or {}))
    return training.fit_length_hypotheses(tuple(rows), roster, protocol_for(*kinds, qc=qc), metadata)


def test_exact_affine_a_fit_materializes_strict_model_bounds_and_frozen_baseline():
    model_contract = import_module("vntyper.scripts.length_model")
    rows = tuple(training_row(str(index), x, 10 + 50 * x) for index, x in enumerate((1, 2, 3), start=1))
    result = fit(rows, "affine-A", "physical-A")
    assert [outcome.candidate_id for outcome in result.outcomes] == ["affine-a", "physical-a"]
    affine = result.outcomes[0]
    assert affine.status == "fitted"
    assert affine.reasons == ()
    assert affine.model.intercept == pytest.approx(10)
    assert affine.model.coefficients == pytest.approx((50,))
    assert affine.model.feature_bounds["A"].minimum == pytest.approx(0.8)
    assert affine.model.feature_bounds["A"].maximum == pytest.approx(3.2)
    assert model_contract.decode_length_model(model_contract.encode_length_model(affine.model)) == affine.model
    assert result.outcomes[1].status == "ineligible"
    assert result.outcomes[1].reasons == ("physical_model_geometry_evidence_unsupported",)
    assert result.baseline.mean_total_repeat_count == 110
    assert result.baseline.independent_group_count == 3
    assert result.baseline.training_roster_sha256 == roster_for(rows).sha256
    assert not hasattr(result.baseline, "member_keys")
    assert not hasattr(result, "selected_candidate")
    assert isinstance(result.outcomes, tuple)
    assert not hasattr(rows[0], "allele_1_truth")
    with pytest.raises(FrozenInstanceError):
        result.baseline.mean_total_repeat_count = 0


def test_training_labels_change_fit_but_nontraining_rows_are_refused_at_source_boundary():
    training = import_module("vntyper.scripts.calibration_length")
    base = tuple(training_row(str(index), x, 10 + 50 * x) for index, x in enumerate((1, 2, 3), start=1))
    changed = (*base[:-1], replace(base[-1], total_truth_repeat_count=200))
    first = fit(base, "affine-A").outcomes[0].model
    second = fit(changed, "affine-A").outcomes[0].model
    assert first.coefficients != second.coefficients
    heldout = (replace(base[0], role="validation"), *base[1:])
    roster = roster_for(heldout)
    with pytest.raises(ValueError, match="training role"):
        training.fit_length_hypotheses(heldout, roster, protocol_for("affine-A"), metadata_for(heldout, roster))


def test_exact_affine_f_fit_uses_f_without_joint_feature_expansion():
    rows = tuple(training_row(str(index), x, 10 + 50 * x, feature="F") for index, x in enumerate((1, 2, 3), start=1))
    outcome = fit(rows, "affine-F").outcomes[0]
    assert outcome.status == "fitted"
    assert outcome.model.feature_order == ("F",)
    assert outcome.model.intercept == pytest.approx(10)
    assert outcome.model.coefficients == pytest.approx((50,))


def test_roster_exactly_binds_keys_and_groups_and_duplicate_groups_are_rejected():
    training = import_module("vntyper.scripts.calibration_length")
    rows = tuple(training_row(str(index), x, 10 + 50 * x) for index, x in enumerate((1, 2, 3), start=1))
    roster = roster_for(rows)
    metadata = metadata_for(rows, roster)
    protocol = protocol_for("affine-A")
    with pytest.raises(ValueError, match="roster exactly"):
        training.fit_length_hypotheses(rows[:-1], roster, protocol, metadata)
    swapped = (replace(rows[0], group_key=rows[1].group_key), rows[1], rows[2])
    with pytest.raises(ValueError, match="duplicate.*group"):
        training.fit_length_hypotheses(swapped, roster, protocol, metadata)
    swapped_features = (replace(rows[0], features=rows[1].features), rows[1], rows[2])
    with pytest.raises(ValueError, match="manifest key"):
        training.fit_length_hypotheses(swapped_features, roster, protocol, metadata)
    with pytest.raises(ValueError, match="digest"):
        training.fit_length_hypotheses(rows, roster, protocol, replace(metadata, training_roster_sha256="0" * 64))
    other_roster = roster_for((replace(rows[0], group_key="other-group"), *rows[1:]))
    with pytest.raises(ValueError, match="supplied roster"):
        training.decode_length_training_metadata(metadata_document(rows, roster), other_roster)


def test_any_missing_feature_or_denominator_qc_failure_fails_whole_candidate_without_dropping_row():
    rows = [training_row("1", 1, 60), training_row("2", 2, 110), training_row("3", 3, 160)]
    rows[1] = replace(rows[1], features=features(depths=(100, 200, 200, 0, 0, 100), manifest_key=rows[1].key))
    result = fit(rows, "affine-A", "affine-F")
    outcomes = {outcome.model_kind: outcome for outcome in result.outcomes}
    assert outcomes["affine-A"].status == "ineligible"
    assert outcomes["affine-A"].reasons == ("missing_A",)
    assert outcomes["affine-A"].model is None
    assert outcomes["affine-F"].status == "fitted"
    assert result.baseline.independent_group_count == 3

    qc_rows = tuple(training_row(str(index), x, 10 + 50 * x) for index, x in enumerate((1, 2, 3), start=1))
    qc_result = fit(
        qc_rows,
        "affine-A",
        metadata_changes={
            "qc": {
                "minimum_denominator_mean_depth": 101,
                "minimum_denominator_covered_fraction": 0.9,
                "minimum_denominator_supporting_fragments": 100,
            }
        },
        qc={"minimum_denominator_mean_depth": 101},
    )
    assert qc_result.outcomes[0].reasons == ("low_invariant_mean_depth",)


def test_zero_range_rank_and_condition_failures_are_explicit_candidate_outcomes():
    one = (training_row("1", 1, 60),)
    assert fit(one, "affine-A").outcomes[0].reasons == ("rank_deficient_design",)
    constant = tuple(training_row(str(index), 2, truth) for index, truth in enumerate((100, 110), start=1))
    assert fit(constant, "affine-A").outcomes[0].reasons == ("zero_training_feature_range",)
    shifted = (
        training_row("1", 10_000_000, 100),
        training_row("2", 10_000_001, 110),
        training_row("3", 10_000_002, 120),
    )
    outcome = fit(shifted, "affine-A", metadata_changes={"maximum_condition_number": 1e6}).outcomes[0]
    assert outcome.reasons == ("ill_conditioned_design",)


def test_noisy_finite_data_fit_without_rounding_or_postfit_clipping():
    rows = tuple(training_row(str(index), x, y) for index, (x, y) in enumerate(((1, 61), (2, 108), (3, 161)), start=1))
    outcome = fit(rows, "affine-A").outcomes[0]
    assert outcome.status == "fitted"
    assert outcome.model.intercept == pytest.approx(10)
    assert outcome.model.coefficients[0] == pytest.approx(50)


def test_large_finite_truth_does_not_overflow_the_training_mean():
    rows = (training_row("1", 1, 1e308), training_row("2", 2, 1e308))
    result = fit(rows, "affine-A")
    assert result.baseline.mean_total_repeat_count == 1e308
    assert result.outcomes[0].status in {"fitted", "ineligible"}


@pytest.mark.parametrize(
    "change,match",
    [
        ({"schema_version": "length-training-metadata-v0"}, "schema_version"),
        ({"study_sha256": "bad"}, "SHA256"),
        ({"maximum_condition_number": True}, "condition"),
        ({"maximum_condition_number": float("inf")}, "condition"),
    ],
)
def test_training_metadata_is_closed_hash_bound_and_strict(change, match):
    training = import_module("vntyper.scripts.calibration_length")
    rows = (training_row("1", 1, 60), training_row("2", 2, 110))
    roster = roster_for(rows)
    raw = metadata_document(rows, roster, **change)
    with pytest.raises(ValueError, match=match):
        training.decode_length_training_metadata(raw, roster)
    valid = metadata_for(rows, roster)
    with pytest.raises(ValueError, match="digest"):
        training.length_training_metadata_document(replace(valid, sha256="0" * 64))


def test_training_metadata_closes_root_and_qc_and_hashes_condition_limit():
    training = import_module("vntyper.scripts.calibration_length")
    rows = (training_row("1", 1, 60), training_row("2", 2, 110))
    roster = roster_for(rows)
    raw = metadata_document(rows, roster)
    raw["extra"] = True
    with pytest.raises(ValueError, match="fields"):
        training.decode_length_training_metadata(raw, roster)
    raw = metadata_document(rows, roster)
    raw["qc"]["extra"] = True
    with pytest.raises(ValueError, match="fields"):
        training.decode_length_training_metadata(raw, roster)
    first = metadata_for(rows, roster, maximum_condition_number=1e6)
    second = metadata_for(rows, roster, maximum_condition_number=1e7)
    assert first.sha256 != second.sha256


def test_training_metadata_and_protocol_qc_must_match_exactly():
    training = import_module("vntyper.scripts.calibration_length")
    rows = (training_row("1", 1, 60), training_row("2", 2, 110))
    roster = roster_for(rows)
    metadata = metadata_for(rows, roster)
    with pytest.raises(ValueError, match="QC"):
        training.fit_length_hypotheses(
            rows,
            roster,
            protocol_for("affine-A", qc={"minimum_denominator_mean_depth": 11}),
            metadata,
        )


def test_baseline_encoding_is_deterministic_and_bound_to_training_metadata_and_roster():
    training = import_module("vntyper.scripts.calibration_length")
    rows = (training_row("1", 1, 60), training_row("2", 2, 110), training_row("3", 3, 160))
    result = fit(rows, "affine-A")
    raw = training.length_baseline_document(result.baseline)
    assert raw["baseline_kind"] == "training-mean-v1"
    assert raw["mean_total_repeat_count"] == 110
    assert training.decode_length_baseline(raw) == result.baseline
    with pytest.raises(ValueError, match="digest"):
        training.length_baseline_document(replace(result.baseline, sha256="0" * 64))
    changed = dict(raw)
    changed["training_evidence_sha256"] = "f" * 64
    assert training.decode_length_baseline(changed).sha256 != result.baseline.sha256
