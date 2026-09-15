"""Captured Kestrel rows replay through the production decision functions."""

from dataclasses import replace

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
from vntyper.scripts.calibration_kestrel_capture import KESTREL_RAW_COLUMNS, build_kestrel_capture
from vntyper.scripts.calibration_kestrel_replay import (
    decode_kestrel_replay_result,
    kestrel_replay_document,
    kestrel_replay_prefilter_frame,
    kestrel_replay_selected_frame,
    replay_kestrel_capture,
)
from vntyper.scripts.identity_candidates import translation_component_from_config
from vntyper.scripts.kestrel_genotyping import _resolve_selection, process_kmer_results
from vntyper.scripts.nomenclature import nomenclature_config

pytestmark = pytest.mark.unit

_GG = "/components/kestrel/alt_filtering/gg_depth_score_threshold"
_FLOOR = "/components/kestrel/confidence_assignment/reporting_floor"
_ACTIVE = "/components/kestrel/confidence_assignment/var_active_region_threshold"
_LOW = "/components/kestrel/confidence_assignment/depth_score_thresholds/low"
_HIGH = "/components/kestrel/confidence_assignment/depth_score_thresholds/high"
_ALT_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/low"
_ALT_MID_LOW = "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low"
_ALT_MID_HIGH = "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high"


def _policy(config: dict[str, object]):
    confidence = config["confidence_assignment"]
    alt_filter = config["alt_filtering"]
    assert isinstance(confidence, dict) and isinstance(alt_filter, dict)
    depth = confidence["depth_score_thresholds"]
    alt = confidence["alt_depth_thresholds"]
    assert isinstance(depth, dict) and isinstance(alt, dict)
    return decode_caller_policy_values(
        {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": ["kestrel"],
            "values": {
                _GG: alt_filter["gg_depth_score_threshold"],
                _FLOOR: confidence["reporting_floor"],
                _ACTIVE: confidence["var_active_region_threshold"],
                _LOW: depth["low"],
                _HIGH: depth["high"],
                _ALT_LOW: alt["low"],
                _ALT_MID_LOW: alt["mid_low"],
                _ALT_MID_HIGH: alt["mid_high"],
            },
        }
    )


def _identity_inputs() -> tuple[str, str, pd.DataFrame]:
    motifs = nomenclature_config["motifs"]
    assert isinstance(motifs, dict)
    pair = motifs["C"] + motifs["S"]
    return "S-C", pair, pd.DataFrame({"Motif": ["S"], "Motif_sequence": [motifs["S"]]})


def _capture(frame: pd.DataFrame, config: dict[str, object], policy=None):
    _, _, motifs = _identity_inputs()
    return build_kestrel_capture(
        frame,
        motifs,
        kestrel_config=config,
        baseline_policy=_policy(config) if policy is None else policy,
        selection=_resolve_selection(config),
        identity_component=translation_component_from_config(nomenclature_config),
        decision_profile_sha256="1" * 64,
        reference_file_bytes=b"invented-reference",
        motif_reference_file_bytes=b"invented-motif-reference",
        kestrel_jar_bytes=b"invented-jar",
        capture_policy={"schema_version": "synthetic-recruitment-v1", "kmer_sizes": [20]},
    )


def _raw(*, depth_alt: int = 7, depth_region: int = 500, ref: str = "G", alt: str = "GG") -> pd.DataFrame:
    pair_name, pair_sequence, _ = _identity_inputs()
    return kestrel_stage_frame(
        "raw",
        motifs=pair_name,
        motif_sequence=pair_sequence,
        depth_alt=depth_alt,
        depth_region=depth_region,
        ref=ref,
        alt=alt,
    )


def _candidate(policy, **changes):
    document = {
        "schema_version": "calibration-caller-policy-values-v1",
        "required_callers": list(policy.required_callers),
        "values": dict(policy.values),
    }
    values = document["values"]
    assert isinstance(values, dict)
    values.update(changes)
    return decode_caller_policy_values(document)


def _replay_rows(config: dict[str, object], depths: list[tuple[int, int]]) -> pd.DataFrame:
    rows = []
    for ordinal, (alternate, active) in enumerate(depths):
        row = _raw(depth_alt=alternate, depth_region=active)
        row.loc[0, "POS"] = 67 + ordinal
        rows.append(row)
    frame = pd.concat(rows, ignore_index=True)
    capture = _capture(frame, config)
    result = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    return kestrel_replay_prefilter_frame(result)


def test_lower_declared_gg_threshold_rescues_a_complete_below_baseline_candidate() -> None:
    config = kestrel_config(
        **{
            "confidence_assignment.reporting_floor": 0.003,
            "confidence_assignment.depth_score_thresholds.low": 0.003,
        }
    )
    capture = _capture(_raw(depth_alt=4, depth_region=1000), config)
    baseline = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    lowered = _candidate(capture.baseline_policy, **{_GG: 0.003})
    rescued = replay_kestrel_capture(
        capture,
        lowered,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )

    assert baseline.disposition == "no-call"
    assert bool(kestrel_replay_prefilter_frame(baseline).iloc[0]["alt_filter_pass"]) is False
    assert rescued.disposition == "called"
    assert bool(kestrel_replay_prefilter_frame(rescued).iloc[0]["alt_filter_pass"]) is True


def test_replay_uses_inclusive_gg_and_reporting_floor_boundaries() -> None:
    gg_config = kestrel_config(
        **{
            "alt_filtering.gg_depth_score_threshold": 0.004,
            "confidence_assignment.reporting_floor": 0.0,
            "confidence_assignment.depth_score_thresholds.low": 0.0,
            "confidence_assignment.depth_score_thresholds.high": 1.0,
        }
    )
    gg = _replay_rows(gg_config, [(3, 1000), (4, 1000), (5, 1000)])
    assert gg["alt_filter_pass"].tolist() == [False, True, True]

    floor_config = kestrel_config(
        **{
            "alt_filtering.gg_depth_score_threshold": 0.0,
            "confidence_assignment.reporting_floor": 0.004,
            "confidence_assignment.depth_score_thresholds.low": 0.004,
            "confidence_assignment.depth_score_thresholds.high": 1.0,
        }
    )
    floor = _replay_rows(floor_config, [(3, 1000), (4, 1000), (5, 1000)])
    assert floor["depth_confidence_pass"].tolist() == [False, True, True]


def test_replay_uses_closed_low_high_depth_score_band_and_ordered_high_boundary() -> None:
    low_config = kestrel_config(
        **{
            "alt_filtering.gg_depth_score_threshold": 0.0,
            "confidence_assignment.reporting_floor": 0.0,
            "confidence_assignment.depth_score_thresholds.low": 0.004,
            "confidence_assignment.depth_score_thresholds.high": 1.0,
        }
    )
    low = _replay_rows(low_config, [(3, 1000), (4, 1000), (5, 1000)])
    assert low["Confidence"].tolist() == ["Negative", "Low_Precision", "Low_Precision"]

    high_config = kestrel_config(
        **{
            "alt_filtering.gg_depth_score_threshold": 0.0,
            "confidence_assignment.reporting_floor": 0.0,
            "confidence_assignment.depth_score_thresholds.low": 0.0,
            "confidence_assignment.depth_score_thresholds.high": 0.004,
            "confidence_assignment.alt_depth_thresholds.low": 0,
            "confidence_assignment.alt_depth_thresholds.mid_low": 1,
            "confidence_assignment.alt_depth_thresholds.mid_high": 2,
        }
    )
    high = _replay_rows(high_config, [(3, 1000), (4, 1000), (5, 1000)])
    assert high["Confidence"].tolist() == ["Low_Precision", "Low_Precision", "High_Precision*"]


def test_replay_uses_discrete_alt_partition_boundaries() -> None:
    config = kestrel_config(
        **{
            "alt_filtering.gg_depth_score_threshold": 0.0,
            "confidence_assignment.reporting_floor": 0.0,
            "confidence_assignment.depth_score_thresholds.low": 0.0,
            "confidence_assignment.depth_score_thresholds.high": 0.001,
        }
    )
    low_and_mid_low = _replay_rows(config, [(19, 1900), (20, 2000), (21, 2100)])
    assert low_and_mid_low["Confidence"].tolist() == ["Low_Precision", "Low_Precision", "High_Precision"]
    mid_high = _replay_rows(config, [(99, 9900), (100, 10000), (101, 10100)])
    assert mid_high["Confidence"].tolist() == ["High_Precision", "High_Precision*", "High_Precision*"]


def test_active_region_threshold_is_inert_under_the_mandated_discrete_partition() -> None:
    config = kestrel_config(
        **{
            "alt_filtering.gg_depth_score_threshold": 0.0,
            "confidence_assignment.reporting_floor": 0.0,
            "confidence_assignment.depth_score_thresholds.low": 0.0,
            "confidence_assignment.depth_score_thresholds.high": 0.001,
            "confidence_assignment.var_active_region_threshold": 200,
        }
    )
    # The only adjacent integer depths are <= low or >= mid_low. The production
    # region predicate requires low < alt < mid_low, so its open interval is empty.
    frame = _replay_rows(config, [(20, 199), (20, 200), (20, 201)])
    assert frame["Confidence"].tolist() == ["Low_Precision", "Low_Precision", "Low_Precision"]


def test_unchanged_policy_replay_exactly_matches_production_including_suppressed_and_invalid_rows(tmp_path) -> None:
    frame = pd.concat(
        [
            _raw(depth_alt=7, depth_region=500),
            _raw(depth_alt=120, depth_region=12000, ref="C", alt="CGGCA"),
            _raw(depth_alt=120, depth_region=12000, ref="G", alt="GGG"),
        ],
        ignore_index=True,
    )[list(KESTREL_RAW_COLUMNS)]
    config = kestrel_config()
    capture = _capture(frame, config)
    replay = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    _, _, motifs = _identity_inputs()
    production = process_kmer_results(
        frame,
        motifs,
        str(tmp_path),
        config,
        identity_component=capture.identity_component,
    )
    prefilter = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t", keep_default_na=False, dtype=str)
    replay_prefilter = kestrel_replay_prefilter_frame(replay).drop(columns=["__Calibration_Source_Row_Ordinal"])
    replay_prefilter = replay_prefilter.fillna("").astype(str)

    pd.testing.assert_frame_equal(replay_prefilter.reset_index(drop=True), prefilter)
    replay_selected = kestrel_replay_selected_frame(replay).drop(columns=["__Calibration_Source_Row_Ordinal"])
    pd.testing.assert_frame_equal(replay_selected, production, check_dtype=False)
    assert len(replay.source_dispositions) == 3
    assert replay.source_dispositions[1].blocking_gates == ("flag_filter_pass",)
    assert "is_valid_frameshift" in replay.source_dispositions[2].blocking_gates


def test_zero_active_depth_is_retained_as_an_explained_no_call() -> None:
    config = kestrel_config()
    capture = _capture(_raw(depth_alt=4, depth_region=0), config)
    result = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    prefilter = kestrel_replay_prefilter_frame(result)
    assert result.disposition == "no-call"
    assert prefilter.iloc[0]["Depth_Score"] is None or pd.isna(prefilter.iloc[0]["Depth_Score"])
    assert result.source_dispositions[0].blocking_gates == ("depth_confidence_pass", "alt_filter_pass")


def test_flags_are_applied_before_selection_so_a_weaker_unflagged_row_wins() -> None:
    artifact = _raw(depth_alt=120, depth_region=12000, ref="C", alt="CGGCA")
    valid = _raw(depth_alt=7, depth_region=500)
    valid.loc[0, "POS"] = 68
    frame = pd.concat([artifact, valid], ignore_index=True)
    capture = _capture(frame, kestrel_config())
    result = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    selected = kestrel_replay_selected_frame(result)
    assert selected.iloc[0]["POS"] == 68
    assert result.source_dispositions[0].blocking_gates == ("flag_filter_pass",)
    assert result.source_dispositions[1].selected is True


def test_equal_evidence_tie_uses_the_frozen_production_position_order() -> None:
    later = _raw(depth_alt=7, depth_region=500)
    later.loc[0, "POS"] = 68
    earlier = _raw(depth_alt=7, depth_region=500)
    earlier.loc[0, "POS"] = 67
    capture = _capture(pd.concat([later, earlier], ignore_index=True), kestrel_config())
    result = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    selected = kestrel_replay_selected_frame(result)
    assert selected.iloc[0]["POS"] == 67
    assert result.source_dispositions[1].selected is True


def test_empty_complete_capture_is_unassessable_and_has_no_synthetic_rows() -> None:
    config = kestrel_config()
    capture = _capture(_raw().iloc[0:0], config)
    result = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    assert result.disposition == "unassessable-no-candidates"
    assert result.source_dispositions == ()
    assert kestrel_replay_prefilter_frame(result).empty


def test_changed_recruitment_commitment_requires_recapture_before_evaluation() -> None:
    capture = _capture(_raw(), kestrel_config())
    with pytest.raises(ValueError, match="recapture"):
        replay_kestrel_capture(capture, capture.baseline_policy, capture_policy_sha256="0" * 64)


def test_replay_result_roundtrip_and_public_boundaries_reject_forged_digests() -> None:
    capture = _capture(_raw(), kestrel_config())
    result = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    assert decode_kestrel_replay_result(kestrel_replay_document(result)) == result
    forged = replace(result, sha256="0" * 64)
    with pytest.raises(ValueError, match="canonical|digest"):
        kestrel_replay_document(forged)
    with pytest.raises(ValueError, match="canonical|digest"):
        kestrel_replay_selected_frame(forged)


@pytest.mark.parametrize(
    "mutation",
    [
        "wrong_disposition",
        "missing_prefilter_row",
        "missing_selected_row",
        "wrong_selected_marker",
        "wrong_blocking_gates",
    ],
)
def test_replay_decoder_rejects_internally_inconsistent_evidence(mutation: str) -> None:
    capture = _capture(_raw(), kestrel_config())
    result = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    document = kestrel_replay_document(result)
    prefilter = document["prefilter"]
    selected = document["selected"]
    dispositions = document["source_dispositions"]
    assert isinstance(prefilter, dict) and isinstance(selected, dict) and isinstance(dispositions, list)
    assert isinstance(prefilter["rows"], list) and isinstance(selected["rows"], list)
    assert isinstance(dispositions[0], dict)
    if mutation == "wrong_disposition":
        document["disposition"] = "no-call"
    elif mutation == "missing_prefilter_row":
        prefilter["rows"] = []
    elif mutation == "missing_selected_row":
        selected["rows"] = []
    elif mutation == "wrong_selected_marker":
        dispositions[0]["selected"] = False
    else:
        dispositions[0]["blocking_gates"] = ["is_frameshift"]

    with pytest.raises(ValueError, match="prefilter|selected|disposition"):
        decode_kestrel_replay_result(document)
