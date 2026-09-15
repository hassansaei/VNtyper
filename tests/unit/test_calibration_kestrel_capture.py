"""Kestrel replay captures are complete, strict and hash bound."""

from copy import deepcopy
from dataclasses import FrozenInstanceError, replace

import pandas as pd
import pytest
from Bio.Seq import Seq

from tests.builders import kestrel_config, kestrel_stage_frame
from tests.unit.test_calibration_caller_policy import policy_document
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
from vntyper.scripts.calibration_kestrel_capture import (
    KESTREL_RAW_COLUMNS,
    build_kestrel_capture,
    decode_kestrel_capture,
    kestrel_capture_document,
)
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.identity_candidates import translation_component_from_config
from vntyper.scripts.kestrel_genotyping import _resolve_selection
from vntyper.scripts.nomenclature import nomenclature_config

pytestmark = pytest.mark.unit

_SHA = "1" * 64


def _baseline():
    document = policy_document(include_advntr=False)
    values = document["values"]
    assert isinstance(values, dict)
    values["/components/kestrel/confidence_assignment/reporting_floor"] = 0.00469
    return decode_caller_policy_values(document)


def _motifs() -> pd.DataFrame:
    return pd.DataFrame({"Motif": ["5", "X"], "Motif_sequence": ["SEQ2", "SEQ1"]})


def _capture(frame: pd.DataFrame | None = None):
    return build_kestrel_capture(
        kestrel_stage_frame("raw") if frame is None else frame,
        _motifs(),
        kestrel_config=kestrel_config(),
        baseline_policy=_baseline(),
        selection=_resolve_selection(kestrel_config()),
        identity_component=translation_component_from_config(nomenclature_config),
        decision_profile_sha256=_SHA,
        reference_file_bytes=b">X-5\nSEQ1SEQ2\n",
        motif_reference_file_bytes=b">X\nSEQ1\n>5\nSEQ2\n",
        kestrel_jar_bytes=b"invented-kestrel-jar",
        capture_policy={"schema_version": "synthetic-kestrel-recruitment-v1", "kmer_sizes": [20]},
    )


def test_capture_roundtrip_preserves_ordered_duplicate_raw_rows_and_distinct_commitments() -> None:
    frame = kestrel_stage_frame("raw", rows=2)
    frame.loc[1, "POS"] = frame.loc[0, "POS"]
    capture = _capture(frame)
    document = kestrel_capture_document(capture)
    decoded = decode_kestrel_capture(deepcopy(document))

    assert decoded == capture
    assert [row.source_row_ordinal for row in capture.rows] == [0, 1]
    assert capture.rows[0].position == capture.rows[1].position
    assert document["raw_columns"] == list(KESTREL_RAW_COLUMNS)
    provenance = document["provenance"]
    assert isinstance(provenance, dict)
    assert provenance["reference_file_sha256"] != provenance["parsed_motif_table_sha256"]
    assert provenance["motif_reference_file_sha256"] != provenance["identity_table_sha256"]
    assert capture.sha256 == canonical_sha256(document)
    with pytest.raises(FrozenInstanceError):
        capture.sha256 = "2" * 64


def test_capture_builder_normalizes_canonical_vcf_position_text() -> None:
    frame = kestrel_stage_frame("raw")
    frame["POS"] = frame["POS"].astype(str)

    capture = _capture(frame)

    assert type(capture.rows[0].position) is int
    assert capture.rows[0].position == 67


def test_capture_builder_normalizes_production_biopython_motif_sequence() -> None:
    frame = kestrel_stage_frame("raw")
    frame["Motif_sequence"] = frame["Motif_sequence"].map(Seq)

    capture = _capture(frame)

    assert type(capture.rows[0].motif_sequence) is str
    assert capture.rows[0].motif_sequence == "SEQ1"


@pytest.mark.parametrize("stage", ["scored", "confidence", "flagged", "final", "named"])
def test_final_or_derived_rows_refuse_capture(stage: str) -> None:
    with pytest.raises(ValueError, match="raw columns"):
        _capture(kestrel_stage_frame(stage))


@pytest.mark.parametrize(
    "column,value",
    [
        ("POS", True),
        ("POS", 0),
        ("Sample", "Del:1"),
        ("Sample", "Del:true:20"),
        ("REF", ""),
        ("Variant", "Substitution"),
    ],
)
def test_raw_cells_are_strict_and_replayable(column: str, value: object) -> None:
    frame = kestrel_stage_frame("raw")
    if column == "POS":
        frame["POS"] = frame["POS"].astype(object)
    frame.loc[0, column] = value
    with pytest.raises(ValueError, match=column):
        _capture(frame)


def test_capture_rejects_config_that_differs_from_frozen_baseline_policy() -> None:
    config = kestrel_config()
    config["confidence_assignment"]["reporting_floor"] = 0.2
    with pytest.raises(ValueError, match="baseline policy"):
        build_kestrel_capture(
            kestrel_stage_frame("raw"),
            _motifs(),
            kestrel_config=config,
            baseline_policy=_baseline(),
            selection=_resolve_selection(config),
            identity_component=translation_component_from_config(nomenclature_config),
            decision_profile_sha256=_SHA,
            reference_file_bytes=b"reference",
            motif_reference_file_bytes=b"motifs",
            kestrel_jar_bytes=b"jar",
            capture_policy={"schema_version": "synthetic-v1"},
        )


def test_document_rejects_unknown_fields_bool_hashes_and_changed_derived_digests() -> None:
    document = kestrel_capture_document(_capture())
    document["extra"] = True
    with pytest.raises(ValueError, match="fields"):
        decode_kestrel_capture(document)
    document = kestrel_capture_document(_capture())
    provenance = document["provenance"]
    assert isinstance(provenance, dict)
    provenance["parsed_motif_table_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="motif table"):
        decode_kestrel_capture(document)
    document = kestrel_capture_document(_capture())
    document["source_row_count"] = True
    with pytest.raises(ValueError, match="source_row_count"):
        decode_kestrel_capture(document)


def test_empty_complete_capture_is_distinct_from_final_only_evidence() -> None:
    empty = kestrel_stage_frame("raw").iloc[0:0]
    capture = _capture(empty)
    assert capture.rows == ()
    assert kestrel_capture_document(capture)["source_row_count"] == 0


def test_public_projection_revalidates_forged_typed_content() -> None:
    capture = _capture()
    with pytest.raises(ValueError, match="canonical|digest"):
        kestrel_capture_document(replace(capture, sha256="0" * 64))


def test_capture_builder_rejects_an_unvalidated_selection_at_its_public_boundary() -> None:
    with pytest.raises(ValueError, match="selection"):
        build_kestrel_capture(
            kestrel_stage_frame("raw"),
            _motifs(),
            kestrel_config=kestrel_config(),
            baseline_policy=_baseline(),
            selection=object(),  # type: ignore[arg-type]
            identity_component=translation_component_from_config(nomenclature_config),
            decision_profile_sha256=_SHA,
            reference_file_bytes=b"reference",
            motif_reference_file_bytes=b"motifs",
            kestrel_jar_bytes=b"jar",
            capture_policy={"schema_version": "synthetic-v1"},
        )


@pytest.mark.parametrize("table", ["raw", "motif"])
def test_capture_rejects_missing_sequence_cells_instead_of_stringifying_them(table: str) -> None:
    raw = kestrel_stage_frame("raw")
    motifs = _motifs()
    if table == "raw":
        raw.loc[0, "Motif_sequence"] = float("nan")
    else:
        motifs.loc[0, "Motif_sequence"] = float("nan")

    with pytest.raises(ValueError, match="Motif_sequence"):
        _capture(raw) if table == "raw" else build_kestrel_capture(
            raw,
            motifs,
            kestrel_config=kestrel_config(),
            baseline_policy=_baseline(),
            selection=_resolve_selection(kestrel_config()),
            identity_component=translation_component_from_config(nomenclature_config),
            decision_profile_sha256=_SHA,
            reference_file_bytes=b"reference",
            motif_reference_file_bytes=b"motifs",
            kestrel_jar_bytes=b"jar",
            capture_policy={"schema_version": "synthetic-v1"},
        )
