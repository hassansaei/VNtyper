"""The pure flattening of a pipeline summary into table cells (#119).

``pipeline_summary.<csv|tsv>`` used to embed a step's result rows as ``"; "``-joined
``json.dumps`` text and carried nothing about the run that produced the table. These
tests pin the replacement rules on a summary shaped exactly as ``summary.start_summary``,
``record_step``, ``pipeline.py`` and ``pipeline_kestrel.py`` write it at schema 3.
"""

from __future__ import annotations

from typing import Any

import pytest

from vntyper.scripts import summary_steps
from vntyper.scripts.summary_flattening import STEP_COLUMNS, run_columns, step_rows

pytestmark = pytest.mark.unit

PROFILE_SHA256 = "5d41402abc4b2a76b9719d911017c592" * 2
ADVNTR_EVIDENCE_DIGEST = "a3f2" * 16
MODEL_SHA256 = "0c9d" * 16

KESTREL_COMMENTS = [
    "VNtyper Kestrel result",
    "VNtyper Version: 2.0.27",
    "Analysis date: 2026-09-05 10:15:42",
    "Reference file: reference/MUC1_motifs.fa",
]

#: One final Kestrel row, column names as ``tests/builders.py`` ``STAGE_COLUMNS["named"]``
#: spells them; every value is a string because ``parse_tsv`` reads text.
KESTREL_ROW = {
    "Motifs": "X-X",
    "POS": "67",
    "REF": "G",
    "ALT": "GG",
    "Motif": "X",
    "Motif_sequence": "GGCCGGGGCCACGGTGTCCCC",
    "Estimated_Depth_AlternateVariant": "10",
    "Estimated_Depth_Variant_ActiveRegion": "1000",
    "Depth_Score": "0.01",
    "Confidence": "High_Precision",
    "Flag": "Not flagged",
    "Nomenclature": "59dupC",
    "Nomenclature_Tier": "A",
}

#: Two adVNTR rows with the ``advntr_genotyping.py`` ``final_columns`` names.
ADVNTR_ROWS = [
    {
        "VID": "25561",
        "Variant": "I22_2_G_LEN1",
        "NumberOfSupportingReads": "7",
        "MeanCoverage": "42",
        "Pvalue": "0.001",
        "RU": "2",
        "POS": "22",
        "REF": "-",
        "ALT": "G",
        "Flag": "Not flagged",
    },
    {
        "VID": "25561",
        "Variant": "D58_2&D59_2",
        "NumberOfSupportingReads": "3",
        "MeanCoverage": "42",
        "Pvalue": "0.04",
        "RU": "2",
        "POS": "58",
        "REF": "CC",
        "ALT": "-",
        "Flag": "Not flagged",
    },
]

COVERAGE_ERROR = (
    "Error parsing TSV file: [Errno 2] No such file or directory: '/work/results/coverage/coverage_summary.tsv'"
)


def _step(
    name: str,
    result_file: str,
    file_type: str,
    command: str,
    parsed_result: Any,
    *,
    md5sum: str | None = "d41d8cd98f00b204e9800998ecf8427e",
    missing: bool = False,
) -> dict[str, Any]:
    """One step record as ``summary.record_step`` writes it."""
    record: dict[str, Any] = {
        "step": name,
        "start": "2026-09-05T10:15:01.000000",
        "end": "2026-09-05T10:15:02.000000",
        "command": command,
        "result_file": result_file,
        "file_type": file_type,
        "md5sum": md5sum,
        "parsed_result": parsed_result,
    }
    if missing:
        record["result_file_missing"] = True
    return record


STEPS: list[dict[str, Any]] = [
    _step(
        summary_steps.STEP_BAM_HEADER,
        "/work/results/fastq_bam_processing/pipeline_info.json",
        "json",
        "parse_header_pipeline_info(...)",
        {"assembly_text": "hg19", "assembly_contig": "chr1", "alignment_pipeline": "bwa"},
    ),
    _step(
        summary_steps.STEP_COVERAGE,
        "/work/results/coverage/coverage_summary.tsv",
        "tsv",
        "calculate_vntr_coverage(...)",
        {"comments": [COVERAGE_ERROR], "data": []},
        md5sum=None,
        missing=True,
    ),
    _step(
        summary_steps.STEP_KESTREL,
        "/work/results/kestrel/kestrel_result.tsv",
        "tsv",
        "run_kestrel(...)",
        {"comments": KESTREL_COMMENTS, "data": [KESTREL_ROW]},
    ),
    _step(
        summary_steps.STEP_ADVNTR,
        "/work/results/advntr/output_adVNTR_result.tsv",
        "tsv",
        "run_advntr(...)",
        {"comments": [], "data": ADVNTR_ROWS},
    ),
]


def _summary(steps: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """A schema-3 summary, keys in the order the run writes them."""
    return {
        "schema_version": 3,
        "decision_policy": "legacy-selection-v1",
        "advntr_evidence_digest": ADVNTR_EVIDENCE_DIGEST,
        "decision_profile_id": "vntyper-default",
        "decision_profile_revision": "1",
        "decision_profile_kind": "packaged",
        "decision_profile_source": "package",
        "decision_profile_sha256": PROFILE_SHA256,
        "decision_profile_snapshot": "provenance/decision_profile.json",
        "pipeline_start": "2026-09-05T10:15:00.000001",
        "version": "2.0.27",
        "input_files": {"bam": "/data/sample.bam"},
        "sample_name": "sample",
        "sample_name_is_explicit": False,
        "reference_assembly_requested": "hg19",
        "reference_key_used": None,
        "reference_path": None,
        "reference_source_effective": "not-required",
        "region_resolved": "chr1:155158000-155163000",
        "steps": STEPS if steps is None else steps,
        "kestrel_counting_mode": "split",
        "advntr_model": {
            "path": "reference/vntr_db_advntr/vntr_data/hg19_genic_VNTRs.db",
            "sha256": MODEL_SHA256,
            "schema_version": "v2",
            "vid": "25561",
            "genomic_interval": "chr1:155160500-155161340",
            "window_bp": 840,
            "n_segments": 14,
            "n_distinct_segments": 9,
            "max_segment_len": 60,
        },
        "pipeline_end": "2026-09-05T10:18:30.000002",
    }


# ---------------------------------------------------------------------------
# run_columns
# ---------------------------------------------------------------------------


def test_every_top_level_scalar_becomes_a_run_prefixed_text_cell() -> None:
    cells = run_columns(_summary())

    assert cells["run_schema_version"] == "3"
    assert cells["run_decision_policy"] == "legacy-selection-v1"
    assert cells["run_decision_profile_sha256"] == PROFILE_SHA256
    assert cells["run_sample_name_is_explicit"] == "False"
    assert cells["run_kestrel_counting_mode"] == "split"
    assert cells["run_pipeline_end"] == "2026-09-05T10:18:30.000002"


def test_a_null_value_is_a_blank_cell_not_the_word_none() -> None:
    cells = run_columns(_summary())

    assert "run_reference_key_used" in cells
    assert cells["run_reference_key_used"] == ""
    assert cells["run_reference_path"] == ""


def test_the_input_files_and_advntr_model_mappings_flatten_with_an_underscore() -> None:
    cells = run_columns(_summary())

    assert cells["run_input_files_bam"] == "/data/sample.bam"
    assert cells["run_advntr_model_sha256"] == MODEL_SHA256
    assert cells["run_advntr_model_window_bp"] == "840"
    assert cells["run_advntr_model_genomic_interval"] == "chr1:155160500-155161340"


def test_the_steps_list_never_becomes_a_run_column() -> None:
    assert [column for column in run_columns(_summary()) if column.startswith("run_steps")] == []
    assert [column for column in run_columns(_summary(steps=[])) if column.startswith("run_steps")] == []


def test_every_recorded_provenance_field_is_covered() -> None:
    """The spec's list: nothing the run records about itself is left out of the table."""
    cells = run_columns(_summary())
    expected = {
        "run_schema_version",
        "run_decision_policy",
        "run_advntr_evidence_digest",
        "run_decision_profile_id",
        "run_decision_profile_revision",
        "run_decision_profile_kind",
        "run_decision_profile_source",
        "run_decision_profile_sha256",
        "run_decision_profile_snapshot",
        "run_pipeline_start",
        "run_pipeline_end",
        "run_version",
        "run_input_files_bam",
        "run_sample_name",
        "run_sample_name_is_explicit",
        "run_reference_assembly_requested",
        "run_reference_key_used",
        "run_reference_path",
        "run_reference_source_effective",
        "run_region_resolved",
        "run_kestrel_counting_mode",
        "run_advntr_model_sha256",
    }

    assert expected <= set(cells)


def test_run_columns_keep_the_summary_s_own_key_order() -> None:
    columns = list(run_columns(_summary()))

    assert columns[:3] == ["run_schema_version", "run_decision_policy", "run_advntr_evidence_digest"]
    assert columns[-1] == "run_pipeline_end"


@pytest.mark.parametrize(
    ("mapping", "expected"),
    [
        pytest.param({"labels": ["a", "b"]}, {"run_labels": "a; b"}, id="scalar-list-joins"),
        pytest.param({"labels": []}, {"run_labels": ""}, id="empty-list-is-blank"),
        pytest.param({"nested": {"inner": {"leaf": 1.5}}}, {"run_nested_inner_leaf": "1.5"}, id="nested-mapping"),
        pytest.param({"rows": [{"a": 1}]}, {}, id="list-of-mappings-omitted"),
        pytest.param({"mixed": ["a", {"b": 1}]}, {}, id="mixed-list-omitted"),
        pytest.param({"flag": True, "count": 0}, {"run_flag": "True", "run_count": "0"}, id="bool-and-zero"),
    ],
)
def test_the_flattening_rules(mapping: dict[str, Any], expected: dict[str, str]) -> None:
    assert run_columns(mapping) == expected


def test_no_run_cell_is_embedded_json() -> None:
    assert [cell for cell in run_columns(_summary()).values() if '{"' in cell or cell.startswith("[")] == []


# ---------------------------------------------------------------------------
# step_rows
# ---------------------------------------------------------------------------


def _row_of(step_name: str) -> dict[str, str]:
    rows = [row for row in step_rows(_summary()) if row["step"] == step_name]
    assert len(rows) == 1, step_name
    return rows[0]


def test_one_row_per_step_in_recording_order() -> None:
    assert [row["step"] for row in step_rows(_summary())] == [
        summary_steps.STEP_BAM_HEADER,
        summary_steps.STEP_COVERAGE,
        summary_steps.STEP_KESTREL,
        summary_steps.STEP_ADVNTR,
    ]


def test_the_fixed_step_columns_carry_the_record_s_values_in_order() -> None:
    row = _row_of(summary_steps.STEP_KESTREL)

    assert list(row)[: len(STEP_COLUMNS)] == list(STEP_COLUMNS)
    assert row["start"] == "2026-09-05T10:15:01.000000"
    assert row["end"] == "2026-09-05T10:15:02.000000"
    assert row["command"] == "run_kestrel(...)"
    assert row["result_file"] == "/work/results/kestrel/kestrel_result.tsv"
    assert row["file_type"] == "tsv"
    assert row["md5sum"] == "d41d8cd98f00b204e9800998ecf8427e"
    assert row["result_file_missing"] == "False"


def test_a_missing_result_file_is_reported_true_and_its_md5_is_blank() -> None:
    row = _row_of(summary_steps.STEP_COVERAGE)

    assert row["result_file_missing"] == "True"
    assert row["md5sum"] == ""


def test_a_single_row_result_explodes_into_data_columns() -> None:
    row = _row_of(summary_steps.STEP_KESTREL)

    assert row["parsed_result_data_Motifs"] == "X-X"
    assert row["parsed_result_data_Depth_Score"] == "0.01"
    assert row["parsed_result_data_Nomenclature"] == "59dupC"
    assert {column for column in row if column.startswith("parsed_result_data_")} == {
        f"parsed_result_data_{field}" for field in KESTREL_ROW
    }
    assert "parsed_result_n_rows" not in row
    assert "parsed_result_data" not in row


def test_comment_lines_join_with_a_pipe() -> None:
    assert _row_of(summary_steps.STEP_KESTREL)["parsed_result_comments"] == (
        "VNtyper Kestrel result | VNtyper Version: 2.0.27 | "
        "Analysis date: 2026-09-05 10:15:42 | Reference file: reference/MUC1_motifs.fa"
    )
    assert _row_of(summary_steps.STEP_ADVNTR)["parsed_result_comments"] == ""


def test_a_multi_row_result_records_only_its_row_count() -> None:
    row = _row_of(summary_steps.STEP_ADVNTR)

    assert row["parsed_result_n_rows"] == "2"
    assert [column for column in row if column.startswith("parsed_result_data")] == []


def test_a_comments_only_result_records_zero_rows() -> None:
    row = _row_of(summary_steps.STEP_COVERAGE)

    assert row["parsed_result_n_rows"] == "0"
    assert row["parsed_result_comments"] == COVERAGE_ERROR


def test_a_json_typed_result_flattens_with_an_underscore() -> None:
    row = _row_of(summary_steps.STEP_BAM_HEADER)

    assert row["parsed_result_assembly_text"] == "hg19"
    assert row["parsed_result_assembly_contig"] == "chr1"
    assert row["parsed_result_alignment_pipeline"] == "bwa"
    assert "parsed_result_comments" not in row
    assert "parsed_result_n_rows" not in row


def test_the_shark_json_step_shape_flattens_to_its_four_fields() -> None:
    """The four-key payload ``tests/unit/test_pipeline_shark_step.py`` pins."""
    shark = _step(
        "SHARK Filtering",
        "/work/results/fastq_bam_processing/sample_shark_step.json",
        "json",
        "shark(...)",
        {
            "filtered_fastq_1": "/work/results/fastq_bam_processing/sample_shark_R1.fastq",
            "filtered_fastq_2": "/work/results/fastq_bam_processing/sample_shark_R2.fastq",
            "kept_reads_r1": "3",
            "kept_reads_r2": "3",
        },
    )

    (row,) = step_rows(_summary(steps=[shark]))

    assert row["parsed_result_filtered_fastq_1"] == "/work/results/fastq_bam_processing/sample_shark_R1.fastq"
    assert row["parsed_result_kept_reads_r1"] == "3"
    assert row["parsed_result_kept_reads_r2"] == "3"
    assert set(row) == set(STEP_COLUMNS) | {
        "parsed_result_filtered_fastq_1",
        "parsed_result_filtered_fastq_2",
        "parsed_result_kept_reads_r1",
        "parsed_result_kept_reads_r2",
    }


def test_an_error_result_keeps_its_message() -> None:
    failed = _step(
        summary_steps.STEP_BAM_HEADER,
        "/work/results/fastq_bam_processing/pipeline_info.xml",
        "xml",
        "parse_header_pipeline_info(...)",
        {"error": "Unsupported file type for result parsing: xml"},
    )

    (row,) = step_rows(_summary(steps=[failed]))

    assert row["parsed_result_error"] == "Unsupported file type for result parsing: xml"


def test_an_unparsed_step_has_only_the_fixed_columns() -> None:
    unparsed = _step(summary_steps.STEP_COVERAGE, "/work/results/coverage/coverage_summary.tsv", "tsv", "x(...)", None)

    (row,) = step_rows(_summary(steps=[unparsed]))

    assert set(row) == set(STEP_COLUMNS)


def test_a_summary_with_no_steps_has_no_rows() -> None:
    assert step_rows(_summary(steps=[])) == []
    assert step_rows({"version": "2.0.27"}) == []


def test_no_step_cell_is_embedded_json() -> None:
    assert [cell for row in step_rows(_summary()) for cell in row.values() if '{"' in cell] == []
