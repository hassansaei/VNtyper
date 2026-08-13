import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from golden_cohort import artifacts, launcher  # noqa: E402

pytestmark = pytest.mark.unit


def test_artifact_and_launcher_modules_share_one_package_identity() -> None:
    assert artifacts.__package__ == launcher.__package__ == "golden_cohort"


def test_missing_artifacts_return_none(tmp_path: Path) -> None:
    missing = tmp_path / "missing"
    assert artifacts.read_tsv(missing, []) is None
    assert artifacts.read_json(missing) is None
    assert artifacts.read_delimited(missing, "\t", []) is None
    assert artifacts.read_report(missing, []) is None


def test_malformed_json_returns_none_and_logs_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    path = tmp_path / "result.json"
    path.write_text("{broken", encoding="utf-8")
    with caplog.at_level(logging.WARNING, logger=artifacts.__name__):
        assert artifacts.read_json(path) is None
    assert f"Could not parse {path}:" in caplog.text


def test_tsv_rows_are_normalized_to_the_header_width(tmp_path: Path) -> None:
    path = tmp_path / "rows.tsv"
    path.write_text("## generated\n\na\tb\n1\n2\t3\t4\n", encoding="utf-8")
    assert artifacts.read_tsv(path, []) == {
        "columns": ["a", "b"],
        "rows": [{"a": "1", "b": ""}, {"a": "2", "b": "3"}],
        "provenance": ["## generated"],
    }


def test_delimited_rows_are_normalized_to_the_header_width(tmp_path: Path) -> None:
    path = tmp_path / "rows.tsv"
    path.write_text("a\tb\n1\n2\t3\t4\n", encoding="utf-8")
    assert artifacts.read_delimited(path, "\t", []) == {
        "columns": ["a", "b"],
        "rows": [{"a": "1", "b": ""}, {"a": "2", "b": "3"}],
    }


def test_empty_delimited_file_has_no_columns_or_rows(tmp_path: Path) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("\n", encoding="utf-8")
    assert artifacts.read_delimited(path, "\t", []) == {"columns": [], "rows": []}


def test_delimited_sort_key_orders_rows_deterministically(tmp_path: Path) -> None:
    path = tmp_path / "rows.csv"
    path.write_text("Sample,value\nb,2\na,3\na,1\n", encoding="utf-8")
    result = artifacts.read_delimited(path, ",", [], sort_key="Sample")
    assert result == {
        "columns": ["Sample", "value"],
        "rows": [
            {"Sample": "a", "value": "1"},
            {"Sample": "a", "value": "3"},
            {"Sample": "b", "value": "2"},
        ],
    }


def test_a_report_written_before_the_masthead_still_reads_as_two_boxes(tmp_path: Path) -> None:
    """The gate compares one release's output against another's, so the older shape has
    to keep reading correctly - the screening box first, the cross-match box second."""
    path = tmp_path / "summary_report.html"
    path.write_text(
        """
        <p class="summary-box summary-positive">Screening positive</p>
        <p class="summary-box">Cross-match neutral</p>
        <table>
          <tr><th>Call</th><th>Value</th></tr>
          <tr><td>Kestrel</td><td>Positive</td></tr>
        </table>
        """,
        encoding="utf-8",
    )
    assert artifacts.read_report(path, []) == {
        "screening": {"text": "Screening positive", "is_positive": True},
        "cross_match": {"text": "Cross-match neutral", "is_positive": False},
        "box_count": 2,
        "tables": [{"header": ["Call", "Value"], "rows": [["Kestrel", "Positive"]]}],
    }


def test_a_report_with_a_masthead_reads_its_state_rather_than_the_cross_match_box(tmp_path: Path) -> None:
    """#242 moved the screening summary into the masthead and left one `summary-box`.
    Reading that remaining box as the screening summary would report a delta on every
    sample that is entirely the reader's, so the masthead is what is read."""
    path = tmp_path / "summary_report.html"
    path.write_text(
        """
        <header class="masthead" data-state="finding">
          <p class="headline">Kestrel detected something.</p>
          <p class="detail">And a second part.</p>
        </header>
        <p class="summary-box">Cross-match neutral</p>
        """,
        encoding="utf-8",
    )
    assert artifacts.read_report(path, []) == {
        "screening": {
            "text": "Kestrel detected something. And a second part.",
            "is_positive": True,
            "state": "finding",
        },
        "cross_match": {"text": "Cross-match neutral", "is_positive": False},
        "box_count": 1,
        "tables": [],
    }


def test_a_masthead_state_the_boolean_could_not_express_is_not_positive(tmp_path: Path) -> None:
    """Three states where the box had two, and "the run did not establish this" is the
    one the boolean could not carry. It is not a finding, and it is not a negative."""
    path = tmp_path / "summary_report.html"
    path.write_text(
        '<header class="masthead" data-state="indeterminate"><p class="headline">No summary available.</p></header>',
        encoding="utf-8",
    )
    report = artifacts.read_report(path, [])

    assert report is not None
    assert report["screening"] == {
        "text": "No summary available.",
        "is_positive": False,
        "state": "indeterminate",
    }
    assert report["cross_match"] is None


def test_donut_totals_are_numeric_strings_in_document_order() -> None:
    html = '"text": "\\uZZZZ", "text": "<b>12</b>", "text": "ignored", "text": "<b>3.5</b>"'
    assert artifacts._donut_totals(html) == ["12", "3.5"]


def test_non_list_records_pass_through_unchanged() -> None:
    records = {"Sample": "one"}
    assert artifacts._sorted_records(records, []) is records


def test_list_records_are_sorted_deterministically() -> None:
    records = [{"Sample": "b", "value": 2}, {"value": 1, "Sample": "a"}]
    assert artifacts._sorted_records(records, []) == [
        {"value": 1, "Sample": "a"},
        {"Sample": "b", "value": 2},
    ]


def test_pipeline_case_collects_result_steps_commands_and_missing_artifacts(tmp_path: Path) -> None:
    output = tmp_path / "output"
    log = tmp_path / "log"
    output.mkdir()
    log.mkdir()
    (log / "result.json").write_text(
        '{"exit_code": 0, "launch_line": "GATE-LAUNCH", "unmapped_read_set": "reads"}', encoding="utf-8"
    )
    (log / "commands.jsonl").write_text('{"command": "samtools view"}\n\n', encoding="utf-8")
    (output / "pipeline_summary.json").write_text(
        '{"steps": [{"step": "Kestrel Genotyping", "output": "volatile"}]}', encoding="utf-8"
    )

    result = artifacts.read_pipeline_case(output, log, [])

    assert result["exit_code"] == 0
    assert result["launch_line"] == "GATE-LAUNCH"
    assert result["unmapped_read_set"] == "reads"
    assert result["pipeline_steps"] == ["Kestrel Genotyping"]
    assert result["executed_commands"] == ["samtools view"]
    assert result["kestrel_result"] is None
    assert result["screening_summary"] is None


def test_pipeline_case_recovers_prior_commands_from_a_trailing_partial_record(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    output = tmp_path / "output"
    log = tmp_path / "log"
    output.mkdir()
    log.mkdir()
    commands_log = log / "commands.jsonl"
    commands_log.write_text('{"command": "first"}\n{"command": "second"}\n{"command":', encoding="utf-8")

    with caplog.at_level(logging.WARNING, logger=artifacts.__name__):
        result = artifacts.read_pipeline_case(output, log, [])

    assert result["executed_commands"] == ["first", "second"]
    assert f"Ignoring trailing partial command record in {commands_log}" in caplog.text


def test_pipeline_case_does_not_hide_a_complete_command_record_without_a_command(tmp_path: Path) -> None:
    output = tmp_path / "output"
    log = tmp_path / "log"
    output.mkdir()
    log.mkdir()
    (log / "commands.jsonl").write_text('{"other": "complete JSON"}', encoding="utf-8")

    with pytest.raises(KeyError, match="command"):
        artifacts.read_pipeline_case(output, log, [])


def test_pipeline_case_does_not_hide_a_newline_terminated_malformed_record(tmp_path: Path) -> None:
    output = tmp_path / "output"
    log = tmp_path / "log"
    output.mkdir()
    log.mkdir()
    (log / "commands.jsonl").write_text('{"command":\n', encoding="utf-8")

    with pytest.raises(ValueError):
        artifacts.read_pipeline_case(output, log, [])


def test_cohort_case_collects_tables_donuts_exports_and_listing(tmp_path: Path) -> None:
    output = tmp_path / "output"
    log = tmp_path / "log"
    output.mkdir()
    log.mkdir()
    (log / "result.json").write_text('{"exit_code": 1, "launch_line": "GATE-LAUNCH"}', encoding="utf-8")
    (output / "cohort_summary.html").write_text(
        """
        <table>
          <tr><th>Sample</th><th>Call</th></tr>
          <tr><td>b</td><td>negative</td></tr>
          <tr><td>a</td><td>positive</td></tr>
        </table>
        <script>{"values": [2, 1], "text": "<b>3</b>"}</script>
        """,
        encoding="utf-8",
    )
    (output / "cohort_kestrel.csv").write_text("Sample,Call\nb,negative\na,positive\n", encoding="utf-8")
    (output / "cohort_kestrel.json").write_text('[{"Sample": "b"}, {"Sample": "a"}]', encoding="utf-8")
    (output / "pseudonymization_table.tsv").write_text("Pseudonym\np2\np1\n", encoding="utf-8")

    result = artifacts.read_cohort_case(output, log, [])

    assert result["exit_code"] == 1
    assert result["cohort_sample_order_raw"] == ["b", "a"]
    assert result["cohort_tables"] == [{"header": ["Sample", "Call"], "rows": [["a", "positive"], ["b", "negative"]]}]
    assert result["cohort_category_counts"] == ["2, 1"]
    assert result["cohort_category_totals"] == ["3"]
    assert result["cohort_kestrel_csv"]["rows"] == [
        {"Sample": "a", "Call": "positive"},
        {"Sample": "b", "Call": "negative"},
    ]
    assert result["cohort_kestrel_json"] == [{"Sample": "a"}, {"Sample": "b"}]
    assert result["pseudonymization_table"]["rows"] == [{"Pseudonym": "p1"}, {"Pseudonym": "p2"}]
    assert result["cohort_output_files"] == [
        "cohort_kestrel.csv",
        "cohort_kestrel.json",
        "cohort_summary.html",
        "pseudonymization_table.tsv",
    ]
