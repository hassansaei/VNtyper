import json
import zipfile

import pytest

from vntyper.scripts.report_integrity import anchor_pipeline_summary, verify_report_integrity

pytestmark = pytest.mark.unit


def test_generate_anchor(tmp_path):
    kestrel_file = tmp_path / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content\r\n")

    cov_dir = tmp_path / "coverage"
    cov_dir.mkdir()
    cov_file = cov_dir / "coverage_summary.tsv"
    cov_file.write_text("coverage content")

    summary = {
        "version": "2.0",
        "sample_name": "sample1",
        "decision_profile_id": "profile1",
        "decision_profile_digest": "digest1",
    }

    anchor = anchor_pipeline_summary(summary, tmp_path)
    assert "run_id" in anchor
    assert "report_integrity_version" in anchor
    assert "decision_files_digest" in anchor
    assert "report_integrity_digest" in anchor
    assert summary["report_integrity"] == anchor


def test_verify_integrity_untouched(tmp_path):
    kestrel_file = tmp_path / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content")

    cov_dir = tmp_path / "coverage"
    cov_dir.mkdir()
    cov_file = cov_dir / "coverage_summary.tsv"
    cov_file.write_text("coverage content")

    summary = {
        "version": "2.0",
        "sample_name": "sample1",
        "decision_profile_id": "profile1",
        "decision_profile_digest": "digest1",
    }

    anchor_pipeline_summary(summary, tmp_path)

    summary_path = tmp_path / "pipeline_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f)

    result = verify_report_integrity(tmp_path)
    assert result["valid"] is True
    assert result["run_id"] == summary["report_integrity"]["run_id"]


def test_verify_tamper_detection(tmp_path):
    kestrel_file = tmp_path / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content")

    cov_file = tmp_path / "coverage_summary.tsv"
    cov_file.write_text("coverage content")

    summary = {
        "version": "2.0",
        "sample_name": "sample1",
        "decision_profile_id": "profile1",
        "decision_profile_digest": "digest1",
    }

    anchor_pipeline_summary(summary, tmp_path)
    summary_path = tmp_path / "pipeline_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f)

    # Tamper kestrel file
    kestrel_file.write_text("kestrel content modified")

    result = verify_report_integrity(tmp_path)
    assert result["valid"] is False
    assert result["error"] == "decision_files_digest mismatch"


def test_verify_tamper_run_id(tmp_path):
    kestrel_file = tmp_path / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content")

    cov_file = tmp_path / "coverage_summary.tsv"
    cov_file.write_text("coverage content")

    summary = {
        "version": "2.0",
        "sample_name": "sample1",
        "decision_profile_id": "profile1",
        "decision_profile_digest": "digest1",
    }

    anchor_pipeline_summary(summary, tmp_path)
    summary["report_integrity"]["run_id"] = "fake-run-id"

    summary_path = tmp_path / "pipeline_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f)

    result = verify_report_integrity(tmp_path)
    assert result["valid"] is False
    assert "report_integrity_digest mismatch" in result["error"]


def test_verify_zip_archive(tmp_path):
    base_dir = tmp_path / "run_dir"
    base_dir.mkdir()

    kestrel_file = base_dir / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content")

    cov_file = base_dir / "coverage_summary.tsv"
    cov_file.write_text("coverage content")

    summary = {"version": "2.0", "sample_name": "sample1"}

    anchor_pipeline_summary(summary, base_dir)

    summary_path = base_dir / "pipeline_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f)

    zip_path = tmp_path / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zip_ref:
        zip_ref.write(kestrel_file, "run_dir/kestrel_result.tsv")
        zip_ref.write(cov_file, "run_dir/coverage_summary.tsv")
        zip_ref.write(summary_path, "run_dir/pipeline_summary.json")

    result = verify_report_integrity(zip_path)
    assert result["valid"] is True


def test_hmac_verification(tmp_path):
    kestrel_file = tmp_path / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content")

    cov_file = tmp_path / "coverage_summary.tsv"
    cov_file.write_text("coverage content")

    summary = {"version": "2.0", "sample_name": "sample1"}

    anchor_pipeline_summary(summary, tmp_path, secret_key="my-secret-key")

    summary_path = tmp_path / "pipeline_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f)

    result = verify_report_integrity(tmp_path, secret_key="my-secret-key")
    assert result["valid"] is True
    assert result["signed"] is True

    result_no_key = verify_report_integrity(tmp_path)
    assert result_no_key["valid"] is False

    result_wrong_key = verify_report_integrity(tmp_path, secret_key="wrong-key")
    assert result_wrong_key["valid"] is False
