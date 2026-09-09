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


def test_canonical_hashing_preserves_trailing_tabs(tmp_path):
    """Trailing tabs in TSV files represent empty fields and must not be stripped."""
    tsv_with_tabs = b"col1\tcol2\tcol3\nval1\tval2\t\n"
    tsv_without_tabs = b"col1\tcol2\tcol3\nval1\tval2\n"

    from vntyper.scripts.report_integrity import _hash_file_content

    hash_with_tabs = _hash_file_content(tsv_with_tabs, "test.tsv")
    hash_without_tabs = _hash_file_content(tsv_without_tabs, "test.tsv")
    assert hash_with_tabs != hash_without_tabs


def test_multistage_decision_files_included(tmp_path):
    """Optional adVNTR and provenance files are hashed into the decision digest when present."""
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    (kestrel_dir / "kestrel_result.tsv").write_text("kestrel content\n")

    cov_dir = tmp_path / "coverage"
    cov_dir.mkdir()
    (cov_dir / "coverage_summary.tsv").write_text("coverage content\n")

    advntr_dir = tmp_path / "advntr"
    advntr_dir.mkdir()
    (advntr_dir / "output_adVNTR_result.tsv").write_text("advntr content\n")
    (advntr_dir / "cross_match_results.tsv").write_text("cross match content\n")

    prov_dir = tmp_path / "provenance"
    prov_dir.mkdir()
    (prov_dir / "decision_profile.json").write_text('{"profile": "v1"}\n')

    summary = {
        "version": "2.0.26",
        "sample_name": "sample1",
        "decision_profile_id": "profile1",
        "decision_profile_sha256": "0" * 64,
    }

    anchor_pipeline_summary(summary, tmp_path)
    (tmp_path / "pipeline_summary.json").write_text(json.dumps(summary))

    result = verify_report_integrity(tmp_path)
    assert result["valid"] is True

    # Tamper with adVNTR cross_match
    (advntr_dir / "cross_match_results.tsv").write_text("tampered cross match\n")
    result_tampered = verify_report_integrity(tmp_path)
    assert result_tampered["valid"] is False
    assert result_tampered["error"] == "decision_files_digest mismatch"


def test_summary_body_tamper_detected(tmp_path):
    """Tampering with summary content (e.g. sample_name or steps) is caught by pre_anchor_summary_digest."""
    kestrel_file = tmp_path / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content\n")

    cov_file = tmp_path / "coverage_summary.tsv"
    cov_file.write_text("coverage content\n")

    summary = {
        "version": "2.0.26",
        "sample_name": "patient_001",
        "steps": [{"step": "kestrel", "status": "completed"}],
    }

    anchor_pipeline_summary(summary, tmp_path)
    # Modify sample_name in summary
    summary["sample_name"] = "patient_tampered"
    (tmp_path / "pipeline_summary.json").write_text(json.dumps(summary))

    result = verify_report_integrity(tmp_path)
    assert result["valid"] is False
    assert "tampered" in str(result["error"])


def test_verify_legacy_v1_archive_directory(tmp_path):
    """Legacy v1.0 archives without pre_anchor_summary_digest and with 7-field payload verify successfully."""
    import hashlib

    from vntyper.scripts.report_integrity import _compute_dir_decision_digest_v1

    kestrel_file = tmp_path / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content\r\n")
    cov_file = tmp_path / "coverage_summary.tsv"
    cov_file.write_text("coverage content\n")

    decision_digest = _compute_dir_decision_digest_v1(tmp_path)

    run_id = "legacy-run-123"
    version = "1.0"
    tool_version = "2.0.26"
    sample_name = "sample_legacy"
    decision_profile_id = "profile1"
    decision_profile_digest = "digest1"

    payload = f"{run_id}:{version}:{tool_version}:{sample_name}:{decision_digest}:{decision_profile_id}:{decision_profile_digest}"
    integrity_digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    summary = {
        "version": tool_version,
        "sample_name": sample_name,
        "decision_profile_id": decision_profile_id,
        "decision_profile_digest": decision_profile_digest,
        "report_integrity": {
            "report_integrity_version": "1.0",
            "run_id": run_id,
            "decision_files_digest": decision_digest,
            "report_integrity_digest": integrity_digest,
        },
    }
    (tmp_path / "pipeline_summary.json").write_text(json.dumps(summary))

    result = verify_report_integrity(tmp_path)
    assert result["valid"] is True
    assert result["run_id"] == run_id
    assert result["decision_files_digest"] == decision_digest

    # Tampering with file in v1.0 fails
    kestrel_file.write_text("tampered")
    result_tampered = verify_report_integrity(tmp_path)
    assert result_tampered["valid"] is False
    assert result_tampered["error"] == "decision_files_digest mismatch"


def test_verify_legacy_v1_zip_archive(tmp_path):
    """Legacy v1.0 zip archives verify correctly."""
    import hashlib

    from vntyper.scripts.report_integrity import _compute_dir_decision_digest_v1

    base_dir = tmp_path / "run_dir"
    base_dir.mkdir()
    kestrel_file = base_dir / "kestrel_result.tsv"
    kestrel_file.write_text("kestrel content\n")
    cov_file = base_dir / "coverage_summary.tsv"
    cov_file.write_text("coverage content\n")

    decision_digest = _compute_dir_decision_digest_v1(base_dir)

    run_id = "legacy-zip-run"
    version = "1.0"
    tool_version = "2.0.26"
    sample_name = "sample_legacy_zip"
    payload = f"{run_id}:{version}:{tool_version}:{sample_name}:{decision_digest}::"
    integrity_digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    summary = {
        "version": tool_version,
        "sample_name": sample_name,
        "report_integrity": {
            "report_integrity_version": "1.0",
            "run_id": run_id,
            "decision_files_digest": decision_digest,
            "report_integrity_digest": integrity_digest,
        },
    }
    summary_file = base_dir / "pipeline_summary.json"
    summary_file.write_text(json.dumps(summary))

    zip_path = tmp_path / "legacy.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.write(kestrel_file, "run_dir/kestrel_result.tsv")
        zf.write(cov_file, "run_dir/coverage_summary.tsv")
        zf.write(summary_file, "run_dir/pipeline_summary.json")

    result = verify_report_integrity(zip_path)
    assert result["valid"] is True


def test_verify_unsupported_version(tmp_path):
    """An unknown integrity version returns an informative error."""
    summary = {
        "version": "2.0.26",
        "sample_name": "sample1",
        "report_integrity": {
            "report_integrity_version": "99.0",
            "run_id": "test",
            "decision_files_digest": "0" * 64,
            "report_integrity_digest": "0" * 64,
        },
    }
    (tmp_path / "pipeline_summary.json").write_text(json.dumps(summary))
    result = verify_report_integrity(tmp_path)
    assert result["valid"] is False
    assert "unsupported report integrity version" in result["error"]
