"""Unit tests for the Data Donation endpoints (Issue #93)."""

import io
import json
import zipfile
from pathlib import Path

import pytest
from app.config import settings
from app.donations import DonationRepository
from starlette.testclient import TestClient

from vntyper.scripts.report_integrity import anchor_pipeline_summary

pytestmark = pytest.mark.unit


def _create_valid_donation_zip(sample_name="test_sample", is_positive=False) -> bytes:
    """Build an in-memory valid, stripped VNtyper result ZIP with a valid integrity anchor."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        kestrel_data = "Region\tCoverage\tMotif\nkestrel_row_1\n"
        coverage_data = "metric\tvalue\nmean\t120.5\nvntr_flank_mean_depth\t115.2\n"

        zf.writestr("kestrel_result.tsv", kestrel_data)
        zf.writestr("coverage/coverage_summary.tsv", coverage_data)

        # Temporary directory to compute genuine anchor
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            p_kd = Path(td) / "kestrel_result.tsv"
            p_kd.write_text(kestrel_data)
            p_cov_dir = Path(td) / "coverage"
            p_cov_dir.mkdir()
            p_cov = p_cov_dir / "coverage_summary.tsv"
            p_cov.write_text(coverage_data)

            summary = {
                "version": "2.0.29",
                "sample_name": sample_name,
                "coverage": {
                    "mean": 120.5,
                    "vntr_flank_mean_depth": 115.2,
                },
                "screening": {
                    "is_positive": is_positive,
                    "kestrel_result": "positive" if is_positive else "negative",
                },
            }
            anchor_pipeline_summary(summary, td)
            zf.writestr("pipeline_summary.json", json.dumps(summary))

    return buf.getvalue()


def test_donation_status_endpoint(client: TestClient) -> None:
    """Status endpoint returns feature flag status and vocabularies."""
    res = client.get("/donations/status/")
    assert res.status_code == 200
    data = res.json()
    assert "enabled" in data
    assert "controlled_vocabularies" in data
    assert "kits" in data["controlled_vocabularies"]
    assert "platforms" in data["controlled_vocabularies"]
    assert "confirmations" in data["controlled_vocabularies"]


def test_donation_submission_disabled_returns_503(client: TestClient, monkeypatch) -> None:
    """When ENABLE_DONATIONS is false, submission returns 503."""
    monkeypatch.setattr(settings, "ENABLE_DONATIONS", False)
    zip_bytes = _create_valid_donation_zip()
    res = client.post(
        "/donations/",
        files={"archive": ("result.zip", zip_bytes, "application/zip")},
        data={
            "consent": "true",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "false",
        },
    )
    assert res.status_code == 503
    assert "disabled" in res.json()["detail"].lower()


def test_donation_submission_validation_errors(client: TestClient, monkeypatch) -> None:
    """Metadata validation enforces consent, coarse dates, HPO syntax, and confirmation method."""
    monkeypatch.setattr(settings, "ENABLE_DONATIONS", True)
    zip_bytes = _create_valid_donation_zip()

    # 1. Missing / False consent
    res = client.post(
        "/donations/",
        files={"archive": ("result.zip", zip_bytes, "application/zip")},
        data={
            "consent": "false",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "false",
        },
    )
    assert res.status_code == 422

    # 2. Positive call without confirmation method
    res = client.post(
        "/donations/",
        files={"archive": ("result.zip", zip_bytes, "application/zip")},
        data={
            "consent": "true",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "true",
            "confirmation_method": "",
        },
    )
    assert res.status_code == 422

    # 3. Non-coarse date (contains day)
    res = client.post(
        "/donations/",
        files={"archive": ("result.zip", zip_bytes, "application/zip")},
        data={
            "consent": "true",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "false",
            "collection_month": "2024-05-12",
        },
    )
    assert res.status_code == 422

    # 4. Invalid HPO term
    res = client.post(
        "/donations/",
        files={"archive": ("result.zip", zip_bytes, "application/zip")},
        data={
            "consent": "true",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "false",
            "phenotype_hpo": "HP:123,INVALID_TERM",
        },
    )
    assert res.status_code == 422


def test_donation_rejects_unstripped_raw_sequence(client: TestClient, monkeypatch) -> None:
    """Server-side re-validation rejects archives containing raw sequence files."""
    monkeypatch.setattr(settings, "ENABLE_DONATIONS", True)

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("pipeline_summary.json", "{}")
        zf.writestr("sample.bam", b"BAM raw sequence data")

    res = client.post(
        "/donations/",
        files={"archive": ("unstripped.zip", buf.getvalue(), "application/zip")},
        data={
            "consent": "true",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "false",
        },
    )
    assert res.status_code == 400
    assert "forbidden extension" in res.json()["detail"].lower()


def test_donation_rejects_tampered_archive(client: TestClient, monkeypatch) -> None:
    """Tampering with decision-bearing files causes integrity verification to fail."""
    monkeypatch.setattr(settings, "ENABLE_DONATIONS", True)
    valid_bytes = _create_valid_donation_zip()

    # Modify kestrel_result.tsv inside zip
    buf = io.BytesIO(valid_bytes)
    with zipfile.ZipFile(buf, "r") as zf_in:
        names = zf_in.namelist()
        contents = {n: zf_in.read(n) for n in names}

    contents["kestrel_result.tsv"] = b"TAMPERED DATA\n"
    tampered_buf = io.BytesIO()
    with zipfile.ZipFile(tampered_buf, "w") as zf_out:
        for n, c in contents.items():
            zf_out.writestr(n, c)

    res = client.post(
        "/donations/",
        files={"archive": ("tampered.zip", tampered_buf.getvalue(), "application/zip")},
        data={
            "consent": "true",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "false",
        },
    )
    assert res.status_code == 400
    assert "integrity verification failed" in res.json()["detail"].lower()


def test_valid_donation_submission_and_aggregates(client: TestClient, monkeypatch, tmp_path) -> None:
    """Valid donation is accepted, stored, and aggregated with privacy protections."""
    monkeypatch.setattr(settings, "ENABLE_DONATIONS", True)

    # Use a temporary SQLite DB for test isolation
    test_repo = DonationRepository(db_url=str(tmp_path / "test_donations.db"))
    monkeypatch.setattr("app.main.get_donation_repo", lambda: test_repo)

    # Submit 1 genuine donation
    zip_bytes = _create_valid_donation_zip(sample_name="sample_alpha", is_positive=False)
    res = client.post(
        "/donations/",
        files={"archive": ("donation.zip", zip_bytes, "application/zip")},
        data={
            "consent": "true",
            "kit": "Twist Exome 2.0",
            "sequencing_platform": "Illumina NovaSeq 6000",
            "positive_call": "false",
            "phenotype_hpo": "HP:0000112, HP:0000083",
            "sex": "XX",
            "collection_month": "2024-05",
        },
    )
    assert res.status_code == 200, res.text
    body = res.json()
    assert body["status"] == "accepted"
    assert "donation_id" in body
    assert "run_id" in body

    # Check aggregates: with 1 sample (< 5), individual kit is suppressed into "Other (<5 samples)"
    agg_res = client.get("/donations/aggregates/")
    assert agg_res.status_code == 200
    agg = agg_res.json()
    assert agg["total_donations"] == 1
    assert agg["positive_count"] == 0
    assert agg["negative_count"] == 1
    # Minimum cell size suppresses individual Twist Exome 2.0 stats
    assert "Twist Exome 2.0" not in agg["by_kit"]
    assert "Other (<5 samples)" in agg["by_kit"]
    assert agg["by_kit"]["Other (<5 samples)"]["negative_power_sufficient"] is False

    # Seed 5 more negative samples for Twist Exome 2.0
    for i in range(5):
        test_repo.save_donation(
            {
                "run_id": f"run-{i}",
                "tool_version": "2.0.29",
                "sample_hash": f"hash-{i}",
                "decision_files_digest": "dig",
                "report_integrity_digest": "idig",
                "kit": "Twist Exome 2.0",
                "sequencing_platform": "Illumina NovaSeq 6000",
                "phenotype_hpo": ["HP:0000112"],
                "positive_call": False,
                "depth_counting_policy": "vntr_flank_mean_depth",
                "mean_coverage": 100.0,
                "flank_mean_depth": 95.0,
            }
        )

    agg_res2 = client.get("/donations/aggregates/")
    agg2 = agg_res2.json()
    assert agg2["total_donations"] == 6
    # Now Twist Exome 2.0 has 6 samples (>= 5), so it is not suppressed
    assert "Twist Exome 2.0" in agg2["by_kit"]
    twist_stats = agg2["by_kit"]["Twist Exome 2.0"]
    assert twist_stats["total_samples"] == 6
    assert twist_stats["negative_power_sufficient"] is True
    assert twist_stats["mean_coverage"] is not None
