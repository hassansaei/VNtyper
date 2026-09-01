"""Calibration evidence remains readable with scripting disabled and offline."""

from pathlib import Path

import pytest
from playwright.sync_api import Browser

from vntyper.scripts.calibration_report import decode_calibration_report, write_calibration_report

pytestmark = pytest.mark.browser


def test_calibration_report_is_complete_without_javascript(browser: Browser, tmp_path: Path) -> None:
    report = decode_calibration_report(
        {
            "schema_version": "calibration-report-v1",
            "phase": "held-out",
            "profile_sha256": "a" * 64,
            "protocol_sha256": "b" * 64,
            "evidence_sha256": "c" * 64,
            "objective": "lexicographic-safety-v1",
            "tier_metrics": [{"tier": tier, "displayed": 1, "exact": 1, "wrong": 0} for tier in ("A", "B", "C")],
            "abstentions": [{"split": "locked-heldout", "reason": "record-tie", "count": 1, "rate": "1/10"}],
            "provenance": {
                "software_versions": ["VNtyper test"],
                "reference_versions": ["reference test"],
                "sample_composition": ["total=10"],
                "assays": ["capture-short-read"],
                "depths": ["30x"],
                "read_lengths": ["150 bp"],
                "independent_array_size": ["not recorded"],
                "mutation_classes": ["duplication"],
                "manifest_hashes": ["partition=" + "d" * 64],
                "access_attempts": ["locked-heldout: one"],
                "boundary_coverage": ["share=1/2"],
                "seeds": ["bootstrap seed=295"],
            },
            "statistics": {
                "intervals": ["two-sided 95% interval"],
                "roc_rows": ["threshold=1/2"],
                "pr_rows": ["precision=1 recall=1"],
                "joint_surface_rows": ["share=1/2"],
            },
            "limitations": [
                "Local append-only custody guards are not proof of independent custody.",
                "Small-n evidence.",
                "Reporting an interval is not a clinical safety claim.",
            ],
        }
    )
    destination = tmp_path / "calibration.html"
    write_calibration_report(destination, report)

    context = browser.new_context(java_script_enabled=False)
    page = context.new_page()
    page.goto(destination.as_uri())
    try:
        text = " ".join(page.locator("main").inner_text().split())
        assert "held-out" in text
        assert "Tier Displayed Exact Wrong" in text
        assert "record-tie" in text
        assert "Reporting an interval is not a clinical safety claim." in text
        assert page.locator("script").count() == 0
    finally:
        context.close()
