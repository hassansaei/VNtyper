"""The real igv-reports sidecar generates and opens without a network."""

from __future__ import annotations

import subprocess
from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import external_requests
from vntyper.scripts.igv_report import run_igv_report
from vntyper.scripts.report_assets import REPORT_IGV_SIDECAR

pytestmark = pytest.mark.browser


def test_the_real_sidecar_generates_and_reads_without_external_requests(
    tmp_path: Path,
    monkeypatch,
    open_report: Callable[..., Page],
) -> None:
    """Exercise create_report, the shipped template and the vendored library together.

    The generation subprocess inherits closed-loopback proxies, and the browser aborts
    every non-file route while recording attempts. The row and cleared startup state
    prove this is a working generated viewer, not merely an inert HTML shell that makes
    no requests.
    """
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + ("A" * 6000) + "\n", encoding="utf-8")
    subprocess.run(["samtools", "faidx", str(reference)], check=True)
    sites = tmp_path / "sites.bed"
    sites.write_text("chr1\t2500\t2501\n", encoding="utf-8")
    sidecar = tmp_path / "igv_report.html"

    for variable in ("HTTPS_PROXY", "HTTP_PROXY", "ALL_PROXY"):
        monkeypatch.setenv(variable, "http://127.0.0.1:9")
    monkeypatch.setenv("NO_PROXY", "")

    run_igv_report(
        sites,
        None,
        reference,
        sidecar,
        report_igv=REPORT_IGV_SIDECAR,
    )

    page = open_report(sidecar, offline=True)
    page.wait_for_function(
        "() => typeof window.igv !== 'undefined' && "
        "document.querySelectorAll('#variant_table tbody tr').length === 1 && "
        "document.getElementById('igvState') === null",
        timeout=15_000,
    )

    assert external_requests(page) == []
    assert page.locator("#variant_table tbody tr").count() == 1
    button = page.locator("#variant_table tbody button")
    assert button.inner_text() == "Showing"
    assert "chr1" in (button.get_attribute("aria-label") or "")
