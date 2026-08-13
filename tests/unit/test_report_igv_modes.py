"""``--report-igv``: what each mode writes into the report, and what it writes beside it.

Three modes, and they differ in what a reader ends up holding:

* ``embedded`` -- the vendored igv.js travels inside the document, gzipped and
  base64-encoded. The archived file is a complete alignment browser needing neither a
  second file nor a network.
* ``sidecar`` -- the report carries no library and points at ``igv_report.html``, which
  ``create_report --standalone`` writes beside it.
* ``off`` -- no alignment browser at all, and ``create_report`` is not run.

Two invariants cut across all three, and both are asserted here rather than assumed.

**A report with no alignments carries no library.** The payload is written only from
the branch that has a session dictionary, so the common case -- a FASTQ run with no
BAM, or any run whose BED never existed -- is a ~30 KB file rather than a ~530 KB one.

**What is embedded is what was verified.** ``report_assets`` checks the decompressed
bytes against the pinned digest before encoding, and the digest it checked is printed
in the report's Provenance section. The test that matters is the round trip *out of the
rendered document*: decode the constant the template wrote, gunzip it, and digest it.
Anything that corrupted the payload on its way through Jinja2 -- autoescaping, a stray
newline, a truncation -- fails there rather than in a browser months later.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pytest

import vntyper
from vntyper.cli import load_config
from vntyper.scripts import generate_report, report_assets, summary_steps
from vntyper.scripts.generate_report import generate_summary_report

pytestmark = pytest.mark.unit

TEMPLATE_DIR = Path(vntyper.__file__).resolve().parent / "templates"

COVERAGE_ROW = {
    "mean": 250.0,
    "median": 248.0,
    "stdev": 12.5,
    "min": 100,
    "max": 400,
    "region_length": 1000,
    "uncovered_bases": 5,
    "percent_uncovered": 0.5,
}

#: The two literals igv-reports writes into its page, and which
#: ``report_formatting.extract_igv_fragments`` lifts back out. Shaped exactly like the
#: real thing so the extraction path under test is the real one: a fake that skipped it
#: would prove nothing about the branch that decides whether a payload is written.
STUB_TABLE_JSON = {
    "headers": ["unique_id", "CHROM", "POS", "REF", "ALT"],
    "rows": [["0", "chr1", "155188205", "G", "GG"]],
}
STUB_SESSION_DICTIONARY = {"0": "session0.json"}

#: The digest of a payload as a reader's browser would compute it: decode, gunzip, hash.
_PAYLOAD = re.compile(r'const IGV_GZ_B64 = "([A-Za-z0-9+/=]+)";')

#: The declaration that puts the library in the document. The checks below look for
#: **this**, not for the bare identifier ``IGV_GZ_B64``: the loader tests
#: ``typeof IGV_GZ_B64 === "undefined"`` to decide which authored state to show, so the
#: name is in every report whether or not the payload is, and asserting on the name
#: alone would fail against a correct 30 KB document. The declaration is the fact that
#: matters - it is what 497 KB of base64 is attached to.
PAYLOAD_DECLARATION = 'const IGV_GZ_B64 = "'


def _igv_reports_page() -> str:
    """A minimal page in the shape ``create_report`` writes.

    Returns:
        str: Markup carrying the container marker and both JavaScript literals.
    """
    return (
        "<html><body>\n"
        '<div id="container">\n<div id="igvDiv"></div>\n</div>\n'
        "<script>\n"
        f"const tableJson = {json.dumps(STUB_TABLE_JSON)}\n"
        f"const sessionDictionary = {json.dumps(STUB_SESSION_DICTIONARY)}\n"
        "</script>\n"
        "</body></html>\n"
    )


@pytest.fixture
def run_with_alignments(tmp_path: Path, monkeypatch) -> Path:
    """A finished run whose BED exists, with ``create_report`` stubbed out.

    igv-reports needs a real BAM and a real reference, which the unit tier has neither
    of. The stub writes the page ``create_report`` would have written, so everything
    downstream of it -- the fragment extraction, the "is there a session" decision, the
    payload -- is the shipped code path.

    Args:
        tmp_path: Pytest's per-test temporary directory.
        monkeypatch: Pytest's patcher.

    Returns:
        Path: The run directory, with a BED beside it.
    """
    (tmp_path / "pipeline_summary.json").write_text(
        json.dumps(
            {
                "version": "9.9.9",
                "input_files": {"bam": "sample.bam"},
                "steps": [
                    {
                        "step": summary_steps.STEP_COVERAGE,
                        "parsed_result": {"comments": [], "data": [COVERAGE_ROW]},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t155188100\t155188300\n", encoding="utf-8")

    def _stub(bed_file, bam_file, fasta_file, output_html, **kwargs) -> None:
        Path(output_html).write_text(_igv_reports_page(), encoding="utf-8")

    monkeypatch.setattr(generate_report, "run_igv_report", _stub)
    return tmp_path


def _render(output_dir: Path, **kwargs: Any) -> str:
    """Render the report and return the HTML.

    Args:
        output_dir: The run directory.
        **kwargs: Forwarded to ``generate_summary_report``.

    Returns:
        str: The rendered document.
    """
    bed = output_dir / "regions.bed"
    if bed.is_file():
        kwargs.setdefault("bed_file", str(bed))
    generate_summary_report(
        output_dir=str(output_dir),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=None,
        config=load_config(None),
        **kwargs,
    )
    return (output_dir / "summary_report.html").read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# The payload is written only where it is needed
# ---------------------------------------------------------------------------


def test_a_report_without_alignments_embeds_no_library(tmp_path: Path) -> None:
    """The commonest report in any cohort must not pay 497 KB for nothing.

    No BED, so no ``create_report``, so no session -- and therefore no payload, whatever
    the mode says. This is also what keeps a report with no alignment data at roughly
    30 KB rather than roughly 530 KB.
    """
    (tmp_path / "pipeline_summary.json").write_text(
        json.dumps({"version": "9.9.9", "input_files": {"bam": "s.bam"}, "steps": []}), encoding="utf-8"
    )

    html = _render(tmp_path)

    assert PAYLOAD_DECLARATION not in html
    assert len(html) < 120_000, f"a report with no alignments is {len(html):,} bytes"


def test_a_report_with_alignments_embeds_the_library(run_with_alignments: Path) -> None:
    """The other half: the branch that has a session does write it."""
    html = _render(run_with_alignments)

    assert PAYLOAD_DECLARATION in html
    assert 450_000 < len(html) < 700_000, f"a report with alignments is {len(html):,} bytes"


@pytest.mark.parametrize("mode", [report_assets.REPORT_IGV_SIDECAR, report_assets.REPORT_IGV_OFF])
def test_the_other_two_modes_embed_no_library(run_with_alignments: Path, mode: str) -> None:
    """A run *with* alignments still carries no payload unless the mode asked for one."""
    html = _render(run_with_alignments, report_igv=mode)

    assert PAYLOAD_DECLARATION not in html


def test_the_embedded_payload_survives_the_round_trip_into_the_document(run_with_alignments: Path) -> None:
    """Decode what the template wrote, gunzip it, and check it against the pinned digest.

    The payload is interpolated into a double-quoted JavaScript string under Jinja2's
    autoescaping rather than marked ``|safe``. Base64's alphabet contains none of the
    five characters autoescaping rewrites, so escaping is a no-op over it -- but that is
    a fact worth measuring rather than reasoning about, because the failure mode is a
    library that silently will not parse in a reader's browser and nowhere else.
    """
    html = _render(run_with_alignments)

    matched = _PAYLOAD.search(html)
    assert matched, "the rendered report carries no IGV_GZ_B64 constant to check"

    payload = matched.group(1)
    expanded = gzip.decompress(base64.b64decode(payload))

    assert hashlib.sha256(expanded).hexdigest() == report_assets.IGV_SHA256
    assert payload == report_assets.igv_payload(report_assets.REPORT_IGV_EMBEDDED)


# ---------------------------------------------------------------------------
# What the reader is told, in each mode
# ---------------------------------------------------------------------------


def test_the_provenance_section_names_the_library_and_its_digest(run_with_alignments: Path) -> None:
    """The reader can check the embedded library against upstream without a browser.

    In full, not truncated: a shortened digest is something a reader can glance at and
    cannot verify, which is the wrong half of the trade for a line whose only job is
    verification.
    """
    html = _render(run_with_alignments)

    assert f"igv.js {report_assets.IGV_VERSION}" in html
    assert report_assets.IGV_SHA256 in html
    assert "embedded in this file" in html


def test_a_report_with_no_library_does_not_claim_one_by_digest(run_with_alignments: Path) -> None:
    """``off`` must not print a digest for something the file does not carry."""
    html = _render(run_with_alignments, report_igv=report_assets.REPORT_IGV_OFF)

    assert report_assets.IGV_SHA256 not in html
    assert "--report-igv off" in html


@pytest.mark.parametrize(
    ("mode", "expected"),
    [
        (report_assets.REPORT_IGV_EMBEDDED, "is embedded in this file, compressed"),
        (report_assets.REPORT_IGV_SIDECAR, "The alignment browser is a separate file."),
        (report_assets.REPORT_IGV_OFF, "Alignment visualisation was switched off for this run."),
    ],
)
def test_each_mode_authors_its_own_state_in_the_markup(run_with_alignments: Path, mode: str, expected: str) -> None:
    """Each mode's alignment frame says which mode it is, before any script runs.

    Written into the markup rather than by the script, because a reader with scripting
    off gets the markup and nothing else - and because a frame that is empty until a
    script fills it is a panel-shaped hole for as long as the script takes.
    """
    html = _render(run_with_alignments, report_igv=mode)

    assert expected in html


def test_a_run_with_no_alignment_session_says_so_without_a_script(tmp_path: Path) -> None:
    """The fourth state, and the one every FASTQ-only run lands in."""
    (tmp_path / "pipeline_summary.json").write_text(
        json.dumps({"version": "9.9.9", "input_files": {"bam": "s.bam"}, "steps": []}), encoding="utf-8"
    )

    html = _render(tmp_path)

    assert "No alignment visualisation is available for this sample." in html
    assert "is embedded in this file, compressed" not in html


# ---------------------------------------------------------------------------
# The sidecar, and the second instance of the same defect
# ---------------------------------------------------------------------------


def test_the_sidecar_igv_report_is_self_contained(tmp_path: Path, monkeypatch) -> None:
    """``create_report`` without ``--standalone`` writes a page that loads igv.js from a CDN.

    That made ``igv_report.html`` - the file ``--report-igv sidecar`` hands the reader as
    *the* alignment browser, and which every mode that runs it leaves in the output
    directory - carry the identical offline defect the summary report was being fixed
    for. It went unnoticed because nobody had opened the sidecar.
    """
    captured: list[list[str]] = []

    def _capture(cmd, **kwargs):
        captured.append([str(part) for part in cmd])
        Path(cmd[cmd.index("--output") + 1]).write_text(_igv_reports_page(), encoding="utf-8")
        return None

    monkeypatch.setattr(generate_report.subprocess, "run", _capture)

    generate_report.run_igv_report(
        bed_file=tmp_path / "regions.bed",
        bam_file=None,
        fasta_file=tmp_path / "ref.fa",
        output_html=tmp_path / "igv_report.html",
        vcf_file=None,
    )

    assert captured, "create_report was never invoked"
    assert "--standalone" in captured[0], (
        f"the sidecar alignment browser is built without --standalone, so it loads igv.js from a CDN: {captured[0]}"
    )


@pytest.mark.parametrize(
    ("mode", "expect_sidecar"),
    [
        (report_assets.REPORT_IGV_EMBEDDED, True),
        (report_assets.REPORT_IGV_SIDECAR, True),
        (report_assets.REPORT_IGV_OFF, False),
    ],
)
def test_what_each_mode_leaves_in_the_output_directory(
    run_with_alignments: Path, mode: str, expect_sidecar: bool
) -> None:
    """The archive's contents, per mode.

    ``off`` does not run ``create_report`` at all - a mode that says "no alignment
    browser" and then spends the time producing one, and archives it, is not off.
    """
    _render(run_with_alignments, report_igv=mode)

    written = sorted(path.name for path in run_with_alignments.iterdir() if path.is_file())

    assert "summary_report.html" in written
    assert ("igv_report.html" in written) is expect_sidecar, f"mode={mode} left {written}"


def test_an_unknown_mode_stops_the_render(run_with_alignments: Path) -> None:
    """A typo must not silently produce a report with no alignment browser in it."""
    with pytest.raises(ValueError, match="Unknown --report-igv mode 'sidecart'"):
        _render(run_with_alignments, report_igv="sidecart")


def test_the_default_mode_is_embedded() -> None:
    """The signature's default, and the one every existing call site gets."""
    import inspect

    signature = inspect.signature(generate_summary_report)

    assert signature.parameters["report_igv"].default == report_assets.REPORT_IGV_EMBEDDED
