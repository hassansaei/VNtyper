"""The embedded alignment view, built by the real generator from a real BAM, opened offline.

Every other test of the summary report's alignment view stubs ``create_report`` out
(``conftest.IGV_REPORTS_PAGE``), and ``test_real_sidecar.py`` runs the real generator
without a BAM. So nothing drove the path a Kestrel call actually takes: a sorted,
indexed ``output.bam`` and a VCF through igv-reports, the fragments lifted into
``summary_report.html``, and igv.js expanding in the reader's browser and loading those
records. These tests do, and assert what a reader would see: the alignment track holds
the haplotype records, the one carrying the insertion is among them, and the view is
centred on the call. A sample with no call gets the authored panel and no library.
"""

from __future__ import annotations

import json
import subprocess
from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

import vntyper
from tests.browser.conftest import external_requests
from vntyper.cli import load_config
from vntyper.scripts import summary_steps
from vntyper.scripts.generate_report import generate_summary_report

pytestmark = pytest.mark.browser

TEMPLATE_DIR = Path(vntyper.__file__).resolve().parent / "templates"
PAIR = "X-5"
#: 120 bp, non-repetitive enough that igv.js draws mismatches where they are.
PAIR_SEQUENCE = ("ACGTTGCAGGCTAACCGTAGCTTGACCATGGCAATCGGTACCTTAGGCAAGTCCGATCGA" * 2)[:120]
CALL_POSITION = 67  # 1-based; the base after it is inserted.

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
CALL_ROW = {
    "Motifs": PAIR,
    "Motif": "5",
    "Variant": "Insertion",
    "POS": CALL_POSITION,
    "REF": PAIR_SEQUENCE[CALL_POSITION - 1],
    "ALT": PAIR_SEQUENCE[CALL_POSITION - 1] + "G",
    "Motif_sequence": PAIR_SEQUENCE[CALL_POSITION - 1 : CALL_POSITION + 12],
    "Estimated_Depth_AlternateVariant": 120,
    "Estimated_Depth_Variant_ActiveRegion": 12000,
    "Depth_Score": 0.01,
    "Confidence": "High_Precision",
    "Flag": "Not flagged",
}
PLACEHOLDER_ROW = {key: "None" for key in CALL_ROW if key not in {"Motifs", "Flag"}} | {"Confidence": "Negative"}

#: igv.js 3: count the records the alignment track loaded for the visible window.
LOADED_RECORDS = """async ([chrom, start, end]) => {
    const track = igvBrowser.findTracks(t => t.type === 'alignment')[0];
    const container = await track.getFeatures(chrom, start, end, 1);
    let records = 0, withInsertion = 0;
    for (const group of container.packedGroups.values())
        for (const row of group.rows)
            for (const alignment of row.alignments) {
                records += 1;
                if ((alignment.insertions || []).length) withInsertion += 1;
            }
    return {records, withInsertion};
}"""


def _run(*command: str) -> None:
    subprocess.run(command, check=True, capture_output=True)


def _write_kestrel_stage(run: Path) -> None:
    """Write a real reference, BAM, VCF and BED as the Kestrel stage leaves them."""
    kestrel = run / "kestrel"
    kestrel.mkdir(parents=True)
    reference = run / "pairs.fa"
    reference.write_text(f">{PAIR}\n{PAIR_SEQUENCE}\n", encoding="utf-8")
    _run("samtools", "faidx", str(reference))

    # Four reference haplotypes and two carrying the insertion after CALL_POSITION.
    start, length = 41, 50
    plain = PAIR_SEQUENCE[start - 1 : start - 1 + length]
    left = CALL_POSITION - start + 1
    inserted = plain[:left] + "G" + plain[left:]
    records = [(f"ref{i}", plain, f"{length}M") for i in range(4)]
    records += [(f"ins{i}", inserted, f"{left}M1I{length - left}M") for i in range(2)]
    sam = run / "output.sam"
    sam.write_text(
        f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{PAIR}\tLN:{len(PAIR_SEQUENCE)}\n"
        + "".join(
            f"{name}\t0\t{PAIR}\t{start}\t60\t{cigar}\t*\t0\t0\t{seq}\t{'I' * len(seq)}\n"
            for name, seq, cigar in records
        ),
        encoding="utf-8",
    )
    _run("samtools", "sort", "-o", str(kestrel / "output.bam"), str(sam))
    _run("samtools", "index", str(kestrel / "output.bam"))

    vcf = kestrel / "output_indel.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={PAIR},length={len(PAIR_SEQUENCE)}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"{PAIR}\t{CALL_POSITION}\t.\t{CALL_ROW['REF']}\t{CALL_ROW['ALT']}\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    (kestrel / "output.bed").write_text(f"{PAIR}\t{CALL_POSITION - 1}\t{CALL_POSITION}\n", encoding="utf-8")


def _render(run: Path, kestrel_row: dict, monkeypatch: pytest.MonkeyPatch) -> Path:
    for variable in ("HTTPS_PROXY", "HTTP_PROXY", "ALL_PROXY"):
        monkeypatch.setenv(variable, "http://127.0.0.1:9")
    monkeypatch.setenv("NO_PROXY", "")
    (run / "pipeline_summary.json").write_text(
        json.dumps(
            {
                "version": "9.9.9",
                "input_files": {"cram": "sample.cram"},
                "steps": [
                    {"step": summary_steps.STEP_COVERAGE, "parsed_result": {"comments": [], "data": [COVERAGE_ROW]}},
                    {"step": summary_steps.STEP_KESTREL, "parsed_result": {"comments": [], "data": [kestrel_row]}},
                ],
            }
        ),
        encoding="utf-8",
    )
    log = run / "pipeline.log"
    log.write_text("Pipeline execution started.\n", encoding="utf-8")
    generate_summary_report(
        output_dir=str(run),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=str(log),
        bed_file=str(run / "kestrel" / "output.bed"),
        bam_file=str(run / "kestrel" / "output.bam"),
        fasta_file=str(run / "pairs.fa"),
        vcf_file=str(run / "kestrel" / "output_indel.vcf"),
        config=load_config(None),
    )
    return run / "summary_report.html"


def test_a_real_call_opens_an_alignment_view_holding_its_records(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    open_report: Callable[..., Page],
    page_errors: dict[Page, list[str]],
) -> None:
    _write_kestrel_stage(tmp_path)
    page = open_report(_render(tmp_path, CALL_ROW, monkeypatch), offline=True)

    page.wait_for_function(
        "() => typeof igvBrowser !== 'undefined' && igvBrowser !== null "
        "&& document.getElementById('igvState') === null",
        timeout=30_000,
    )

    tracks = page.evaluate("() => igvBrowser.findTracks(() => true).map(t => [t.type, t.name])")
    # igv.js labels a track by file name minus its last extension: `output_indel` here,
    # `output_indel.vcf` for the pipeline's bgzipped file.
    assert any(kind == "variant" and name.startswith("output_indel") for kind, name in tracks)
    assert ["alignment", "output"] in tracks
    # One locus comes back as a string, several as an array.
    locus = page.evaluate("() => igvBrowser.currentLoci()")
    assert isinstance(locus, str)
    chrom, span = locus.split(":")
    first, last = (int(value.replace(",", "")) for value in span.split("-"))
    assert chrom == PAIR
    assert first <= CALL_POSITION <= last

    loaded = page.evaluate(LOADED_RECORDS, [PAIR, first, last])
    assert loaded == {"records": 6, "withInsertion": 2}
    assert page.locator("#variant_table tbody button").inner_text() == "Showing"
    assert external_requests(page) == []
    assert page_errors[page] == []


def test_a_real_no_call_report_shows_the_reason_and_no_library(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    open_report: Callable[..., Page],
    page_errors: dict[Page, list[str]],
) -> None:
    """The same stage directory minus the BED, as the Kestrel stage leaves a no-call run."""
    _write_kestrel_stage(tmp_path)
    (tmp_path / "kestrel" / "output.bed").unlink()
    page = open_report(_render(tmp_path, PLACEHOLDER_ROW, monkeypatch), offline=True)

    panel = " ".join(page.locator("#igvState").inner_text().split())
    assert panel.startswith("No alignment view: Kestrel called no variant in this sample.")
    assert "not a sign of a failed step or an incomplete installation" in panel
    assert page.evaluate("() => typeof window.igv") == "undefined"
    assert "igv.js not included (this report has no alignment view)" in page.inner_text("body")
    assert external_requests(page) == []
    assert page_errors[page] == []
