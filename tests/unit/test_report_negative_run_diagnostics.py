"""What a clean negative run says about itself: the log, the report and the region file.

An external user screening a 3,500-exome cohort asked whether their installation was
broken. Every negative sample logged three WARNINGs -- ``BED file does not exist or not
provided``, ``No Kestrel data found in pipeline summary`` and ``fastp output file not
found`` -- and its report said only that "this run produced no alignment session", while
its Provenance section claimed igv.js was "embedded in this file". None of that was a
fault: the alignment view is drawn around a Kestrel call, a negative sample has none, and
fastp runs only on FASTQ input. These tests pin that a negative run says so, at INFO,
in words that name the cause, and that the report's Provenance tells the truth.

They also pin the one real defect underneath: ``kestrel/output.bed`` left by an earlier
positive run into the same output directory survived a re-run that called nothing, and
the report then drew an alignment view at a position the current run never called.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import pandas as pd
import pytest

import vntyper
from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.cli import load_config
from vntyper.scripts import generate_report, igv_absence, report_assets, summary_steps
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.kestrel_genotyping import process_kmer_results
from vntyper.scripts.nomenclature import nomenclature_config

pytestmark = pytest.mark.unit

TEMPLATE_DIR = Path(vntyper.__file__).resolve().parent / "templates"

COVERAGE_ROW = {
    "mean": 187.2,
    "median": 80.0,
    "stdev": 285.81,
    "min": 0,
    "max": 1267,
    "region_length": 4501,
    "uncovered_bases": 255,
    "percent_uncovered": 5.67,
}

#: The row ``output_empty_result`` writes for a sample with no call.
PLACEHOLDER_ROW = {
    "Motif": "None",
    "Variant": "None",
    "POS": "None",
    "REF": "None",
    "ALT": "None",
    "Motif_sequence": "None",
    "Estimated_Depth_AlternateVariant": "None",
    "Estimated_Depth_Variant_ActiveRegion": "None",
    "Depth_Score": "None",
    "Confidence": "Negative",
}

CALL_ROW = {
    "Motifs": "X-5",
    "Motif": "5",
    "Variant": "Insertion",
    "POS": 67,
    "REF": "G",
    "ALT": "GG",
    "Motif_sequence": "GGCCACCACCCTG",
    "Estimated_Depth_AlternateVariant": 120,
    "Estimated_Depth_Variant_ActiveRegion": 12000,
    "Depth_Score": 0.01,
    "Confidence": "High_Precision",
    "Flag": "Not flagged",
}


def _write_run(output_dir: Path, kestrel_rows: list[dict] | None, *, fastq_qc: bool = False) -> Path:
    """Write a finished run's ``pipeline_summary.json`` and its Kestrel stage directory.

    Args:
        output_dir: The run directory.
        kestrel_rows: The Kestrel step's rows, or None for a run with no Kestrel step.
        fastq_qc: Whether the run recorded a FASTQ quality-control (fastp) step.

    Returns:
        Path: ``output_dir``.
    """
    steps: list[dict] = [
        {"step": summary_steps.STEP_COVERAGE, "parsed_result": {"comments": [], "data": [COVERAGE_ROW]}}
    ]
    if kestrel_rows is not None:
        steps.append({"step": summary_steps.STEP_KESTREL, "parsed_result": {"comments": [], "data": kestrel_rows}})
    if fastq_qc:
        steps.append({"step": summary_steps.STEP_FASTQ_QC, "parsed_result": {}})
    (output_dir / "pipeline_summary.json").write_text(
        json.dumps({"version": "9.9.9", "input_files": {"cram": "sample.cram"}, "steps": steps}),
        encoding="utf-8",
    )
    (output_dir / "kestrel").mkdir(exist_ok=True)
    return output_dir


def _render(output_dir: Path, **kwargs) -> str:
    """Render the report the way the pipeline does, and return its HTML."""
    log_file = output_dir / "pipeline.log"
    log_file.write_text("Pipeline execution started.\n", encoding="utf-8")
    generate_summary_report(
        output_dir=str(output_dir),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=str(log_file),
        bed_file=str(output_dir / "kestrel" / "output.bed"),
        bam_file=str(output_dir / "kestrel" / "output.bam"),
        config=load_config(None),
        **kwargs,
    )
    return (output_dir / "summary_report.html").read_text(encoding="utf-8")


def _visible_text(document: str) -> str:
    document = re.sub(r"<(script|style)\b.*?</\1>", " ", document, flags=re.DOTALL)
    document = re.sub(r"<!--.*?-->", " ", document, flags=re.DOTALL)
    return " ".join(re.sub(r"<[^>]+>", " ", document).split())


def _report_records(caplog: pytest.LogCaptureFixture, level: int) -> list[str]:
    return [
        record.getMessage()
        for record in caplog.records
        if record.name == generate_report.logger.name and record.levelno == level
    ]


@pytest.fixture
def no_igv_generator(monkeypatch: pytest.MonkeyPatch) -> list[tuple]:
    """Record every igv-reports invocation instead of running it."""
    calls: list[tuple] = []
    monkeypatch.setattr(generate_report, "run_igv_report", lambda *args, **kwargs: calls.append(args))
    return calls


# ---------------------------------------------------------------------------
# The reported case: a negative sample, a fresh output directory
# ---------------------------------------------------------------------------


def test_a_negative_run_logs_no_warning_about_its_missing_alignment_view(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    """The three WARNINGs a user grepping 3,500 logs would read as an installation fault."""
    caplog.set_level(logging.DEBUG)
    _render(_write_run(tmp_path, [PLACEHOLDER_ROW]))

    assert _report_records(caplog, logging.WARNING) == []
    assert no_igv_generator == []


def test_a_negative_run_says_at_info_why_there_is_no_alignment_view(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    caplog.set_level(logging.INFO)
    _render(_write_run(tmp_path, [PLACEHOLDER_ROW]))

    [message] = [m for m in _report_records(caplog, logging.INFO) if "alignment view" in m]
    assert "Kestrel called no variant" in message
    assert "expected" in message
    assert "not a sign of a failed step or an incomplete installation" in message


def test_a_negative_report_names_the_cause_of_its_empty_alignment_panel(
    tmp_path: Path, no_igv_generator: list[tuple]
) -> None:
    text = _visible_text(_render(_write_run(tmp_path, [PLACEHOLDER_ROW])))

    assert "Kestrel called no variant in this sample" in text
    assert "not a sign of a failed step or an incomplete installation" in text
    assert "This run produced no alignment session" not in text


def test_a_report_without_an_alignment_view_does_not_claim_igv_is_embedded(
    tmp_path: Path, no_igv_generator: list[tuple]
) -> None:
    """Provenance used to print the embedded mode's line whether or not a payload was written."""
    html = _render(_write_run(tmp_path, [PLACEHOLDER_ROW]))

    assert 'const IGV_GZ_B64 = "' not in html
    text = _visible_text(html)
    assert "embedded in this file" not in text
    assert "Alignment browser: igv.js not included (this report has no alignment view)" in text


def test_a_negative_run_logs_its_kestrel_result_as_a_result_not_as_missing_data(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    caplog.set_level(logging.INFO)
    _render(_write_run(tmp_path, [PLACEHOLDER_ROW]))

    messages = _report_records(caplog, logging.INFO)
    assert any("Kestrel called no variant for this sample" in m for m in messages)
    assert not any("No Kestrel data found" in m for m in _report_records(caplog, logging.WARNING))


def test_a_run_with_no_kestrel_rows_at_all_still_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    """Zero rows is not a negative result: a negative run writes the placeholder row."""
    caplog.set_level(logging.WARNING)
    _render(_write_run(tmp_path, []))

    assert any("No Kestrel rows" in m for m in _report_records(caplog, logging.WARNING))


# ---------------------------------------------------------------------------
# fastp: FASTQ-only, so its absence means something only for a FASTQ run
# ---------------------------------------------------------------------------


def test_a_cram_run_does_not_warn_that_fastp_output_is_missing(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    caplog.set_level(logging.INFO)
    _render(_write_run(tmp_path, [PLACEHOLDER_ROW]))

    assert not any("fastp" in m for m in _report_records(caplog, logging.WARNING))
    assert any("fastp runs only on FASTQ input" in m for m in _report_records(caplog, logging.INFO))


def test_a_fastq_run_that_lost_its_fastp_output_still_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    caplog.set_level(logging.WARNING)
    _render(_write_run(tmp_path, [PLACEHOLDER_ROW], fastq_qc=True))

    assert any("fastp output file not found" in m for m in _report_records(caplog, logging.WARNING))


def test_load_fastp_output_reports_an_unexpected_absence_below_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.DEBUG)
    assert generate_report.load_fastp_output(tmp_path / "output.json", expected=False) == {}
    assert _report_records(caplog, logging.WARNING) == []


# ---------------------------------------------------------------------------
# The states that are not normal still say so -- and say what to do
# ---------------------------------------------------------------------------


def test_a_call_whose_region_file_is_missing_warns_and_names_the_file(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    caplog.set_level(logging.WARNING)
    text = _visible_text(_render(_write_run(tmp_path, [CALL_ROW])))

    [message] = [m for m in _report_records(caplog, logging.WARNING) if "alignment view" in m]
    assert "Kestrel called a variant" in message
    assert str(tmp_path / "kestrel" / "output.bed") in message
    assert "--bed-file" in message
    assert "Kestrel called a variant, but the region file" in text


def test_a_region_file_left_by_an_earlier_run_is_not_drawn_for_a_negative_sample(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, no_igv_generator: list[tuple]
) -> None:
    """Defence in depth for directories written before the Kestrel stage owned its BED."""
    run = _write_run(tmp_path, [PLACEHOLDER_ROW])
    (run / "kestrel" / "output.bed").write_text("X-5\t66\t67\n", encoding="utf-8")
    caplog.set_level(logging.WARNING)

    text = _visible_text(_render(run))

    assert no_igv_generator == []
    [message] = [m for m in _report_records(caplog, logging.WARNING) if "alignment view" in m]
    assert "earlier run" in message
    assert "was not used" in text


def test_a_call_with_its_region_file_still_runs_the_generator(tmp_path: Path, no_igv_generator: list[tuple]) -> None:
    run = _write_run(tmp_path, [CALL_ROW])
    (run / "kestrel" / "output.bed").write_text("X-5\t66\t67\n", encoding="utf-8")

    with pytest.raises(ValueError, match="without writing its expected report"):
        _render(run)
    assert len(no_igv_generator) == 1


def test_off_mode_keeps_its_own_reason(tmp_path: Path, no_igv_generator: list[tuple]) -> None:
    text = _visible_text(_render(_write_run(tmp_path, [PLACEHOLDER_ROW]), report_igv=report_assets.REPORT_IGV_OFF))

    assert "Alignment visualisation was switched off for this run." in text
    assert "Kestrel called no variant in this sample" not in text


# ---------------------------------------------------------------------------
# The decision itself
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("mode", "bed_exists", "state", "has_call", "expected"),
    [
        (report_assets.REPORT_IGV_OFF, True, summary_steps.STEP_READ, True, igv_absence.IGV_SWITCHED_OFF),
        (report_assets.REPORT_IGV_EMBEDDED, False, summary_steps.STEP_READ, False, igv_absence.IGV_NO_KESTREL_CALL),
        (report_assets.REPORT_IGV_EMBEDDED, True, summary_steps.STEP_READ, False, igv_absence.IGV_STALE_REGION_FILE),
        (report_assets.REPORT_IGV_EMBEDDED, False, summary_steps.STEP_READ, True, igv_absence.IGV_REGION_FILE_MISSING),
        (report_assets.REPORT_IGV_EMBEDDED, False, summary_steps.STEP_ABSENT, False, igv_absence.IGV_NO_KESTREL_RESULT),
        (report_assets.REPORT_IGV_EMBEDDED, True, summary_steps.STEP_READ, True, None),
        # A summary with no Kestrel step and an explicit BED: `vntyper report --bed-file`.
        (report_assets.REPORT_IGV_SIDECAR, True, summary_steps.STEP_ABSENT, False, None),
    ],
)
def test_decide_igv_absence(tmp_path: Path, mode, bed_exists, state, has_call, expected) -> None:
    (tmp_path / "kestrel").mkdir()
    bed = tmp_path / "kestrel" / "output.bed"
    if bed_exists:
        bed.write_text("X-5\t66\t67\n", encoding="utf-8")

    decision = igv_absence.decide_igv_absence(
        report_igv=mode, bed_file=str(bed), kestrel_state=state, kestrel_has_call=has_call
    )

    assert decision.reason == expected
    assert decision.build_view is (expected is None)


def test_an_operator_region_file_is_drawn_even_for_a_negative_sample(tmp_path: Path) -> None:
    """`vntyper report --bed-file my_regions.bed` is a choice, not a leftover."""
    bed = tmp_path / "my_regions.bed"
    bed.write_text("X-5\t66\t67\n", encoding="utf-8")

    decision = igv_absence.decide_igv_absence(
        report_igv=report_assets.REPORT_IGV_EMBEDDED,
        bed_file=str(bed),
        kestrel_state=summary_steps.STEP_READ,
        kestrel_has_call=False,
    )

    assert decision.build_view


@pytest.mark.parametrize(
    ("path", "expected"),
    [("run/kestrel/output.bed", True), ("run/output.bed", False), ("run/kestrel/regions.bed", False), (None, False)],
)
def test_is_kestrel_region_file(path, expected) -> None:
    assert igv_absence.is_kestrel_region_file(path) is expected


def test_decide_igv_absence_treats_no_bed_path_as_missing() -> None:
    decision = igv_absence.decide_igv_absence(
        report_igv=report_assets.REPORT_IGV_EMBEDDED,
        bed_file=None,
        kestrel_state=summary_steps.STEP_READ,
        kestrel_has_call=True,
    )
    assert decision.reason == igv_absence.IGV_REGION_FILE_MISSING
    assert decision.level == logging.WARNING
    assert "no region file was given" in decision.message


def test_every_normal_absence_is_logged_below_warning() -> None:
    assert igv_absence.LOG_LEVELS[igv_absence.IGV_SWITCHED_OFF] == logging.INFO
    assert igv_absence.LOG_LEVELS[igv_absence.IGV_NO_KESTREL_CALL] == logging.INFO
    for reason in (
        igv_absence.IGV_STALE_REGION_FILE,
        igv_absence.IGV_REGION_FILE_MISSING,
        igv_absence.IGV_NO_KESTREL_RESULT,
    ):
        assert igv_absence.LOG_LEVELS[reason] == logging.WARNING


def test_provenance_without_a_view_says_nothing_is_embedded() -> None:
    for mode in (report_assets.REPORT_IGV_EMBEDDED, report_assets.REPORT_IGV_SIDECAR):
        assert (
            report_assets.igv_provenance(mode, view_available=False)
            == "not included (this report has no alignment view)"
        )
    assert "embedded in this file" in report_assets.igv_provenance(report_assets.REPORT_IGV_EMBEDDED)


# ---------------------------------------------------------------------------
# The Kestrel stage owns kestrel/output.bed
# ---------------------------------------------------------------------------


def _kestrel_frame(depth_alt: int, depth_region: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    motifs = nomenclature_config["motifs"]
    assert isinstance(motifs, dict)
    frame = kestrel_stage_frame(
        "raw",
        motifs="S-C",
        motif_sequence=motifs["C"] + motifs["S"],
        depth_alt=depth_alt,
        depth_region=depth_region,
        ref="G",
        alt="GG",
    )
    return frame, pd.DataFrame({"Motif": ["S"], "Motif_sequence": [motifs["S"]]})


def test_a_rerun_that_calls_nothing_removes_the_previous_runs_region_file(tmp_path: Path) -> None:
    """Reproduces the stale BED: a positive run, then a negative one into the same directory."""
    frame, motifs = _kestrel_frame(depth_alt=120, depth_region=12000)
    assert len(process_kmer_results(frame, motifs, str(tmp_path), kestrel_config())) == 1
    assert (tmp_path / "output.bed").is_file()

    frame, motifs = _kestrel_frame(depth_alt=1, depth_region=100000)
    assert process_kmer_results(frame, motifs, str(tmp_path), kestrel_config()).empty

    assert not (tmp_path / "output.bed").exists()


def test_a_rerun_that_calls_again_rewrites_the_region_file(tmp_path: Path) -> None:
    (tmp_path / "output.bed").write_text("stale-pair\t1\t2\n", encoding="utf-8")
    frame, motifs = _kestrel_frame(depth_alt=120, depth_region=12000)

    process_kmer_results(frame, motifs, str(tmp_path), kestrel_config())

    assert (tmp_path / "output.bed").read_text(encoding="utf-8") == "S-C\t66\t67\n"


def test_an_unremovable_stale_region_file_stops_the_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "output.bed").write_text("X-5\t66\t67\n", encoding="utf-8")
    original_unlink = Path.unlink

    def refuse(self: Path, *args, **kwargs) -> None:
        if self.name == "output.bed":
            raise PermissionError("read-only")
        original_unlink(self, *args, **kwargs)

    monkeypatch.setattr(Path, "unlink", refuse)
    frame, motifs = _kestrel_frame(depth_alt=1, depth_region=100000)

    with pytest.raises(RuntimeError, match="could not be removed"):
        process_kmer_results(frame, motifs, str(tmp_path), kestrel_config())


def test_a_rerun_whose_vcf_has_no_indel_removes_the_previous_runs_region_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """The stage's earliest no-call return never reaches ``process_kmer_results``."""
    from vntyper.scripts import kestrel_genotyping as kg

    (tmp_path / "output.bed").write_text("X-5\t66\t67\n", encoding="utf-8")
    vcf = tmp_path / "output.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")
    monkeypatch.setattr(kg, "generate_header", lambda reference_vntr: ["##fileformat=VCFv4.2"])
    caplog.set_level(logging.INFO, logger=kg.logger.name)

    assert kg.process_kestrel_output(str(tmp_path), vcf, "ref.fa", {}, {}) is None

    assert not (tmp_path / "output.bed").exists()
    assert "Negative" in (tmp_path / "kestrel_result.tsv").read_text(encoding="utf-8")
    # A sample with no call is the normal state: it is said at INFO, never WARNING. Other
    # warnings are the environment's (CI has no bcftools), so only the no-call ones count.
    no_call_warnings = [
        r.getMessage()
        for r in caplog.records
        if r.levelno >= logging.WARNING
        and any(text in r.getMessage() for text in ("insertion", "deletion", "empty", "output.bed", "no variant"))
    ]
    assert no_call_warnings == []
    assert any("Kestrel called no variant" in r.getMessage() for r in caplog.records)
