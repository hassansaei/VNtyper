"""``vntyper report`` must actually be callable (AGENTS.md trap 11).

The subcommand was dead: :func:`vntyper.scripts.cli_report.handle_report` passed
three keyword arguments ``generate_summary_report()`` does not accept
(``input_files``, ``pipeline_version``, ``mean_vntr_coverage``) and omitted the
required positional ``log_file``. Nothing caught it, because no test ever reached
the call and the only signal was a ``TypeError`` at the moment a user ran it.

The defence is a spy that **binds the real signature**. A spy accepting
``*args, **kwargs`` records a broken call just as happily as a correct one; one
that runs ``inspect.signature(generate_summary_report).bind(...)`` fails exactly
where the real call would, while still keeping the report generator itself out of
the test. The tests below drive the real ``cli.main(["report", ...])``, so the
argparse plumbing, the config fallbacks and the call site are all exercised
together.

``setup_logging()`` clears every root handler, pytest's included, so it is stubbed
out wherever ``main()`` runs here.
"""

import argparse
import inspect
import logging
from pathlib import Path

import pytest

from vntyper import cli
from vntyper.scripts import cli_report
from vntyper.scripts.generate_report import generate_summary_report

pytestmark = pytest.mark.unit


@pytest.fixture(autouse=True)
def _no_logging_reconfiguration(monkeypatch):
    """Keep ``setup_logging`` from tearing down pytest's log capture."""
    monkeypatch.setattr(cli, "setup_logging", lambda log_level=logging.INFO, log_file=None: None)


@pytest.fixture
def bound_calls(monkeypatch):
    """Replace the report generator with a spy that binds its real signature.

    Returns:
        list[inspect.BoundArguments]: Appended to once per call. Binding happens
        inside the spy, so a call the real function would reject raises
        ``TypeError`` from the spy with the same message.
    """
    signature = inspect.signature(generate_summary_report)
    calls: list[inspect.BoundArguments] = []

    def _spy(*args, **kwargs):
        calls.append(signature.bind(*args, **kwargs))

    monkeypatch.setattr(cli_report, "generate_summary_report", _spy)
    return calls


def test_the_spy_rejects_the_arguments_the_real_function_rejects(bound_calls) -> None:
    """Guard the guard: a spy that accepts anything proves nothing."""
    with pytest.raises(TypeError):
        # Deliberately wrong, which is the whole point: mypy flags it here for the
        # same reason it flagged the production call site before #179 fixed it.
        cli_report.generate_summary_report(  # type: ignore[call-arg]
            output_dir=".", template_dir=".", report_file="r.html", nonexistent=1
        )
    assert bound_calls == []


def test_report_subcommand_calls_the_generator_with_arguments_it_accepts(tmp_path, bound_calls) -> None:
    """The whole defect, end to end: ``vntyper report`` must not raise."""
    cli.main(["report", "-o", str(tmp_path)])

    assert len(bound_calls) == 1, "the report handler did not reach generate_summary_report"
    bound = bound_calls[0]
    bound.apply_defaults()
    assert Path(bound.arguments["output_dir"]) == tmp_path
    assert bound.arguments["report_file"] == "summary_report.html"
    assert bound.arguments["config"], "the handler passed no configuration; the generator raises ValueError on None"


def test_the_removed_keywords_are_not_passed(tmp_path, bound_calls) -> None:
    """``input_files``, ``pipeline_version`` and ``mean_vntr_coverage`` are read
    out of ``pipeline_summary.json`` by the generator itself, never passed in."""
    cli.main(["report", "-o", str(tmp_path)])

    passed = set(bound_calls[0].arguments)
    assert not passed & {"input_files", "pipeline_version", "mean_vntr_coverage"}


def test_the_log_file_is_passed_through(tmp_path, bound_calls) -> None:
    """``log_file`` is required and was never supplied, so the log section of
    every report generated this way would have been missing even after the
    TypeError was gone."""
    log_file = tmp_path / "elsewhere.log"
    log_file.write_text("pipeline log\n", encoding="utf-8")

    cli.main(["--log-file", str(log_file), "report", "-o", str(tmp_path)])

    assert bound_calls[0].arguments["log_file"] == str(log_file)


def test_an_existing_pipeline_log_in_the_output_dir_wins(tmp_path, bound_calls) -> None:
    """``vntyper report -o out/`` regenerates a report from a finished run, whose
    log is ``out/pipeline.log`` -- not whatever relative default the CLI resolved."""
    pipeline_log = tmp_path / "pipeline.log"
    pipeline_log.write_text("pipeline log\n", encoding="utf-8")

    cli.main(["report", "-o", str(tmp_path)])

    assert bound_calls[0].arguments["log_file"] == str(pipeline_log)


def test_an_explicit_log_file_beats_the_output_dir_default(tmp_path, bound_calls) -> None:
    """An explicit ``--log-file`` is never overridden by the probe above."""
    (tmp_path / "pipeline.log").write_text("wrong log\n", encoding="utf-8")
    explicit = tmp_path / "explicit.log"
    explicit.write_text("right log\n", encoding="utf-8")

    cli.main(["--log-file", str(explicit), "report", "-o", str(tmp_path)])

    assert bound_calls[0].arguments["log_file"] == str(explicit)


def test_the_igv_inputs_are_discovered_under_the_input_dir(tmp_path, bound_calls) -> None:
    """The handler's existing ``kestrel/output.{bam,bed}`` probe, pinned."""
    kestrel = tmp_path / "kestrel"
    kestrel.mkdir()
    (kestrel / "output.bam").write_bytes(b"")
    (kestrel / "output.bed").write_text("", encoding="utf-8")

    cli.main(["report", "-o", str(tmp_path), "--input-dir", str(tmp_path)])

    bound = bound_calls[0]
    assert bound.arguments["bam_file"] == kestrel / "output.bam"
    assert bound.arguments["bed_file"] == kestrel / "output.bed"


def test_explicit_igv_inputs_are_not_overwritten(tmp_path, bound_calls) -> None:
    """A user-supplied ``--bam-file`` must survive the probe."""
    kestrel = tmp_path / "kestrel"
    kestrel.mkdir()
    (kestrel / "output.bam").write_bytes(b"")
    chosen = tmp_path / "chosen.bam"
    chosen.write_bytes(b"")

    cli.main(["report", "-o", str(tmp_path), "--input-dir", str(tmp_path), "--bam-file", str(chosen)])

    assert bound_calls[0].arguments["bam_file"] == chosen


def test_flanking_falls_back_to_the_configured_default(tmp_path, bound_calls) -> None:
    cli.main(["report", "-o", str(tmp_path)])
    assert bound_calls[0].arguments["flanking"] == 50


def test_an_explicit_flanking_value_is_used(tmp_path, bound_calls) -> None:
    cli.main(["report", "-o", str(tmp_path), "--flanking", "120"])
    assert bound_calls[0].arguments["flanking"] == 120


def test_the_handler_no_longer_carries_a_call_arg_ignore() -> None:
    """The ``# type: ignore[call-arg]`` existed only to keep mypy green while the
    call was known-wrong. Leaving it behind would hide the next such defect."""
    source = inspect.getsource(cli_report)
    assert "generate_summary_report(" in source, "the call site vanished; this assertion would be vacuous"
    assert "type: ignore[call-arg]" not in source


def test_handle_report_is_callable_directly(tmp_path, bound_calls) -> None:
    """The handler keeps the uniform dispatch signature ``main()`` calls it with."""
    args = argparse.Namespace(
        output_dir=str(tmp_path),
        input_dir=None,
        report_file=None,
        bed_file=None,
        bam_file=None,
        reference_fasta=None,
        flanking=None,
    )
    cli_report.handle_report(
        args,
        config={"default_values": {}},
        parser=argparse.ArgumentParser(),
        log_level_value=20,
        log_file_str=None,
    )

    assert len(bound_calls) == 1
    assert bound_calls[0].arguments["log_file"] is None


# ---------------------------------------------------------------------------
# The IGV panel resolves BAM, BED and VCF from one effective run directory (#167)
# ---------------------------------------------------------------------------
#
# SPECIFICATION (REVIEW.md finding 1, HIGH -- supersedes the plan's original
# "resolve only the VCF" text): resolving only the VCF is not enough.
# `generate_report.py` only runs IGV generation `if bed_file and
# os.path.exists(bed_file)`, and before this fix `--bam-file`/`--bed-file`
# were discovered ONLY when `--input-dir` was given. So `vntyper report -o
# <run>` (no `--input-dir`) resolved no BED at all, skipped IGV generation
# entirely, and a newly-resolved VCF would have changed nothing for that
# invocation -- exactly the case `--output-dir`'s own help text describes
# ("Directory containing pipeline results"). BAM, BED and VCF are now
# resolved from one effective run directory, in the same precedence for each:
# the explicit flag, then `--input-dir`, then `--output-dir`. This mirrors
# `resolve_log_file`, which already reads `<output-dir>/pipeline.log` on
# exactly that basis.


def _kestrel_run(run_dir: Path, vcf_name: str = "output_indel.vcf.gz") -> tuple[Path, Path, Path]:
    """Lay out a standard ``kestrel/`` directory with a BAM, a BED and a VCF."""
    kestrel = run_dir / "kestrel"
    kestrel.mkdir(parents=True, exist_ok=True)
    bam = kestrel / "output.bam"
    bam.write_bytes(b"")
    bed = kestrel / "output.bed"
    bed.write_text("chr1\t1\t2\n", encoding="utf-8")
    vcf = kestrel / vcf_name
    vcf.write_bytes(b"")
    return bam, bed, vcf


def _kestrel_vcf(run_dir: Path, name: str = "output_indel.vcf.gz") -> Path:
    kestrel = run_dir / "kestrel"
    kestrel.mkdir(parents=True, exist_ok=True)
    vcf = kestrel / name
    vcf.write_bytes(b"")
    return vcf


def test_the_vcf_is_discovered_under_the_output_dir(tmp_path, bound_calls) -> None:
    vcf = _kestrel_vcf(tmp_path)
    cli.main(["report", "-o", str(tmp_path)])
    assert bound_calls[0].arguments["vcf_file"] == str(vcf)


def test_the_vcf_is_discovered_under_the_input_dir_when_one_is_given(tmp_path, bound_calls) -> None:
    run = tmp_path / "run"
    vcf = _kestrel_vcf(run)
    cli.main(["report", "-o", str(tmp_path / "out"), "--input-dir", str(run)])
    assert bound_calls[0].arguments["vcf_file"] == str(vcf)


def test_the_compressed_vcf_wins_over_the_uncompressed_one(tmp_path, bound_calls) -> None:
    """``select_best_vcf_file``'s existing preference, reached from the report path."""
    _kestrel_vcf(tmp_path, "output_indel.vcf")
    gz = _kestrel_vcf(tmp_path, "output_indel.vcf.gz")
    cli.main(["report", "-o", str(tmp_path)])
    assert bound_calls[0].arguments["vcf_file"] == str(gz)


def test_an_explicit_vcf_file_beats_the_discovered_one(tmp_path, bound_calls) -> None:
    _kestrel_vcf(tmp_path)
    explicit = tmp_path / "mine.vcf.gz"
    explicit.write_bytes(b"")
    cli.main(["report", "-o", str(tmp_path), "--vcf-file", str(explicit)])
    assert bound_calls[0].arguments["vcf_file"] == str(explicit)


def test_no_vcf_anywhere_stays_none_and_does_not_raise(tmp_path, bound_calls) -> None:
    """bcftools is optional; a missing VCF is a warning, never an error."""
    cli.main(["report", "-o", str(tmp_path)])
    assert bound_calls[0].arguments["vcf_file"] is None


def test_an_explicit_vcf_file_that_does_not_exist_is_still_passed_through(tmp_path, bound_calls) -> None:
    """The user asked for it by name; ``run_igv_report`` warns and skips the track."""
    cli.main(["report", "-o", str(tmp_path), "--vcf-file", str(tmp_path / "absent.vcf.gz")])
    assert bound_calls[0].arguments["vcf_file"] == str(tmp_path / "absent.vcf.gz")


def test_the_bam_and_bed_are_discovered_under_the_output_dir_alone(tmp_path, bound_calls) -> None:
    """The redesign's core claim: ``-o`` alone (no ``--input-dir``) now resolves the panel."""
    bam, bed, _ = _kestrel_run(tmp_path)
    cli.main(["report", "-o", str(tmp_path)])
    bound = bound_calls[0]
    assert bound.arguments["bam_file"] == bam
    assert bound.arguments["bed_file"] == bed


def test_the_input_dir_wins_over_the_output_dir_for_bam_and_bed(tmp_path, bound_calls) -> None:
    run = tmp_path / "run"
    bam, bed, _ = _kestrel_run(run)
    # A decoy kestrel/ directly under --output-dir, which must lose to --input-dir.
    _kestrel_run(tmp_path, vcf_name="output_indel.vcf")
    cli.main(["report", "-o", str(tmp_path), "--input-dir", str(run)])
    bound = bound_calls[0]
    assert bound.arguments["bam_file"] == bam
    assert bound.arguments["bed_file"] == bed


def test_the_input_dir_wins_over_the_output_dir_for_the_vcf(tmp_path, bound_calls) -> None:
    run = tmp_path / "run"
    vcf = _kestrel_vcf(run)
    _kestrel_vcf(tmp_path, "output_indel.vcf")  # a decoy under --output-dir, must lose
    cli.main(["report", "-o", str(tmp_path), "--input-dir", str(run)])
    assert bound_calls[0].arguments["vcf_file"] == str(vcf)


def test_explicit_bam_and_bed_win_over_output_dir_discovery(tmp_path, bound_calls) -> None:
    _kestrel_run(tmp_path)
    chosen_bam = tmp_path / "chosen.bam"
    chosen_bam.write_bytes(b"")
    chosen_bed = tmp_path / "chosen.bed"
    chosen_bed.write_text("", encoding="utf-8")
    cli.main(
        [
            "report",
            "-o",
            str(tmp_path),
            "--bam-file",
            str(chosen_bam),
            "--bed-file",
            str(chosen_bed),
        ]
    )
    bound = bound_calls[0]
    assert bound.arguments["bam_file"] == chosen_bam
    assert bound.arguments["bed_file"] == chosen_bed


def test_nothing_anywhere_leaves_bam_bed_and_vcf_none_without_raising(tmp_path, bound_calls) -> None:
    """``-o`` alone, no kestrel directory anywhere: the handler must not raise."""
    cli.main(["report", "-o", str(tmp_path)])
    bound = bound_calls[0]
    assert bound.arguments["bam_file"] is None
    assert bound.arguments["bed_file"] is None
    assert bound.arguments["vcf_file"] is None


def test_vntyper_report_reaches_run_igv_report_with_all_three_tracks(tmp_path, monkeypatch) -> None:
    """The review's specific instruction: prove the panel, not just the argument.

    Unlike the ``bound_calls``-based tests above, which replace
    ``generate_summary_report`` entirely, this drives the real generator so the
    assertion lands on the actual call into ``run_igv_report`` -- the function
    that decides whether the IGV panel is built at all, and the one
    ``generate_report.py:492`` gates on ``bed_file`` alone.
    """
    from vntyper.scripts import generate_report

    bam, bed, vcf = _kestrel_run(tmp_path)

    calls: list[tuple[tuple, dict]] = []
    monkeypatch.setattr(generate_report, "run_igv_report", lambda *a, **k: calls.append((a, k)))

    cli.main(["report", "-o", str(tmp_path)])

    assert len(calls) == 1, "run_igv_report must be invoked once the bed file is resolved"
    call_args, call_kwargs = calls[0]
    assert call_args[0] == bed
    assert call_args[1] == bam
    assert call_kwargs["vcf_file"] == str(vcf)
