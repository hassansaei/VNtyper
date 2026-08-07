"""Unit tests for ``run_kestrel``'s handling of a pre-existing ``output.vcf`` (#212).

The skip was unconditional, gated on no flag, undocumented and untested. Because it
``return``ed rather than falling through, it also skipped the two statements that turn a
VCF into a result -- so re-running into a directory left by an interrupted run could
report a confident negative for a sample carrying a pathogenic variant.

Two things are pinned here:

1. A pre-existing VCF no longer short-circuits the stage. It is removed, with a WARNING
   naming it, and Kestrel runs. Deliberate reuse belongs behind the ``--resume`` flag
   proposed in #20, not behind an unflagged ``is_file()`` check.
2. A run where no k-mer size produced a VCF raises instead of returning ``None``. The
   ``break`` is reachable only inside ``if vcf_path.is_file()``, so a Kestrel invocation
   that exits 0 without writing a VCF used to fall out of the loop silently.

``run_command``, ``convert_sam_to_bam_and_index`` and ``process_kestrel_output`` are all
module attributes of ``kestrel_genotyping`` (the first is imported into it, the other two
are defined in it) and are called unqualified, so ``monkeypatch.setattr`` on the module
substitutes them for the calls ``run_kestrel`` makes.
"""

import logging

import pytest

from vntyper.scripts import kestrel_genotyping as kg

pytestmark = pytest.mark.unit


def _config():
    """Build the smallest pipeline config ``run_kestrel`` reads.

    Returns:
        dict: A config exposing only ``tools.java_path``, the one key the function
        dereferences without a default.
    """
    return {"tools": {"java_path": "java"}}


def _run(vcf, tmp_path):
    """Invoke ``run_kestrel`` with placeholder inputs.

    Args:
        vcf (pathlib.Path): The VCF path Kestrel is told to write.
        tmp_path (pathlib.Path): The output directory.

    Returns:
        None
    """
    kg.run_kestrel(
        vcf,
        str(tmp_path),
        "r1.fq",
        "r2.fq",
        "ref.fa",
        "kestrel.jar",
        _config(),
        "sample",
    )


def test_a_stale_vcf_does_not_skip_the_kestrel_run(tmp_path, monkeypatch):
    """The regression test #212 asks for: assert Kestrel ran, never that it silently didn't."""
    vcf = tmp_path / "output.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    executed = []

    def fake_run_command(command, log_file=None, **kwargs):
        executed.append(command)
        vcf.write_text("fresh\n", encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)

    _run(vcf, tmp_path)

    assert any("kestrel" in command.lower() for command in executed), "Kestrel was skipped"


def test_a_stale_vcf_is_removed_before_the_run(tmp_path, monkeypatch, caplog):
    """Unlinking makes the reason legible in the log rather than surfacing as a Java error."""
    vcf = tmp_path / "output.vcf"
    vcf.write_text("stale\n", encoding="utf-8")
    seen = {}

    def fake_run_command(command, log_file=None, **kwargs):
        seen.setdefault("existed_at_launch", vcf.is_file())
        vcf.write_text("fresh\n", encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)

    with caplog.at_level(logging.WARNING):
        _run(vcf, tmp_path)

    assert "existed_at_launch" in seen, "Kestrel never ran"
    assert seen["existed_at_launch"] is False, "the stale VCF was still present when Kestrel ran"
    warnings = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
    assert any("output.vcf" in message for message in warnings), warnings


def test_a_run_that_produces_no_vcf_raises_rather_than_returning_silently(tmp_path, monkeypatch, caplog):
    """Outcome 2 of #212, closed at the source.

    A Kestrel invocation that exits 0 without writing a VCF used to fall out of the loop
    and return None. ``record_step`` then wrote md5sum=None and an empty data list, which
    both the HTML report and cohort mode render as a negative.
    """
    monkeypatch.setattr(kg, "run_command", lambda *a, **kw: True)

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="produced no VCF"):
        _run(tmp_path / "output.vcf", tmp_path)

    errors = [record.getMessage() for record in caplog.records if record.levelno == logging.ERROR]
    assert any("produced no VCF" in message for message in errors), errors


def test_post_processing_runs_after_a_stale_vcf_is_replaced(tmp_path, monkeypatch):
    """The two statements the old ``return`` skipped are what produce kestrel_result.tsv."""
    vcf = tmp_path / "output.vcf"
    vcf.write_text("stale\n", encoding="utf-8")
    calls = []

    def fake_run_command(command, log_file=None, **kwargs):
        vcf.write_text("fresh\n", encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: calls.append("bam"))
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: calls.append("post"))

    _run(vcf, tmp_path)

    assert calls == ["bam", "post"]


def test_the_vcf_a_later_kmer_size_writes_is_not_removed(tmp_path, monkeypatch):
    """The removal must not eat a VCF this same call just produced.

    With more than one configured k-mer size the loop can iterate, so the guard has to be
    positioned where it cannot fire on a VCF written by an earlier iteration. Here the
    first k-mer size exits 0 without a VCF and the second writes one; post-processing must
    see that file and it must survive the call.
    """
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20, 25]}})
    vcf = tmp_path / "output.vcf"
    vcf.write_text("stale\n", encoding="utf-8")
    launches = []
    post_processing_saw = []

    def fake_run_command(command, log_file=None, **kwargs):
        launches.append(command)
        if len(launches) == 2:
            vcf.write_text("fresh\n", encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: post_processing_saw.append(vcf.is_file()))

    _run(vcf, tmp_path)

    assert len(launches) == 2, "the second k-mer size never ran"
    assert post_processing_saw == [True], "post-processing did not see the VCF the run produced"
    assert vcf.is_file(), "the VCF the successful k-mer size wrote was removed"
