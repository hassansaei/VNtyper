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
from pathlib import Path

import pytest

from vntyper.scripts import kestrel_genotyping as kg

pytestmark = pytest.mark.unit

#: What a successful Kestrel run writes, minimally: a header the records can be read against.
#: The five sites below used to write ``"fresh\n"``, which #223 now correctly classifies as
#: unusable -- so they would exercise the failure path instead of the success path they mean.
USABLE_VCF = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


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
        vcf_path=vcf,
        output_dir=str(tmp_path),
        fastq_files=("r1.fq", "r2.fq"),
        reference_vntr="ref.fa",
        kestrel_path="kestrel.jar",
        config=_config(),
        sample_name="sample",
    )


def test_a_stale_vcf_does_not_skip_the_kestrel_run(tmp_path, monkeypatch):
    """The regression test #212 asks for: assert Kestrel ran, never that it silently didn't."""
    vcf = tmp_path / "output.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    executed = []

    def fake_run_command(command, log_file=None, **kwargs):
        executed.append(command)
        vcf.write_text(USABLE_VCF, encoding="utf-8")
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
        vcf.write_text(USABLE_VCF, encoding="utf-8")
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

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="produced no usable VCF"):
        _run(tmp_path / "output.vcf", tmp_path)

    errors = [record.getMessage() for record in caplog.records if record.levelno == logging.ERROR]
    assert any("produced no usable VCF" in message for message in errors), errors


def test_post_processing_runs_after_a_stale_vcf_is_replaced(tmp_path, monkeypatch):
    """The two statements the old ``return`` skipped are what produce kestrel_result.tsv."""
    vcf = tmp_path / "output.vcf"
    vcf.write_text("stale\n", encoding="utf-8")
    calls = []

    def fake_run_command(command, log_file=None, **kwargs):
        vcf.write_text(USABLE_VCF, encoding="utf-8")
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
            vcf.write_text(USABLE_VCF, encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: post_processing_saw.append(vcf.is_file()))

    _run(vcf, tmp_path)

    assert len(launches) == 2, "the second k-mer size never ran"
    assert post_processing_saw == [True], "post-processing did not see the VCF the run produced"
    assert vcf.is_file(), "the VCF the successful k-mer size wrote was removed"


def test_a_kestrel_invocation_that_exits_non_zero_aborts_the_run(tmp_path, monkeypatch, caplog):
    """A failed Kestrel command must raise, not fall through to the next k-mer size.

    This is the sibling of the no-VCF path above and the sharper of the two: an exit
    status the pipeline ignored would let a sample with a pathogenic variant reach the
    report as a negative, which is the failure mode #212 exists to close. The error
    names the log file, because that is the only place the Java stack trace survives.
    """
    monkeypatch.setattr(kg, "run_command", lambda *args, **kwargs: False)

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="Kestrel failed for kmer size 20"):
        _run(tmp_path / "output.vcf", tmp_path)

    assert "kestrel_kmer_20.log" in caplog.text, "the error must name the log that holds the failure"


def test_runner_receives_the_planned_command_log_criticality_and_cwd(tmp_path, monkeypatch):
    """Extracting invocation planning must not change the execution boundary."""
    vcf = tmp_path / "output.vcf"
    calls = []

    def fake_run_command(command, log_file=None, **kwargs):
        calls.append((command, log_file, kwargs))
        vcf.write_text(USABLE_VCF, encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)

    kg.run_kestrel(
        vcf_path=vcf,
        output_dir=str(tmp_path),
        fastq_files=("r1.fq", "r2.fq", "single.fq"),
        reference_vntr="ref.fa",
        kestrel_path="kestrel.jar",
        config=_config(),
        sample_name="sample",
        log_level=logging.DEBUG,
        cwd="/project/root",
    )

    assert len(calls) == 1
    command, log_file, kwargs = calls[0]
    assert "--loglevel DEBUG" in command
    assert log_file == str(tmp_path / "kestrel_kmer_20.log")
    assert kwargs == {"critical": True, "cwd": "/project/root"}
    assert "-ssample r1.fq r2.fq single.fq --hapfmt" in command


# --------------------------------------------------------------------------------------
# #223 -- a VCF that cannot be parsed must never become a confident Negative
# --------------------------------------------------------------------------------------


def test_a_headerless_vcf_with_records_does_not_manufacture_a_negative(tmp_path, monkeypatch):
    """The defect is a *manufactured Negative*, not an unraised exception.

    Kestrel exits 0 and writes a VCF that has lost its ``#CHROM`` header but still carries two
    real indel records. ``filter_vcf`` passes the records through, the derived files inherit
    the missing header, ``read_vcf_without_comments`` returns two empty frames, and
    ``output_empty_result`` writes the ``Negative`` placeholder -- with no ERROR logged at all.
    A positive converted to a negative.

    Assertion 1 is the defect and is what fails before the fix. The raise is asserted last, and
    on its own would be satisfiable by an unrelated error, so it is never the primary
    assertion -- this repository has already been bitten once by a RED test that passed for an
    unrelated reason.
    """
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20, 25]}})
    vcf = tmp_path / "output.vcf"
    records = "chr1\t155160000\t.\tC\tCG\t.\t.\tDP=100\n" * 2

    def fake_run_command(command, log_file=None, **kwargs):
        vcf.write_text("##fileformat=VCFv4.2\n" + records, encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)

    with pytest.raises(RuntimeError) as excinfo:
        _run(vcf, tmp_path)

    # 1. THE DEFECT: no manufactured result of any kind was written.
    written = sorted(path.name for path in tmp_path.rglob("*.tsv"))
    assert written == [], f"a result file was written for an unparsable VCF: {written}"

    # 2. The run failed loudly, and the message says the output was unusable rather than absent.
    assert "no usable VCF" in str(excinfo.value)


def test_every_configured_kmer_size_is_tried_before_giving_up(tmp_path, monkeypatch):
    """An unusable VCF is the same kind of event as no VCF: warn, remove it, try the next size.

    Aborting on the first one would contradict the deliberate fall-through for an exit-zero run
    that wrote nothing, which ``test_the_vcf_a_later_kmer_size_writes_is_not_removed`` already
    pins. Retrying is safe because the terminal raise still catches the case where every
    configured size fails, so no path reaches a manufactured negative either way.
    """
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20, 25]}})
    vcf = tmp_path / "output.vcf"
    launches = []
    post_processing_saw = []

    def fake_run_command(command, log_file=None, **kwargs):
        launches.append(command)
        if len(launches) == 1:
            vcf.write_text("##fileformat=VCFv4.2\nchr1\t1\t.\tC\tCG\t.\t.\tDP=1\n", encoding="utf-8")
        else:
            vcf.write_text(USABLE_VCF, encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: post_processing_saw.append(True))

    _run(vcf, tmp_path)

    assert len(launches) == 2, "the unusable VCF from k=20 did not fall through to k=25"
    assert post_processing_saw == [True], "the usable VCF from k=25 was not post-processed"


def test_an_unusable_vcf_is_removed_before_the_next_kmer_size_runs(tmp_path, monkeypatch):
    """``vcf_path`` is the same path for every k-mer size.

    A stale unusable file left in place would be re-examined by the next iteration and reported
    against the wrong k-mer size, and if that iteration wrote nothing it would look like a
    result. The pre-loop unlink exists for the same reason.
    """
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20, 25]}})
    vcf = tmp_path / "output.vcf"
    existed_at_launch = []

    def fake_run_command(command, log_file=None, **kwargs):
        existed_at_launch.append(vcf.is_file())
        vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)

    with pytest.raises(RuntimeError):
        _run(vcf, tmp_path)

    assert existed_at_launch == [False, False], "an unusable VCF survived into the next iteration"


def test_the_warning_for_an_unusable_vcf_names_the_reason(tmp_path, monkeypatch, caplog):
    """A silent negative had no ERROR at all. Whatever replaces it has to say what happened."""
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20]}})
    vcf = tmp_path / "output.vcf"

    def fake_run_command(command, log_file=None, **kwargs):
        vcf.write_text("##fileformat=VCFv4.2\nchr1\t1\t.\tC\tCG\t.\t.\tDP=1\n", encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)

    with caplog.at_level(logging.DEBUG), pytest.raises(RuntimeError):
        _run(vcf, tmp_path)

    warnings = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
    assert any("no #CHROM header line" in message for message in warnings), warnings
    assert any("1 data line(s)" in message for message in warnings), warnings


def test_a_failed_removal_still_raises_the_terminal_error_type(tmp_path, monkeypatch):
    """The removal must not become a new way to leave the loop.

    ``run_kestrel``'s contract is that it either post-processes a usable VCF or raises
    ``RuntimeError``. An unhandled ``OSError`` from the unlink would escape as a bare
    filesystem error, skipping the terminal raise entirely and handing the caller an exception
    type it has no reason to expect.
    """
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20, 25]}})
    vcf = tmp_path / "output.vcf"

    def fake_run_command(command, log_file=None, **kwargs):
        vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return True

    def refuse_unlink(self, *args, **kwargs):
        raise PermissionError("read-only filesystem")

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(Path, "unlink", refuse_unlink)

    with pytest.raises(RuntimeError, match="could not be removed"):
        _run(vcf, tmp_path)


# --------------------------------------------------------------------------------------
# #255 -- retry attempts must not leave artefacts for a later attempt to pick up
# --------------------------------------------------------------------------------------


def test_a_discarded_attempt_leaves_no_sam_for_a_later_one_to_convert(tmp_path, monkeypatch):
    """Attempt isolation: a discarded attempt must leave no artefact for a later one.

    Scoped honestly. The pinned Kestrel opens its haplotype output with Java's truncating
    `FileOutputStream`, so an attempt that reaches SAM initialisation overwrites the file --
    the common case is self-cleaning. The gap is an attempt that exits 0 *before* reaching
    it: the predecessor's SAM survives, and a later successful attempt converts whatever
    occupies that path into `output.bam`, the alignment the report's IGV track shows.

    The genotype comes from the VCF and is unaffected either way, which is why this is easy
    to miss. This test pins the isolation, not a claim about how often it matters.
    """
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20, 25]}})
    vcf = tmp_path / "output.vcf"
    sam = tmp_path / "output.sam"
    launches = []
    converted = []

    def fake_run_command(command, log_file=None, **kwargs):
        launches.append(command)
        if len(launches) == 1:
            # exits 0, writes an unusable VCF and a SAM -- the attempt that gets discarded
            vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
            sam.write_text("@HD\tVN:1.6\nFIRST-ATTEMPT\n", encoding="utf-8")
        else:
            vcf.write_text(USABLE_VCF, encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda s, o: converted.append(Path(s).exists()))
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)

    _run(vcf, tmp_path)

    assert len(launches) == 2, "the unusable attempt did not fall through to the next k-mer size"
    assert converted == [False], "the discarded attempt's SAM survived into the successful one"


def test_an_absent_sam_is_not_an_error_when_an_attempt_is_discarded(tmp_path, monkeypatch):
    """Kestrel need not have written a SAM at all. Absent is the desired state, so removing
    it must not turn a normal discard into a RuntimeError."""
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {"kmer_sizes": [20, 25]}})
    vcf = tmp_path / "output.vcf"
    launches = []

    def fake_run_command(command, log_file=None, **kwargs):
        launches.append(command)
        vcf.write_text("##fileformat=VCFv4.2\n" if len(launches) == 1 else USABLE_VCF, encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)

    _run(vcf, tmp_path)  # must not raise

    assert len(launches) == 2


# --------------------------------------------------------------------------------------
# #255 -- samtools' exit status is a result, not something to infer from a path test
# --------------------------------------------------------------------------------------


def test_a_failed_sam_to_bam_conversion_raises_instead_of_being_inferred(tmp_path, monkeypatch):
    """Success used to be inferred from `os.path.exists` on the two outputs. A failed
    `samtools view` still leaves a truncated BAM, so the SAM was deleted, the truncated BAM
    kept, and the only record was a log file nobody reads."""
    sam = tmp_path / "output.sam"
    sam.write_text("@HD\tVN:1.6\n", encoding="utf-8")

    def failing_view(command, log_file=None, **kwargs):
        (tmp_path / "output.bam").write_text("truncated", encoding="utf-8")
        return False

    monkeypatch.setattr(kg, "run_command", failing_view)

    with pytest.raises(RuntimeError, match="SAM to BAM"):
        kg.convert_sam_to_bam_and_index(str(sam), str(tmp_path))

    assert sam.exists(), "the SAM must survive a failed conversion - it is the only copy of the data"


def test_a_failed_bam_index_raises(tmp_path, monkeypatch):
    sam = tmp_path / "output.sam"
    sam.write_text("@HD\tVN:1.6\n", encoding="utf-8")
    calls = []

    def run(command, log_file=None, **kwargs):
        calls.append(command)
        (tmp_path / "output.bam").write_text("bam", encoding="utf-8")
        # Keyed off the samtools subcommand, not a substring of the whole line: pytest's
        # tmp_path embeds the test name, which contains "index", so a naive `in` test made
        # the *view* command look like the index one.
        return not command.startswith("samtools index")

    monkeypatch.setattr(kg, "run_command", run)

    with pytest.raises(RuntimeError, match="Indexing"):
        kg.convert_sam_to_bam_and_index(str(sam), str(tmp_path))

    assert sam.exists()


def test_a_successful_conversion_still_deletes_the_sam(tmp_path, monkeypatch):
    """The behaviour that must not regress: on success the SAM is cleaned up."""
    sam = tmp_path / "output.sam"
    sam.write_text("@HD\tVN:1.6\n", encoding="utf-8")

    def run(command, log_file=None, **kwargs):
        (tmp_path / "output.bam").write_text("bam", encoding="utf-8")
        (tmp_path / "output.bam.bai").write_text("bai", encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", run)

    assert kg.convert_sam_to_bam_and_index(str(sam), str(tmp_path)) == str(tmp_path / "output.bam")
    assert not sam.exists()
