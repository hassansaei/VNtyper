"""Counting as its own KAnalyze step, and the guard that keeps it call-neutral (#262).

Kestrel 1.0.1 never configures the ``CountModule`` it builds -- ``getCountModule`` calls
``setKSize``, ``setTempDirName``, ``addPostCountFilterDefinition`` and ``setFreeSegment``
and none of ``setKmerThreadCount``, ``setSplitThreadCount`` or ``setThreads`` -- so
counting runs one k-mer and one split thread while being roughly 88% of its work.

Handing it a pre-built IKC is supported: ``IkcCountMap.preModuleRun`` adopts the sample's
single source when the format is ``ikc``, sets ``rmLastTemp = false`` so Kestrel never
deletes a file it was handed, and returns false so the count module never runs. That last
point is what makes the allowlist necessary: once the IKC is supplied, a count-affecting
option in ``additional_settings`` applies to Kestrel's expectations and not to the
counting that actually happened, with no error anywhere.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from vntyper.scripts import kestrel_genotyping as kg
from vntyper.scripts.kestrel_command import CALL_ONLY_OPTIONS, construct_kestrel_command
from vntyper.scripts.kestrel_counting import attempt_directory, ikc_path

pytestmark = pytest.mark.unit

USABLE_VCF = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _config():
    """Return the smallest pipeline config ``run_kestrel`` reads."""
    return {"tools": {"java_path": "java", "kanalyze": "kanalyze.jar"}}


def _settings(**overrides):
    """Return a ``kestrel_settings`` mapping with the shipped defaults plus overrides."""
    settings = {"kmer_sizes": [20], "java_memory": "12g"}
    settings.update(overrides)
    return {"kestrel_settings": settings}


def _run_kestrel(tmp_path, **overrides):
    """Invoke ``run_kestrel`` against ``tmp_path`` with placeholder inputs."""
    kwargs = {
        "vcf_path": tmp_path / "output.vcf",
        "output_dir": str(tmp_path),
        "fastq_files": ("r1.fq", "r2.fq"),
        "reference_vntr": "ref.fa",
        "kestrel_path": "kestrel.jar",
        "config": _config(),
        "sample_name": "sample",
        "threads": 8,
    }
    kwargs.update(overrides)
    return kg.run_kestrel(**kwargs)


def _recorder(tmp_path, *, count_ok=True, call_ok=True, write_vcf=True):
    """Return a ``run_command`` stand-in and the list it records commands into."""
    commands: list[str] = []

    def fake_run_command(command, log_file=None, **kwargs):
        commands.append(command)
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        Path(log_file).write_text("", encoding="utf-8")
        if " count " in command:
            # KAnalyze writes into the attempt directory the planner derived.
            ikc = ikc_path(tmp_path, 20)
            ikc.parent.mkdir(parents=True, exist_ok=True)
            ikc.write_text("ikc", encoding="utf-8")
            (ikc.parent / "segment_0.bin").write_text("seg", encoding="utf-8")
            return count_ok
        if write_vcf:
            (tmp_path / "output.vcf").write_text(USABLE_VCF, encoding="utf-8")
        return call_ok

    return fake_run_command, commands


def _stub_postprocessing(monkeypatch):
    """Replace the two stages that turn a VCF into a result."""
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)


# ---------------------------------------------------------------------------
# The two-step invocation
# ---------------------------------------------------------------------------


def test_counting_runs_first_and_kestrel_is_handed_the_ikc(tmp_path, monkeypatch):
    """The ordering is the contract: Kestrel cannot adopt a file that is not there yet."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, commands = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert len(commands) == 2
    assert " count " in commands[0] and "-m ikc" in commands[0]
    assert " count " not in commands[1]
    assert f"-f ikc {ikc_path(tmp_path, 20)}" in commands[1]


def test_kestrel_is_no_longer_given_the_fastq_operands(tmp_path, monkeypatch):
    """A supplied IKC must be the sample's sole source, or preModuleRun will not adopt it."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, commands = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert "r1.fq" not in commands[1]
    assert "r2.fq" not in commands[1]
    assert "r1.fq" in commands[0] and "r2.fq" in commands[0]


def test_the_thread_budget_reaches_the_count_command(tmp_path, monkeypatch):
    """--threads 8 allocates -d4 -l2 -t2; the call step stays single-threaded."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, commands = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path, threads=8)

    assert "-d 4 -l 2 -t 2" in commands[0]


def test_the_call_step_takes_serial_gc_and_the_count_step_does_not(tmp_path, monkeypatch):
    """Measured, and they differ on purpose.

    The count step is genuinely parallel once the budget reaches it, so it keeps G1. The
    call step is single-threaded, where G1 spawns a GC worker per core against one
    application thread.
    """
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, commands = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert "-XX:+UseSerialGC" not in commands[0]
    assert commands[1].startswith("java -Xmx12g -XX:+UseSerialGC -jar ")


def test_the_two_steps_log_to_different_files(tmp_path, monkeypatch):
    """run_command opens its log with "w", so one file would erase the count diagnostics."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    logs: list[str] = []

    def fake_run_command(command, log_file=None, **kwargs):
        logs.append(str(log_file))
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        Path(log_file).write_text("", encoding="utf-8")
        if " count " not in command:
            (tmp_path / "output.vcf").write_text(USABLE_VCF, encoding="utf-8")
        return True

    monkeypatch.setattr(kg, "run_command", fake_run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert logs == [str(tmp_path / "kanalyze_count_kmer_20.log"), str(tmp_path / "kestrel_kmer_20.log")]


# ---------------------------------------------------------------------------
# The IKC lifecycle
# ---------------------------------------------------------------------------


def test_the_attempt_directory_is_removed_after_a_successful_run(tmp_path, monkeypatch):
    """Kestrel never deletes an IKC it was handed, so nothing else will."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, _ = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert not attempt_directory(tmp_path, 20).exists()


def test_a_partial_ikc_from_a_failed_count_is_removed(tmp_path, monkeypatch):
    """A half-written IKC must never be offered to a later attempt or a later run."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, _ = _recorder(tmp_path, count_ok=False)
    monkeypatch.setattr(kg, "run_command", run_command)

    with pytest.raises(RuntimeError, match="KAnalyze"):
        _run_kestrel(tmp_path)

    assert not attempt_directory(tmp_path, 20).exists()


def test_the_ikc_is_removed_when_the_call_step_fails(tmp_path, monkeypatch):
    """The failing path matters most: nothing else runs to clean up after it."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, _ = _recorder(tmp_path, call_ok=False)
    monkeypatch.setattr(kg, "run_command", run_command)

    with pytest.raises(RuntimeError, match="Kestrel failed"):
        _run_kestrel(tmp_path)

    assert not attempt_directory(tmp_path, 20).exists()


def test_the_offloaded_segment_files_go_with_it(tmp_path, monkeypatch):
    """Removing the directory takes KAnalyze's segments too, whatever it names them."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, _ = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert list(tmp_path.glob("**/segment_*.bin")) == []


def test_a_failed_count_reports_the_count_failure_not_a_cleanup_failure(tmp_path, monkeypatch):
    """The operator must see why counting failed, not why a temporary file survived."""
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, _ = _recorder(tmp_path, count_ok=False)
    monkeypatch.setattr(kg, "run_command", run_command)

    def refuse(*args, **kwargs):
        raise OSError("cleanup denied")

    monkeypatch.setattr("vntyper.scripts.kestrel_counting.shutil.rmtree", refuse)

    with pytest.raises(RuntimeError, match="KAnalyze"):
        _run_kestrel(tmp_path)


def test_a_failed_count_never_falls_back_to_stock_kestrel(tmp_path, monkeypatch):
    """A broken counting step must fail loudly.

    An automatic fallback would turn a broken step into a silently slower run, which is
    the failure mode this codebase is least able to notice. The kill switch is
    operator-selected and nothing selects it automatically.
    """
    monkeypatch.setattr(kg, "kestrel_config", _settings())
    run_command, commands = _recorder(tmp_path, count_ok=False)
    monkeypatch.setattr(kg, "run_command", run_command)

    with pytest.raises(RuntimeError, match="KAnalyze"):
        _run_kestrel(tmp_path)

    assert len(commands) == 1, "Kestrel must not run after the count step failed"


# ---------------------------------------------------------------------------
# The operator kill switch and IKC retention
# ---------------------------------------------------------------------------


def test_split_counting_can_be_disabled_by_an_operator(tmp_path, monkeypatch):
    """The legacy single-command path, restored exactly."""
    monkeypatch.setattr(kg, "kestrel_config", _settings(split_counting=False))
    run_command, commands = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert len(commands) == 1
    assert " count " not in commands[0]
    assert "-f ikc" not in commands[0]
    assert "r1.fq r2.fq" in commands[0]


def test_keep_ikc_retains_the_attempt_directory_for_diagnosis(tmp_path, monkeypatch):
    """A successful run deletes the IKC that would be needed to diagnose a bad result."""
    monkeypatch.setattr(kg, "kestrel_config", _settings(keep_ikc=True))
    run_command, _ = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path)

    assert ikc_path(tmp_path, 20).is_file()


def test_a_replacement_config_without_any_of_the_new_keys_still_runs(tmp_path, monkeypatch):
    """--config-path replaces the whole config (trap 2), so every key needs a default."""
    monkeypatch.setattr(kg, "kestrel_config", {"kestrel_settings": {}})
    run_command, commands = _recorder(tmp_path)
    monkeypatch.setattr(kg, "run_command", run_command)
    _stub_postprocessing(monkeypatch)

    _run_kestrel(tmp_path, config={"tools": {"java_path": "java"}})

    assert len(commands) == 2, "the shipped default is split counting"
    assert "kanalyze.jar" in commands[0], "the kanalyze path needs a shipped default too"


# ---------------------------------------------------------------------------
# The additional_settings allowlist
# ---------------------------------------------------------------------------

COUNT_AFFECTING = [
    "--mincount 3",
    "--mincount=3",
    "--minsize 11",
    "--minsize=11",
    "--minmask 4",
    "--seqfilter quality:20",
    "--quality quality:20",
    "--charset Latin-1",
    "--format fastq",
    "--memcount",
    "--nomemcount",
    "--noseqfilter",
    "--temploc /tmp/x",
    "--rmikc",
    "--normikc",
    "--free",
    "--nofree",
    "-k 21",
    "--ksize 21",
    "-k21",
    "--lib /x.jar",
    "--liburl file:///x.jar",
    "-f ikc",
    "-fikc",
    "extra_sample.fastq.gz",
    "-- extra_sample.fastq.gz",
    "--ambivar --memcount",
    "--flank",
    "--ambivar=yes",
    "--maxrepeat 4",
    "--noautoflank",
]


def _kestrel_command(**overrides):
    """Build a Kestrel call command with a supplied IKC by default."""
    kwargs = {
        "kmer_size": 20,
        "kestrel_path": "k.jar",
        "reference_vntr": "r.fa",
        "output_dir": "/o",
        "fastq_files": ["/a.fastq.gz"],
        "vcf_out": "/o/out.vcf",
        "java_path": "java",
        "java_memory": "12g",
        "max_align_states": 40,
        "max_hap_states": 40,
        "log_level": "INFO",
        "sample_name": "s",
        "ikc_path": "/o/kmer_20/kestrel_kmers.ikc",
    }
    kwargs.update(overrides)
    return construct_kestrel_command(**kwargs)


@pytest.mark.parametrize("setting", COUNT_AFFECTING)
def test_settings_that_could_desynchronise_the_two_steps_are_rejected(setting):
    """Each of these would apply to Kestrel but not to the separate count step."""
    with pytest.raises(ValueError):
        _kestrel_command(additional_settings=setting)


@pytest.mark.parametrize(
    "setting",
    ["--noambivar", "--flank 5", "--flank=5", "--maxalignstates 60", "--mindiff -5", "--weight -10,1,1,1,0"],
)
def test_a_call_only_setting_is_still_accepted(setting):
    """The guard must not make the allowlisted options unusable."""
    assert _kestrel_command(additional_settings=setting).endswith(setting)


def test_a_repeated_allowlisted_option_is_accepted():
    """Kestrel takes the last one; neither reaches the counter."""
    assert _kestrel_command(additional_settings="--flank 5 --flank 7").endswith("--flank 5 --flank 7")


def test_a_forbidden_option_hidden_behind_an_allowlisted_flag_is_rejected():
    """A boolean state machine over a flag list is what makes this fail open."""
    with pytest.raises(ValueError, match="requires a value"):
        _kestrel_command(additional_settings="--flank --memcount")


def test_the_allowlist_does_not_apply_without_a_supplied_ikc():
    """Without the split there is no second command to desynchronise from.

    Enforcing it anyway would reject configurations that are perfectly valid for stock
    Kestrel, which is exactly why the guard and the split are inseparable rather than
    merely committed together.
    """
    command = _kestrel_command(ikc_path=None, additional_settings="--mincount 3")

    assert command.endswith("--mincount 3")


def test_the_sample_and_grouping_guard_still_applies_without_an_ikc():
    """The older, narrower guard is unconditional and must stay that way."""
    with pytest.raises(ValueError, match="sample or input grouping"):
        _kestrel_command(ikc_path=None, additional_settings="--filespersample 2")


@pytest.mark.parametrize("absent", ["--noautoflank", "--maxrepeat", "--mincount", "--memcount", "--lib"])
def test_the_allowlist_contains_no_option_that_should_not_be_there(absent):
    """Two of these are not Kestrel options at all; three reach the counter."""
    assert absent not in CALL_ONLY_OPTIONS


@pytest.mark.parametrize("present", ["--countrev", "--nocountrev", "--peakscan", "--flank"])
def test_genuine_call_side_options_are_allowed(present):
    """--countrev never reaches the counter: getCountModule does not read it."""
    assert present in CALL_ONLY_OPTIONS


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def test_logstderr_is_gone_from_the_call_command():
    """OptLogStderr and OptLogStdout set the same field, so --logstdout always won."""
    command = _kestrel_command()

    assert "--logstderr" not in command
    assert "--logstdout" in command


def test_the_kestrel_stage_records_which_counting_path_produced_the_result(tmp_path, monkeypatch):
    """A result must be attributable to the code that produced it."""
    from vntyper.scripts import pipeline_kestrel

    monkeypatch.setattr(kg, "kestrel_config", _settings())
    monkeypatch.setattr(pipeline_kestrel, "record_step", lambda *a, **k: None)
    summary: dict = {"steps": []}

    pipeline_kestrel.run_kestrel_stage(
        fastq_files=("r1.fq",),
        dirs={"kestrel": str(tmp_path / "kestrel")},
        config={"tools": {"kestrel": "k.jar"}, "reference_data": {"muc1_reference_vntr": "/refs/muc1.fa"}},
        sample_name="s",
        log_level=logging.INFO,
        cwd="/project",
        summary=summary,
        summary_file_path=str(tmp_path / "summary.json"),
        runner=lambda **kwargs: None,
    )

    assert summary["kestrel_counting_mode"] == "split"


def test_the_recorded_mode_follows_the_kill_switch(tmp_path, monkeypatch):
    """Otherwise the provenance would say "split" for a run that did not split."""
    from vntyper.scripts import pipeline_kestrel

    monkeypatch.setattr(kg, "kestrel_config", _settings(split_counting=False))
    monkeypatch.setattr(pipeline_kestrel, "record_step", lambda *a, **k: None)
    summary: dict = {"steps": []}

    pipeline_kestrel.run_kestrel_stage(
        fastq_files=("r1.fq",),
        dirs={"kestrel": str(tmp_path / "kestrel")},
        config={"tools": {"kestrel": "k.jar"}, "reference_data": {"muc1_reference_vntr": "/refs/muc1.fa"}},
        sample_name="s",
        log_level=logging.INFO,
        cwd="/project",
        summary=summary,
        summary_file_path=str(tmp_path / "summary.json"),
        runner=lambda **kwargs: None,
    )

    assert summary["kestrel_counting_mode"] == "internal"
