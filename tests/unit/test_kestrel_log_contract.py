"""Kestrel exits 0 after a fatal error; its log is what says so (#338).

The ERROR lines below are copied from the pinned ``kestrel.jar`` run against each
failure: a missing motif reference, and a missing IKC. The first reached users as "no
usable VCF" on every sample. The second writes a header-only VCF, which the VCF contract
correctly accepts as an empty result, so without this check it would become a negative.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pytest

from vntyper.scripts import kestrel_genotyping as kg
from vntyper.scripts.kestrel_counting import attempt_directory, ikc_path
from vntyper.scripts.kestrel_log_contract import (
    MAX_QUOTED_ERRORS,
    describe_kestrel_errors,
    find_kestrel_errors,
    read_kestrel_errors,
)

pytestmark = pytest.mark.unit

REFERENCE_ERROR = (
    "23:22:07 [KestrelRunner] ERROR e.g.kestrel.runner.KestrelRunner - Error reading reference "
    "sequence(s): Error reading reference sequence: IO error reading sequence source missing.fa"
)
IKC_ERRORS = [
    "23:22:07 [KestrelRunner] ERROR edu.gatech.kestrel.counter.CountMap - Unexpected error setting "
    "sample S: Error setting k-mer counts (in map post-run) for sample S: Cannot open indexed k-mer count "
    "(IKC) file: File does not exist: /x/kestrel_kmers.ikc",
    "23:22:07 [KestrelRunner] ERROR e.g.kestrel.runner.KestrelRunner - null: Error setting k-mer counts "
    "(in map post-run) for sample S",
]
HEALTHY = [
    "23:21:25 [main] INFO  edu.gatech.kestrel.clui.Main - Log level: INFO",
    "23:21:25 [KestrelRunner] INFO  e.g.kestrel.runner.KestrelRunner - Reading references",
    "23:21:25 [KestrelRunner] INFO  e.g.kestrel.runner.KestrelRunner - Processing sample: example",
]
HEADER_ONLY_VCF = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


# ---------------------------------------------------------------------------
# find_kestrel_errors: which lines count
# ---------------------------------------------------------------------------


def test_a_healthy_log_has_no_errors():
    assert find_kestrel_errors(HEALTHY) == []


def test_the_missing_reference_error_is_found():
    assert find_kestrel_errors([*HEALTHY[:2], REFERENCE_ERROR + "\n"]) == [REFERENCE_ERROR]


def test_every_error_line_is_returned_in_order():
    assert find_kestrel_errors([HEALTHY[0], *IKC_ERRORS, HEALTHY[2]]) == IKC_ERRORS


def test_a_warning_is_not_an_error():
    """Kestrel logs WARN before the ERROR for a missing IKC; only the ERROR decides."""
    warn = "23:22:07 [KestrelRunner] WARN  edu.gatech.kestrel.counter.CountMap - Error in post-module run"
    assert find_kestrel_errors([warn]) == []


@pytest.mark.parametrize(
    "line",
    [
        "23:21:25 [KestrelRunner] INFO  e.g.kestrel.runner.KestrelRunner - Processing sample: ERROR",
        "ERROR something printed without logback formatting",
        "[KestrelRunner] ERROR e.g.kestrel.runner.KestrelRunner - no timestamp",
        "23:21:25 KestrelRunner ERROR no thread bracket",
    ],
)
def test_the_word_error_elsewhere_does_not_match(line):
    """Anchoring on timestamp and thread keeps a sample or path named ERROR from failing a run."""
    assert find_kestrel_errors([line]) == []


def test_a_millisecond_timestamp_is_accepted():
    """A logback pattern with milliseconds must not silently disable the check."""
    line = "23:22:07.123 [KestrelRunner] ERROR e.g.kestrel.runner.KestrelRunner - boom"
    assert find_kestrel_errors([line]) == [line]


# ---------------------------------------------------------------------------
# read_kestrel_errors: the file boundary
# ---------------------------------------------------------------------------


def test_reading_returns_the_errors_in_the_file(tmp_path):
    log = tmp_path / "kestrel_kmer_20.log"
    log.write_text("\n".join([*HEALTHY, REFERENCE_ERROR]) + "\n", encoding="utf-8")
    assert read_kestrel_errors(log) == [REFERENCE_ERROR]


def test_a_missing_log_reports_no_errors(tmp_path):
    assert read_kestrel_errors(tmp_path / "absent.log") == []


def test_undecodable_bytes_do_not_hide_an_error(tmp_path):
    log = tmp_path / "kestrel.log"
    log.write_bytes(b"\xff\xfe garbage\n" + REFERENCE_ERROR.encode() + b"\n")
    assert read_kestrel_errors(log) == [REFERENCE_ERROR]


def test_an_unreadable_log_fails_closed(tmp_path, caplog):
    """A log that exists but cannot be read cannot vouch for an exit 0."""
    unreadable = tmp_path / "is_a_directory.log"
    unreadable.mkdir()
    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError, match="Could not read the Kestrel log"):
        read_kestrel_errors(unreadable)
    assert "Could not read the Kestrel log" in caplog.text


# ---------------------------------------------------------------------------
# describe_kestrel_errors: the operator-facing message
# ---------------------------------------------------------------------------


def test_the_message_names_the_cause_the_kmer_and_the_log():
    msg = describe_kestrel_errors([REFERENCE_ERROR], kmer_size=20, log_file="/out/kestrel/kestrel_kmer_20.log")
    assert "1 error(s) for k-mer size 20 but exited 0" in msg
    assert "Error reading reference sequence(s)" in msg
    assert msg.endswith("Full log: /out/kestrel/kestrel_kmer_20.log")
    assert "more" not in msg


def test_the_message_quotes_at_most_the_cap_and_counts_the_rest():
    errors = [f"23:00:0{i} [T] ERROR x - cause {i}" for i in range(MAX_QUOTED_ERRORS + 2)]
    msg = describe_kestrel_errors(errors, kmer_size=25, log_file="k.log")
    assert f"{MAX_QUOTED_ERRORS + 2} error(s)" in msg
    for i in range(MAX_QUOTED_ERRORS):
        assert f"cause {i}" in msg
    assert f"cause {MAX_QUOTED_ERRORS}" not in msg
    assert "... and 2 more" in msg


def test_exactly_the_cap_has_no_remainder_line():
    errors = [f"23:00:0{i} [T] ERROR x - cause {i}" for i in range(MAX_QUOTED_ERRORS)]
    assert "more" not in describe_kestrel_errors(errors, kmer_size=20, log_file="k.log")


def test_describing_no_errors_is_a_caller_bug():
    with pytest.raises(ValueError, match="at least one error line"):
        describe_kestrel_errors([], kmer_size=20, log_file="k.log")


# ---------------------------------------------------------------------------
# run_kestrel: an exit 0 with a logged ERROR never becomes a result
# ---------------------------------------------------------------------------


def _config():
    return {"tools": {"java_path": "java", "kanalyze": "kanalyze.jar", "samtools": "samtools"}}


def _fake_run_command(tmp_path, *, call_log_lines, vcf_text):
    """KAnalyze succeeds; Kestrel exits 0, writes ``call_log_lines`` and maybe a VCF."""
    calls: list[str] = []

    def run_command(command, log_file=None, **kwargs):
        calls.append(command)
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        if " count " in command:
            Path(log_file).write_text("", encoding="utf-8")
            kmer_size = int(re.search(r" -k (\d+) ", command).group(1))
            ikc_path(tmp_path, kmer_size).write_text("ikc", encoding="utf-8")
            return True
        Path(log_file).write_text("\n".join(call_log_lines) + "\n", encoding="utf-8")
        if vcf_text is not None:
            (tmp_path / "output.vcf").write_text(vcf_text, encoding="utf-8")
        return True

    return run_command, calls


def _run(tmp_path, monkeypatch, *, call_log_lines, vcf_text, kmer_sizes=(20,)):
    run_command, calls = _fake_run_command(tmp_path, call_log_lines=call_log_lines, vcf_text=vcf_text)
    processed: list[object] = []
    monkeypatch.setattr(kg, "run_command", run_command)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: processed.append(a))
    kwargs = {
        "vcf_path": tmp_path / "output.vcf",
        "output_dir": str(tmp_path),
        "fastq_files": ("r1.fq", "r2.fq"),
        "reference_vntr": "ref.fa",
        "kestrel_path": "kestrel.jar",
        "config": _config(),
        "sample_name": "sample",
        "threads": 4,
        "runtime_component": {"kestrel_settings": {"kmer_sizes": list(kmer_sizes), "java_memory": "12g"}},
    }
    return kwargs, calls, processed


def test_a_header_only_vcf_with_a_logged_error_is_not_reported_as_negative(tmp_path, monkeypatch):
    """The silent-negative path: missing IKC, exit 0, valid empty VCF."""
    kwargs, _calls, processed = _run(
        tmp_path, monkeypatch, call_log_lines=[*HEALTHY[:2], *IKC_ERRORS], vcf_text=HEADER_ONLY_VCF
    )
    with pytest.raises(RuntimeError, match="Cannot open indexed k-mer count"):
        kg.run_kestrel(**kwargs)
    assert processed == [], "a failed Kestrel run reached post-processing"


def test_a_missing_reference_raises_with_the_cause_not_no_usable_vcf(tmp_path, monkeypatch):
    """The reported bug: the operator must see the reference error, not a data-shaped message."""
    kwargs, _calls, processed = _run(tmp_path, monkeypatch, call_log_lines=[REFERENCE_ERROR], vcf_text=None)
    with pytest.raises(RuntimeError) as excinfo:
        kg.run_kestrel(**kwargs)
    message = str(excinfo.value)
    assert "Error reading reference sequence(s)" in message
    assert "kestrel_kmer_20.log" in message
    assert "no usable VCF" not in message
    assert processed == []


def test_a_logged_error_stops_at_the_first_kmer_size(tmp_path, monkeypatch):
    """Every measured cause is the same at every k, so retrying would only repeat it."""
    kwargs, calls, _processed = _run(
        tmp_path, monkeypatch, call_log_lines=[REFERENCE_ERROR], vcf_text=None, kmer_sizes=(20, 25)
    )
    with pytest.raises(RuntimeError, match="k-mer size 20"):
        kg.run_kestrel(**kwargs)
    assert len(calls) == 2, "the second k-mer size ran after a logged error"


def test_the_attempt_directory_is_still_cleaned_up(tmp_path, monkeypatch):
    kwargs, _calls, _processed = _run(tmp_path, monkeypatch, call_log_lines=[REFERENCE_ERROR], vcf_text=None)
    with pytest.raises(RuntimeError):
        kg.run_kestrel(**kwargs)
    assert not attempt_directory(tmp_path, 20).exists()


def test_a_clean_log_still_produces_a_result(tmp_path, monkeypatch):
    """The check must not turn a genuine empty result into a failure."""
    kwargs, _calls, processed = _run(tmp_path, monkeypatch, call_log_lines=HEALTHY, vcf_text=HEADER_ONLY_VCF)
    kg.run_kestrel(**kwargs)
    assert len(processed) == 1


def test_no_vcf_on_every_kmer_size_names_the_logs(tmp_path, monkeypatch):
    """Without a logged error the terminal message must still point at where to look."""
    kwargs, _calls, _processed = _run(tmp_path, monkeypatch, call_log_lines=HEALTHY, vcf_text=None, kmer_sizes=(20, 25))
    with pytest.raises(RuntimeError, match="no usable VCF") as excinfo:
        kg.run_kestrel(**kwargs)
    message = str(excinfo.value)
    assert "kestrel_kmer_20.log" in message and "kestrel_kmer_25.log" in message
    assert "#338" in message
