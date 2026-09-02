"""Count and route converted FASTQs before any downstream consumer runs."""

from __future__ import annotations

import gzip
import logging
import threading
import zlib
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts import pipeline_read_routing
from vntyper.scripts.pipeline_read_routing import count_fastq_records, route_converted_fastqs

pytestmark = pytest.mark.unit


def _write_fastq(path: Path, records: int, *, lines_per_record: int = 4) -> None:
    if lines_per_record == 4:
        text = "".join(f"@read-{index}\nACGT\n+\n!!!!\n" for index in range(records))
    else:
        text = "".join(f"line-{index}\n" for index in range(records * lines_per_record))
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


def _paths(
    tmp_path: Path,
    counts: tuple[int, int, int, int],
    *,
    lines_per_record: int = 4,
) -> tuple[str, str, str, str]:
    paths = tuple(tmp_path / name for name in ("r1.fastq.gz", "r2.fastq.gz", "other.fastq.gz", "single.fastq.gz"))
    for path, count in zip(paths, counts, strict=True):
        _write_fastq(path, count, lines_per_record=lines_per_record)
    r1, r2, other, single = paths
    return str(r1), str(r2), str(other), str(single)


def test_a_zero_byte_file_takes_the_size_fast_path_without_opening_it(tmp_path: Path, monkeypatch) -> None:
    fastq = tmp_path / "empty.fastq.gz"
    fastq.touch()

    def fail_open(*args, **kwargs):
        raise AssertionError("a zero-byte FASTQ must not be decompressed")

    monkeypatch.setattr(pipeline_read_routing.gzip, "open", fail_open)

    assert count_fastq_records(fastq, lines_per_record=4) == 0


def test_plain_fastq_records_are_counted_from_complete_line_groups(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq, 2)

    assert count_fastq_records(fastq, lines_per_record=4) == 2


def test_gzipped_fastq_records_are_counted_after_decompression(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fastq.gz"
    _write_fastq(fastq, 3)

    assert count_fastq_records(fastq, lines_per_record=4) == 3


def test_an_incomplete_record_fails_closed_naming_the_file_and_line_counts(tmp_path: Path) -> None:
    fastq = tmp_path / "truncated.fastq"
    fastq.write_text("one\ntwo\nthree\nfour\nfive\n")

    with pytest.raises(ValueError, match=rf"{fastq}.*5.*4"):
        count_fastq_records(fastq, lines_per_record=4)


def test_a_missing_fastq_fails_closed_naming_the_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fastq.gz"

    with pytest.raises(ValueError, match=str(missing)):
        count_fastq_records(missing, lines_per_record=4)


def test_counting_rejects_a_non_positive_record_width_before_division(tmp_path: Path) -> None:
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq, 1)

    with pytest.raises(ValueError, match="lines_per_record"):
        count_fastq_records(fastq, lines_per_record=0)


def test_paired_outputs_are_returned_in_mate_order_and_the_layout_is_logged(tmp_path: Path, caplog) -> None:
    produced = _paths(tmp_path, (2, 2, 0, 0))

    with caplog.at_level(logging.INFO, logger=pipeline_read_routing.logger.name):
        routed = route_converted_fastqs(produced, config={})

    assert routed == (produced[0], produced[1])
    assert any("paired" in record.getMessage() for record in caplog.records)


def test_single_end_other_output_is_the_only_downstream_fastq(tmp_path: Path, caplog) -> None:
    produced = _paths(tmp_path, (0, 0, 5, 0))

    with caplog.at_level(logging.INFO, logger=pipeline_read_routing.logger.name):
        routed = route_converted_fastqs(produced, config={})

    assert routed == (produced[2],)
    assert any("single" in record.getMessage() for record in caplog.records)


def test_replacement_configs_use_the_shipped_four_line_default(tmp_path: Path) -> None:
    produced = _paths(tmp_path, (0, 0, 2, 0))

    assert route_converted_fastqs(produced, config={"utils": {}}) == (produced[2],)


def test_the_shipped_four_line_setting_counts_fastq_records(tmp_path: Path) -> None:
    produced = _paths(tmp_path, (0, 0, 3, 0), lines_per_record=4)

    assert route_converted_fastqs(produced, config={"utils": {"fastq_validation_lines": 4}}) == (produced[2],)


def test_fastq_record_structure_cannot_be_configured_away_from_four_lines(tmp_path: Path) -> None:
    """FASTQ has exactly four lines per record; another value would miscount reads."""
    produced = _paths(tmp_path, (0, 0, 1, 0), lines_per_record=2)

    with pytest.raises(ValueError, match="must be 4"):
        route_converted_fastqs(produced, config={"utils": {"fastq_validation_lines": 2}})


@pytest.mark.parametrize("bad_value", [0, -1, True, "4"])
def test_invalid_lines_per_record_configuration_fails_closed(tmp_path: Path, bad_value: object) -> None:
    produced = _paths(tmp_path, (0, 0, 1, 0))

    with pytest.raises(ValueError, match="fastq_validation_lines"):
        route_converted_fastqs(produced, config={"utils": {"fastq_validation_lines": bad_value}})


@pytest.mark.parametrize(
    ("counts", "selected_indices", "expected_record"),
    [
        (
            (3, 3, 0, 1),
            (0, 1, 3),
            'READ_SET_ROUTING {"counts":{"other":0,"r1":3,"r2":3,"single":1},"layout":"mixed",'
            '"selected":["r1.fastq.gz","r2.fastq.gz","single.fastq.gz"]}',
        ),
        (
            (3, 3, 2, 1),
            (0, 1, 2, 3),
            'READ_SET_ROUTING {"counts":{"other":2,"r1":3,"r2":3,"single":1},"layout":"mixed",'
            '"selected":["r1.fastq.gz","r2.fastq.gz","other.fastq.gz","single.fastq.gz"]}',
        ),
        (
            (0, 0, 2, 1),
            (2, 3),
            'READ_SET_ROUTING {"counts":{"other":2,"r1":0,"r2":0,"single":1},"layout":"mixed",'
            '"selected":["other.fastq.gz","single.fastq.gz"]}',
        ),
    ],
)
def test_counted_mixed_layout_routes_every_file_and_emits_one_canonical_record(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    counts: tuple[int, int, int, int],
    selected_indices: tuple[int, ...],
    expected_record: str,
) -> None:
    produced = _paths(tmp_path, counts)

    with caplog.at_level(logging.INFO, logger=pipeline_read_routing.logger.name):
        routed = route_converted_fastqs(produced, config={})

    assert routed == tuple(produced[index] for index in selected_indices)
    records = [record.getMessage() for record in caplog.records if record.getMessage().startswith("READ_SET_ROUTING ")]
    assert records == [expected_record]


@pytest.mark.parametrize("counts", [(3, 2, 1, 0), (1, 0, 0, 0), (0, 1, 0, 0)])
def test_invalid_layout_failure_names_every_file_and_its_record_count(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    counts: tuple[int, int, int, int],
) -> None:
    produced = _paths(tmp_path, counts)

    with (
        caplog.at_level(logging.ERROR, logger=pipeline_read_routing.logger.name),
        pytest.raises(ValueError) as exc_info,
    ):
        route_converted_fastqs(produced, config={})

    message = str(exc_info.value)
    assert "invalid" in message
    assert "mate outputs are inconsistent" in message
    assert "cannot be consumed without dropping reads" not in message
    for path, count in zip(produced, counts, strict=True):
        assert path in message
        assert f"{count} records" in message
    assert [record.getMessage() for record in caplog.records if record.levelno == logging.ERROR] == [message]


def test_an_undecompressible_fastq_fails_closed_naming_the_file(tmp_path: Path) -> None:
    corrupt = tmp_path / "corrupt.fastq.gz"
    corrupt.write_text("not gzip", encoding="utf-8")

    with pytest.raises(ValueError, match=str(corrupt)):
        count_fastq_records(corrupt, lines_per_record=4)


def test_a_gzip_with_a_corrupt_deflate_body_fails_closed_with_path_cause_and_log(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    corrupt = tmp_path / "body-corrupt.fastq.gz"
    compressed = bytearray(gzip.compress(b"@read-0\nACGT\n+\n!!!!\n", mtime=0))
    compressed[10] = (compressed[10] & ~0x07) | 0x07
    corrupt.write_bytes(compressed)

    with (
        caplog.at_level(logging.ERROR, logger=pipeline_read_routing.logger.name),
        pytest.raises(ValueError) as exc_info,
    ):
        count_fastq_records(corrupt, lines_per_record=4)

    message = str(exc_info.value)
    assert message.startswith(f"Could not count FASTQ records in {corrupt}:")
    assert "invalid block type" in message
    assert isinstance(exc_info.value.__cause__, zlib.error)
    assert str(exc_info.value.__cause__) in message
    assert [record.getMessage() for record in caplog.records if record.levelno == logging.ERROR] == [message]


def test_empty_layout_failure_names_every_file_with_zero_records(tmp_path: Path) -> None:
    produced = _paths(tmp_path, (0, 0, 0, 0))

    with pytest.raises(ValueError) as exc_info:
        route_converted_fastqs(produced, config={})

    message = str(exc_info.value)
    assert "empty" in message
    assert "cannot be consumed without dropping reads" in message
    assert "mate outputs are inconsistent" not in message
    for path in produced:
        assert path in message
    assert message.count("0 records") == 4


def test_a_stranded_path_from_the_router_is_never_ignored(tmp_path: Path, monkeypatch) -> None:
    produced = _paths(tmp_path, (1, 1, 0, 0))
    monkeypatch.setattr(
        pipeline_read_routing,
        "route_fastqs",
        lambda layout, paths, counts: ((paths["r1"], paths["r2"]), (paths["other"],)),
    )

    with pytest.raises(ValueError) as exc_info:
        route_converted_fastqs(produced, config={})

    message = str(exc_info.value)
    assert produced[2] in message
    assert "0 records" in message
    assert "cannot be consumed without dropping reads" in message
    assert "mate outputs are inconsistent" not in message


# ---------------------------------------------------------------------------
# Binary counting on a thread pool (#262)
# ---------------------------------------------------------------------------


def test_counting_never_decodes_the_file_into_text(tmp_path: Path, monkeypatch) -> None:
    """The red test for #262: decoding 3.85M lines into str objects cost 1044ms.

    Classifying the read layout and logging one line needs the number of newlines and
    nothing else. This fails against the text iterator, which opens in "rt".
    """
    fastq = tmp_path / "reads.fastq.gz"
    _write_fastq(fastq, 1)
    real_open = gzip.open

    def binary_only(path, mode="rb", *args, **kwargs):
        assert "t" not in mode, f"count_fastq_records must open in binary mode, got {mode!r}"
        return real_open(path, mode, *args, **kwargs)

    monkeypatch.setattr(pipeline_read_routing.gzip, "open", binary_only)

    assert count_fastq_records(fastq, lines_per_record=4) == 1


def test_a_plain_fastq_is_also_opened_in_binary(tmp_path: Path, monkeypatch) -> None:
    """The uncompressed branch must not keep decoding text either."""
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq, 2)
    real_open = Path.open

    def binary_only(self, mode="r", *args, **kwargs):
        assert "b" in mode, f"count_fastq_records must open in binary mode, got {mode!r}"
        return real_open(self, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "open", binary_only)

    assert count_fastq_records(fastq, lines_per_record=4) == 2


@pytest.mark.parametrize("records", [1, 2, 17])
def test_binary_counting_returns_the_same_integers_as_text_counting(tmp_path: Path, records: int) -> None:
    """The contract is the integer, and it must not move for any input the text path took."""
    plain = tmp_path / f"plain-{records}.fastq"
    gzipped = tmp_path / f"gz-{records}.fastq.gz"
    _write_fastq(plain, records)
    _write_fastq(gzipped, records)

    assert count_fastq_records(plain, lines_per_record=4) == records
    assert count_fastq_records(gzipped, lines_per_record=4) == records


def test_a_final_line_without_a_terminator_is_still_counted(tmp_path: Path) -> None:
    """Counting b"\\n" would lose it, and the text iterator did not.

    A FASTQ whose last quality line has no trailing newline is a complete record, and
    dropping that line would make the line count indivisible and fail a valid file.
    """
    ragged = tmp_path / "ragged.fastq"
    ragged.write_text("@r1\nACGT\n+\n!!!!", encoding="utf-8")

    assert count_fastq_records(ragged, lines_per_record=4) == 1


def test_a_non_utf8_fastq_is_still_rejected(tmp_path: Path) -> None:
    """The text reader raised UnicodeError, which this function turned into ValueError.

    Binary counting cannot get that for free, so validation is restored explicitly
    rather than silently dropped -- this is a deliberate contract that is kept, not a
    behaviour that survived by accident.
    """
    bad = tmp_path / "latin1.fastq"
    bad.write_bytes(b"@r1\n\xff\xfe\n+\n!!\n")

    with pytest.raises(ValueError, match=str(bad)):
        count_fastq_records(bad, lines_per_record=4)


def test_a_multibyte_character_spanning_a_chunk_boundary_is_accepted(tmp_path: Path, monkeypatch) -> None:
    """Strict-decoding each chunk independently rejects valid files.

    A two-byte character straddling the read boundary raises "unexpected end of data"
    when each chunk is decoded on its own. An incremental decoder carries the partial
    character across. Every FASTQ fixture in this repository is ASCII, so without this
    test the regression is invisible.
    """
    fastq = tmp_path / "multibyte.fastq"
    # "é" is 0xC3 0xA9; put the read boundary between its two bytes.
    fastq.write_bytes("@ré\nACGT\n+\n!!!!\n".encode())
    monkeypatch.setattr(pipeline_read_routing, "_COUNT_CHUNK_BYTES", 4)

    assert count_fastq_records(fastq, lines_per_record=4) == 1


def test_a_truncated_multibyte_character_at_end_of_file_is_rejected(tmp_path: Path) -> None:
    """final=True is what catches a file whose last character is half-written."""
    bad = tmp_path / "truncated-char.fastq"
    bad.write_bytes(b"@r1\nACGT\n+\n!!!\xc3")

    with pytest.raises(ValueError, match=str(bad)):
        count_fastq_records(bad, lines_per_record=4)


def test_the_four_counts_run_concurrently(tmp_path: Path) -> None:
    """Serial counting spent 1044ms decoding; zlib releases the GIL, so they overlap.

    Asserting on overlap rather than on elapsed time: a timing assertion is flaky on a
    loaded machine. Counting distinct worker threads is racy too -- with two-record
    fixtures the first worker can finish before the second submit, so the executor
    reuses it and only one thread identity is ever seen, which is what happened under
    coverage tracing. A barrier inside the count is deterministic: it can
    only be passed while every count is in flight at once, and a serial implementation
    trips its timeout instead of passing by luck.
    """
    produced = _paths(tmp_path, (2, 2, 0, 0))
    threads: set[int] = set()
    # One party per produced file: every count must be in flight at once, which is
    # exactly what a pool sized to the number of files provides.
    overlap = threading.Barrier(len(produced), timeout=10)
    real_count = pipeline_read_routing.count_fastq_records

    def record_thread(path, *, lines_per_record):
        threads.add(threading.get_ident())
        overlap.wait()
        return real_count(path, lines_per_record=lines_per_record)

    with mock.patch.object(pipeline_read_routing, "count_fastq_records", record_thread):
        assert route_converted_fastqs(produced, config={}) == (produced[0], produced[1])

    assert not overlap.broken
    assert len(threads) > 1
    assert threading.get_ident() not in threads


def test_the_first_corrupt_file_reported_is_r1_not_whichever_thread_finished(tmp_path: Path) -> None:
    """Parallel counting must not make the error message depend on scheduling.

    R1 and R2 are both invalid. The text path resolved them in R1/R2/other/single order
    through a dict comprehension, so R1 was always the reported one; resolving futures
    as they complete would make the message a race.
    """
    produced = _paths(tmp_path, (1, 1, 0, 0))
    for index in (0, 1):
        Path(produced[index]).write_bytes(b"not gzip at all")

    with pytest.raises(ValueError, match="r1.fastq.gz"):
        route_converted_fastqs(produced, config={})


def test_a_counting_failure_still_raises_through_the_thread_pool(tmp_path: Path, monkeypatch) -> None:
    """future.result() re-raises, so the per-file ValueError contract is preserved."""
    produced = _paths(tmp_path, (2, 2, 0, 0))

    def explode(*args, **kwargs):
        raise ValueError("boom")

    monkeypatch.setattr(pipeline_read_routing, "count_fastq_records", explode)

    with pytest.raises(ValueError, match="boom"):
        route_converted_fastqs(produced, config={})
