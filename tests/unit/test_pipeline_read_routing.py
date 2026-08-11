"""Count and route converted FASTQs before any downstream consumer runs."""

from __future__ import annotations

import gzip
import logging
import zlib
from pathlib import Path

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
