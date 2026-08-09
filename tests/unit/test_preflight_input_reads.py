"""Nonblocking and bounded local-file reads before alignment probes."""

from __future__ import annotations

import json
import logging
import multiprocessing
import os
import time
from collections.abc import Callable
from multiprocessing.connection import Connection
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from vntyper.scripts import alignment_preflight
from vntyper.scripts.pipeline_alignment import build_alignment_preflight_kwargs, prepare_alignment_target

pytestmark = pytest.mark.unit


def _reference_reason_worker(path: str, sender: Connection) -> None:
    try:
        sender.send(("returned", alignment_preflight._reference_unavailable_reason(path)))
    except BaseException as error:
        sender.send(("raised", type(error).__name__, str(error)))
    finally:
        sender.close()


def _bed_read_worker(path: str, sender: Connection) -> None:
    try:
        build_alignment_preflight_kwargs(
            in_path="/input/sample.bam",
            output_dir="/output",
            output_name="input",
            file_format="bam",
            config={},
            threads=1,
            bed_file=Path(path),
            reference_assembly="hg19",
            fast_mode=False,
        )
        sender.send(("returned",))
    except BaseException as error:
        sender.send(("raised", type(error).__name__, str(error)))
    finally:
        sender.close()


def _run_with_outer_deadline(worker: Callable[[str, Connection], None], path: Path) -> tuple[Any, ...]:
    context = multiprocessing.get_context("spawn")
    receiver, sender = context.Pipe(duplex=False)
    process = context.Process(target=worker, args=(str(path), sender))
    started = time.monotonic()
    process.start()
    sender.close()
    process.join(2)
    elapsed = time.monotonic() - started
    blocked = process.is_alive()
    if blocked:
        process.terminate()
        process.join(1)
    try:
        assert not blocked, "pre-probe local-file inspection waited for a FIFO writer"
        assert not process.is_alive(), "pre-probe local-file inspection did not terminate"
        assert elapsed < 3, f"pre-probe local-file inspection blocked for {elapsed:.2f} s"
        assert receiver.poll(), "pre-probe local-file inspection returned no result"
        return receiver.recv()
    finally:
        receiver.close()


def test_reference_fifo_is_rejected_without_waiting_for_a_writer(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    os.mkfifo(reference)

    result = _run_with_outer_deadline(_reference_reason_worker, reference)

    assert result == ("returned", "reference FASTA is not a regular file")


def test_supplied_bed_fifo_is_rejected_without_waiting_for_a_writer(tmp_path: Path) -> None:
    bed_file = tmp_path / "target.bed"
    os.mkfifo(bed_file)

    result = _run_with_outer_deadline(_bed_read_worker, bed_file)

    assert result == ("raised", "ValueError", f"alignment target BED is not a regular file: {bed_file}")


def test_target_preparation_rejects_a_supplied_fifo_before_returning_it(tmp_path: Path) -> None:
    bed_file = tmp_path / "target.bed"
    os.mkfifo(bed_file)

    with pytest.raises(ValueError) as raised:
        prepare_alignment_target(
            input_type="BAM",
            bam="/input/sample.bam",
            cram=None,
            output_dir=tmp_path / "output",
            reference_assembly="hg19",
            config={},
            bed_file=bed_file,
            custom_regions=None,
        )

    assert str(raised.value) == f"alignment target BED is not a regular file: {bed_file}"


def test_direct_bam_preflight_rejects_a_bed_fifo_before_view_or_probe_work(tmp_path: Path) -> None:
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"BAM")
    bed_file = tmp_path / "target.bed"
    os.mkfifo(bed_file)

    with (
        patch("vntyper.scripts.alignment_preflight.build_alignment_view") as build_view,
        patch("vntyper.scripts.alignment_preflight.capture_command") as capture,
        pytest.raises(ValueError) as raised,
    ):
        alignment_preflight.run_preflight(
            str(alignment),
            str(tmp_path / "output"),
            "sample",
            "bam",
            {},
            1,
            bed_file=bed_file,
        )

    assert str(raised.value) == f"alignment target BED is not a regular file: {bed_file}"
    build_view.assert_not_called()
    capture.assert_not_called()


def test_oversized_bed_is_rejected_before_preflight_work(tmp_path: Path) -> None:
    bed_file = tmp_path / "target.bed"
    bed_file.write_text("#" * 65, encoding="utf-8")

    with pytest.raises(ValueError) as raised:
        build_alignment_preflight_kwargs(
            in_path="/input/sample.bam",
            output_dir=tmp_path / "output",
            output_name="input",
            file_format="bam",
            config={"utils": {"preflight_text_max_bytes": 64}},
            threads=1,
            bed_file=bed_file,
            reference_assembly="hg19",
            fast_mode=False,
        )

    assert str(raised.value) == f"alignment target BED exceeds configured 64-byte preflight limit: {bed_file}"


def test_oversized_fai_keeps_the_successful_probe_but_skips_coverage_inference(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(f"{reference}.fai").write_text("chr1\t4\t6\t4\t5\n" + "x" * 65, encoding="utf-8")

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", return_value=(True, "decoded")),
        caplog.at_level("WARNING"),
    ):
        resolved, source, uncovered, _binding = alignment_preflight.resolve_reference(
            "/run/view.cram",
            (("cli", str(reference)),),
            "chr1:1-2",
            None,
            {"utils": {"preflight_text_max_bytes": 64}},
            1,
            str(tmp_path),
            "sample",
            ("chr1", "chr2"),
            "abc",
        )

    assert resolved is not None
    assert Path(resolved).read_bytes() == reference.read_bytes()
    assert (source, uncovered) == ("cli", ())
    assert "coverage unavailable" in caplog.text.lower()
    assert _binding is not None
    _binding.close()
    assert not any(record.levelno >= logging.ERROR for record in caplog.records)


@pytest.mark.parametrize("invalid", [0, -1, True, "64", None, 1.5])
def test_invalid_preflight_text_limit_fails_before_reading_the_bed(tmp_path: Path, invalid: object) -> None:
    bed_file = tmp_path / "target.bed"
    bed_file.write_text("chr1\t1\t2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="utils.preflight_text_max_bytes must be a positive integer"):
        build_alignment_preflight_kwargs(
            in_path="/input/sample.bam",
            output_dir=tmp_path / "output",
            output_name="input",
            file_format="bam",
            config={"utils": {"preflight_text_max_bytes": invalid}},
            threads=1,
            bed_file=bed_file,
            reference_assembly="hg19",
            fast_mode=False,
        )


def test_shipped_config_declares_the_default_preflight_text_limit() -> None:
    config = json.loads((Path(__file__).parents[2] / "vntyper" / "config.json").read_text(encoding="utf-8"))

    assert config["utils"]["preflight_text_max_bytes"] == 1048576
