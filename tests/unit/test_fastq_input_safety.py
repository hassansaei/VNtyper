"""FASTQ quality-control destinations must never overwrite operator inputs."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from vntyper.scripts import fastq_bam_processing

pytestmark = pytest.mark.unit

CONFIG = {
    "tools": {"fastp": "fastp"},
    "bam_processing": {
        "compression_level": 6,
        "disable_adapter_trimming": True,
        "deduplication": True,
        "dup_calc_accuracy": 3,
        "length_required": 50,
        "qualified_quality_phred": 20,
    },
}

FASTP_DESTINATION_NAMES = (
    "output_R1.fastq.gz",
    "output_R2.fastq.gz",
    "output.html",
    "output.json",
    "output_fastp.log",
)


def _assert_fastq_alias_rejected_before_work(
    tmp_path: Path,
    *,
    input_slot: str,
    fastq_path: Path,
) -> None:
    mate = tmp_path / "input" / "other.fastq.gz"
    mate.parent.mkdir(exist_ok=True)
    mate.write_bytes(b"other-fastq")
    fastq_1 = fastq_path if input_slot == "fastq1" else mate
    fastq_2 = fastq_path if input_slot == "fastq2" else mate
    command_runner = MagicMock(return_value=True)

    with (
        patch.object(fastq_bam_processing, "run_command", command_runner),
        pytest.raises(ValueError, match="aliases protected input"),
    ):
        fastq_bam_processing.process_fastq(fastq_1, fastq_2, 4, str(tmp_path / "output"), "output", CONFIG)

    command_runner.assert_not_called()


@pytest.mark.parametrize("input_slot", ["fastq1", "fastq2"])
@pytest.mark.parametrize("destination_name", FASTP_DESTINATION_NAMES)
def test_direct_fastq_alias_of_any_fastp_destination_is_rejected_before_work(
    tmp_path: Path,
    input_slot: str,
    destination_name: str,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    fastq_path = output / destination_name
    fastq_path.write_bytes(b"operator-fastq")

    _assert_fastq_alias_rejected_before_work(tmp_path, input_slot=input_slot, fastq_path=fastq_path)

    assert fastq_path.read_bytes() == b"operator-fastq"


@pytest.mark.parametrize("input_slot", ["fastq1", "fastq2"])
@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_fastq_filesystem_alias_of_a_fastp_destination_is_rejected_before_work(
    tmp_path: Path,
    input_slot: str,
    alias_kind: str,
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    destination = output / "output.json"
    destination.write_bytes(b"operator-fastq")
    fastq_path = tmp_path / "input" / "protected.fastq.gz"
    fastq_path.parent.mkdir()
    if alias_kind == "symlink":
        fastq_path.symlink_to(destination)
    else:
        fastq_path.hardlink_to(destination)

    _assert_fastq_alias_rejected_before_work(tmp_path, input_slot=input_slot, fastq_path=fastq_path)

    assert fastq_path.read_bytes() == b"operator-fastq"


def test_single_link_regular_fastp_artifacts_remain_replaceable(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    for destination_name in FASTP_DESTINATION_NAMES:
        (output / destination_name).write_bytes(b"stale-run-artifact")
    input_root = tmp_path / "input"
    input_root.mkdir()
    fastq_1 = input_root / "r1.fastq.gz"
    fastq_2 = input_root / "r2.fastq.gz"
    fastq_1.write_bytes(b"r1")
    fastq_2.write_bytes(b"r2")
    command_runner = MagicMock(return_value=True)

    with patch.object(fastq_bam_processing, "run_command", command_runner):
        fastq_bam_processing.process_fastq(fastq_1, fastq_2, 4, output, "output", CONFIG)

    command_runner.assert_called_once()
