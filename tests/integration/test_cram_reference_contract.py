"""A-209 reference-dependent and reference-free CRAM pipeline contracts."""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from scripts.make_cram_fixtures import build_reference_dependent_fixture

logger = logging.getLogger(__name__)

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_DATA_CONFIG = json.loads((REPO_ROOT / "tests" / "test_data_config.json").read_text(encoding="utf-8"))
PURPOSE_FIXTURES = TEST_DATA_CONFIG["purpose_fixtures"]["cram_reference_contract"]


def _require_purpose_fixture(name: str) -> Path:
    """Resolve one registered fixture or name the ordinary derivation command."""
    path = REPO_ROOT / PURPOSE_FIXTURES[name]
    assert path.is_file(), f"Missing purpose fixture {path}; run `make cram-fixtures`."
    return path


def _local_config(tmp_path: Path) -> Path:
    """Pin reference lookup and coverage to the purpose fixture's compact contig."""
    payload = json.loads((REPO_ROOT / "vntyper" / "config.json").read_text(encoding="utf-8"))
    payload["cram"]["local_ref_path"] = str(tmp_path / "local-ref" / "%2s" / "%2s" / "%s")
    payload["bam_processing"]["assemblies"]["GRCh37"]["vntr_region_coords"] = "1-10000"
    config = tmp_path / "config.json"
    config.write_text(json.dumps(payload), encoding="utf-8")
    return config


def _no_ref_config(tmp_path: Path) -> Path:
    """Remove every explicit hg19 reference while keeping terminal resolution last."""
    config = _local_config(tmp_path)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["reference_data"]["cram_reference_hg19"] = None
    payload["reference_data"]["bwa_reference_hg19"] = None
    payload["cram"]["reference_candidate_order"] = [
        "cli",
        "config_cram_reference",
        "config_bwa_reference",
        "htslib_resolved",
    ]
    config.write_text(json.dumps(payload), encoding="utf-8")
    return config


def _run_cram(
    cram: Path,
    output_dir: Path,
    config: Path,
    *,
    reference: Path | None = None,
    fast_mode: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run the purpose CRAM over its own compact target interval."""
    command = [
        sys.executable,
        "-m",
        "vntyper.cli",
        "--config-path",
        str(config),
        "-l",
        "DEBUG",
        "pipeline",
        "--cram",
        str(cram),
        "--reference-assembly",
        "hg19",
        "--custom-regions",
        "chr1:1-10000",
        "--threads",
        "2",
        "--keep-intermediates",
        "--output-dir",
        str(output_dir),
    ]
    if fast_mode:
        command.append("--fast-mode")
    if reference is not None:
        command.extend(["--reference-fasta", str(reference)])
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    logger.info("CRAM contract stdout:\n%s", result.stdout)
    logger.info("CRAM contract stderr:\n%s", result.stderr)
    return result


def _assert_no_reference_flags_in_cram_commands(result: subprocess.CompletedProcess[str]) -> None:
    """Require every logged CRAM consumer to use the terminal no-reference form."""
    markers = (
        "Running captured command: samtools view ",
        "Executing region slicing with command: samtools view ",
        "Executing filtering with command: set -o pipefail; samtools view ",
        "Calculating VNTR coverage with command: samtools depth ",
    )
    for marker in markers:
        commands = [line for line in result.stderr.splitlines() if marker in line]
        assert commands, f"No logged CRAM command matched {marker!r}"
        assert all(" -T " not in command and " --reference " not in command for command in commands), commands


def test_a209_1_missing_reference_names_the_digest_and_candidates_before_stages(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """The actual header ``UR:`` target is renamed before the failing decode."""
    _require_purpose_fixture("reference_dependent_cram")
    _require_purpose_fixture("reference_fasta")
    fixture = build_reference_dependent_fixture(tmp_path / "purpose")
    with pysam.AlignmentFile(str(fixture.cram), "rc", reference_filename=str(fixture.reference)) as alignment:
        sequence = alignment.header.to_dict()["SQ"][0]
    ur_target = Path(sequence["UR"])
    expected_m5 = hashlib.md5(("A" * 10_000).encode()).hexdigest()
    assert ur_target.samefile(fixture.reference)
    assert sequence["M5"] == expected_m5

    hidden_target = ur_target.with_name(f"{ur_target.name}.a209-missing")
    ur_target.rename(hidden_target)
    output = tmp_path / "missing-output"
    try:
        config = _local_config(tmp_path)
        result = _run_cram(fixture.cram, output, config)
    finally:
        hidden_target.rename(ur_target)

    diagnostic = f"{result.stdout}\n{result.stderr}"
    assert result.returncode != 0
    assert "contig=chr1" in diagnostic
    assert f"M5={expected_m5}" in diagnostic
    for source in (
        "source=cli",
        "source=config_cram_reference",
        "source=config_bwa_reference",
        "source=htslib-resolved",
    ):
        assert source in diagnostic
    assert not (output / "kestrel").exists()
    assert not (output / "coverage").exists()
    assert not (output / "pipeline_summary.json").exists()


def test_a209_2_explicit_reference_completes_the_reference_dependent_cram(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    cram = _require_purpose_fixture("reference_dependent_cram")
    reference = _require_purpose_fixture("reference_fasta")
    output = tmp_path / "explicit-reference-output"

    result = _run_cram(cram, output, _local_config(tmp_path), reference=reference)

    assert result.returncode == 0
    assert (output / "pipeline_summary.json").is_file()
    assert (output / "kestrel" / "kestrel_result.tsv").is_file()


def test_a209_3_no_ref_cram_completes_without_an_explicit_reference(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    cram = _require_purpose_fixture("no_ref_cram")
    output = tmp_path / "no-ref-output"

    result = _run_cram(cram, output, _no_ref_config(tmp_path), fast_mode=False)

    assert result.returncode == 0
    assert "Resolved CRAM reference through htslib-resolved" in result.stderr
    _assert_no_reference_flags_in_cram_commands(result)
    assert (output / "pipeline_summary.json").is_file()
    assert (output / "kestrel" / "kestrel_result.tsv").is_file()
