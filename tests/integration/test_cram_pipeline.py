"""Real CRAM pipeline cases using the shared transport-independent oracle."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from tests.integration.test_pipeline_integration import _fresh_output_dir, _run_local_pipeline
from tests.parametrization import get_cram_test_cases, get_cram_test_ids
from tests.support.orchestration import run_cram_test_case

pytestmark = pytest.mark.integration

REPOSITORY = Path(__file__).resolve().parents[2]
REFERENCE_COMPRESSED_CRAM = REPOSITORY / "tests/data/cram/reference-compressed/example_b178_hg19_subset.cram"


def _reheader_reference_compressed_cram(source: Path, destination: Path, reference: Path | None) -> None:
    """Replace every local reference URI without decoding the CRAM records."""
    header = subprocess.check_output(["samtools", "view", "-H", "--no-PG", str(source)], text=True)
    rewritten: list[str] = []
    sequence_lines = 0
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            fields = [field for field in line.split("\t") if not field.startswith("UR:")]
            if reference is not None:
                fields.append(f"UR:{reference}")
            line = "\t".join(fields)
            sequence_lines += 1
        rewritten.append(line)
    assert sequence_lines > 0
    header_path = destination.with_suffix(".header.sam")
    header_path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    with destination.open("wb") as output:
        subprocess.run(
            ["samtools", "reheader", "-P", str(header_path), str(source)],
            stdout=output,
            check=True,
        )


def _run_reference_failure(tmp_path: Path, reference: Path | None) -> subprocess.CompletedProcess[str]:
    """Run the real reference-compressed CRAM with every valid fallback disabled."""
    input_root = tmp_path / "operator-input"
    input_root.mkdir()
    input_cram = input_root / "reference-failure.cram"
    _reheader_reference_compressed_cram(REFERENCE_COMPRESSED_CRAM, input_cram, reference)

    config_payload = json.loads((REPOSITORY / "vntyper/config.json").read_text(encoding="utf-8"))
    config_payload["reference_data"]["cram_reference_hg19"] = None
    config_payload["reference_data"]["bwa_reference_hg19"] = None
    config_payload["cram"]["local_ref_path"] = str(tmp_path / "empty-reference-cache" / "%2s/%2s/%s")
    config = tmp_path / "config.json"
    config.write_text(json.dumps(config_payload), encoding="utf-8")

    return subprocess.run(
        [
            sys.executable,
            "-m",
            "vntyper.cli",
            "--config-path",
            str(config),
            "-l",
            "DEBUG",
            "pipeline",
            "--cram",
            str(input_cram),
            "--threads",
            "2",
            "--reference-assembly",
            "hg19",
            "--output-dir",
            str(tmp_path / "output"),
            "--fast-mode",
            "--keep-intermediates",
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def _assert_causal_decode_failure(
    result: subprocess.CompletedProcess[str],
    output: Path,
    cause: str,
) -> None:
    """Require a real htslib reference decode failure before downstream stages."""
    diagnostic = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == 1, diagnostic
    assert "cram_decode_slice" in diagnostic
    assert cause in diagnostic
    assert not (output / "coverage").exists()
    assert not (output / "kestrel").exists()


@pytest.mark.parametrize("cram_case", get_cram_test_cases(), ids=get_cram_test_ids())
def test_cram_input(tmp_path: Path, ensure_test_data: None, cram_case: dict) -> None:
    """Run one declared CRAM through the same strict oracle as every other transport.

    Args:
        tmp_path: Pytest-managed output root.
        ensure_test_data: Session fixture proving the registered data are available.
        cram_case: One entry of ``integration_tests.cram_tests``.
    """
    output_dir = _fresh_output_dir(tmp_path, cram_case["test_name"])
    run_cram_test_case(cram_case, _run_local_pipeline, output_dir)


def test_reference_compressed_cram_fails_causally_when_reference_is_missing(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """Prove the real fixture cannot decode from an empty pinned M5 cache."""
    result = _run_reference_failure(tmp_path, None)

    _assert_causal_decode_failure(result, tmp_path / "output", "Unable to fetch reference")


def test_reference_compressed_cram_fails_causally_with_wrong_reference(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """Prove an existing but incorrect FASTA cannot decode the real fixture."""
    wrong_reference = tmp_path / "wrong.fa"
    wrong_reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="ascii")
    subprocess.run(["samtools", "faidx", str(wrong_reference)], check=True)

    result = _run_reference_failure(tmp_path, wrong_reference)

    _assert_causal_decode_failure(result, tmp_path / "output", "MD5 checksum reference mismatch")
