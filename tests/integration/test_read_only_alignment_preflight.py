"""Real read-only alignment-view proofs for milestone-4 index handling."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from vntyper.scripts.alignment_contract import index_candidate_names
from vntyper.scripts.alignment_preflight import run_preflight
from vntyper.scripts.command_builders import build_cram_unmapped_filter_command, build_samtools_slice_command

pytestmark = pytest.mark.integration


def _tree_digest(root: Path) -> dict[str, str]:
    """Return content digests for every regular file beneath ``root``."""
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


@pytest.mark.parametrize("file_format", ["bam", "cram"])
def test_unindexed_read_only_alignment_uses_a_run_local_index_for_real_slice(
    tmp_path: Path,
    ensure_test_data: None,
    monkeypatch: pytest.MonkeyPatch,
    file_format: str,
) -> None:
    """BAM and CRAM inputs stay byte-identical while htslib finds the view index."""
    repository = Path(__file__).parents[2]
    input_root = tmp_path / "patient-input"
    output_root = tmp_path / "run"
    input_root.mkdir()
    output_root.mkdir()
    if file_format == "bam":
        source = repository / "tests/data/example_b178_hg19_subset.bam"
        alignment = input_root / "sample.bam"
        region = "chr1:155158000-155163000"
        reference: Path | None = None
    else:
        fixture_root = repository / "tests/data/cram/reference-dependent"
        source = fixture_root / "reference-dependent.cram"
        alignment = input_root / "sample.cram"
        region = "chr1:1-10000"
        reference = fixture_root / "reference.fa"
    shutil.copyfile(source, alignment)
    alignment.chmod(0o444)
    input_root.chmod(0o555)
    before = _tree_digest(input_root)
    monkeypatch.setenv("REF_PATH", str(tmp_path / "no-reference-cache" / "%2s/%2s/%s"))
    config = {
        "tools": {"samtools": "samtools"},
        "bam": {"unmapped_scan": "auto"},
        "cram": {
            "allow_ambient_reference_resolution": False,
            "local_ref_path": str(tmp_path / "no-reference-cache" / "%2s/%2s/%s"),
            "reference_probe_timeout_seconds": 10,
            "unmapped_scan": "auto",
            "reference_candidate_order": ["cli", "htslib_resolved"],
        },
    }

    try:
        plan = run_preflight(
            str(alignment),
            str(output_root),
            "input",
            file_format,
            config,
            2,
            region=region,
            reference_fasta=str(reference) if reference is not None else None,
            header_contigs=("chr1",),
        )
        slice_path = output_root / "slice.bam"
        command = build_samtools_slice_command(
            samtools_path="samtools",
            in_bam=plan.view_path,
            output_bam=slice_path,
            region=region,
            reference_path=plan.reference_path,
            threads=2,
        )
        completed = subprocess.run(
            command,
            shell=True,
            executable="/bin/bash",
            capture_output=True,
            text=True,
            check=False,
        )
    finally:
        input_root.chmod(0o755)
        alignment.chmod(0o644)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert int(subprocess.check_output(["samtools", "view", "-c", str(slice_path)], text=True)) > 0
    assert Path(plan.index_path).is_file()
    assert not Path(plan.index_path).is_symlink()
    assert all(not os.path.lexists(candidate) for candidate in index_candidate_names(str(alignment), file_format))
    assert _tree_digest(input_root) == before


def test_a_wrong_sample_source_bai_is_ignored_before_real_bam_slice(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """A valid wrong BAI cannot authorize the empty slice reproduced in round 3."""
    repository = Path(__file__).parents[2]
    input_root = tmp_path / "patient-input"
    output_root = tmp_path / "run"
    input_root.mkdir()
    output_root.mkdir()
    alignment = input_root / "sample.bam"
    supplied_index = input_root / "sample.bam.bai"
    shutil.copyfile(repository / "tests/data/example_b178_hg19_subset.bam", alignment)
    shutil.copyfile(repository / "tests/data/example_40cf_hg38_subset.bam.bai", supplied_index)
    before = _tree_digest(input_root)

    plan = run_preflight(
        str(alignment),
        str(output_root),
        "input",
        "bam",
        {"tools": {"samtools": "samtools"}, "bam": {"unmapped_scan": "auto"}},
        2,
        region="chr1:155158000-155163000",
    )
    slice_path = output_root / "slice.bam"
    command = build_samtools_slice_command(
        samtools_path="samtools",
        in_bam=plan.view_path,
        output_bam=slice_path,
        region="chr1:155158000-155163000",
        threads=2,
    )
    completed = subprocess.run(command, shell=True, executable="/bin/bash", check=False)

    assert completed.returncode == 0
    assert int(subprocess.check_output(["samtools", "view", "-c", str(slice_path)], text=True)) == 29736
    assert not Path(plan.index_path).is_symlink()
    assert supplied_index.samefile(input_root / "sample.bam.bai")
    assert _tree_digest(input_root) == before


@pytest.mark.parametrize("file_format", ["bam", "cram"])
def test_single_end_and_placed_unmapped_reads_survive_the_complete_scan(
    tmp_path: Path,
    ensure_test_data: None,
    monkeypatch: pytest.MonkeyPatch,
    file_format: str,
) -> None:
    """Flag-4 extraction retains records that flag 12 and the BAI tail both lose."""
    repository = Path(__file__).parents[2]
    source_bam = repository / "tests/data/derived/single_end/example_b178_hg19_single_end.bam"
    reference = repository / "reference/alignment/chr1.hg19.fa"
    alignment = tmp_path / f"single-end.{file_format}"
    if file_format == "bam":
        shutil.copyfile(source_bam, alignment)
    else:
        subprocess.run(
            ["samtools", "view", "-C", "-T", str(reference), "-o", str(alignment), str(source_bam)],
            check=True,
        )
    monkeypatch.setenv("REF_PATH", str(tmp_path / "no-reference-cache" / "%2s/%2s/%s"))
    config = {
        "tools": {"samtools": "samtools"},
        "bam": {"unmapped_scan": "auto"},
        "cram": {
            "local_ref_path": str(tmp_path / "no-reference-cache" / "%2s/%2s/%s"),
            "reference_probe_timeout_seconds": 10,
            "unmapped_scan": "auto",
            "reference_candidate_order": ["cli", "htslib_resolved"],
        },
    }
    plan = run_preflight(
        str(alignment),
        str(tmp_path / "run"),
        "input",
        file_format,
        config,
        2,
        region="chr1:155158000-155163000",
        reference_fasta=str(reference) if file_format == "cram" else None,
        header_contigs=("chr1",),
    )
    output = tmp_path / "unmapped.bam"
    command = build_cram_unmapped_filter_command(
        samtools_path="samtools",
        in_bam=plan.view_path,
        unmapped_bam=output,
        threads=2,
        reference_path=plan.reference_path,
    )
    completed = subprocess.run(
        command,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
        check=False,
    )

    expected = int(subprocess.check_output(["samtools", "view", "-c", "-f", "4", str(source_bam)], text=True))
    recovered = int(subprocess.check_output(["samtools", "view", "-c", str(output)], text=True))
    assert plan.unmapped_scan == "stream"
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert expected == 4807
    assert recovered == expected
