"""Real read-only alignment-view proofs for milestone-4 index handling."""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from vntyper.scripts.alignment_contract import AlignmentPlan, index_candidate_names
from vntyper.scripts.alignment_preflight import run_preflight
from vntyper.scripts.command_builders import (
    build_cram_unmapped_filter_command,
    build_samtools_depth_command,
    build_samtools_slice_command,
)
from vntyper.scripts.fastq_bam_processing import process_bam_to_fastq

pytestmark = pytest.mark.integration


def _tree_digest(root: Path) -> dict[str, str]:
    """Return content digests for every regular file beneath ``root``."""
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _write_all_unplaced_bam(path: Path) -> tuple[str, ...]:
    """Write five literal flag-4 records beneath one declared sequence."""
    names = ("unplaced-0", "unplaced-1", "unplaced-2", "unplaced-3", "unplaced-4")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as alignment:
        for name in names:
            read = pysam.AlignedSegment(alignment.header)
            read.query_name = name
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.flag = 4
            read.reference_id = -1
            read.reference_start = -1
            alignment.write(read)
    return names


def _bam_qnames(path: Path) -> tuple[str, ...]:
    """Read a BAM's sorted QNAME multiset through real samtools."""
    output = subprocess.check_output(["samtools", "view", str(path)], text=True)
    return tuple(sorted(line.partition("\t")[0] for line in output.splitlines()))


def test_indexed_bam_recovery_handles_a_nonempty_all_unplaced_alignment(tmp_path: Path) -> None:
    """MED4: indexed htslib recovery returns all five literal QNAMEs exactly once."""
    input_root = tmp_path / "patient-input"
    preflight_root = tmp_path / "preflight-output"
    conversion_root = tmp_path / "run-output"
    input_root.mkdir()
    preflight_root.mkdir()
    conversion_root.mkdir()
    source = input_root / "all-unplaced.bam"
    expected_names = tuple(sorted(_write_all_unplaced_bam(source)))
    config = {"tools": {"samtools": "samtools"}, "bam": {"unmapped_scan": "auto"}}
    plan = run_preflight(
        str(source),
        str(preflight_root),
        "input",
        "bam",
        config,
        2,
        region="chr1:1-10",
        header_contigs=("chr1",),
    )
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t0\t10\n", encoding="utf-8")

    process_bam_to_fastq(
        output=str(conversion_root),
        output_name="output",
        threads=2,
        config=config,
        plan=plan,
        fast_mode=False,
        keep_intermediates=True,
        delete_intermediates=False,
        bed_file=bed,
    )
    plan.close()

    recovered = conversion_root / "output_unmapped.bam"
    final_bam = conversion_root / "output_sliced.bam"
    assert plan.unmapped_scan == "indexed"
    assert _bam_qnames(source) == expected_names
    assert _bam_qnames(recovered) == expected_names
    assert _bam_qnames(final_bam) == expected_names


@pytest.mark.parametrize("file_format", ["bam", "cram"])
def test_real_pipeline_preserves_a_read_only_alignment_tree(
    tmp_path: Path,
    ensure_test_data: None,
    file_format: str,
) -> None:
    """A-225-3 runs the CLI through conversion while patient bytes remain unchanged."""
    repository = Path(__file__).parents[2]
    input_root = tmp_path / "patient-input"
    input_root.mkdir()
    output_root = tmp_path / "run"
    alignment = input_root / f"sample.{file_format}"
    if file_format == "bam":
        source = repository / "tests/data/derived/single_end/example_b178_hg19_single_end.bam"
        assert source.is_file(), f"Missing derived fixture {source}; run `make cram-fixtures`."
        reference_args: list[str] = []
        target_args: list[str] = []
        config_path = repository / "vntyper/config.json"
    else:
        fixture_root = repository / "tests/data/cram/reference-dependent"
        source = fixture_root / "reference-dependent.cram"
        reference = fixture_root / "reference.fa"
        assert source.is_file() and reference.is_file(), "Missing purpose fixtures; run `make cram-fixtures`."
        reference_args = ["--reference-fasta", str(reference)]
        target_args = ["--custom-regions", "chr1:1-10000"]
        config_payload = json.loads((repository / "vntyper/config.json").read_text(encoding="utf-8"))
        config_payload["bam_processing"]["assemblies"]["GRCh37"]["vntr_region_coords"] = "1-10000"
        config_path = tmp_path / "cram-config.json"
        config_path.write_text(json.dumps(config_payload), encoding="utf-8")
    shutil.copyfile(source, alignment)
    alignment.chmod(0o444)
    input_root.chmod(0o555)
    before = _tree_digest(input_root)
    command = [
        sys.executable,
        "-m",
        "vntyper.cli",
        "--config-path",
        str(config_path),
        "pipeline",
        f"--{file_format}",
        str(alignment),
        "--reference-assembly",
        "hg19",
        "--threads",
        "2",
        "--keep-intermediates",
        "--output-dir",
        str(output_root),
        *reference_args,
        *target_args,
    ]

    try:
        completed = subprocess.run(command, cwd=repository, capture_output=True, text=True, check=False)
    finally:
        input_root.chmod(0o755)
        alignment.chmod(0o644)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert (output_root / "fastq_bam_processing/output_sliced.bam").is_file()
    assert all(not os.path.lexists(candidate) for candidate in index_candidate_names(str(alignment), file_format))
    assert _tree_digest(input_root) == before


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

    plan: AlignmentPlan | None = None
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
            index_path=plan.stable_index_path,
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
        if plan is not None:
            plan.close()
        input_root.chmod(0o755)
        alignment.chmod(0o644)

    assert plan is not None
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert int(subprocess.check_output(["samtools", "view", "-c", str(slice_path)], text=True)) > 0
    assert Path(plan.index_path).is_file()
    assert not Path(plan.index_path).is_symlink()
    assert all(not os.path.lexists(candidate) for candidate in index_candidate_names(str(alignment), file_format))
    assert _tree_digest(input_root) == before


def test_unindexed_read_only_reference_uses_a_run_local_index_for_every_consumer(
    tmp_path: Path,
    ensure_test_data: None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unindexed operator FASTA stays unchanged while probe and slice decode it."""
    repository = Path(__file__).parents[2]
    fixture_root = repository / "tests/data/cram/reference-dependent"
    operator_root = tmp_path / "operator-reference"
    output_root = tmp_path / "run"
    operator_root.mkdir()
    output_root.mkdir()
    reference = operator_root / "reference.fa"
    alignment = tmp_path / "input.cram"
    shutil.copyfile(fixture_root / "reference.fa", reference)
    shutil.copyfile(fixture_root / "reference-dependent.cram", alignment)
    before = _tree_digest(operator_root)
    reference.chmod(0o444)
    operator_root.chmod(0o555)
    monkeypatch.setenv("REF_PATH", str(tmp_path / "no-reference-cache" / "%2s/%2s/%s"))
    config = {
        "tools": {"samtools": "samtools"},
        "cram": {
            "allow_ambient_reference_resolution": False,
            "local_ref_path": str(tmp_path / "no-reference-cache" / "%2s/%2s/%s"),
            "reference_probe_timeout_seconds": 10,
            "unmapped_scan": "auto",
            "reference_candidate_order": ["cli", "htslib_resolved"],
        },
    }

    plan: AlignmentPlan | None = None
    try:
        plan = run_preflight(
            str(alignment),
            str(output_root),
            "input",
            "cram",
            config,
            2,
            region="chr1:1-10000",
            reference_fasta=str(reference),
            header_contigs=("chr1",),
        )
        slice_path = output_root / "slice.bam"
        command = build_samtools_slice_command(
            samtools_path="samtools",
            in_bam=plan.view_path,
            index_path=plan.stable_index_path,
            output_bam=slice_path,
            region="chr1:1-10000",
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
        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert int(subprocess.check_output(["samtools", "view", "-c", str(slice_path)], text=True)) == 50
        assert plan.reference_source == "cli"
        assert plan.reference_path is not None
        assert plan.reference_path != str(reference)
        assert Path(plan.reference_path).is_symlink()
        # The generated index is retained under its own name, not replaced by a link to
        # this process's descriptor for it, which would point the entry at itself (#238).
        generated_index = Path(f"{plan.reference_path}.fai")
        assert not generated_index.is_symlink()
        assert generated_index.stat(follow_symlinks=False).st_nlink == 1
        assert generated_index.read_text(encoding="utf-8").startswith("chr1\t")
        assert not Path(f"{reference}.fai").exists()
    finally:
        if plan is not None:
            plan.close()
        operator_root.chmod(0o755)
        reference.chmod(0o644)

    assert _tree_digest(operator_root) == before


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
    plan.close()

    assert completed.returncode == 0
    assert int(subprocess.check_output(["samtools", "view", "-c", str(slice_path)], text=True)) == 29736
    assert not Path(plan.index_path).is_symlink()
    assert supplied_index.samefile(input_root / "sample.bam.bai")
    assert _tree_digest(input_root) == before


def test_atomic_bai_replacement_after_preflight_still_slices_the_proven_index_inode(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """Post-preflight pathname replacement must not reopen a wrong valid BAI."""
    repository = Path(__file__).parents[2]
    input_root = tmp_path / "patient-input"
    preflight_root = tmp_path / "preflight"
    conversion_root = tmp_path / "conversion"
    input_root.mkdir()
    preflight_root.mkdir()
    conversion_root.mkdir()
    alignment = input_root / "sample.bam"
    shutil.copyfile(repository / "tests/data/example_b178_hg19_subset.bam", alignment)
    plan = run_preflight(
        str(alignment),
        str(preflight_root),
        "input",
        "bam",
        {"tools": {"samtools": "samtools"}, "bam": {"unmapped_scan": "auto"}},
        2,
        region="chr1:155158000-155163000",
    )
    wrong_index = tmp_path / "wrong-40cf.bai"
    shutil.copyfile(repository / "tests/data/example_40cf_hg38_subset.bam.bai", wrong_index)
    wrong_index.replace(plan.index_path)
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t155157999\t155163000\n", encoding="utf-8")

    try:
        process_bam_to_fastq(
            output=conversion_root,
            output_name="output",
            threads=2,
            config={"tools": {"samtools": "samtools"}},
            plan=plan,
            fast_mode=True,
            keep_intermediates=True,
            bed_file=bed,
        )
        depth_output = tmp_path / "depth.tsv"
        depth_command = build_samtools_depth_command(
            samtools_path="samtools",
            threads=2,
            region="chr1:155158000-155163000",
            bam_file=plan.view_path,
            index_path=plan.stable_index_path,
            coverage_output=depth_output,
        )
        subprocess.run(depth_command, shell=True, executable="/bin/bash", check=True)
    finally:
        plan.close()

    sliced = conversion_root / "output_sliced.bam"
    assert int(subprocess.check_output(["samtools", "view", "-c", str(sliced)], text=True)) == 29_736
    assert len(depth_output.read_text(encoding="utf-8").splitlines()) == 5_001


def test_atomic_bam_replacement_after_preflight_slices_the_bound_original(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """A fresh BAI and later real slice must consume the same opened BAM inode."""
    repository = Path(__file__).parents[2]
    original = repository / "tests/data/example_b178_hg19_subset.bam"
    replacement_source = repository / "tests/data/example_40cf_hg38_subset.bam"
    original_digest = hashlib.sha256(original.read_bytes()).hexdigest()
    replacement_digest = hashlib.sha256(replacement_source.read_bytes()).hexdigest()
    input_root = tmp_path / "patient-input"
    output_root = tmp_path / "run"
    input_root.mkdir()
    output_root.mkdir()
    alignment = input_root / "sample.bam"
    replacement = input_root / "replacement.bam"
    shutil.copyfile(original, alignment)
    shutil.copyfile(replacement_source, replacement)
    plan = run_preflight(
        str(alignment),
        str(output_root),
        "input",
        "bam",
        {"tools": {"samtools": "samtools"}, "bam": {"unmapped_scan": "auto"}},
        2,
        region="chr1:155158000-155163000",
    )
    replacement.replace(alignment)
    slice_path = output_root / "slice.bam"
    command = build_samtools_slice_command(
        samtools_path="samtools",
        in_bam=plan.view_path,
        output_bam=slice_path,
        region="chr1:155158000-155163000",
        threads=2,
    )

    try:
        completed = subprocess.run(
            command,
            shell=True,
            executable="/bin/bash",
            capture_output=True,
            text=True,
            check=False,
        )
    finally:
        plan.close()

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert int(subprocess.check_output(["samtools", "view", "-c", str(slice_path)], text=True)) == 29736
    assert hashlib.sha256(alignment.read_bytes()).hexdigest() == replacement_digest
    assert hashlib.sha256(original.read_bytes()).hexdigest() == original_digest


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
    plan.close()

    expected = int(subprocess.check_output(["samtools", "view", "-c", "-f", "4", str(source_bam)], text=True))
    recovered = int(subprocess.check_output(["samtools", "view", "-c", str(output)], text=True))
    assert plan.unmapped_scan == "stream"
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert expected == 4807
    assert recovered == expected


def test_zero_placed_cram_with_incomplete_literal_star_fetch_selects_stream(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """A real multi-slice CRAM disproves zero placed count as sufficient proof."""
    repository = Path(__file__).parents[2]
    source = repository / "tests/data/remapped/bwa/hg38_ensembl/example_7a61_hg38_ensembl_bwa.bam"
    assert source.is_file(), f"Missing registered fixture {source}; run `make download-test-data`."
    alignment = tmp_path / "zero-placed.cram"
    quoted_source = shlex.quote(str(source))
    quoted_alignment = shlex.quote(str(alignment))
    derivation = (
        "set -o pipefail; "
        f"{{ samtools view -H {quoted_source}; samtools view -F 4 {quoted_source}; "
        f"samtools view -f 4 {quoted_source} '*'; }} | "
        f"samtools view -C --output-fmt-option no_ref=1 -o {quoted_alignment} -"
    )
    subprocess.run(derivation, shell=True, executable="/bin/bash", check=True)
    subprocess.run(["samtools", "index", str(alignment)], check=True)

    idxstats = subprocess.check_output(["samtools", "idxstats", str(alignment)], text=True)
    rows = [line.split("\t") for line in idxstats.splitlines()]
    placed = sum(int(fields[3]) for fields in rows if fields[0] != "*")
    unplaced = int(next(fields[3] for fields in rows if fields[0] == "*"))
    whole = int(subprocess.check_output(["samtools", "view", "-c", "-f", "4", str(alignment)], text=True))
    literal_star = int(subprocess.check_output(["samtools", "view", "-c", "-f", "4", str(alignment), "*"], text=True))
    assert (placed, unplaced, whole, literal_star) == (0, 622_690, 622_690, 2_690)

    plan = run_preflight(
        str(alignment),
        str(tmp_path / "run"),
        "input",
        "cram",
        {
            "tools": {"samtools": "samtools"},
            "cram": {
                "local_ref_path": str(tmp_path / "no-reference-cache" / "%2s/%2s/%s"),
                "reference_probe_timeout_seconds": 10,
                "unmapped_scan": "auto",
                "reference_candidate_order": ["htslib_resolved"],
            },
        },
        2,
        region="1:1-10",
        header_contigs=("1",),
    )
    try:
        assert plan.unmapped_scan == "stream"
    finally:
        plan.close()


def test_nonfast_slice_and_complete_recovery_merge_is_a_disjoint_flag_four_union(
    tmp_path: Path,
    ensure_test_data: None,
) -> None:
    """HIGH1: the final BAM must contain every source flag-4 QNAME exactly once."""
    repository = Path(__file__).parents[2]
    source = repository / "tests/data/derived/single_end/example_b178_hg19_single_end.bam"
    assert source.is_file(), f"Missing derived fixture {source}; run `make cram-fixtures`."
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t155157999\t155163000\n", encoding="utf-8")
    output = tmp_path / "conversion"
    output.mkdir()
    process_bam_to_fastq(
        output=str(output),
        output_name="output",
        threads=2,
        config={"tools": {"samtools": "samtools"}},
        plan=AlignmentPlan(
            input_path=str(source),
            view_path=str(source),
            file_format="bam",
            index_path=str(source) + ".bai",
            reference_path=None,
            reference_source="not-required",
            uncovered_contigs=(),
            unmapped_scan="stream",
        ),
        fast_mode=False,
        keep_intermediates=True,
        bed_file=bed,
    )

    final_bam = output / "output_sliced.bam"
    source_names = sorted(
        line.partition("\t")[0]
        for line in subprocess.check_output(["samtools", "view", "-f", "4", str(source)], text=True).splitlines()
    )
    final_names = sorted(
        line.partition("\t")[0]
        for line in subprocess.check_output(["samtools", "view", "-f", "4", str(final_bam)], text=True).splitlines()
    )

    assert len(source_names) == 4_807
    assert len(final_names) == 4_807
    assert final_names == source_names
