"""Execute Kestrel 1.0.1 against every non-empty BAM-conversion FASTQ."""

from __future__ import annotations

import hashlib
import shlex
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.fastq_bam_processing import process_bam_to_fastq
from vntyper.scripts.kestrel_execution import KestrelCommandArguments, plan_kestrel_invocations
from vntyper.scripts.kestrel_genotyping import kestrel_config, process_kestrel_output
from vntyper.scripts.pipeline_read_routing import count_fastq_records
from vntyper.scripts.utils import load_config

pytestmark = pytest.mark.integration

REPOSITORY = Path(__file__).resolve().parents[2]
KESTREL_SHA256 = "ceca2682578d641323405b3e605c85eebbe181574d056f0fabdf0ea2bf0bd1f1"


def _registered_bam(test_config: dict[str, Any], filename: str) -> Path:
    """Resolve one BAM only through the canonical downloaded-data registry."""
    matches = [
        REPOSITORY / resource["local_path"] / resource["filename"]
        for resource in test_config["file_resources"]
        if resource["filename"] == filename
    ]
    assert len(matches) == 1, f"Expected exactly one registered resource named {filename}, found {matches}"
    return matches[0]


def _convert_to_four_fastqs(source: Path, assembly: str, output: Path, config: dict[str, Any]) -> tuple[Path, ...]:
    """Run the production fast-mode BAM conversion and return all four destinations."""
    output.mkdir(parents=True)
    paths = process_bam_to_fastq(
        output=output,
        output_name="output",
        threads=2,
        config=config,
        plan=AlignmentPlan(
            input_path=str(source),
            view_path=str(source),
            file_format="bam",
            index_path=str(source) + ".bai",
            reference_path=None,
            reference_source="not-required",
            uncovered_contigs=(),
            unmapped_scan="indexed",
        ),
        reference_assembly=assembly,
        fast_mode=True,
        delete_intermediates=False,
    )
    return tuple(Path(path) for path in paths)


def _vcf_records(path: Path) -> Counter[tuple[str, ...]]:
    """Parse stable VCF tuples, excluding headers and sample-depth serialization."""
    records: Counter[tuple[str, ...]] = Counter()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            fields = line.split("\t")
            records[tuple(fields[:9])] += 1
    return records


def _run_kestrel(
    *,
    fastqs: tuple[Path, ...],
    output: Path,
    sample_name: str,
    config: dict[str, Any],
) -> tuple[Counter[tuple[str, ...]], str]:
    """Build and execute one calibrated real Kestrel invocation."""
    output.mkdir(parents=True)
    settings = kestrel_config["kestrel_settings"]
    reference_vntr = REPOSITORY / config["reference_data"]["muc1_reference_vntr"]
    kestrel_jar = REPOSITORY / config["tools"]["kestrel"]
    vcf = output / "output.vcf"
    invocation = plan_kestrel_invocations(
        fastq_files=fastqs,
        kmer_sizes=settings["kmer_sizes"],
        output_dir=output,
        command_arguments=KestrelCommandArguments(
            kestrel_path=kestrel_jar,
            reference_vntr=reference_vntr,
            vcf_out=vcf,
            java_path=config["tools"]["java_path"],
            java_memory=settings["java_memory"],
            max_align_states=settings["max_align_states"],
            max_hap_states=settings["max_hap_states"],
            log_level="INFO",
            sample_name=sample_name,
            additional_settings=settings["additional_settings"],
        ),
    )[0]

    command_words = shlex.split(invocation.command)
    assert sum(word.startswith("-s") and not word.startswith("--") for word in command_words) == 1
    assert all(command_words.count(str(path)) == 1 for path in fastqs)
    completed = subprocess.run(
        invocation.command,
        cwd=REPOSITORY,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    sam = output / "output.sam"
    assert vcf.is_file() and vcf.stat().st_size > 0
    assert sam.is_file() and sam.stat().st_size > 0

    vcf_text = vcf.read_text(encoding="utf-8")
    assert "##source=Kestrel1.0.1" in vcf_text
    process_kestrel_output(output, vcf, reference_vntr, kestrel_config, config)
    result = pd.read_csv(output / "kestrel_result.tsv", sep="\t", comment="#")
    assert len(result) == 1
    return _vcf_records(vcf), str(result.iloc[0]["Confidence"])


@pytest.mark.parametrize(
    ("filename", "assembly", "expected_pair_records", "expected_additional_records", "expected_confidence"),
    [
        ("example_b178_hg19_subset.bam", "hg19", 6_280, 0, "High_Precision*"),
        ("example_40cf_hg38_subset.bam", "hg38", 9_069, 9, "Negative"),
    ],
)
def test_every_nonempty_conversion_fastq_is_one_kestrel_sample(
    tmp_path: Path,
    test_config: dict[str, Any],
    ensure_test_data: None,
    filename: str,
    assembly: str,
    expected_pair_records: int,
    expected_additional_records: int,
    expected_confidence: str,
) -> None:
    """Real single-read FASTQs are retained without changing calibrated calls."""
    del ensure_test_data
    config = load_config(REPOSITORY / "vntyper/config.json")
    source = _registered_bam(test_config, filename)
    assert source.is_file() and Path(str(source) + ".bai").is_file()
    fastqs = _convert_to_four_fastqs(source, assembly, tmp_path / "conversion", config)
    record_counts = tuple(count_fastq_records(path, lines_per_record=4) for path in fastqs)
    assert len(fastqs) == 4
    assert record_counts[0] > 0 and record_counts[1] > 0
    assert record_counts[2] == 0 and record_counts[3] > 0

    pair_fastqs = fastqs[:2]
    all_nonempty = tuple(path for path, count in zip(fastqs, record_counts, strict=True) if count > 0)
    assert all_nonempty == (fastqs[0], fastqs[1], fastqs[3])
    jar = REPOSITORY / config["tools"]["kestrel"]
    assert hashlib.sha256(jar.read_bytes()).hexdigest() == KESTREL_SHA256

    pair_records, pair_confidence = _run_kestrel(
        fastqs=pair_fastqs,
        output=tmp_path / "pair-only",
        sample_name=source.stem,
        config=config,
    )
    all_records, all_confidence = _run_kestrel(
        fastqs=all_nonempty,
        output=tmp_path / "all-nonempty",
        sample_name=source.stem,
        config=config,
    )

    assert sum(pair_records.values()) == expected_pair_records
    assert pair_records <= all_records
    assert sum((all_records - pair_records).values()) == expected_additional_records
    if expected_additional_records == 0:
        assert all_records == pair_records
    assert pair_confidence == expected_confidence
    assert all_confidence == expected_confidence
