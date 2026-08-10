"""Typed planning for ordered Kestrel invocations."""

from pathlib import Path

import pytest

from vntyper.scripts.kestrel_execution import (
    KestrelCommandArguments,
    KestrelInvocation,
    plan_kestrel_invocations,
)

pytestmark = pytest.mark.unit


def _arguments() -> KestrelCommandArguments:
    return KestrelCommandArguments(
        kestrel_path="/opt/kestrel.jar",
        reference_vntr=Path("/refs/muc1.fa"),
        vcf_out=Path("/out/output.vcf"),
        java_path="java",
        java_memory="12g",
        max_align_states=40,
        max_hap_states=40,
        log_level="INFO",
        sample_name="patient-7",
        additional_settings="--flank 5",
    )


def test_planner_returns_one_immutable_invocation_per_kmer_in_stable_order(tmp_path: Path) -> None:
    fastqs: list[str | Path] = [
        Path("/reads/R1.fastq.gz"),
        "/reads/R2.fastq.gz",
        "/reads/single.fastq.gz",
    ]

    planned = plan_kestrel_invocations(
        fastq_files=fastqs,
        kmer_sizes=[20, 25],
        output_dir=tmp_path,
        command_arguments=_arguments(),
    )

    assert planned == (
        KestrelInvocation(20, planned[0].command, tmp_path / "kestrel_kmer_20.log"),
        KestrelInvocation(25, planned[1].command, tmp_path / "kestrel_kmer_25.log"),
    )
    assert all("/reads/R1.fastq.gz /reads/R2.fastq.gz /reads/single.fastq.gz" in item.command for item in planned)
    fastqs.append("/reads/late.fastq.gz")
    assert all("late.fastq.gz" not in item.command for item in planned)


@pytest.mark.parametrize("kmer_sizes", [(), [], [0], [-1], [True], ["20"], [20, None]])
def test_planner_rejects_empty_boolean_nonpositive_or_unsupported_kmers(kmer_sizes: object, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="k-mer"):
        plan_kestrel_invocations(
            fastq_files=("reads.fastq.gz",),
            kmer_sizes=kmer_sizes,  # type: ignore[arg-type]
            output_dir=tmp_path,
            command_arguments=_arguments(),
        )


@pytest.mark.parametrize("fastq_files", ["reads.fastq.gz", Path("reads.fastq.gz")])
def test_planner_rejects_a_scalar_fastq_container(fastq_files: object, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="sequence"):
        plan_kestrel_invocations(
            fastq_files=fastq_files,  # type: ignore[arg-type]
            kmer_sizes=(20,),
            output_dir=tmp_path,
            command_arguments=_arguments(),
        )


@pytest.mark.parametrize("field", ["output_dir", "command_arguments"])
def test_planner_rejects_unsupported_planner_argument_types(field: str, tmp_path: Path) -> None:
    values = {
        "fastq_files": ("reads.fastq.gz",),
        "kmer_sizes": (20,),
        "output_dir": tmp_path,
        "command_arguments": _arguments(),
    }
    values[field] = str(tmp_path) if field == "output_dir" else {}

    with pytest.raises(ValueError, match=field):
        plan_kestrel_invocations(**values)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("java_path", 7),
        ("max_align_states", True),
        ("max_hap_states", "40"),
        ("sample_name", None),
        ("additional_settings", Path("--flank")),
    ],
)
def test_command_arguments_reject_unsupported_field_types(field: str, value: object) -> None:
    values = {
        "kestrel_path": "/opt/kestrel.jar",
        "reference_vntr": Path("/refs/muc1.fa"),
        "vcf_out": Path("/out/output.vcf"),
        "java_path": "java",
        "java_memory": "12g",
        "max_align_states": 40,
        "max_hap_states": 40,
        "log_level": "INFO",
        "sample_name": "patient-7",
        "additional_settings": "",
    }
    values[field] = value

    with pytest.raises(ValueError, match=field):
        KestrelCommandArguments(**values)  # type: ignore[arg-type]
