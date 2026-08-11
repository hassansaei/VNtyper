"""Plan typed, ordered invocations of the vendored Kestrel genotyper."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.kestrel_command import construct_kestrel_command


@dataclass(frozen=True)
class KestrelInvocation:
    """One rendered Kestrel command and its dedicated log destination."""

    kmer_size: int
    command: str
    log_file: Path


@dataclass(frozen=True)
class KestrelCommandArguments:
    """Kestrel command values shared by every configured k-mer size."""

    kestrel_path: str | Path
    reference_vntr: str | Path
    vcf_out: str | Path
    java_path: str
    java_memory: str
    max_align_states: int
    max_hap_states: int
    log_level: str
    sample_name: str
    additional_settings: str

    def __post_init__(self) -> None:
        """Reject unsupported values before command rendering."""
        for field_name in ("kestrel_path", "reference_vntr", "vcf_out"):
            value = getattr(self, field_name)
            if not isinstance(value, (str, Path)) or not str(value):
                raise ValueError(f"{field_name} must be a non-empty string or Path.")
        for field_name in ("java_path", "java_memory", "log_level", "sample_name"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{field_name} must be a non-empty string.")
        if not isinstance(self.additional_settings, str):
            raise ValueError("additional_settings must be a string.")
        for field_name in ("max_align_states", "max_hap_states"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{field_name} must be a positive integer.")


def plan_kestrel_invocations(
    *,
    fastq_files: Sequence[str | Path],
    kmer_sizes: Sequence[int],
    output_dir: Path,
    command_arguments: KestrelCommandArguments,
) -> tuple[KestrelInvocation, ...]:
    """Render one immutable invocation per configured k-mer size.

    Args:
        fastq_files: Ordered FASTQ operands belonging to one sample.
        kmer_sizes: Positive integer sizes to try in stable order.
        output_dir: Kestrel output directory and log root.
        command_arguments: Shared calibrated command settings.

    Returns:
        Immutable invocations in the same order as ``kmer_sizes``.

    Raises:
        ValueError: If the k-mer sequence or planner argument types are invalid.
    """
    if isinstance(kmer_sizes, (str, bytes)) or not isinstance(kmer_sizes, Sequence) or not kmer_sizes:
        raise ValueError("k-mer sizes must be a non-empty sequence of positive integers.")
    normalized_kmers = tuple(kmer_sizes)
    if any(isinstance(kmer, bool) or not isinstance(kmer, int) or kmer <= 0 for kmer in normalized_kmers):
        raise ValueError("k-mer sizes must be positive integers.")
    if not isinstance(output_dir, Path):
        raise ValueError("output_dir must be a Path.")
    if not isinstance(command_arguments, KestrelCommandArguments):
        raise ValueError("command_arguments must be KestrelCommandArguments.")
    if isinstance(fastq_files, (str, Path)) or not isinstance(fastq_files, Sequence):
        raise ValueError("fastq_files must be a non-scalar sequence of FASTQ paths.")

    immutable_fastqs = tuple(fastq_files)
    return tuple(
        KestrelInvocation(
            kmer_size=kmer_size,
            command=construct_kestrel_command(
                kmer_size=kmer_size,
                kestrel_path=command_arguments.kestrel_path,
                reference_vntr=command_arguments.reference_vntr,
                output_dir=output_dir,
                fastq_files=immutable_fastqs,
                vcf_out=command_arguments.vcf_out,
                java_path=command_arguments.java_path,
                java_memory=command_arguments.java_memory,
                max_align_states=command_arguments.max_align_states,
                max_hap_states=command_arguments.max_hap_states,
                log_level=command_arguments.log_level,
                sample_name=command_arguments.sample_name,
                additional_settings=command_arguments.additional_settings,
            ),
            log_file=output_dir / f"kestrel_kmer_{kmer_size}.log",
        )
        for kmer_size in normalized_kmers
    )
