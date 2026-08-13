"""Plan typed, ordered invocations of the vendored Kestrel genotyper."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.kestrel_command import construct_kestrel_command
from vntyper.scripts.kestrel_count import construct_kanalyze_count_command
from vntyper.scripts.kestrel_counting import attempt_directory, count_log_path, ikc_path


@dataclass(frozen=True)
class KestrelInvocation:
    """One rendered Kestrel attempt: what to run, where to log it, what to clean up.

    Attributes:
        kmer_size: The configured k-mer size this attempt is for.
        command: The Kestrel call command.
        log_file: Where the call's output goes.
        count_command: The KAnalyze count command run first, or None when split
            counting is disabled and Kestrel counts internally as it always did.
        count_log: Where the count step's output goes. Separate from ``log_file``
            because ``run_command`` opens its log with ``"w"``.
        attempt_dir: The directory this attempt owns, removed on every exit path.
    """

    kmer_size: int
    command: str
    log_file: Path
    count_command: str | None = None
    count_log: Path | None = None
    attempt_dir: Path | None = None


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
    #: Extra JVM options for the Kestrel call. The call step is single-threaded, so
    #: SerialGC avoids G1 spawning a GC worker per core against one application thread.
    java_opts_call: str = ""
    #: Extra JVM options for the KAnalyze count step, empty to keep the G1 default:
    #: that step is genuinely parallel once the thread budget is allocated to it.
    java_opts_count: str = ""
    #: ``kanalyze.jar``. Empty disables the split, as does ``split_counting``.
    kanalyze_path: str | Path = ""
    #: The run's total thread budget for the count step.
    threads: int = 4
    #: The operator kill switch. False keeps Kestrel counting internally, exactly as it
    #: did before the split, and is never selected automatically.
    split_counting: bool = True

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
    # The planner is the single place that knows attempt layout, so the loop in
    # `run_kestrel` never derives a path of its own and cannot drift from cleanup.
    split = command_arguments.split_counting and bool(command_arguments.kanalyze_path)
    return tuple(
        _plan_one(kmer_size, immutable_fastqs, output_dir, command_arguments, split) for kmer_size in normalized_kmers
    )


def _plan_one(
    kmer_size: int,
    fastq_files: tuple[str | Path, ...],
    output_dir: Path,
    command_arguments: KestrelCommandArguments,
    split: bool,
) -> KestrelInvocation:
    """Render one attempt's commands, logs and owned directory.

    Args:
        kmer_size: The k-mer size for this attempt.
        fastq_files: The sample's ordered FASTQ operands.
        output_dir: Kestrel output directory and log root.
        command_arguments: Shared calibrated command settings.
        split: Whether counting runs as its own KAnalyze step.

    Returns:
        KestrelInvocation: The immutable record for this attempt.
    """
    attempt_dir = attempt_directory(output_dir, kmer_size)
    supplied_ikc = ikc_path(output_dir, kmer_size)
    count_command = (
        construct_kanalyze_count_command(
            kmer_size=kmer_size,
            kanalyze_path=command_arguments.kanalyze_path,
            output_dir=attempt_dir,
            ikc_out=supplied_ikc,
            fastq_files=fastq_files,
            java_path=command_arguments.java_path,
            java_memory=command_arguments.java_memory,
            java_opts=command_arguments.java_opts_count,
            threads=command_arguments.threads,
        )
        if split
        else None
    )
    return KestrelInvocation(
        kmer_size=kmer_size,
        command=construct_kestrel_command(
            kmer_size=kmer_size,
            kestrel_path=command_arguments.kestrel_path,
            reference_vntr=command_arguments.reference_vntr,
            output_dir=output_dir,
            fastq_files=fastq_files,
            vcf_out=command_arguments.vcf_out,
            java_path=command_arguments.java_path,
            java_memory=command_arguments.java_memory,
            max_align_states=command_arguments.max_align_states,
            max_hap_states=command_arguments.max_hap_states,
            log_level=command_arguments.log_level,
            sample_name=command_arguments.sample_name,
            additional_settings=command_arguments.additional_settings,
            java_opts=command_arguments.java_opts_call,
            ikc_path=supplied_ikc if split else None,
        ),
        log_file=output_dir / f"kestrel_kmer_{kmer_size}.log",
        count_command=count_command,
        count_log=count_log_path(output_dir, kmer_size) if split else None,
        attempt_dir=attempt_dir if split else None,
    )
