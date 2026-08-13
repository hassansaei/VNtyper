"""The extracted Kestrel command builder remains the public implementation."""

from __future__ import annotations

import pytest

from vntyper.scripts import kestrel_genotyping
from vntyper.scripts.kestrel_command import construct_kestrel_command

pytestmark = pytest.mark.unit


def test_kestrel_genotyping_reexports_the_extracted_builder() -> None:
    assert kestrel_genotyping.construct_kestrel_command is construct_kestrel_command


def test_the_paired_command_without_a_supplied_ikc_is_pinned_byte_for_byte() -> None:
    """The legacy shape, still reachable through ``kestrel_settings.split_counting``.

    ``--logstderr`` is gone from all three expectations in this module and that is a
    clarity fix rather than a behaviour change: ``OptLogStderr`` and ``OptLogStdout``
    call the same setter on a single ``logFile`` field, so last-one-wins and
    ``--logstdout`` -- emitted second -- has always won. ``run_command`` merges stderr
    into stdout regardless.
    """
    command = construct_kestrel_command(
        kmer_size=20,
        kestrel_path="/opt/kestrel.jar",
        reference_vntr="/ref/muc1.fa",
        output_dir="/out",
        fastq_files=("/in/R1.fastq.gz", "/in/R2.fastq.gz"),
        vcf_out="/out/output.vcf",
        java_path="/usr/bin/java",
        java_memory="12g",
        max_align_states=40,
        max_hap_states=40,
        log_level="info",
        sample_name="SAMPLE1",
    )

    assert command == (
        "/usr/bin/java -Xmx12g -jar /opt/kestrel.jar -k 20 "
        "--maxalignstates 40 --maxhapstates 40 "
        "-r /ref/muc1.fa -o /out/output.vcf "
        "-sSAMPLE1 /in/R1.fastq.gz /in/R2.fastq.gz "
        "--hapfmt sam -p /out/output.sam --logstdout "
        "--loglevel INFO --temploc /out"
    )


@pytest.mark.parametrize(
    ("output_dir", "expected_sam", "expected_temp"),
    [
        (".", "./output.sam", "."),
        ("/out/", "/out//output.sam", "/out/"),
    ],
)
def test_output_directory_spelling_is_preserved_byte_for_byte(
    output_dir: str,
    expected_sam: str,
    expected_temp: str,
) -> None:
    command = construct_kestrel_command(
        kmer_size=20,
        kestrel_path="/opt/kestrel.jar",
        reference_vntr="/ref/muc1.fa",
        output_dir=output_dir,
        fastq_files=("/in/R1.fastq.gz", "/in/R2.fastq.gz"),
        vcf_out="/out/output.vcf",
        java_path="/usr/bin/java",
        java_memory="12g",
        max_align_states=40,
        max_hap_states=40,
        log_level="info",
        sample_name="SAMPLE1",
    )

    assert command == (
        "/usr/bin/java -Xmx12g -jar /opt/kestrel.jar -k 20 "
        "--maxalignstates 40 --maxhapstates 40 "
        "-r /ref/muc1.fa -o /out/output.vcf "
        "-sSAMPLE1 /in/R1.fastq.gz /in/R2.fastq.gz "
        f"--hapfmt sam -p {expected_sam} --logstdout "
        f"--loglevel INFO --temploc {expected_temp}"
    )


@pytest.mark.parametrize(
    ("fastq_files", "expected_inputs"),
    [
        (("/in/R1.fastq.gz",), "/in/R1.fastq.gz"),
        (
            ("/in/R1.fastq.gz", "/in/R2.fastq.gz", "/in/other.fastq.gz"),
            "/in/R1.fastq.gz /in/R2.fastq.gz /in/other.fastq.gz",
        ),
        (
            ("/in/R1.fastq.gz", "/in/R2.fastq.gz", "/in/other.fastq.gz", "/in/single.fastq.gz"),
            "/in/R1.fastq.gz /in/R2.fastq.gz /in/other.fastq.gz /in/single.fastq.gz",
        ),
    ],
)
def test_one_three_and_four_fastq_commands_preserve_exact_order(
    fastq_files: tuple[str, ...], expected_inputs: str
) -> None:
    command = construct_kestrel_command(
        kmer_size=20,
        kestrel_path="/opt/kestrel.jar",
        reference_vntr="/ref/muc1.fa",
        output_dir="/out",
        fastq_files=fastq_files,
        vcf_out="/out/output.vcf",
        java_path="/usr/bin/java",
        java_memory="12g",
        max_align_states=40,
        max_hap_states=40,
        log_level="info",
        sample_name="SAMPLE1",
    )

    assert command == (
        "/usr/bin/java -Xmx12g -jar /opt/kestrel.jar -k 20 "
        "--maxalignstates 40 --maxhapstates 40 "
        "-r /ref/muc1.fa -o /out/output.vcf "
        f"-sSAMPLE1 {expected_inputs} "
        "--hapfmt sam -p /out/output.sam --logstdout "
        "--loglevel INFO --temploc /out"
    )
    assert command.count("-sSAMPLE1") == 1
