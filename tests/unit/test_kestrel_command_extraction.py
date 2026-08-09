"""The extracted Kestrel command builder remains the public implementation."""

from __future__ import annotations

import pytest

from vntyper.scripts import kestrel_genotyping
from vntyper.scripts.kestrel_command import construct_kestrel_command

pytestmark = pytest.mark.unit


def test_kestrel_genotyping_reexports_the_extracted_builder() -> None:
    assert kestrel_genotyping.construct_kestrel_command is construct_kestrel_command


def test_extraction_preserves_the_existing_paired_command_byte_for_byte() -> None:
    command = construct_kestrel_command(
        kmer_size=20,
        kestrel_path="/opt/kestrel.jar",
        reference_vntr="/ref/muc1.fa",
        output_dir="/out",
        fastq_1="/in/R1.fastq.gz",
        fastq_2="/in/R2.fastq.gz",
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
        "--hapfmt sam -p /out/output.sam --logstderr --logstdout "
        "--loglevel INFO --temploc /out"
    )
