#!/usr/bin/env python3
# vntyper/scripts/file_processing.py

import json
import importlib.resources as pkg_resources


def filter_vcf(input_path, output_path):
    """
    Filter a VCF file to extract indels (insertions and deletions) and write them to a new file.

    An indel is a type of variant where either the reference allele (REF) or the alternate allele (ALT)
    has a different length (e.g., one is a single base, and the other is multiple bases).

    Args:
        input_path (str): Path to the input VCF file.
        output_path (str): Path to the output VCF file containing only indels.
    """
    with pkg_resources.open_text("vntyper", "config.json") as f:
        config_data = json.load(f)
    snv_length = config_data.get("file_processing", {}).get("snv_length", 1)

    with open(input_path, "r") as vcf_file, open(output_path, "w") as indel_file:
        for line in vcf_file:
            if line.startswith("##"):
                indel_file.write(line)
            elif line.startswith("#CHROM"):
                indel_file.write(line)
            else:
                _, _, _, ref, alt, *_ = line.split("\t")
                if (len(ref) == snv_length and len(alt) != snv_length) or (
                    len(ref) != snv_length and len(alt) == snv_length
                ):
                    indel_file.write(line)


def filter_indel_vcf(indel_vcf, output_ins, output_del):
    """
    Further filter the indel VCF file to separate insertions and deletions into two separate files.

    Insertions and deletions are differentiated based on the length of the reference allele (REF)
    and the alternate allele (ALT):
        - An insertion occurs when the ALT allele is longer than the REF allele.
        - A deletion occurs when the REF allele is longer than or equal to the ALT allele.

    Args:
        indel_vcf (str): Path to the input VCF file containing indels.
        output_ins (str): Path to the output VCF file for insertions.
        output_del (str): Path to the output VCF file for deletions.
    """
    with pkg_resources.open_text("vntyper", "config.json") as f:
        config_data = json.load(f)
    snv_length = config_data.get("file_processing", {}).get("snv_length", 1)

    with open(indel_vcf, "r") as vcf_file, open(
        output_ins, "w"
    ) as insertion_file, open(output_del, "w") as deletion_file:
        for line in vcf_file:
            if line.startswith("##"):
                insertion_file.write(line)
                deletion_file.write(line)
            elif line.startswith("#CHROM"):
                insertion_file.write(line)
                deletion_file.write(line)
            else:
                _, _, _, ref, alt, *_ = line.split("\t")

                if len(ref) == snv_length and len(alt) > snv_length:
                    insertion_file.write(line)
                else:
                    deletion_file.write(line)
