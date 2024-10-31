import pandas as pd


def filter_vcf(input_path, output_path):
    """
    Filter a VCF file to extract indels (insertions and deletions) and write them to a new file.

    An indel is a type of variant where either the reference allele (REF) or the alternate allele (ALT)
    has a different length (e.g., one is a single base, and the other is multiple bases).

    Args:
        input_path (str): Path to the input VCF file.
        output_path (str): Path to the output VCF file containing only indels.
    """
    with open(input_path, "r") as vcf_file, open(output_path, "w") as indel_file:
        for line in vcf_file:
            if line.startswith("##"):
                # Metadata lines that begin with "##" are written as-is to the output file
                indel_file.write(line)
            elif line.startswith("#CHROM"):
                # The header line that starts with "#CHROM" contains the column names.
                # This line must be written to the output file.
                indel_file.write(line)
            else:
                # Split the line into columns using tab as the delimiter.
                # The VCF format is tab-delimited, and the columns are:
                # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
                # We are particularly interested in the REF (reference allele) and ALT (alternate allele) columns.
                _, _, _, ref, alt, *_ = line.split("\t")

                # Determine if the variant is an indel:
                # An indel occurs when either REF or ALT is of different length than 1 base.
                if (len(ref) == 1 and len(alt) != 1) or (len(ref) != 1 and len(alt) == 1):
                    # Write lines containing indels to the output file.
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
    with open(indel_vcf, "r") as vcf_file, \
            open(output_ins, "w") as insertion_file, \
            open(output_del, "w") as deletion_file:
        for line in vcf_file:
            if line.startswith("##"):
                # Metadata lines that begin with "##" are written as-is to both output files
                insertion_file.write(line)
                deletion_file.write(line)
            elif line.startswith("#CHROM"):
                # The header line that starts with "#CHROM" contains the column names.
                # This line must be written to both output files.
                insertion_file.write(line)
                deletion_file.write(line)
            else:
                # Split the line into columns using tab as the delimiter.
                # Extract the REF (reference allele) and ALT (alternate allele) columns.
                _, _, _, ref, alt, *_ = line.split("\t")

                # Determine if the variant is an insertion or deletion:
                if len(ref) == 1 and len(alt) > 1:
                    # If ALT is longer than REF, it's an insertion.
                    insertion_file.write(line)
                else:
                    # Otherwise, it's a deletion.
                    deletion_file.write(line)
