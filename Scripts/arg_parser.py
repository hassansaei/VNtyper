import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Given raw fastq files, this pipeline genotypes MUC1-VNTR using Kestrel and Code-adVNTR methods')
    parser.add_argument('-ref', '--reference_file', type=str, metavar='Reference', help='FASTA-formatted reference file and indexes', required=True)
    parser.add_argument('-r1', '--fastq1', type=str, default=None, help='Fastq file first pair', required=False)
    parser.add_argument('-r2', '--fastq2', type=str, default=None, help='Fastq file second pair', required=False)
    parser.add_argument('-o', '--output', type=str, default=None, help='Output file name', required=True)
    parser.add_argument('-ref_VNTR', '--reference_VNTR', type=str, metavar='Reference', help='MUC1-specific reference file', required=True)
    parser.add_argument('-t', '--threads', type=str, default=None, help='Number of threads (CPU)')
    parser.add_argument('-p', '--tools_path', type=str, default=None, help='Path to the VNtyper directory', required=True)
    parser.add_argument('-w', '--working_dir', type=str, default=None, help='The path to the output', required=True)
    parser.add_argument('-m', '--reference_vntr', type=str, default=None, help='adVNTR reference VNTR database', required=False)
    parser.add_argument("--ignore_advntr", action="store_true", help="Skip adVNTR genotyping of MUC1-VNTR")
    parser.add_argument("--bam", action="store_true", help="BAM file as an input")
    parser.add_argument("--fastq", action="store_true", help="Paired-end fastq files as an input")
    parser.add_argument('-a', '--alignment', type=str, default=None, help='Alignment File (with an index file .bai)', required=False)

    return parser.parse_args()
