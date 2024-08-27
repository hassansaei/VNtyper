import subprocess as sp
import logging
import os

def process_fastq(fastq_1, fastq_2, threads, output, output_name, config):
    # Extract paths and parameters from the config
    fastp_path = config["tools"]["fastp"]
    compression_level = config["bam_processing"]["compression_level"]
    disable_adapter_trimming = config["bam_processing"]["disable_adapter_trimming"]
    deduplication = config["bam_processing"]["deduplication"]
    dup_calc_accuracy = config["bam_processing"]["dup_calc_accuracy"]
    length_required = config["bam_processing"]["length_required"]
    temp_directory = config["temp_directory"]
    
    # Construct the fastp command using parameters from the config
    QC_command = f"{fastp_path} --thread {threads} --in1 {fastq_1} --in2 {fastq_2} --out1 {output}{output_name}_R1.fastq.gz --out2 {output}{output_name}_R2.fastq.gz --compression {compression_level}"
    
    # Add optional parameters based on config
    if disable_adapter_trimming:
        QC_command += " --disable_adapter_trimming"
    if deduplication:
        QC_command += " --dedup"
    QC_command += f" --dup_calc_accuracy {dup_calc_accuracy} --length_required {length_required} --html {output}{temp_directory}{output_name}.html"
    
    # Start FASTQ quality control
    logging.info('Starting quality control for FASTQ files...')
    process = sp.Popen(QC_command, shell=True)
    process.wait()
    logging.info('Quality control passed for FASTQ files.')

def process_bam_to_fastq(in_bam, output, output_name, threads, config):
    # Extract paths and parameters from the config
    samtools_path = config["tools"]["samtools"]
    bam_region = config["bam_processing"]["bam_region"]

    # Use samtools to slice the BAM region (e.g., extracting a specific region of the genome)
    command_slice = f"{samtools_path} view -b {in_bam} {bam_region} -o {output}{output_name}_chr1.bam"
    
    # Execute BAM slicing
    logging.info('Starting BAM region slicing...')
    process = sp.Popen(command_slice, shell=True)
    process.wait()
    logging.info('BAM region slicing completed.')

    # Use tee with process substitution to pipe the BAM to multiple samtools view processes concurrently
    command_filter = (
        f"{samtools_path} view -@ {threads} -h {in_bam} | tee "
        f">(samtools view -b -f 4 -F 264 -@ {threads} - -o {output}{output_name}_unmapped1.bam) "
        f">(samtools view -b -f 8 -F 260 -@ {threads} - -o {output}{output_name}_unmapped2.bam) "
        f">(samtools view -b -f 12 -F 256 -@ {threads} - -o {output}{output_name}_unmapped3.bam) "
        f"> /dev/null"
    )

    logging.info('Starting BAM filtering for unmapped and partially mapped reads...')
    process = sp.Popen(command_filter, shell=True, executable='/bin/bash')  # Use bash to interpret process substitution
    process.wait()
    logging.info('BAM filtering completed.')

    # Use samtools to merge the sliced BAM file and the filtered unmapped BAM files into one
    command_merge = f"{samtools_path} merge -@ {threads} {output}{output_name}_vntyper.bam {output}{output_name}_chr1.bam " \
                    f"{output}{output_name}_unmapped1.bam {output}{output_name}_unmapped2.bam {output}{output_name}_unmapped3.bam"

    logging.info('Merging BAM files...')
    process = sp.Popen(command_merge, shell=True)
    process.wait()
    logging.info('BAM merging completed.')

    # Sort the merged BAM file by read name (-n) and then convert it to paired-end FASTQ files
    command_sort_fastq = f"{samtools_path} sort -n -@ {threads} {output}{output_name}_vntyper.bam -o {output}{output_name}_VN.bam && " \
                         f"{samtools_path} fastq -@ {threads} {output}{output_name}_VN.bam -1 {output}{output_name}_R1.fastq.gz -2 {output}{output_name}_R2.fastq.gz"

    logging.info('Sorting and converting BAM to FASTQ...')
    process = sp.Popen(command_sort_fastq, shell=True)
    process.wait()
    logging.info('BAM to FASTQ conversion completed.')

    # Remove intermediate BAM files after successful conversion
    logging.info('Removing intermediate BAM files...')
    for ext in ["*.bam", "*.bai"]:
        files_to_remove = f"{output}{output_name}{ext}"
        command_rm = f"rm -f {files_to_remove}"  # Using -f flag to force remove and ignore nonexistent files
        process = sp.Popen(command_rm, shell=True)
        process.wait()
        logging.info(f"Removed {ext} files (if they existed).")
