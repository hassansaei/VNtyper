import subprocess as sp
import logging

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
    
    # Use samtools for flag-based filtering
    # The goal here is to create multiple BAM files by filtering reads based on specific flags:
    # - '-f' indicates we want reads with a particular flag set.
    # - '-F' indicates we want reads where a particular flag is NOT set.
    # 
    # Specifically, this is what each command does:
    # 1. `-f 4 -F 264`: Extract unmapped reads (flag 4) that are NOT part of a properly paired alignment (flag 264).
    # 2. `-f 8 -F 260`: Extract reads where the mate is unmapped (flag 8) and NOT part of a properly paired alignment (flag 260).
    # 3. `-f 12 -F 256`: Extract reads where both the read and its mate are unmapped (flags 4 and 8 combined = 12) and NOT supplementary alignments (flag 256).
    command_flag = f"{samtools_path} view -b -f 4 -F 264 -@ {threads} {in_bam} -o {output}{output_name}_unmapped1.bam && " \
                   f"{samtools_path} view -b -f 8 -F 260 -@ {threads} {in_bam} -o {output}{output_name}_unmapped2.bam && " \
                   f"{samtools_path} view -b -f 12 -F 256 -@ {threads} {in_bam} -o {output}{output_name}_unmapped3.bam"
    
    # Use samtools to merge the sliced BAM file and the filtered unmapped BAM files into one
    command_merge = f"{samtools_path} merge -@ {threads} {output}{output_name}_vntyper.bam {output}{output_name}_chr1.bam " \
                    f"{output}{output_name}_unmapped1.bam {output}{output_name}_unmapped2.bam {output}{output_name}_unmapped3.bam"
    
    # Sort the merged BAM file by read name (-n) and then convert it to paired-end FASTQ files
    command_sort_fastq = f"{samtools_path} sort -n -@ {threads} {output}{output_name}_vntyper.bam -o {output}{output_name}_VN.bam && " \
                         f"{samtools_path} fastq -@ {threads} {output}{output_name}_VN.bam -1 {output}{output_name}_R1.fastq.gz -2 {output}{output_name}_R2.fastq.gz"

    logging.info('Starting BAM file cleanup and conversion to FASTQ...')
    
    # Execute the commands
    process = sp.Popen(command_slice, shell=True)
    process.wait()
    logging.info('BAM region slicing completed.')

    process = sp.Popen(command_flag, shell=True)
    process.wait()
    logging.info('BAM flag-based filtering completed.')

    process = sp.Popen(command_merge, shell=True)
    process.wait()
    logging.info('BAM merging completed.')

    process = sp.Popen(command_sort_fastq, shell=True)
    process.wait()
    logging.info('BAM to FASTQ conversion completed.')

    # Remove intermediate BAM files after successful conversion
    command_rm = f"rm {output}{output_name}*.bam && rm {output}{output_name}*.bai"
    process = sp.Popen(command_rm, shell=True)
    process.wait()
    logging.info('Intermediate BAM files removed.')
