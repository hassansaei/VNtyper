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

def process_bam_to_fastq(in_bam, output, output_name, threads, config, reference_assembly="hg19", fast_mode=False, delete_intermediates=True, keep_intermediates=False):
    samtools_path = config["tools"]["samtools"]
    
    if reference_assembly == "hg38":
        bam_region = config["bam_processing"]["bam_region_hg38"]
    else:
        bam_region = config["bam_processing"]["bam_region_hg19"]

    # Slicing BAM region
    final_bam = f"{output}{output_name}_sliced.bam"
    if keep_intermediates and os.path.exists(final_bam):
        logging.info(f"Reusing existing BAM slice: {final_bam}")
    else:
        command_slice = f"{samtools_path} view -P -b {in_bam} {bam_region} -o {final_bam}"
        logging.info(f'Starting BAM region slicing for {reference_assembly}...')
        if sp.call(command_slice, shell=True) != 0:
            logging.error('BAM region slicing failed.')
            return None, None
        logging.info('BAM region slicing completed.')

    if fast_mode:
        logging.info('Fast mode enabled: Skipping filtering of unmapped and partially mapped reads.')
    else:
        # Filtering step
        command_filter = (
            f"{samtools_path} view -@ {threads} -h {in_bam} | tee "
            f">(samtools view -b -f 4 -F 264 -@ {threads} - -o {output}{output_name}_unmapped1.bam) "
            f">(samtools view -b -f 8 -F 260 -@ {threads} - -o {output}{output_name}_unmapped2.bam) "
            f">(samtools view -b -f 12 -F 256 -@ {threads} - -o {output}{output_name}_unmapped3.bam) "
            f"> /dev/null"
        )
        logging.info('Starting BAM filtering for unmapped and partially mapped reads...')
        if sp.call(command_filter, shell=True, executable='/bin/bash') != 0:
            logging.error('BAM filtering failed.')
            return None, None
        command_merge = f"{samtools_path} merge -@ {threads} {output}{output_name}_sliced_unmapped.bam {final_bam} " \
                        f"{output}{output_name}_unmapped1.bam {output}{output_name}_unmapped2.bam {output}{output_name}_unmapped3.bam"
        logging.info('BAM filtering completed.')

        # Merging BAM files
        logging.info('Merging BAM files...')
        if sp.call(command_merge, shell=True) != 0:
            logging.error('BAM merging failed.')
            return None, None
        logging.info('BAM merging completed.')

        final_bam = f"{output}{output_name}_sliced_unmapped.bam"

    # Sorting and converting BAM to FASTQ using pipes
    final_fastq_1 = f"{output}{output_name}_R1.fastq.gz"
    final_fastq_2 = f"{output}{output_name}_R2.fastq.gz"
    
    if keep_intermediates and os.path.exists(final_fastq_1) and os.path.exists(final_fastq_2):
        logging.info(f"Reusing existing FASTQ files: {final_fastq_1} and {final_fastq_2}")
    else:
        command_sort_fastq = (
            f"{samtools_path} sort -n -@ {threads} {final_bam} | "
            f"{samtools_path} fastq -@ {threads} - -1 {final_fastq_1} -2 {final_fastq_2} -0 /dev/null -s /dev/null"
        )
        logging.info('Sorting and converting BAM to FASTQ using pipes...')
        if sp.call(command_sort_fastq, shell=True) != 0:
            logging.error('BAM to FASTQ conversion failed.')
            return None, None
        logging.info('BAM to FASTQ conversion completed.')

    # Remove intermediate BAM files if required
    if delete_intermediates and not keep_intermediates:
        logging.info('Removing intermediate BAM files...')
        intermediate_files = [
            f"{output}{output_name}_unmapped1.bam",
            f"{output}{output_name}_unmapped2.bam",
            f"{output}{output_name}_unmapped3.bam",
        ]
        for file in intermediate_files:
            if os.path.exists(file):
                os.remove(file)
        logging.info('Intermediate files removed.')

    # Return the paths to the generated FASTQ files
    return final_fastq_1, final_fastq_2
