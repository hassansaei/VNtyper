import subprocess as sp
import logging

def process_fastq(fastq_1, fastq_2, threads, output, output_name):
    QC_command = f"fastp --thread {threads} --in1 {fastq_1} --in2 {fastq_2} --out1 {output}{output_name}_R1.fastq.gz --out2 {output}{output_name}_R2.fastq.gz --compression 6 --disable_adapter_trimming --dedup --dup_calc_accuracy 3 --length_required 40 --html {output}temp/{output_name}.html"
    logging.info('Fastq file quality control: deduplication, low quality read removal, and size correction...')
    print("Quality control step!\n")
    process = sp.Popen(QC_command, shell=True)
    process.wait()
    print("QC passed!\n")
    logging.info('Quality control passed..')

def process_bam_to_fastq(in_bam, output, output_name, threads):
    # Commands for BAM cleanup and conversion to FASTQ
    command_slice = f"/SOFT/./sambamba-0.6.8 slice {in_bam} chr1:155158000-155163000 -o {output}{output_name}_chr1.bam"
    command_flag = f"samtools view -u -f 4 -F264 -@ {threads} {in_bam} > {output}{output_name}_unmapped1.bam && " \
                   f"samtools view -u -f 8 -F260 -@ {threads} {in_bam} > {output}{output_name}_unmapped2.bam && " \
                   f"samtools view -u -f 12 -F256 -@ {threads} {in_bam} > {output}{output_name}_unmapped3.bam"
    command_merge = f"/SOFT/./sambamba-0.6.8 merge -t {threads} {output}{output_name}_vntyper.bam {output}{output_name}_chr1.bam " \
                    f"{output}{output_name}_unmapped1.bam {output}{output_name}_unmapped2.bam {output}{output_name}_unmapped3.bam"
    command_sort_fastq = f"samtools sort -n -@ {threads} {output}{output_name}_vntyper.bam -o {output}{output_name}_VN.bam && " \
                         f"samtools fastq -@ {threads} {output}{output_name}_VN.bam -1 {output}{output_name}_R1.fastq.gz -2 {output}{output_name}_R2.fastq.gz"

    logging.info('BAM file cleanup and converting to fastq...')
    print('BAM cleanup and converting to fastq...')
    process = sp.Popen(command_slice, shell=True)
    process.wait()
    process = sp.Popen(command_flag, shell=True)
    process.wait()
    process = sp.Popen(command_merge, shell=True)
    process.wait()
    process = sp.Popen(command_sort_fastq, shell=True)
    process.wait()
    logging.info('BAM to Fastq conversion finished!')
    print('BAM to Fastq conversion finished!')

    # Remove intermediate BAM files
    command_rm = f"rm {output}{output_name}*.bam && rm {output}{output_name}*.bai"
    process = sp.Popen(command_rm, shell=True)
    process.wait()
