# ——————————————————>>>>>>>>>>
# Project Name: RNA-Seq analysis pipeline
# Author: Hua-nan ZHAO
# E-mail: hermanzhaozzzz@gmail.com
# ——————————————————
# Update log:
#     2022-08-22: start project
# ——————————————————
# pipeline:
# 1. cutadapt
# 2. STAR alingment
# 3. samtools addreplacerg - > RG
# 4. samtools sort by position
# 5. sambamba mark duplicates and build bam index
# ——————————————————>>>>>>>>>>
import os
import json
# ------------------------------------------------------------------->>>>>>>>>>
# FUNCTIONS
# ------------------------------------------------------------------->>>>>>>>>>
def print_head(SAMPLES, MODE):
    print('----------\nSAMPLES:')
    [print('\t' + i) for i in SAMPLES]
    print('----------\nMODE:')
    print('\t' + MODE)
    print('----------\n\n')

def check_cmd(x):
    return any(
        os.access(os.path.join(path, x), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )

def check_read(x):
    if x == "PE":
        read = ['R1', 'R2']
    elif x == "SE":
        read = ['SE']
    else:
        raise ValueError()
    return read
# ------------------------------------------------------------------->>>>>>>>>>
# SAMPLE INFO
# ------------------------------------------------------------------->>>>>>>>>>
with open('./samples.json') as f:
    dt = json.loads(f.read())

SAMPLES = dt['samples']
MODE = dt['seq_mode']
THREAD = dt['thread']
READ = check_read(MODE)

print_head(SAMPLES, MODE)
print(READ)
# ------------------------------------------------------------------->>>>>>>>>>
# RUN INFO
# ------------------------------------------------------------------->>>>>>>>>>
# THREAD = os.cpu_count() - 1
# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
GENOME = dt["genome"]
STAR_INDEX = "star_index"
ANNOTATION_GTF = "annotation_gtf"
# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# check if cmd exists
assert check_cmd("fastp")
assert check_cmd("samtools")
assert check_cmd("STAR")
assert check_cmd("sambamba")

# manually set cmd path
FASTP = "fastp"
SAMTOOLS = "samtools"
STAR = "STAR"
SAMBAMBA = "sambamba"

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        expand("../bam/{sample}_Aligned_sort_rmdup.bam",sample=SAMPLES),
        expand("../bam/{sample}_Aligned_sort_rmdup.bam.bai",sample=SAMPLES),
        '../featureCounts/all_feature.txt'
# ------------------------------------------------------------------->>>>>>>>>>
# trim adaptor
# ------------------------------------------------------------------->>>>>>>>>>
rule fastp_trim_adaptor:
    input:
        fwd="../fastq/{sample}_R1.fastq.gz",
        rev="../fastq/{sample}_R2.fastq.gz"
    output:
        fwd=temp("../fix.fastq/{sample}_R1_cutadapt.fastq.gz"),
        rev=temp("../fix.fastq/{sample}_R2_cutadapt.fastq.gz"),
        html="../fix.fastq/{sample}.html",
        json="../fix.fastq/{sample}.json"
    log:
        "../fix.fastq/{sample}.log"
    message:
        "trim_adaptor {input}"
    shell:
        """
        {FASTP} -w {THREAD} -h {output[html]} -j  {output[json]} \
            -i {input[fwd]} -I {input[rev]} -o {output[fwd]} \
            -O {output[rev]} 2> {log}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# STAR mapping
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule STAR_mapping:
    input:
        fq1 = "../fix.fastq/{sample}_R1_cutadapt.fastq.gz",
        fq2 = "../fix.fastq/{sample}_R2_cutadapt.fastq.gz"
    output:
        temp("../bam/{sample}_Aligned.out.bam")
    log:
        "../bam/{sample}_Aligned.out.log"
    params:
        "../bam/{sample}_"
    shell:
        """
        {STAR} \
        --genomeDir {STAR_INDEX} \
        --runThreadN {THREAD} \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params} \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical > {log} 2>&1
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# add @RG tag (mostly for GATK SNP/SNV calling)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule add_RG_tag:
    input:
        "../bam/{sample}_Aligned.out.bam"
    output:
        temp("../bam/{sample}_Aligned.out.fix_RG.bam")
    params:
        tag = "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
    shell:
        """
        {SAMTOOLS} addreplacerg -r {params.tag} -@ {THREAD} -O BAM -o {output} --reference {GENOME} {input}
        """
rule filter_bam:
    input:
        "../bam/{sample}_Aligned.out.fix_RG.bam"
    output:
        temp("../bam/{sample}_Aligned.out.fix_RG_filter.bam")
    shell:
        """{SAMTOOLS} view -@ {THREAD} -F 4 -F 8 -hb {input} -o {output}"""
# 12
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position(not sort by name)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:
        "../bam/{sample}_Aligned.out.fix_RG_filter.bam"
    output:
        "../bam/{sample}_Aligned_sort.bam"
    shell:
        """
        {SAMTOOLS} sort \
        -O BAM \
        -o {output} \
        -T {output}.temp \
        -@ {THREAD}  \
        {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# sambamba rmdup and build bam index
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule sambamba_rmdup_and_build_index:
    input:
        "../bam/{sample}_Aligned_sort.bam"
    output:
        "../bam/{sample}_Aligned_sort_rmdup.bam",
        "../bam/{sample}_Aligned_sort_rmdup.bam.bai"
    log:
        "../bam/{sample}_Aligned_sort_rmdup.log"
    shell:
        """
        {SAMBAMBA} markdup \
                --remove-duplicates \
                --nthreads={THREAD} \
                --show-progress \
                --sort-buffer-size 8192 \
                {input} {output[0]} > {log} 2>&1
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# featureCounts
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# https://www.jianshu.com/p/9cc4e8657d62
# -T 使用的线程数
# -p 如果是paird end 就用, 只能用在paired-end的情况中，会统计fragment而不统计read
# -B 在-p选择的条件下，只有两端read都比对上的fragment才会被统计
# -t 将exon作为一个feature
# -g 将gene_id作为一个feature
# -a 参考的gtf/gff
# -o 输出文件
# 最后加上bam文件，有几个就加几个
rule featureCounts:
    input:
        expand("../bam/{sample}_Aligned_sort_rmdup.bam",sample=SAMPLES)
    output:
        '../featureCounts/all_feature.txt'
    log:
        '../featureCounts/run_FC.log'
    shell:
        """
        {FEATURECOUNTS} \
            -T {THREADS} \
            -p \
            -t exon \
            -g gene_id \
            -a {ANNOTATION_GTF} \
            -o {output} \
            {input} 2>{log}
        """