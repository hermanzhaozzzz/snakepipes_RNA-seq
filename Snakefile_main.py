# _*_ coding: UTF-8 _*_

########################################################################
# ZHAO Huanan
# 2021-02-01
# pcif1 RNA-Seq analysis pipeline
######################################################################## 
# before this, make sure you have done fastqc + multiqc 
# https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc
# and add trim rule to make a better trim
# --------------------------------------------------------------->>>>>>>
# pipeline
# --------------------------------------------------------------->>>>>>>
# 1. cutadapt
# 2. STAR alingment 
# 3. samtools addreplacerg - > RG
# 4. samtools sort by position
# 5. picard mark duplicates
# 6. samtools build bam index
# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
# make sure the fastqqc and the multiqc are in you PATH
# get the application path
with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
    print('PATH cutadapt:', CUTADAPT)
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
    print('PATH samtools:', SAMTOOLS)    
with os.popen("which STAR") as path:
    STAR = path.read().strip()
    print('PATH star:', STAR)   
with os.popen("which java") as path:
    JAVA = path.read().strip()
    print('PATH java:', JAVA)

# picard version 2.23
PICARD = "/home/zhaohuanan/zhaohn_HD/1.apps/picard/picard.jar"


# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
HG38_FA_DICT = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
STAR_HG38_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/star_hg38"
HG38_GTF = "/home/zhaohuanan/zhaohn_HD/2.database/annotate_hg38/20200714_ComprehensiveGeneAnnotation_Chr_gencode.v29.annotation.gtf"

# 这里是小鼠的, 没改变量名
# HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_gencode_GRCm38.p6/GRCm38.p6.genome.fa"
# HG38_FA_DICT = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_gencode_GRCm38.p6/GRCm38.p6.genome.fa"
# STAR_HG38_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_gencode_GRCm38.p6/star_index_150bp"
# HG38_GTF = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_annotation/genome_gencode_GRCm38.p6/gencode.vM25.annotation.gtf"


# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
    "pcif1-WT-1",
    "pcif1-WT-2",
    "pcif1-WT-3",
    "pcif1-WT-4",
    "pcif1-KO-1",
    "pcif1-KO-2",
    "pcif1-KO-3",
    "pcif1-KO-4"
#    "test"
]


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        expand("../fastq/{sample}_R1.fq.gz",sample=SAMPLES),
        expand("../fastq/{sample}_R2.fq.gz",sample=SAMPLES),
        expand("../fix.fastq/293T-RNASeq-{sample}_R1_cutadapt.fq.gz",sample=SAMPLES),
        expand("../fix.fastq/293T-RNASeq-{sample}_R2_cutadapt.fq.gz",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned.out.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned.out.fix_RG.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned_sort.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{sample}_Aligned_sort.bam.bai",sample=SAMPLES),
#         expand("../fpkm/{sample}",sample=SAMPLES)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# cutadapter
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule TruSeq_cutadapt:
    input:
        "../fix.fastq/293T-RNASeq-{sample}_R1.fastq.gz",
        "../fix.fastq/293T-RNASeq-{sample}_R2.fastq.gz"
    output:
        "../fix.fastq/293T-RNASeq-{sample}_R1_cutadapt.fq.gz",
        "../fix.fastq/293T-RNASeq-{sample}_R2_cutadapt.fq.gz"
    log:
        "../fix.fastq/293T-RNASeq-{sample}_cutadapt.log"
    shell:# using illumina universal adaptor
        """
        srun -T 24 -c 24 \
        {CUTADAPT} -j 24 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 55 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
        """
# rule NexteraSeq_cutadapt:
#     input:
#         "../fastq/{sample}_R1.fq.gz",
#         "../fastq/{sample}_R2.fq.gz"
#     output:
#         "../fix.fastq/293T-RNASeq-{sample}_R1_cutadapt.fq.gz",
#         "../fix.fastq/293T-RNASeq-{sample}_R2_cutadapt.fq.gz"
#     log:
#         "../fix.fastq/293T-RNASeq-{sample}_cutadapt.log"
#     shell:# using illumina NexteraSeq adaptor
#         """
#         {CUTADAPT} -j 24 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
#         -m 55 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
#         -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
#         -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1
#         """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# STAR mapping
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule STAR_mapping:
    input:
        fq1 = "../fix.fastq/293T-RNASeq-{sample}_R1_cutadapt.fq.gz",
        fq2 = "../fix.fastq/293T-RNASeq-{sample}_R2_cutadapt.fq.gz"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned.out.bam"
    log:
        "../bam/293T-RNASeq-{sample}_Aligned.out.log"
    params:
        "../bam/293T-RNASeq-{sample}_"
    shell:
        """
        {STAR} \
        --genomeDir {STAR_HG38_INDEX} \
        --runThreadN 24 \
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
        "../bam/293T-RNASeq-{sample}_Aligned.out.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned.out.fix_RG.bam"
    params:
        tag = "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
    shell:
        """
        {SAMTOOLS} addreplacerg -r {params.tag} -@ 24 -O BAM -o {output} --reference {HG38_FA_DICT} {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position(not sort by name)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned.out.fix_RG.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned_sort.bam"
    shell:
        """
        {SAMTOOLS} sort \
        -O BAM \
        -o {output} \
        -T {output}.temp \
        -@ 24  \
        {input}
        """
#         -m 4G
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# picard mark duplicate
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_mark_duplicate:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned_sort.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.bam",
        "../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.matrix"
    log:
        "../bam/293T-RNASeq-{sample}_Aligned_sort_MarkDup.log"
    shell:
        """
        {JAVA} -Xms100g -Xmx100g -XX:ParallelGCThreads=24 \
        -jar {PICARD} MarkDuplicates \
        I={input} \
        O={output[0]} \
        M={output[1]} \
        ASO=coordinate 2>{log}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools build bam index
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_index:
    input:
        "../bam/293T-RNASeq-{sample}_Aligned_sort.bam"
    output:
        "../bam/293T-RNASeq-{sample}_Aligned_sort.bam.bai"
    shell:
        """
        {SAMTOOLS} index -@ 24 \
        {input} \
        {output}
        """
