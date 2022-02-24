# _*_ coding: UTF-8 _*_

########################################################################
# ZHAO Huanan
# 2022-02-23
# RNA-Seq analysis pipeline
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
# 5. sambamba mark duplicates and build bam index
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
    
with os.popen("which sambamba") as path:
    SAMBAMBA = path.read().strip()
    print('PATH sambamba:', SAMBAMBA)


# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
DB_PATH = "/lustre1/chengqiyi_pkuhpc/zhaohn/1.database"

# HG38
# GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index/hg38_only_chromosome.fa"
# STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index"
# ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index/201902-RefSeq_gene.from_UCSC.hg38.rm_XM_XR.sorted.gtf"

# 这里是小鼠的, 没改变量名
# GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_gencode_GRCm38.p6/GRCm38.p6.genome.fa"
# STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_gencode_GRCm38.p6/star_index_150bp"
# ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_annotation/genome_gencode_GRCm38.p6/gencode.vM25.annotation.gtf"


# 这里是拟南芥的, 没改变量名
GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_ensemblgenomes_tair10.28/genome.fa"
STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_ensemblgenomes_tair10.28/star_index_50bp"
ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_annotation/genome_ensemblgenomes_tair10.28/Arabidopsis_thaliana.TAIR10.28.gtf"

THREADS = '20'

# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
    'Col0-1',
    'Col0-2',
    'sgs3-1',
    'sgs3-2'
]


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        expand("../fix.fastq/{sample}_R1_cutadapt.fastq.gz",sample=SAMPLES),
        expand("../bam/{sample}_Aligned.out.bam",sample=SAMPLES),
        expand("../bam/{sample}_Aligned.out.fix_RG.bam",sample=SAMPLES),
        expand("../bam/{sample}_Aligned_sort.bam",sample=SAMPLES),
        expand("../bam/{sample}_Aligned_sort_rmdup.bam",sample=SAMPLES),
        expand("../bam/{sample}_Aligned_sort_rmdup.bam.bai",sample=SAMPLES),
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# cutadapter
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule TruSeq_cutadapt:
    input:
        "../fastq/{sample}_R1.fastq.gz",
    output:
        "../fix.fastq/{sample}_R1_cutadapt.fastq.gz",
    log:
        "../fix.fastq/{sample}_cutadapt.log"
    shell:# using illumina universal adaptor
        """
        {CUTADAPT} -j {THREADS} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -o {output} {input} > {log} 2>&1
        """
        # -m指定最小长度
# rule NexteraSeq_cutadapt:
#     input:
        # "../fastq/{sample}_R1.fastq.gz",
#     output:
        # "../fix.fastq/{sample}_R1_cutadapt.fastq.gz",
#     log:
        # "../fix.fastq/{sample}_cutadapt.log"
#     shell:# using illumina NexteraSeq adaptor
#         """
#         {CUTADAPT} -j {THREADS} --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
#         -m 55 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
#         -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
#         -o {output[0]} -p {output} {input} > {log} 2>&1
#         """
        # -m指定最小长度


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# STAR mapping
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule STAR_mapping:
    input:
        fq1 = "../fix.fastq/{sample}_R1_cutadapt.fastq.gz",
    output:
        "../bam/{sample}_Aligned.out.bam"
    log:
        "../bam/{sample}_Aligned.out.log"
    params:
        "../bam/{sample}_"
    shell:
        """
        {STAR} \
        --genomeDir {STAR_INDEX} \
        --runThreadN {THREADS} \
        --readFilesIn {input.fq1} \
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
        "../bam/{sample}_Aligned.out.fix_RG.bam"
    params:
        tag = "'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'"
    shell:
        """
        {SAMTOOLS} addreplacerg -r {params.tag} -@ {THREADS} -O BAM -o {output} --reference {GENOME} {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position(not sort by name)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:
        "../bam/{sample}_Aligned.out.fix_RG.bam"
    output:
        "../bam/{sample}_Aligned_sort.bam"
    shell:
        """
        {SAMTOOLS} sort \
        -O BAM \
        -o {output} \
        -T {output}.temp \
        -@ {THREADS}  \
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
                --nthreads={THREADS} \
                --show-progress \
                --sort-buffer-size 8192 \
                {input} {output[0]} > {log} 2>&1
        """