CUFFDIFF = "/home/zhaohuanan/anaconda3/envs/snakepipes_cutadapt-STARmapping-FPKM-sortBAM/bin/cuffdiff"
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
HG38_GTF = "/home/zhaohuanan/zhaohn_HD/2.database/annotate_hg38/20200714_ComprehensiveGeneAnnotation_Chr_gencode.v29.annotation.gtf"

# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
    'BE4-0706-rep1'
]
CTRL = [
    '2Mock-1_combined'
]

# SAMPLES = [
#     '2Mock-1_combined',# 前面试bg后面是exp
#     '2BE4-All-1_combined',
#     '2Vector-1_combined',
#     'BE3-1_combined',
#     'BE3-2_combined',
#     'BE4-1_combined',
#     'BE4-2_combined',
#     'EM-1_combined',
#     'EM-2_combined',
#     'BE4-0706-rep1',
#     'M2-0706-rep1'
# ]

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        expand("../bam/293T-RNASeq-{sample}_Aligned_sort.bam",sample=SAMPLES),
        expand("../bam/293T-RNASeq-{ctrl}_Aligned_sort.bam",ctrl=CTRL),
        expand("../cuffdiff/{sample}_vs_{ctrl}",sample=SAMPLES,ctrl=CTRL)


rule cuffdiff:
    # 摘要：如何解决cuffdiff运行结果中表达量为0的情况？ cuffdiff -o cdiffout -b ref.fasta -u -p 15 --library-type fr-firststrand \ -L H_3,O_3,F_3 cuffmergeout/merged.gtf \ tophato 阅读全文
    input:
        "../bam/293T-RNASeq-{sample}_Aligned_sort.bam",
        "../bam/293T-RNASeq-{ctrl}_Aligned_sort.bam"
    output:
        directory("../cuffdiff/{sample}_vs_{ctrl}")
#     params:
#         l1 = '{sample}',
#         l2 = '{ctrl}'
    log:
        "../cuffdiff/{sample}_vs_{ctrl}.log"
    shell:
        """
        {CUFFDIFF} \
        -o {output} \
        -b {HG38_FA} \
        -p 24 \
        --FDR 0.01 \
        -u {HG38_GTF} \
        {input[0]} \
        {input[1]} \
        --no-update-check 2>{log}
        """
