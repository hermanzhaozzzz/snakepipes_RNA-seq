with os.popen("which featureCounts") as path:
    FEATURECOUNTS = path.read().strip()
    print('PATH cutadapt:', FEATURECOUNTS)
    
    
    
# human 
# HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/bwa_hg38/hg38_only_chromosome.fa"
# HG38_GTF = "/home/zhaohuanan/zhaohn_HD/2.database/annotate_hg38/20200714_ComprehensiveGeneAnnotation_Chr_gencode.v29.annotation.gtf"
# mouse
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_gencode_GRCm38.p6/GRCm38.p6.genome.fa"
HG38_GTF = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_annotation/genome_gencode_GRCm38.p6/gencode.vM25.annotation.gtf"
# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>


# 前面WT后面是TREAT
SAMPLES = [
    "pcif1-WT-1",
    "pcif1-WT-2",
    "pcif1-WT-3",
    "pcif1-WT-4",
    "pcif1-KO-1",
    "pcif1-KO-2",
    "pcif1-KO-3",
    "pcif1-KO-4"
]


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        '../featureCounts/all_feature.txt'

rule featureCounts:
    params:
        BAM = ' '.join(expand("../bam/293T-RNASeq-{sample}_Aligned_sort.bam",sample=SAMPLES))
    output:
        '../featureCounts/all_feature.txt'
    log:
        '../featureCounts/run_FC.log'
    shell:
        """
        {FEATURECOUNTS} \
        -T 24 \
        -p \
        -t exon \
        -g gene_id \
        -a {HG38_GTF} \
        -o {output} \
        {params.BAM} 2>{log}
        """

# -T 使用的线程数
# -p 如果是paird end 就用
# -t 将exon作为一个feature
# -g 将gene_id作为一个feature
# -a 参考的gtf/gff
# -o 输出文件
# 最后加上bam文件，有几个就加几个
