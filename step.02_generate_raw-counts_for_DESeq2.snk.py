with os.popen("which featureCounts") as path:
    FEATURECOUNTS = path.read().strip()
    print('PATH cutadapt:', FEATURECOUNTS)
    
    

# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
DB_PATH = "/lustre1/chengqiyi_pkuhpc/zhaohn/1.database"

# HG38
# GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index/hg38_only_chromosome.fa"
# STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index"
# ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_annotation/genome_ucsc_hg38/ucsc_hg38_genes-and-gene-predictions_NCBI-refseq_refFlat.gtf"

# 这里是小鼠的, 没改变量名
GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_gencode_GRCm38.p6/GRCm38.p6.genome.fa"
STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_gencode_GRCm38.p6/star_index_150bp"
# ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_annotation/genome_gencode_GRCm38.p6/gencode.vM25.annotation.gtf"
ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_annotation/genome_gencode_GRCm38.p6/ucsc_mm10_genes-and-gene-predictions_NCBI-refseq_refFlat.gtf"


# 这里是拟南芥的, 没改变量名
# GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_ensemblgenomes_tair10.28/genome.fa"
# STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_ensemblgenomes_tair10.28/star_index_50bp"
# ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_annotation/genome_ensemblgenomes_tair10.28/Arabidopsis_thaliana.TAIR10.28.gtf"

THREADS = '20'
# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>


# 前面WT后面是TREAT
# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
    'LI-rep1',
    'LI-rep2',
    'mLN-rep1',
    'mLN-rep2',
    'pLN-rep1',
    'pLN-rep2',
    'SI-rep1',
    'SI-rep2',
    'spleen-rep1',
    'spleen-rep2',
]




# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# https://www.jianshu.com/p/9cc4e8657d62
# -T 使用的线程数
# -p 如果是paird end 就用, 只能用在paired-end的情况中，会统计fragment而不统计read
# -B 在-p选择的条件下，只有两端read都比对上的fragment才会被统计
# -t 将exon作为一个feature
# -g 将gene_id作为一个feature
# -a 参考的gtf/gff
# -o 输出文件
# 最后加上bam文件，有几个就加几个

# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule all:
    input:
        '../featureCounts/all_feature.txt'

rule featureCounts:
    params:
        BAM = ' '.join(expand("../bam/{sample}_Aligned_sort_rmdup.bam",sample=SAMPLES))
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
        {params.BAM} 2>{log}
        """
