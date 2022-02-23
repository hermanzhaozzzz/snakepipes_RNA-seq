with os.popen("which featureCounts") as path:
    FEATURECOUNTS = path.read().strip()
    print('PATH cutadapt:', FEATURECOUNTS)
    
    

# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
DB_PATH = "/lustre1/chengqiyi_pkuhpc/zhaohn/1.database"

# HG38
GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index/hg38_only_chromosome.fa"
STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index"
ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.star_index/201902-RefSeq_gene.from_UCSC.hg38.rm_XM_XR.sorted.gtf"

# 这里是小鼠的, 没改变量名
# GENOME = f"{DB_PATH}/db_genomes/genome_fa/genome_gencode_GRCm38.p6/GRCm38.p6.genome.fa"
# STAR_INDEX = f"{DB_PATH}/db_genomes/genome_fa/genome_gencode_GRCm38.p6/star_index_150bp"
# ANNOTATION_GTF = f"{DB_PATH}/db_genomes/genome_annotation/genome_gencode_GRCm38.p6/gencode.vM25.annotation.gtf"


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
    'KD_NC-rep1',
    'KD_NC-rep2',
    'KD_NC-rep3',
    'KD_CCNB1-rep1',
    'KD_CCNB1-rep2',
    'KD_CCNB1-rep3',
    'KD_PUS7-rep1',
    'KD_PUS7-rep2',
    'KD_PUS7-rep3'
]



# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
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

# -T 使用的线程数
# -p 如果是paird end 就用
# -t 将exon作为一个feature
# -g 将gene_id作为一个feature
# -a 参考的gtf/gff
# -o 输出文件
# 最后加上bam文件，有几个就加几个
