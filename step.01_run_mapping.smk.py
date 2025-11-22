# snakemake --profile polaris.cnshort -j 10 -s step.01_run_mapping.smk.py --rerun-incomplete --restart-times 0 --keep-going -n
# snakemake --profile polaris.cnshort -j 10 --cores 20 -s step.01_run_mapping.smk.py -n
configfile: "config.yaml"


treatments = config["treatments"]
controls = config["controls"]
samples = treatments + controls

input_dir = config["input_dir"]
output_base = config["output_base"]

# featureCounts 链特异性（0/1/2），以及是否用 rmdup BAM 做计数
FC_STRAND = config.get("featurecounts_strand", 0)
COUNT_USE_RMDUP = config.get("count_use_rmdup", False)

# 标准做法（bulk RNA-seq）是不要去重。重复读数大多反映真实高表达，非 UMI 数据无法可靠区分 PCR duplicate 与生物学重复。
BAMS_FOR_FC = (
    expand(f"{output_base}/bam/{{sample}}_Aligned_sort_rmdup.bam", sample=samples)
    if COUNT_USE_RMDUP
    else
    expand(f"{output_base}/bam/{{sample}}_Aligned_sort.bam", sample=samples)
)


# 读取 mpileup 配置
MPI_CFG = config.get("mpileup", {})
MPI_ENABLE = MPI_CFG.get("enable", False)
MPI_USE_RMDUP = MPI_CFG.get("use_rmdup", True)
MPI_MAX_DEPTH = MPI_CFG.get("max_depth", 10000)
MPI_MIN_MAPQ  = MPI_CFG.get("min_mapq", 20)
MPI_MIN_BQ = MPI_CFG.get("min_bq", 20)
MPI_MUT_THRES = MPI_CFG.get("mut_threshold", 0)

# 选择做 mutation 用的 BAM（去重 或 未去重）
BAM_FOR_MUT = (
    lambda s: f"{output_base}/bam/{s}_Aligned_sort_rmdup.bam" if MPI_USE_RMDUP
              else f"{output_base}/bam/{s}_Aligned_sort.bam"
)

# 统一收集 rule all 目标
ALL_TARGETS = [
    # RNA-seq mapping + raw counts：
    *expand(f"{output_base}/bam/{{sample}}_Aligned_sort_rmdup.bam", sample=samples),
    *expand(f"{output_base}/bam/{{sample}}_Aligned_sort_rmdup.bam.bai", sample=samples),
    f"{output_base}/featureCounts/all_feature.txt",
]

# 如果开启 mpileup，则追加 mpileup 及其解析结果到 all 目标
if MPI_ENABLE:
    ALL_TARGETS += [
        *expand(f"{output_base}/mut_tables/{{sample}}.base_mut.tsv.gz", sample=samples),
        *expand(f"{output_base}/mut_tables/{{sample}}.base_mut.filtered.tsv", sample=samples),
    ]
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# rule all
# ------------------------------------------------------------------------------------------>>>>>>>>>>
localrules: filter_mut_table
rule all:
    input: ALL_TARGETS
# ------------------------------------------------------------------->>>>>>>>>>
# trim adaptor
# ------------------------------------------------------------------->>>>>>>>>>
rule fastp_trim_adaptor:
    input:
        fwd=f"{input_dir}/{{sample}}_R1.fastq.gz",
        rev=f"{input_dir}/{{sample}}_R2.fastq.gz"
    output:
        fwd=temp(f"{output_base}/fix.fastq/{{sample}}_R1_cutadapt.fastq.gz"),
        rev=temp(f"{output_base}/fix.fastq/{{sample}}_R2_cutadapt.fastq.gz"),
        html=f"{output_base}/fix.fastq/{{sample}}.html",
        json=f"{output_base}/fix.fastq/{{sample}}.json"
    threads: config["threads"]
    log: f"{output_base}/fix.fastq/{{sample}}.log"
    shell:
        r"""
        fastp -w {threads} -h {output[html]} -j {output[json]} \
              -i {input[fwd]} -I {input[rev]} -o {output[fwd]} -O {output[rev]} 2> {log}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# STAR mapping
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule STAR_mapping:
    input:
        fq1 = f"{output_base}/fix.fastq/{{sample}}_R1_cutadapt.fastq.gz",
        fq2 = f"{output_base}/fix.fastq/{{sample}}_R2_cutadapt.fastq.gz"
    output:
        temp(f"{output_base}/bam/{{sample}}_Aligned.out.bam")
    threads: config["threads"]
    log: f"{output_base}/bam/{{sample}}_Aligned.out.log"
    params:
        prefix = f"{output_base}/bam/{{sample}}_",
        star_idx = config["star_index"]
    shell:
        r"""
        STAR \
          --genomeDir {params.star_idx} \
          --runThreadN {threads} \
          --readFilesIn {input.fq1} {input.fq2} \
          --readFilesCommand zcat \
          --outFileNamePrefix {params.prefix} \
          --outSAMtype BAM Unsorted \
          --outSAMstrandField intronMotif \
          --outFilterIntronMotifs RemoveNoncanonical > {log} 2>&1
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# add @RG tag (mostly for GATK SNP/SNV calling)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule add_RG_tag:
    input:  f"{output_base}/bam/{{sample}}_Aligned.out.bam"
    output: temp(f"{output_base}/bam/{{sample}}_Aligned.out.fix_RG.bam")
    params:
        tag=lambda w: f"'@RG\\tID:{w.sample}\\tSM:{w.sample}\\tPL:ILLUMINA'",
        genome = config["genome"]
    threads: config["threads"]
    shell:
        "samtools addreplacerg -r {params.tag} -@ {threads} -O BAM -o {output} --reference {params.genome} {input}"

rule filter_bam:
    input:  f"{output_base}/bam/{{sample}}_Aligned.out.fix_RG.bam"
    output: temp(f"{output_base}/bam/{{sample}}_Aligned.out.fix_RG_filter.bam")
    threads: config["threads"]
    shell:
        "samtools view -@ {threads} -F 4 -F 8 -hb {input} -o {output}"
# 12
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position(not sort by name)
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:  f"{output_base}/bam/{{sample}}_Aligned.out.fix_RG_filter.bam"
    output: f"{output_base}/bam/{{sample}}_Aligned_sort.bam"
    threads: config["threads"]
    shell:
        r"""
        samtools sort -O BAM -o {output} -T {output}.temp -@ {threads} {input}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# sambamba rmdup and build bam index
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule sambamba_rmdup_and_build_index:
    input:  f"{output_base}/bam/{{sample}}_Aligned_sort.bam"
    output:
        f"{output_base}/bam/{{sample}}_Aligned_sort_rmdup.bam",
        f"{output_base}/bam/{{sample}}_Aligned_sort_rmdup.bam.bai"
    threads: config["threads"]
    log: f"{output_base}/bam/{{sample}}_Aligned_sort_rmdup.log"
    shell:
        r"""
        sambamba markdup --remove-duplicates --nthreads={threads} \
                         --show-progress --sort-buffer-size 8192 \
                         {input} {output[0]} > {log} 2>&1
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# featureCounts
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# https://www.jianshu.com/p/e20a3b73dcd0
# https://www.jianshu.com/p/9cc4e8657d62
# -T 使用的线程数
# -p 如果是paird end 就用, 只能用在paired-end的情况中，会统计fragment而不统计read
# -B 在-p选择的条件下，只有两端read都比对上的fragment才会被统计
# -t 将exon作为一个feature
# -g 将gene_id作为一个feature
# -a 参考的gtf/gff
# -o 输出文件
# 最后加上bam文件，有几个就加几个

# 如果是 single-end，把 -p -B 去掉

# 标准做法（bulk RNA-seq）是不要去重。重复读数大多反映真实高表达，非 UMI 数据无法可靠区分 PCR duplicate 与生物学重复。
# 只有 UMI 流程（或明确的 PCR 过度扩增异常）才做 UMI 去重；ChIP/ATAC 常做去重，但 RNA-seq 计数通常不做。

rule featureCounts:
    input:  BAMS_FOR_FC
    output: f"{output_base}/featureCounts/all_feature.txt"
    params: anno_gtf = config["annotation_gtf"]
    threads: config["threads"]
    log:    f"{output_base}/featureCounts/run_FC.log"
    shell:
        r"""
        featureCounts \
          -T {threads} \
          -p -B \
          -t exon \
          -g gene_id \
          -s {FC_STRAND} \
          -a {params.anno_gtf} \
          -o {output} \
          {input} 2>{log}
        """
rule spike_in_mpileup:
    input:
        bam = lambda wc: BAM_FOR_MUT(wc.sample),
        bai = lambda wc: BAM_FOR_MUT(wc.sample) + ".bai",
        ref = config["genome"]
    output:
        mpileup = temp(f"{output_base}/mpileup/{{sample}}.mpileup")
    threads: config["threads"]
    shell:
        r"""
        samtools mpileup {input.bam} \
            --reference {input.ref} \
            --max-depth {MPI_MAX_DEPTH} \
            -q {MPI_MIN_MAPQ} -Q {MPI_MIN_BQ} \
            > {output.mpileup}
        """
rule parse_mpileup:
    input:
        mpileup = rules.spike_in_mpileup.output.mpileup
    output:
        table = f"{output_base}/mut_tables/{{sample}}.base_mut.tsv.gz"
    params:
        temp_dir = f"{output_base}/mut_tables/{{sample}}.tmp",
        mut_thres = MPI_MUT_THRES
    threads: config["threads"]
    shell:
        r"""
        bioat bam mpileup2table {input.mpileup} \
            -o {output.table} \
            --mutation_number_threshold {params.mut_thres} \
            --temp_dir {params.temp_dir} \
            --threads {threads}
        """
rule filter_mut_table:
    input:
        table = rules.parse_mpileup.output.table  # .base_mut.tsv.gz
    output:
        filtered = f"{output_base}/mut_tables/{{sample}}.base_mut.filtered.tsv"
    threads: config["threads"]
    shell:
        r"""
        set -eu
        export LC_ALL=C

        in="{input.table}"
        out="{output.filtered}"

        echo "[INFO] mawk 流式过滤启动" >&2
        echo "[INFO] Input : $in" >&2
        echo "[INFO] Output: $out" >&2

        # 解压(.gz) -> 删除 4 列 -> 按 mut_num>0 过滤 -> 输出明文 TSV
        gzip -dc "$in" | \
        mawk -F'\t' -v OFS='\t' '
            NR==1 {{
                # 选择保留列，并找到 mut_num 列
                for (i = 1; i <= NF; i++) {{
                    if ($i == "mut_num") col_mut = i
                    if ($i != "ambiguous_count" && $i != "deletion" && $i != "insertion" && $i != "ambiguous")
                        keep[++k] = i
                }}
                if (!col_mut) {{
                    print "WARNING: mut_num column not found, skip mut_num filter" > "/dev/stderr"
                    col_mut = -1
                }}
                # 输出表头（去掉 4 个列）
                for (j = 1; j <= k; j++)
                    printf "%s%s", $(keep[j]), (j < k ? OFS : ORS)
                next
            }}
            {{
                # 如果没有 mut_num，只做删列（容错）
                if (col_mut < 0) {{
                    for (j = 1; j <= k; j++) {{
                        field = (keep[j] <= NF ? $(keep[j]) : ".")
                        printf "%s%s", field, (j < k ? OFS : ORS)
                    }}
                    next
                }}

                # 列数异常，直接跳过该行
                if (NF < col_mut) next

                # 读取 mut_num，非正数丢弃；非数字会按 0 处理，也会被丢弃
                val = $(col_mut) + 0
                if (val <= 0) next

                # 正常输出保留列，缺失列补 "."
                for (j = 1; j <= k; j++) {{
                    field = (keep[j] <= NF ? $(keep[j]) : ".")
                    printf "%s%s", field, (j < k ? OFS : ORS)
                }}
            }}
        ' > "$out"

        echo "[INFO] 完成: $out" >&2
        """