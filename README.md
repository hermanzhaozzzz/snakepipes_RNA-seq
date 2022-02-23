1. https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc
先跑质控流程，主要是看看需不需要trim 5‘ 端的不稳定序列和3'端质量不够高的reads
2. 跑mapping流程，如果是star mapping，就用这个流程即可
3. 执行Snakefile进行 cutadapt、Star mapping、bam sort等处理
4. 后续1，可以选择GATK SNP calling和差异表达分析
5. 后续2，差异表达基因
    - 可以使用Snakefile.fc.deseq2.py计算差异表达基因的featureCounts
    - 并且使用format_all.featureCounts.form.ipynb整理表格数据
    - 最后使用DESeq2_plots.ipynb绘图