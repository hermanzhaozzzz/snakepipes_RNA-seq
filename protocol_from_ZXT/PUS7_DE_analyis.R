## 01.library
rm(list =ls() )
setwd("/Users/xiaotingzhang/Desktop/tmp/xt_RNA-seq/NC_vs_CCNB1")
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
## 02.read in data
PUS7_deg <- read.csv("deg_result.csv",header=T)
PUS7_deg_down <- dplyr::filter(PUS7_deg,sig == "down")
PUS7_deg_up <- dplyr::filter(PUS7_deg,sig == "up")

PUS7_deg_down_ID_change = bitr(PUS7_deg_down$X, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
PUS7_deg_up_ID_change = bitr(PUS7_deg_up$X, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")


ego_down <- enrichGO(gene = PUS7_deg_down_ID_change$ENTREZID, 
                   #universe = names(geneList), #背景基因集
                   OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                   ont = "BP", #也可以是 CC  BP  MF中的一种
                   pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                   pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                   qvalueCutoff = 1,
                   readable = TRUE)

ego_up <- enrichGO(gene = PUS7_deg_up_ID_change$ENTREZID, 
                     #universe = names(geneList), #背景基因集
                     OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                     ont = "BP", #也可以是 CC  BP  MF中的一种
                     pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                     pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                     qvalueCutoff = 1,
                     readable = TRUE)
##simple
# egoUP_2 <- simplify(ego_up, cutoff=0.7, by="p.adjust",select_fun=min)
# barplot(egoUP_2, showCategory=50,title="EnrichmentGO_UP")

dotplot(ego_up, showCategory=30,title="EnrichmentGO_UP")
dotplot(ego_down, showCategory=30,title="EnrichmentGO_down")

write.csv(ego_down,"ego_down_0417.txt")
write.csv(ego_up,"ego_up_0417.txt")

##KEGG
kk_up <- enrichKEGG(gene = PUS7_deg_up_ID_change$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
barplot(kk_up,title="Enrichment KEGG UP")


kk_down <- enrichKEGG(gene = PUS7_deg_down_ID_change$ENTREZID,
                    organism = 'hsa', #KEGG可以用organism = 'hsa'
                    pvalueCutoff = 1)
barplot(kk_down,title="Enrichment KEGG down")
