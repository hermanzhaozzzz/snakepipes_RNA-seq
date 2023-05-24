rm(list = ls())
#setwd("/Users/xiaotingzhang/Documents/YiLab/data/140_Proj_tsRNASeq/25-tsRNA-rawdata-20220501")
setwd("/Users/zhaohuanan/Desktop/tsRNA")

# import ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(pheatmap)
library("reshape2")
library(vegan)
library(purrrlyr)


# readin tRNA data -------------------------------------------------------------
title <- c("PUS10_R1","PUS10_R2","PUS1_R1","PUS1_R2","PUS3_R1","PUS3_R2",
           "PUS7_R1","PUS7_R2","PUS7L_R1","PUS7L_R2","PUSL1_R1","PUSL1_R2",
           "RD1_R1","RD1_R2","RD2_R1","RD2_R2","RD3_R1","RD3_R2","TRUB1_R1",
           "TRUB1_R2","WT_R1","WT_R2")
data <- read.csv("20220428_tsRNA_15_50_tidy_RPM.tRNA.csv",header=F,row.names = 1)

names(data) <- title

##计算每一列的最大值
data %>% apply(1,max) -> data$max_row
# hist(log2(data$max_row + 1))

data_filter <- dplyr::filter(data,max_row >= 10)

dim(data)
dim(data_filter)
head(data_filter)
##将最后一列max_row去除，这个也是第一次使用mutate
data_filter <- data_filter %>% 
  mutate(max_row=NULL)

head(data_filter)
# pheatmap(t(log2(data_filter+1)),show_colnames = F)


# readin spikein ----------------------------------------------------------

spikein <- read.csv("20220428_tsRNA_15_50_tidy_RPM.spikein.csv",header=F)
names(spikein) <- c("id",title)
head(spikein)
spikein <- dplyr::filter(spikein,(id == "spikein_34nt_1-34" )|(id == "spikein_50nt_1-50" ))
head(spikein)
spikein34 <- dplyr::filter(spikein,id == "spikein_34nt_1-34" )
spikein50 <- dplyr::filter(spikein,id == "spikein_50nt_1-50" )
rownames(spikein34) <- spikein34$id
rownames(spikein50) <- spikein50$id
spikein34 <- spikein34[,-1]
spikein50 <- spikein50[,-1]
# rownames(spikein) <- spikein$id
# spikein <- spikein[,-1]
# 
# head(spikein)
# pheatmap(spikein)
# 
# spikein <- as.data.frame(t(spikein))
# spikein$ratio <- spikein$`spikein_50nt_1-50`/spikein$`spikein_34nt_1-34`


# 1.RPS of 34nt -----------------------------------------------------------

data_filter_spikein34 <- rbind(spikein34,data_filter)
data_filter_spikein34 <- as.data.frame(apply(data_filter_spikein34,2,function(x){x/x[1]}*10^15))
data_filter_spikein34 <- data_filter_spikein34[-1,]
head(data_filter_spikein34)
# 
# pheatmap(t(log2(data_filter_spikein34+1)),
#          show_colnames = F,
#          main = "data_filter_spikein34_normlize log2",
#          clustering_distance_rows = "correlation")

pheatmap(t(log2(data_filter_spikein34+1)),
         show_colnames = F,
         main = "data_filter_spikein34_normlize log2",
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2")

pheatmap(t(log2(data_filter_spikein34+1)),
         show_colnames = F,
         main = "data_filter_spikein34_normlize log2",
         clustering_method = "ward.D2")

pheatmap(t(log2(data_filter_spikein34+1)),
         show_colnames = F,
         main = "data_filter_spikein34_normlize log2",
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D")



pheatmap(t(log2(data_filter_spikein34+1)),
         show_colnames = F,
         main = "data_filter_spikein34_normlize log2",
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2")

pheatmap(t(log2(data_filter_spikein34+1)),
         show_colnames = F,
         main = "data_filter_spikein34_normlize log2")


pheatmap(t(log2(data_filter_spikein34+1)),
         show_colnames = F,
         main = "data_filter_spikein34_normlize log2",
         clustering_distance_rows = "manhattan")

# pheatmap(t(data_filter_spikein34),
#           show_colnames = F,
#          main = "data_filter_spikein34_normlize ")
# 

# 2.RPS for spikein_50nt --------------------------------------------------

data_filter_spikein50 <- rbind(spikein50,data_filter)
data_filter_spikein50 <- as.data.frame(apply(data_filter_spikein50,2,function(x){x/x[1]}*10^6))
data_filter_spikein50 <- data_filter_spikein50[-1,]
head(data_filter_spikein50)

pheatmap(t(log2(data_filter_spikein50+1)),
         show_colnames = F,
         main = "data_filter_spikein50_normlize log2")
pheatmap(t(data_filter_spikein50),
         show_colnames = F,
         main = "data_filter_spikein50_normlize ")
