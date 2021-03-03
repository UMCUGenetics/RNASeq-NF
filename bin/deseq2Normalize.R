#!/usr/bin/env Rscript
library("DESeq2")
# Specifies the raw read_counts (first argument)  
# Specifies the output file name (second argument)

args = commandArgs(trailingOnly=TRUE)
counts_raw <- read.delim(args[1], comment.char="#", row.names = 1)
cds = as.matrix(counts_raw[,7:ncol(counts_raw)])
colnames(cds) <- colnames(counts_raw)[7:ncol(counts_raw)]
#Get column data and set dummy condition
(colData <- data.frame(row.names=colnames(cds),condition=factor(c(rep("exp", ncol(cds))))))
#Create DEseq2 dataset
dds <- DESeqDataSetFromMatrix(countData = cds[, rownames(colData)], colData = colData, design = ~ 1)
#Estimate size factors
size.factors <- estimateSizeFactors(dds)
#Normalize counts
cds.norm <- counts(size.factors, normalized = T)

write.table(
  cds.norm,
  file = paste0(args[2], "_featureCounts_deseq2.txt"),
  row.names = T,
  quote = F,
  sep = "   "
)
