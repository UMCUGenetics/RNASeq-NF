#!/usr/bin/env Rscript

library(edgeR)
# Specifies the raw read_counts (first argument)  
# Specifies the output file name (second argument)

args = commandArgs(trailingOnly=TRUE)

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}


fc.df <- read.delim(args[1], comment.char="#")
# Normalize 
counts <- data.frame(fc.df[,8:ncol(fc.df)])
colnames(counts) <- colnames(fc.df)[8:ncol(fc.df)]
rownames(counts) <- fc.df$Geneid
x <- DGEList(counts=counts, genes=fc.df[,c("Geneid","Length")] )
# RPKM/CPM normalization
df.rpkm <- rpkm(x,x$genes$Length, normalized.lib.sizes=F)
df.cpm <- cpm(x)
df.tpm = tpm3(counts, x$genes$Length)

write.table(df.rpkm, file=paste0(args[2],"_featureCounts_RPKM.txt"), row.names=T, quote=F,sep="   ")
write.table(df.cpm, file=paste0(args[2],"_featureCounts_CPM.txt"), row.names=T, quote=F,sep="   ")
write.table(df.tpm, file=paste0(args[2],"_featureCounts_TPM.txt"), row.names=T, quote=F,sep="   ")
