#!/usr/bin/env Rscript

library(edgeR)
# Specifies the raw read_counts (first argument)  
# Specifies the exon_gene_sizes (second argument)
# Specifies the output file name (third argument)

args = commandArgs(trailingOnly=TRUE)

outname <- args[1]
raw_read_counts <- read.table(args[2],sep="\t",header=T,row.names=1)
exon_gene_sizes <- read.table(args[3],sep="\t",header=F,row.names=1)


nrsamples <- ncol(raw_read_counts)
nrrows <- nrow(raw_read_counts)

tab <- matrix(data=NA, nrow=nrrows, ncol=nrsamples)

for (j in 1:nrsamples){
  RPKM = rpkm(raw_read_counts[j], exon_gene_sizes, normalized.lib.sizes=F, log=F)
  for(i in 1:nrrows){
    tab[i,j] = RPKM$V2[i]
  }
}

df <- data.frame(tab)
colnames(df) <- colnames(raw_read_counts)
rownames(df) <- rownames(raw_read_counts)
outfile=paste0(args[1],"_readCounts_RPKM.txt")
write.table(df, file=outfile, row.names=T,col.names=NA, quote=F,sep="   ")
