#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
path <- as.character(args[1])
out_name <-as.character(args[2])

quant.files <- list.files(path, full.names = T, pattern = '*.sf', recursive = T)

raw_counts <- list()
tpm_counts <- list()

for (sample in quant.files){
  sample_id <- basename(dirname(sample))
  data <- read.table(sample, sep = "\t", header=T, stringsAsFactors=FALSE)
  raw_counts[[sample_id]] <- data[,c(1,5)]
  tpm_counts[[sample_id]] <- data[,c(1,4)]
  colnames(raw_counts[[sample_id]]) <- c("ENSEMBL_GeneID", sample_id)
  colnames(tpm_counts[[sample_id]]) <- c("ENSEMBL_GeneID", sample_id)
}
# Merge dataframe by ENSEMBL_ID
df.raw <-Reduce(function(x,y) merge(x = x, y = y, by ="ENSEMBL_GeneID"), raw_counts)
df.tpm <-Reduce(function(x,y) merge(x = x, y = y, by ="ENSEMBL_GeneID"), tpm_counts)

# Write merged count tables
write.table(df.raw, paste(out_name, "_Salmon_counts_raw", ".txt",sep=""), sep="\t", quote= F, row.names = F)
write.table(df.tpm, paste(out_name, "_Salmon_counts_tpm", ".txt",sep=""), sep="\t", quote= F, row.names = F)


