# /usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
GTF=args[1]
EXON_OUT=args[2]

library(GenomicFeatures)
txdb <- makeTxDbFromGFF(GTF, format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
write.table(exonic.gene.sizes, file=EXON_OUT, quote=FALSE, sep='\t', col.names = FALSE)
