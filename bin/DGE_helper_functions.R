library("dplyr")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("ggrepel")


#library("gplots")
#library('ReportingTools')
#library("GO.db")
#library("tibble")

#library("tidyverse")

#initializes for use of specific annotation database. Returns a list object with anno_db, gene_id_column and gene_name_column objects.
init_annotation_db <- function (use_biomart=T, species=NULL){
  if(use_biomart){
    gene_id_column = "ensembl_gene_id"
    gene_name_column = "external_gene_name"
    
    library("biomaRt")
    
    #listDatasets(useMart("ensembl"))
    #listEnsembl("GRCh=37", mart=useMart("ensembl"))
    if(species=="human_grch37"){
      grch37 = useEnsembl(biomart="ensembl",GRCh=37)
      anno_db = useDataset("hsapiens_gene_ensembl", grch37)
    }
    else if(species=="human_grch38"){
      anno_db = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    }
    else if(species=="mouse_grcm38") {
      mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="http://aug2020.archive.ensembl.org")
      anno_db = useDataset("mmusculus_gene_ensembl", mart)
    }
    else if(species=="mouse_grcm39"){
      anno_db = useDataset("mmusculus_gene_ensembl", useMart("ensembl"))  
    }
    else if(species=="dog_canfam3.1"){
      anno_db = useDataset("clfamiliaris_gene_ensembl", useMart("ensembl"))
    }
    else if(species=="cat_gFelis_catus_9.0"){
      anno_db = useDataset("fcatus_gene_ensembl", useMart("ensembl"))
    }
  } else {
    
    gene_id_column = "ENSEMBL"
    gene_name_column= "SYMOL"
    
    if(species=="human_grch37"){
      library("org.Hs.eg.db")
      anno_db <-org.Hs.eg.db
    }
    else if(species=="human_grch38"){
      anno_db = NULL
    }
    else if(species=="mouse_grcm39"){
      library("org.Mm.eg.db")
      anno_db <- org.Mm.eg.db 
    }
    else if(species=="dog_canfam3.1"){
      anno_db = NULL
    }
    else if(species=="cat_gFelis_catus_9.0"){
      anno_db = NULL
    }
  }
  
  result <- list("gene_id_column"=gene_id_column, "gene_name_column"=gene_name_column, "anno_db"=anno_db)
  
  return(result)  
}

setCores <- function(nrcores=2){
  #set appropriate multicore using path as indicator for unix/windows
  library("BiocParallel")
  if(.Platform$OS.type == "unix") {
    register(MulticoreParam(nrcores))  
  } else{
    register(SnowParam(nrcores))
  }
}

#Correct sampleids for proper processing in R, samples starting with a digit get an 'X' prepended. Also '-' are replaced with a '.'
correctIds <- function(x){
  x = as.character(x)
  if(is.list(x) || is.vector(x)){
    for(i in 1:length(x)){
      c = substring(x[i], 1, 1)
      if(c>="0" && c<="9"){
        x[i] = paste0("X", x[i])
      }
    }
  }else{
    c = substring(x, 1, 1)
    if(c>="0" && c<="9"){
      x = paste0("X", x)
    }
  }
  gsub("\\-", ".", x)
}

readMetadataFile <- function(metafile= NULL){
  md <-
    data.frame(read.table(
      #paste0(project, "_MetaData.txt"),
      metafile,
      sep = '\t',
      header = T,
      comment.char = ""
    ))
  
  return (md)
}

getMinimalGroupSize <- function(metadata=NULL){
  min_groupsize=100
  groups=unique(metadata$condition)
  for(group in groups){
    min_groupsize=min(min_groupsize, nrow(metadata[metadata$condition==group,]))
  }
  groupsize=max(min_groupsize, 2)
  
  return(groupsize)
}

readCountsFile <- function(countfile = NULL, sep='\t'){
  rawcounts <-
    data.frame(read.table(
      countfile,
      sep = sep,
      header = T,
      row.names = 1
    ))
  return(rawcounts)
}

readResultsFile <- function(resultfile = NULL){
  rawcounts <-
    data.frame(read.table(
      resultfile,
      sep = ',',
      header = T
      
    ))
  return(rawcounts)
}

mergeCountFiles <- function(countfile1=NULL, countfile2=NULL) {
  rawcounts1 <-readCountsFile(countfile1)
  rawcounts2 <-readCountsFile(countfile2)
  
  rawcounts = rawcounts1 %>% rownames_to_column("geneid") %>% full_join(rawcounts2) %>% column_to_rownames("geneid")
  return(rawcounts)
}

getDataForContrast <- function(dds=NULL, contrast=NULL, fcval=1.5, fdrval=0.05, anno_db=NULL) {
  #Generate comparision of interest
  title <- paste0("DESeq2_Volcano_",contrast[1],"_",contrast[2],"_versus_",contrast[3])
  
  #Filter hits according
  #results_contrast are all results for contrast
  results_contrast <- results(dds, contrast = contrast)
  
  #Filter adjusted p-val and L2FC threshold
  results_contrast$filter <- 0
  results_contrast$filter[results_contrast$padj < fdrval & results_contrast$log2FoldChange > log2(fcval)] <- 1
  results_contrast$filter[results_contrast$padj < fdrval & results_contrast$log2FoldChange < -log2(fcval)] <- -1
  results_contrast$filter <- as.factor(results_contrast$filter)
  #filtered only contains the results matching tresholds
  filtered <- results_contrast
  
  #Annotate gene names
  if (!is.null(anno_db)) {
    select <- AnnotationDbi::select
    
    genesymbols <- select(anno_db,  rownames(rawcounts), c(gene_id_column, gene_name_column), gene_id_column)
    symbols <- genesymbols[!duplicated(genesymbols[[gene_name_column]]),]
    symbols.subset <- symbols[symbols[[gene_id_column]] %in% rownames(filtered),]
    rownames(symbols.subset) <- symbols.subset[[gene_id_column]]
    symbols.subset.2 <- subset(symbols.subset,select=c(gene_name_column))
    final <- merge(symbols.subset.2 ,as(filtered, "data.frame"), by='row.names', all=T)
    rownames(final) <- final[,1]
    # Write Filtered genes to output file
  } else {
    final <- filtered
  } 
  final <- final[order(final$padj),]
  
  #write only genes passing filter
  filtered = final[final$filter!=0,]
  filtered <- subset(filtered, select=-c( `filter`))
  write.csv(data.frame("GENEID"=rownames(filtered),filtered), file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],"_padj_",fdrval,".csv"), quote = F, row.names = F)
  
  #write full set
  final <- subset(final, select=-c(`filter`))
  write.csv(data.frame("GENEID"=rownames(final),final), file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],".csv"), quote = F, row.names = F)
  
  return(results_contrast) 
}

plotVolcano <- function(results_contrast=NULL, contrast=NULL, fcval=1.5, fdrval=0.05){
  #
  # -------- QC figures --------------#
  #pdf(paste0("Volcano_",contrast[1],"_",contrast[2],"_versus_",contrast[3],".pdf"))
  par(mfrow=c(1,1))
  # 1) Volcano plot displaying up/down regulated genes with padj_cutoff, l2fc_cutoff
  with(results_contrast, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste0("DE_Volcano_",contrast[1],"_",contrast[2],"_versus_",contrast[3]),xlim = c(-6, 6)))
  with(subset(results_contrast, padj<fdrval & log2FoldChange < -log2(fcval)), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(results_contrast, padj<fdrval & log2FoldChange > log2(fcval)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  abline(v = -log2(fcval), col="blue", lwd=3, lty=2)
  abline(v = log2(fcval), col="red", lwd=3, lty=2)
  legend("bottomright", xjust=1, yjust=1, legend=c("Up","Down"), pch=20, col=c("red","blue"))
}

plotMA <- function(results_contrast=NULL, contrast=NULL,fdrval=0.05){
  # MA Plot showing 
  DESeq2::plotMA(results_contrast, ylim = c(-1,1), alpha = fdrval, main=paste0("DE_MA_",contrast[1],"_",contrast[2],"_versus_",contrast[3]))  
}

#Method used in previous versions of DGE analysis scripts, kept for completeness
generateReportObsolete <- function(dds=NULL, contrast=NULL, fcval=1.5, fdrval=0.05, anno_db=NULL, gene_names=NULL) {
  
  #Generate comparision of interest
  title <- paste0("DESeq2_Volcano_",contrast[1],"_",contrast[2],"_versus_",contrast[3])
  
  #Filter hits according
  res <- results(dds, contrast = contrast)
  
  #Filter adjusted p-val and L2FC threshold
  res$filter <- 0
  res$filter[res$padj < fdrval & res$log2FoldChange > log2(fcval)] <- 1
  res$filter[res$padj < fdrval & res$log2FoldChange < -log2(fcval)] <- -1
  res$filter <- as.factor(res$filter)
  filtered <- res
  #Annotate gene names
  if (!is.null(anno_db)) {
    select <- AnnotationDbi::select    
    
    if(FALSE){
      genesymbols <- select(anno_db, gene_names, "SYMBOL", "ENSEMBL")
      symbols <- genesymbols[!duplicated(genesymbols$ENSEMBL),]
      colnames(symbols) <- c("ENSEMBL","SYMBOL")
      symbols.subset <- symbols[symbols$ENSEMBL %in% rownames(filtered),]
      rownames(symbols.subset) <- symbols.subset$ENSEMBL
      symbols.subset.2 <- subset(symbols.subset,select=c("SYMBOL"))
      # Merge annotation with res data frame 
    }else{
      #genesymbols <- select(anno_db, gene_names, gene_name_column, gene_id_column)
      genesymbols <- select(anno_db,  rownames(rawcounts), c(gene_id_column, gene_name_column), gene_id_column)
      symbols <- genesymbols[!duplicated(genesymbols[[gene_name_column]]),]
      #colnames(symbols) <- c(gene_id_column, gene_name_column)
      symbols.subset <- symbols[symbols[[gene_id_column]] %in% rownames(filtered),]
      rownames(symbols.subset) <- symbols.subset[[gene_id_column]]
      symbols.subset.2 <- subset(symbols.subset,select=c(gene_name_column))
    }
    
    #final <- merge(symbols.subset.2 ,filtered, by='row.names', all=T)
    final <- merge(symbols.subset.2 ,as(filtered, "data.frame"), by='row.names', all=T)
    rownames(final) <- final[,1]
    # Write Filtered genes to output file
  } else {
    final <- filtered
  } 
  final <- final[order(final$padj),]
  
  
  #write only genes passing filter
  #write.csv(final[final$filter!=0,], file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],"_padj_",fdrval,".csv"), quote = F, row.names = T)
  filtered = final[final$filter!=0,]
  
  #filtered <- subset(filtered, select=-c(`Row.names`, `filter`))
  filtered <- subset(filtered, select=-c( `filter`))
  
  #filtered <- subset(filtered, select=-c(`row.names`, `filter`))
  write.csv(data.frame("GENEID"=rownames(filtered),filtered), file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],"_padj_",fdrval,".csv"), quote = F, row.names = F)
  #write full set
  #write.csv(final, file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],".csv"), quote = F, row.names = T)
  #final <- subset(final, select=-c(`Row.names`, `filter`))
  final <- subset(final, select=-c(`filter`))
  write.csv(data.frame("GENEID"=rownames(final),final), file = paste0(contrast[1],"_",contrast[2],"_vs_",contrast[3],".csv"), quote = F, row.names = F)
  
  #
  # -------- QC figures --------------#
  pdf(paste0("Volcano_",contrast[1],"_",contrast[2],"_versus_",contrast[3],".pdf"))
  par(mfrow=c(1,1))
  # 1) Volcano plot displaying up/down regulated genes with padj_cutoff, l2fc_cutoff
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste0("DE_Volcano_",contrast[1],"_",contrast[2],"_versus_",contrast[3]),xlim = c(-6, 6)))
  with(subset(res, padj<fdrval & log2FoldChange < -log2(fcval)), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(res, padj<fdrval & log2FoldChange > log2(fcval)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  abline(v = -fcval, col="blue", lwd=3, lty=2)
  abline(v = fcval, col="red", lwd=3, lty=2)
  legend("bottomright", xjust=1, yjust=1, legend=c("Up","Down"), pch=20, col=c("red","blue"))
  # MA Plot showing 
  DESeq2::plotMA(res, ylim = c(-1,1), alpha = fdrval, main=paste0("DE_MA_",contrast[1],"_",contrast[2],"_versus_",contrast[3]))
  # Write Filtered genes to output file
}


generatePCARepel <- function(pcadata,
                             color = "condition",
                             shape = NULL,
                             title = NULL
) {
  if(is.null(shape)){
    print(
      ggplot(pcadata, aes(x=PC1, y=PC2, label = name)) + 
        geom_point(aes_string(color = color), size = 3, alpha = .8) + 
        #geom_text(aes(label = name), vjust = 2, size = 2) + 
        geom_text_repel(size=2) +
        ggtitle(paste0("DE-PCA ", title)) + 
        theme(plot.title = element_text(hjust = 0.5, face ="bold"))
    )
  }
  else{
    print(
      ggplot(pcadata, aes(x=PC1, y=PC2, label = name)) + 
        geom_point(aes_string(shape = shape, color = color), size = 3, alpha = .8) + 
        #geom_text(aes(label = name), vjust = 2, size = 2) + 
        geom_text_repel(size=2) +
        ggtitle(paste0("DE-PCA ", title)) + 
        theme(plot.title = element_text(hjust = 0.5, face ="bold"))
    )
  }
}

generateHeatmap <- function(dds, title=NULL, rowlabeling="condition", filename=NA) {
  #Heatmap
  rld <- vst(dds, blind=FALSE)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix), rld[[rowlabeling]], sep = " | ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           main = title,
           filename=filename)
}

read_in_feature_counts_custom <- function(file){
  cnt<- read.table(file, sep="\t", header=T, comment.char="#")
  exclude = c("Chr", "Start", "End", "Strand", "Length", "gene_name")
  for (excl_col in exclude){
    if(excl_col %in% colnames(cnt)){
      cnt<- cnt %>% dplyr::select(-all_of(excl_col))
    }
  }
  return(cnt)
}

add_label_to_columnnames <- function(dataframe, label){
  for(i in 2:length(colnames(dataframe))){
    colnames(dataframe)[i] = paste0(colnames(dataframe)[i],"_", label)
  }
  return(dataframe)
}

merge_counts <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "Geneid")
}

merge_md <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "sampleid")
}


