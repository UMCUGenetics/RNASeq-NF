#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
path <- as.character(args[1])
myoutname <-as.character(args[2])
 
print(paste0("############################"))
##Read files names
files <- list.files(path=path, pattern="*.txt")
print(sprintf("## Files to be merged are: ##"))
print(files)
print(paste0("############################"))
 
# using perl to manpulate file names by trimming file extension
labs <- paste("", gsub("\\.txt", "", files, perl=TRUE), sep="")
 
##Load all files to list object, use paste to return the trimpping parts to file name
print(sprintf("######### file read START ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
 
cov <- list()
for (i in labs) {
print(i)
filepath <- file.path(path,paste(i,".txt",sep=""))
cov[[i]] <- read.table(filepath,sep = "\t", header=F, stringsAsFactors=FALSE)
colnames(cov[[i]]) <- c("ENSEMBL_GeneID", i)
}
print(sprintf("######### file read END ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
 
## construct one data frame from list of data.frames using reduce function
print(sprintf("######### merge START ######### %s", 
format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
df <-Reduce(function(x,y) merge(x = x, y = y, by ="ENSEMBL_GeneID"), cov)
 
write.table(df,paste(myoutname, "_counts_merged", ".txt",sep=""), sep="\t", quote= F, row.names = F)
 
print(sprintf("######### MERGE FUNCTION COMPLETE ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
