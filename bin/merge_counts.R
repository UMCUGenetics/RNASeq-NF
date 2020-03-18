args <- commandArgs(TRUE)
path <- as.character(args[1])
myoutname <-as.character(args[2])
 
print(paste0("############################"))
##Read files names
files <- list.files(path=path, pattern="*.txt")
print(sprintf("## Files to be merged are: ##"))
print(files)
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
filepath <- file.path(path,paste(i,".cnt",sep=""))
cov[[i]] <- read.table(filepath,sep = "\t", header=F, stringsAsFactors=FALSE)
colnames(cov[[i]]) <- c("ENSEMBL_GeneID", i)
}
print(sprintf("######### file read END ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
 
## construct one data frame from list of data.frames using reduce function
print(sprintf("######### merge START ######### %s", 
format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
df <-Reduce(function(x,y) merge(x = x, y = y, by ="ENSEMBL_GeneID"), cov)
print(sprintf("######### merge END ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
 
print(sprintf("Exported merged table within work directory in txt and Rdata format with file name merged_%s_%s", make.names(format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")), myoutname))
 
write.table(df,paste(path, myoutname,".txt",sep=""), sep="\t", quote= F, row.names = F)
 
save.image(paste(path, myoutname,".Rdata",sep=""))
 
print(sprintf("######### MERGE FUNCTION COMPLETE ######### %s", format(Sys.time(),"%b_%d_%Y_%H_%M_%S_%Z")))
