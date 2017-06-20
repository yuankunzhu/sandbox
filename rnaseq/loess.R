dir = '~/Documents/chopwork/mateusz/htg/6_samples/Chop/RSEM'
setwd(dir)
filelist = list.files(pattern = "genes")
bgi <- readRDS('/Users/zhuy/Documents/chopwork/mateusz/htg/NanoString-HTG-Illumina-miRNA/jim/Nanostring/BGI_TPM_Normalized5.rds')
colClasses = c("character", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"); colClasses[5] = "numeric"
datalist = lapply(filelist, function(x)read.table(x, header=T, row.names = 1, colClasses = colClasses))
exp = do.call("cbind", datalist)
logexp = voom(log(do.call("cbind", datalist)+1,2))
logexp2 = voom(do.call("cbind", datalist))
logexp = voom(log(do.call("cbind", datalist)+1,2), normalize.method=)
gene = "ALK"; ensembl = map[ which(map$Symbol==gene), ]$Ensembl; logexp$E[ensembl,]

rsemDir <- '~/Documents/chopwork/mateusz/htg/6_samples/Chop/RSEM'
setwd(rsemDir)
filelist = list.files(pattern = "genes")
do.call(rbind, strsplit(filelist, "\\."))[,1]
datalist2 = lapply(
  filelist, function(x) read.table(
    x, header=T, row.names = 1, colClasses = colClasses
  )
)
exp2 = do.call("cbind", datalist2)
colnames(exp2) <- do.call(rbind, strsplit(filelist, "\\."))[,1]

NormLoess(log2(1+exp2[rowSums(exp2) > 0, , drop=FALSE]))['ENSG00000183454',]