library(DESeq)
library(waveslim)
library(brainwaver)
argument<-commandArgs(trailingOnly = TRUE)

FileName<-substr(argument[1],11,nchar(argument[1]))  # -filename
fdr<-substr(argument[2],6,nchar(argument[2]))  #-pvalue
Output1<-substr(argument[3],10,nchar(argument[3]))  # -output1
Output2<-substr(argument[4],10,nchar(argument[4]))  # -output2
Condition1<-substr(argument[5],7,nchar(argument[5])) # -con1
Condition2<-substr(argument[6],7,nchar(argument[6])) # -con2
fdr<-as.numeric(fdr)
countsTable<-read.delim(FileName,header=TRUE,stringsAsFactors=TRUE)
rownames(countsTable)<-countsTable$gene
countsTable<-countsTable[,-1]

fields<-colnames(countsTable)
#n1<-length(fields[fields=grep('Con1',fields)]) 
#n2<-length(fields[fields=grep('Con2',fields)])
n1<-length(fields[fields=grep(Condition1,fields)])
n2<-length(fields[fields=grep(Condition2,fields)])

conds<-rep(c("c1", "c2"), times = c(n1,n2))


cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )


# estimateVarianceFunctions deprecated, replaced with estimateDispersions to use the params that are closest to what was used for estimateVarianceFunctions
#if (n1==1 & n2==1) cds <- estimateVarianceFunctions( cds,method="blind" ) else  cds <- estimateVarianceFunctions(cds)
if (n1==1 & n2==1) cds<-estimateDispersions(cds,method="blind",sharingMode="fit-only") else  cds<-estimateDispersions(cds,method="per-condition",sharingMode="fit-only",fitType="local")


res <- nbinomTest( cds, "c1", "c2")

temp<- res[!(is.na(res$pval)), ]
thresh<-compute.FDR(temp$pval,fdr)
thresh
resSig<-temp[temp$pval<=thresh, ]
#write.table(data.frame(temp$id,temp$pval),file=paste("DE_all_WithReplicate",FileName,sep='_'),sep="\t")
write.table(data.frame(temp),file=Output1,sep="\t")

resSigDataFrame=data.frame(resSig$id,resSig$pval)
write.table(paste("FDR=",fdr),file=Output2,row.names = FALSE,col.names = FALSE)
#write.table(resSigDataFrame,file=paste("DE_Significant_WithReplicate",FileName,sep='_'),append = TRUE,sep="\t")
write.table(data.frame(resSig),file=Output2,append = TRUE,sep="\t")


