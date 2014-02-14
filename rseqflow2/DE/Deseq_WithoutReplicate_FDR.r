library(DESeq)
library(waveslim)
library(brainwaver)
fisher <- function(x)
{
   x <- x[!is.na(x)]
   q <- -2*sum(log(x))
   d <- 2*length(x)
   1-pchisq(q, df=d)
}

argument<-commandArgs(trailingOnly = TRUE)
FileName<-substr(argument[1],11,nchar(argument[1]))  # -filename
fdr<-substr(argument[2],6,nchar(argument[2]))  #-fdr
Output1<-substr(argument[3],10,nchar(argument[3]))  # -output1
Output2<-substr(argument[4],10,nchar(argument[4]))  # -output2

fdr<-as.numeric(fdr)
countsTable<-read.delim(FileName,header=TRUE,stringsAsFactors=TRUE)
countsRead<-countsTable[,3:4]
rownames(countsRead)<-paste(countsTable$gene,countsTable$exon,sep = ":")
conds<-rep(c("c1", "c2"))
cds <- newCountDataSet( countsRead, conds )
cds <- estimateSizeFactors( cds )
sizeFactors(cds)

# estimateVarianceFunctions deprecated, replaced with estimateDispersions
#cds <- estimateVarianceFunctions( cds,method="blind" ) 
cds <- estimateDispersions( cds,method="blind",sharingMode="fit-only")
res <- nbinomTest( cds, "c1", "c2")
temp<- res[!(is.na(res$pval)), ]

#resSig<-temp[temp$pval<=pvalue, ]
#resSigDataFrame=data.frame(resSig$id,resSig$pval)
#write.table(resSigDataFrame,file="DE_Output_Significant_exon.txt",sep="\t")



q.value<-data.frame(gene=countsTable[,1],p=rep(0,nrow(countsTable)))
q.value[,2]<- res[,7]

genename<-unique(countsTable[,1])

p.combined<-data.frame(gene=genename,pval=numeric(length(genename)))
for(i in 1:length(genename)){
  trans <- q.value[q.value[,1]==genename[i],]  
  p.combined[i,2]<- fisher(trans[,2])
}

thresh<-compute.FDR(p.combined$pval,fdr)
#resSig<-p.combined[p.combined$pval<=thresh, ]



write.table(p.combined,file=Output1,sep="\t")
write.table(paste("FDR=",fdr),,file=Output2,row.names = FALSE,col.names = FALSE)
write.table(p.combined[p.combined$pval<=thresh, ],file=Output2,,append = TRUE,sep="\t")

