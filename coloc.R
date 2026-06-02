#!/usr/bin/env Rscript
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
data1<-fread(exp_path,header=T)
exp_path <- args[2]
data2<-fread(exp_path,header=T)
library(data.table)
N1 <- as.numeric(args[3])
N2<- as.numeric(args[4])
N3 <- as.numeric(args[5])
N4 <- as.numeric(args[6])
colnames(data1)=c("SNP","hg18","bq","A2","A1", "beta","pval", "se", "eaf", "n", "z" )
colnames(data2)=c("SNP","hg18","bq","A2","A1", "beta","pval", "se", "eaf", "n", "z" )
#data1$varbeta <- (data1$SE)^2
#data2$varbeta <- (data2$SE)^2
input <- merge(data1, data2, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
library("coloc")
result <- coloc.abf(dataset1=list(pvalues=input$pval_eqtl, snp=input$SNP, type="cc", s=N1, N=N2,MAF=input$eaf_eqtl), dataset2=list(pvalues=input$pval_gwas,MAF=input$eaf_gwas, snp=input$SNP, type="quant", N=N4))
final<-data.frame(t(data.frame(print(result[[1]]))))
write.csv(final,"tmp.csv",quote=F,row.names=F)





