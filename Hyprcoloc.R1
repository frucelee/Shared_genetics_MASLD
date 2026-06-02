#!/bin/bash
library(data.table)
data21<-fread("/scratch/users/s/h/shifang/data/GCST90728570_modified.tsv",header=T)
data21$snpid<-paste(data21$chromosome,"_",data21$base_pair_location)
data22<-fread("/scratch/users/s/h/shifang/data/GCST90627749_modified.tsv",header=T)
data22$snpid<-paste(data22$chrom,"_",data22$pos)
names(data22)[names(data22)=="rsids"]="rsid"
path<-setwd("/scratch/users/s/h/shifang/coloc/tenk10k/tenk10K_eQTL_sel")
fileNames = list.files(path="/scratch/users/s/h/shifang/coloc/tenk10k/tenk10K_eQTL_sel",pattern=".txt", full.names = TRUE)
coloc_sum <- c()
for (j in c(1:length(fileNames))){
 data2<-data21
data4<-data22
  data1<-fread(fileNames[j],header=T)
data1$snpid<-paste(data1$CHR,"_",data1$POS)
bb<-intersect(data1$snpid,data2$snpid)
data3<-data2[data2$snpid%in%bb,]
bb<-intersect(data3$snpid,data4$snpid)
data5<-data4[data4$snpid%in%bb,]
data3<-data3[,c("snpid","beta","standard_error","rsid")]
data5<-data5[,c("snpid","beta","sebeta","rsid")]
data1<-data1[,c("snpid","BETA","SE")]
data3<-subset(data3,standard_error>0&standard_error<1)
data5<-subset(data5,sebeta>0&sebeta<1)
data1<-subset(data1,SE>0&SE<1)
input <- merge(data3, data5, by="snpid", all=FALSE, suffixes=c("_eqtl","_gwas"))
input1 <- merge(input, data1, by="snpid", all=FALSE, suffixes=c("_eqtl","_gwas"))
input1=input1[!duplicated(input1$rsid_eqtl),]
rownames(input1)<-input1$rsid_eqtl
betas<-input1[,c(2,5,8)]
ses<-input1[,c(3,6,9)]
colnames(betas)<- c("GCST90728570","GCST90627749",j)
colnames(ses)<- c("GCST90728570","GCST90627749",j)
betas <- as.data.frame(betas)
ses <- as.data.frame(ses)
rownames(betas)<-input1$rsid_eqtl
rownames(ses)<-input1$rsid_eqtl

library(hyprcoloc)
traits <- c("GCST90728570","GCST90627749",j)
rsid <- input1$rsid_eqtl
betas <- as.matrix(betas)
ses <- as.matrix(ses)
result <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
result<-data.frame(data.frame(result$results))
if (!is.null(result) && nrow(result) > 0) {
    result$ID <- fileNames[j]
coloc_sum <- bind_rows(coloc_sum, result)
write.csv(coloc_sum,"hycoloc.csv")
  }
}
