## Code used for the colocalization analysis
#Colocalization analysis between hypothyroidism and MASLD
##step1. Create test.R2 
##test.R2
library(data.table)
#library(magrittr) # needs to be run every time you start R and want to use %>%
#library(dplyr)
#library(tidyr)
#library("tibble")
#args <- commandArgs(trailingOnly=TRUE)
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
data1$varbeta <- (data1$SE)^2
data2$varbeta <- (data2$SE)^2
input <- merge(data1, data2, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
library("coloc")
#names(input )[names(input )=="V5"]="rs_id"
result <- coloc.abf(dataset1=list(pvalues=input$pval_eqtl, snp=input$SNP, type="cc", s=N1, N=N2,MAF=input$eaf_eqtl), dataset2=list(pvalues=input$pval_gwas,MAF=input$eaf_gwas, snp=input$SNP, type="quant", N=N4))
#library(dplyr)
dd<-data.frame(t(data.frame(print(result[[1]]))))
write.csv(dd,"123.csv",quote=F,row.names=F)

##step2. Perform the Colocalization analysis 
#!/bin/bash
#
#SBATCH --output=66_mpi.txt
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000

while read -r P1 P2 P3 P4 P5 P6 P7 P8 P9; do
    echo "Processing: $P1 $P2 $P3 $P4"

    OUT_PREFIX="${P1}_${P2}_${P3}"

    #
    START=$((P5 - 500000))
    END=$((P5 + 500000))

    #
    awk -v chr="$P4" -v start="$START" -v end="$END" '$2 == chr && $3 >= start && $3 <= end' \
        "/scratch/users/s/h/shifang/ldsc/data/used/$P2" > 123.txt

    awk -v chr="$P4" -v start="$START" -v end="$END" '$2 == chr && $3 >= start && $3 <= end' \
        "/scratch/users/s/h/shifang/ldsc/data/used/$P3" > 456.txt

    # fine-mapping
    Rscript --vanilla test.R2 123.txt 456.txt "$P6" "$P7" "$P8" "$P9"

    #
    if [[ -f "123.csv" ]]; then
        mv 123.csv "${OUT_PREFIX}.csv"
    else
        echo "Warning: 123.csv not found for $OUT_PREFIX"
    fi

    #
    rm -f 123.txt 456.txt

done <ID


#Colocalization analysis between hypothyroidism and MASLD
##step1. Extract the sc-eQTL from the tenk10K based on the genes
#!/bin/bash
#
#SBATCH --output=66_mpi.txt
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

set -euo pipefail

GENE_LIST="cc"
INPUT_PATTERN="*tsv"
OUTPUT_DIR="eqtl1"
GENE_COL=18

if [[ ! -f "$GENE_LIST" ]]; then
    echo "ERROR: gene list file not found: $GENE_LIST" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
files=( $INPUT_PATTERN )
shopt -u nullglob

if [[ ${#files[@]} -eq 0 ]]; then
    echo "ERROR: no files matched pattern: $INPUT_PATTERN" >&2
    exit 1
fi

TMP_GENE_LIST=$(mktemp)
awk 'NF>0 {gsub(/\r/, "", $1); print $1}' "$GENE_LIST" | sort -u > "$TMP_GENE_LIST"

for f in "${files[@]}"; do
    base=$(basename "$f")
    echo "Processing $base ..."

    LC_ALL=C awk -v genes_file="$TMP_GENE_LIST" -v outdir="$OUTPUT_DIR" -v base="$base" -v col="$GENE_COL" '
        BEGIN {
            FS=OFS="\t"
            while ((getline line < genes_file) > 0) {
                gene[line] = 1
            }
            close(genes_file)
        }
        FNR == 1 {
            header = $0
            next
        }
        ($col in gene) {
            out = outdir "/" $col "_" base ".txt"
            if (!(out in header_written)) {
                print header > out
                header_written[out] = 1
            }
            print >> out
        }
    ' "$f"
done

rm -f "$TMP_GENE_LIST"

echo "Done."
echo "Output files are in: $OUTPUT_DIR"

##step2. Perform the multi-trait colocalization analysis using hyprcoloc
#library(susieR)
library(data.table)
#library(magrittr) # needs to be run every time you start R and want to use %>%
#library(dplyr)
#library(tidyr)
#library("tibble")
#args <- commandArgs(trailingOnly=TRUE)
data21<-fread("/scratch/users/s/h/shifang/ldsc/data/raw/EBV/GCST90809825.h.tsv",header=T)
data21$snpid<-paste(data21$chromosome,"_",data21$base_pair_location)
data22<-fread("/scratch/users/s/h/shifang/ldsc/data/raw/EBV/finngen_R12_J10_ASTHMA_EXMORE",header=T)
#names(data22)[names(data22)=="#chrom"]="chr"
data22$snpid<-paste(data22$chrom,"_",data22$pos)
names(data22)[names(data22)=="rsids"]="rsid"
path<-setwd("/scratch/users/s/h/shifang/ldsc/data/coloc/tenk10k/eqtl1/sel")
fileNames = list.files(path="/scratch/users/s/h/shifang/ldsc/data/coloc/tenk10k/eqtl1/sel",pattern=".txt", full.names = TRUE)
coloc_sum <- c()
for (j in c(1:length(fileNames))){
 data2<-data21
data4<-data22
  data1<-fread(fileNames[j],header=T)
#colnames(data1)<-c("CELL_ID","CELL_TYPE","snpid","SNPID","GENE","GENE_ID","CHR","POS","A1","A2","A2_FREQ_ONEK1K","A2_FREQ_HRC","SPEARMANS_RHO","S_STATISTICS","P_VALUE","Q_VALUE","FDR","RSQUARE","GENOTYPED","ROUND")
#data1<-subset(data1,data1$ROUND=="1")
data1$snpid<-paste(data1$CHR,"_",data1$POS)
#names(data2)[names(data2)=="SNP"]="snpid"
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
# 方法A：将列表转换为矩阵
betas <- as.matrix(betas)
ses <- as.matrix(ses)

result <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
library(dplyr)
dd<-data.frame(t(data.frame( result$results)))
dd$ID<-fileNames[j]
if (dim(dd)[1] > 0) {
coloc_sum<-rbind(dd,coloc_sum)
write.csv(coloc_sum,"hycoloc.csv")
}
}
write.csv(coloc_sum,"hycoloc.csv")

