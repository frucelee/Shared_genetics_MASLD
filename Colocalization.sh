## Code used for the colocalization analysis
#Colocalization analysis between hypothyroidism and MASLD
##step1. Create test.R2 
##test.R2
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
data1$varbeta <- (data1$SE)^2
data2$varbeta <- (data2$SE)^2
input <- merge(data1, data2, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
library("coloc")
result <- coloc.abf(dataset1=list(pvalues=input$pval_eqtl, snp=input$SNP, type="cc", s=N1, N=N2,MAF=input$eaf_eqtl), dataset2=list(pvalues=input$pval_gwas,MAF=input$eaf_gwas, snp=input$SNP, type="quant", N=N4))
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

