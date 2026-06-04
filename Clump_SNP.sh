#!/bin/bash
awk '$2<=5e-08' CPASSOC_gwas_GCST90627749_GCST90728570.tsv >output.txt
sed  -i '1i SNP P' output.txt
plink --bfile /scratch/users/s/h/shifang/ldsc/MAGMA/g1000_eur --clump output.txt  --clump-r2 0.2 --out tmp --clump-kb 1000 --clump-p1 5e-8 --clump-p2 1e-5 --threads 8
awk '{print $1,$3}' tmp.clumped > CPASSOC_GCST90627749_GCST90728570_clup_SNP.txt
sed '1d' CPASSOC_GCST90627749_GCST90728570_clup_SNP.txt |awk '{print $2}' | head -n -2 > CPASSOC_GCST90627749_GCST90728570_SNP_only.txt
rm -r output.txt tmp* 
