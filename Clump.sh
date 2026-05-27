#!/bin/bash
awk '$2<=5e-08' ${id} >output1.txt
sed  -i '1i SNP P' output1.txt
plink --bfile /scratch/users/s/h/shifang/ldsc/MAGMA/g1000_eur --clump output1.txt  --clump-r2 0.2 --out T123 --clump-kb 1000 --clump-p1 5e-8 --clump-p2 1e-5 --threads 8
awk '{print $1,$3}' T123.clumped > ${id}_SNP.txt
sed '1d' ${id} |awk '{print $2}' | head -n -2 > output.txt
