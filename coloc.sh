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
