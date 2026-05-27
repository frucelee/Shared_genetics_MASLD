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
