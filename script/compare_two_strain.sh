#!/bin/bash
## Synopsis: bash compare_two_strain.sh $ASM1<FASTA> $ASM2<FASTA>
## Requirements: minimap2, last.
## It compare two strains, output summary of the mutation rate (in,del,mism) with the length of the alignments (\sqrt(len * len)).
set -uex

ASM1=$1
ASM2=$2
OUTPUT_PREFIX=$3
TEMP=$RANDOM

lastdb "$TEMP"_db "$ASM1" 
last-train "$TEMP"_db "$ASM2" > "$TEMP".par
lastal -p "$TEMP".par "$TEMP"_db "$ASM2" | last-split | last-split -r | tee "$OUTPUT_PREFIX".maf | maf-convert tab > "$OUTPUT_PREFIX".tsv
cargo run --release --bin annotate_break_points -- --alignments "$OUTPUT_PREFIX".tsv | tail -n+2 > "$OUTPUT_PREFIX".brk.tsv
cargo run --release --bin summarize_maf -- 2000 < "$OUTPUT_PREFIX".maf |  tr -d '\n'
echo -ne "\t"
tail -n+1 "$OUTPUT_PREFIX".brk.tsv | wc -l 
rm "$TEMP"_db* "$TEMP".par
rm "$OUTPUT_PREFIX".tsv