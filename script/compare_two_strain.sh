#!/bin/bash
## Synopsis: bash compare_two_strain.sh $ASM1<FASTA> $ASM2<FASTA>
## Requirements: minimap2, last.
## It compare two strains, output summary of the mutation rate (in,del,mism) with the length of the alignments (\sqrt(len * len)).
set -uex

ASM1=$1
ASM2=$2
TEMP=$RANDOM

lastdb "$TEMP"_db "$ASM1" 
last-train "$TEMP"_db "$ASM2" > "$TEMP".par
lastal -p "$TEMP".par "$TEMP"_db "$ASM2" | last-split | last-split -r | tee temp.maf | cargo run --release --bin summarize_maf -- 2000
rm "$TEMP"_db* "$TEMP".par