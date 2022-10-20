#!/bin/bash
set -uex
## Synopsis: <REFERENCE:fasta> <Annotation:GFF3> <PREFIX:Dirname> <CONTIGS:[fasta]>
## Require: minimap2, liftoff, last,
REFERENCE=$1
ANNOTATION=$2
PREFIX=$3
CONTIGS="$4"

mkdir -p "$PREFIX"


### 0. Prepare contigs.
NUM=0
for contig in $CONTIGS
do
    strain=$( basename "$contig" )
    awk --assign ctg="$NUM" 'BEGIN{tig=0} ($0 ~ />/){print(">" ctg "_" tig); tig+=1; next}{print $0}' "$contig" > "$PREFIX"/"$strain"
    NUM=$(( NUM + 1 ))
done
awk '($0 ~ />/){print(">ref");next}{print $0}' "$REFERENCE" > "$PREFIX"/reference.fa

### 1. liftover alignments
function liftover_by_liftoff () {
    GFFFILE=$1
    TARGET_CTG=$2
    REF_CTG=$3
    OUT_GFF=${TARGET_CTG%.fa}.gff3
    OUT_TSV=${TARGET_CTG%.fa}.tsv
    TEMP="$RANDOM"
    liftoff -o "$OUT_GFF" -dir "$TEMP" -g "$GFFFILE" "$TARGET_CTG" "$REF_CTG"
    grep -v "^#" "$OUT_GFF" |\
        awk 'BEGIN{OFS="\t"}($3 ~/gene/){if ($7 == "+"){ dir = "1"} else { dir = "-1"};print($1,$4,$5,dir,"gene")}' > "$OUT_TSV"
    rm -r "$TEMP"
}

awk '($0 ~ /^#/){print $0;next} ($0 ~ /^NC_037304.1/){print $0}' "$ANNOTATION" | sed -e 's/NC_037304.1/ref/g'> "$PREFIX"/reference.gff3
grep -v "^#" "$PREFIX"/reference.gff3 |\
    awk 'BEGIN{OFS="\t"}($3 ~ /gene/){if ($7 == "+"){ dir = "1"} else { dir = "-1"};print($1,$4,$5,dir,"gene")}' > "$PREFIX"/reference.tsv
for contig in $CONTIGS
do
    strain=$( basename "$contig" )    
    TARGET="$PREFIX"/"$strain"
    liftover_by_liftoff "$PREFIX"/reference.gff3 "$TARGET" "$PREFIX"/reference.fa
done

### 2. All vs All alignments
cat "$PREFIX"/*.fa | sed -e 's/>\(.*\)/>db_\1/g'> "$PREFIX"/seqs.fasta
lastdb "$PREFIX"/seqs.db "$PREFIX"/seqs.fasta
for seq in "$PREFIX"/*.fa
do
    outfile=${seq%.fa}.maf
    echo "##maf version=1 scoring=lastz.v1.03.73" > "$outfile"
    lastal -E0.0001 "$PREFIX"/seqs.db "$seq" | last-split -r >> "$outfile"
done