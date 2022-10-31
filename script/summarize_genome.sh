#!/bin/bash
## Synopsis: <REFERENCE:fasta> <Annotation:GFF3> <CONTIGS:[fasta]> <OUTPUT_DIR> <OUTPUT_PREFIX>
## Require: rust, liftoff, last
set -uex

cargo --version >&2
lastal --version >&2 
liftoff --version >&2

REFERENCE=$1
ANNOTATION=$2
CONTIGS=$3
OUTPUT_DIR=$4
PREFIX=$5

## 0. Prepare outputs
mkdir -p "$OUTPUT_DIR"
OUT_GFF="$OUTPUT_DIR"/"$PREFIX".gff3
LAST_DB="$OUTPUT_DIR"/"$PREFIX"_lastdb
LAST_PAR="$OUTPUT_DIR"/"$PREFIX".par
LAST_MAF="$OUTPUT_DIR"/"$PREFIX".maf
VAR_INFO="$OUTPUT_DIR"/"$PREFIX".tsv

## 1. Annotation
TEMP="$RANDOM"
liftoff -o "$OUT_GFF" -dir "$TEMP" -g "$ANNOTATION" "$CONTIGS" "$REFERENCE" >&2
rm -r "$TEMP"

## 2. Diff between references.
lastdb "$LAST_DB" "$REFERENCE"
last-train "$LAST_DB" "$CONTIGS" > "$LAST_PAR"
lastal -p "$LAST_PAR" "$LAST_DB" "$CONTIGS" | last-split | last-split -r > "$LAST_MAF"
rm "$LAST_DB"*

## Summarize.
### Output to stdout.
cargo run --release --bin summarize_annotation -- --contigs "$CONTIGS" --gff "$ANNOTATION" --maf "$LAST_MAF" --output "$VAR_INFO"