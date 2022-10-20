#!/bin/bash
# Synopsis bash repetitive_analysis.sh $REFERENCE $ASSEMBLIES $PREFIX
# Requires: last, rust.
set -uex

REFERENCE=$1
ASSEMBLIES=$2
PREFIX=$3
mkdir -p "$PREFIX"

cargo run --release --bin count_kmers -- --genomes "$ASSEMBLIES" > "$PREFIX"/kmer_counts.tsv
cargo run --release --bin annotate_repetitive_kmers -- --genomes "$ASSEMBLIES" > "$PREFIX"/repeat_annotations.tsv

lastdb "$PREFIX"/ref_genome_db "$REFERENCE"
last-train "$PREFIX"/ref_genome_db "$ASSEMBLIES" > "$PREFIX"/param.par
lastal -p "$PREFIX"/param.par "$PREFIX"/ref_genome_db "$ASSEMBLIES" | last-split | maf-convert tab > "$PREFIX"/aln.tsv
cargo run --release --bin annotate_break_points -- --alignments "$PREFIX"/aln.tsv > "$PREFIX"/break_points.tsv

cargo run --release --bin compare_annotations -- \
    --genomes "$ASSEMBLIES" --repeat-annotation "$PREFIX"/repeat_annotations.tsv \
    --break-point-annotation "$PREFIX"/break_points.tsv --PREFIX-prefix "$PREFIX" > "$PREFIX"/p_value.txt
