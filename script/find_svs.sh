#!/bin/bash
## Synopsis: bash find_svs $INPUT_BAM $OUTPUT_PREFIX
## Requirements: rust and samtools
## It extracts SVs found in the given BAM file.
set -uex 
INPUT_BAM=$1
OUTPUT_PREFIX=$2

samtools view -h "$INPUT_BAM" | cargo run --release --bin filter_large_indel -- --alignments - | samtools sort -O BAM > "$OUTPUT_PREFIX".indel.bam
samtools view -h "$INPUT_BAM" | cargo run --release --bin filter_break_points -- --alignments - | samtools sort -O BAM > "$OUTPUT_PREFIX".stst.bam