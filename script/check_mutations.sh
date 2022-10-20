#!/bin/bash
# Synopsis: check_mutations.sh <REFERENCE>(fasta) <STARAINS>(fasta) <OUTPUT_PREFIX>
# Require: rust, minimap2, mafft.
REFERENCE=$1
STRAINS=$2
PREFIX=$3

mkdir -p "$PREFIX"
## 1. Aligns all the sequence into the reference.
minimap2 -x asm5 --eqx --secondary=no -c "$REFERENCE" "$STRAINS" > "$PREFIX"/alns.paf
cargo run --release --bin diff_aln -- --alignments "$PREFIX"/alns.paf --prefix "$PREFIX" > "$PREFIX"/alns.diff

## 2. Aggregate these sequences.
mkdir -p "$PREFIX"/subseqs/
cargo run --release --bin focus_on_mismatches -- --reference "$REFERENCE" --contigs "$STRAINS" --alignments "$PREFIX"/alns.paf --prefix "$PREFIX"/subseqs/

## 3. MAFFT 
mkdir -p "$PREFIX"/mafft/
for subseq in "$PREFIX"/subseqs/*.fa
do
    filename=${subseq#"$PREFIX"/subseqs/}
    filename=${filename%.fa}
    mafft --preservecase --auto "$subseq" 2> /dev/null > "$PREFIX"/mafft/"$filename".mul.fa
done

## 4. Filtering only supports from two or more reads.
mkdir -p "$PREFIX"/filtered/
touch "$PREFIX"/variant_nums.tsv
rm "$PREFIX"/variant_nums.tsv
cargo build --release
for subseq in "$PREFIX"/mafft/*.mul.fa
do
    if "$PWD"/target/release/is_supported_by --pileup "$subseq" >> "$PREFIX"/variant_nums.tsv
    then
        cp "$subseq" "$PREFIX"/filtered/
    fi
done