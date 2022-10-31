#!/bin/bash
echo "---"
echo "  genomes:"

for strain in an1 c24 col0 cvi eri kyo ler pacbio sha reference
do
    fasta="$PWD"/result/annotations/"$strain".fa
    annotation="$PWD"/result/annotations/"$strain".tsv
    echo "   -"
    echo "    name: ${strain}"
    echo "    sequence_files:"
    echo "       - ${fasta}"
    echo "    feature_files:"
    echo "      gene:"
    echo "       - ${annotation}"
done


echo "  alignment:"
echo "      program: importer"
echo "      parameter:"
for strain in an1 c24 col0 cvi eri kyo ler pacbio sha reference
do
    maf="$PWD"/result/annotations/"$strain".maf
    echo "       - ${maf}"
done
