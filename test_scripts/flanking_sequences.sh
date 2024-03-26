#!/bin/bash

input_directory=$1
output_directory=$2
resources_directory=$3

# Get the flanking sequences of the merged files and store in a new directory.
printf "Retrieving flanking sequences\n"

printf "${total_vcf} VCF files to process\n"

for file in $input_directory;
do
    fill-fs -r ~/.vep/homo_sapiens/110_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz -l 15 $file > ${output_directory}/fs_$(basename "$file");
done
