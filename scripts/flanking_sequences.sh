#!/bin/bash

output_directory=$1
resources_directory=$2

# Get the flanking sequences of the merged files and store in a new directory.
printf "Retrieving flanking sequences\n"

total_vcf=$(ls ${output_directory}/1_Filtering/Merged/*.vcf| wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${output_directory}/2_FlankingSequences
for file in ${output_directory}/1_Filtering/Merged/*.vcf;
do
    progressBar $current_vcf $total_vcf

    fill-fs -r ~/.vep/homo_sapiens/110_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz -l 15 $file > ${output_directory}/2_FlankingSequences/fs_$(basename "$file");

    current_vcf=$(($current_vcf+1))
    progressBar $current_vcf $total_vcf
done
