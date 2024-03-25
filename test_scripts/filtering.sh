#!/bin/bash

# Script to perform the Hard Filtering using GATK. The files must be divided since the values are different for indels and SNVs. The filter is applied to both and the files are merged together again.

# Get directories passed as arguments
data_directory=$1
output_directory=$2
snp_filters=$3
indel_filters=$4

# Create the new directories to store the results.
mkdir ${output_directory}/Filtering
mkdir ${output_directory}/Filtering/SNVs
mkdir ${output_directory}/Filtering/Indels

# Separate the indels from the snvs using the option SelectVariants from GATK.
printf "Separating indels and SNVs \n"

total_vcf=$(ls ${data_directory}/*.vcf | wc -l)

printf "${total_vcf} VCF files to process\n"

for file in ${data_directory}/*.vcf; do

  gatk SelectVariants -V $file -select-type SNP -O ${output_directory}/Filtering/SNVs/snvs_${file##*/} 2>>${output_directory}/Filtering/SNVs/snvs_error;
  gatk SelectVariants -V $file -select-type INDEL -select-type MIXED -O ${output_directory}/Filtering/Indels/indels_${file##*/} 2>>${output_directory}/Filtering/Indels/indels_error;

done

# Filter SNVs using the values sugested by the GATK best practices. With the exception of SOR, which we use > 4.0 instead of 3.0. Save the results in a new directory.
printf "Filtering SNVs\n"

total_vcf=$(ls ${output_directory}/Filtering/SNVs/*.vcf | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${output_directory}/Filtering/FilteredSNVs
for file in ${output_directory}/Filtering/SNVs/snvs_*.vcf; do

  gatk VariantFiltration -V $file \
  --arguments_file ${snp_filters} \
  -O ${output_directory}/Filtering/FilteredSNVs/filtered_${file##*/} 2>>${output_directory}/Filtering/FilteredSNVs/snvs_filter_error;

 done

# Filter Indels using the values sugested by the GATK best practices. Save the results in a new directory.
printf "Filtering indels\n"

total_vcf=$(ls ${output_directory}/Filtering/Indels/*.vcf | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${output_directory}/Filtering/FilteredIndels
for file in ${output_directory}/Filtering/Indels/indels_*.vcf; do

  gatk VariantFiltration -V $file \
  --arguments_file ${indel_filters} \
  -O ${output_directory}/Filtering/FilteredIndels/filtered_${file##*/} 2>>${output_directory}/Filtering/FilteredIndels/indels_filter_error;

done

# Merge the Indels and SNVs together and store it in a new directory.
printf "Merging SNVs and indels\n"

total_vcf=$(ls ${output_directory}/Filtering/FilteredIndels/*.vcf | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${output_directory}/Filtering/Merged
for file in ${output_directory}/Filtering/FilteredIndels/*.vcf; do

  fname=${file#*indels_}; gatk MergeVcfs -I $file -I ${output_directory}/Filtering/FilteredSNVs/filtered_snvs_$fname \
  -O ${output_directory}/Filtering/Merged/filtered_$fname 1>> ${output_directory}/Filtering/exit_code 2>>${output_directory}/Filtering/Merged/merge_error;

done
