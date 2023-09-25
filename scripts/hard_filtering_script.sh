#!/bin/bash

# Script to perform the Hard Filtering using GATK. The files must be divided since the values are different for indels and SNVs. The filter is applied to both and the files are merged together again.

# Get directories passed as arguments
data_directory=$1
root_directory=$2

# Create the new directories to store the results.
mkdir ${root_directory}/1_Filtering
cd ${root_directory}/1_Filtering
mkdir SNVs
mkdir Indels

# Separate the indels from the snvs using the option SelectVariants from GATK.
printf "Separating indels and SNVs \n"

total_vcf=$(ls ${data_directory}/*.vcf | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

for file in ${data_directory}/*.vcf; do
  
  progressBar $current_vcf $total_vcf

  gatk SelectVariants -V $file -select-type SNP -O SNVs/snvs_${file##*/} 2>>SNVs/snvs_error;
  gatk SelectVariants -V $file -select-type INDEL -select-type MIXED -O Indels/indels_${file##*/} 2>>Indels/indels_error;

  current_vcf=$(($current_vcf+1))
  progressBar $current_vcf $total_vcf
done

# Filter SNVs using the values sugested by the GATK best practices. With the exception of SOR, which we use > 4.0 instead of 3.0. Save the results in a new directory.
printf "Filtering SNVs\n"

total_vcf=$(ls SNVs/*.vcf | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir FilteredSNVs
for file in SNVs/snvs_*.vcf; do
  progressBar $current_vcf $total_vcf

  gatk VariantFiltration -V $file \
  -filter "QD < 2.0 " --filter-name QD2S \
  -filter "MQ < 40.0" --filter-name MQ40 \
  -filter "FS > 60.0" --filter-name FS60 \
  -filter "SOR > 4.0" --filter-name SOR4 \
  -filter "MQRankSum < -12.5" --filter-name MQRS-12.5 \
  -filter "ReadPosRankSum < -8.0" --filter-name RPRS-8 \
  -O FilteredSNVs/filtered_${file##*/} 2>>FilteredSNVs/snvs_filter_error;

  current_vcf=$(($current_vcf+1))
  progressBar $current_vcf $total_vcf
 done


# Filter Indels using the values sugested by the GATK best practices. Save the results in a new directory.
printf "Filtering indels\n"

total_vcf=$(ls Indels/*.vcf | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir FilteredIndels
for file in Indels/indels_*.vcf; do
  progressBar $current_vcf $total_vcf

  gatk VariantFiltration -V $file \
  -filter "QD < 2.0 " --filter-name QD2I  \
  -filter "FS > 200.0" --filter-name FS200 \
  -filter "SOR > 10.0" --filter-name SOR10 \
  -filter "ReadPosRankSum < -20.0" --filter-name RPRS-20 \
  -O FilteredIndels/filtered_${file##*/} 2>>FilteredIndels/indels_filter_error;

  current_vcf=$(($current_vcf+1))
  progressBar $current_vcf $total_vcf
done

# Merge the Indels and SNVs together and store it in a new directory.
printf "Merging SNVs and indels\n"

total_vcf=$(ls FilteredIndels/*.vcf | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir Merged
for file in FilteredIndels/*.vcf; do
  progressBar $current_vcf $total_vcf

  fname=${file#*indels_}; gatk MergeVcfs -I $file -I FilteredSNVs/filtered_snvs_$fname -O Merged/filtered_$fname 1>> exit_code 2>>Merged/merge_error;
  
  current_vcf=$(($current_vcf+1))
  progressBar $current_vcf $total_vcf
done
