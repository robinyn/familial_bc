#!/bin/bash

# Script to perform the Hard Filtering using GATK. The files must be divided since the values are different for indels and SNVs. The filter is applied to both and the files are merged together again.

# Get directories passed as arguments
file=$2
output_directory=$3
snp_filters=$4
indel_filters=$5

function separateVariants {
  # Separate the indels from the snvs using the option SelectVariants from GATK.
  gatk SelectVariants -V $file -select-type SNP \
  -O ${output_directory}/SNVs/snvs_${file##*/} 2>>${output_directory}/SNVs/snvs_error;
  gatk SelectVariants -V $file -select-type INDEL -select-type MIXED \
  -O ${output_directory}/Indels/indels_${file##*/} 2>>${output_directory}/Indels/indels_error;
}

function filterSNVs {
  # Filter SNVs using the values sugested by the GATK best practices. With the exception of SOR, which we use > 4.0 instead of 3.0. Save the results in a new directory.
  gatk VariantFiltration -V $file \
  --arguments_file ${snp_filters} \
  -O ${output_directory}/FilteredSNVs/filtered_${file##*/} 2>>${output_directory}/FilteredSNVs/snvs_filter_error;
}

function filterIndels {
  # Filter Indels using the values sugested by the GATK best practices. Save the results in a new directory.
  gatk VariantFiltration -V $file \
  --arguments_file ${indel_filters} \
  -O ${output_directory}/FilteredIndels/filtered_${file##*/} 2>>${output_directory}/FilteredIndels/indels_filter_error;
}

function mergeVariants {
  # Merge the Indels and SNVs together and store it in a new directory.
  fname=${file#*indels_}; gatk MergeVcfs -I $file -I ${output_directory}/FilteredSNVs/filtered_snvs_$fname \
  -O ${output_directory}/Merged/filtered_$fname 1>> ${output_directory}/exit_code 2>>${output_directory}/Merged/merge_error;
}

"$@"