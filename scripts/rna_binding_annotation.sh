#!/bin/bash

output_directory=$1
resources_directory=$2

file_dir=${output_directory}/3_BaseAnnotation/

total_vcf=$(find ${file_dir} -type f -name "*.vcf" | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${output_directory}/4_EncodeRBSAnnotation
rsync -a -f"+ */" -f"- *" ${file_dir} ${output_directory}/4_EncodeRBSAnnotation

find ${file_dir} -type f -name "*.vcf" | while read file
do
    progressBar $current_vcf $total_vcf

    path_name=$(echo ${file%/*} | sed "s|3_BaseAnnotation|4_EncodeRBSAnnotation|")

    vep -i $file -o ${path_name}/encode_$(basename "$file") \
    --vcf \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --custom ${resources_directory}/EncodeFiles/encode_rna_binding_1.bed.gz,Encode,bed,overlap \
    --fields "Feature","Encode" \
    --keep_csq \
    --vcf_info_field Encode \
    --no_stats 1>>${output_directory}/4_EncodeRBSAnnotation/vep.log 2>>${output_directory}/4_EncodeRBSAnnotation/vep.error

    current_vcf=$(($current_vcf+1))
    progressBar $current_vcf $total_vcf
done
