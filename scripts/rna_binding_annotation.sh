#!/bin/bash

root_directory=$1
resources_directory=$2

total_vcf=$(ls ${root_directory}/3_BaseAnnotation/*.vcf;| wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${root_directory}/4_EncodeRBSAnnotation

for f in ${root_directory}/3_BaseAnnotation/*.vcf;
do
    progressBar $current_vcf $total_vcf

    perl -X

    vep -i $f -o ${root_directory}/4_EncodeRBSAnnotation/encode_$(basename "$f") \
    --vcf \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --custom ${resources_directory}/EncodeFiles/encode_rna_binding_1.bed.gz,Encode,bed,overlap \
    --fields "Feature","Encode" \
    --keep_csq \
    --vcf_info_field Encode \
    --no_stats;

    current_vcf=$(($current_vcf+1))
done
