#!/bin/bash

file=$1
output_directory=$2

temp_file=temp_${file##*/}

bcftools annotate -c FISHER:=INFO/FS $file > ${output_directory}/temp/${temp_file}

fill-fs -r ~/.vep/homo_sapiens/111_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz -l 15 \
    ${output_directory}/temp/${temp_file} > ${output_directory}/temp/temp_fs_$(basename "$file");

bcftools annotate -c FSEQ:=INFO/FS ${output_directory}/temp/temp_fs_$(basename "$file") > ${output_directory}/temp/fs_$(basename "$file")
bcftools annotate -c FS:=INFO/FISHER ${output_directory}/temp/fs_$(basename "$file") > ${output_directory}/fs_$(basename "$file")

rm -r ${output_directory}/temp/*
