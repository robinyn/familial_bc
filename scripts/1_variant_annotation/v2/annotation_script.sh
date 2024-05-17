#!/bin/bash

# Script to retrieve the flanking sequences of each variant using VCFTools, and perform the annotation in the vcf files using VEP. Installation of all necessary software is detailed in the README file.

# REMEMBER to activate conda env e_vep.

output_directory=$1
resources_directory=$2
data_type=$3
data_directory=$4

# Annotate the files using VEP and store in a new directory.
# Custom annotation to retrieve allele frequency from SweGen database.
# Custom annotation for conservation values from PhyloP and GERP.
printf "Generating base annotations\n"

if [[ ${data_type} == "swea" ]]
then
    file_dir=${output_directory}/2_FlankingSequences/
else
    file_dir=${data_directory}/
fi

total_vcf=$(find ${file_dir} -type f -name "*.vcf" | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${output_directory}/3_BaseAnnotation
rsync -a -f"+ */" -f"- *" ${file_dir} ${output_directory}/3_BaseAnnotation/

find ${file_dir} -type f -name "*.vcf" | while read file
do
    progressBar $current_vcf $total_vcf

    if [[ ${data_type} == "swea" ]]
    then
        output_path=$(echo ${file%/*} | sed "s|2_FlankingSequences|3_BaseAnnotation|")
    else
        output_path=$(echo ${file%/*} | sed "s|${data_directory}|${output_directory}/3_BaseAnnotation|")
    fi

    vep -i $file -o ${output_path}/vep_$(basename "$file") \
    --cache \
    --vcf \
    --fork 4 \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --custom ${resources_directory}/SweGen/anon-SweGen_STR_NSPHS_1000samples_freq_hg19.vcf.gz,SweGen,vcf,exact,0,AF \
    --af \
    --af_1kg \
    --af_gnomade \
    --check_existing \
    --custom ${resources_directory}/All_hg19_RS.bw,GERP,bigwig,exact \
    --custom ${resources_directory}/hg19.100way.phyloP100way.bw,PhyloP,bigwig,exact \
    --fields "Gene","SYMBOL","Feature","Existing_variation","STRAND","EXON","INTRON","Consequence","Codons","AF","EUR_AF","SweGen_AF","gnomADe_AF","gnomADe_NFE_AF","PhyloP","GERP" \
    --no_stats 1>>${output_directory}/3_BaseAnnotation/vep.log 2>>${output_directory}/3_BaseAnnotation/vep.error

    current_vcf=$(($current_vcf+1))

    progressBar $current_vcf $total_vcf
done
