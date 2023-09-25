#!/bin/bash

# Script to retrieve the flanking sequences of each variant using VCFTools, and perform the annotation in the vcf files using VEP. Installation of all necessary software is detailed in the README file.

# REMEMBER to activate conda env e_vep.

root_directory=$1
resources_directory=$2

# Annotate the files using VEP and store in a new directory.
# Custom annotation to retrieve allele frequency from SweGen database.
# Custom annotation for conservation values from PhyloP and GERP.
printf "Generating base annotations\n"

total_vcf=$(ls ${root_directory}/2_FlankingSequences/*.vcf| wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${root_directory}/3_BaseAnnotation
for file in ${root_directory}/2_FlankingSequences/*.vcf
do
    progressBar $current_vcf $total_vcf

    perl -X

    vep -i $file -o ${root_directory}/3_BaseAnnotation/$(basename "$file" .vcf)_VEP.vcf \
    --cache \
    --vcf \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --custom ${resources_directory}/SweGen/anon-SweGen_STR_NSPHS_1000samples_freq_hg19.vcf.gz,SweGen,vcf,exact,0,AF \
    --af \
    --af_1kg \
    --af_gnomad \
    --check_existing \
    --custom ${resources_directory}/All_hg19_RS.bw,GERP,bigwig,exact \
    --custom ${resources_directory}/hg19.100way.phyloP100way.bw,PhyloP,bigwig,exact \
    --fields "Feature","Existing_variation","STRAND","EXON","INTRON","Consequence","Codons","AF","EUR_AF","SweGen_AF","gnomAD_AF","gnomAD_NFE_AF","PhyloP","GERP", \
    --no_stats

    current_vcf=$(($current_vcf+1))
 done
