root_directory=$1

total_vcf=$(ls ${root_directory}/4_EncodeRBSAnnotation/*.vcf| wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${root_directory}/5_GeneNameAnnotation

for f in ${root_directory}/4_EncodeRBSAnnotation/*.vcf; do
    progressBar $current_vcf $total_vcf

    vep -i $f -o ${root_directory}/5_GeneNameAnnotation/gene_name_$(basename "$f") \
    --cache \
    --vcf \
    --offline \
    --assembly GRCh37 \
    --no_stats \
    --fields "Gene","SYMBOL" \
    --force \
    --keep_csq \
    --vcf_info_field Gene 1>>${root_directory}/5_GeneNameAnnotation/vep.log 2>>${root_directory}/5_GeneNameAnnotation/vep.error

    current_vcf=$(($current_vcf+1))
    progressBar $current_vcf $total_vcf
done
