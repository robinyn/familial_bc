
output_directory=$1
resources_directory=$2

file_dir=${output_directory}/4_EncodeRBSAnnotation/

total_vcf=$(find ${file_dir} -type f -name "*.vcf" | wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${output_directory}/5_TargetScanAnnotation
rsync -a -f"+ */" -f"- *" ${file_dir} ${output_directory}/5_TargetScanAnnotation

find ${file_dir} -type f -name "*.vcf" | while read file
do
    progressBar $current_vcf $total_vcf

    path_name=$(echo ${file%/*} | sed "s|4_EncodeRBSAnnotation|5_TargetScanAnnotation|")

    vep -i $file -o ${path_name}/targetscan_$(basename "$file") \
    --vcf \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --custom ${resources_directory}/targetscan/targetscan_miRNA_sorted.bed.gz,TargetScan,bed,overlap \
    --fields "Feature","TargetScan" \
    --keep_csq \
    --vcf_info_field TargetScan \
    --no_stats 1>>${output_directory}/5_TargetScanAnnotation/vep.log 2>>${output_directory}/5_TargetScanAnnotation/vep.error

    current_vcf=$(($current_vcf+1))
    progressBar $current_vcf $total_vcf
done
