
root_directory=$1
resources_directory=$2

total_vcf=$(ls ${root_directory}/5_GeneNameAnnotation/*.vcf| wc -l)
current_vcf=0

printf "${total_vcf} VCF files to process\n"

mkdir ${root_directory}/6_TargetScanAnnotation

for f in ${root_directory}/5_GeneNameAnnotation/*.vcf;
do
    progressBar $current_vcf $total_vcf

    perl -X

    vep -i $f -o ${root_directory}/6_TargetScanAnnotation/targetscan_$(basename "$f") \
    --vcf \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --custom ${resources_directory}/targetscan/targetscan_miRNA_sorted.bed.gz,TargetScan,bed,overlap \
    --fields "Feature","TargetScan" \
    --keep_csq \
    --vcf_info_field TargetScan \
    --no_stats;

    current_vcf=$(($current_vcf+1))
done