while getopts i:o:r: flag
do
    case "${flag}" in
        i) file=${OPTARG};;
		o) output_directory=${OPTARG};;
        r) resources_directory=${OPTARG};;
    esac
done

vep -i $file -o ${output_directory}/vep_$(basename "$file") \
    --cache \
    --vcf \
    --fork 4 \
    --offline \
    --force \
    --assembly GRCh37 \
    --distance 0 \
    --check_existing \
    --af \
    --af_1kg \
    --af_gnomade \
    --custom file=${resources_directory}/SweGen/anon-SweGen_STR_NSPHS_1000samples_freq_hg19.vcf.gz,short_name=SweGen,format=vcf,type=exact,fields=AF \
    --custom file=${resources_directory}/All_hg19_RS.bw,short_name=GERP,format=bigwig,type=exact \
    --custom file=${resources_directory}/hg19.100way.phyloP100way.bw,short_name=PhyloP,format=bigwig,type=exact \
    --custom file=${resources_directory}/EncodeFiles/encode_rna_binding_1.bed.gz,short_name=Encode,format=bed,type=overlap \
    --custom file=${resources_directory}/targetscan/targetscan_miRNA_sorted.bed.gz,short_name=TargetScan,format=bed,type=overlap \
    --custom file=${resources_directory}/clinvar_20240312.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields="ALLELEID%CLNDISDB%\
    CLNDN%CLNHGVS%CLNREVSTAT%CLNSIG%CLNSIGCONF%CLNSIGINCL%CLNVC%CLNVCSO%CLNVI%DBVARID%GENEINFO%MC%ORIGIN" \
    --fields "Allele,Gene,SYMBOL,Feature,BIOTYPE,Existing_variation,STRAND,EXON,INTRON,Consequence,Codons,AF,EUR_AF,SweGen_AF,gnomADe_AF,gnomADe_NFE_AF,PhyloP,GERP,\
    Encode,TargetScan,ClinVar,ClinVar_ALLELEID,ClinVar_CLNDISDB,ClinVar_CLNDN,ClinVar_CLNHGVS,ClinVar_CLNREVSTAT,ClinVar_CLNSIG,ClinVar_CLNSIGCONF,ClinVar_CLNSIGINCL,\
    ClinVar_CLNVC,ClinVar_CLNVCSO,ClinVar_CLNVI,ClinVar_DBVARID,ClinVar_GENEINFO,ClinVar_MC,ClinVar_ORIGIN" \
    --no_stats 1>>${output_directory}/vep.log 2>>${output_directory}/vep.error
