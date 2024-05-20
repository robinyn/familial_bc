base_dir=$1

cd ${base_dir}

## PhyloP
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw

## GERP
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw

# Download latest ClinVar VCF annotations for GRCh37 (version used for annotation)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi

# Download and process ENCODE RBS data
mkdir ${base_dir}/EncodeFiles
cd ${base_dir}/EncodeFiles

wget -O encode_urls.txt https://www.encodeproject.org/batch_download/?control_type!=*&status=released&perturbed=false&assay_slims=RNA+binding&assembly=hg19&files.file_type=bed+narrowPeak&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=eCLIP&type=Experiment
wget -O encode_metadata.tsv https://www.encodeproject.org/metadata/?control_type%21=%2A&status=released&perturbed=false&assay_slims=RNA+binding&assembly=hg19&files.file_type=bed+narrowPeak&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=eCLIP&type=Experiment

grep "GRCh37" encode_metadata.tsv | grep "1, 2" > files_to_download.tsv
cat ~/files_to_download.tsv | cut -f1 | while read line; do file=$(grep "$line" encode_urls.txt); wget $file ; done

## Decompress the bed files.
for file in *; do filename=${file%.gz}; gunzip -c $file > $filename; done

## Join only column 1-8 from all files, and remove chr from the chromosome column and _IDR from prot column.
for f in *.bed; do
     while read line; do
         feature=$(echo "$line" | cut -f4-8 | tr $'\t' ':');
         pos=$(echo "$line" | cut -f1-3);
         echo "$pos"$'\t'"$feature" | sed 's/_IDR//';
     done < "$f" >> encode_rna_binding.bed;
 done

## Sort the file and compress
cat encode_rna_binding.bed | sort -k1,1 -k2,2n -k3,3n | bgzip -c > encode_rna_binding_1.bed.gz

## Index with tabix
tabix -p bed encode_rna_binding_1.bed.gz

cd ${base_dir}/

# Download VEP cache files
# mkdir ${base_dir}/vep_cache
# cd ${base_dir}/vep_cache

# Download cache files for VEP (VEP v.110)
# wget https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh37.tar.gz
# tar xzf homo_sapiens_vep_110_GRCh38.tar.gz

# wget https://ftp.ensembl.org/pub/grch37/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
# samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
