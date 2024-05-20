# Identification of candidate pathogenic variants in familial breast cancer

## Analysis environment
### Software installation
A conda environment can be created using the provided [YAML file](/gatk-vep.yaml) to install all the required softwares and dependencies for the pipeline.

If the conda environment fails to resolve, the required dependencies must be installed manually. The installation and setup instructions are almost identical to the original version of the [pipeline](https://github.com/aboffelli/variant_annotation).

#### Genome Analysis Toolkit (GATK)

```sh
# Create a bin folder in the root directory if it doesn't exist already
cd ~
mkdir bin

# Download GATK from Broad Institute
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip

# Unzip the downloaded file and remove the zip file
unzip gatk-4.4.0.0.zip
rm gatk-4.4.0.0.zip

# Add the directory to PATH
export PATH=$PATH:~/bin/gatk-4.4.0.0/
```

#### Ensembl Variant Effect Predictor (VEP)

VEP requires the following dependencies before it can be installed:
* gcc
* g++
* Perl (>=5.10)
* Perl libraries [Archive::Zip](https://metacpan.org/pod/Archive::Zip) and [DBI](https://metacpan.org/pod/DBI)

```sh
# Clone the VEP repository from Ensembl
git clone https://github.com/Ensembl/ensembl-vep.git

# Move to the cloned repository
cd ensembl-vep

# Install remaining dependencies. The GRCh37 cache for Homo sapiens can and should also be installed using this script
perl INSTALL.pl
```

The full installation instructions and required dependencies can be found [here](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)

#### PLINK
The PLINK binary can be downloaded from [here](https://www.cog-genomics.org/plink/). Unzip the downloaded file, rename the folder to "plink" and move the folder to the bin folder created earlier.

```sh
# Add the pline folder to PATH
export PATH=$PATH:~/bin/plink
```

#### R packages
While the pipeline does not require R, all of the downstream analyses were performed using R and the following packages are required to reproduce all of the plots and tables from the project.

### Downloading resources

#### VEP Cache files
If the required VEP cache files could not be installed using the "INSTALL.pl" file, they can be downloaded from the Ensembl FTP server.

```sh
# Download and decompress the cache files for VEP (VEP v.110)
wget https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh37.tar.gz
tar xzf homo_sapiens_vep_110_GRCh37.tar.gz

# Download and index the GRCh37 primary assembly fasta file
wget https://ftp.ensembl.org/pub/grch37/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# Move primary assembly file and the index to the downloaded cache directory and move to the VEP cache directory
mv Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz homo_sapiens/110_GRCh37/
mv Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.fai homo_sapiens/110_GRCh37/
mv Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.gzi homo_sapiens/110_GRCh37/

mv homo_sapiens ~/.vep
```

#### PhyloP and GERP

```sh
## PhyloP
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw

## GERP
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw
```

#### ClinVar

```sh
# Download the tab delimited ClinVar annotations
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# Retrieve only GRCh37 lines and save in a file.
zcat variant_summary.txt.gz | grep "GRCh37" > variant_summary_GRCh37.txt
```

#### ENCODE RBP binding sites
```sh
# Download list of files
wget -O encode_urls.txt https://www.encodeproject.org/batch_download/?control_type!=*&status=released&perturbed=false&assay_slims=RNA+binding&assembly=hg19&files.file_type=bed+narrowPeak&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=eCLIP&type=Experiment

# Download the metadata file
wget -O encode_metadata.tsv https://www.encodeproject.org/metadata/?control_type%21=%2A&status=released&perturbed=false&assay_slims=RNA+binding&assembly=hg19&files.file_type=bed+narrowPeak&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=eCLIP&type=Experiment

# Grep files with GRCh37 assembly and Biological replicates 1,2
grep "GRCh37" encode_metadata.tsv | grep "1, 2" > files_to_download.tsv

# Download the selected files
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
```
#### SweGen
The SweGen resources are not publically available and were retrieved from local storage.

## Running the pipeline
The pipeline can now be run using the wrapper script.

```sh
bash scripts/pipeline.sh -d input/data -s resources/directory -f start_flag -t data/type -o output/directory
```
* Start flag determines which step of the pipeline the wrapper will start from:
  * 1 - Filtering
  * 2 - Flanking sequence retrieval
  * 3 - VEP annotations
  * 4 - ENCODE RBS annotations
  * 5 - miRNA target site annotations
  * 6 - ClinVar annotations
* Data type can either be "swea" or "bridges"

For the SWEA data, the pipeline was run from the filtering step (1), and for the BRIDGES data, the pipeline was run from the VEP annotations "3".

## Parsing results
Before any of the downstream analyses could begin, the annotations were parsed into 3 TSV files using the 2 Python parser scripts ([parse_annot.py](./parse_annot.py), [parse_results.py](./parse_results.py)).

```sh
# Run the first parser to combine all the annotations into one large TSV file
python scripts/parse_annot.py -i input/files -o output/directory -t data/type

# Run the second parser to separate them into per_transcript, per_variant, and per_sample files
python scripts/parse_results.py -i input/file -o output/directory -t data/type
```

## Downstream analyses
### AF calculations
The allele frequencies of the variants were also calculated using a Python script ([af_calc.py](./af_calc.py)).
```sh
bash scripts/af_calc.sh list/of/samples/to/exclude path/to/per/sample output/file
```

* The first argument to the script is a list of sampleIDs to exclde. This was added as a way to remove HEBCS samples from the calculations.
* The per_sample_summary.tsv file from the parser is required for the allele frequency calculations.

### Analysis of synonymous variants
The analysis of synonymous variants and the the pseudo-ranking scores were calculated using the [annotations_combined_analysis.R](../../3_results/annotation_combined_analysis.R) script.

The script outputs a ranked table of synonymous variants as final output.

### Association analysis
The input files for the association analyses performed with the BRIDGES data were created using the [association_analysis.R](../../2_association_analysis/association_analysis.R) script.

The analysis was then run using PLINK.

```sh
# Run PLINK
plink --file input/file --no-fid --no-parents --assoc --compound-genotypes --adjust

# Reformat output files
cat plink.assoc | tr -s ' ' | sed 's/^[[:space:]]//' | sed 's/[[:space:]]$//' | tr ' ' '\t' > output_reformatted.txt

cat plink.assoc.adjusted | tr -s ' ' | sed 's/^[[:space:]]//' | sed 's/[[:space:]]$//' | tr ' ' '\t' > output_adjusted_reformatted.txt
```

* PLINK only needs the basename of the two input files (e.g. if the input files were named assoc.ped, assoc.map, then the input file name is assoc).

The results from the association analyses were then read into R using [assoc_analysis_results.R](../../3_results/assoc_analysis_results.R) and combined with the annotations from the pipeline.

### Summary statistics and plotting
All of the summary statistics and plots found in the report can be generating using the following scripts:
* [summary_stats.R](../../3_results/summary_stats.R)
* [result_plots.R](../../3_results/result_plots.R)