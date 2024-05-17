# Identification of candidate pathogenic variants in familial breast cancer

## Analysis environment
### Software installation
A conda environment can be created using the provided [YAML file](/gatk-vep.yaml) to install all the required softwares and dependencies for the pipeline.

If the conda environment fails to resolve, the required dependencies must be installed manually.

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

## Running the pipeline

## Downstream analysis