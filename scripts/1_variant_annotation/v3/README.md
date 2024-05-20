# Variant Annotation Pipeline V3

This is the most up-to-date version of the pipeline. Several features have been added including a setup script which will download all publically available resources required for the variant annotation, a Docker/Singularity image with all of the required dependencies and scripts, as well as the ability to use multiple cores for parallel processing.

## Installation

### Downloading resource files

Use the [setup script](./setup.sh) to download all of the required resource files. If installing the pipeline manually, uncomment lines 46-54 from the setup script.

Filter files containing filter names and values to be used to filter SNPs and indels are required for the filtering step. Follow the GATK arguments file format. Or download and use the default files [snp_filters.txt](./filter_files/snp_filters.txt) and [indel_filters.txt](./filter_files/indel_filters.txt).

Please note that SweGen data is not publically available, but is required for the current version of the pipeline to function.

```sh
# Run the setup script
bash setup.sh /directory/for/resources/files

# Move the filter files into the resources directory
mv snp_filters.txt /directory/for/resources/files
mv indel_filters.txt /directory/for/resources/files
```

### Using Docker

```sh
# Pull docker image from Docker Hub and build a container
docker pull robinyn/variant_annotation_pipeline_v3
```

### Using Singularity

```sh
# Pull docker image from Docker Hub and build a container
singularity pull docker://robinyn/variant_annotation_pipeline_v3
```

### Manual installation using Conda

Create a new Conda environment using the provided [YAML file](./gatk-vep.yaml).

```sh
# Create Conda environment
conda env create -f gatk-vep.yaml

# Move VEP cache files downloaded with the setup script to VEP directory
mv Homo_sapiens.GRCh37.37.dna.primary_assembly.fa.gz homo_sapiens/111_GRCh37
mv Homo_sapiens.GRCh37.37.dna.primary_assembly.fa.gz.fai homo_sapiens/111_GRCh37
mv Homo_sapiens.GRCh37.37.dna.primary_assembly.fa.gz.gzi homo_sapiens/111_GRCh37
mv homo_sapiens ./vep
```

## Running the pipeline

### Using Docker

```sh
# Start an interactive shell within Docker container and attach the resources, input data, output directories as volumes
docker run -v /resources/directory/on/host:/home/resources -v /input/data/directory/on/host:/home/data -v /output/directory/on/host:/home/output -it robinyn/variant_annotation_pipeline_v3

# Activate the environment
conda activate gatk-vep

# Run the pipeline inside the Docker container using the wrapper script
python wrapper.py annotate -d /input/data/directory -r /resources/directory -o /output/directory -c number_of_cores
```

### Using Singularity

```sh
# Start an interactive shell within Singularity container and attach the resources, input data, output directories as volumes
singularity shell --bind /resources/directory/on/host:/home/resources -bind /input/data/directory/on/host:/home/data -bind /output/directory/on/host:/home/output variant_annotation_pipeline_v3_latest.sif

# Activate the environment
conda activate gatk-vep

# Run the pipeline inside the Singularity container using the wrapper script
python wrapper.py annotate -d /input/data/directory -r /resources/directory -o /output/directory -c number_of_cores
```

### Manual installation using Conda

```sh
# Activate the environment
conda activate gatk-vep

# Run the pipeline using the wrapper script
python wrapper.py annotate -d /input/data/directory -r /resources/directory -o /output/directory -c number_of_cores
```

The wrapper script takes the following command-line arguments as input:
* -d / --data-directory: input data directory
* -o / --output-directory: output data directory
* -r / --resources-directory: resources directory
* -c / --cores: number of cores to use
* --indel-filter: directory and name of the text file containing filter names and values to be used to filter indels. Follow the GATK arguments file format.
* --snp-filter: directory and name of the text file containing filter names and values to be used to filter SNPs. Follow the GATK arguments file format.
* -q / --quiet: supress the ouput of progress bars
* --filt: perform the hard filtering step
* --fs: perform the flanking sequences retrieval
* --vep: perform the VEP annotation step
* --custom: performed the custom annotation step

If no options to perform specific parts of the pipeline are provided, the wrapper will run through the entire pipeline from the filtering step until the custom annotations step.