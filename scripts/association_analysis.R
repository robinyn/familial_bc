## Title: association_analysis.R
## Author: Euisuk Robin Han
## Description: A script for the case-control association analysis of BRIDGES data 
## Dependencies: Tested and stable with: 
##                R version 4.2.3
##                tidyverse (version 2.0.0)

library(tidyverse)

# Set working directory
setwd("~/Desktop/thesis/association_analysis/")

variants = read_tsv("~/Desktop/thesis/bridges_results/bridges_synonymous_list.tsv") %>% 
  mutate(variant=str_remove(variant, "chr"))

pathogenic_variants = read_tsv("~/Desktop/thesis/bridges_results/bridges_pathogenic_list.tsv") %>% 
  mutate(variant=str_remove(variant, "chr"))

# Filter to remove non-synonymous variants, check to make sure they are SNPs 
variant_list = variants %>% 
  mutate(variant=str_remove(variant, "chr")) %>% 
  mutate(ref_allele=str_split(variant, "-", simplify=TRUE)[,3]) %>% 
  mutate(alt_allele=str_split(variant, "-", simplify=TRUE)[,4]) %>% 
  filter(str_length(ref_allele)==1 & str_length(alt_allele)==1) %>% 
  dplyr::select(c("variant", "samples", "ref_allele", "alt_allele")) %>% 
  unique()

# Read phenotype annotations for the samples
control_phenotypes = read_tsv("bridges_annotation/controls_phenotypes.txt") %>% 
  filter(study!="HEBCS") # HEBCS withdrew from BCAC
case_phenotypes = read_tsv("bridges_annotation/cases_phenotypes.txt") %>% 
  filter(study!="HEBCS") # HEBCS withdrew from BCAC

# Filter control samples to exclude HEBCS samples and filter by population
control_samples = control_phenotypes %>% 
  #filter(famHist==1) %>% # Select samples with family history
  #filter(ethnicityClass==1) %>% # Select population (1 - European, 5 - East Asian)
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only includes female samples
  mutate(status=1)

# Filter case samples to exclude HEBCS samples and filter by population
case_samples = case_phenotypes %>% 
  filter(fhnumber>1 & fhnumber!=888 & fhnumber!=777) %>% # Select samples with family history
  #filter(AgeDiagIndex<50) %>% 
  #filter(ethnicityClass==1) %>%  # Select population (1 - European, 5 - East Asian)
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only includes female samples
  mutate(status=2) 

# Create a sample list to create the PED files with
sample_list = control_samples %>% rbind(case_samples)

# Read the per_sample summary file to generate the PED files from
per_sample = read_tsv("~/Desktop/thesis/bridges_results/per_sample_summary.tsv") %>% 
  mutate(variant=str_remove(variant, "chr"))

samples_with_pathogenic=per_sample %>% 
  filter(variant %in% pathogenic_variants$variant) 

# Filter variants and samples to make sure they passed the past filtering steps
per_sample = per_sample %>% 
  filter(!sampleID %in% (samples_with_pathogenic$sampleID %>% unique())) %>% 
  filter(sampleID %in% sample_list$BRIDGES_ID) %>% 
  filter(variant %in% variant_list$variant)

# Pivot the per_sample table to a wider format to match the PED format
ped_out = per_sample %>% pivot_wider(names_from = variant, values_from = genotype)

# Join the sample_list with ped_out by the sampleIDs to create the final PED file
ped_out = left_join(sample_list, ped_out, by=join_by("BRIDGES_ID"=="sampleID"))

# Create a map file from the PED file 
map_out = colnames(ped_out)[-c(1,2,3)] %>% 
  as.data.frame() %>% 
  setNames("variantID") %>% 
  mutate(chromosome=str_split(variantID, "-", simplify = TRUE)[,1]) %>% 
  mutate(bp_position=str_split(variantID, "-", simplify = TRUE)[,2]) %>% 
  dplyr::select(c(chromosome, variantID, bp_position))

# Replace NA with homogeneous reference genotype
ped_out = ped_out %>% 
  mutate(across(colnames(ped_out)[-c(1,2,3)], 
                ~replace_na(., sprintf("%1$s%1$s", str_split(cur_column(), "-", simplify = TRUE)[[3]]))))

# Write the generated PED/MAP files to disk
write_tsv(ped_out, "assoc_analysis_no_pathogenic_fam_hist_1st_dg.ped", na="00", col_names = FALSE)
write_tsv(map_out, "assoc_analysis_no_pathogenic_fam_hist_1st_dg.map", na="00", col_names = FALSE)




