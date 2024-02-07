## Title: association_analysis.R
## Author: Euisuk Robin Han
## Last modified: 02 Feb 24
## Description: A script for the case-control association analysis of BRIDGES data 
## Dependencies: Tested and stable with: 
##                R version 4.2.3
##                tidyverse (version 2.0.0)

library(tidyverse)

# Set working directory
setwd("~/Desktop/thesis/association_analysis/")

# Read per_transcript/per_variant summary files into memory
per_transcript=read_tsv("~/Desktop/thesis/bridges_results/per_transcript_summary.tsv") %>% 
  mutate(variant=str_remove(variant, "chr"))
per_variant=read_tsv("~/Desktop/thesis/bridges_results/per_variant_summary.tsv") %>%
  dplyr::select(-c(gnome_AD_AF, gnome_AD_NFE_AF)) %>% 
  mutate(variant=str_remove(variant, "chr"))

# Read targeted gene list for BRIDGES cohort
gene_list=read_tsv("~/Desktop/thesis/bridges_gene_list.txt", col_names = F)

# Remove any variants mapping to non-targeted genes
per_transcript = per_transcript %>%
  filter(gene %in% gene_list$X1)

# Remove any variants mapping to non-targeted genes
per_variant = per_variant %>%
  separate_longer_delim(cols=gene, delim="|") %>%
  filter(gene %in% gene_list$X1) %>%
  mutate(p=10**-abs(phyloP)) %>% # Create a new column and convert phyloP scores to p-values
  mutate(FDR=p.adjust(p, method="BH", n=2861327195)) # Adjust for multiple testing (n=number of non N bases in hg19)

# Select relevant columns from per_transcript summary and join with per_variant summary
variants = per_transcript %>%
  dplyr::select(c(variant, transcript_id, variant_type, codon, encode, rscu, 
                  ref_ese, alt_ese, ref_ess, alt_ess, miRNA_target)) %>% 
  left_join(per_variant %>% dplyr::select(-type), by = join_by(variant==variant)) 

# Create new columns and extract information from ClinVar annotation
variants = variants %>% 
  mutate(type=variant_type) %>%
  dplyr::select(-variant_type) %>% 
  mutate(pathogenicity=if_else(!is.na(ClinVar), str_split(ClinVar, "\\|", simplify = TRUE)[,7], NA)) %>% 
  mutate(clinvar_phenotype=if_else(!is.na(ClinVar),str_split(ClinVar, "\\|", simplify=TRUE)[,14], NA)) %>% 
  mutate(phenotype_BOC=if_else(str_detect(clinvar_phenotype, 
                                          "Breast|breast|Ovarian|ovarian|Hereditary_cancer-predisposing_syndrome"), TRUE, FALSE)) %>% 
  mutate(clinvar_id=str_split(ClinVar, "\\|", simplify=TRUE)[,3]) %>% 
  group_by(across(-c(transcript_id, codon, type))) %>% 
  summarise(transcript_id=paste(transcript_id, collapse="|"),
            codon=paste(codon, collapse="|"),
            type=paste(type, collapse="|")) %>% 
  relocate(any_of(c("gene", "known_variation", "type", "AF", "EUR_AF", 
                    "SweAF", "rscu", "phyloP", "phyloP_sig", "ClinVar")), 
           .after=variant) %>% 
  ungroup()

# Filter to remove non-synonymous variants, check to make sure they are SNPs 
variant_list = variants %>% 
  filter(str_detect(type, "synonymous")) %>% 
  mutate(ref_allele=str_split(variant, "-", simplify=TRUE)[,3]) %>% 
  mutate(alt_allele=str_split(variant, "-", simplify=TRUE)[,4]) %>% 
  filter(str_length(ref_allele)==1 & str_length(alt_allele)==1) %>% 
  dplyr::select(c("variant", "samples", "ref_allele", "alt_allele")) %>% 
  unique()

#synonymous_variants = variants %>% filter(variant %in% variant_list$variant) %>% unique()

# Read phenotype annotations for the samples
control_phenotypes = read_tsv("bridges_annotation/controls_phenotypes.txt")
case_phenotypes = read_tsv("bridges_annotation/cases_phenotypes.txt")

# Filter control samples to exclude HEBCS samples and filter by population
control_samples = control_phenotypes %>% 
  filter(study!="HEBCS") %>% # HEBCS withdrew from BCAC
  #filter(famHist==1) %>% # Select samples with family history
  filter(ethnicityClass==1) %>% # Select population (1 - European, 5 - East Asian)
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only includes female samples
  mutate(status=1)

# Filter case samples to exclude HEBCS samples and filter by population
case_samples = case_phenotypes %>% 
  filter(study!="HEBCS") %>% # HEBCS withdrew from BCAC
  filter(famHist==1) %>% # Select samples with family history
  filter(ethnicityClass==1) %>%  # Select population (1 - European, 5 - East Asian)
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only includes female samples
  mutate(status=2) 

# Create a sample list to create the PED files with
sample_list = control_samples %>% rbind(case_samples)

# Read the per_sample summary file to generate the PED files from
per_sample = read_tsv("~/Desktop/thesis/bridges_results/per_sample_summary.tsv")

# Filter variants and samples to make sure they passed the past filtering steps
per_sample = per_sample %>% 
  mutate(variant=str_remove(variant, "chr")) %>% 
  filter(variant %in% variant_list$variant) %>% 
  filter(sampleID %in% sample_list$BRIDGES_ID)

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

######## REPLACE WITH DATA FROM combined_analysis.R FOR CONSISTENCY ##########
# Identify pathogenic variants to remove samples with them from the analysis
samples_pathogenic = per_sample %>% 
  left_join(variants %>% dplyr::select(variant, pathogenicity) %>% unique(), by=join_by(variant==variant))

samples_pathogenic = samples_pathogenic %>% 
  filter(pathogenicity=="Pathogenic" | pathogenicity=="Likely_pathogenic" | pathogenicity=="Pathogenic/Likely_pathogenic") %>% 
  dplyr::select(sampleID) %>% 
  unique()

# Remove samples with pathogenic variants
ped_out = ped_out %>% 
  filter(!(BRIDGES_ID %in% samples_pathogenic$sampleID))

# Replace NA with homogeneous reference genotype
ped_out = ped_out %>% 
  mutate(across(colnames(ped_out)[-c(1,2,3)], 
                ~replace_na(., sprintf("%1$s%1$s", str_split(cur_column(), "-", simplify = TRUE)[[3]]))))

# Write the generated PED/MAP files to disk
write_tsv(ped_out, "assoc_analysis_no_pathogenic_fam_hist_european.ped", na="00", col_names = FALSE)
write_tsv(map_out, "assoc_analysis_no_pathogenic_fam_hist_european.map", na="00", col_names = FALSE)




