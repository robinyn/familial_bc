## Title: combined_analysis.R
## Author: Euisuk Robin Han
## Description: A script for the analysis of the annotated SWEA/BRIDGES data
## Dependencies: Tested and stable with: 
##                R version 4.2.3
##                tidyverse (version 2.0.0)
##                ggplot2 (version 3.4.4)

library(tidyverse)
library(ggplot2)
library(org.Hs.eg.db)
library(DBI)

setwd("~/Desktop/thesis")

################################## SELECT CORRECT DATA TYPE ##################################
data_type = "bridges"
##############################################################################################

if(data_type=="swea"){
  # Read variant summary files - remove gnome_AD columns since the annotation no longer works
  per_transcript=read_tsv("swea_results/per_transcript_summary.tsv")
  per_variant=read_tsv("swea_results/per_variant_summary.tsv") %>%
    dplyr::select(-c(gnome_AD_AF, gnome_AD_NFE_AF))
  
  per_sample=read_tsv("swea_results/per_sample_summary.tsv") 
  
  # Read the list of genes targeted for the analysis
  gene_list=read_tsv("swea_gene_list.txt", col_names = F)
}

if(data_type=="bridges"){
  # Read sample phenotype annotations 
  case_phenotypes = read_tsv("association_analysis/bridges_annotation/cases_phenotypes.txt") %>% 
    filter(study!="HEBCS") # HEBCS withdrew from BCAC
  control_phenotypes = read_tsv("association_analysis/bridges_annotation/controls_phenotypes.txt") %>% 
    filter(study!="HEBCS") # HEBCS withdrew from BCAC
  
  # Filter control samples to exclude HEBCS samples
  control_samples = control_phenotypes %>%
    dplyr::select(c(BRIDGES_ID))

  # Filter case samples to exclude HEBCS samples
  case_samples = case_phenotypes %>%
    dplyr::select(c(BRIDGES_ID))
  
  # Create a sample list
  sample_list = control_samples %>% rbind(case_samples)
  
  # Read variant summary files - remove gnome_AD columns since the annotation no longer works
  per_transcript=read_tsv("bridges_results/per_transcript_summary.tsv")
  per_variant=read_tsv("bridges_results/per_variant_summary.tsv") %>%
    dplyr::select(-c(gnome_AD_AF, gnome_AD_NFE_AF))
  
  per_sample=read_tsv("bridges_results/per_sample_summary.tsv") 
  
  variants_to_keep = per_sample %>% 
    filter(sampleID %in% sample_list$BRIDGES_ID) %>% 
    dplyr::select(variant) %>% 
    unique()
  
  # Read the list of genes targeted for the analysis
  gene_list=read_tsv("bridges_gene_list.txt", col_names = F)
}

# GET SYNONYMOUS GENE NAMES
dbCon=org.Hs.eg_dbconn()
# write your SQL query
sqlQuery='SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol=dbGetQuery(dbCon, sqlQuery)
# subset to get your results
gene_list=aliasSymbol %>% dplyr::select(alias_symbol, symbol) %>% filter(symbol %in% gene_list$X1)

# Run the script that identifies all protein coding transcripts of the targeted genes
source("~/Dev/familial_bc/scripts/get_canonical_transcript.R")

# Read the list of known breast cancer genes
bc_genes=read_tsv("known_BC_genes.txt", col_names = F)

# Filter the results to only include transcripts if the gene was targeted 
transcripts = per_transcript %>%
  filter(gene %in% gene_list$alias_symbol)

# Filter the results and calculate p-values and FDR from phyloP scores.
variants = per_variant %>%
  separate_longer_delim(cols=gene, delim="|") %>%
  filter(gene %in% gene_list$alias_symbol) %>%
  mutate(p=10**-abs(phyloP)) %>%
  mutate(FDR=p.adjust(p, method="BH", n=2861327195)) # Correcting for all non-N bases in hg19

# Subset only the synonymous variants, determine if phyloP scores are significant
synonymous_variants = variants %>%
  filter(str_detect(type, "synonymous")) %>%
  mutate(phyloP_adj_sig = if_else(FDR<=0.05, TRUE, FALSE)) %>% 
  mutate(phyloP_raw_sig = if_else(p<0.05, TRUE, FALSE))

# Subset only the synonymous variants, remove transcript if NMD transcript variant, filter non-protein coding transcripts
synonymous_transcripts = transcripts %>%
  filter(str_detect(variant_type, "synonymous")) %>%
  filter(!str_detect(variant_type, "NMD_transcript_variant")) %>% 
  filter(transcript_id %in% canonical_transcripts$ensembl_transcript_id)

# Select wanted columns from the synonymous transcript dataframe 
synonymous_transcripts = synonymous_transcripts %>%
  dplyr::select(c(variant, transcript_id, variant_type, codon, encode, rscu, ref_ese, alt_ese, ref_ess, alt_ess, miRNA_target))

# Join the per_transcript data to per_variant data - reformat and parse ClinVar info
synonymous_variants = synonymous_transcripts %>%
  left_join(synonymous_variants %>% dplyr::select(-type), by = join_by(variant==variant)) %>%
  mutate(type=variant_type) %>%
  filter(!is.na(phyloP)) %>% # Remove variants that are not SNVs
  mutate(pathogenicity=if_else(!is.na(ClinVar), str_split(ClinVar, "\\|", simplify = TRUE)[,7], "Unreported")) %>% 
  mutate(clinvar_phenotype=if_else(!is.na(ClinVar),str_split(ClinVar, "\\|", simplify=TRUE)[,14], "Unreported")) %>% 
  mutate(clinvar_review_status=if_else(!is.na(ClinVar),str_split(ClinVar, "\\|", simplify=TRUE)[,25], "Unreported")) %>% 
  mutate(phenotype_BOC=if_else(str_detect(clinvar_phenotype, 
                              "Breast|breast|Ovarian|ovarian|Hereditary_cancer-predisposing_syndrome"), TRUE, FALSE)) %>% 
  mutate(clinvar_id=if_else(!is.na(ClinVar), str_split(ClinVar, "\\|", simplify=TRUE)[,3], "Unreported"))

# Collapse the variant list, so that redundant entries are removed
synonymous_variants_collapsed = synonymous_variants %>% 
  group_by(across(-c(transcript_id, codon, type))) %>% 
  summarise(transcript_id=paste(transcript_id, collapse="|"),
            codon=paste(codon, collapse="|"),
            type=paste(type, collapse="|")) %>% 
            relocate(any_of(c("gene", "known_variation", "type", "AF", "EUR_AF", 
                              "SweAF", "rscu", "phyloP", "phyloP_sig", "ClinVar")), 
                     .after=variant) %>% 
  ungroup()

# Calculate the quantiles for the dRSCU values and select wanted quantiles
rscu_quantiles= quantile(synonymous_variants_collapsed$rscu, probs=seq(0,1,0.05)) %>% 
  as.data.frame() %>% 
  rownames_to_column()

colnames(rscu_quantiles)=c("quantiles", "rscu")

rscu_quantiles = rscu_quantiles %>% 
  filter(quantiles %in% c("5%", "10%", "20%", "30%", "50%", "70%", "80%", "90%", "95%"))
  
# Generate a density plot of the dRSCU values and show the quantiles
ggplot(synonymous_variants_collapsed) + 
  geom_density(aes(x=rscu)) + 
  geom_vline(data=rscu_quantiles, aes(xintercept=rscu), linetype=2, color="grey")

synonymous_variants_final=synonymous_variants_collapsed 

# Select pathogenic variants from list of all synonymous variants - pathogenic in breast/ovarian cancer
pathogenic_BOC_variants = synonymous_variants_final %>% 
  filter((pathogenicity=="Pathogenic" & phenotype_BOC) | (pathogenicity=="Likely_pathogenic" & phenotype_BOC) | (pathogenicity=="Pathogenic/Likely_pathogenic" & phenotype_BOC))

samples_with_pathogenic = per_sample %>% 
  filter(variant %in% pathogenic_BOC_variants$variant) %>% 
  group_by(sampleID) %>% 
  summarise(n=n())

# Select non-pathogenic variants and VUS from list of all synonymous variants - non-pathogenic in breast/ovarian cancer
synonymous_variants_final = synonymous_variants_final %>%
  filter(!variant %in% pathogenic_BOC_variants$variant)

synonymous_variants_final = synonymous_variants_final %>% 
  mutate(ESS=case_when(is.na(ref_ess) & !is.na(alt_ess)~"Create",
                       !is.na(ref_ess) & is.na(alt_ess)~"Remove",
                       !is.na(ref_ess) & !is.na(alt_ess)~"Alter",
                       is.na(ref_ess) & is.na(alt_ess)~NA)) %>% 
  mutate(ESE=case_when(is.na(ref_ese) & !is.na(alt_ese)~"Create",
                       !is.na(ref_ese) & is.na(alt_ese)~"Remove",
                       !is.na(ref_ese) & !is.na(alt_ese)~"Alter",
                       is.na(ref_ese) & is.na(alt_ese)~NA)) %>% 
  mutate(ESS_ESE_overlap=if_else(!is.na(ESS) | !is.na(ESE), TRUE, FALSE)) %>% 
  mutate(known=if_else(!is.na(known_variation), TRUE, FALSE)) %>% 
  mutate(AF_reported=if_else(!is.na(AF) & !is.na(EUR_AF) & !is.na(Swe_AF), TRUE, FALSE)) %>% 
  mutate(splice_region_overlap=if_else(str_detect(type, "splice_region_variant"), TRUE, FALSE)) %>% 
  mutate(RBS_overlap=if_else(!is.na(encode), TRUE, FALSE))

# Rank variants based on chosen thresholds 
synonymous_variants_final = synonymous_variants_final %>% 
  mutate(clinvar_review_score=case_when(clinvar_review_status=="Unreported"~0,
                                        clinvar_review_status=="practice guideline"~4,
                                        clinvar_review_status=="reviewed_by_expert_panel"~3,
                                        clinvar_review_status=="criteria_provided,_multiple_submitters,_no_conflicts"~2,
                                        clinvar_review_status=="criteria_provided,_conflicting_interpretations"~1,
                                        clinvar_review_status=="criteria_provided,_single_submitter"~1,
                                        clinvar_review_status=="no_assertion_criteria_provided"~0,
                                        clinvar_review_status=="no_assertion_provided"~0)) %>% 
  mutate(clinvar_rank=case_when(pathogenicity=="Unreported"~0,
                                pathogenicity=="not_provided"~0,
                                pathogenicity=="Conflicting_interpretations_of_pathogenicity"~1,
                                pathogenicity=="Uncertain_significance"~1,
                                (pathogenicity=="Pathogenic" | pathogenicity=="Likely_pathogenic" | pathogenicity=="Pathogenic/Likely_pathogenic")~2,
                                (pathogenicity=="Benign" | pathogenicity=="Likely_benign" | pathogenicity=="Benign/Likely_benign") & !phenotype_BOC~3,
                                (pathogenicity=="Benign" | pathogenicity=="Likely_benign" | pathogenicity=="Benign/Likely_benign") & phenotype_BOC~4)) %>% 
  mutate(rscu_rank=case_when(rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="5%"]~1,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="10%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="5%"]~2,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="20%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="10%"]~3,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="30%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="20%"]~4,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="50%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="30%"]~5,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="70%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="50%"]~5,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="80%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="70%"]~4,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="90%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="80%"]~3,
                             rscu<=rscu_quantiles$rscu[rscu_quantiles$quantiles=="95%"] & rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="90%"]~2,
                             rscu>rscu_quantiles$rscu[rscu_quantiles$quantiles=="95%"]~1)) %>% 
  mutate(risk_gene_rank=if_else(gene %in% bc_genes$X1, -3, 0)) %>% 
  mutate(ESS_ESE_rank=if_else(ESS_ESE_overlap, 0, 1)) %>% 
  mutate(known_rank=if_else(known, 1, 0)) %>% 
  mutate(AF_rank=if_else(AF_reported, 1+AF+EUR_AF+Swe_AF, 0)) %>% 
  mutate(splice_region_rank=if_else(splice_region_overlap, 0, 3)) %>% 
  mutate(RBS_rank=if_else(RBS_overlap, 0, 1)) %>%
  mutate(phyloP_rank=if_else(phyloP_raw_sig & phyloP>0, 0, 1))

if(data_type=="swea"){
  synonymous_variants_final = synonymous_variants_final %>% 
    separate_longer_delim(samples, "|") %>% 
    mutate(with_pathogenic=if_else(samples %in% samples_with_pathogenic$sampleID, 1, 0)) %>% 
    group_by(across(-c(samples, with_pathogenic))) %>% 
    summarise(samples=paste(samples, collapse="|"),
              with_pathogenic=sum(with_pathogenic)) %>% 
    ungroup() %>% 
    mutate(n_samples = str_count(samples, "\\|") + 1) %>% 
    unique()
  
  synonymous_variants_final = synonymous_variants_final %>% 
    mutate(rank=ESS_ESE_rank+known_rank+rscu_rank+AF_rank+splice_region_rank+RBS_rank+phyloP_rank+(clinvar_rank*clinvar_review_score)+risk_gene_rank+with_pathogenic) %>% 
    arrange(rank)
    
  # Remove problematic variant -> custom annotation scripts cannot handle multiple alternate alleles at one location
  # only one such synonymous variant was found in the SWEA data, none in the BRIDGES data.
  synonymous_variants_final = synonymous_variants_final %>%
    filter(!(variant=="6-35423662-A-C,T" & codon=="ccA/ccC"))

  synonymous_variants_final$variant[synonymous_variants_final$variant=="6-35423662-A-C,T"]="6-35423662-A-T"
  # Export list of variatns as TSV file
  # write_tsv(synonymous_variants_final, file="swea_synonymous_list.tsv")
  # write_tsv(pathogenic_BOC_variants, file="swea_pathogenic_list.tsv")
}
if(data_type=="bridges"){
  synonymous_variants_final = synonymous_variants_final %>%
    separate_longer_delim(samples, "|") %>% 
    mutate(n_cases=if_else(samples %in% case_phenotypes$BRIDGES_ID, 1, 0)) %>% 
    mutate(n_controls=if_else(samples %in% control_phenotypes$BRIDGES_ID, 1, 0)) %>% 
    mutate(with_pathogenic=if_else(samples %in% samples_with_pathogenic$sampleID, 1, 0)) %>%
    group_by(across(-c(samples, n_cases, n_controls, with_pathogenic))) %>% 
    summarise(samples=paste(samples, collapse="|"),
              n_cases=sum(n_cases),
              n_controls=sum(n_controls),
              with_pathogenic=sum(with_pathogenic)) %>% 
    ungroup() %>% 
    mutate(n_samples=n_cases + n_controls) %>% 
    unique()
  
  synonymous_variants_final = synonymous_variants_final %>% 
    filter(variant %in% variants_to_keep$variant) %>% 
    mutate(case_rank = ((n_cases-n_controls)/n_samples)*100) %>% 
    mutate(rank=ESS_ESE_rank+known_rank+rscu_rank+AF_rank+splice_region_rank+RBS_rank+phyloP_rank+(clinvar_rank*clinvar_review_score)+risk_gene_rank-case_rank+with_pathogenic) %>% 
    arrange(rank)

  # Export list of variatns as TSV file
  # write_tsv(synonymous_variants_final, file="bridges_synonymous_list.tsv")
  # write_tsv(pathogenic_BOC_variants, file="bridges_pathogenic_list.tsv")
}
