library(tidyverse)
library(ggplot2)
library(VennDiagram)

setwd("~/Desktop/RESULTS/BRIDGES/association/basic_association/output")

# Read phenotype annotations for the samples
control_phenotypes = read_tsv("../../../bridges_annotation/controls_phenotypes.txt") %>% 
  filter(study!="HEBCS")
case_phenotypes = read_tsv("../../../bridges_annotation/cases_phenotypes.txt") %>% 
  filter(study!="HEBCS")

# Read variant annotations
synonymous_table = read_tsv("~/Desktop/thesis/bridges_results/bridges_synonymous_list.tsv") %>% 
  mutate(variant = str_remove(variant, "chr"))

# # TEST
# a1 = read_tsv("../../input/output_reformatted.txt")
# a2 = read_tsv("../../input/output_reformatted_adjusted.txt")
# 
# TEST = a1 %>% 
#   filter(TEST=="ADD") %>% 
#   merge(a2 %>% dplyr::select(-1), by = "SNP")
# 
# TEST_sig = TEST %>% filter(FDR_BH<0.05 | FDR_BH=='Inf')
# 
# TEST_sig = TEST_sig %>% 
#   dplyr::select(c(SNP, A1, OR, UNADJ, FDR_BH)) %>% 
#   left_join(synonymous_table, by=join_by(SNP==variant))
# 
# #######################################################

# General breast cancer risk - No subsetting for fam hist
a1 = read_tsv("no_pathogenic_reccurent.tsv")
a2 = read_tsv("no_pathogenic_reccurent_adjusted.tsv")

no_path_recur = a1 %>% 
  merge(a2 %>% dplyr::select(-1), by = "SNP")

no_path_recur_significant = no_path_recur %>% filter(FDR_BH<0.05 | FDR_BH=='Inf')

no_path_recur_significant = no_path_recur_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, UNADJ, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

#####################

# History of BOC in 1st degree family member
a1 = read_tsv("fam_hist_1dg_no_pathogenic_recurrent.tsv")
a2 = read_tsv("fam_hist_1dg_no_pathogenic_recurrent_adjusted.tsv")

fdg_no_path_recur = a1 %>% 
  merge(a2 %>% dplyr::select(-1), by = "SNP")

fdg_no_path_recur_significant = fdg_no_path_recur %>% filter(FDR_BH<0.05 | FDR_BH=='Inf') 

fdg_no_path_recur_significant = fdg_no_path_recur_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, UNADJ, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

#######################

# History of BOC in 1st degree family member AND diagnosis u50
a1 = read_tsv("u50_fam_hist_1dg_no_pathogenic_recurrent.tsv")
a2 = read_tsv("u50_fam_hist_1dg_no_pathogenic_recurrent_adjusted.tsv")

u50_fdg_no_path_recur = a1 %>% 
  merge(a2 %>% dplyr::select(-1), by = "SNP")

u50_fdg_no_path_recur_significant = u50_fdg_no_path_recur %>% filter(FDR_BH<0.05 | FDR_BH=='Inf') 

u50_fdg_no_path_recur_significant = u50_fdg_no_path_recur_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, UNADJ, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

#######################

# History of BOC in 1st degree family member AND multiple family members affected AND diagnosis u50
a1 = read_tsv("u50_multi_fam_hist_1dg_no_pathogenic_recurrent.tsv")
a2 = read_tsv("u50_multi_fam_hist_1dg_no_pathogenic_recurrent_adjusted.tsv")

u50_multi_fdg_no_path_recur = a1 %>% 
  merge(a2 %>% dplyr::select(-1), by = "SNP")

u50_multi_fdg_no_path_recur_significant = u50_multi_fdg_no_path_recur %>% filter(FDR_BH<0.05 | FDR_BH=='Inf') 

u50_multi_fdg_no_path_recur_significant = u50_multi_fdg_no_path_recur_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, UNADJ, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

# 
# variants = u50_multi_fdg_no_path_recur_significant %>%
#   filter(OR > 1) %>% filter(phyloP_sig, phyloP>0)
# 
# sampleIDs = synonymous_table %>%
#   filter(variant %in% variants$SNP) %>%
#   separate_longer_delim(samples, "|") %>%
#   dplyr::select(samples) %>%
#   unique()
# 
# samples = case_phenotypes %>%
#   filter(BRIDGES_ID %in% sampleIDs$samples)
# 
# samples %>%
#   filter(ER_statusIndex!=888 & PR_statusIndex!=888 & HER2_statusIndex!=888) %>%
#   group_by(ER_statusIndex, PR_statusIndex, HER2_statusIndex) %>%
#   summarise(n=n())
# 
# samples %>%
#   filter(ER_statusIndex==888 | PR_statusIndex==888 | HER2_statusIndex==888) %>% nrow()
# 
# samples %>% group_by(MorphologygroupIndex_corr) %>% summarise(n=n())
# 
# samples %>% group_by(StageIndex) %>% summarise(n=n())
# 
# samples %>% group_by(GradeIndex) %>% summarise(n=n())
