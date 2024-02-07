library(tidyverse)

setwd("~/Desktop/thesis/association_analysis/")

no_subset_1 = read_tsv("no_path_fam_hist_results.txt")
no_subset_2 = read_tsv("no_path_fam_hist_results_adjusted.txt")

no_subset = no_subset_1 %>% 
  merge(no_subset_2 %>% dplyr::select(-1), by = "SNP")

no_subset_significant = no_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf')    

european_subset_1 = read_tsv("no_path_fam_hist_european_results.txt")
european_subset_2 = read_tsv("no_path_fam_hist_european_results_adjusted.txt")

european_subset = european_subset_1 %>%
  merge(european_subset_2 %>% dplyr::select(-1), by = "SNP")

european_subset_significant = european_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf')




# 
# easian_subset_1 = read_tsv("assoc_no_pathogenic_easian_reformatted.tsv")
# easian_subset_2 = read_tsv("assoc_no_pathogenic_easian_adjusted_reformatted.tsv")
# 
# easian_subset = easian_subset_1 %>% 
#   merge(easian_subset_2 %>% select(-1), by = "SNP")
# 
# easian_subset_significant = easian_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf')  