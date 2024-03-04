library(tidyverse)
library(ggplot2)
library(VennDiagram)

setwd("~/Desktop/thesis/association_analysis/")

# Read phenotype annotations for the samples
control_phenotypes = read_tsv("bridges_annotation/controls_phenotypes.txt") %>% 
  filter(study!="HEBCS")
case_phenotypes = read_tsv("bridges_annotation/cases_phenotypes.txt") %>% 
  filter(study!="HEBCS")


# No subset
no_subset_1 = read_tsv("assoc_no_pathogenic_fam_hist.tsv")
no_subset_2 = read_tsv("assoc_no_pathogenic_fam_hist_adjusted.tsv")

no_subset = no_subset_1 %>% 
  merge(no_subset_2 %>% dplyr::select(-1), by = "SNP")

rm(no_subset_1)
rm(no_subset_2)

no_subset_significant = no_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf')    

# European subset
european_subset_1 = read_tsv("assoc_european_no_pathogenic_fam_hist.tsv")
european_subset_2 = read_tsv("assoc_european_no_pathogenic_fam_hist_adjusted.tsv")

european_subset = european_subset_1 %>% 
  merge(european_subset_2 %>% dplyr::select(-1), by = "SNP")

rm(european_subset_1)
rm(european_subset_2)

european_subset_significant = european_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf')    

# East/South Asian subset
asian_subset_1 = read_tsv("assoc_es_asian_no_pathogenic_fam_hist.tsv")
asian_subset_2 = read_tsv("assoc_es_asian_no_pathogenic_fam_hist_adjusted.tsv")

asian_subset = asian_subset_1 %>% 
  merge(asian_subset_2 %>% dplyr::select(-1), by = "SNP")

rm(asian_subset_1)
rm(asian_subset_2)

asian_subset_significant = asian_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf') 

# Under 50
u50_subset_1 = read_tsv("assoc_u50_no_pathogenic_fam_hist.tsv")
u50_subset_2 = read_tsv("assoc_u50_no_pathogenic_fam_hist_adjusted.tsv")

u50_subset = u50_subset_1 %>% 
  merge(u50_subset_2 %>% dplyr::select(-1), by = "SNP")

rm(u50_subset_1)
rm(u50_subset_2)

u50_subset_significant = u50_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf') 

# At least 1 first dg family
fdg_subset_1 = read_tsv("assoc_no_pathogenic_fam_hist_1st_dg.tsv")
fdg_subset_2 = read_tsv("assoc_no_pathogenic_fam_hist_1st_dg_adjusted.tsv")

fdg_subset = fdg_subset_1 %>% 
  merge(fdg_subset_2 %>% dplyr::select(-1), by = "SNP")

rm(fdg_subset_1)
rm(fdg_subset_2)

fdg_subset_significant = fdg_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf') 

# At least 1 first dg family AND u50
u50_fdg_subset_1 = read_tsv("assoc_no_pathogenic_fam_hist_u50_1st_dg.tsv")
u50_fdg_subset_2 = read_tsv("assoc_no_pathogenic_fam_hist_u50_1st_dg_adjusted.tsv")

u50_fdg_subset = u50_fdg_subset_1 %>% 
  merge(u50_fdg_subset_2 %>% dplyr::select(-1), by = "SNP")

rm(u50_fdg_subset_1)
rm(u50_fdg_subset_2)

u50_fdg_subset_significant = u50_fdg_subset %>% filter(FDR_BH<0.05 | FDR_BH=='Inf') 

synonymous_table = read_tsv("~/Desktop/thesis/bridges_results/bridges_synonymous_list.tsv") %>% 
  mutate(variant = str_remove(variant, "chr"))

no_subset_significant = no_subset_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

european_subset_significant = european_subset_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

asian_subset_significant = asian_subset_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

fdg_subset_significant = fdg_subset_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

u50_subset_significant = u50_subset_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

u50_fdg_subset_significant = u50_fdg_subset_significant %>% 
  dplyr::select(c(SNP, A1, A2, F_A, F_U, OR, FDR_BH)) %>% 
  left_join(synonymous_table, by=join_by(SNP==variant))

t1 = no_subset_significant #%>% filter(OR>1)
t2 = european_subset_significant #%>% filter(OR>1)
t3 = asian_subset_significant #%>% filter(OR>1)
t4 = u50_subset_significant
t5 = fdg_subset_significant
t6 = u50_fdg_subset_significant

venn_list=list(No_sub=t1$SNP, Eur_sub=t2$SNP, Asian_sub=t3$SNP)
x = venn.diagram(venn_list, filename=NULL)
grid.draw(x)

all_three=Reduce(intersect, list(t1$SNP,t6$SNP))

# 
# samples= case_phenotypes %>% filter(brca12==1)
# 
# samples %>% 
#   filter(ER_statusIndex!=888 & PR_statusIndex!=888 & HER2_statusIndex!=888) %>% 
#   group_by(ER_statusIndex, PR_statusIndex, HER2_statusIndex) %>% 
#   summarise(n=n())
# 
# samples %>% 
#   filter(ER_statusIndex==888 & PR_statusIndex==888 &HER2_statusIndex==888) %>% nrow()
# 
# samples %>% group_by(MorphologygroupIndex_corr) %>% summarise(n=n())
# 
# samples %>% group_by(StageIndex) %>% summarise(n=n())
# 
# samples %>% group_by(GradeIndex) %>% summarise(n=n())
