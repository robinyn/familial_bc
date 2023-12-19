library(tidyverse)
library(UpSetR)

setwd("~/Desktop/thesis")

swea_results = read_tsv("swea_results/swea_synonymous_list.tsv")
bridges_results = read_tsv("bridges_results/bridges_synonymous_list.tsv")

hp_genes = read_tsv("known_BC_genes.txt", col_names = F)

bridges_results$variant = str_remove_all(bridges_results$variant, "chr")

variants_both_cohort = bridges_results$variant[bridges_results$variant %in% swea_results$variant] 

results_both = bridges_results %>% 
  filter(variant %in% variants_both_cohort)

bridges_splice_region = bridges_results %>% 
  filter(str_detect(type, "splice"))

bridges_ESS_ESE = bridges_results %>% 
  filter(ESS_ESE_overlap)

bridges_clinvar_reported = bridges_results %>% 
  filter(pathogenicity!="Unreported")

bridges_RBS_overlap = bridges_results %>% 
  filter(RBS_overlap)

bridges_conserved = bridges_results %>% 
  filter(phyloP_raw_sig & phyloP>0)

bridges_dRSCU = bridges_results %>% 
  filter(rscu_rank==1)

bridges_hp = bridges_results %>% 
  filter(gene %in% hp_genes$X1)

bridges_matrix = list(conserved = bridges_conserved$variant %>% unique(),
                      splice_region = bridges_splice_region$variant %>% unique(), 
                      ESS_ESE_overlap = bridges_ESS_ESE$variant %>% unique(), 
                      RBS_overlap = bridges_RBS_overlap$variant %>% unique(),
                      dRSCU_high = bridges_dRSCU$variant %>% unique(),
                      ClinVar_reported = bridges_clinvar_reported$variant %>% unique())

upset(fromList(bridges_matrix),
      sets = c("conserved", "splice_region", "dRSCU_high", "ESS_ESE_overlap", "RBS_overlap", "ClinVar_reported"),
      keep.order=TRUE,
      group.by="degree",
      set_size.show=TRUE,
      set_size.scale_max = 8000,
      sets.x.label = "Number of variants",
      mainbar.y.label = "Number of variants",
      nintersects = 50,
      main.bar.color="gray40")


swea_splice_region = swea_results %>% 
  filter(str_detect(type, "splice"))

swea_ESS_ESE = swea_results %>% 
  filter(ESS_ESE_overlap)

swea_clinvar_reported = swea_results %>% 
  filter(pathogenicity!="Unreported")

swea_RBS_overlap = swea_results %>% 
  filter(RBS_overlap)

swea_conserved = swea_results %>% 
  filter(phyloP_raw_sig & phyloP>0)

swea_dRSCU = swea_results %>% 
  filter(rscu_rank==1)

swea_hp = swea_results %>% 
  filter(gene %in% hp_genes$X1)

swea_matrix = list(conserved = swea_conserved$variant %>% unique(),
                      splice_region = swea_splice_region$variant %>% unique(), 
                      ESS_ESE_overlap = swea_ESS_ESE$variant %>% unique(), 
                      RBS_overlap = swea_RBS_overlap$variant %>% unique(),
                      dRSCU_high = swea_dRSCU$variant %>% unique(),
                      ClinVar_reported = swea_clinvar_reported$variant %>% unique())

upset(fromList(swea_matrix),
      sets = c("conserved", "splice_region", "dRSCU_high", "ESS_ESE_overlap", "RBS_overlap", "ClinVar_reported"),
      keep.order=TRUE,
      group.by="degree",
      set_size.show=TRUE,
      set_size.scale_max = 700,
      sets.x.label = "Number of variants",
      mainbar.y.label = "Number of variants",
      nintersects=50,
      main.bar.color="gray40")

swea_rscu_quantiles= quantile(swea_results$rscu, probs=seq(0,1,0.05)) %>% 
  as.data.frame() %>% 
  rownames_to_column()

colnames(swea_rscu_quantiles)=c("quantiles", "rscu")

swea_rscu_quantiles = swea_rscu_quantiles %>% 
  filter(quantiles %in% c("5%", "95%"))

bridges_rscu_quantiles= quantile(bridges_results$rscu, probs=seq(0,1,0.05)) %>% 
  as.data.frame() %>% 
  rownames_to_column()

colnames(bridges_rscu_quantiles)=c("quantiles", "rscu")

bridges_rscu_quantiles = bridges_rscu_quantiles %>% 
  filter(quantiles %in% c("5%", "95%"))

ggplot(swea_results) + 
  geom_density(aes(x=rscu)) + 
  geom_vline(aes(xintercept=mean(rscu)), color="grey") +
  geom_vline(data=swea_rscu_quantiles, aes(xintercept=rscu), color="grey", linetype=2) +
  annotate("text", x=0.1, y=0.8, label="Mean") +
  annotate("text", x=1.3, y=0.8, label="95th percentile") +
  annotate("text", x=-1.05, y=0.8, label="5th percentile") +
  xlab("dRSCU") + 
  ylab("Density")

ggplot(bridges_results) + 
  geom_density(aes(x=rscu)) + 
  geom_vline(aes(xintercept=mean(rscu)), color="grey") +
  geom_vline(data=bridges_rscu_quantiles, aes(xintercept=rscu), color="grey", linetype=2) +
  annotate("text", x=0.125, y=1.15, label="Mean") +
  annotate("text", x=1.3, y=1.15, label="95th percentile") +
  annotate("text", x=-0.85, y=1.15, label="5th percentile") +
  xlab("dRSCU") +
  ylab("Density")


swea_df = data.frame(raw=c(46322, 951, 170, 58, 29, 333, 231, 96))

swea_df = swea_df %>%
  mutate(total_perc=raw/46322*100) %>%
  mutate(synonymous_perc=raw/951*100)


bridges_df = data.frame(raw=c(103653, 11919, 3494, 5190, 448, 4556, 3579, 1158))

bridges_df = bridges_df %>%
  mutate(total_perc=raw/103653*100) %>%
  mutate(synonymous_perc=raw/11919*100)


