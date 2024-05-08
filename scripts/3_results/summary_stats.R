library(tidyverse)
library(UpSetR)

setwd("~/Desktop/RESULTS")

swea_results = read_tsv("SWEA/swea_synonymous_list.tsv")
swea_pathogenic = read_tsv("SWEA/swea_pathogenic_list.tsv")
swea_af = read_tsv("SWEA/swea_af.tsv")
bridges_results = read_tsv("BRIDGES/bridges_synonymous_list.tsv") 
bridges_pathogenic = read_tsv("BRIDGES/bridges_pathogenic_list.tsv")
bridges_AF = read_tsv("BRIDGES/bridges_af.tsv")

swea_results=swea_results %>% 
  filter(!(variant=="6-35423662-A-C,T" & codon=="ccA/ccC"))

swea_results$variant[swea_results$variant=="6-35423662-A-C,T"]="6-35423662-A-T"

synonymous = (swea_results$variant %>% unique() %>% length()) + swea_pathogenic %>% nrow()
pathogenic = swea_pathogenic %>% nrow()
conserved = swea_results %>% filter(phyloP_rank==0) %>% nrow()
unreported = swea_results %>% filter(!AF_reported, !known, is.na(ClinVar)) %>% nrow()
splice_region = swea_results %>% filter(str_detect(type, "splice")) %>% nrow()
ESS_ESE = swea_results %>% filter(ESS_ESE_overlap) %>% nrow()
RBS = swea_results %>% filter(RBS_overlap) %>% nrow()
RSCU = swea_results %>% filter(rscu_rank==1) %>% nrow()

swea_table = data.frame(category=c("synonymous", "pathogenic", "conserved", "unreported", 
                                   "splice_region", "ESS_ESE", "RBS", "RSCU"),
                        values=c(synonymous, pathogenic, conserved, unreported, splice_region, ESS_ESE, RBS, RSCU))

swea_per_gene = swea_results %>% group_by(gene) %>% summarise(n=n())

synonymous = (bridges_results$variant %>% unique() %>% length()) + bridges_pathogenic %>% nrow()
pathogenic = bridges_pathogenic %>% nrow()
conserved = bridges_results %>% filter(phyloP_rank==0) %>% nrow()
unreported = bridges_results %>% filter(!AF_reported, !known, is.na(ClinVar)) %>% nrow()
splice_region = bridges_results %>% filter(str_detect(type, "splice")) %>% nrow()
ESS_ESE = bridges_results %>% filter(ESS_ESE_overlap) %>% nrow()
RBS = bridges_results %>% filter(RBS_overlap) %>% nrow()
RSCU = bridges_results %>% filter(rscu_rank==1) %>% nrow()

bridges_table = data.frame(category=c("synonymous", "pathogenic", "conserved", "unreported", 
                                   "splice_region", "ESS_ESE", "RBS", "RSCU"),
                        values=c(synonymous, pathogenic, conserved, unreported, splice_region, ESS_ESE, RBS, RSCU))

bridges_per_gene = bridges_results %>% group_by(gene) %>% summarise(n=n())

both = intersect(swea_results$variant, bridges_results$variant %>% str_remove("chr"))

both = bridges_results %>% mutate(variant=str_remove(variant, "chr")) %>% filter(variant %in% both)

synonymous = (both %>% filter(pathogenicity!="Pathogenic/Likely_pathogenic") %>% nrow())
conserved = both %>% filter(phyloP_rank==0) %>% nrow()
unreported = both %>% filter(!AF_reported, !known, is.na(ClinVar)) %>% nrow()
splice_region = both %>% filter(str_detect(type, "splice")) %>% nrow()
ESS_ESE = both %>% filter(ESS_ESE_overlap) %>% nrow()
RBS = both %>% filter(RBS_overlap) %>% nrow()
RSCU = both %>% filter(rscu_rank==1) %>% nrow()

both_table = data.frame(category=c("total", "conserved", "unreported", 
                                      "splice_region", "ESS_ESE", "RBS", "RSCU"),
                           values=c(synonymous, conserved, unreported, splice_region, ESS_ESE, RBS, RSCU))