library(tidyverse)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/thesis")

per_transcript=read_tsv("bridges_results/per_transcript_summary.tsv")
per_variant=read_tsv("bridges_results/per_variant_summary.tsv") %>%
  dplyr::select(-c(gnome_AD_AF, gnome_AD_NFE_AF))

gene_list=read_tsv("bridges_gene_list.txt", col_names = F)

bc_genes=read_tsv("known_BC_genes.txt", col_names = F)

source("~/Dev/familial_bc/scripts/get_canonical_transcript.R")

transcripts = per_transcript %>%
  filter(gene %in% gene_list$X1)

variants = per_variant %>%
  separate_longer_delim(cols=gene, delim="|") %>%
  filter(gene %in% gene_list$X1) %>%
  mutate(p=10**-abs(phyloP)) %>%
  mutate(FDR=p.adjust(p, method="BH", n=2861327195))

synonymous_variants = variants %>%
  filter(str_detect(type, "synonymous")) %>%
  filter(!is.na(phyloP)) %>% # Remove synonymous variants that are not SNVs
  mutate(phyloP_adj_sig = if_else(FDR<=0.05, TRUE, FALSE)) %>% 
  mutate(phyloP_raw_sig = if_else(p<0.05, TRUE, FALSE))

synonymous_transcripts = transcripts %>%
  filter(str_detect(variant_type, "synonymous")) %>%
  filter(!str_detect(variant_type, "NMD_transcript_variant")) %>% 
  filter(transcript_id %in% protein_coding$ensembl_transcript_id)

synonymous_transcripts = synonymous_transcripts %>%
  dplyr::select(c(variant, transcript_id, variant_type, codon, encode, rscu, ref_ese, alt_ese, ref_ess, alt_ess, miRNA_target))

synonymous_variants = synonymous_transcripts %>%
  left_join(synonymous_variants %>% dplyr::select(-type), by = join_by(variant==variant)) %>%
  mutate(type=variant_type) %>%
  dplyr::select(-variant_type) %>% 
  mutate(pathogenicity=if_else(!is.na(ClinVar), str_split(ClinVar, "\\|", simplify = TRUE)[,7], NA)) %>% 
  mutate(clinvar_phenotype=str_split(ClinVar, "\\|", simplify=TRUE)[,14]) %>% 
  mutate(phenotype_BOC=if_else(str_detect(clinvar_phenotype, 
                              "Breast|breast|Ovarian|ovarian|Hereditary_cancer-predisposing_syndrome"), TRUE, FALSE)) %>% 
  mutate(clinvar_id=str_split(ClinVar, "\\|", simplify=TRUE)[,3])

synonymous_variants_collapsed = synonymous_variants %>% 
  group_by(across(-c(transcript_id, codon, type))) %>% 
  summarise(transcript_id=paste(transcript_id, collapse="|"),
            codon=paste(codon, collapse="|"),
            type=paste(type, collapse="|")) %>% 
            relocate(any_of(c("gene", "known_variation", "type", "AF", "EUR_AF", 
                              "SweAF", "rscu", "phyloP", "phyloP_sig", "ClinVar")), 
                     .after=variant) %>% 
  ungroup()

rscu_quantiles= quantile(synonymous_variants_collapsed$rscu, probs=seq(0,1,0.05)) %>% 
  as.data.frame() %>% 
  rownames_to_column()

colnames(rscu_quantiles)=c("quantiles", "rscu")

rscu_quantiles = rscu_quantiles %>% 
  filter(quantiles %in% c("5%", "10%", "20%", "30%", "50%", "70%", "80%", "90%", "95%"))
  
ggplot(synonymous_variants_collapsed) + 
  geom_density(aes(x=rscu)) + 
  geom_vline(data=rscu_quantiles, aes(xintercept=rscu), linetype=2, color="grey")

synonymous_variants_final=synonymous_variants_collapsed 

synonymous_variants_final$pathogenicity[is.na(synonymous_variants_final$pathogenicity)]="Unreported"

pathogenic_BOC_variants = synonymous_variants_final %>% 
  filter((pathogenicity=="Pathogenic" & phenotype_BOC) | 
           (pathogenicity=="Likely_pathogenic" & phenotype_BOC) | 
           (pathogenicity=="Likely_pathogenic" & phenotype_BOC))

synonymous_variants_final = synonymous_variants_final %>% 
  filter(!(pathogenicity=="Pathogenic" & phenotype_BOC)) %>% 
  filter(!(pathogenicity=="Likely_pathogenic" & phenotype_BOC)) %>% 
  filter(!(pathogenicity=="Pathogenic/Likely_pathogenic" & phenotype_BOC)) 

synonymous_variants_final = synonymous_variants_final %>% 
  mutate(ESS=case_when(is.na(ref_ess) & !is.na(alt_ess)~"Create",
                       !is.na(ref_ess) & is.na(alt_ess)~"Remove",
                       is.na(ref_ess) & is.na(alt_ess)~NA)) %>% 
  mutate(ESE=case_when(is.na(ref_ese) & !is.na(alt_ese)~"Create",
                       !is.na(ref_ese) & is.na(alt_ese)~"Remove",
                       is.na(ref_ese) & is.na(alt_ese)~NA)) %>% 
  mutate(ESS_ESE_overlap=if_else(!is.na(ESS) | !is.na(ESE), TRUE, FALSE)) %>% 
  mutate(known=if_else(!is.na(known_variation), TRUE, FALSE)) %>% 
  mutate(AF_reported=if_else(!is.na(AF) & !is.na(EUR_AF) & !is.na(Swe_AF), TRUE, FALSE)) %>% 
  mutate(splice_region_overlap=if_else(str_detect(type, "splice_region_variant"), TRUE, FALSE)) %>% 
  mutate(RBS_overlap=if_else(!is.na(encode), TRUE, FALSE)) %>% 
  mutate(n_samples=str_count(samples, "\\|")+1)

synonymous_variants_final = synonymous_variants_final %>% 
  mutate(clinvar_rank=case_when(pathogenicity=="Pathogenic" | pathogenicity=="Likely_pathogenic"~1,
                                pathogenicity=="Conflicting_interpretations_of_pathogenicity"~2,
                                pathogenicity=="Unreported"~2,
                                pathogenicity=="Uncertain_significance"~2,
                                pathogenicity=="Benign" | pathogenicity=="Likely_benign" | pathogenicity=="Benign/Likely_benign"~3)) %>% 
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
  mutate(ESS_ESE_rank=if_else(ESS_ESE_overlap, 0, 1)) %>% 
  mutate(known_rank=if_else(known, 1, 0)) %>% 
  mutate(AF_rank=if_else(AF_reported, 1+AF+EUR_AF+Swe_AF, 0)) %>% 
  mutate(splice_region_rank=if_else(splice_region_overlap, 0, 1)) %>% 
  mutate(RBS_rank=if_else(RBS_overlap, 0, 1)) %>%
  #mutate(phyloP_rank=if_else(phyloP_adj_sig & phyloP>0, -phyloP, abs(phyloP))) %>% 
  mutate(phyloP_rank=if_else(phyloP_raw_sig & phyloP>0, 0, 1)) %>% 
  mutate(rank=ESS_ESE_rank+known_rank+rscu_rank+AF_rank+splice_region_rank+RBS_rank+phyloP_rank) %>% 
  arrange(rank)

synonymous_variants_known_genes = synonymous_variants_final %>% filter(gene %in% bc_genes$X1)

# write_tsv(synonymous_variants_final, file="bridges_synonymous_list.tsv")
# write_tsv(synonymous_variants_known_genes, file="bridges_synonymous_known_genes_list.tsv")
