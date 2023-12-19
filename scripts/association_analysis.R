library(tidyverse)

selectVariants = function(synonymous=T){
  per_transcript=read_tsv("~/Desktop/thesis/bridges_results/per_transcript_summary.tsv")
  per_variant=read_tsv("~/Desktop/thesis/bridges_results/per_variant_summary.tsv") %>%
    dplyr::select(-c(gnome_AD_AF, gnome_AD_NFE_AF))
  
  gene_list=read_tsv("~/Desktop/thesis/bridges_gene_list.txt", col_names = F)
  
  per_transcript = per_transcript %>%
    filter(gene %in% gene_list$X1)
  
  per_variant = per_variant %>%
    separate_longer_delim(cols=gene, delim="|") %>%
    filter(gene %in% gene_list$X1) %>%
    mutate(p=10**-abs(phyloP)) %>%
    mutate(FDR=p.adjust(p, method="BH", n=2861327195))
  
  variants = per_transcript %>%
    dplyr::select(c(variant, transcript_id, variant_type, codon, encode, rscu, 
                    ref_ese, alt_ese, ref_ess, alt_ess, miRNA_target)) %>% 
    left_join(per_variant %>% dplyr::select(-type), by = join_by(variant==variant)) %>% 
    mutate(type=variant_type) %>%
    dplyr::select(-variant_type) %>% 
    mutate(pathogenicity=if_else(!is.na(ClinVar), str_split(ClinVar, "\\|", simplify = TRUE)[,7], NA)) %>% 
    mutate(clinvar_phenotype=str_split(ClinVar, "\\|", simplify=TRUE)[,14]) %>% 
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
  
  variant_list = variants %>% 
    filter(str_detect(type, "synonymous")) %>% 
    dplyr::select(c("variant", "samples"))
}

setwd("~/Desktop/association_analysis/")

control_phenotypes = read_tsv("bridges_annotation/controls_phenotypes.txt")
case_phenotypes = read_tsv("bridges_annotation/cases_phenotypes.txt")

control_samples = control_phenotypes %>% 
  filter(study!="HEBCS") %>% # HEBCS withdrew from BCAC
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only include female samples
  mutate(status=1)

case_samples = case_phenotypes %>% 
  filter(study!="HEBCS") %>% # HEBCS withdrew from BCAC
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only include female samples
  mutate(status=2) 

