library(tidyverse)

setwd("~/Desktop/association_analysis/")

per_transcript=read_tsv("~/Desktop/thesis/bridges_results/per_transcript_summary.tsv") %>% 
  mutate(variant=str_remove(variant, "chr"))
per_variant=read_tsv("~/Desktop/thesis/bridges_results/per_variant_summary.tsv") %>%
  dplyr::select(-c(gnome_AD_AF, gnome_AD_NFE_AF)) %>% 
  mutate(variant=str_remove(variant, "chr"))

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

variant_list = variants %>% 
  filter(str_detect(type, "synonymous")) %>% 
  mutate(ref_allele=str_split(variant, "-", simplify=TRUE)[,3]) %>% 
  mutate(alt_allele=str_split(variant, "-", simplify=TRUE)[,4]) %>% 
  filter(str_length(ref_allele)==1 & str_length(alt_allele)==1) %>% 
  dplyr::select(c("variant", "samples", "ref_allele", "alt_allele")) %>% 
  unique()

synonymous_variants = variants %>% filter(variant %in% variant_list$variant) %>% unique()

control_phenotypes = read_tsv("bridges_annotation/controls_phenotypes.txt")
case_phenotypes = read_tsv("bridges_annotation/cases_phenotypes.txt")

control_samples = control_phenotypes %>% 
  filter(study!="HEBCS") %>% # HEBCS withdrew from BCAC
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only includes female samples
  mutate(status=1)

case_samples = case_phenotypes %>% 
  filter(study!="HEBCS") %>% # HEBCS withdrew from BCAC
  dplyr::select(c(BRIDGES_ID)) %>% 
  mutate(sex=2) %>% # Data only includes female samples
  mutate(status=2) 

sample_list = control_samples %>% rbind(case_samples)

per_sample = read_tsv("~/Desktop/thesis/bridges_results/per_sample_summary.tsv")
per_sample = per_sample %>% 
  mutate(variant=str_remove(variant, "chr")) %>% 
  filter(variant %in% variant_list$variant) %>% 
  filter(sampleID %in% sample_list$BRIDGES_ID)

ped_out = per_sample %>% pivot_wider(names_from = variant, values_from = genotype)

ped_out = left_join(sample_list, ped_out, by=join_by("BRIDGES_ID"=="sampleID"))

map_out = colnames(ped_out)[-c(1,2,3)] %>% 
  as.data.frame() %>% 
  setNames("variantID") %>% 
  mutate(chromosome=str_split(variantID, "-", simplify = TRUE)[,1]) %>% 
  mutate(bp_position=str_split(variantID, "-", simplify = TRUE)[,2]) %>% 
  dplyr::select(c(chromosome, variantID, bp_position))

samples_pathogenic = per_sample %>% 
  left_join(variants %>% select(variant, pathogenicity) %>% unique(), by=join_by(variant==variant))

samples_pathogenic = samples_pathogenic %>% 
  filter(pathogenicity=="Pathogenic" | pathogenicity=="Likely_pathogenic" | pathogenicity=="Pathogenic/Likely_pathogenic") %>% 
  select(sampleID) %>% 
  unique()

ped_out = ped_out %>% 
  filter(!(BRIDGES_ID %in% samples_pathogenic$sampleID))

# Replace NA with homogeneous reference genotype
ped_out = ped_out %>% 
  mutate(across(colnames(ped_out)[-c(1,2,3)], 
                ~replace_na(., sprintf("%1$s%1$s", str_split(cur_column(), "-", simplify = TRUE)[[3]]))))

write_tsv(ped_out, "assoc_analysis.ped", na="00", col_names = FALSE)
write_tsv(map_out, "assoc_analysis.map", na="00", col_names = FALSE)




