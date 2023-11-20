library(tidyverse)
library(ggplot2)

setwd("~/Desktop/thesis")

per_variant=read_tsv("swea_results/per_variant_summary.tsv")
per_transcript=read_tsv("swea_results/per_transcript_summary.tsv")

gene_list=read_tsv("swea_gene_list.tsv", col_names = F)

synonymous = per_variant %>% 
  #separate_longer_delim(type, delim="|") %>% 
  filter(str_detect(type, "synonymous_variant")) %>% 
  filter(gene %in% gene_list$X1) %>% 
  mutate(n_sample=lengths(str_split(samples, "\\|"))) %>% 
  select(-samples) %>% 
  mutate(pathogenicity=str_split(ClinVar, "\\|", simplify = TRUE)[,7])

synonymous$known_variation[!is.na(synonymous$known_variation)]=FALSE
synonymous$known_variation[is.na(synonymous$known_variation)]=TRUE

synonymous$pathogenicity[synonymous$pathogenicity==""]=NA
synonymous$pathogenicity[synonymous$pathogenicity=="Benign" | synonymous$pathogenicity=="Likely_benign" ]="Benign/Likely_benign" 
synonymous$pathogenicity[synonymous$pathogenicity=="Pathogenic"]="Pathogenic/Likely_pathogenic"

synonymous_unknown = synonymous %>% 
  filter(is.na(AF)&is.na(EUR_AF)&is.na(Swe_AF)&is.na(known_variation)&is.na(ClinVar))

synonymous_selected = synonymous %>% 
  filter(AF<0.1 & phyloP>=-log10(0.05))

synonymous_pathogenic = synonymous_selected %>% 
  filter(str_detect(ClinVar, "pathogenic"))

synonymous_not_benign = synonymous_selected %>% 
  filter(str_detect(ClinVar, "Benign|benign", negate=TRUE) | is.na(ClinVar))

synonymous_gene = synonymous %>% 
  group_by(gene) %>% 
  summarise(n=n())

synonymous_transcripts = per_transcript %>% 
  filter(variant %in% synonymous_not_benign$variant) %>% 
  filter(str_detect(variant_type,"synonymous_variant"))

p = ggplot() +
  geom_point(data=synonymous, aes(x=phyloP, y=gerp, color="blue", size=AF)) +
  geom_point(data=synonymous_selected, aes(x=phyloP, y=gerp, color="red", size=AF)) +
  scale_x_continuous(limits=c(-10, 10)) +
  scale_y_continuous(limits=c(-10, 10))
print(p)

p = ggplot() +
  geom_point(data=synonymous, aes(x=phyloP, y=gerp, color="blue", size=AF)) +
  geom_point(data=synonymous_pathogenic, aes(x=phyloP, y=gerp, color="red", size=AF)) +
  scale_x_continuous(limits=c(-10, 10)) +
  scale_y_continuous(limits=c(-10, 10))
print(p)

