library(tidyverse)
library(ggplot2)

setwd("~/Desktop/thesis")

per_variant=read_tsv("swea_results/per_variant_summary.tsv")
per_transcript=read_tsv("swea_results/reformatted_transcript_summary.tsv")

gene_list=read_tsv("swea_gene_list.tsv", col_names = F)

cnames = colnames(per_variant) 
cnames = cnames[-12]

colnames(per_variant) = cnames

per_variant = per_variant[,1:13]

synonymous = per_variant %>% 
  #separate_longer_delim(type, delim="|") %>% 
  filter(str_detect(type, "synonymous_variant")) %>% 
  filter(gene %in% gene_list$X1) %>% 
  mutate(n_sample=lengths(str_split(samples, "\\|"))) %>% 
  select(-samples)

synonymous$known_variation[!is.na(synonymous$known_variation)]=FALSE
synonymous$known_variation[is.na(synonymous$known_variation)]=TRUE

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
  filter(variant %in% synonymous_selected$variant) %>% 
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

