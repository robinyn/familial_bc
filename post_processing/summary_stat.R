library(tidyverse)
library(biomaRt)
library(ggplot2)
library(GenomeInfoDb)

setwd("~/Desktop/thesis_resources")

#results_vcf=read_tsv("./9_ParsedResults/output.tsv")

summary_file=read_tsv("per_variant_summary.tsv")

gene_list=read_tsv("swea_gene_list.tsv", col_names = F)

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org")

# Retrieve chromosomal coordinates of genes
gene_map = getBM(attributes=c('external_gene_name', 'chromosome_name', 'start_position',
                              'end_position', '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end'),
                 filters='hgnc_symbol',
                 values=unique(gene_list$X1),
                 mart=mart)

# Remove haplotypes and non-standard chromosomes
gene_map = gene_map %>%
  filter(chromosome_name=="X" | !is.na(as.numeric(chromosome_name))) %>%
  arrange(as.numeric(chromosome_name)) %>%
  mutate(width=end_position-start_position)

summary_file = summary_file %>% 
  mutate(num_samples=lengths(str_split(samples, "\\|")))

consequence_types = paste(summary_file$type, collapse = "|") %>% 
  str_split("\\|") %>% 
  unlist() %>%  
  unique()

consequence_summary = c()

for(csq in consequence_types){
  consequence_summary[csq]=summary_file %>% 
    filter(str_detect(type, csq)) %>% 
    nrow()
}

consequence_summary = consequence_summary %>% 
  as.data.frame() %>% 
  rownames_to_column()

colnames(consequence_summary) = c("variant_type", "count")

p = ggplot(data=consequence_summary, aes(x=count, y=variant_type, fill=variant_type)) +
  geom_bar(stat="identity") +
  labs(x="Count", y="Consequence") +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        legend.position = "none")

print(p)

# p = ggplot(data=consequence_summary, aes(x="", y=count, fill=variant_type)) +
#   geom_col() +
#   coord_polar(theta = "y")
# 
# print(p)

total_var = summary_file %>% nrow()

summary_file = summary_file %>% 
  separate_longer_delim(gene, delim="|") %>% 
  filter(gene %in% gene_list$X1)

by_gene = summary_file %>% group_by(gene) %>% summarize(n=n()) %>% 
  arrange(n)

p = ggplot(data=by_gene, aes(x=n, y=gene)) +
  geom_bar(stat="identity")

print(p)

p = ggplot(data=summary_file, aes(AF)) +
  geom_density(aes(y = after_stat(count)))
print(p)

p = ggplot(data=summary_file, aes(EUR_AF)) +
  geom_density(aes(y = after_stat(count)))
print(p)

p = ggplot(data=summary_file, aes(Swe_AF)) +
  geom_density(aes(y = after_stat(count)))
print(p)

filtered_var = summary_file %>% nrow()
known_var = summary_file %>% filter(!is.na(known_variation)) %>% nrow()
unknown_var = summary_file %>% filter(is.na(known_variation)) %>% nrow()
miRNA = summary_file %>% filter(!is.na(miRNA)) %>% nrow()
clinvar = summary_file %>% filter(!is.na(ClinVar)) %>% nrow()

# sprintf("Number of variants found in screened genes: %d (%.2f%%)", filtered_var, (filtered_var/total_var)*100)
# sprintf("Number of known variants: %d (%.2f%%)", known_var, (known_var/filtered_var)*100)
# sprintf("Number of unknown variants: %d (%.2f%%)", unknown_var, (unknown_var/filtered_var)*100)
# sprintf("Number of variants overlapping predicted miRNA target sites: %d (%.2f%%)", miRNA, (miRNA/filtered_var)*100)
# sprintf("Number of variants reported in ClinVar: %d (%.2f%%)", clinvar, (clinvar/filtered_var)*100)

summary_table = data.frame(category=c("total", "filtered", "known", "miRNA_TS_overlap", "ClinVar_reported"), 
                 counts=c(total_var, filtered_var, known_var, miRNA, clinvar),
                 percentage=c(100, (filtered_var/total_var)*100, (known_var/filtered_var)*100, (miRNA/filtered_var)*100, (clinvar/filtered_var)*100))

summary_plot = summary_table %>% 
  filter(category!="total") %>% 
  dplyr::select(category, counts) %>% 
  mutate(percentage=counts/filtered_var*100)

p = ggplot(data=summary_plot, aes(x=category, y=percentage)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c("filtered", "known", "miRNA_TS_overlap", "ClinVar_reported"))

print(p)

print(summary_table)