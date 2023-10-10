library(tidyverse)
library(biomaRt)
library(ggplot2)
library(GenomeInfoDb)

setwd("~/Desktop/")

results_vcf=read_tsv("./9_ParsedResults/output.tsv")

gene_list=read_tsv("swea_gene_list.tsv", col_names = F)

genome_info = Seqinfo(genome="hg19") %>% 
  as.data.frame() %>% 
  rownames_to_column("chromosome") %>% 
  dplyr::mutate(start=0) %>% 
  dplyr::mutate(end=seqlengths) %>% 
  dplyr::select(chromosome, start, end) %>% 
  dplyr::filter(str_detect(chromosome, "_", negate=TRUE))

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org")

# Retrieve chromosomal coordinates of genes
gene_map = getBM(attributes=c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
              filters='hgnc_symbol',
              values=unique(gene_list$X1),
              mart=mart)

# Remove haplotypes and non-standard chromosomes
gene_map = gene_map %>% 
  filter(chromosome_name=="X" | !is.na(as.numeric(chromosome_name))) %>% 
  arrange(as.numeric(chromosome_name)) %>% 
  mutate(width=end_position-start_position)

p=ggplot(data=gene_map) +
  facet_grid(rows=vars(chromosome_name))

print(p)
