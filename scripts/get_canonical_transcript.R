library(biomaRt)
library(tidyverse)

ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh=37)

protein_coding=getBM(attributes=c("hgnc_symbol", "ensembl_transcript_id", "transcript_biotype", 
                                  "transcript_gencode_basic"), filters="hgnc_symbol", 
                             values=gene_list$alias_symbol, mart=ensembl) %>% 
  filter(transcript_biotype=="protein_coding")

cds_lengths=getBM(attributes=c("ensembl_transcript_id", "cds_length"), filters="ensembl_transcript_id", 
                     values=protein_coding$ensembl_transcript_id, mart=ensembl)

canonical_transcripts=protein_coding  %>% 
  filter(transcript_gencode_basic==1) %>% 
  left_join(cds_lengths, by=join_by(ensembl_transcript_id==ensembl_transcript_id)) %>% 
  group_by(hgnc_symbol) %>% 
  top_n(1, cds_length)
