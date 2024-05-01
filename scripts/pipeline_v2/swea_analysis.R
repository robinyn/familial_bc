library(tidyverse)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/thesis")

per_transcript=read_tsv("swea_results/per_transcript_summary.tsv")
per_variant=read_tsv("swea_results/per_variant_summary.tsv") %>% 
  dplyr::select(-c(gnome_AD_AF, gnome_AD_NFE_AF))

gene_list=read_tsv("swea_gene_list.tsv", col_names = F)

transcripts = per_transcript %>% 
  filter(gene %in% gene_list$X1)

variants = per_variant %>% 
  separate_longer_delim(cols=gene, delim="|") %>% 
  filter(gene %in% gene_list$X1) %>% 
  mutate(p=10**-abs(phyloP)) %>% 
  mutate(FDR=p.adjust(p, method="BH", n=2861327195))

synonymous_variants = variants %>% 
  filter(str_detect(type, "synonymous")) %>% 
  mutate(phyloP_sig = if_else(FDR<=0.05, TRUE, FALSE))

synonymous_transcripts = transcripts %>% 
  filter(str_detect(variant_type, "synonymous"))

synonymous_transcripts = synonymous_transcripts %>% 
  dplyr::select(c(variant, codon, encode, rscu, ref_ese, alt_ese, ref_ess, alt_ess, miRNA_target))

synonymous_variants = synonymous_variants %>% 
  left_join(synonymous_transcripts, by = join_by(variant==variant))

synonymous_significant = synonymous_variants %>% 
  filter(FDR<=0.05)

geneset1 = gene_list$X1[1:(length(gene_list$X1)/2)]
geneset2 = gene_list$X1[((length(gene_list$X1)/2)+1):length(gene_list$X1)]

ggplot() + 
  geom_violin(data=synonymous_variants %>% 
                filter(gene %in% geneset1) %>% 
                filter(str_detect(type, "synonymous")), 
              aes(x=factor(gene), y=phyloP)) +
  geom_jitter(data=synonymous_variants %>% filter(gene %in% geneset1), 
             aes(x=factor(gene), y=phyloP, color=phyloP_sig)) +
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot() + 
  geom_violin(data=synonymous_variants %>% 
                filter(gene %in% geneset2) %>% 
                filter(str_detect(type, "synonymous")), 
              aes(x=factor(gene), y=phyloP)) +
  geom_jitter(data=synonymous_variants %>% filter(gene %in% geneset2), 
              aes(x=factor(gene), y=phyloP, color=phyloP_sig)) +
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot() + 
  geom_violin(data=synonymous_variants %>% 
                filter(gene %in% geneset1), 
              aes(x=factor(gene),y=rscu)) +
  stat_summary(data=synonymous_variants %>% 
                 filter(gene %in% geneset1),
               aes(group = gene, x=gene, y=rscu), 
               fun = mean, geom = 'point', size=1.5, alpha=0.9) +
  geom_hline(yintercept = 0, linetype=2) +
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot() + 
  geom_violin(data=synonymous_variants %>% 
                filter(gene %in% geneset2), 
              aes(x=factor(gene),y=rscu)) +
  stat_summary(data=synonymous_variants %>% 
                 filter(gene %in% geneset2),
               aes(group = gene, x=gene, y=rscu), 
               fun = mean, geom = 'point', size=1.5, alpha=0.9) +
  geom_hline(yintercept = 0, linetype=2) +
  theme(axis.text.x = element_text(angle=90, hjust=1))


known_variation_data = synonymous_variants %>% 
  dplyr::select(variant, known_variation, gene, AF, EUR_AF, Swe_AF, ClinVar) %>% 
  unique() %>% 
  mutate(known=if_else(!is.na(known_variation), T, F)) %>% 
  mutate(pathogenicity=if_else(!is.na(ClinVar), str_split(ClinVar, "\\|", simplify = TRUE)[,7], NA))

known_variation_data$pathogenicity[known_variation_data$pathogenicity=="Benign" | known_variation_data$pathogenicity=="Likely_benign"] = "Benign/Likely_benign"
known_variation_data$pathogenicity[known_variation_data$pathogenicity=="Pathogenic"] = "Pathogenic/Likely_pathogenic"
known_variation_data$pathogenicity[is.na(known_variation_data$pathogenicity)] = "Unreported"

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

ggplot(data=known_variation_data, aes(x="", fill=known)) +
  geom_bar() +
  coord_polar(theta="y") +
  blank_theme +
  guides(fill=guide_legend("Variant reported previously")) +
  geom_text(stat="count", 
            aes(x=1.6, 
                label=after_stat(scales::percent(..count../sum(..count..), 1))),
            position = position_stack(vjust=0.5),
            size=4)

ggplot(data=known_variation_data, aes(x="", fill=pathogenicity)) +
  geom_bar() +
  coord_polar(theta="y") +
  blank_theme +
  guides(fill=guide_legend("Pathogenicity")) +
  geom_text(stat="count", 
            aes(x=1.6, 
                label=after_stat(scales::percent(..count../sum(..count..), 1))),
            position = position_stack(vjust=0.5),
            size=3)

synonymous_variants_original = synonymous_variants
synonymous_variants = synonymous_variants %>% unique()

ggplot() + 
  geom_boxplot(data=synonymous_variants %>% 
                filter(gene %in% geneset1) %>% 
                filter(str_detect(type, "synonymous")), 
              aes(x=reorder(gene, phyloP, median), y=phyloP),
              outlier.shape = NA) +
  geom_jitter(data=synonymous_variants %>% filter(gene %in% geneset1), 
              aes(x=factor(gene), y=phyloP, color=phyloP_sig, alpha=0.1)) +
  geom_hline(yintercept = -log10(0.05), linetype=2)+
  geom_hline(yintercept = log10(0.05), linetype=2)+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "none") +
  labs(x="Gene")

ggplot() + 
  geom_boxplot(data=synonymous_variants %>% 
                filter(gene %in% geneset2) %>% 
                filter(str_detect(type, "synonymous")), 
              aes(x=reorder(gene, phyloP, median), y=phyloP),
              outlier.shape = NA) +
  geom_jitter(data=synonymous_variants %>% filter(gene %in% geneset2), 
              aes(x=factor(gene), y=phyloP, color=phyloP_sig, alpha=0.1)) +
  geom_hline(yintercept = -log10(0.05), linetype=2)+
  geom_hline(yintercept = log10(0.05), linetype=2)+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "none")+
  labs(x="Gene")

ggplot() + 
  geom_boxplot(data=synonymous_variants %>% 
                filter(gene %in% geneset1), 
              aes(x=reorder(gene,rscu, median),y=rscu)) +
  stat_summary(data=synonymous_variants %>%
                filter(gene %in% geneset1),
              aes(group = gene, x=gene, y=rscu),
              fun = mean, geom = 'point', size=1.5, alpha=0.9) +
  geom_hline(yintercept = 0, linetype=2) +
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot() + 
  geom_boxplot(data=synonymous_variants %>% 
                filter(gene %in% geneset2), 
              aes(x=reorder(gene,rscu,median),y=rscu)) +
  stat_summary(data=synonymous_variants %>%
                filter(gene %in% geneset2),
              aes(group = gene, x=gene, y=rscu),
              fun = mean, geom = 'point', size=1.5, alpha=0.9) +
  geom_hline(yintercept = 0, linetype=2) +
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot(data=synonymous_variants) + 
  geom_freqpoly(aes(x=AF, color="black"), bins=100) + 
  geom_freqpoly(aes(x=EUR_AF, color="red"), bins=100) + 
  geom_freqpoly(aes(x=Swe_AF, color="blue"), bins=100) +
  theme(legend.title=element_blank()) +
  scale_color_discrete(labels=c("AF", "EUR_AF", "Swe_AF"))

synonymous_variants = synonymous_variants %>% 
  mutate(overlaps_ESE_ESS = if_else(!is.na(ref_ese) | !is.na(ref_ess) | !is.na(alt_ese) | !is.na(alt_ess),
                                    TRUE, FALSE)) %>% 
  mutate(overlaps_RBS = if_else(!is.na(encode), TRUE, FALSE))

# Select variants reported at pathogenic/likely pathogenic in ClinVar
pathogenic_variants = per_variant %>% 
  filter(str_detect(ClinVar, "Pathogenic|Likely_pathogenic")) %>% 
  mutate(pathogenicity=str_split(ClinVar, "\\|", simplify = TRUE)[,7]) %>% 
  mutate(clinvar_phenotype=str_split(ClinVar, "\\|", simplify=TRUE)[,14]) %>% 
  mutate(phenotype_BOC=if_else(str_detect(clinvar_phenotype, "Breast|breast|Ovarian|ovarian"), TRUE, FALSE))

pathogenic_variants_BOC = pathogenic_variants %>% filter(phenotype_BOC)

search_string = pathogenic_variants_BOC$samples %>% 
  paste(collapse="\n")

# Generate a list of sample IDs
samples = per_variant %>%
  dplyr::select(samples) %>%
  separate_longer_delim(cols=samples, delim="|") %>%
  unique() %>%
  mutate(p_var=0)

# Count number of pathogenic variants found in each sample
for(sample in samples$samples){
  samples$p_var[samples$samples==sample]=length(grep(sample,
                                                     unlist(str_split(paste(pathogenic_variants_BOC$samples,
                                                                            collapse = "|"), "\\|"))))
}

samples_stat = samples %>% group_by(p_var) %>% summarise(n=n()) %>% 
  mutate(percentage=n/3403*100)

stat_df = data.frame(category = c("ESE_ESS", "ESE_ESS", "RBS", "RBS"), 
                     type=c("total", "values", "total", "values"), 
                     values=c(984, 499, 984, 237))

ggplot(data=stat_df) +
  geom_bar(aes(x=category, y=values, fill=type),
           stat="identity",
           position="identity")

ggplot(samples_stat, aes(label=round(percentage,2))) + 
  geom_bar(aes(x=p_var, y=n, fill=as.factor(p_var)), stat="identity") +
  theme(legend.position="none") +
  geom_text(aes(x=p_var, y=n/2))


synonymous_unknown=synonymous_variants %>% 
  filter(is.na(AF)&is.na(EUR_AF)&is.na(Swe_AF)) %>% 
  filter(is.na(ClinVar)) %>% 
  filter(is.na(known_variation))

splice_region_synonymous_variants = synonymous_variants_original %>% 
  filter(str_detect(type, "splice")) %>% 
  mutate(pathogenicity=str_split(ClinVar, "\\|", simplify=TRUE)[,7])

