library(tidyverse)
library(ggplot2)
library(FSA)

######################## RSCU ###########################

# Generate a density plot of the dRSCU values and show the quantiles
ggplot(synonymous) + 
  geom_density(aes(x=rscu), bw=0.15) + 
  geom_vline(data=rscu_quantiles, aes(xintercept=rscu), linetype=2, color="grey") +
  geom_vline(aes(xintercept=mean(rscu)), color="red", linetype=2)

ggplot(synonymous) + 
  geom_histogram(aes(x=rscu), binwidth = 0.5)

ggplot(synonymous) +
  geom_boxplot(aes(x=reorder(pathogenicity, rscu, mean), y=rscu)) +
  geom_hline(data=rscu_quantiles, aes(yintercept=rscu), linetype=2, color="grey") +
  coord_flip()

ggplot(synonymous) + 
  geom_point(aes(x=rscu, y=phyloP)) + 
  geom_hline(aes(yintercept=-log10(0.05))) +
  geom_hline(aes(yintercept=log10(0.05))) 

kruskal.test(phyloP~pathogenicity, data=synonymous)
dunnTest(phyloP~pathogenicity, data=synonymous, method="bh")

################## Variant types/consequences ########################
# Counts of different variant consequences
variant_types = variants_final %>% 
  separate_longer_delim(cols=type, delim="&") %>% 
  group_by(type) %>% 
  summarise(n=n()) %>% 
  arrange(-n)

# Bar plot of variant types
ggplot(variants_final %>% separate_longer_delim(cols=type, delim="&")) +
  geom_bar(aes(y=fct_rev(fct_infreq(type)), fill=type), stat="count", width=0.9) +
  theme(legend.position = "none") +
  labs(x="Count", y="Variant type")

######################## SAMPLES WITH PATHOGENIC #######################

samples_num_pathogenic = samples_with_pathogenic %>% 
  group_by(n) %>% 
  summarise(num_path=n())

# BRIDGES
samples_num_pathogenic = samples_with_pathogenic %>% 
  filter(sampleID %in% case_samples$BRIDGES_ID) %>% 
  group_by(n) %>% 
  summarise(num_path=n())

######################### CLINVAR ALL VARIANTS #########################

# Generate counts of ClinVar significance 
pathogenicity_table = variants_final %>% 
  group_by(pathogenicity) %>% 
  summarise(n=n())

# Calculate percentage
clinvar_distribution = pathogenicity_table %>% 
  filter(pathogenicity!="Unreported") %>% 
  mutate(percent = n/sum(n)*100)

clinvar_review = variants_final %>% 
  group_by(pathogenicity, clinvar_review_status) %>% 
  summarise(n=n())

# Bar plot of ClinVar clinical significance
ggplot() + 
  geom_bar(data=clinvar_distribution, 
           aes(y=fct_reorder(pathogenicity, percent), x=percent, fill=as.factor(pathogenicity)), stat="identity") + 
  scale_fill_viridis_d() +
  theme(legend.position = "none") +
  labs(y="Clinical significance", x="Percentage")

######################### CLINVAR SYNONYMOUS VARIANTS #########################

# Generate counts of ClinVar significance 
pathogenicity_table = synonymous %>% 
  group_by(pathogenicity) %>% 
  summarise(n=n())

# Calculate percentage
clinvar_distribution = pathogenicity_table %>% 
  filter(pathogenicity!="Unreported") %>% 
  mutate(percent = n/sum(n)*100)

clinvar_review = synonymous %>% 
  group_by(pathogenicity, clinvar_review_status) %>% 
  summarise(n=n())

# Bar plot of ClinVar clinical significance
ggplot() + 
  geom_bar(data=clinvar_distribution, 
           aes(y=fct_reorder(pathogenicity, percent), x=percent, fill=as.factor(pathogenicity)), stat="identity") + 
  scale_fill_viridis_d() +
  theme(legend.position = "none") +
  labs(y="Clinical significance", x="Percentage")

######################### phyloP - variant type #########################

phyloP_mean_per_variant_type = variants_final %>% 
  separate_longer_delim(type, delim="&") %>% 
  group_by(type) %>% 
  summarise(mean_phyloP=mean(phyloP, na.rm=TRUE))

ggplot(variants_final %>% separate_longer_delim(type, delim="&") %>% filter(!is.na(phyloP)),
       aes(x=reorder(type, phyloP, median), y=phyloP)) +
  geom_boxplot(outlier.shape = NA, na.rm=TRUE) +
  geom_hline(yintercept=0, linetype=2, alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5, color="blue") +
  geom_hline(yintercept=log10(0.05), linetype=2, alpha=0.5, color="blue") +
  stat_summary(fun="mean", geom="point", shape=1, size=1, color="red") +
  coord_flip()

mean(variants_final$phyloP[str_detect(variants_final$type, "synonymous")])

phyloP_sig_variants = variants_final %>% filter(str_detect(type, "synonymous") & phyloP_sig)


######################### phyloP - ClinVar #########################

ggplot(variants_final %>% filter(!is.na(phyloP)),
       aes(x=reorder(pathogenicity, phyloP, median), y=phyloP)) +
  geom_boxplot(outlier.shape = NA, na.rm=TRUE) +
  geom_hline(yintercept=0, linetype=2, alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5, color="blue") +
  geom_hline(yintercept=log10(0.05), linetype=2, alpha=0.5, color="blue") +
  stat_summary(fun="mean", geom="point", shape=1, size=1, color="red") +
  coord_flip()

kruskal.test()

# t.test(variants_final$phyloP[variants_final$pathogenicity=="Unreported"], 
#        variants_final$phyloP[variants_final$pathogenicity=="Other"])


######################### RSCU - ClinVar #########################

ggplot(variants_final %>% filter(!is.na(rscu)),
       aes(x=reorder(pathogenicity, rscu, median), y=rscu)) +
  geom_boxplot(outlier.shape = NA, na.rm=TRUE) +
  geom_hline(yintercept=0, linetype=2, alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5, color="blue") +
  geom_hline(yintercept=log10(0.05), linetype=2, alpha=0.5, color="blue") +
  stat_summary(fun="mean", geom="point", shape=1, size=1, color="red") +
  coord_flip()

######################### ESE ESS #########################

ggplot(variants_final %>% filter(!is.na(ESE)), aes(x=ESE)) +
  geom_bar(aes(fill=ESE),stat="count")

ggplot(variants_final %>% filter(!is.na(ESS)), aes(x=ESS)) +
  geom_bar(aes(fill=ESS),stat="count")

ggplot(variants_final %>% filter(!is.na(ESE)), aes(x=ESE)) + 
  geom_bar(aes(fill=pathogenicity), stat="count", position="dodge")

ggplot(variants_final %>% filter(!is.na(ESS)), aes(x=ESS)) + 
  geom_bar(aes(fill=pathogenicity), stat="count", position="dodge")

############## ESE ESS for synonymous #####################

ggplot(synonymous %>% filter(!is.na(ESE)), aes(x=ESE)) +
  geom_bar(aes(fill=ESE),stat="count")

ggplot(synonymous %>% filter(!is.na(ESS)), aes(x=ESS)) +
  geom_bar(aes(fill=ESS),stat="count")

ggplot(synonymous %>% filter(!is.na(ESE)), aes(x=ESE)) + 
  geom_bar(aes(fill=pathogenicity), stat="count", position="dodge")

ggplot(synonymous %>% filter(!is.na(ESS)), aes(x=ESS)) + 
  geom_bar(aes(fill=pathogenicity), stat="count", position="dodge")

######################### RBP binding site #########################

RBP_per_variant = synonymous %>% mutate(n_RBP=str_count(encode, "&")+1) %>% 
  group_by(n_RBP) %>% 
  summarise(n=n())

RBP_list = synonymous %>% separate_longer_delim(encode, "&") %>% 
  mutate(encode=str_split(encode, "_", simplify=TRUE)[,1]) %>% 
  dplyr::select(encode) %>% 
  group_by(encode) %>% 
  summarise(n=n())

RBP_pathogenicity = synonymous %>% separate_longer_delim(encode, "&") %>% 
  mutate(encode=str_split(encode, "_", simplify=TRUE)[,1]) %>% 
  dplyr::select(c(encode, pathogenicity)) %>% 
  group_by(encode, pathogenicity) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from=pathogenicity, values_from = n) %>% 
  ungroup() %>% 
  mutate(across(-encode, ~replace_na(.x, 0))) %>% 
  mutate(percent_pathogenic = `Pathogenic/Likely_pathogenic`/rowSums(across(!encode))*100)

RBP_Clinvar = synonymous %>% separate_longer_delim(encode, "&") %>% 
  mutate(encode=str_split(encode, "_", simplify=TRUE)[,1]) %>% 
  dplyr::select(pathogenicity) %>% 
  group_by(pathogenicity) %>% 
  summarise(n=n()) 

ggplot(variants_final%>% separate_longer_delim(encode, "&") %>% filter(!is.na(encode)) %>% 
         mutate(encode=str_split(encode, "_", simplify=TRUE)[,1])) +
  geom_bar(aes(x=encode, fill=pathogenicity), stat="count", position="dodge") +
  coord_flip()


test = bridges

test = test %>% mutate(ESE_ESS_RBS = if_else(ESS_ESE_overlap & RBS_overlap, TRUE, FALSE))

ggplot(test %>% filter(pathogenicity=="Pathogenic/Likely_pathogenic")) +
  geom_bar(aes(x=RBS_overlap), stat="count", position="dodge")

