## Title: result_plots.R
## Author: Euisuk Robin Han
## Description: A script to generate the result plots in the thesis report

library(tidyverse)
library(ggplot2)
library(FSA)
library(tidytext)
library(ggh4x)
library(ggridges)

############################ Read SWEA data ############################
swea = read_tsv("~/Desktop/RESULTS/SWEA/swea_all_variants.tsv")
swea_rscu = read_tsv("~/Desktop/RESULTS/SWEA/swea_rscu_quantiles.tsv")
swea_samples_pathogenic = read_tsv("~/Desktop/RESULTS/SWEA/swea_samples_pathogenic.tsv")
swea_synonymous = read_tsv("~/Desktop/RESULTS/SWEA/swea_synonymous_list.tsv")

########################## Read BRIDGES data ###########################
bridges = read_tsv("~/Desktop/RESULTS/BRIDGES/bridges_all_variants.tsv")
bridges_rscu = read_tsv("~/Desktop/RESULTS/BRIDGES/bridges_rscu_quantiles.tsv")
bridges_samples_pathogenic = read_tsv("~/Desktop/RESULTS/BRIDGES/bridges_samples_pathogenic.tsv")
bridges_synonymous = read_tsv("~/Desktop/RESULTS/BRIDGES/bridges_synonymous_list.tsv")

case_phenotypes = read_tsv("~/Desktop/RESULTS/BRIDGES/bridges_annotation/cases_phenotypes.txt") %>% 
  filter(study!="HEBCS") # HEBCS withdrew from BCAC
control_phenotypes = read_tsv("~/Desktop/RESULTS/BRIDGES/bridges_annotation/controls_phenotypes.txt") %>% 
  filter(study!="HEBCS") # HEBCS withdrew from BCAC

########################### Set plot theme #############################
plot_theme = theme(panel.background = element_blank(),
                   axis.line = element_line(color = "black", linewidth = 1),
                   axis.ticks = element_line(linewidth = 1),
                   axis.ticks.length = unit(0.2, "cm"),
                   axis.text.x = element_text(size = 18),
                   axis.text.y = element_text(size = 18),
                   axis.title = element_text(size = 20))

#################### Variant types/consequences ########################
swea_types = swea %>% 
  separate_longer_delim(cols=type, delim="&") %>% 
  group_by(type) %>% 
  summarise(n=n()) %>% 
  arrange(-n) %>% 
  mutate(percent = (n/sum(n))*100)

bridges_types = bridges %>% 
  separate_longer_delim(cols=type, delim="&") %>% 
  group_by(type) %>% 
  summarise(n=n()) %>% 
  arrange(-n) %>% 
  mutate(percent = (n/sum(n))*100)

swea_types = swea_types %>% 
  mutate(type = str_to_sentence(type)) %>% 
  mutate(type = str_replace(type, "utr", "UTR")) %>% 
  mutate(type = str_replace(type, "polypyrimidine", "PT")) %>% 
  add_row(type="Other", 
          n=sum(swea_types$n[swea_types$percent<1]), 
          percent=sum(swea_types$percent[swea_types$percent<1])) %>% 
  filter(percent>1) %>% 
  mutate(data = "SWEA")

bridges_types = bridges_types %>% 
  mutate(type = str_to_sentence(type)) %>% 
  mutate(type = str_replace(type, "utr", "UTR")) %>% 
  mutate(type = str_replace(type, "polypyrimidine", "PT")) %>% 
  add_row(type="Other", 
          n=sum(bridges_types$n[bridges_types$percent<1]), 
          percent=sum(bridges_types$percent[bridges_types$percent<1])) %>% 
  filter(percent>1) %>% 
  mutate(data = "BRIDGES")

df = swea_types %>% rbind(bridges_types) %>% 
  mutate(data = factor(data, c("SWEA","BRIDGES")))

ggplot(df) +
  geom_bar(aes(x=percent, y=reorder(type, percent), fill=type), stat="identity") +
  theme(legend.position = "none") +
  labs(x="Percentage", y="Variant type") +
  facet_wrap("data", nrow=1, scales = "free_y") +
  plot_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=16, hjust=0),
        panel.spacing = unit(1, "lines")) +
  scale_fill_viridis_d()


######################## SAMPLES WITH PATHOGENIC #######################

# SWEA
swea_samples_num_pathogenic = swea_samples_pathogenic %>% 
  group_by(n) %>% 
  summarise(num_path=n())

# BRIDGES
bridges_samples_num_pathogenic = bridges_samples_pathogenic %>% 
  filter(sampleID %in% case_phenotypes$BRIDGES_ID) %>% 
  group_by(n) %>% 
  summarise(num_path=n())


######################### CLINVAR ALL VARIANTS #########################

swea_pathogenicity = swea %>% 
  group_by(pathogenicity) %>% 
  summarise(n=n())%>% 
  filter(pathogenicity!="Unreported") %>% 
  mutate(percent = n/sum(n)*100) %>% 
  mutate(data = "SWEA")

bridges_pathogenicity = bridges %>% 
  group_by(pathogenicity) %>% 
  summarise(n=n())%>% 
  filter(pathogenicity!="Unreported") %>% 
  mutate(percent = n/sum(n)*100) %>% 
  mutate(data = "BRIDGES")

pathogenicity_table = swea_pathogenicity %>% 
  rbind(bridges_pathogenicity) %>% 
  mutate(data = factor(data, levels=c("SWEA", "BRIDGES"))) %>% 
  mutate(pathogenicity = str_replace(pathogenicity, "_", " ")) %>% 
  mutate(type = "All variants")

swea_synonymous_pathogenicity = swea %>%
  filter(str_detect(type, "synonymous")) %>% 
  group_by(pathogenicity) %>% 
  summarise(n=n())%>% 
  mutate(percent = n/sum(n)*100) %>% 
  mutate(data = "SWEA")

bridges_synonymous_pathogenicity = bridges %>% 
  filter(str_detect(type, "synonymous")) %>% 
  group_by(pathogenicity) %>% 
  summarise(n=n())%>% 
  mutate(percent = n/sum(n)*100) %>% 
  mutate(data = "BRIDGES")

synonymous_pathogenicity_table = swea_synonymous_pathogenicity %>% 
  rbind(bridges_synonymous_pathogenicity) %>% 
  mutate(data = factor(data, levels=c("SWEA", "BRIDGES"))) %>% 
  mutate(pathogenicity = str_replace(pathogenicity, "_", " ")) %>% 
  mutate(type = "Synonymous variants")

pathogenicity_table = pathogenicity_table %>% 
  rbind(synonymous_pathogenicity_table)

rm(list=c("swea_pathogenicity", "bridges_pathogenicity", 
          "swea_synonymous_pathogenicity", "bridges_synonymous_pathogenicity", 
          "synonymous_pathogenicity_table"))

ggplot() + 
  geom_bar(data=pathogenicity_table, 
           aes(y=reorder_within(pathogenicity, percent, list(type, data)), x=percent, fill=as.factor(pathogenicity)), stat="identity") + 
  scale_fill_viridis_d() +
  scale_y_reordered() +
  theme(legend.position = "none") +
  labs(y="Clinical significance", x="Percentage") +
  facet_grid2(type~data, scales = "free", independent = "y") +
  plot_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=18),
        strip.text.x = element_text(hjust=0),
        panel.spacing = unit(1, "line"))


######################### phyloP  #########################

df = swea %>% dplyr::select(variant, type, pathogenicity, phyloP) %>% 
  separate_longer_delim(type, delim="&") %>% 
  mutate(data="SWEA") %>% 
  rbind(bridges %>% 
          dplyr::select(variant, type, pathogenicity, phyloP) %>% 
          separate_longer_delim(type, delim="&") %>% 
          mutate(data = "BRIDGES"))

df_summary = df %>% 
  group_by(data, type) %>% 
  summarise(n=n()) %>% 
  filter(n<10)

df = df %>% 
  mutate(type = if_else(type %in% df_summary$type, "Other", type)) %>% 
  mutate(data = factor(data, levels=c("SWEA", "BRIDGES")))
  
ggplot(df %>% 
         filter(!is.na(phyloP)) %>% 
         mutate(pathogenicity=str_replace(pathogenicity, "_", " ")),
       aes(y=reorder(pathogenicity, phyloP, mean), x=phyloP, fill=after_stat(x))) +
  #geom_boxplot(outlier.shape = NA, na.rm = TRUE, position = position_nudge(x=-0.2)) +
  #geom_jitter(alpha=0.5, width=0.2) +
  geom_density_ridges_gradient(quantile_lines=TRUE, quantile_fun=mean) +
  #stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom="pointrange", color="black") +
  scale_fill_steps2(breaks=c(log10(0.05), -log10(0.05)), low="blue", mid="darkgrey", high="red") +
  plot_theme +
  labs(x="phyloP", y="Clinical significance") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text=element_text(size=18, hjust=0)) +
  facet_wrap("data")

ggplot(df %>% 
         filter(!is.na(phyloP)) %>% 
         mutate(pathogenicity=str_replace(type, "_", " ")) %>% 
         separate_longer_delim(type, delim="&") %>% 
         mutate(type=str_to_sentence(type)) %>% 
         mutate(type = str_replace(type, "polypyrimidine", "PT")) %>% 
         mutate(type = str_replace(type, "utr", "UTR")),
       aes(y=reorder(type, phyloP, mean), x=phyloP, fill=after_stat(x))) +
  #geom_boxplot(outlier.shape = NA, na.rm = TRUE, position = position_nudge(x=-0.2)) +
  #geom_jitter(alpha=0.5, width=0.2) +
  geom_density_ridges_gradient(quantile_lines=TRUE, quantile_fun=mean) +
  #stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom="pointrange", color="black") +
  scale_fill_steps2(breaks=c(log10(0.05), -log10(0.05)), low="blue", mid="darkgrey", high="red") +
  plot_theme +
  labs(x="phyloP", y="Clinical significance") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text=element_text(size=18, hjust=0)) +
  facet_wrap("data")

######################### RSCU - ClinVar #########################

df = swea %>% dplyr::select(variant, type, pathogenicity, rscu) %>% 
  separate_longer_delim(type, delim="&") %>% 
  mutate(data="SWEA") %>% 
  rbind(bridges %>% 
          dplyr::select(variant, type, pathogenicity, rscu) %>% 
          separate_longer_delim(type, delim="&") %>% 
          mutate(data = "BRIDGES")) %>% 
  mutate(data = factor(data, levels=c("SWEA", "BRIDGES"))) %>% 
  mutate(pathogenicity=str_replace(pathogenicity, "_", " "))

ggplot(df %>% filter(!is.na(rscu)),
       aes(y=reorder(pathogenicity, rscu, mean), x=rscu, fill=pathogenicity)) +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=mean) +
  plot_theme +
  labs(x="dRSCU", y="Clinical significance") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text=element_text(size=18, hjust=0)) +
  facet_wrap("data", scales="free_y", nrow=2) +
  geom_vline(data=df %>% filter(data == "SWEA"), aes(xintercept=swea_rscu$rscu[swea_rscu$quantiles=="5%"]), linetype=2) +
  geom_vline(data=df %>% filter(data == "SWEA"), aes(xintercept=swea_rscu$rscu[swea_rscu$quantiles=="95%"]), linetype=2) +
  geom_vline(data=df %>% filter(data == "BRIDGES"), aes(xintercept=bridges_rscu$rscu[swea_rscu$quantiles=="5%"]), linetype=2) +
  geom_vline(data=df %>% filter(data == "BRIDGES"), aes(xintercept=bridges_rscu$rscu[swea_rscu$quantiles=="95%"]), linetype=2) +
  scale_fill_viridis_d()

ggplot(df %>% filter(!is.na(rscu)), aes(x=rscu)) +
  geom_density() +
  facet_wrap("data") +
  plot_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text=element_text(size=18, hjust=0)) +
  geom_vline(data=df %>% filter(data == "SWEA"), aes(xintercept=swea_rscu$rscu[swea_rscu$quantiles=="5%"]), linetype=2) +
  geom_vline(data=df %>% filter(data == "SWEA"), aes(xintercept=swea_rscu$rscu[swea_rscu$quantiles=="95%"]), linetype=2) +
  geom_vline(data=df %>% filter(data == "BRIDGES"), aes(xintercept=bridges_rscu$rscu[swea_rscu$quantiles=="5%"]), linetype=2) +
  geom_vline(data=df %>% filter(data == "BRIDGES"), aes(xintercept=bridges_rscu$rscu[swea_rscu$quantiles=="95%"]), linetype=2) 

