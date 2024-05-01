## Title: pca_population_strat.R
## Author: Euisuk Robin Han
## Last modified: 02 Feb 24
## Description: Script for the PCA of the BRIDGES genotype data to control for population stratification
## Dependencies: Tested and stable with: 
##                R version 4.2.3
##                tidyverse (version 2.0.0)
##                ggplot2 (version 3.4.4)
##                flashpcaR (version 2.1)

library(tidyverse)
library(ggplot2)
library(flashpcaR)

setwd("~/Desktop/RESULTS/BRIDGES/")

# Read sample phenotype files 
control_phenotypes = read_tsv("bridges_annotation/controls_phenotypes.txt")
case_phenotypes = read_tsv("bridges_annotation/cases_phenotypes.txt")

# eigenvec = read_table("test.eigenvec", col_names=FALSE)
# eigenval = scan("test.eigenval")
# 
# eigenvec = eigenvec[,-1]
# names(eigenvec)[1] = "sampleID"
# names(eigenvec)[2:ncol(eigenvec)] <- paste0("PC", 1:(ncol(eigenvec)-1))
# 
# pve = data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# 
# p = ggplot(pve, aes(PC, pve)) +
#   geom_bar(stat = "identity") +
#   ylab("Percentage variance explained") +
#   theme_light()
# 
# print(p)
# 
# p = ggplot(eigenvec, aes(PC1, PC3)) + geom_point(size = 1)
# 
# print(p)

# Read PLINK genotype data (PED/MAP/FAM files)
plink_files = "/Users/robinhan/Desktop/RESULTS/BRIDGES/association/input/plink"

# Run PCA with flashpcaR
pca_res = flashpca(plink_files, ndim=10)

plot_data = pca_res$projection %>% as.data.frame()

plot_data=plot_data %>%
  rownames_to_column(var="sampleID") %>%
  mutate(sampleID=str_remove(sampleID, ":.*"))

phenotypes = case_phenotypes %>%
  dplyr::select(c(BRIDGES_ID, ethnicityClass, ethnicitySubClass)) %>%
  rbind(control_phenotypes %>%
          dplyr::select(c(BRIDGES_ID, ethnicityClass, ethnicitySubClass)))

final_data = plot_data %>%
  left_join(phenotypes, by=join_by(sampleID==BRIDGES_ID))

#Inter01_Lib26_TSP0256_1
#UCAM_Brdg01_Lib47_TSP0107_3

ggplot(final_data %>% filter(ethnicityClass!=888)) +
  geom_point(aes(x=V2, y=V1, color=as.factor(ethnicityClass)))

plot_data=plot_data %>% 
  mutate(IID=sampleID) %>% 
  mutate(FID=sampleID) %>% 
  dplyr::select(c(FID, IID, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10))

write_tsv(plot_data, "~/Desktop/RESULTS/BRIDGES/association/input/pca.tsv", col_names = TRUE)

