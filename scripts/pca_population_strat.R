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

setwd("~/Desktop/thesis/association_analysis/")

# Read sample phenotype files 
control_phenotypes = read_tsv("bridges_annotation/controls_phenotypes.txt")
case_phenotypes = read_tsv("bridges_annotation/cases_phenotypes.txt")

# Read PLINK genotype data (PED/MAP/FAM files)
plink_files = "/Users/robinhan/Desktop/thesis/association_analysis/test"

# Run PCA with flashpcaR
pca_res = flashpca(plink_files, ndim=10, stand = "binom2")

plot_data = pca_res$projection %>% as.data.frame()

plot_data=plot_data %>% 
  rownames_to_column(var="sampleID") %>% 
  mutate(sampleID=str_remove(sampleID, ":.*"))

phenotypes = case_phenotypes %>% 
  dplyr::select(c(BRIDGES_ID, ethnicityClass)) %>% 
  rbind(control_phenotypes %>% 
          dplyr::select(c(BRIDGES_ID, ethnicityClass)))

final_data = plot_data %>% 
  left_join(phenotypes, by=join_by(sampleID==BRIDGES_ID))

#Inter01_Lib26_TSP0256_1
#UCAM_Brdg01_Lib47_TSP0107_3

ggplot(final_data) + 
  geom_point(aes(x=PC1, y=PC2, color=as.factor(ethnicityClass))) + xlim(0, 0.001) + ylim(-0.0001, 0.0005)  
  
