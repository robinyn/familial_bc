library(tidyverse)

setwd("~/server_mnt/H730/noncoding/home/eu3337ha-s/familial")

# targetscan_combined = read_tsv("resources/targetscan/targetscan_combined.bed", col_names=F)
# 
# # Remove unwanted columns and calculate length of predicted target sites
# targetscan_combined = targetscan_combined[,1:8] %>% 
#   mutate(width=X3-X2)
# 
# # Filter target sites with 6 > length > 8
# targetscan_combined = targetscan_combined %>% 
#   filter(width<9) %>% 
#   filter(width>5)

conserved_sites=read_tsv("~/Downloads/Conserved_Site_Context_Scores.txt")
nonconserved_sites=read_tsv("~/Downloads/Nonconserved_Site_Context_Scores.txt")

combined_sites = conserved_sites %>% 
  rbind(nonconserved_sites)

combined_sites = combined_sites %>% 
  filter()