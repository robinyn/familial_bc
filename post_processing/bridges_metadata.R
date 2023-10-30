library(tidyverse)
library(ggplot2)

setwd("~/Desktop/thesis_resources")

summary_file=read_tsv("bridges_annotation/concept_699_kvist_pheno_v13.txt")

num_unknown_age = summary_file %>% filter(AgeDiagIndex==888) %>% nrow()

age_diagnosis = summary_file %>% 
  select(BRIDGES_ID, AgeDiagIndex, fhnumber, fhscore) %>% 
  filter(AgeDiagIndex!=888)

age_diagnosis$fhnumber[age_diagnosis$fhnumber==888]=NA
age_diagnosis$fhnumber[age_diagnosis$fhnumber==777]=NA

grp1 = age_diagnosis %>% filter(fhnumber==0)
grp2 = age_diagnosis %>% filter(fhnumber>0)

p = ggplot(age_diagnosis) +
  geom_freqpoly(data=grp1, stat="count", aes(x=AgeDiagIndex), size=1) +
  geom_freqpoly(data=grp2, stat="count", aes(x=AgeDiagIndex), size=1) +
  geom_vline(aes(xintercept=mean(AgeDiagIndex),col="red"), linetype=2) +
  geom_vline(aes(xintercept=median(AgeDiagIndex), col="blue"), linetype=2) +
  scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  labs(title="Age at diagnosis of index tumour", x="Age", y="Count") +
  theme(legend.position = "none")

print(p)
