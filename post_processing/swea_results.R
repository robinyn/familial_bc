library(tidyverse)
library(ggplot2)

setwd("~/Desktop/thesis")

per_variant=read_tsv("swea_results/per_variant_summary.tsv")
per_transcript=read_tsv("swea_results/reformatted_transcript_summary.tsv")
