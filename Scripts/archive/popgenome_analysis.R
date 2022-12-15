library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(valr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_qtl_peaks.RData")
load("Processed_Data/Ce_Genome-wide_Neutrality_stats.Rda")

df_popgenome <- neutrality_df %>%
  dplyr::rename(chrom=CHROM, start=startWindow, end=endWindow)

QTL_popgenome <- valr::bed_intersect(df_qtl_bed, df_popgenome)

df_popgenome %>%
  dplyr::filter(statistic == "Tajima.D") %>%
  ggplot(.) +
  geom_density(alpha = 0.5, fill='red') +
  aes(x=value) +
  geom_density(data=dplyr::filter(QTL_popgenome, statistic.y == "Tajima.D"), aes(x=value.y), fill='blue', alpha=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(-3,5)

QTL_popgenome_summary <- QTL_popgenome %>%
  na.omit() %>%
  dplyr::group_by(chrom, start.x, end.x, statistic.y) %>%
  dplyr::summarise(value = mean(value.y))
