library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_GWA.RData")
load("Processed_Data/FileS4_GWA_mapping.RData")

df_ratio_wide <- read_tsv("Processed_Data/ascr_ratio_YA_OutlierRemoval.tsv", col_names = T)

ascr5_res <- process_mapping %>%
  distinct(strain, value) %>%
  na.omit()

ascr5_pro <- df_GWA_fraction_YA %>%
  dplyr::select(strain, pro=ascr.5)

ascr5_comp <- ascr5_pro %>%
  dplyr::left_join(., ascr5_res, by='strain') %>%
  na.omit()

cor(ascr5_comp$pro, ascr5_comp$value)

ascr5_comp %>%
  ggplot(.) +
  geom_point() +
  aes(x=pro, y=value) +
  theme_bw() +
  labs(x="Relative quantity", y="ascr#5 response")

