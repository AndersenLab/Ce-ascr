library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggcorrplot)
library(picante)
library(qgraph)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_ascr_data_norm_h2.RData")

test_result_YA <- cor.table(df_YA_pheno_wide, cor.method="spearman")[[1]]
test_result_YA_h2_filtered <- cor.table(df_YA_pheno_wide_h2_filtered, cor.method="spearman")[[1]]

## corr heatmap

ggcorrplot(test_result_YA, hc.order = TRUE, tl.srt = 90, type = "lower", legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=9, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=9, color='black'),
        legend.position = c(0.5,0.8),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.2,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr.pdf", width = 7.5, height = 7.5, units = "in")

ggcorrplot(test_result_YA_h2_filtered, hc.order = TRUE, tl.srt = 90, type = "lower", legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=9, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=9, color='black'),
        legend.position = c(0.5,0.8),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.2,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr_h2_filtered.pdf", width = 7.5, height = 7.5, units = "in")

## scatter plot for important inverse correlation

## network analysis
 #cor
Graph_cor_YA_h2_filtered <- qgraph(test_result_YA_h2_filtered, graph = "cor", layout = "spring", 
                        labels = colnames(test_result_YA_h2_filtered), minimum = 0.1)  #  Unregularized partial correlation network
qgraph(Graph_cor_YA_h2_filtered, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_cor_min01_h2_filtered")
 #pcor
Graph_pcor_YA_h2_filtered <- qgraph(test_result_YA_h2_filtered, graph = "pcor", layout = "spring", 
                        labels = colnames(test_result_YA_h2_filtered), minimum = 0.1)  #  Unregularized partial correlation network
qgraph(Graph_pcor_YA_h2_filtered, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_pcor_min01_h2_filtered")


### analysis for each pathway ###

 #1. Biosynthetic pathway for simple ascarosides (a)

df_simple1 <- test_result_YA %>%
  as.data.frame %>%
  dplyr::select("bhas#16", "ascr#12", "bhas#12", "ascr#11","bhas#11") %>%
  t(.) %>%
  as.data.frame %>%
  dplyr::select("bhas#16", "ascr#12", "bhas#12", "ascr#11","bhas#11") %>% 
  t(.)

## corr heatmap

ggcorrplot(df_simple1, hc.order = FALSE, tl.srt = 90, legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.2,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr_simple1.pdf", width = 7.5, height = 7.5, units = "in")

## network analysis

Graph_cor_YA_simple1 <- qgraph(df_simple1, graph = "cor", layout = "spring", 
                        labels = colnames(df_simple1), minimum = 0.1)  #  Unregularized partial correlation network
qgraph(Graph_cor_YA_simple1, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_cor_min01_simple1")

Graph_pcor_YA_simple1 <- qgraph(df_simple1, graph = "pcor", layout = "spring", 
                               labels = colnames(df_simple1), minimum = 0.1)  #  Unregularized partial correlation network
qgraph(Graph_pcor_YA_simple1, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_pcor_min01_simple1")

## scatter plot for important edge

plot_simple1_scatter1 <- df_GWA_fraction_YA %>%
  as.data.frame %>%
  dplyr::select(ascr.11, bhas.11, bhas.12) %>%
  ggplot(.) + 
  geom_point(alpha=0.8) +
  aes(x=bhas.11, y=bhas.12) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="bhas#11", y="bhas#12")

plot_simple1_scatter1

plot_simple1_scatter2 <- df_GWA_fraction_YA %>%
  as.data.frame %>%
  dplyr::select(ascr.11, bhas.11, bhas.12) %>%
  ggplot(.) + 
  geom_point(alpha=0.8) +
  aes(x=ascr.11, y=bhas.12) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="ascr#11", y="bhas#12")

plot_simple1_scatter2

ggsave(plot_simple1_scatter1, file="Plots/plot_simple1_scatter1.pdf", width = 5, height = 4, units = "in")
ggsave(plot_simple1_scatter2, file="Plots/plot_simple1_scatter2.pdf", width = 5, height = 4, units = "in")


#2. Biosynthetic pathway for simple ascarosides (b)

df_simple2 <- test_result_YA %>%
  as.data.frame %>%
  dplyr::select("bhos#22", "oscr#18", "bhos#18", "oscr#10", "oscr#3", "bhos#10", "oscr#1", "oscr#7", "ascr#5") %>%
  t(.) %>%
  as.data.frame %>%
  dplyr::select("bhos#22", "oscr#18", "bhos#18", "oscr#10", "oscr#3", "bhos#10", "oscr#1", "oscr#7", "ascr#5") %>%
  t(.)

## corr heatmap

ggcorrplot(df_simple2, hc.order = FALSE, tl.srt = 90, legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.2,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr_simple2.pdf", width = 7.5, height = 7.5, units = "in")

## network analysis

Graph_cor_YA_simple2 <- qgraph(df_simple2, graph = "cor", layout = "spring", 
                                labels = colnames(df_simple2), minimum = 0.1)  #  Unregularized partial correlation network
qgraph(Graph_cor_YA_simple2, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_cor_min01_simple2")

Graph_pcor_YA_simple2 <- qgraph(df_simple2, graph = "pcor", layout = "spring", 
                               labels = colnames(df_simple2), minimum = 0.1)  #  Unregularized partial correlation network
qgraph(Graph_pcor_YA_simple2, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_pcor_min01_simple2")

## scatter plot for important edge

plot_simple2_scatter1 <- df_GWA_fraction_YA %>%
  as.data.frame %>%
  dplyr::select(oscr.10, oscr.18, ascr.5, oscr.7) %>%
  ggplot(.) + 
  geom_point(alpha=0.8) +
  aes(x=oscr.10, y=oscr.18) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="oscr#10", y="oscr#18")

plot_simple2_scatter1

plot_simple2_scatter2 <- df_GWA_fraction_YA %>%
  as.data.frame %>%
  dplyr::select(oscr.10, oscr.18, ascr.5, oscr.7) %>%
  ggplot(.) + 
  geom_point(alpha=0.8) +
  aes(x=ascr.5, y=oscr.7) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="ascr#5", y="oscr#7")

plot_simple2_scatter2

ggsave(plot_simple2_scatter1, file="Plots/plot_simple2_scatter1.pdf", width = 5, height = 4, units = "in")
ggsave(plot_simple2_scatter2, file="Plots/plot_simple2_scatter2.pdf", width = 5, height = 4, units = "in")


#3. Biosynthetic pathway for modified ascarosides (b)

df_modified <- test_result_YA %>%
  as.data.frame %>%
  dplyr::select("bhas#26", "bhas#22", "ascr#18", "bhas#18", "ascr#10", "ascr#3", "bhas#10", 
                "ascr#1", "ascr#7", "ascr#9", "glas#3", "anglas#7", "icas#3", "icas#9", "icas#10", 
                "osas#2", "osas#9", "osas#10",  
                "iglas#2", "ascr#81") %>%
  t(.) %>%
  as.data.frame %>%
  dplyr::select("bhas#26", "bhas#22", "ascr#18", "bhas#18", "ascr#10", "ascr#3", "bhas#10", 
                "ascr#1", "ascr#7", "ascr#9", "glas#3", "anglas#7", "icas#3", "icas#9", "icas#10", 
                "osas#2", "osas#9", "osas#10",  
                "iglas#2", "ascr#81") %>%
  t(.)

## corr heatmap

ggcorrplot(df_modified, hc.order = TRUE, tl.srt = 90, legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.2,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr_modified.pdf", width = 7.5, height = 7.5, units = "in")

## network analysis

Graph_cor_YA_modified <- qgraph(df_modified, graph = "cor", layout = "spring", 
                                labels = colnames(df_modified), minimum = 0.2)  #  Unregularized partial correlation network
qgraph(Graph_cor_YA_modified, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_cor_min02_modified")

Graph_pcor_YA_modified <- qgraph(df_modified, graph = "pcor", layout = "spring", 
                                labels = colnames(df_modified), minimum = 0.2)  #  Unregularized partial correlation network
qgraph(Graph_pcor_YA_modified, filetype = "pdf", height = 7.5, width = 7.5, filename = "Plots/YA_cor_network_pcor_min02_modified")

## scatter plot for important edge

plot_modified_scatter1 <- df_GWA_fraction_YA %>%
  as.data.frame %>%
  dplyr::select(ascr.3, icas.10) %>%
  ggplot(.) + 
  geom_point(alpha=0.8) +
  aes(x=ascr.3, y=icas.10) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="ascr#3", y="icas#10")

plot_modified_scatter1

  plot_modified_scatter2 <- df_GWA_fraction_YA %>%
  as.data.frame %>%
  dplyr::select(oscr.10, oscr.18, ascr.5, oscr.7) %>%
  ggplot(.) + 
  geom_point(alpha=0.8) +
  aes(x=ascr.5, y=oscr.7) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="ascr#5", y="oscr#7")

plot_modified_scatter2

ggsave(plot_modified_scatter1, file="Plots/plot_modified_scatter1.pdf", width = 5, height = 4, units = "in")
ggsave(plot_modified_scatter2, file="Plots/plot_modified_scatter2.pdf", width = 5, height = 4, units = "in")


### comparisons among simple 1, 2 and modified

df_comp <- test_result_YA %>%
  as.data.frame %>%
  dplyr::select("bhas#16", "ascr#12", "bhas#12", "ascr#11","bhas#11",
                "bhos#22", "oscr#18", "bhos#18", "oscr#10", "oscr#3", "bhos#10", "oscr#1", "oscr#7", "ascr#5",
    "bhas#26", "bhas#22", "ascr#18", "bhas#18", "ascr#10", "ascr#3", "bhas#10", 
                "ascr#1", "ascr#7", "ascr#9", "glas#3", "anglas#7", "icas#3", "icas#9", "icas#10", 
                "osas#2", "osas#9", "osas#10",  
                "iglas#2", "ascr#81") %>%
  t(.) %>%
  as.data.frame %>%
  dplyr::select("bhas#16", "ascr#12", "bhas#12", "ascr#11","bhas#11",
                "bhos#22", "oscr#18", "bhos#18", "oscr#10", "oscr#3", "bhos#10", "oscr#1", "oscr#7", "ascr#5",
                "bhas#26", "bhas#22", "ascr#18", "bhas#18", "ascr#10", "ascr#3", "bhas#10", 
                "ascr#1", "ascr#7", "ascr#9", "glas#3", "anglas#7", "icas#3", "icas#9", "icas#10", 
                "osas#2", "osas#9", "osas#10",  
                "iglas#2", "ascr#81") %>%
  t(.)

## corr heatmap

# ordered by pahtway
ggcorrplot(df_comp, hc.order = F, tl.srt = 90, legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.2,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr_comp.pdf", width = 10, height = 10, units = "in")

# ordered by similarity
ggcorrplot(df_comp, hc.order = T, tl.srt = 90, legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.2,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr_comp_clustered.pdf", width = 10, height = 10, units = "in")

