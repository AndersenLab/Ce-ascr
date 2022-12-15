library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggrepel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

ascr_YA <- read.csv(file="Raw/20190807_ascr_raw_YA.csv") %>%
  dplyr::select(-ascr_total, -iglu.1, -iglu.2, -oscr.3, -oscr.7) %>%
  tidyr::separate(strain, sep = " ", into = c("strain", "batch"), remove = F)  %>%
  dplyr::filter(!strain %in% c("JU238", "JU363", "JU1516", "QG537", "QG538", "CX11400", "QG1", "LSJ1")) %>%
  dplyr::mutate(strain = ifelse(strain == "JU491", "JU1491", 
                                ifelse(strain=="JU1580", "JU1793", strain))) %>% ## fix strain name
  dplyr::mutate(ascr_sum = dplyr::select(., anglas.3:uglas.11) %>% rowSums()) %>%
  dplyr::ungroup()

ascr_YA_summary <- ascr_YA %>%
  dplyr::mutate(norm_ascr_sum=ascr_sum/bpTIC) %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(norm_ascr_sum= mean(norm_ascr_sum), ascr_sum=mean(ascr_sum)) %>%
  dplyr::ungroup()

ascr_YA_summary %>%
  ggplot(.) +
  geom_col() +
  aes(x=reorder(strain, -norm_ascr_sum), y=norm_ascr_sum) +
  theme_bw() +
  theme(axis.text.x = element_text(size= 8, angle = 90, vjust = 0.5, hjust=1), axis.title.x=element_blank(),
        panel.grid = element_blank()) +
  labs(y="Normalized ascaroside abundance")

ggsave(file = "Plots/plot_bpTIC_norm.pdf", width = 7.5, height = 4)

df_GWA_norm_abundance <- ascr_YA_summary %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::select(strain, norm_ascr_sum)

write_tsv(df_GWA_norm_abundance, path="Processed_Data/ascr_abundance_YA_bpTICnorm.tsv", col_names = T)

daf22_gt <- read.table(file="Processed_Data/daf_22_Glu266Ala.tsv") %>%
  dplyr::mutate(tested = ifelse(V4 %in% unique(ascr_YA_summary$strain), "Yes", "No")) %>%
  dplyr::mutate(V4=ifelse(V4=="ECA259", "PB306", V4)) %>%
  dplyr::mutate(daf22_gt=ifelse(V5=="0/0", "REF", 
                                  ifelse(V5=="./.", "Del", "Q266A")))

ascr_YA_summary %>%
  dplyr::left_join(dplyr::select(daf22_gt, strain=V4, daf22_gt)) %>%
  ggplot(.) +
  geom_col() +
  aes(x=reorder(strain, -norm_ascr_sum), y=norm_ascr_sum, fill=factor(daf22_gt, levels=c("REF","Q266A", "Del"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size= 8, angle = 90, vjust = 0.5, hjust=1), axis.title.x=element_blank(),
        panel.grid = element_blank(),
        legend.title=element_blank(),
        legend.position = c(0.8,0.8)) +
  labs(y="Normalized ascaroside abundance") +
  scale_fill_manual(values = c("grey70","blue","red"))

ggsave(file = "Plots/plot_bpTIC_norm_daf22_Q266A.png", width = 7.5, height = 4)

ascr_YA_summary %>%
  dplyr::left_join(dplyr::select(daf22_gt, strain=V4, daf22_gt)) %>%
  ggplot(.) +
  geom_jitter(width=0.2) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.8) +
  aes(x=factor(daf22_gt, levels=c("REF","Q266A", "Del")), y=norm_ascr_sum, fill=factor(daf22_gt, levels=c("REF","Q266A", "Del"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size= 8, angle = 90, vjust = 0.5, hjust=1), axis.title.x=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none') +
  labs(y="Normalized ascaroside abundance") +
  scale_fill_manual(values = c("grey70","blue","red"))

ggsave(file = "Plots/plot_bpTIC_norm_daf22_Q266A_boxplot.png", width = 4, height = 4)

ascr_YA_summary %>%
  ggplot(.) +
  geom_col() +
  aes(x=reorder(strain, -norm_ascr_sum), y=log10(norm_ascr_sum)) +
  theme_bw() +
  theme(axis.text.x = element_text(size= 7.5, angle = 90, vjust = 0.5, hjust=1, color='black'), axis.title.x=element_blank(),
        panel.grid = element_blank()) +
  labs(y="log10 (Normalized ascaroside abundance)")
  
ggsave(file = "Plots/plot_bpTIC_norm_log.pdf", width = 7.5, height = 4)

### Calculating ascaroside fractions

ascr_YA_frac <- ascr_YA %>%
  tidyr::gather(-strain, -batch, -bpTIC, -ascr_sum, key=feature, value=abundance) %>%
  dplyr::group_by(strain, feature) %>%
  dplyr::mutate(abundance = mean(abundance), ascr_sum = mean(ascr_sum)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, feature, abundance, .keep_all = T) %>%
  dplyr::select(-batch) %>%
  dplyr::mutate(ascr_fraction = abundance/ascr_sum)

ascr_YA_frac$feature <- factor(ascr_YA_frac$feature, levels = c("ascr.1", "ascr.3", "ascr.5", "ascr.7", "ascr.8", "ascr.9", "ascr.10", "ascr.11", "ascr.12", "ascr.15", "ascr.18", "ascr.22",  "ascr.81", "anglas.3", "anglas.7", "bhas.10", "bhas.11", "bhas.12", "bhas.16", "bhas.18", "bhas.22", "bhas.26", "bhos.10", "bhos.18", "bhos.22", "glas.1", "glas.3", "icas.3", "icas.9", "icas.10", "iglas.2", "osas.2", "osas.9", "osas.10", "oscr.1", "oscr.9", "oscr.10", "oscr.11", "oscr.12", "oscr.14", "oscr.18", "oscr.19", "oscr.21", "uglas.11"),
                                      labels = c("ascr#1", "ascr#3", "ascr#5", "ascr#7", "ascr#8", "ascr#9", "ascr#10", "ascr#11", "ascr#12", "ascr#15", "ascr#18", "ascr#22",  "ascr#81", "anglas#3", "anglas#7", "bhas#10", "bhas#11", "bhas#12", "bhas#16", "bhas#18", "bhas#22", "bhas#26", "bhos#10", "bhos#18", "bhos#22", "glas#1", "glas#3", "icas#3", "icas#9", "icas#10", "iglas#2", "osas#2", "osas#9", "osas#10", "oscr#1", "oscr#9", "oscr#10", "oscr#11", "oscr#12", "oscr#14", "oscr#18", "oscr#19", "oscr#21", "uglas#11"))

for (f in unique(ascr_YA_frac$feature)) {
  
  ascr_YA_frac %>%
    dplyr::filter(feature == f) %>%
    dplyr::filter(strain != "JU1400") %>%
    ggplot(.) +
    geom_col() +
    aes(x=reorder(strain, -ascr_fraction), y=ascr_fraction*100) +
    theme_bw() +
    theme(axis.text.x = element_text(size= 9, angle = 90, vjust = 0.5, hjust=1), 
          axis.title.y = element_text(size= 10), 
          axis.title.x=element_blank(),
          panel.grid = element_blank()) +
    labs(y="Relative abundance")
  
  ggsave(glue::glue("Plots/ascr_fraction/{f}_fraction.png"), width = 10, height = 4)
  
}

### df wide spread ###

df_YA_pheno <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::select(strain, feature, ascr_fraction) %>%
  tidyr::spread(key=strain, value=ascr_fraction)

df_YA_pheno_wide <- t(df_YA_pheno[2:ncol(df_YA_pheno)])
colnames(df_YA_pheno_wide) <- df_YA_pheno[[1]]

# grouping of ascarosides

ascr_short = c("ascr#3", "ascr#1", "ascr#12", "ascr#10", "ascr#9", "ascr#7", "ascr#11")
ascr_long = c("ascr#18", "ascr#22", "ascr#15")
oscr_short = c("ascr#5", "oscr#1", "oscr#12", "oscr#10", "oscr#11", "oscr#9", "oscr#14")
oscr_long = c("oscr#21", "oscr#18", "oscr#19")
bhos = c("bhos#10", "bhos#18", "bhos#22")
four_modified = c("icas#3", "icas#9", "icas#10", "osas#2", "osas#9", "osas#10")
C_modified = c("iglas#2", "glas#1", "glas#3", "anglas#3", "anglas#7", "ascr#81", "ascr#8", "uglas#11")
bhas_short = c("bhas#10", "bhas#11", "bhas#12")
bhas_long = c("bhas#18", "bhas#22", "bhas#16", "bhas#20", "bhas#26")

df_ascr_group <- data.frame(feature = c(as.character(unique(df_YA_pheno$feature)))) %>%
  dplyr::mutate(group1 = ifelse(feature %in% ascr_short, "ascr(s)",
                                ifelse(feature %in% ascr_long, "ascr(l)",
                                       ifelse(feature %in% oscr_short, "oscr(s)",
                                              ifelse(feature %in% oscr_long, "oscr(l)",
                                                     ifelse(feature %in% bhos, "bhos",
                                                            ifelse(feature %in% four_modified, "modified(4')",
                                                                   ifelse(feature %in% C_modified, "modified(C)", 
                                                                          ifelse(feature %in% bhas_short, "bhas(s)",
                                                                                 ifelse(feature %in% bhas_long, "bhas(l)", "etc")))))))))) %>%
  dplyr::group_by(group1) %>%
  dplyr::mutate(n_member_group1=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group2 = ifelse(feature %in% c(ascr_short, ascr_long), "ascr",
                                ifelse(feature %in% c(oscr_short, oscr_long), "oscr",
                                       ifelse(feature %in% c(bhas_long, bhas_short), "bhas", group1))))

df_ascr_group$group1 <- factor(df_ascr_group$group1, levels = c("ascr(s)", "ascr(l)", "bhas(s)", "bhas(l)", "modified(C)", "modified(4')", "oscr(s)", "oscr(l)", "bhos"))

group_pal <- c("#9E0142", "#F46D43", "#FDAE61", "#FEE08B", "#A6761D", "#666666", "#3288BD", "#5E4FA2", "#66C2A5")

### average ascaroside fraction and cv across all wild isolates ###

plot_ascr_pop_log <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::left_join(., df_ascr_group, by = 'feature') %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 0.5) +
  aes(x=reorder(feature, -ascr_fraction, FUN=median), y=log10(ascr_fraction*100), fill= group1) +
  theme_bw() +
  theme(axis.text.x = element_text(size= 9.5, angle = 90, vjust = 0.5, hjust=1, color='black'), 
        axis.text.y = element_text(size=9.5, color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color='black'),
        panel.grid = element_blank(),
        legend.position ='right',
        legend.title = element_blank(),
        legend.text = element_text(size=9, color='black'),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.2,"in")) +
  labs(y="log10(Relative quantity (%))") +
  #scale_y_continuous(trans='log10', breaks = c(10, 0.1, 0.001), labels = comma) +
  scale_fill_manual(values = group_pal) +
  guides(fill= guide_legend(ncol=1))

plot_ascr_pop_log

ggsave(plot_ascr_pop_log, file = "Plots/plot_ascr_pop_log.png", width = 7.5, height = 3)

plot_ascr_pop <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::left_join(., df_ascr_group, by = 'feature') %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 0.5) +
  aes(x=reorder(feature, -ascr_fraction, FUN=median), y=ascr_fraction*100, fill= group1) +
  theme_bw() +
  theme(axis.text.x = element_text(size= 9.5, angle = 90, vjust = 0.5, hjust=1, color='black'), 
        axis.text.y = element_text(size=9.5, color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color='black'),
        panel.grid = element_blank(),
        legend.position ='right',
        legend.title = element_blank(),
        legend.text = element_text(size=9, color='black'),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.2,"in")) +
  labs(y="log10(Relative quantity (%))") +
  #scale_y_continuous(trans='log10', breaks = c(10, 0.1, 0.001), labels = comma) +
  scale_fill_manual(values = group_pal) +
  guides(fill= guide_legend(ncol=1))

plot_ascr_pop

ggsave(plot_ascr_pop, file = "Plots/plot_ascr_pop.png", width = 7.5, height = 3)

df_ascr_YA_frac_summary <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(median = median(ascr_fraction), cv=sd(ascr_fraction)/mean(ascr_fraction)) %>%
  dplyr::ungroup()

### heat map for ascr_fraction ###

plot_heatmap_ascrfrac <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::left_join(., df_ascr_group, by = 'feature') %>%
  ggplot(.) +
  geom_tile() +
  aes(x=feature, y=strain, fill=log2(ascr_fraction)) +
  theme_bw() +
  theme(axis.text.x = element_text(size= 9, angle = 90, vjust = 0.5, hjust=1, color='black'), 
        axis.text.y = element_text(size= 9, color='black'), 
        axis.title=element_blank(), 
        legend.position = 'bottom',
        legend.key.width = unit(0.3,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=8, color='black'),
        panel.spacing = unit(0.07, 'line'),
        strip.text = element_text(size=8, color='black')) +
  scale_fill_gradient(low = "black", high = "turquoise") +
  facet_grid(~group1, scale='free', space='free') +
  labs(fill="log2 (Relative quantity)")

plot_heatmap_ascrfrac

ggsave(plot_heatmap_ascrfrac, file = "Plots/plot_heatmap_ascrfrac.pdf", width = 7.5, height = 10.5)


### PCA and outlier detection ###

strain_YA<-row.names(df_YA_pheno_wide)

df_YA_ascrfrac_PC_trim.fit_scaled <- prcomp(df_YA_pheno_wide, scale. = TRUE)
df_YA_PCA_summary <- data.frame(t(summary(df_YA_ascrfrac_PC_trim.fit_scaled)[[6]]), PC=c(1:44))

plot_YA_PCA_summary <- df_YA_PCA_summary %>%
  ggplot(.) +
  geom_col() +
  aes(x=PC, y=Proportion.of.Variance) +
  theme_bw()

plot_YA_PCA_summary_cum <- df_YA_PCA_summary %>%
  ggplot(.) +
  geom_col() +
  aes(x=PC, y=Cumulative.Proportion) +
  theme_bw()

ggsave(plot_YA_PCA_summary, file = "Plots/plot_YA_PCA_summary.pdf", width = 7.5, height = 3.5)
ggsave(plot_YA_PCA_summary_cum, file = "Plots/plot_YA_PCA_summary_cum.pdf", width = 7.5, height = 3.5)


#df_ascr_group$group1 <- factor(df_ascr_group$group, levels = c("ascr_short", "oscr_short", "bhas_short", "four_modified", "C_modified", "bhos", "ascr_long", "oscr_long","bhas_long", "etc"))

## plot for loadings

df_YA_ascrfrac_PC_loading<-data.frame(feature=rownames(df_YA_ascrfrac_PC_trim.fit_scaled[[2]]), df_YA_ascrfrac_PC_trim.fit_scaled[[2]]) %>%
  tidyr::gather(key="PC", value="loading", -feature)

# group1

plot_YA_PC1_loading1 <- df_YA_ascrfrac_PC_loading %>%
  dplyr::filter(PC == "PC1") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='none') +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC1")

plot_YA_PC1_loading1

plot_YA_PC2_loading1 <- df_YA_ascrfrac_PC_loading %>%
  dplyr::filter(PC == "PC2") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='none') +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC2")

plot_YA_PC3_loading1 <- df_YA_ascrfrac_PC_loading %>%
  dplyr::filter(PC == "PC3") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='bottom',
        legend.title = element_blank()) +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC3") +
  guides(fill= guide_legend(nrow=2))

## group2

plot_YA_PC1_loading2 <- df_YA_ascrfrac_PC_loading %>%
  dplyr::filter(PC == "PC1") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='none') +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC1")

plot_YA_PC1_loading2

plot_YA_PC2_loading2 <- df_YA_ascrfrac_PC_loading %>%
  dplyr::filter(PC == "PC2") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='none') +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC2")

plot_YA_PC3_loading2 <- df_YA_ascrfrac_PC_loading %>%
  dplyr::filter(PC == "PC3") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='bottom',
        legend.title = element_blank()) +
  ggtitle("PC3") +
  scale_fill_manual(values = group_pal) +
  guides(fill= guide_legend(nrow=1))


plot_YA_PC1_loading_group1 <- cowplot::plot_grid(plot_YA_PC1_loading1, plot_YA_PC2_loading1, plot_YA_PC3_loading1, 
                                                 ncol = 1, rel_heights = c(1,1,1.2))
#plot_YA_PC1_loading_group2 <- cowplot::plot_grid(plot_YA_PC1_loading2, plot_YA_PC2_loading2, plot_YA_PC3_loading2, 
#                                                 ncol = 1, rel_heights = c(1,1,1.2))

ggsave(plot_YA_PC1_loading_group1, file = "Plots/plot_YA_PC1_loading_group1.pdf", width = 7.5, height = 9)
#ggsave(plot_YA_PC1_loading_group2, file = "Plots/plot_YA_PC1_loading_group2.pdf", width = 7.5, height = 9)

df_YA_ascrfrac_PC_traits <-df_YA_ascrfrac_PC_trim.fit_scaled[[5]]
#pc_traits <-data.frame(pcDF.fit[[5]])
df_YA_ascrfrac_PC_traits <- data.frame(strain=strain_YA, df_YA_ascrfrac_PC_traits)

plot_PC12_YA <- df_YA_ascrfrac_PC_traits %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, PC1 > mean(df_YA_ascrfrac_PC_traits$PC1)+6*sd(df_YA_ascrfrac_PC_traits$PC1) |
                                          PC1 < mean(df_YA_ascrfrac_PC_traits$PC1)-6*sd(df_YA_ascrfrac_PC_traits$PC1)), aes(label = strain), size = 2.5) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, PC2 > mean(df_YA_ascrfrac_PC_traits$PC2)+6*sd(df_YA_ascrfrac_PC_traits$PC2) |
                                          PC2 < mean(df_YA_ascrfrac_PC_traits$PC2)-6*sd(df_YA_ascrfrac_PC_traits$PC2)), aes(label = strain), size = 2.5) +
  #geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, strain %in% c("ED3052", "JU1242")), aes(label = strain), size = 2.5) +
  aes(x=PC1, y=PC2) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC1 (19.7%)", y="PC2 (15.7%)")

plot_PC12_YA

plot_PC13_YA <- df_YA_ascrfrac_PC_traits %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, PC1 > mean(df_YA_ascrfrac_PC_traits$PC1)+6*sd(df_YA_ascrfrac_PC_traits$PC1) |
                                          PC1 < mean(df_YA_ascrfrac_PC_traits$PC1)-6*sd(df_YA_ascrfrac_PC_traits$PC1)), aes(label = strain), size = 2.5) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, PC3 > mean(df_YA_ascrfrac_PC_traits$PC3)+6*sd(df_YA_ascrfrac_PC_traits$PC3) |
                                          PC3 < mean(df_YA_ascrfrac_PC_traits$PC3)-6*sd(df_YA_ascrfrac_PC_traits$PC3)), aes(label = strain), size = 2.5) +
  #geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, strain %in% c("ED3052", "JU1242")), aes(label = strain), size = 2.5) +
  aes(x=PC1, y=PC3) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC1 (19.7%)", y="PC3 (8.6%)")

plot_PC13_YA

plot_PC23_YA <- df_YA_ascrfrac_PC_traits %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, PC2 > mean(df_YA_ascrfrac_PC_traits$PC2)+6*sd(df_YA_ascrfrac_PC_traits$PC2) |
                                          PC2 < mean(df_YA_ascrfrac_PC_traits$PC2)-6*sd(df_YA_ascrfrac_PC_traits$PC2)), aes(label = strain), size = 2.5) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, PC3 > mean(df_YA_ascrfrac_PC_traits$PC3)+6*sd(df_YA_ascrfrac_PC_traits$PC3) |
                                          PC3 < mean(df_YA_ascrfrac_PC_traits$PC3)-6*sd(df_YA_ascrfrac_PC_traits$PC3)), aes(label = strain), size = 2.5) +
  #geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, strain %in% c("ED3052", "JU1242")), aes(label = strain), size = 2.5) +
  aes(x=PC2, y=PC3) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC2 (15.7%)", y="PC3 (8.6%)")

plot_PC23_YA

plot_PC123_YA <- cowplot::plot_grid(plot_PC12_YA, plot_PC13_YA, plot_PC23_YA, nrow = 1)
plot_PC123_YA

ggsave(plot_PC123_YA, file = "Plots/plot_PC123_YA.pdf", width = 7.5, height = 2.7)

plot_PC123_YA_v <- cowplot::plot_grid(plot_PC12_YA, plot_PC13_YA, plot_PC23_YA, ncol = 1)
ggsave(plot_PC123_YA_v, file = "Plots/plot_PC123_YA_v.pdf", width = 2.5, height = 7)

### PCA after outlier (JU1242) removal ###

df_YA_pheno_or <- ascr_YA_frac %>%
  dplyr::filter(!strain %in% c("JU1400", "JU1242")) %>%
  dplyr::select(strain, feature, ascr_fraction) %>%
  tidyr::spread(key=strain, value=ascr_fraction)

df_YA_pheno_wide_or <- t(df_YA_pheno_or[2:ncol(df_YA_pheno_or)])
colnames(df_YA_pheno_wide_or) <- df_YA_pheno_or[[1]]

### PCA and outlier detection ###

strain_YA_or<-row.names(df_YA_pheno_wide_or)

df_YA_ascrfrac_PC_trim.fit_scaled_or <- prcomp(df_YA_pheno_wide_or, scale. = TRUE)
df_YA_PCA_summary_or <- data.frame(t(summary(df_YA_ascrfrac_PC_trim.fit_scaled_or)[[6]]), PC=c(1:44))

plot_YA_PCA_summary_or <- df_YA_PCA_summary_or %>%
  ggplot(.) +
  geom_col() +
  aes(x=PC, y=Proportion.of.Variance) +
  theme_bw()

plot_YA_PCA_summary_cum_or <- df_YA_PCA_summary_or %>%
  ggplot(.) +
  geom_col() +
  aes(x=PC, y=Cumulative.Proportion) +
  theme_bw()

ggsave(plot_YA_PCA_summary, file = "Plots/plot_YA_PCA_summary_or.pdf", width = 7.5, height = 3.5)
ggsave(plot_YA_PCA_summary_cum, file = "Plots/plot_YA_PCA_summary_cum_or.pdf", width = 7.5, height = 3.5)

## plot for loadings

df_YA_ascrfrac_PC_loading_or<-data.frame(feature=rownames(df_YA_ascrfrac_PC_trim.fit_scaled_or[[2]]), df_YA_ascrfrac_PC_trim.fit_scaled_or[[2]]) %>%
  tidyr::gather(key="PC", value="loading", -feature)

# group1

plot_YA_PC1_loading1_or <- df_YA_ascrfrac_PC_loading_or %>%
  dplyr::filter(PC == "PC1") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='none') +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC1")

plot_YA_PC1_loading1_or

plot_YA_PC2_loading1_or <- df_YA_ascrfrac_PC_loading_or %>%
  dplyr::filter(PC == "PC2") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='none') +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC2")

plot_YA_PC3_loading1_or <- df_YA_ascrfrac_PC_loading_or %>%
  dplyr::filter(PC == "PC3") %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color='black', size=0.1) +
  aes(x=reorder(feature, loading), y=loading, fill=group1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        legend.position='bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.15,"in"),
        legend.text = element_text(size=9, color='black')) +
  scale_fill_manual(values = group_pal) +
  ggtitle("PC3") +
  guides(fill= guide_legend(nrow=2))

plot_YA_PC1_loading_group1_or <- cowplot::plot_grid(plot_YA_PC1_loading1_or, plot_YA_PC2_loading1_or, plot_YA_PC3_loading1_or, 
                                                 ncol = 1, rel_heights = c(1,1,1.3))

ggsave(plot_YA_PC1_loading_group1_or, file = "Plots/plot_YA_PC1_loading_group1_or.pdf", width = 5, height = 7)

df_YA_ascrfrac_PC_traits_or <-df_YA_ascrfrac_PC_trim.fit_scaled_or[[5]]
df_YA_ascrfrac_PC_traits_or <- data.frame(strain=strain_YA_or, df_YA_ascrfrac_PC_traits_or)

plot_PC12_YA_or <- df_YA_ascrfrac_PC_traits_or %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_or, PC1 > mean(df_YA_ascrfrac_PC_traits_or$PC1)+3*sd(df_YA_ascrfrac_PC_traits_or$PC1) |
                                          PC1 < mean(df_YA_ascrfrac_PC_traits_or$PC1)-3*sd(df_YA_ascrfrac_PC_traits_or$PC1)), aes(label = strain), size = 2.5, nudge_x = -1) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_or, PC2 > mean(df_YA_ascrfrac_PC_traits_or$PC2)+3*sd(df_YA_ascrfrac_PC_traits_or$PC2) |
                                          PC2 < mean(df_YA_ascrfrac_PC_traits_or$PC2)-3*sd(df_YA_ascrfrac_PC_traits_or$PC2)), aes(label = strain), size = 2.5, nudge_x = -1) +
  #geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_or, strain %in% c("ED3052", "JU1242")), aes(label = strain), size = 2.5) +
  aes(x=PC1, y=PC2) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC1 (20.6%)", y="PC2 (11.7%)")

plot_PC12_YA_or

plot_PC13_YA_or <- df_YA_ascrfrac_PC_traits_or %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_or, PC1 > mean(df_YA_ascrfrac_PC_traits_or$PC1)+3*sd(df_YA_ascrfrac_PC_traits_or$PC1) |
                                          PC1 < mean(df_YA_ascrfrac_PC_traits_or$PC1)-3*sd(df_YA_ascrfrac_PC_traits_or$PC1)), aes(label = strain), size = 2.5, nudge_x = -1) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_or, PC3 > mean(df_YA_ascrfrac_PC_traits_or$PC3)+3*sd(df_YA_ascrfrac_PC_traits_or$PC3) |
                                          PC3 < mean(df_YA_ascrfrac_PC_traits_or$PC3)-3*sd(df_YA_ascrfrac_PC_traits_or$PC3)), aes(label = strain), size = 2.5, nudge_x = -1) +
  #geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_or, strain %in% c("ED3052", "JU1242")), aes(label = strain), size = 2.5) +
  aes(x=PC1, y=PC3) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC1 (20.6%)", y="PC3 (6.5%)") +
  scale_y_continuous(breaks=c(-6,-3,0,3))

plot_PC13_YA_or

plot_PC23_YA_or <- df_YA_ascrfrac_PC_traits_or %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_or, PC2 > mean(df_YA_ascrfrac_PC_traits_or$PC2)+3*sd(df_YA_ascrfrac_PC_traits_or$PC2) |
                                          PC2 < mean(df_YA_ascrfrac_PC_traits_or$PC2)-3*sd(df_YA_ascrfrac_PC_traits_or$PC2) |
                                          PC3 > mean(df_YA_ascrfrac_PC_traits_or$PC3)+3*sd(df_YA_ascrfrac_PC_traits_or$PC3) |
                                          PC3 < mean(df_YA_ascrfrac_PC_traits_or$PC3)-3*sd(df_YA_ascrfrac_PC_traits_or$PC3)), aes(label = strain), size = 2.5, nudge_x = -1) +
  aes(x=PC2, y=PC3) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC2 (11.7%)", y="PC3 (6.5%)")+
  scale_y_continuous(breaks=c(-6,-3,0,3))

plot_PC23_YA_or

plot_PC123_YA_or <- cowplot::plot_grid(plot_PC12_YA_or, plot_PC13_YA_or, plot_PC23_YA_or, nrow = 1)
plot_PC123_YA_or

ggsave(plot_PC123_YA_or, file = "Plots/plot_PC123_YA_or.pdf", width = 7.5, height = 2.7)

plot_PC123_YA_or_v <- cowplot::plot_grid(plot_PC12_YA_or, plot_PC13_YA_or, plot_PC23_YA_or, ncol = 1)
ggsave(plot_PC123_YA_or_v, file = "Plots/plot_PC123_YA_or_v.pdf", width = 2.5, height = 7)

### Pie chart for outliers and N2 ###

ascr_YA_frac_select <- ascr_YA_frac %>%
  dplyr::filter(strain %in% c("N2","CB4856", "JU1242", "ED3052", "WN2001")) %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  dplyr::group_by(group1, strain) %>%
  dplyr::summarise(group_fraction = sum(ascr_fraction)) %>%
  dplyr::ungroup()

ascr_YA_frac_select$strain <- factor(ascr_YA_frac_select$strain, levels = c("N2","CB4856", "JU1242", "ED3052", "WN2001"))

plot_5strains_pie <- ascr_YA_frac_select %>%
  ggplot(.) +
  geom_col(color= 'black', size = 0.1) +
  aes(x="", y=group_fraction, fill=group1) +
  #geom_label(data = dplyr::filter(ascr_YA_frac_select, ascr_fraction >= 0.05), aes(label = group2), size = 2) +
  coord_polar("y", start=0) +
  facet_grid(~strain) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.key.size = unit(0.1,"in"),
        legend.text = element_text(size=8, color='black'),
        panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = group_pal) +
  guides(fill= guide_legend(nrow=2))

plot_5strains_pie

ggsave(plot_5strains_pie, file = "Plots/plot_5strains_pie_label.pdf", width = 7.5, height = 2)

### Pie chart for all strains ###

ascr_YA_frac_select_others <- ascr_YA_frac %>%
  dplyr::filter(!strain %in% c("N2","JU1242", "ED3052", "CB4856", "WN2001", "JU1400")) %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  dplyr::group_by(group1, strain) %>%
  dplyr::summarise(group_fraction = sum(ascr_fraction)) %>%
  dplyr::ungroup()

plot_otherstrains_pie <- ascr_YA_frac_select_others %>%
  ggplot(.) +
  geom_col(color= 'black', size = 0.1) +
  aes(x="", y=group_fraction, fill=group1) +
  #geom_label(data = dplyr::filter(ascr_YA_frac_select, ascr_fraction >= 0.05), aes(label = group2), size = 2) +
  coord_polar("y", start=0) +
  facet_wrap(~strain, ncol=8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.1,"in"),
        legend.text = element_text(size=8, color='black'),
        panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = group_pal) +
  guides(fill= guide_legend(nrow=2))

plot_otherstrains_pie

ggsave(plot_otherstrains_pie, file = "Plots/plot_otherstrains_pie.pdf", width = 7.5, height = 15)

### pie chart for interesting strains

ascr_YA_frac_select_interest <- ascr_YA_frac %>%
  dplyr::filter(strain %in% c("ECA36", "JU258", "JU775", "PB303", "QG557")) %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  dplyr::group_by(group1, strain) %>%
  dplyr::summarise(group_fraction = sum(ascr_fraction)) %>%
  dplyr::ungroup()

plot_intstrains_pie <- ascr_YA_frac_select_interest %>%
  ggplot(.) +
  geom_col(color= 'black', size = 0.1) +
  aes(x="", y=group_fraction, fill=group1) +
  #geom_label(data = dplyr::filter(ascr_YA_frac_select, ascr_fraction >= 0.05), aes(label = group2), size = 2) +
  coord_polar("y", start=0) +
  facet_wrap(~strain, ncol=8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.1,"in"),
        legend.text = element_text(size=8, color='black'),
        panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = group_pal) +
  guides(fill= guide_legend(nrow=2))

plot_intstrains_pie

ggsave(plot_intstrains_pie, file = "Plots/plot_intstrains_pie.pdf", width = 7.5, height = 5)

### Table for cegwas - remove JU1400

df_GWA_fraction_YA <- data.frame(strain = rownames(df_YA_pheno_wide), df_YA_pheno_wide) %>%
  dplyr::filter(!strain %in% c("JU1400"))

pheno_strains_YA <- unique(df_GWA_fraction_YA$strain)

write_tsv(df_GWA_fraction_YA, path="Processed_Data/ascr_fraction_YA_OutlierRemoval.tsv", col_names = T)

save(pheno_strains_YA, ascr_YA_frac, df_ascr_group, df_YA_pheno_wide, df_GWA_fraction_YA, file="Processed_Data/df_GWA.RData")


View(gsub("PB306", "ECA259", df_GWA_fraction_YA))
