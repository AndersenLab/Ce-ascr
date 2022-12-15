library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggrepel)
library(sommer)

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

ascr_YA_frac$feature <- gsub("bhas#11", "bhos#11", ascr_YA_frac$feature) ## ascr name fix

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
    labs(y="Relative abundance (%)")
  
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
bhos = c("bhos#10", "bhos#11", "bhos#18", "bhos#22")
four_modified = c("icas#3", "icas#9", "icas#10", "osas#2", "osas#9", "osas#10")
C_modified = c("iglas#2", "glas#1", "glas#3", "anglas#3", "anglas#7", "ascr#81", "ascr#8", "uglas#11")
bhas = c("bhas#10", "bhas#12","bhas#18", "bhas#22", "bhas#16", "bhas#20", "bhas#26")

df_ascr_group <- data.frame(feature = c(as.character(unique(df_YA_pheno$feature)))) %>%
  dplyr::mutate(group1 = ifelse(feature %in% ascr_short, "ascr(s)",
                                ifelse(feature %in% ascr_long, "ascr(l)",
                                       ifelse(feature %in% oscr_short, "oscr(s)",
                                              ifelse(feature %in% oscr_long, "oscr(l)",
                                                     ifelse(feature %in% bhos, "bhos",
                                                            ifelse(feature %in% four_modified, "modified(4')",
                                                                   ifelse(feature %in% C_modified, "modified(C)", 
                                                                          ifelse(feature %in% bhas, "bhas", "etc"))))))))) %>%
  dplyr::group_by(group1) %>%
  dplyr::mutate(n_member_group1=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group2 = ifelse(feature %in% c(ascr_short, ascr_long), "ascr",
                                ifelse(feature %in% c(oscr_short, oscr_long), "oscr",group1)))

df_ascr_group$group1 <- factor(df_ascr_group$group1, levels = c("ascr(s)", "ascr(l)", "bhas","modified(C)", "modified(4')", "oscr(s)", "oscr(l)", "bhos"))

group_pal <- c("#9E0142", "#F46D43", "#FDAE61", "#A6761D", "#666666", "#3288BD", "#5E4FA2", "#66C2A5")

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

### Table for cegwas - remove JU1400

df_GWA_fraction_YA <- data.frame(strain = rownames(df_YA_pheno_wide), df_YA_pheno_wide*1e6)

pheno_strains_YA <- unique(df_GWA_fraction_YA$strain)

write_tsv(df_GWA_fraction_YA, path="Processed_Data/ascr_fraction_YA_rescaled.tsv", col_names = T)

save(pheno_strains_YA, ascr_YA_frac, df_ascr_group, group_pal, df_YA_pheno_wide, df_GWA_fraction_YA, file="Processed_Data/df_GWA.RData")

df_a3a5_ratios <- df_GWA_fraction_YA %>%
  dplyr::mutate(a3a5=ascr.3/ascr.5, a5a3=ascr.5/ascr.3) %>%
  dplyr::select(strain, a3a5, a5a3)

write_tsv(df_a3a5_ratios, path="Processed_Data/ascr35_ratio_YA.tsv", col_names = T)

df_bhas <- df_GWA_fraction_YA %>%
  dplyr::select(strain, "bhas.12", "bhas.10", "bhas.18", "bhas.22") %>%
  dplyr::mutate(bhas.12=bhas.12*1e5, bhas.10=bhas.10*1e5, bhas.18=bhas.18*1e5, bhas.22=bhas.22*1e5)

write_tsv(df_bhas, path="Processed_Data/bhas_YA.tsv", col_names = T)

min(df_YA_pheno_wide * 1e6)

### Heritability ###

geno_matrix <- read.table(file = "Raw/Analysis_Results-20211018_EIGEN/Genotype_Matrix/Genotype_Matrix.tsv", header = T)
A <- sommer::A.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains_YA]))

df_h2_YA <- data.frame(feature = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_GWA_fraction_YA)-1)) {
  
  feature <- colnames(df_GWA_fraction_YA)[i+1]
  
  df_y <- df_GWA_fraction_YA %>%
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=feature) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(feature, h2, h2_SE)
  
  df_h2_YA <- rbind(df_h2_YA, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

df_h2_YA$feature <- gsub("bhas#11", "bhos#11", df_h2_YA$feature) ## ascr name fix

View(df_h2_YA)

plot_frac_h2_YA_histo <- df_h2_YA %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_frac_h2_YA_histo

ggsave(file="Plots/plot_frac_YA_h2_histo.png", plot=plot_frac_h2_YA_histo, width = 7.5, height = 5, units = "in")

df_ascr_CN <- read.csv(file="Raw/ascr_n_carbon.csv")

df_h2_YA$feature <- gsub("\\.", "#", df_h2_YA$feature)

df_ascr_YA_frac_summary <- ascr_YA_frac %>%
  dplyr::left_join(., df_ascr_CN, by='feature') %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::group_by(feature, carbons) %>%
  dplyr::summarise(median = median(ascr_fraction)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., df_h2_YA, by="feature")

plot_cn_h2 <- df_ascr_YA_frac_summary %>%
  dplyr::left_join(., df_ascr_group, by = 'feature') %>%
  ggplot(.) +
  geom_point(alpha=0.7, shape=21, size=2) +
  aes(x=carbons, y=as.numeric(h2)*100, fill=group1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=9, color='black'),
        axis.title = element_text(size=10, color='black'),
        legend.title = element_blank(),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.position=c(0.8,0.75),
        legend.background = element_blank(),
        legend.text = element_text(size=8.5, color='black')) +
  geom_hline(yintercept = 10, color='blue', linetype = 2, alpha=0.7) +
  labs(x="FA length (Cn)", y="Narrow-sense heritability (%)") +
  scale_fill_manual(values = group_pal)

plot_cn_h2

ggsave(file="Plots/plot_cv_h2.png", plot=plot_cv_h2, width = 5, height = 5, units = "in")

df_ascr_YA_frac_summary_high_h2 <- df_ascr_YA_frac_summary %>%
  dplyr::filter(h2>0.1)

### Table for cegwas - remove JU1400 / h2 filtered

df_YA_pheno_h2_filtered <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::filter(feature %in% df_ascr_YA_frac_summary_high_h2$feature) %>%
  dplyr::select(strain, feature, ascr_fraction) %>%
  tidyr::spread(key=strain, value=ascr_fraction) 

df_YA_pheno_wide_h2_filtered <- t(df_YA_pheno_h2_filtered[2:ncol(df_YA_pheno_h2_filtered)])
colnames(df_YA_pheno_wide_h2_filtered) <- df_YA_pheno_h2_filtered[[1]]

df_GWA_fraction_YA_h2_filtered <- data.frame(strain = rownames(df_YA_pheno_wide_h2_filtered), df_YA_pheno_wide_h2_filtered)

write_tsv(df_GWA_fraction_YA_h2_filtered, path="Processed_Data/ascr_fraction_YA_h2_filtered.tsv", col_names = T)

save(pheno_strains_YA, ascr_YA_frac, df_ascr_YA_frac_summary_high_h2, 
     df_ascr_group, df_YA_pheno_wide, df_YA_pheno_wide_h2_filtered, 
     df_GWA_fraction_YA, df_GWA_fraction_YA_h2_filtered, 
     df_h2_YA, df_ascr_YA_frac_summary, group_pal, 
     file="Processed_Data/df_ascr_data_norm_h2.RData")

### ascr#3 vs ascr#5

cor(df_YA_pheno_wide_h2_filtered[,2],df_YA_pheno_wide_h2_filtered[,3], method = 'spearman') ## ascr#3 vs ascr#5

plot_ascr3_ascr5 <- data.frame(ascr.3=df_YA_pheno_wide_h2_filtered[,2], ascr.5=df_YA_pheno_wide_h2_filtered[,3]) %>%
  ggplot(.) +
  geom_point() +
  aes(x=ascr.3, y=ascr.5) +
  theme_bw() +
  labs(x="ascr#3", y="ascr#5")

ggsave(file="Plots/plot_ascr3_vs_ascr_5.png", plot=plot_ascr3_ascr5, width = 5, height = 5, units = "in")

## heritabiliy

df_y_a3a5 <- data.frame(strain=rownames(df_YA_pheno_wide_h2_filtered),
                        a3_a5=df_YA_pheno_wide_h2_filtered[,2]/df_YA_pheno_wide_h2_filtered[,3]) %>%
  dplyr::arrange(strain) %>%
  dplyr::select(strain, value=a3_a5) %>%
  dplyr::mutate(strain = as.character(strain))

h2_res_a3a5 <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y_a3a5)

h2_a3a5 <- as.numeric(pin(h2_res_a3a5, h2_a3a5 ~ (V1) / (V1+V2))[[1]][1])
h2_SE_a3a5 <- pin(h2_res_a3a5, h2_a3a5 ~ (V1) / (V1+V2))[[2]]
h2_a3a5


### PCA and outlier detection ###

ascr_YA_frac_h2_filtered <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::filter(feature %in% df_ascr_YA_frac_summary_high_h2$feature) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(ascr_sum_h2_filtered = sum(abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, feature, abundance, .keep_all = T) %>%
  dplyr::mutate(ascr_fraction_h2_filtered = abundance/ascr_sum_h2_filtered)

df_YA_pheno_h2_filtered <- ascr_YA_frac_h2_filtered %>%
  dplyr::select(strain, feature, ascr_fraction_h2_filtered) %>%
  tidyr::spread(key=strain, value=ascr_fraction_h2_filtered)

df_YA_pheno_wide_h2_filtered <- t(df_YA_pheno_h2_filtered[2:ncol(df_YA_pheno_h2_filtered)])
colnames(df_YA_pheno_wide_h2_filtered) <- df_YA_pheno_h2_filtered[[1]]

strain_YA_h2_filtered<-row.names(df_YA_pheno_wide_h2_filtered)

df_YA_ascrfrac_PC_trim.fit_scaled_h2_filtered <- prcomp(df_YA_pheno_wide_h2_filtered, scale. = TRUE)
df_YA_PCA_summary_h2_filtered <- data.frame(t(summary(df_YA_ascrfrac_PC_trim.fit_scaled_h2_filtered)[[6]]), PC=c(1:23))

plot_YA_PCA_summary_h2_filtered <- df_YA_PCA_summary_h2_filtered %>%
  ggplot(.) +
  geom_col() +
  aes(x=PC, y=Proportion.of.Variance) +
  theme_bw()

plot_YA_PCA_summary_cum_h2_filtered <- df_YA_PCA_summary_h2_filtered %>%
  ggplot(.) +
  geom_col() +
  aes(x=PC, y=Cumulative.Proportion) +
  theme_bw()

ggsave(plot_YA_PCA_summary_h2_filtered, file = "Plots/plot_YA_PCA_summary_h2_filtered.png", width = 7.5, height = 3.5)
ggsave(plot_YA_PCA_summary_cum_h2_filtered, file = "Plots/plot_YA_PCA_summary_cum_h2_filtered.png", width = 7.5, height = 3.5)


#df_ascr_group$group1 <- factor(df_ascr_group$group, levels = c("ascr_short", "oscr_short", "bhas_short", "four_modified", "C_modified", "bhos", "ascr_long", "oscr_long","bhas_long", "etc"))

## plot for loadings

df_YA_ascrfrac_PC_loading_h2_filtered<-data.frame(feature=rownames(df_YA_ascrfrac_PC_trim.fit_scaled_h2_filtered[[2]]), 
                                                  df_YA_ascrfrac_PC_trim.fit_scaled_h2_filtered[[2]]) %>%
  tidyr::gather(key="PC", value="loading", -feature)

# group1

plot_YA_PC1_loading1_h2_filtered <- df_YA_ascrfrac_PC_loading_h2_filtered %>%
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

plot_YA_PC1_loading1_h2_filtered

plot_YA_PC2_loading1_h2_filtered <- df_YA_ascrfrac_PC_loading_h2_filtered %>%
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

plot_YA_PC2_loading1_h2_filtered

plot_YA_PC3_loading1_h2_filtered <- df_YA_ascrfrac_PC_loading_h2_filtered %>%
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

plot_YA_PC3_loading1_h2_filtered

## group2

plot_YA_PC1_loading2_h2_filtered <- df_YA_ascrfrac_PC_loading_h2_filtered %>%
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

plot_YA_PC1_loading2_h2_filtered

plot_YA_PC2_loading2_h2_filtered <- df_YA_ascrfrac_PC_loading_h2_filtered %>%
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

plot_YA_PC2_loading2_h2_filtered

plot_YA_PC3_loading2_h2_filtered <- df_YA_ascrfrac_PC_loading_h2_filtered %>%
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

plot_YA_PC3_loading2_h2_filtered

plot_YA_PC1_loading_group1_h2_filtered <- cowplot::plot_grid(plot_YA_PC1_loading1_h2_filtered, 
                                                             plot_YA_PC2_loading1_h2_filtered, 
                                                             plot_YA_PC3_loading1_h2_filtered, 
                                                             ncol = 1, rel_heights = c(1,1,1.2))
plot_YA_PC1_loading_group2_h2_filtered <- cowplot::plot_grid(plot_YA_PC1_loading2_h2_filtered, 
                                                             plot_YA_PC2_loading2_h2_filtered, 
                                                             plot_YA_PC3_loading2_h2_filtered, 
                                                             ncol = 1, rel_heights = c(1,1,1.2))

#ggsave(plot_YA_PC1_loading_group1_h2_filtered, file = "Plots/plot_YA_PC1_loading_group1_h2_filtered.pdf", width = 7.5, height = 9)
ggsave(plot_YA_PC1_loading_group2_h2_filtered, file = "Plots/plot_YA_PC1_loading_group2_h2_filtered.png", width = 7.5, height = 9)

df_YA_ascrfrac_PC_traits_h2_filtered <-df_YA_ascrfrac_PC_trim.fit_scaled_h2_filtered[[5]]
#pc_traits <-data.frame(pcDF.fit[[5]])
df_YA_ascrfrac_PC_traits_h2_filtered <- data.frame(strain=strain_YA_h2_filtered, df_YA_ascrfrac_PC_traits_h2_filtered)

plot_PC12_YA_h2_filtered <- df_YA_ascrfrac_PC_traits_h2_filtered %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_h2_filtered, strain %in% c("JU1242")), 
                   aes(label = strain), size = 2.5) +
  aes(x=PC1, y=PC2) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC1 (40.3%)", y="PC2 (12.1%)")

plot_PC12_YA_h2_filtered

plot_PC13_YA_h2_filtered <- df_YA_ascrfrac_PC_traits_h2_filtered %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_h2_filtered, PC3 > 3| PC3 < -4), aes(label = strain), size = 2.5) +
  #geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits, strain %in% c("ED3052", "JU1242")), aes(label = strain), size = 2.5) +
  aes(x=PC1, y=PC3) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC1 (19.7%)", y="PC3 (8.6%)")

plot_PC13_YA_h2_filtered

plot_PC23_YA_h2_filtered <- df_YA_ascrfrac_PC_traits_h2_filtered %>%
  ggplot(.) +
  geom_point(size = .8, alpha = 0.7) +
  geom_label_repel(data = dplyr::filter(df_YA_ascrfrac_PC_traits_h2_filtered, PC3 > 3), aes(label = strain), size = 2.5) +
  aes(x=PC2, y=PC3) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        panel.grid = element_blank()) +
  labs(x="PC2 (15.7%)", y="PC3 (8.6%)")

plot_PC23_YA_h2_filtered

plot_PC123_YA_h2_filtered <- cowplot::plot_grid(plot_PC12_YA_h2_filtered, plot_PC13_YA_h2_filtered, plot_PC23_YA_h2_filtered, nrow = 1)
plot_PC123_YA_h2_filtered

ggsave(plot_PC123_YA_h2_filtered, file = "Plots/plot_PC123_YA_h2_filtered.png", width = 7.5, height = 2.7)

plot_PC123_YA_v_h2_filtered <- cowplot::plot_grid(plot_PC12_YA_h2_filtered, plot_PC13_YA_h2_filtered, plot_PC23_YA_h2_filtered, ncol = 1)
ggsave(plot_PC123_YA_v_h2_filtered, file = "Plots/plot_PC123_YA_v_h2_filtered.png", width = 2.5, height = 7)


save(df_YA_PCA_summary_h2_filtered, df_YA_ascrfrac_PC_loading_h2_filtered, df_YA_ascrfrac_PC_traits_h2_filtered,
     file="Processed_Data/df_PCA.RData")


## Piechart for all strains

plot_allstrains_h2_filtered_pie <- df_ascr_YA_frac_h2_filtered %>%
  ggplot(.) +
  geom_col(color= 'black', size = 0.1) +
  aes(x="", y=group_fraction, fill=group2) +
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
  guides(fill= guide_legend(nrow=1))

plot_allstrains_h2_filtered_pie

ggsave(plot_allstrains_h2_filtered_pie, file = "Plots/plot_allstrains_pie_h2_filtered.png", width = 7.5, height = 15)

ascr_YA_frac_select <- ascr_YA_frac %>%
  dplyr::filter(strain %in% c("N2","ED3052", "EG4349","JU1242")) %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  dplyr::group_by(group1, strain) %>%
  dplyr::summarise(group_fraction = sum(ascr_fraction)) %>%
  dplyr::ungroup()

ascr_YA_frac_select$strain <- factor(ascr_YA_frac_select$strain, levels = c("N2","ED3052", "EG4349","JU1242"))

plot_4strains_pie <- ascr_YA_frac_select %>%
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
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.1,"in"),
        legend.text = element_text(size=8, color='black'),
        panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = group_pal) +
  guides(fill= guide_legend(nrow=2))

plot_4strains_pie

ascr_YA_frac_h2_filtered <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::filter(feature %in% df_ascr_YA_frac_summary_high_h2$feature) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(ascr_sum_h2_filtered = sum(abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, feature, abundance, .keep_all = T) %>%
  dplyr::mutate(ascr_fraction_h2_filtered = abundance/ascr_sum_h2_filtered)

plot_4strains_pie_h2_filtered <- ascr_YA_frac_h2_filtered %>%
  dplyr::filter(strain %in% c("N2","ED3052", "EG4349","JU1242")) %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  ggplot(.) +
  geom_col(color= 'black', size = 0.1) +
  aes(x="", y=ascr_fraction_h2_filtered, fill=group1) +
  #geom_label(data = dplyr::filter(ascr_YA_frac_select, ascr_fraction >= 0.05), aes(label = group2), size = 2) +
  coord_polar("y", start=0) +
  facet_grid(~factor(strain, levels=c("N2","ED3052", "EG4349","JU1242"))) +
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
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FEE08B", "#A6761D", "#666666", "#3288BD", "#5E4FA2")) +
  guides(fill= guide_legend(nrow=2))

plot_4strains_pie_h2_filtered
