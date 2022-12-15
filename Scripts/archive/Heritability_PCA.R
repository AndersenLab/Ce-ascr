library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(sommer)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_GWA.RData")

### Heritability ###

geno_matrix <- read.table(file = "Processed_Data/Analysis_Results-20210225_total_abundance_bpTIC_normalized/Genotype_Matrix/Genotype_Matrix.tsv", header = T)
A <- sommer::A.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains_YA]))
E <- sommer::E.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains_YA]))

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

View(df_h2_YA)

plot_frac_h2_YA_histo <- df_h2_YA %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_frac_h2_YA_histo

ggsave(file="Plots/plot_frac_YA_h2_histo.png", plot=plot_frac_h2_YA_histo, width = 7.5, height = 5, units = "in")

df_h2_YA$feature <- gsub("\\.", "#", df_h2_YA$feature)

df_ascr_YA_frac_summary <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(median = median(ascr_fraction), cv=sd(ascr_fraction)/mean(ascr_fraction)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., df_h2_YA, by="feature")

plot_cv_h2 <- df_ascr_YA_frac_summary %>%
  ggplot(.) +
  geom_point() +
  aes(x=cv, y=as.numeric(h2)) +
  theme_bw() +
  labs(x="Coeffecient of variation", y="Narrow-sense heritability")

plot_cv_h2

ggsave(file="Plots/plot_cv_h2.png", plot=plot_cv_h2, width = 5, height = 5, units = "in")

cor(df_ascr_YA_frac_summary$cv, as.numeric(df_ascr_YA_frac_summary$h2), method='spearman')

plot_median_h2 <- df_ascr_YA_frac_summary %>%
  ggplot(.) +
  geom_point() +
  aes(x=log10(median), y=as.numeric(h2)) +
  theme_bw() +
  labs(x="log10(Median relative abundance)", y="Narrow-sense heritability")

plot_median_h2

ggsave(file="Plots/plot_median_h2.png", plot=plot_median_h2, width = 5, height = 5, units = "in")

cor(df_ascr_YA_frac_summary$median, as.numeric(df_ascr_YA_frac_summary$h2), method='spearman')

df_ascr_YA_frac_summary_high_h2 <- df_ascr_YA_frac_summary %>%
  dplyr::filter(h2>0.2)


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
df_YA_PCA_summary_h2_filtered <- data.frame(t(summary(df_YA_ascrfrac_PC_trim.fit_scaled_h2_filtered)[[6]]), PC=c(1:16))

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
  labs(x="PC1 (19.7%)", y="PC2 (15.7%)")

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

df_ascr_YA_frac_h2_filtered <- ascr_YA_frac_h2_filtered %>%
  dplyr::left_join(., df_ascr_group, by='feature') %>%
  dplyr::group_by(group2, strain) %>%
  dplyr::summarise(group_fraction = sum(ascr_fraction_h2_filtered)) %>%
  dplyr::ungroup()

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


### pairwise-ratio ###

ratio <- data.frame(combn(as.character(unique(ascr_YA_frac$feature)),m=2))

ratio_list_YA<-list()

for (i in 1:ncol(ratio)) {
  as1 <- dplyr::filter(ascr_YA_frac,feature==as.character(ratio[1,i]))%>%
    rename(as1 = feature,amount1 = abundance)
  as2 <- dplyr::filter(ascr_YA_frac,feature==as.character(ratio[2,i]))%>%
    rename(as2 = feature,amount2 = abundance)
  comb.as <- left_join(as1,as2,by="strain")%>%
    dplyr::mutate(ratio1= amount1/amount2,
                  ratio2= amount2/amount1)%>%
    dplyr::select(strain, ratio1, ratio2)
  
  colnames(comb.as) <- c("strain",
                         paste(as.character(ratio[1,i]),as.character(ratio[2,i]),sep="/"),
                         paste(as.character(ratio[2,i]),as.character(ratio[1,i]),sep="/"))
  as.ratio.long <- tidyr::gather(comb.as,trait,value,-strain)
  ratio_list_YA[[i]] <- as.ratio.long
}

all_ratio_YA <- bind_rows(ratio_list_YA)


### pair-wise ratio h2 ###

df_ratio_wide <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  tidyr::spread(trait,value)

write_tsv(df_ratio_wide, path="Processed_Data/ascr_ratio_YA_OutlierRemoval.tsv", col_names = T)

df_ratio_h2_YA <- data.frame(trait = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_ratio_wide)-1)) {
  
  trait <- colnames(df_ratio_wide)[i+1]
  
  df_y <- df_ratio_wide %>%
    dplyr::filter(strain != "JU1516") %>%  ### same isotype with JU1581
    dplyr::filter(strain != "JU1400") %>% ### Exclude extreme outlier
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=trait) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(trait, h2, h2_SE)
  
  df_ratio_h2_YA <- rbind(df_ratio_h2_YA, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

View(df_ratio_h2_YA)

plot_ratio_h2_YA_histo <- df_ratio_h2_YA %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_ratio_h2_YA_histo

ggsave(file="Plots/plot_ratio_YA_h2_histo.pdf", plot=plot_ratio_h2_YA_histo, width = 7.5, height = 5, units = "in")

### heritability heatmap ###

df_ratio_h2_wide <- df_ratio_h2_YA %>%
  dplyr::select(-h2_SE) %>%
  dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=h2)

rownames(df_ratio_h2_wide) <- df_ratio_h2_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_h2_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
           legend.title = "Narrow-sense\nheritability\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.4,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Narrow-sense\nheritability\n")

ggsave("Plots/plot_ratio_h2_heatmap.pdf", width = 7.5, height = 7.5, units = "in")

### example plot for high heritability to check what's happening ###

ECA36_relatedness <- A %>%
  as.data.frame %>%
  dplyr::select(relateness_ECA36=ECA36) %>%
  dplyr::add_rownames(var='strain')

df_example_ratio <- all_ratio_YA %>%
  dplyr::filter(trait == "ascr#10/icas#3") %>%
  dplyr::left_join(., ECA36_relatedness, by='strain') %>%
  na.omit()
  
df_example_ratio %>%
  ggplot(.) +
  geom_point() +
  aes(x=relateness_ECA36, y=value) +
  theme_bw()

### after removal of ECA36 and DL238 ###

### pair-wise ratio h2 ###

df_ratio_wide2 <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  dplyr::filter(!strain %in% c("ECA36", "DL238")) %>%
  tidyr::spread(trait,value)

df_ratio_h2_YA2 <- data.frame(trait = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_ratio_wide2)-1)) {
  
  trait <- colnames(df_ratio_wide2)[i+1]
  
  df_y <- df_ratio_wide2 %>%
    dplyr::filter(strain != "JU1516") %>%  ### same isotype with JU1581
    dplyr::filter(strain != "JU1400") %>% ### Exclude extreme outlier
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=trait) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(trait, h2, h2_SE)
  
  df_ratio_h2_YA2 <- rbind(df_ratio_h2_YA2, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

View(df_ratio_h2_YA2)

save(df_h2_YA, df_ratio_h2_YA, df_ratio_h2_YA2, file="Processed_Data/df_heritability.RData")

plot_ratio_h2_YA_histo2 <- df_ratio_h2_YA2 %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_ratio_h2_YA_histo2

ggsave(file="Plots/plot_ratio_YA_h2_histo_outlier_removed.pdf", plot=plot_ratio_h2_YA_histo2, width = 7.5, height = 5, units = "in")

### heritability heatmap ###

df_ratio_h2_wide2 <- df_ratio_h2_YA2 %>%
  dplyr::select(-h2_SE) %>%
  dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=h2)

rownames(df_ratio_h2_wide2) <- df_ratio_h2_wide2$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_h2_wide2, -trait2)), hc.order = FALSE, tl.srt = 90, 
           legend.title = "Narrow-sense\nheritability\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.4,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Narrow-sense\nheritability\n")

ggsave("Plots/plot_ratio_h2_heatmap_outlier_removed.pdf", width = 7.5, height = 7.5, units = "in")


### example plot for high heritability to check what's happening ###

N2_relatedness <- A %>%
  as.data.frame %>%
  dplyr::select(N2_relatedness=N2) %>%
  dplyr::add_rownames(var='strain')

df_example_ratio2 <- all_ratio_YA %>%
  dplyr::filter(trait == "ascr#3/icas#3", !strain %in% c("DL238", "ECA36")) %>%
  dplyr::left_join(., N2_relatedness, by='strain') %>%
  na.omit()

df_example_ratio2 %>%
  ggplot(.) +
  geom_point() +
  aes(x=N2_relatedness, y=value) +
  theme_bw()


### after outlier removal by 3 standard deviation ###

### pair-wise ratio h2 ###

df_ratio_wide3 <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(mean=mean(value), sd=sd(value)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(value>=mean-(3*sd) & value <= mean+(3*sd)) %>%
  dplyr::select(-mean, -sd) %>%
  tidyr::spread(trait,value)

df_ratio_h2_YA3 <- data.frame(trait = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_ratio_wide3)-1)) {
  
  trait <- colnames(df_ratio_wide3)[i+1]
  
  df_y <- df_ratio_wide3 %>%
    dplyr::filter(strain != "JU1516") %>%  ### same isotype with JU1581
    dplyr::filter(strain != "JU1400") %>% ### Exclude extreme outlier
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=trait) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(trait, h2, h2_SE)
  
  df_ratio_h2_YA3 <- rbind(df_ratio_h2_YA3, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

View(df_ratio_h2_YA3)

save(df_h2_YA, df_ratio_h2_YA, df_ratio_h2_YA2, df_ratio_h2_YA3, file="Processed_Data/df_heritability.RData")

plot_ratio_h2_YA_histo3 <- df_ratio_h2_YA3 %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_ratio_h2_YA_histo3

ggsave(file="Plots/plot_ratio_YA_h2_histo_outlier_removed_bySD.pdf", plot=plot_ratio_h2_YA_histo3, width = 7.5, height = 5, units = "in")

### heritability heatmap ###

df_ratio_h2_wide3 <- df_ratio_h2_YA3 %>%
  dplyr::select(-h2_SE) %>%
  dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=h2)

rownames(df_ratio_h2_wide3) <- df_ratio_h2_wide3$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_h2_wide3, -trait2)), hc.order = FALSE, tl.srt = 90, 
           legend.title = "Narrow-sense\nheritability\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.4,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Narrow-sense\nheritability\n")

ggsave("Plots/plot_ratio_h2_heatmap_outlier_removed_bySD.pdf", width = 7.5, height = 7.5, units = "in")


### example plot for high heritability to check what's happening ###

N2_relatedness <- A %>%
  as.data.frame %>%
  dplyr::select(N2_relatedness=N2) %>%
  dplyr::add_rownames(var='strain')

df_example_ratio3 <- all_ratio_YA %>%
  dplyr::filter(trait == "ascr#3/icas#3") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(mean=mean(value), sd=sd(value)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(value>=mean-(3*sd) & value <= mean+(3*sd)) %>%
  dplyr::select(-mean, -sd) %>%
  dplyr::left_join(., N2_relatedness, by='strain') %>%
  na.omit()

df_example_ratio3 %>%
  ggplot(.) +
  geom_point() +
  aes(x=N2_relatedness, y=value) +
  theme_bw()

### heatmap of CV (coeffecient of variation) ###

df_ratio_CV <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(mean=mean(value), sd=sd(value), cv=sd/mean) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait,cv)

df_ratio_CV_wide <- df_ratio_CV %>%
  #dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=cv)

rownames(df_ratio_CV_wide) <- df_ratio_CV_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_CV_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
           legend.title = "Narrow-sense\nheritability\n") +
  theme(axis.text.x = element_text(size=11, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=11, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.4,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right") +
  scale_fill_gradient(limits=c(0,10), name = "Coefficient of variation")

ggsave("Plots/plot_ratio_CV_heatmap.pdf", width = 7.5, height = 7.5, units = "in")


