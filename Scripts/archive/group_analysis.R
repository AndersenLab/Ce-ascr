library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggcorrplot)
library(picante)
library(qgraph)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_GWA.RData")
load("Processed_Data/df_qtl_peaks.RData")

df_chr_length <- data.table::fread("Processed_Data/chr_info.tsv") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::rename(CHROM = V1, stop= V2) %>%
  dplyr::mutate(start=0)

### correlation analysis for each group ###

ascr_YA_frac_grouped <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::left_join(., df_ascr_group, by = 'feature') %>%
  dplyr::group_by(strain, group1) %>%
  dplyr::summarise(sum_group_fraction = sum(ascr_fraction)) %>%
  dplyr::ungroup()

df_YA_group <- ascr_YA_frac_grouped %>%
  dplyr::select(strain, group1, sum_group_fraction) %>%
  tidyr::spread(key=strain, value=sum_group_fraction)

df_YA_group_wide <- t(df_YA_group[2:ncol(df_YA_group)])
colnames(df_YA_group_wide) <- df_YA_group[[1]]

test_result_YA_group <- cor.table(df_YA_group_wide, cor.method="spearman")[[1]]

ggcorrplot(test_result_YA_group, hc.order = FALSE, tl.srt = 90, legend.title = "Correlation\n(Rho)\n") +
  theme(axis.text.x = element_text(size=9, color='black', vjust=0.6, hjust = 1),
        axis.text.y = element_text(size=9, color='black'),
        legend.position = 'bottom',
        legend.key.width = unit(0.2,"in"),
        legend.key.height = unit(0.1,"in"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9, color='black'),
        panel.grid=element_blank()) +
  scale_y_discrete(position = "right")

ggsave("Plots/plot_corr_group.pdf", width = 3, height = 4, units = "in")

plot_ascrs_oscrs <- df_YA_group_wide %>%
  as.data.frame() %>%
  dplyr::rename(ascrs='ascr(s)', oscrs='oscr(s)') %>%
  ggplot(.) +
  geom_point(alpha=0.6) +
  aes(x=ascrs, y=oscrs) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="Relative quantity of short-chain ascr#s", y="Relative quantity of short-chain oscr#s")

plot_ascrs_oscrs

ggsave(plot_ascrs_oscrs, file="Plots/plot_ascrs_oscrs.pdf", width = 4, height = 4, units = "in")

plot_ascrl_oscrl <- df_YA_group_wide %>%
  as.data.frame() %>%
  dplyr::rename(ascrl='ascr(l)', oscrl='oscr(l)') %>%
  ggplot(.) +
  geom_point(alpha=0.6) +
  aes(x=ascrl, y=oscrl) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'))+
  labs(x="Relative quantity of long-chain ascr", y="Relative quantity of long-chain oscr")

plot_ascrl_oscrl

#ggsave(plot_ascrl_oscrl, file="Plots/plot_ascrl_oscrl.pdf", width = 4, height = 4, units = "in")

plot_ascr3_ascr5 <- df_YA_pheno_wide %>%
  as.data.frame() %>%
  dplyr::rename(ascr3='ascr#3', ascr5='ascr#5') %>%
  ggplot(.) +
  geom_point(alpha=0.6) +
  aes(x=ascr3, y=ascr5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'))+
  labs(x="Relative quantity of ascr#3", y="Relative quantity of ascr#5")

plot_ascr3_ascr5

#ggsave(plot_ascr3_ascr5, file="Plots/plot_ascr3_ascr5.pdf", width = 4, height = 4, units = "in")

cor(df_YA_pheno_wide[,2], df_YA_pheno_wide[,3], method='spearman')

ascr_YA_frac_grouped_woascr35 <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400", !feature %in% c("ascr#3", "ascr#5")) %>%
  dplyr::filter(!feature %in% c("ascr#3", "ascr#5")) %>%
  dplyr::left_join(., df_ascr_group, by = 'feature') %>%
  dplyr::group_by(strain, group1) %>%
  dplyr::summarise(sum_group_fraction = sum(ascr_fraction)) %>%
  dplyr::ungroup()

df_YA_group_woascr35 <- ascr_YA_frac_grouped_woascr35 %>%
  dplyr::select(strain, group1, sum_group_fraction) %>%
  tidyr::spread(key=strain, value=sum_group_fraction)

df_YA_group_wide_woascr35 <- t(df_YA_group_woascr35[2:ncol(df_YA_group_woascr35)])
colnames(df_YA_group_wide_woascr35) <- df_YA_group_woascr35[[1]]

plot_ascrs_oscrs_woascr35 <- df_YA_group_wide_woascr35 %>%
  as.data.frame() %>%
  dplyr::rename(ascrs='ascr(s)', oscrs='oscr(s)') %>%
  ggplot(.) +
  geom_point(alpha=0.6) +
  aes(x=ascrs, y=oscrs) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="Relative quantity of short-chain ascr (ex ascr#3)", y="Relative quantity of short-chain oscr (ex ascr#5)")

plot_ascrs_oscrs_woascr35

#ggsave(plot_ascrs_oscrs_woascr35, file="Plots/plot_ascrs_oscrs_woascr35.pdf", width = 4, height = 4, units = "in")

cor(df_YA_group_wide_woascr35[,1], df_YA_group_wide_woascr35[,7], method='spearman')

### data for GWA ###

df_GWA_group_fraction_YA <- data.frame(strain = rownames(df_YA_group_wide), df_YA_group_wide) %>%
  dplyr::filter(!strain %in% c("JU1242","JU1400")) %>%
  dplyr::mutate(strain = ifelse(strain=="PB306","ECA259", strain))

write_tsv(df_GWA_group_fraction_YA, path="Processed_Data/ascr_group_fraction_YA_OutlierRemoval.tsv", col_names = T)

### mapping results ###

### group fraction mapping ###

# load independent tests result
total_independent_tests_group <- read.table("Processed_Data/Analysis_Results-20201029_ascr_group/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_tests_group <- total_independent_tests_group[[1]]

independent_test_cutoff_group <- -log10(0.05/independent_tests_group)

### load mappings

qtl_directory_group <- "Processed_Data/Analysis_Results-20201029_ascr_group/Mappings/Data/"

qtl_files_group <- list.files(qtl_directory_group)

for(qtl in 1:length(qtl_files_group)){
  
  temp_qtl_all <- data.table::fread(glue::glue("{qtl_directory_group}{qtl_files_group[qtl]}"))
  
  if(!exists("qtl_all_group")){
    qtl_all_group <- temp_qtl_all
  } else {
    qtl_all_group <- dplyr::bind_rows(qtl_all_group, temp_qtl_all)
  }
}

df_qtl_bonferroni_group <- qtl_all_group %>%
  dplyr::mutate(eigen_BF_adjusted = independent_test_cutoff_group+log10(9), BF_BF_adjusted = BF+log10(9)) %>% ## same as log10(10^(-independent_test_cutoff)/9) 
  dplyr::mutate(sig_eigen_adjusted = ifelse(log10p >= eigen_BF_adjusted, 1, 0), 
                sig_BF_adjusted = ifelse(log10p >= BF_BF_adjusted, 1, 0)) %>%
  dplyr::filter(sig_eigen_adjusted  == 1)

df_qtl_bonferroni_group_peaks <- df_qtl_bonferroni_group %>%
  dplyr::distinct(trait, CHROM, startPOS, peakPOS,endPOS, var.exp, interval_size, log10p, eigen_BF_adjusted,BF_BF_adjusted) %>%
  na.omit()

save(qtl_all_group, df_qtl_bonferroni_group_peaks, file="Processed_Data/df_qtl_peaks_groups.RData")

plot_group_summary <- df_qtl_bonferroni_group_peaks %>%
  #dplyr::filter(CHROM != "MtDNA") %>%
  ggplot(.) +
  geom_rect(aes(xmin = peakPOS/1e6,    # this is the plot boundary for LD and gene plots
                xmax = peakPOS/1e6,    # this is the plot boundary for LD and gene plots
                ymin = 0, 
                ymax = Inf), color = 'red') +
  geom_rect(aes(xmin = startPOS/1e6,    
                xmax = endPOS/1e6,    
                ymin = 0, 
                ymax = Inf), alpha=0.2) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 0.2, ymax=0.3), 
            color='transparent', fill='transparent', size =0.1) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "lines"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y=element_text(size=9, color = "black", angle = 90),  
        strip.text.y.left = element_text(angle=0),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size=11, color='black')) +
  facet_grid(trait~CHROM, scales = 'free', space = 'free', switch = 'y') +
  xlab("Genomic position (Mb)")

ggsave(plot_group_summary, file="Plots/plot_mapping_summary_groups.pdf", width = 7.5, height = 2, units = "in")

plot_frac_summary <- df_qtl_bonferroni_peaks %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  ggplot(.) +
  geom_rect(aes(xmin = peakPOS/1e6,    # this is the plot boundary for LD and gene plots
                xmax = peakPOS/1e6,    # this is the plot boundary for LD and gene plots
                ymin = 0, 
                ymax = Inf), color = 'red') +
  geom_rect(aes(xmin = startPOS/1e6,    
                xmax = endPOS/1e6,    
                ymin = 0, 
                ymax = Inf), alpha=0.2) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 0.2, ymax=0.3), 
            color='transparent', fill='transparent', size =0.1) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "lines"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y=element_text(size=9, color = "black", angle = 90),  
        strip.text.y.left = element_text(angle=0),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size=11, color='black')) +
  facet_grid(trait~CHROM, scales = 'free', space = 'free', switch = 'y') +
  xlab("Genomic position (Mb)")

plot_fracgroup <- cowplot::plot_grid(plot_group_summary, plot_frac_summary, 
                                     ncol = 1, rel_heights = c(1,2.2), align ='v')

ggsave(plot_fracgroup, file="Plots/plot_mapping_summary_fracgroups.pdf", width = 7.5, height = 5.5, units = "in")


### pxg ###

geno_matrix_group <- read.table(file = "Processed_Data/Analysis_Results-20201029_ascr_group/Genotype_Matrix/Genotype_Matrix.tsv", header = T)

for (i in unique(df_qtl_bonferroni_group_peaks$trait) {
  
  
}

df_gt_qtl_group <- geno_matrix_group %>%
  dplyr::filter(CHROM %in% df_qtl_bonferroni_group_peaks$CHROM, POS %in% df_qtl_bonferroni_group_peaks$peakPOS) %>%
  dplyr::select(-REF, -ALT) %>%
  tidyr::gather(-CHROM, -POS, key="strain", value="genotype") %>%
  dplyr::mutate(peakPOS=POS) %>%
  dplyr::left_join(., dplyr::select(df_qtl_bonferroni_group, peakPOS, trait), by='peakPOS') %>%
  dplyr::select(-POS) %>%
  dplyr::left_join(., df_GWA_fraction_YA, by='strain') %>%
  tidyr::gather(-CHROM, -peakPOS, -strain, -genotype, -trait, key='feature', value = "value") %>%
  dplyr::filter(trait == feature) %>%
  dplyr::distinct()

for (i in unique(df_gt_qtl_group$trait)) {
  
  df_gt_qtl_fraction %>%
    dplyr::filter(trait == i) %>%
    dplyr::mutate(marker = paste(CHROM, peakPOS, sep=":"), genotype = ifelse(genotype==-1, "REF", "ALT")) %>%
    ggplot(.) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
    aes(x=factor(genotype, levels = c("REF", "ALT")), y=value) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=10, color='black'),
          axis.title.y = element_text(size=11, color='black'),
          axis.title.x = element_blank()) +
    facet_grid(trait~marker)
  
  ggsave(file = glue::glue("Plots/pxg/fractions/plot_pxg_fraction_{i}.png"), width=7.5, height = 5)
  
}


