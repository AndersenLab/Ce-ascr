library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggcorrplot)
library(picante)
library(qgraph)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_GWA.RData")
## load("Processed_Data/df_qtl_peaks.RData")

df_chr_length <- data.table::fread("Processed_Data/chr_info.tsv") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::rename(CHROM = V1, stop= V2) %>%
  dplyr::mutate(start=0)

## ratio traits ##

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
                         paste(as.character(ratio[1,i]),as.character(ratio[2,i]),sep="_"),
                         paste(as.character(ratio[2,i]),as.character(ratio[1,i]),sep="_"))
  as.ratio.long <- tidyr::gather(comb.as,trait,value,-strain)
  ratio_list_YA[[i]] <- as.ratio.long
}

all_ratio_YA <- bind_rows(ratio_list_YA)

df_ratio_wide <- all_ratio_YA %>%
  dplyr::filter(!strain %in% c("JU1400", "JU1242")) %>%
  tidyr::spread(trait,value) %>%
  dplyr::mutate(strain = ifelse(strain=="PB306","ECA259", strain))

save(all_ratio_YA, df_ratio_wide, file = "Processed_Data/ratio_trait.RData")

write_tsv(df_ratio_wide, path="Processed_Data/ascr_ratio_YA_OutlierRemoval.tsv", col_names = T)

### analyze mapping results ###

df_chr_length <- data.table::fread("Processed_Data/chr_info.tsv") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::rename(CHROM = V1, stop= V2) %>%
  dplyr::mutate(start=0)

# load independent tests result
total_independent_tests_ratio <- read.table("Processed_Data/Analysis_Results-20201103_ascr_ratio/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_tests_ratio <- total_independent_tests_ratio[[1]]

independent_test_cutoff_ratio <- -log10(0.05/independent_tests_ratio)

### load ratio mapping ###

qtl_directory_ratio <- "Processed_Data/Analysis_Results-20201103_ascr_ratio/Mappings/Data/"

qtl_files_ratio <- list.files(qtl_directory_ratio)[grepl("*mapping.tsv", list.files(qtl_directory_ratio))]

qtl_all_ratio <- NULL

for(qtl in 1:length(qtl_files_ratio)){
  temp_qtl_all <- data.table::fread(glue::glue("{qtl_directory_ratio}{qtl_files_ratio[qtl]}")) %>%
    na.omit()
  
  if(!exists("qtl_all_ratio")){
    qtl_all_ratio <- temp_qtl_all
  } else {
    qtl_all_ratio <- dplyr::bind_rows(qtl_all_ratio, temp_qtl_all)
  }
}

df_qtl_ratio_bonferroni <- qtl_all_ratio %>%
  dplyr::mutate(eigen_BF_adjusted = independent_test_cutoff_ratio+log10(2070), BF_BF_adjusted = BF+log10(207)) %>% ## same as log10(10^(-independent_test_cutoff)/2070) 
  dplyr::mutate(sig_eigen_adjusted = ifelse(log10p >= eigen_BF_adjusted, 1, 0), 
                sig_BF_adjusted = ifelse(log10p >= BF_BF_adjusted, 1, 0)) %>%
  dplyr::filter(sig_eigen_adjusted  == 1)

df_ratio_mapping_peaks <- df_qtl_ratio_bonferroni %>%
  dplyr::distinct(marker, trait, .keep_all=T)
  
save(qtl_all_ratio, df_qtl_ratio_bonferroni, df_ratio_mapping_peaks, file="Processed_Data/df_ratio_mapping.RData")

## histogram ##

df_ratio_mapping_peaks %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  ggplot(.) +
  geom_histogram(binwidth=1, aes(x=POS/1e6) ) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=2), 
            color='transparent', fill='transparent', size =0.1) +
  theme_bw() +
  facet_grid(~CHROM, scale='free_x', space='free_x') +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.1, 'line'),
        axis.text.x = element_text(size=10, color='black')) +
  labs(x='Genomic position (Mb)', y='-log10(P)') +
  ggtitle("Pair-wise ratio QTL")

ggsave("Plots/plot_pairwise_histo.pdf", width = 7.5, height = 3, units = "in")


qtl_all_ratio_MT <- qtl_all_ratio %>%
  dplyr::filter(CHROM=="MtDNA")

qtl_all_ratio_MT %>%
  ggplot(.) +
  geom_jitter(width=0.2, size = 0.8, alpha=0.8) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.6) +
  #geom_text(data=dplyr::filter(qtl_all_ratio_MT, allele==1),  
                  label=dplyr::filter(qtl_all_ratio_MT, allele==1)$strain,
                  size=2, color='red') +
  aes(x=factor(allele), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank())+
  facet_wrap(marker~trait, scales='free')

ggsave("Plots/plot_MtDNA_pxg.pdf", width = 10, height = 9, units = "in")

df_MT_qtl_geno <- qtl_all_ratio_MT %>%
  dplyr::distinct(marker, allele, strain) %>%
  tidyr::spread(key='marker', value='allele')

df_MT_qtl_geno %>%
  dplyr::mutate(geno_comb = paste(MtDNA_262, MtDNA_355, MtDNA_4202, MtDNA_6766, MtDNA_12503, sep="_")) %>%
  na.omit() %>%
  dplyr::group_by(geno_comb) %>%
  dplyr::summarise(n=n())

write.csv(df_MT_qtl_geno, file="MtDNA_QTL_geno.csv")


df_MT_qtl_geno_pheno <- df_MT_qtl_geno %>%
  dplyr::mutate(geno_comb = paste(MtDNA_262, MtDNA_355, MtDNA_4202, MtDNA_6766, MtDNA_12503, sep="_")) %>%
  na.omit() %>%
  dplyr::left_join(., dplyr::distinct(qtl_all_ratio_MT, strain, trait, value), by='strain')

df_MT_qtl_geno_pheno %>%
  na.omit() %>%
  ggplot(.) +
  geom_jitter(width=0.2, size = 0.8, alpha=0.8) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.6) +
  aes(x=geno_comb, y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  facet_wrap(~trait, scales='free', nrow=4)

ggsave("Plots/plot_MtDNA_pxg_genocomb.pdf", width = 9, height = 10, units = "in")

df_MT_qtl_geno_pheno %>%
  dplyr::filter(geno_comb == "1_1_-1_1_1") %>%
  na.omit() %>%
  ggplot(.) +
  geom_col() +
  aes(x=strain, y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90))+
  facet_wrap(~trait, scales='free', nrow=4)

ggsave("Plots/plot_MtDNA_pxg_genocomb_alts.pdf", width = 9, height = 10, units = "in")

unique(df_MT_qtl_geno_pheno$geno_comb)


df_qtl_summary_adjusted %>%
  dplyr::distinct(CHROM, POS, trait, .keep_all=T) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  ggplot(.) +
  geom_point(alpha=0.7, aes(x=POS/1e6, y=log10p), size = 0.5) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 7, ymax=8), 
            color='transparent', fill='transparent', size =0.1) +
  geom_hline(yintercept = eigen_BF_adjusted, color='red') +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.spacing = unit(0.1, "lines")) +
  facet_grid(~CHROM, scale='free_x', space = 'free') +
  labs(x="Genomic position (Mb)", y="-log10(p)")

ggsave("Plots/plot_ratio_summary_manhattan.pdf", width = 7.5, height = 3, units = "in")

### Mitochondrial mapping

df_qtl_summary_adjusted %>%
  dplyr::distinct(CHROM, POS, trait, .keep_all=T) %>%
  dplyr::filter(CHROM == "MtDNA") %>%
  ggplot(.) +
  geom_point(alpha=0.7, aes(x=POS/1e6, y=log10p), size = 0.5) +
  geom_hline(yintercept = eigen_BF_adjusted, color='red') +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.spacing = unit(0.1, "lines")) +
  facet_grid(~CHROM, scale='free_x', space = 'free') +
  labs(x="Genomic position (Mb)", y="-log10(p)")

ggsave("Plots/plot_ratio_summary_manhattan_MtDNA.pdf", width = 5, height = 2.5, units = "in")

### example qtl hospot heatmap ###

#1. V_14773914

ratioqtl_V_14773914 <- df_ratio_mapping %>%
  dplyr::filter(CHROM == "V", POS == 14773914)

df_V_14773914_wide <- ratioqtl_V_14773914 %>%
  dplyr::distinct(trait, log10p) %>%
  tidyr::separate(trait, sep='_', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=log10p)

rownames(df_V_14773914_wide) <- df_V_14773914_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_V_14773914_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,max(ratioqtl_V_14773914$log10p)+1), name = "log10(p)\n", low='white', high='red') +
  ggtitle("V:14773914")

ggsave("Plots/plot_V_14773914_log10p_heatmap.pdf", width = 7.5, height = 7.5, units = "in")

#2. II_1060797

ratioqtl_II_1060797 <- df_ratio_mapping %>%
  dplyr::filter(CHROM == "II", POS == 1060797)

df_II_1060797_wide <- ratioqtl_II_1060797 %>%
  dplyr::distinct(trait, log10p) %>%
  tidyr::separate(trait, sep='_', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=log10p)

rownames(df_II_1060797_wide) <- df_II_1060797_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_II_1060797_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,max(ratioqtl_II_1060797$log10p)+1), name = "log10(p)\n", low='white', high='red') +
  ggtitle("II_1060797")

ggsave("Plots/plot_II_1060797_log10p_heatmap.pdf", width = 7.5, height = 7.5, units = "in")

#3. II_13655631

ratioqtl_II_13655631 <- df_ratio_mapping %>%
  dplyr::filter(CHROM == "II", POS == 13655631)

df_II_13655631_wide <- ratioqtl_II_13655631 %>%
  dplyr::distinct(trait, log10p) %>%
  tidyr::separate(trait, sep='_', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=log10p)

rownames(df_II_13655631_wide) <- df_II_13655631_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_II_13655631_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,max(ratioqtl_II_13655631$log10p)+1), name = "log10(p)\n", low='white', high='red') +
  ggtitle("II_13655631")

ggsave("Plots/plot_II_13655631_log10p_heatmap.pdf", width = 7.5, height = 7.5, units = "in")


#4. II_2232339

ratioqtl_II_2232339 <- df_ratio_mapping %>%
  dplyr::filter(CHROM == "II", POS == 2232339)

df_II_2232339_wide <- ratioqtl_II_2232339 %>%
  dplyr::distinct(trait, log10p) %>%
  tidyr::separate(trait, sep='_', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=log10p)

rownames(df_II_2232339_wide) <- df_II_2232339_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_II_2232339_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,max(ratioqtl_II_2232339$log10p)+1), name = "log10(p)\n", low='white', high='red') +
  ggtitle("II_2232339")

ggsave("Plots/plot_II_2232339_log10p_heatmap.pdf", width = 7.5, height = 7.5, units = "in")



#5. MtDNA_262

ratioqtl_MtDNA_262 <- df_ratio_mapping %>%
  dplyr::filter(CHROM == "MtDNA", POS == 262)

df_MtDNA_262_wide <- ratioqtl_MtDNA_262 %>%
  dplyr::distinct(trait, log10p) %>%
  tidyr::separate(trait, sep='_', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=log10p)

rownames(df_MtDNA_262_wide) <- df_MtDNA_262_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_MtDNA_262_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,max(ratioqtl_MtDNA_262$log10p)+1), name = "log10(p)\n", low='white', high='red') +
  ggtitle("MtDNA_262")

ggsave("Plots/plot_MtDNA_262_log10p_heatmap.pdf", width = 7.5, height = 7.5, units = "in")


### Ratio traits of Interest ###

ratio_interest <- c("ascr.18_ascr.9", "ascr.9_ascr.18", "ascr.18_ascr.1", "ascr.1_ascr.18", "ascr.18_ascr.10", "ascr.10_ascr.18", "oscr.10_ascr.5", "ascr.5_oscr.10", "icas.10_icas.9", "icas.9_icas.10", "osas.10_osas.9", "osas.9_osas.10", "ascr.22_bhas.22", "bhas.22_ascr.22", "ascr.18_bhas.18", "bhas.18_ascr.18", "oscr.10_bhos.10", "bhos.10_oscr.10", "oscr.18_bhos.18", "bhos.18_oscr.18", "oscr.18_ascr.18", "ascr.18_oscr.18", "ascr.5_ascr.9", "ascr.9_ascr.5", "ascr.10_ascr.3", "ascr.3_ascr.10", "ascr.1_ascr.7", "ascr.7_ascr.1", "ascr.3_mbas.3", "mbas.3_ascr.3", "ascr.3_icas.3", "icas.3_ascr.3", "ascr.7_iglas.7", "iglas.7_ascr.7", "ascr.3_iglas.3", "iglas.3_ascr.3", "ascr.1_uglas.11", "uglas.11_ascr.1", "ascr.1_anglas.1", "anglas.1_ascr.1", "ascr.1_glas.1", "glas.1_ascr.1", "ascr.3_glas.3", "glas.3_ascr.3", "ascr.7_ascr.81", "ascr.81_ascr.7")

df_ratio_mapping_interest <- df_ratio_mapping_peaks %>%
  dplyr::filter(trait %in% ratio_interest) %>%
  dplyr::left_join(., df_qtl_summary, by = c('CHROM', 'POS', 'trait'))
  
eigen_BF_adjusted_ratio = independent_test_cutoff_ratio+log10(46) ## same as log10(10^(-independent_test_cutoff)/44) 

df_qtl_summary_adjusted_interest <- df_ratio_mapping_interest %>%
  dplyr::mutate(eigen_BF_adjusted = eigen_BF_adjusted_ratio) %>% 
  dplyr::mutate(sig_eigen_adjusted = ifelse(log10p >= eigen_BF_adjusted, 1, 0)) %>%
  dplyr::filter(sig_eigen_adjusted  == 1) %>%
  dplyr::group_by(CHROM, POS) %>%
  dplyr::mutate(N_qtl_traits=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(marker=paste(CHROM, POS, sep='_')) %>%
  na.omit()

df_qtl_summary_adjusted_interest$trait <- factor(df_qtl_summary_adjusted_interest$trait, levels = c("ascr.18_ascr.9", "ascr.9_ascr.18", "ascr.18_ascr.1", "ascr.1_ascr.18", "ascr.18_ascr.10", "ascr.10_ascr.18", "oscr.10_ascr.5", "ascr.5_oscr.10", "icas.10_icas.9", "icas.9_icas.10", "osas.10_osas.9", "osas.9_osas.10", "ascr.22_bhas.22", "bhas.22_ascr.22", "ascr.18_bhas.18", "bhas.18_ascr.18", "oscr.10_bhos.10", "bhos.10_oscr.10", "oscr.18_bhos.18", "bhos.18_oscr.18", "oscr.18_ascr.18", "ascr.18_oscr.18", "ascr.5_ascr.9", "ascr.9_ascr.5", "ascr.10_ascr.3", "ascr.3_ascr.10", "ascr.1_ascr.7", "ascr.7_ascr.1", "ascr.3_mbas.3", "mbas.3_ascr.3", "ascr.3_icas.3", "icas.3_ascr.3", "ascr.7_iglas.7", "iglas.7_ascr.7", "ascr.3_iglas.3", "iglas.3_ascr.3", "ascr.1_uglas.11", "uglas.11_ascr.1", "ascr.1_anglas.1", "anglas.1_ascr.1", "ascr.1_glas.1", "glas.1_ascr.1", "ascr.3_glas.3", "glas.3_ascr.3", "ascr.7_ascr.81", "ascr.81_ascr.7"))

df_qtl_bonferroni_peaks %>%
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
        strip.background = element_blank()) +
  facet_grid(trait~CHROM, scales = 'free', space = 'free', switch = 'y') +
  xlab("Genomic position (Mb)")

df_qtl_summary_adjusted_interest %>%
  dplyr::rename(peakPOS=POS) %>%
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
        strip.background = element_blank()) +
  facet_grid(trait~CHROM, scales = 'free', space = 'free', switch = 'y') +
  xlab("Genomic position (Mb)")

ggsave("Plots/plot_mapping_summary_ratios_interest.pdf", width = 7.5, height = 6, units = "in")


