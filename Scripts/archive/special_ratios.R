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
load("Processed_Data/df_qtl_peaks_groups.RData")

df_special_ratios <- df_GWA_fraction_YA %>%
  dplyr::mutate(ascr.10_osas.9.10 = ascr.10/(osas.9+osas.10), osas.9.10_ascr.10 = (osas.9+osas.10)/ascr.10,
                ascr.10_icas.9.10 = ascr.10/(icas.9+icas.10), icas.9.10_ascr.10 = (icas.9+icas.10)/ascr.10,
                ascr.7_ascr.81.8=ascr.7/(ascr.81+ascr.8), ascr.81.8_ascr.7 = (ascr.81+ascr.8)/ascr.7,
                oscrs_ascrs = (ascr.5+oscr.9+oscr.1+oscr.10+oscr.18+bhos.10+bhos.18+bhos.22)/(ascr.9+ascr.1+ascr.7+ascr.10+ascr.3 +ascr.18+bhas.10+bhas.18+bhas.22),
                ascrs_oscrs = (ascr.9+ascr.1+ascr.7+ascr.10+ascr.3+ascr.18+bhas.10+bhas.18+bhas.22)/(ascr.5+oscr.9+oscr.1+oscr.10+oscr.18+bhos.10+bhos.18+bhos.22)) %>%
  dplyr::select(strain, ascr.10_osas.9.10, osas.9.10_ascr.10, ascr.10_icas.9.10, icas.9.10_ascr.10, ascr.7_ascr.81.8, ascr.81.8_ascr.7, oscrs_ascrs, ascrs_oscrs)

write_tsv(df_special_ratios, path="Processed_Data/df_special_ratios.tsv", col_names = T)

# geno matrix

geno_matrix_sp <- read.table(file = "Processed_Data/20201104_DL1/Analysis_Results-20201104/Genotype_Matrix/Genotype_Matrix.tsv", header = T)

geno_matrix_sp %>% dplyr::filter(POS == 2066787)

# load independent tests result
total_independent_tests_sp <- read.table("Processed_Data/20201104_DL1/Analysis_Results-20201104/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_tests_sp <- total_independent_tests_sp[[1]]

independent_test_sp_cutoff <- -log10(0.05/independent_tests_sp)

### load mappings

qtl_directory_sp <- "Processed_Data/20201104_DL1/Analysis_Results-20201104/Mappings/Data/"

qtl_files_sp <- list.files(qtl_directory_sp)[grepl("*mapping.tsv", list.files(qtl_directory_sp))]

qtl_all_sp <- NULL

for(qtl in 1:length(qtl_files_sp)){
  temp_qtl_all <- data.table::fread(glue::glue("{qtl_directory_sp}{qtl_files_sp[qtl]}"))
  
  if(!exists("qtl_all_sp")){
    qtl_all_sp <- temp_qtl_all
  } else {
    qtl_all_sp <- dplyr::bind_rows(qtl_all_sp, temp_qtl_all)
  }
}

df_qtl_sp_bonferroni <- qtl_all_sp %>%
  dplyr::mutate(eigen_BF_adjusted = independent_test_sp_cutoff+log10(8), BF_BF_adjusted = BF+log10(8)) %>% ## same as log10(10^(-independent_test_cutoff)/8) 
  dplyr::mutate(sig_eigen_adjusted = ifelse(log10p >= eigen_BF_adjusted, 1, 0), 
                sig_BF_adjusted = ifelse(log10p >= BF_BF_adjusted, 1, 0)) %>%
  dplyr::filter(sig_eigen_adjusted  == 1)

df_qtl_sp_bonferroni_peaks <- df_qtl_sp_bonferroni %>%
  dplyr::distinct(trait, CHROM, POS, startPOS, peakPOS,endPOS, var.exp, interval_size, log10p, eigen_BF_adjusted,BF_BF_adjusted) %>%
  na.omit()

df_qtl_sp_bonferroni_peaks %>%
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

ggsave("Plots/plot_spratio_summary.pdf", width = 7.5, height = 2, units = "in")

plot_spratio_summary <- df_qtl_sp_bonferroni_peaks %>%
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
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank(),
        strip.text.x =element_blank()) +
  facet_grid(trait~CHROM, scales = 'free', space = 'free', switch = 'y') +
  xlab("Genomic position (Mb)")

plot_spratio_summary

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
        axis.title = element_blank()) +
  facet_grid(trait~CHROM, scales = 'free', space = 'free', switch = 'y') +
  xlab("Genomic position (Mb)")

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
        axis.title = element_text(size=11, color='black'),
        strip.text.x = element_blank()) +
  facet_grid(trait~CHROM, scales = 'free', space = 'free', switch = 'y') +
  xlab("Genomic position (Mb)")

plot_groupspfrac <- cowplot::plot_grid(plot_group_summary, plot_spratio_summary, plot_frac_summary, 
                                     ncol = 1, rel_heights = c(1,1,2.3), align ='v')

ggsave(plot_groupspfrac, file="Plots/plot_mapping_summary_groupspfrac.pdf", width = 7.5, height = 6, units = "in")



### oscrs vs ascrs ###

plot_oscr_ascr_man <- qtl_all_sp %>%
  dplyr::filter(trait == "oscrs_ascrs", CHROM != "MtDNA") %>%
  ggplot(.) +
  geom_point(size = 0.5, alpha = 0.6) +
  aes(x=POS/1e6, y=log10p) +
  geom_hline(yintercept = independent_test_sp_cutoff+log10(8), color = 'red', size = 0.4, linetype = 2) +
  geom_hline(yintercept = independent_test_sp_cutoff, color = 'blue', size = 0.4, linetype = 2) +
  theme_bw() +
  facet_grid(~CHROM, scale='free_x', space='free_x') +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.1, 'line')) +
  labs(x='Genomic position (Mb)', y='-log10(P)') +
  ggtitle("oscrs vs ascrs")

plot_oscr_ascr_man

ggsave(plot_oscr_ascr_man, file = "Plots/plot_oscr_ascr_man.png", width=7.5, height = 3)

df_oscr_ascr_mapping <- qtl_all_sp %>%
  dplyr::filter(trait == "oscrs_ascrs", CHROM != "MtDNA") 

plot_pxg_oscrascr <- df_oscr_ascr_mapping %>%
  na.omit() %>%
  ggplot(.) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  aes(x=factor(allele, levels = c(-1,1), labels = c("REF", "ALT")), y=value) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank()) +
  facet_grid(~marker) +
  labs(y="oscr#s : ascr#s")

plot_pxg_oscrascr

ggsave(plot_pxg_oscrascr, file = "Plots/plot_pxg_oscrascr.pdf", width=6, height = 3)

### Fine_mapping ###

df_oscr_ascr_fine <- read.table(file = "Processed_Data/20201104_DL1/Analysis_Results-20201104/Fine_Mappings/Data/oscrs_ascrs_snpeff_genes.tsv", header = T)

## Chr II QTL ##

df_oscr_ascr_fine_II_summary <- df_oscr_ascr_fine %>%
  dplyr::filter(CHROM=="II") %>%
  dplyr::distinct(MARKER, .keep_all=T) %>%
  dplyr::arrange(-VARIANT_LOG10p)

df_oscr_ascr_fine_II_summary %>%
  dplyr::filter(VARIANT_LOG10p > 5) %>%
  dplyr::distinct(WBFeature_ID)

df_mecr1 <- df_oscr_ascr_fine %>%
  dplyr::filter(AMINO_ACID_CHANGE == "p.Gly159Val")

df_mecr1 %>%
  ggplot(.) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  aes(x=factor(STRAIN_GENOTYPE, levels = c("C","A"), labels = c("REF", "ALT")), y=Phenotype_Value) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank()) +
  facet_grid(~AMINO_ACID_CHANGE) +
  labs(y="oscr#s : ascr#s")

ggsave(file = "Plots/plot_pxg_mecr1.pdf", width=2.5, height = 4)

## Chr IV QTL ##

df_oscr_ascr_fine_IV_summary <- df_oscr_ascr_fine %>%
  dplyr::filter(CHROM=="IV") %>%
  dplyr::distinct(MARKER, .keep_all=T) %>%
  dplyr::arrange(-VARIANT_LOG10p)

df_oscr_ascr_fine_IV_summary %>%
  dplyr::filter(VARIANT_LOG10p > 4) %>%
  dplyr::distinct(WBFeature_ID)

df_oscr_ascr_fine %>%
  dplyr::filter(GENE_NAME == "WBGene00021882", VARIANT_LOG10p > 3) %>%
  ggplot(.) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank()) +
  facet_grid(GENE_NAME~MARKER, scales='free') +
  labs(y="oscr#s : ascr#s")

ggsave(file = "Plots/plot_pxg_nprt1.pdf", width=5, height = 4)


## Chr X QTL ##

df_oscr_ascr_fine_X_summary <- df_oscr_ascr_fine %>%
  dplyr::filter(CHROM=="X") %>%
  dplyr::distinct(MARKER, .keep_all=T) %>%
  dplyr::arrange(-VARIANT_LOG10p)

df_oscr_ascr_fine_X_summary %>%
  dplyr::filter(VARIANT_LOG10p > 4) %>%
  dplyr::distinct(WBFeature_ID)

df_oscr_ascr_fine %>%
  dplyr::filter(GENE_NAME == "WBGene00018657", VARIANT_LOG10p > 2) %>%
  ggplot(.) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank()) +
  facet_grid(~MARKER, scales = 'free') +
  labs(y="oscr#s : ascr#s")

ggsave(file = "Plots/plot_pxg_acl4.pdf", width=6, height = 4)


### ascr#10:(icas#9+icas#10) ###

df_ascr10_icas910_fine <- read.table(file = "Processed_Data/20201104_DL1/Analysis_Results-20201104/Fine_Mappings/Data/ascr.10_icas.9.10_snpeff_genes.tsv", header = T)

df_ascr10_icas910_fine_X_summary <- df_ascr10_icas910_fine %>%
  dplyr::filter(CHROM=="X") %>%
  dplyr::distinct(MARKER, .keep_all=T) %>%
  dplyr::arrange(-VARIANT_LOG10p)

df_ascr10_icas910_fine_X_summary %>%
  dplyr::filter(VARIANT_LOG10p > 3) %>%
  dplyr::distinct(WBFeature_ID)

df_ascr10_icas910_fine %>%
  dplyr::filter(GENE_NAME == "WBGene00008451", VARIANT_LOG10p > 2) %>%
  ggplot(.) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank()) +
  facet_grid(~MARKER, scales = 'free') +
  labs(y="ascr#10 : icas#9+icas#10")

ggsave(file = "Plots/plot_pxg_cest3.pdf", width=5, height = 3)


View(df_ascr10_icas910_fine %>% dplyr::distinct(STRAIN, Phenotype_Value))



df_mecr1_regress <- df_mecr1 %>%
  dplyr::select(STRAIN, STRAIN_GENOTYPE, Phenotype_Value) %>%
  dplyr::mutate(resid_mecr1 = residuals(lm(Phenotype_Value~STRAIN_GENOTYPE)))
  
df_mecr1_regress %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=STRAIN_GENOTYPE, y=resid_mecr1)

df_acl4 <- df_oscr_ascr_fine %>%
  dplyr::filter(GENE_NAME == "WBGene00018657")

df_acl4 %>%
  dplyr::left_join(., dplyr::select(df_mecr1_regress, -STRAIN_GENOTYPE), by = c('STRAIN', 'Phenotype_Value')) %>%
  dplyr::filter(VARIANT_LOG10p > 2) %>%
  ggplot(.) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  aes(x=STRAIN_GENOTYPE, y=resid_mecr1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank()) +
  facet_grid(~MARKER, scales = 'free') +
  labs(y="oscr#s : ascr#s (mecr-1 regression)")
