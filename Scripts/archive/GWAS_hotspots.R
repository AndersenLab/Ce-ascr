library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_ascr_data_norm_h2.RData")

df_chr_length <- data.table::fread("Processed_Data/chr_info.tsv") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::rename(CHROM = V1, stop= V2) %>%
  dplyr::mutate(start=0)

## eigen thresholds

# load independent tests result
total_independent_tests <- read.table("Raw/Analysis_Results-20211018_EIGEN/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_tests <- total_independent_tests[[1]]

independent_test_cutoff <- -log10(0.05/independent_tests)

### load mappings

qtl_directory <- "Raw/Analysis_Results-20211018_EIGEN/Mapping/Processed/"

qtl_files <- list.files(qtl_directory)

for(qtl in 1:length(qtl_files)){
  
  temp_qtl_all <- data.table::fread(glue::glue("{qtl_directory}{qtl_files[qtl]}"))
  
  if(!exists("qtl_all")){
    qtl_all <- temp_qtl_all
  } else {
    qtl_all <- dplyr::bind_rows(qtl_all, temp_qtl_all)
  }
}

df_qtl_bonferroni <- qtl_all %>%
  dplyr::mutate(eigen_BF_adjusted = independent_test_cutoff+log10(23), BF_BF_adjusted = BF+log10(23)) %>% ## same as log10(10^(-independent_test_cutoff)/16) 
  dplyr::mutate(sig_eigen_adjusted = ifelse(log10p >= eigen_BF_adjusted, 1, 0), 
                sig_BF_adjusted = ifelse(log10p >= BF_BF_adjusted, 1, 0)) %>%
  dplyr::filter(sig_eigen_adjusted  == 1)

df_qtl_bonferroni_peaks <- df_qtl_bonferroni %>%
  dplyr::distinct(trait, CHROM, startPOS, peakPOS,endPOS, var.exp, interval_size, log10p, eigen_BF_adjusted,BF_BF_adjusted) %>%
  na.omit()

mecr1_locus <- qtl_all %>%
  dplyr::filter(POS==13692928) %>%
  dplyr::distinct(CHROM, POS, log10p, trait)

eigen_BF_adjusted = independent_test_cutoff+log10(23)

qtl_all_a3a5 <- qtl_all %>%
  dplyr::filter(trait %in% c("ascr.3", "ascr.5"))

save(qtl_all_a3a5, df_qtl_bonferroni, df_qtl_bonferroni_peaks, eigen_BF_adjusted, file="Processed_Data/df_qtl_peaks.RData")

heri_features <- gsub("\\#", "\\_", unique(df_ascr_YA_frac_summary_high_h2$feature))
  
df_qtl_bonferroni_peaks %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::filter(trait %in% heri_features) %>%
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

ggsave("Plots/plot_mapping_summary_individuals.pdf", width = 7.5, height = 4, units = "in")




### manplots ###

for (i in unique(df_gt_qtl_fraction$trait)) {
  
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

### group fraction mapping ###

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
  dplyr::mutate(eigen_BF_adjusted = independent_test_cutoff+log10(7), BF_BF_adjusted = BF+log10(7)) %>% ## same as log10(10^(-independent_test_cutoff)/46) 
  dplyr::mutate(sig_eigen_adjusted = ifelse(log10p >= eigen_BF_adjusted, 1, 0), 
                sig_BF_adjusted = ifelse(log10p >= BF_BF_adjusted, 1, 0)) %>%
  dplyr::filter(sig_eigen_adjusted  == 1)

df_qtl_bonferroni_group_peaks <- df_qtl_bonferroni_group %>%
  dplyr::distinct(trait, CHROM, startPOS, peakPOS,endPOS, var.exp, interval_size, log10p, eigen_BF_adjusted,BF_BF_adjusted) %>%
  na.omit()

df_qtl_bonferroni_group_peaks %>%
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

ggsave("Plots/plot_mapping_summary_groups.pdf", width = 7.5, height = 4, units = "in")

### manplots ###

geno_matrix_frac <- read.table(file = "Processed_Data/Analysis_Results-20201029_ascr_fraction/Genotype_Matrix/Genotype_Matrix.tsv", header = T)

df_gt_qtl_fraction <- geno_matrix_frac %>%
  dplyr::filter(CHROM %in% df_qtl_bonferroni_peaks$CHROM, POS %in% df_qtl_bonferroni_peaks$peakPOS) %>%
  dplyr::select(-REF, -ALT) %>%
  tidyr::gather(-CHROM, -POS, key="strain", value="genotype") %>%
  dplyr::mutate(peakPOS=POS) %>%
  dplyr::left_join(., dplyr::select(df_qtl_bonferroni, peakPOS, trait), by='peakPOS') %>%
  dplyr::select(-POS) %>%
  dplyr::left_join(., df_GWA_fraction_YA, by='strain') %>%
  tidyr::gather(-CHROM, -peakPOS, -strain, -genotype, -trait, key='feature', value = "value") %>%
  dplyr::filter(trait == feature) %>%
  dplyr::distinct()

for (i in unique(df_gt_qtl_fraction$trait)) {
  
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




### Fine mapping ###

gene_directory <- "Processed_Data/Analysis_Results-20200928/Fine_Mappings/Data/"

gene_files <- list.files(gene_directory, pattern = "*_genes.tsv")

for(gene in 1:length(gene_files)){
  temp_gene_all <- data.table::fread(glue::glue("{gene_directory}{gene_files[gene]}"))
  
  if(!exists("gene_all")){
    gene_all <- temp_gene_all
  } else {
    gene_all <- dplyr::bind_rows(gene_all, temp_gene_all)
  }
}

gene_all_summary <- gene_all %>%
  distinct(GENE_NAME, WBGeneID, PEAK_MARKER, TRAIT, MARKER, VARIANT_LOG10p) %>%
  na.omit()

### acox ###

gene_all %>%
  dplyr::filter(MARKER %in% c("I_12940048","IV_6246979"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_acox_pxg.pdf", width = 6, height = 4, units = "in")

### top genes ###

df_top_genes <- gene_all %>%
  dplyr::filter(VARIANT_LOG10p > 1) %>%
  dplyr::distinct(GENE_NAME, WBGeneID, TRAIT, PEAK_MARKER,VARIANT_LOG10p, MARKER) %>%
  na.omit() %>%
  dplyr::group_by(WBGeneID) %>%
  dplyr::mutate(n_TRIAT = n_distinct(TRAIT), max_P = max(as.numeric(VARIANT_LOG10p))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(max_P == VARIANT_LOG10p) %>%
  dplyr::distinct(GENE_NAME, WBGeneID, TRAIT, PEAK_MARKER, n_TRIAT, max_P, .keep_all = T)

gene_all %>%
  dplyr::filter(MARKER %in% c("X_11039760"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_ldb1_pxg.pdf", width = 7.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("X_2711739"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_unc2_pxg.pdf", width = 7.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("X_10739667"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_egl36_pxg.pdf", width = 7.5, height = 3.5, units = "in")


gene_all %>%
  dplyr::filter(MARKER %in% c("II_13410182"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_efl2_pxg.pdf", width = 5.5, height = 3.5, units = "in")


gene_all %>%
  dplyr::filter(MARKER %in% c("X_11077830"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_W04G310_pxg.pdf", width = 7.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("IV_2179042"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_Y71G10AL1_pxg.pdf", width = 5.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("V_15939281"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_srw55_pxg.pdf", width = 5.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("X_2017525"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_y102a11a9_pxg.pdf", width = 5.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("IV_13782643"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_lact9_pxg.pdf", width = 5.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("IV_2155200"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_gcy37_pxg.pdf", width = 5.5, height = 3.5, units = "in")

gene_all %>%
  dplyr::filter(MARKER %in% c("IV_13641431"), !is.na(STRAIN_GENOTYPE)) %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.2, alpha = 0.7) +
  aes(x=STRAIN_GENOTYPE, y=Phenotype_Value) +
  theme_bw() +
  theme(axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        panel.grid=element_blank()) +
  facet_wrap(GENE_NAME~TRAIT, scale= 'free') +
  labs(x="Marker genotype", y="Relative quantity")

ggsave("Plots/plot_ugt22_pxg.pdf", width = 3.5, height = 3.5, units = "in")
