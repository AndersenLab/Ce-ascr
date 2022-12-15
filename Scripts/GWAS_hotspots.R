library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggnewscale)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_ascr_data_norm_h2.RData")

df_chr_length <- data.table::fread("Processed_Data/chr_info.tsv") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::rename(CHROM = V1, stop= V2) %>%
  dplyr::mutate(start=0)

## eigen thresholds

# load independent tests result
total_independent_tests <- read.table("Raw/Analysis_Results-20211021_ratio_outliers/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_tests <- total_independent_tests[[1]]

independent_test_cutoff <- -log10(0.05/independent_tests)

eigen_BF_adjusted = independent_test_cutoff+log10(24)

### load mappings

### ascr#3:ascr##5 traits

qtl_a3a5 <- data.table::fread("Raw/Analysis_Results-20211021_ratio_outliers/Mapping/Processed/processed_a3a5_AGGREGATE_mapping.tsv")

plot_man_a3a5 <- qtl_a3a5 %>%
  dplyr::distinct(CHROM, POS, trait, .keep_all=T) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(QTL = ifelse(log10p >=eigen_BF_adjusted, "yes", "no")) %>%
  ggplot(.) +
  geom_rect(data=dplyr::filter(qtl_a3a5,CHROM != "MtDNA", log10p >=eigen_BF_adjusted),
            aes(xmin = startPOS/1e6,    # this is the plot boundary for LD and gene plots
                xmax = endPOS/1e6,    # this is the plot boundary for LD and gene plots
                ymin = 0, 
                ymax = Inf, 
                fill = "blue"), 
            color = "blue",fill = "cyan",linetype = 2, 
            alpha=.2) +  
  geom_point(alpha=0.7, aes(x=POS/1e6, y=log10p, fill=QTL), size = 1, shape=21) +
  geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 7, ymax=8), 
            color='transparent', fill='transparent', size =0.1) +
  geom_hline(yintercept = eigen_BF_adjusted, color='darkgrey', alpha = 0.8, linetype=2) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.spacing = unit(0.1, "lines"),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        axis.text.x = element_blank(),
        title = element_text(size=11, color='black'),
        legend.position = 'none') +
  facet_grid(~CHROM, scale='free_x', space = 'free') +
  labs(x="Genomic position (Mb)", y="-log10(p)") +
  ggtitle("ascr#3:ascr#5") +
  scale_y_continuous(breaks = c(0,3,6,9,12)) +
  scale_fill_manual(values=c("grey","red"))

plot_man_a3a5

strains.of.interest <- c("N2", "ED3052", "NIC256", "JU258", "LKC34","NIC166")

plot_pxg_a3a5 <- qtl_a3a5 %>%
  dplyr::filter(!is.na(peak_id)) %>%
  dplyr::select(CHROM, marker, trait, startPOS, peakPOS, endPOS, AF1, value, strain, allele, peak_id) %>%
  dplyr::distinct() %>%
  dplyr::mutate(startPOS = startPOS/1000000,
                peakPOS = peakPOS/1000000,
                endPOS = endPOS/1000000) %>%
  dplyr::mutate(allele = dplyr::case_when(allele == "-1" ~ "REF",
                                          allele == "1" ~ "ALT",
                                          TRUE ~ "NA"),
                allele = factor(allele, levels = c("REF", "ALT"))) %>%
  dplyr::filter(allele != "NA" | !is.na(allele)) %>%
  dplyr::mutate(SOI = strain %in% strains.of.interest,
                SOI.2 = if_else(SOI == TRUE, true = strain, false = "")) %>%
  droplevels() %>%
  dplyr::arrange(SOI.2) %>%
  ggplot2::ggplot(mapping = aes(x = allele, y = value, text = SOI.2)) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::geom_violin(aes(fill = allele), alpha = 0.8, scale = "count", draw_quantiles = c(0.25, 0.5, 0.75)) +
  ggplot2::scale_fill_manual(values = c("#726E75","#720E07"), guide = FALSE) +
  ggnewscale::new_scale("fill") +
  ggplot2::geom_point(aes(fill = SOI), position = ggbeeswarm::position_beeswarm(), size = 1.5, shape = 21) +
  # geom_point(aes(colour = sweep.share*100), size = 1.1, position = pos) +
  ggplot2::scale_fill_manual(values = c("#9297C4","#D33E43"), guide = FALSE) +
  # scale_colour_gradient(low = "black", high = "violetred", name = "Selective Sweep (% Chromosome)") +
  ggrepel::geom_text_repel(aes(label = SOI.2),
                           colour = "black", position = ggbeeswarm::position_beeswarm()) +
  ggplot2::theme(legend.position = "bottom",
                 axis.text = element_text(size=10, color='black'),
                 axis.title.y = element_text(size=11, color='black'),
                 axis.title.x = element_blank(),
                 title = element_text(size=11, color='black')) +
  ggplot2::labs(y = "Trait Value",
                x = "Genotype") +
  ggplot2::facet_grid(~marker)

plot_pxg_a3a5

### all fraction traits ###

qtl_directory <- "Raw/Analysis_Results-20211026/Mapping/Processed/aggregate_mapping/"

qtl_files <- list.files(qtl_directory)[grepl("AGGREGATE_mapping.tsv", list.files(qtl_directory))]

qtl_all <- NULL

for(qtl in 1:length(qtl_files)){
  
  temp_qtl_all <- data.table::fread(glue::glue("{qtl_directory}{qtl_files[qtl]}"))
  
  if(!exists("qtl_all")){
    qtl_all <- temp_qtl_all
  } else {
    qtl_all <- dplyr::bind_rows(qtl_all, temp_qtl_all)
  }
}

df_qtl_bonferroni <- rbind(qtl_all,qtl_a3a5) %>%
  dplyr::mutate(eigen_BF_adjusted = independent_test_cutoff+log10(24), BF_BF_adjusted = BF+log10(24)) %>% ## same as log10(10^(-independent_test_cutoff)/23) 
  dplyr::mutate(sig_eigen_adjusted = ifelse(log10p >= eigen_BF_adjusted, 1, 0), 
                sig_BF_adjusted = ifelse(log10p >= BF_BF_adjusted, 1, 0)) %>%
  dplyr::filter(sig_eigen_adjusted  == 1)

df_qtl_bonferroni_peaks <- df_qtl_bonferroni %>%
  dplyr::distinct(trait, CHROM, startPOS, peakPOS,endPOS, var.exp, interval_size, log10p, eigen_BF_adjusted,BF_BF_adjusted) %>%
  na.omit()

mecr1_locus <- qtl_all %>%
  dplyr::filter(POS==13692928) %>%
  dplyr::distinct(CHROM, POS, log10p, trait)

eigen_BF_adjusted = independent_test_cutoff+log10(24)

qtl_all_a3a5 <- qtl_all %>%
  dplyr::filter(trait %in% c("ascr.3", "ascr.5"))

save(qtl_a3a5, qtl_all_a3a5, df_qtl_bonferroni, df_qtl_bonferroni_peaks, eigen_BF_adjusted, file="Processed_Data/df_qtl_peaks.RData")

heri_features <- unique(df_ascr_YA_frac_summary_high_h2$feature)
df_qtl_bonferroni_peaks$trait <- gsub("\\_", "\\#", df_qtl_bonferroni_peaks$trait)

df_qtl_bonferroni_peaks_h2 <- df_qtl_bonferroni_peaks %>%
  dplyr::filter(trait %in% heri_features | trait == "a3a5")

length(unique(df_qtl_bonferroni_peaks_h2$trait))


  
df_qtl_bonferroni_peaks %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::filter(trait %in% heri_features | trait == "a3a5") %>%
  dplyr::mutate(trait=ifelse(trait =="a3a5", "ascr#3:\nascr#5", trait)) %>%
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

### manplots ###

geno_matrix_frac <- read.table(file = "Raw/Analysis_Results-20211018_EIGEN/Genotype_Matrix/Genotype_Matrix.tsv", header = T)

df_gt_qtl_fraction <- geno_matrix_frac %>%
  dplyr::filter(CHROM %in% df_qtl_bonferroni_peaks$CHROM, POS %in% df_qtl_bonferroni_peaks$peakPOS) %>%
  dplyr::select(-REF, -ALT) %>%
  tidyr::gather(-CHROM, -POS, key="strain", value="genotype") %>%
  dplyr::mutate(peakPOS=POS) %>%
  dplyr::left_join(., dplyr::select(df_qtl_bonferroni, peakPOS, trait), by='peakPOS') %>%
  dplyr::select(-POS) %>%
  dplyr::left_join(., df_GWA_fraction_YA, by='strain') %>%
  tidyr::gather(-CHROM, -peakPOS, -strain, -genotype, -trait, key='feature', value = "value") %>%
  dplyr::filter(trait ==  gsub("\\.", "\\_", feature)) %>%
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
