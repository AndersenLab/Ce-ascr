library(genetics)

load("Processed_Data/df_qtl_peaks.RData")

## Linkage disequilibrium (LD) between QTL

processed_mapping <- qtl_a3a5 %>%
  dplyr::filter(log10p >=eigen_BF_adjusted) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F)

gm <- read.table("Raw/Analysis_Results-20211021_ratio_outliers/Genotype_Matrix/Genotype_Matrix.tsv", header = T)
snp_df <- processed_mapping %>% na.omit()
ld_snps <- dplyr::filter(gm, CHROM %in% snp_df$CHROM, POS %in% snp_df$POS)

ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS,
                                       sep = "_"), data.frame(ld_snps)[, 5:ncol(ld_snps)])
  
sn <- list()
  
for (i in 1:nrow(ld_snps)) {
    sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T",
                                                    gsub(-1, "A/A", ld_snps[i, 4:ncol(ld_snps)]))))
  }
  
test <- data.frame(sn)
colnames(test) <- (ld_snps$snp_id)
ldcalc <- t(genetics::LD(test)[[4]])^2
diag(ldcalc) <- 1
  
TRAIT <- unique(processed_mapping$trait)

save(ldcalc, file="Processed_Data/LD_a3a5.RData")

ldcalc %>%
    as.data.frame() %>%
    dplyr::mutate(QTL1 = rownames(.),
                  trait = TRAIT) %>%
    tidyr::pivot_longer(cols = -c(QTL1, trait), names_to = "QTL2", values_to = "r2") %>%
    dplyr::filter(!is.na(r2)) %>%
    dplyr::select(QTL1, QTL2, everything()) %>%
    ggplot(., mapping = aes(x = factor(QTL1, levels=c("II_2429716","II_13692928","IV_6958736","X_148173"), labels = c("II:2429716","II:13692928","IV:6958736","X:148173")), 
                                       y = factor(QTL2, levels=c("II_2429716","II_13692928","IV_6958736","X_148173"), labels = c("II:2429716","II:13692928","IV:6958736","X:148173")))) + 
    theme_classic() +
    geom_tile(aes(fill = r2),colour = "black", size = 2) + 
    geom_text(aes(label = round(r2, 4))) + 
    scale_fill_gradient(low="darkgreen", high="red", limits = c(0, 1), name = expression(r^2)) + 
    theme(axis.title = element_blank(),
          axis.text = element_text(colour = "black")) + 
    labs(title = "Linkage Disequilibrium")




```
