library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(valr)
library(clusterProfiler) ## BiocManager::install("clusterProfiler")
library(org.Ce.eg.db) ##BiocManager::install("org.Ce.eg.db")
library(biomaRt) ##BiocManager::install("biomaRt")
library(enrichplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_qtl_peaks.RData")

df_qtl_bed <- df_qtl_bonferroni_peaks %>%
  dplyr::filter(interval_size < 1e6) %>%
  dplyr::distinct(CHROM, startPOS, endPOS) %>%
  dplyr::rename(chrom=CHROM, start=startPOS, end=endPOS) %>%
  valr::bed_merge()

WS270_gene <- read.table("Processed_Data/gene.protein_coding_WS270.bed") %>%
  dplyr::select(V1,V2,V3,V4) %>%
  tidyr::separate(V4, into=c("V4","V5")) %>%
  dplyr::select(chrom=V1, start=V2, end=V3, gene=V5)

### Genes in species-wide divergent regions

QTL_genes <- valr::bed_intersect(df_qtl_bed, WS270_gene)

QTL_genes_only <- unique(QTL_genes$gene.y)

write.table(QTL_genes_only, "Processed_Data/QTL_genes_only.tsv", row.names=F, col.names=F, quote=F)

### GO enrichment

### GO enrichment tests

enrich_go <- function(x = NULL){
  
  # GO term enrichment  
  wb_ids <- dplyr::filter(x) %>% dplyr::pull(WBGeneID)
  
  mf <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  bp <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  return(list(mf, bp))
}

### Go enrichment test for the divergent regions
# set up data base to extract gene ids

df_GO <- enrich_go(data.frame(WBGeneID=QTL_genes_only))

gene_name_Pro_BP <- setReadable(df_GO[[2]], OrgDb = org.Ce.eg.db)
barplot(gene_name_Pro_BP, showCategory=100) 

gene_name_Pro_MF <- setReadable(df_GO[[1]], OrgDb = org.Ce.eg.db)
barplot(gene_name_Pro_MF, showCategory=100)

### not relavant analysis ###