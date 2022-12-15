library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(ape)
library(tibble)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/strain_geo.RData")

## tree generation

# Generate phylogeny files

#system(glue::glue("python Scripts/vcf2phylip-master/vcf2phylip.py -i Processed_Data/ascr94_complete_sites.vcf.gz --fasta --nexus --nexus-binary"))
#WI_phy <- read.phyDat("Processed_Data/ascr94_complete_sites.min4.phy", format = "interleaved")
#save(WI_phy, file = "Processed_Data/WI_phy.RData")

load("Processed_Data/WI_phy.RData")

dm <- dist.ml(WI_phy)
treeNJ <- NJ(dm)

tree_NJ_rooted <- root(treeNJ, "ECA36")

tree_pt<-ggtree(tree_NJ_rooted,
                branch.length="none")
plot_tree <- ggtree(tree_NJ_rooted, branch.length="none")
plot_tree +
  geom_text(aes(label=label), size=2, hjust=-0.1)

df_tree <- na.omit(plot_tree[[1]])

## genome-wide tree with labels for geographic origins

df_node <- plot_tree[[1]] %>%
  dplyr::select(node, label)

isotype_info_geo <- indep_strain_info_geo %>%
  dplyr::select(strain, geo) %>%
  dplyr::select(label=strain, geo)

unique(pheno_strains_YA[!pheno_strains_YA %in% unique(df_node$label)])

df_node_geo <- df_node %>% 
  dplyr::left_join(., isotype_info_geo, by = 'label') %>%
  dplyr::distinct(node, geo) 

df_node_geo$node <- as.numeric(df_node_geo$node)

tbl_tree_geo <- as.tibble(tree_NJ_rooted) %>%
  dplyr::full_join(., df_node_geo, by = 'node') %>%
  dplyr::rename(label = geo, strain = label) %>%
  dplyr::mutate(label = ifelse(is.na(label), "Z_link", as.character(label))) %>%
  dplyr::mutate(branch.length = ifelse(branch.length < 0, 0, branch.length)) ## to control negative branch length

tree_ph_geo <- as.phylo(tbl_tree_geo, use.labels = T)

plot_tree_NJ_geo <- ggtree(tree_ph_geo, branch.length="none") + 
  aes(color = label, size = label) + 
  theme(legend.position = 'bottom', 
        legend.title = element_blank(), 
        legend.text = element_text(size=9, color = "black"),
        legend.direction = "horizontal", 
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.1,"in")) +
  #theme(legend.position = 'none', legend.text = element_text(size = 10), legend.title = element_blank(), plot.margin = margin(b=-8, l=-10, r=-6, t=-10, unit = "in")) +
  scale_color_manual(values = c(geo.colours, "Z_link"="gray")) + 
  guides(col= guide_legend(nrow=3))  +
  scale_size_manual(values=c(.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.2))

plot_tree_NJ_geo

#ggsave("tree_geo.pdf", plot_tree_NJ_geo_rotate, width = 7.5, height = 5)

save(WI_phy, tree_NJ_rooted, plot_tree, df_tree, tbl_tree_geo, tree_ph_geo, file = "Processed_Data/tree.RData")
