library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(Hmisc)

mecr1_alt <- c("ED3052","JU258","JU2838","JU3144","JU3167","JU3169","LKC34","NIC1780","NIC1783","NIC1785","NIC1788",
               "NIC1789","NIC1794","NIC1796","NIC1798","NIC1799","NIC1802","NIC251","NIC256","NIC261","NIC271",
               "NIC274","QG2838","QG2843","QG2846","QG2854","QG2855","QG2857","QG2875","QG2932","WN2001",
               "NIC1782","NIC1790","NIC1792","QG4003","QG4006","QG4008","QG4021","QG4193")

### mecr-1 fine mapping for a3a5 trait ###

df_a3a5_mecr1 <- read_tsv("Raw/Analysis_Results-20211021_ratio_outliers/Fine_Mappings/Data/a3a5_II_12422412-13915254_bcsq_genes.tsv")

df_a3a5_mecr1 %>%
  dplyr::mutate(VARIANT_IMPACT=ifelse(is.na(VARIANT_IMPACT), "Unknown", VARIANT_IMPACT)) %>%
  arrange(desc(VARIANT_IMPACT)) %>%
  ggplot(.) +
  geom_vline(xintercept = 13692928/1e6, color='black', alpha = 0.5, linetype=2) +
  geom_vline(xintercept = 13185987/1e6, color='cyan', alpha = 0.5, linetype=2) +
  geom_vline(xintercept = 13187314/1e6, color='cyan', alpha = 0.5, linetype=2) +
  geom_point(alpha=0.7, aes(x=POS/1e6, y=VARIANT_LOG10p, fill=VARIANT_IMPACT), size = 1.5, shape=21) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        legend.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        title = element_text(size=11, color='black'),
        legend.title=element_blank(),
        legend.background = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values=c("purple", "grey","grey", "grey")) +
  scale_x_continuous(expand=c(0.01,0.01))+
  labs(x="Genomic position (Mb)", y="-log10(p)")

save(df_a3a5_mecr1, file="Processed_Data/mecr1_a3a5_fine.RData")

### pairwise ratio association with mecr-1 ###

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

df_ratio_split <- all_ratio_YA %>%
  tidyr::separate(trait, into=c("feature1", "feature2"), sep="_") %>%
  dplyr::mutate(mecr1_gt = ifelse(strain %in% mecr1_alt, 1, -1)) %>%
  dplyr::filter(feature1 %in% gsub("\\.","#", heri_features) & feature2 %in% gsub("\\.","#", heri_features))

df_NA_ratio <- data.frame(feature1=unique(df_ratio_split$feature1), 
                          feature2=unique(df_ratio_split$feature1),
                          Variance_explained=NA,
                          pvalue=NA, mean=NA)

df_ratio_mecr1_ve <- df_ratio_split %>%  
  dplyr::group_by(feature1, feature2) %>%
  dplyr::summarise(Variance_explained=(rcorr(value, mecr1_gt)$r[1,2])^2, 
                   pvalue=rcorr(value, mecr1_gt)$P[1,2], mean=mean(value)) %>%
  rbind(., df_NA_ratio)

save(df_ratio_mecr1_ve, file="Processed_Data/df_ratio_mecr1_ve.RData")

df_ratio_mecr1_ve %>%
  ggplot(.) +
  geom_tile(aes(x=feature1, y=feature2, fill=Variance_explained*100)) +
  theme_bw() +
  theme(axis.text.y = element_text(size= 10, color='black'), 
        axis.title = element_text(size= 11, color='black'), 
        axis.text.x = element_text(size= 10, angle = 90, vjust = 0.5, hjust=1, color='black'), 
        panel.grid = element_blank()) +
  scale_fill_gradient(low='white', high='blue') +
  labs(fill="Variance\nexplained (%)")




df_a3_fine_mecr1 <- read_tsv("Processed_Data/20201104_DL2/Analysis_Results-20201104/Fine_Mappings/Data/ascr.3_snpeff_genes.tsv") %>%
  dplyr::filter(POS>=12839435 & POS<=13847547, CHROM=="II")

df_a5_fine_mecr1 <- read_tsv("Processed_Data/20201104_DL2/Analysis_Results-20201104/Fine_Mappings/Data/ascr.5_snpeff_genes.tsv") %>%
  dplyr::filter(POS>=12839435 & POS<=13847547, CHROM=="II")

qtl_all_a3a5_mecr <- rbind(df_a3_fine_mecr1, df_a5_fine_mecr1) %>%
  dplyr::distinct(CHROM, POS, TRAIT, VARIANT_LOG10p, GENE_NAME, VARIANT_IMPACT)

save(qtl_all_a3a5_mecr, file="Processed_Data/mecr1_a3a5_fine.RData")

qtl_all_a3a5_mecr %>%
  ggplot(.) +
  geom_vline(xintercept = 13692928/1e6, color='black', alpha = 0.5, linetype=2) +
  geom_vline(xintercept = 13185987/1e6, color='cyan', alpha = 0.5, linetype=2) +
  geom_vline(xintercept = 13187314/1e6, color='cyan', alpha = 0.5, linetype=2) +
  geom_point(alpha=0.7, aes(x=POS/1e6, y=VARIANT_LOG10p, fill=gsub("\\.", "#", TRAIT)), size = 1.5, shape=21) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        legend.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        title = element_text(size=11, color='black'),
        legend.position = c(0.1,0.9),
        legend.title=element_blank(),
        legend.background = element_blank()) +
  scale_fill_manual(values=c("orange","purple")) +
  scale_x_continuous(expand=c(0.01,0.01))+
  labs(x="Genomic position (Mb)", y="-log10(p)")


### mecr-1 edit ###

df_mecr1_edit <- read.csv(file="Raw/mecr1_edit_ascr.csv") %>%
  tidyr::gather(-feature, key="strain", value="abundance") %>%
  dplyr::filter(feature %in% heri_features)

mecr1_ratio <- data.frame(combn(as.character(unique(df_mecr1_edit$feature)),m=2))

mecr1_ratio_list_YA<-list()

for (i in 1:ncol(mecr1_ratio)) {
  as1 <- dplyr::filter(df_mecr1_edit,feature==as.character(mecr1_ratio[1,i]))%>%
    rename(as1 = feature,amount1 = abundance)
  as2 <- dplyr::filter(df_mecr1_edit,feature==as.character(mecr1_ratio[2,i]))%>%
    rename(as2 = feature,amount2 = abundance)
  comb.as <- left_join(as1,as2,by="strain")%>%
    dplyr::mutate(ratio1= amount1/amount2,
                  ratio2= amount2/amount1)%>%
    dplyr::select(strain, ratio1, ratio2)
  
  colnames(comb.as) <- c("strain",
                         paste(as.character(mecr1_ratio[1,i]),as.character(mecr1_ratio[2,i]),sep="_"),
                         paste(as.character(mecr1_ratio[2,i]),as.character(mecr1_ratio[1,i]),sep="_"))
  as.ratio.long <- tidyr::gather(comb.as,trait,value,-strain)
  mecr1_ratio_list_YA[[i]] <- as.ratio.long
}

all_ratio_YA_mecr1 <- bind_rows(mecr1_ratio_list_YA) %>%
  tidyr::separate(strain, into=c("strain","rep"),sep="\\.") %>%
  dplyr::select(-rep)

length(unique(all_ratio_YA_mecr1$trait))

all_ratio_YA_mecr1_E1 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(strain %in% c("N2", "E1")) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(pvalue=t.test(value~strain)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait, pvalue) %>%
  dplyr::mutate(group="E1")

all_ratio_YA_mecr1_E2 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(strain %in% c("N2", "E2"), value!="Inf") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(pvalue=t.test(value~strain)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait, pvalue) %>%
  dplyr::mutate(group="E2")

all_ratio_YA_mecr1_E12 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(strain %in% c("N2", "E1", "E2"), value!="Inf") %>%
  dplyr::mutate(strain=ifelse(strain =="N2", "N2", "E")) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(pvalue=t.test(value~strain)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait, pvalue, padjust=p.adjust(pvalue, "fdr"))

all_ratio_YA_mecr1_E12

t.test(all_ratio_YA_mecr1_E1$strain, all_ratio_YA_mecr1_E1$value)

mecr1_sigs <- all_ratio_YA_mecr1_E12 %>%
  dplyr::filter(padjust < 0.05) %>%
  dplyr::select(trait)

save(all_ratio_YA_mecr1, all_ratio_YA_mecr1_E12, mecr1_sigs, file="Processed_Data/df_ratio_mecr1_edit.RData")



### mecr-1 edit boxplot ###

plot_mecr1_a3a5 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(trait == "ascr#3_ascr#5") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","E1","E2","ED","NIC"), labels=c("N2(159G)","N2(159V)(1)","N2(159V)(2)","ED3052","NIC256")), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#3 : ascr#5")

plot_mecr1_a3a5

plot_mecr1_a5a18 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(trait == "ascr#5_ascr#18") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","E1","E2","ED","NIC"), labels=c("N2(159G)","N2(159V)(1)","N2(159V)(2)","ED3052","NIC256")), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#5 : ascr#18")

plot_mecr1_a5a18

plot_mecr1_o9os10 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(trait == "oscr#9_osas#10") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","E1","E2","ED","NIC"), labels=c("N2(159G)","N2(159V)(1)","N2(159V)(2)","ED3052","NIC256")), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="oscr#9 : osas#10")

plot_mecr1_o9os10

plot_mecr1_a3o9 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(trait == "ascr#3_oscr#9") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","E1","E2","ED","NIC"), labels=c("N2(159G)","N2(159V)(1)","N2(159V)(2)","ED3052","NIC256")), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#3 : oscr#9")

plot_mecr1_a3o9

plot_mecr1_a5o14 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(trait == "ascr#5_oscr#14") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","E1","E2","ED","NIC"), labels=c("N2(159G)","N2(159V)(1)","N2(159V)(2)","ED3052","NIC256")), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#5 : oscr#14")

plot_mecr1_a5o14

plot_mecr1_a18a3 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(trait == "ascr#18_ascr#3") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","E1","E2","ED","NIC"), labels=c("N2(159G)","N2(159V)(1)","N2(159V)(2)","ED3052","NIC256")), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#18 : ascr#3")

plot_mecr1_a18a3

plot_mecr1_a1a3 <- all_ratio_YA_mecr1 %>%
  dplyr::filter(trait == "ascr#1_ascr#3") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","E1","E2","ED","NIC"), labels=c("N2(159G)","N2(159V)(1)","N2(159V)(2)","ED3052","NIC256")), y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#1 : ascr#3")

plot_mecr1_a1a3




###

load("Processed_Data/df_GWA.RData")



df_a5a3_ratio_YA <- df_GWA_fraction_YA %>%
  dplyr::select(strain, ascr.3, ascr.5) %>%
  dplyr::mutate(a5_a3=ascr.5/ascr.3, mecr1_gt = ifelse(strain %in% mecr1_alt, "ALT", "REF")) 

df_a5a3_ratio_YA %>%
  ggplot(.) +
  geom_jitter(width=0.2, alpha=0.8) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  aes(x=factor(mecr1_gt, levels = c("REF", "ALT")), y=a5_a3) +
  theme_bw() +
  labs(x="mecr-1 genotype", y="ascr#5 : ascr#3 (YA)")

ggsave("Plots/mecr1_YA.png", width=3, height=4)

### mixed stage

ascr_mixed <- read.csv(file="~/Dropbox/AndersenLab/LabFolders/PastMembers/Daehan/Manuscripts/Metoblites/Data/20190807_ascr_raw_mixed.csv") %>%
  tidyr::separate(strain, sep = " ", into = c("strain", 
                                              "batch"), remove = F)  %>%
  dplyr::filter(!strain %in% c("JU238", "JU363", "JU1516", "QG537", "QG538", "CX11400", "QG1", "LSJ1")) %>%
  dplyr::mutate(strain = ifelse(strain == "JU491", "JU1491", 
                                ifelse(strain=="JU1580", "JU1793", strain))) ## fix strain name 

ascr_mixed$strain <- gsub("\\.00","", ascr_mixed$strain)
  
df_a5a3_ratio_mix <- ascr_mixed %>%
  dplyr::select(strain, ascr.3, ascr.5) %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(ascr.3=mean(ascr.3), ascr.5=mean(ascr.5)) %>%
  dplyr::mutate(a5_a3=ascr.5/ascr.3, mecr1_gt = ifelse(strain %in% mecr1_alt, "ALT", "REF")) 

df_a5a3_ratio_mix %>%
  ggplot(.) +
  geom_jitter(width=0.2, alpha=0.8) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  aes(x=factor(mecr1_gt, levels = c("REF", "ALT")), y=a5_a3) +
  theme_bw() +
  labs(x="mecr-1 genotype", y="ascr#5 : ascr#3 (mixed)")  

ggsave("Plots/mecr1_mixed.png", width=3, height=4)

