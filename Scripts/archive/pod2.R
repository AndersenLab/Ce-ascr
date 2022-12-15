library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(Hmisc)

pod2_alt <- c("ED3052","JU258","JU2838","JU3144","JU3167","JU3169","LKC34","NIC1780","NIC1783","NIC1785","NIC1788",
               "NIC1789","NIC1794","NIC1796","NIC1798","NIC1799","NIC1802","NIC251","NIC256","NIC261","NIC271",
               "NIC274","QG2838","QG2843","QG2846","QG2854","QG2855","QG2857","QG2875","QG2932","WN2001",
               "NIC1782","NIC1790","NIC1792","QG4003","QG4006","QG4008","QG4021","QG4193")

### pod-2 fine mapping ###

df_a3_fine_pod2 <- read_tsv("Processed_Data/20201104_DL2/Analysis_Results-20201104/Fine_Mappings/Data/ascr.3_snpeff_genes.tsv") %>%
  dplyr::filter(POS>=314739 & POS<=2634728, CHROM=="II")

df_a5_fine_pod2 <- read_tsv("Processed_Data/20201104_DL2/Analysis_Results-20201104/Fine_Mappings/Data/ascr.5_snpeff_genes.tsv") %>%
  dplyr::filter(POS>=314739 & POS<=2634728, CHROM=="II")

qtl_all_a3a5_acot <- rbind(df_a3_fine_pod2, df_a5_fine_pod2) %>%
  dplyr::distinct(CHROM, POS, TRAIT, VARIANT_LOG10p, GENE_NAME, VARIANT_IMPACT)

save(qtl_all_a3a5_acot, file="Processed_Data/pod2_a3a5_fine.RData")

qtl_all_a3a5_acot %>%
  ggplot(.) +
  #geom_vline(xintercept = 13692928/1e6, color='black', alpha = 0.5, linetype=2) +
  #geom_vline(xintercept = 13185987/1e6, color='cyan', alpha = 0.5, linetype=2) +
  #geom_vline(xintercept = 13187314/1e6, color='cyan', alpha = 0.5, linetype=2) +
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

### pairwise ratio association with pod-2 ###

load("Processed_Data/ratio_trait.RData")

df_ratio_split <- all_ratio_YA %>%
  tidyr::separate(trait, into=c("feature1", "feature2"), sep="_") %>%
  dplyr::mutate(pod2_gt = ifelse(strain %in% pod2_alt, 1, -1)) %>%
  dplyr::filter(feature1 %in% gsub("\\.","#", heri_features) & feature2 %in% gsub("\\.","#", heri_features))
  
df_ratio_pod2_ve <- df_ratio_split %>%  
  dplyr::group_by(feature1, feature2) %>%
  dplyr::summarise(Variance_explained=(rcorr(value, pod2_gt)$r[1,2])^2, 
                   pvalue=rcorr(value, pod2_gt)$P[1,2], mean=mean(value)) 

save(df_ratio_pod2_ve, file="Processed_Data/df_ratio_pod2_ve.RData")

df_ratio_pod2_ve %>%
  ggplot(.) +
  geom_tile() +
  aes(x=feature1, y=feature2, fill=Variance_explained*100) +
  theme_bw() +
  theme(axis.text.y = element_text(size= 10, color='black'), 
        axis.title = element_text(size= 11, color='black'), 
        axis.text.x = element_text(size= 10, angle = 90, vjust = 0.5, hjust=1, color='black'), 
        panel.grid = element_blank()) +
  scale_fill_gradient(low='white', high='blue') +
  labs(fill="Variance\nexplained (%)")
  
### pod-2 edit ###

df_pod2_edit <- read.csv(file="Raw/pod2_edit_ascr.csv") %>%
  tidyr::gather(-feature, key="strain", value="abundance")

pod2_ratio <- data.frame(combn(as.character(unique(df_pod2_edit$feature)),m=2))

pod2_ratio_list_YA<-list()

for (i in 1:ncol(pod2_ratio)) {
  as1 <- dplyr::filter(df_pod2_edit,feature==as.character(pod2_ratio[1,i]))%>%
    rename(as1 = feature,amount1 = abundance)
  as2 <- dplyr::filter(df_pod2_edit,feature==as.character(pod2_ratio[2,i]))%>%
    rename(as2 = feature,amount2 = abundance)
  comb.as <- left_join(as1,as2,by="strain")%>%
    dplyr::mutate(ratio1= amount1/amount2,
                  ratio2= amount2/amount1)%>%
    dplyr::select(strain, ratio1, ratio2)
  
  colnames(comb.as) <- c("strain",
                         paste(as.character(pod2_ratio[1,i]),as.character(pod2_ratio[2,i]),sep="_"),
                         paste(as.character(pod2_ratio[2,i]),as.character(pod2_ratio[1,i]),sep="_"))
  as.ratio.long <- tidyr::gather(comb.as,trait,value,-strain)
  pod2_ratio_list_YA[[i]] <- as.ratio.long
}

all_ratio_YA_pod2 <- bind_rows(pod2_ratio_list_YA) %>%
  tidyr::separate(strain, into=c("strain","rep"),sep="\\.") %>%
  dplyr::select(-rep)

length(unique(all_ratio_YA_pod2$trait))

all_ratio_YA_pod2_E1 <- all_ratio_YA_pod2 %>%
  dplyr::filter(strain %in% c("N2", "E1")) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(pvalue=t.test(value~strain)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait, pvalue) %>%
  dplyr::mutate(group="E1")

all_ratio_YA_pod2_E2 <- all_ratio_YA_pod2 %>%
  dplyr::filter(strain %in% c("N2", "E2"), value!="Inf") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(pvalue=t.test(value~strain)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait, pvalue) %>%
  dplyr::mutate(group="E2")

all_ratio_YA_pod2_E12 <- all_ratio_YA_pod2 %>%
  dplyr::filter(strain %in% c("N2", "E1", "E2"), value!="Inf") %>%
  dplyr::mutate(strain=ifelse(strain =="N2", "N2", "E")) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(pvalue=t.test(value~strain)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait, pvalue)

all_ratio_YA_pod2_E12

t.test(all_ratio_YA_pod2_E1$strain, all_ratio_YA_pod2_E1$value)

save(all_ratio_YA_pod2, file="Processed_Data/df_ratio_pod2_edit.RData")

### pod-2 edit boxplot ###

plot_pod2_a3a5 <- all_ratio_YA_pod2 %>%
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

plot_pod2_a3a5

plot_pod2_a5o14 <- all_ratio_YA_pod2 %>%
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

plot_pod2_a5o14

plot_pod2_a18a3 <- all_ratio_YA_pod2 %>%
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

plot_pod2_a18a3

plot_pod2_a1a3 <- all_ratio_YA_pod2 %>%
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

plot_pod2_a1a3




###

load("Processed_Data/df_GWA.RData")



df_a5a3_ratio_YA <- df_GWA_fraction_YA %>%
  dplyr::select(strain, ascr.3, ascr.5) %>%
  dplyr::mutate(a5_a3=ascr.5/ascr.3, pod2_gt = ifelse(strain %in% pod2_alt, "ALT", "REF")) 

df_a5a3_ratio_YA %>%
  ggplot(.) +
  geom_jitter(width=0.2, alpha=0.8) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  aes(x=factor(pod2_gt, levels = c("REF", "ALT")), y=a5_a3) +
  theme_bw() +
  labs(x="pod-2 genotype", y="ascr#5 : ascr#3 (YA)")

ggsave("Plots/pod2_YA.png", width=3, height=4)

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
  dplyr::mutate(a5_a3=ascr.5/ascr.3, pod2_gt = ifelse(strain %in% pod2_alt, "ALT", "REF")) 

df_a5a3_ratio_mix %>%
  ggplot(.) +
  geom_jitter(width=0.2, alpha=0.8) +
  geom_boxplot(outlier.alpha = 0, alpha=0.5) +
  aes(x=factor(pod2_gt, levels = c("REF", "ALT")), y=a5_a3) +
  theme_bw() +
  labs(x="pod-2 genotype", y="ascr#5 : ascr#3 (mixed)")  

ggsave("Plots/pod2_mixed.png", width=3, height=4)

