library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(sommer)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_GWA.RData")

### Heritability ###

geno_matrix <- read.table(file = "Processed_Data/Analysis_Results-20201029_ascr_fraction/Genotype_Matrix/Genotype_Matrix.tsv", header = T)
A <- sommer::A.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains_YA]))
E <- sommer::E.mat(t(geno_matrix[,colnames(geno_matrix) %in% pheno_strains_YA]))

df_h2_YA <- data.frame(feature = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_GWA_fraction_YA)-1)) {
  
  feature <- colnames(df_GWA_fraction_YA)[i+1]
  
  df_y <- df_GWA_fraction_YA %>%
    dplyr::filter(strain != "JU1516") %>%  ### same isotype with JU1581
    dplyr::filter(strain != "JU1400") %>% ### Exclude extreme outlier
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=feature) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(feature, h2, h2_SE)
  
  df_h2_YA <- rbind(df_h2_YA, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

View(df_h2_YA)

plot_frac_h2_YA_histo <- df_h2_YA %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_frac_h2_YA_histo

ggsave(file="Plots/plot_frac_YA_h2_histo.png", plot=plot_frac_h2_YA_histo, width = 7.5, height = 5, units = "in")

df_h2_YA$feature <- gsub("\\.", "#", df_h2_YA$feature)

df_ascr_YA_frac_summary <- ascr_YA_frac %>%
  dplyr::filter(strain != "JU1400") %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(median = median(ascr_fraction), cv=sd(ascr_fraction)/mean(ascr_fraction)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., df_h2_YA, by="feature")

plot_cv_h2 <- df_ascr_YA_frac_summary %>%
  ggplot(.) +
  geom_point() +
  aes(x=cv, y=as.numeric(h2)) +
  theme_bw() +
  labs(x="Coeffecient of variation", y="Narrow-sense heritability")

plot_cv_h2

cor(df_ascr_YA_frac_summary$cv, as.numeric(df_ascr_YA_frac_summary$h2), method='spearman')

save()

### pairwise-ratio ###

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
                         paste(as.character(ratio[1,i]),as.character(ratio[2,i]),sep="/"),
                         paste(as.character(ratio[2,i]),as.character(ratio[1,i]),sep="/"))
  as.ratio.long <- tidyr::gather(comb.as,trait,value,-strain)
  ratio_list_YA[[i]] <- as.ratio.long
}

all_ratio_YA <- bind_rows(ratio_list_YA)


### pair-wise ratio h2 ###

df_ratio_wide <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  tidyr::spread(trait,value)

write_tsv(df_ratio_wide, path="Processed_Data/ascr_ratio_YA_OutlierRemoval.tsv", col_names = T)

df_ratio_h2_YA <- data.frame(trait = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_ratio_wide)-1)) {
  
  trait <- colnames(df_ratio_wide)[i+1]
  
  df_y <- df_ratio_wide %>%
    dplyr::filter(strain != "JU1516") %>%  ### same isotype with JU1581
    dplyr::filter(strain != "JU1400") %>% ### Exclude extreme outlier
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=trait) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(trait, h2, h2_SE)
  
  df_ratio_h2_YA <- rbind(df_ratio_h2_YA, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

View(df_ratio_h2_YA)

plot_ratio_h2_YA_histo <- df_ratio_h2_YA %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_ratio_h2_YA_histo

ggsave(file="Plots/plot_ratio_YA_h2_histo.pdf", plot=plot_ratio_h2_YA_histo, width = 7.5, height = 5, units = "in")

### heritability heatmap ###

df_ratio_h2_wide <- df_ratio_h2_YA %>%
  dplyr::select(-h2_SE) %>%
  dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=h2)

rownames(df_ratio_h2_wide) <- df_ratio_h2_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_h2_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Narrow-sense\nheritability\n")

ggsave("Plots/plot_ratio_h2_heatmap.pdf", width = 7.5, height = 7.5, units = "in")

### example plot for high heritability to check what's happening ###

ECA36_relatedness <- A %>%
  as.data.frame %>%
  dplyr::select(relateness_ECA36=ECA36) %>%
  dplyr::add_rownames(var='strain')

df_example_ratio <- all_ratio_YA %>%
  dplyr::filter(trait == "ascr#10/icas#3") %>%
  dplyr::left_join(., ECA36_relatedness, by='strain') %>%
  na.omit()
  
df_example_ratio %>%
  ggplot(.) +
  geom_point() +
  aes(x=relateness_ECA36, y=value) +
  theme_bw()

### after removal of ECA36 and DL238 ###

### pair-wise ratio h2 ###

df_ratio_wide2 <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  dplyr::filter(!strain %in% c("ECA36", "DL238")) %>%
  tidyr::spread(trait,value)

df_ratio_h2_YA2 <- data.frame(trait = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_ratio_wide2)-1)) {
  
  trait <- colnames(df_ratio_wide2)[i+1]
  
  df_y <- df_ratio_wide2 %>%
    dplyr::filter(strain != "JU1516") %>%  ### same isotype with JU1581
    dplyr::filter(strain != "JU1400") %>% ### Exclude extreme outlier
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=trait) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(trait, h2, h2_SE)
  
  df_ratio_h2_YA2 <- rbind(df_ratio_h2_YA2, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

View(df_ratio_h2_YA2)

save(df_h2_YA, df_ratio_h2_YA, df_ratio_h2_YA2, file="Processed_Data/df_heritability.RData")

plot_ratio_h2_YA_histo2 <- df_ratio_h2_YA2 %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_ratio_h2_YA_histo2

ggsave(file="Plots/plot_ratio_YA_h2_histo_outlier_removed.pdf", plot=plot_ratio_h2_YA_histo2, width = 7.5, height = 5, units = "in")

### heritability heatmap ###

df_ratio_h2_wide2 <- df_ratio_h2_YA2 %>%
  dplyr::select(-h2_SE) %>%
  dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=h2)

rownames(df_ratio_h2_wide2) <- df_ratio_h2_wide2$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_h2_wide2, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Narrow-sense\nheritability\n")

ggsave("Plots/plot_ratio_h2_heatmap_outlier_removed.pdf", width = 7.5, height = 7.5, units = "in")


### example plot for high heritability to check what's happening ###

N2_relatedness <- A %>%
  as.data.frame %>%
  dplyr::select(N2_relatedness=N2) %>%
  dplyr::add_rownames(var='strain')

df_example_ratio2 <- all_ratio_YA %>%
  dplyr::filter(trait == "ascr#3/icas#3", !strain %in% c("DL238", "ECA36")) %>%
  dplyr::left_join(., N2_relatedness, by='strain') %>%
  na.omit()

df_example_ratio2 %>%
  ggplot(.) +
  geom_point() +
  aes(x=N2_relatedness, y=value) +
  theme_bw()


### after outlier removal by 3 standard deviation ###

### pair-wise ratio h2 ###

df_ratio_wide3 <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(mean=mean(value), sd=sd(value)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(value>=mean-(3*sd) & value <= mean+(3*sd)) %>%
  dplyr::select(-mean, -sd) %>%
  tidyr::spread(trait,value)

df_ratio_h2_YA3 <- data.frame(trait = NA, h2 = NA, h2_SE = NA)

for(i in 1:(ncol(df_ratio_wide3)-1)) {
  
  trait <- colnames(df_ratio_wide3)[i+1]
  
  df_y <- df_ratio_wide3 %>%
    dplyr::filter(strain != "JU1516") %>%  ### same isotype with JU1581
    dplyr::filter(strain != "JU1400") %>% ### Exclude extreme outlier
    dplyr::arrange(strain) %>%
    dplyr::select(strain, value=trait) %>%
    dplyr::mutate(strain = as.character(strain))
  
  h2_res <- sommer::mmer(value~1, random=~vs(strain,Gu=A), data=df_y)
  
  h2 <- as.numeric(pin(h2_res, h2 ~ (V1) / (V1+V2))[[1]][1])
  h2_SE <- pin(h2_res, h2 ~ (V1) / (V1+V2))[[2]]
  
  h2_combine = c(trait, h2, h2_SE)
  
  df_ratio_h2_YA3 <- rbind(df_ratio_h2_YA3, h2_combine) %>%
    na.omit() %>%
    dplyr::arrange(desc(h2))
}

View(df_ratio_h2_YA3)

save(df_h2_YA, df_ratio_h2_YA, df_ratio_h2_YA2, df_ratio_h2_YA3, file="Processed_Data/df_heritability.RData")

plot_ratio_h2_YA_histo3 <- df_ratio_h2_YA3 %>%
  ggplot(.) +
  geom_histogram(binwidth = 0.05) +
  aes(x=as.numeric(h2)) +
  theme_bw() +
  labs(x="Narrow-sense heritability (h2)")

plot_ratio_h2_YA_histo3

ggsave(file="Plots/plot_ratio_YA_h2_histo_outlier_removed_bySD.pdf", plot=plot_ratio_h2_YA_histo3, width = 7.5, height = 5, units = "in")

### heritability heatmap ###

df_ratio_h2_wide3 <- df_ratio_h2_YA3 %>%
  dplyr::select(-h2_SE) %>%
  dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=h2)

rownames(df_ratio_h2_wide3) <- df_ratio_h2_wide3$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_h2_wide3, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Narrow-sense\nheritability\n")

ggsave("Plots/plot_ratio_h2_heatmap_outlier_removed_bySD.pdf", width = 7.5, height = 7.5, units = "in")


### example plot for high heritability to check what's happening ###

N2_relatedness <- A %>%
  as.data.frame %>%
  dplyr::select(N2_relatedness=N2) %>%
  dplyr::add_rownames(var='strain')

df_example_ratio3 <- all_ratio_YA %>%
  dplyr::filter(trait == "ascr#3/icas#3") %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(mean=mean(value), sd=sd(value)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(value>=mean-(3*sd) & value <= mean+(3*sd)) %>%
  dplyr::select(-mean, -sd) %>%
  dplyr::left_join(., N2_relatedness, by='strain') %>%
  na.omit()

df_example_ratio3 %>%
  ggplot(.) +
  geom_point() +
  aes(x=N2_relatedness, y=value) +
  theme_bw()

### heatmap of CV (coeffecient of variation) ###

df_ratio_CV <- all_ratio_YA %>%
  dplyr::filter(strain %in% unique(df_GWA_fraction_YA$strain)) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(mean=mean(value), sd=sd(value), cv=sd/mean) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(trait,cv)

df_ratio_CV_wide <- df_ratio_CV %>%
  #dplyr::mutate(h2=as.numeric(h2)) %>%
  tidyr::separate(trait, sep='/', into=c("trait1","trait2")) %>%
  tidyr::spread(key=trait1, value=cv)

rownames(df_ratio_CV_wide) <- df_ratio_CV_wide$trait2

ggcorrplot(as.matrix(dplyr::select(df_ratio_CV_wide, -trait2)), hc.order = FALSE, tl.srt = 90, 
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
  scale_fill_gradient(limits=c(0,10), name = "Coefficient of variation")

ggsave("Plots/plot_ratio_CV_heatmap.pdf", width = 7.5, height = 7.5, units = "in")


