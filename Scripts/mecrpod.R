library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/heri_features.RData")

df_mecrpod_edit <- read.csv(file="Raw/mecrpod.csv") %>%
  dplyr::filter(feature %in% heri_features) %>%
  dplyr::rename(abundance=value, strain=Strain)

## edit strain ratio traits ##

mecrpod_ratio <- data.frame(combn(as.character(unique(df_mecrpod_edit$feature)),m=2))

mecrpod_ratio_list_YA<-NULL

for (i in 1:ncol(mecrpod_ratio)) {
  as1 <- dplyr::filter(df_mecrpod_edit,feature==as.character(mecrpod_ratio[1,i]))%>%
    rename(as1 = feature,amount1 = abundance)
  as2 <- dplyr::filter(df_mecrpod_edit,feature==as.character(mecrpod_ratio[2,i]))%>%
    rename(as2 = feature,amount2 = abundance)
  comb.as <- left_join(as1,as2,by=c("strain","rep"))%>%
    dplyr::mutate(ratio1= amount1/amount2,
                  ratio2= amount2/amount1)%>%
    dplyr::select(strain, ratio1, ratio2)
  
  colnames(comb.as) <- c("strain",
                         paste(as.character(mecrpod_ratio[1,i]),as.character(mecrpod_ratio[2,i]),sep="_"),
                         paste(as.character(mecrpod_ratio[2,i]),as.character(mecrpod_ratio[1,i]),sep="_"))
  as.ratio.long <- tidyr::gather(comb.as,trait,value,-strain)
  mecrpod_ratio_list_YA[[i]] <- as.ratio.long
}

all_ratio_YA_mecrpod <- bind_rows(mecrpod_ratio_list_YA) %>%
  tidyr::separate(strain, into=c("strain","rep"),sep="\\.") %>%
  dplyr::select(-rep)

save(all_ratio_YA_mecrpod, file="Processed_Data/df_ratio_double_edit.RData")

plot_mecrpod_a3a5 <- all_ratio_YA_mecrpod %>%
  dplyr::filter(trait == "ascr#3_ascr#5") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=factor(strain, levels=c("N2","ECA2818","ECA3130","ECA3131","ECA3128","ECA3129"), 
               labels=c("N2","mecr-1(159V)","pod-2(1516Y)","pod-2(1516Y)","pod-2(1516Y) mecr-1(159V)","pod-2(1516Y) mecr-1(159V)")), 
      y=value,
      color=strain) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#3 : ascr#5")

plot_mecrpod_a3a5


all_ratio_YA_mecrpod %>%
  dplyr::filter(trait == "ascr#5_ascr#3") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=strain, y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#5 : ascr#3")

all_ratio_YA_mecrpod %>%
  dplyr::filter(trait == "ascr#1_ascr#3") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=strain, y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#1 : ascr#3")


all_ratio_YA_mecrpod %>%
  dplyr::filter(trait == "ascr#18_ascr#3") %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red", alpha = 0.5) +
  #geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", binwidth=0.00001) +
  aes(x=strain, y=value) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black')) +
  labs(x="", y="ascr#18 : ascr#3")


