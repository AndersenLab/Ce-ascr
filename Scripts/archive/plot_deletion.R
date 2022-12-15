library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/gene_ref_flat.Rda")

gene_del <- c("F57G9.3", "Y57A10C.1", "Y57A10C.11", "Y57A10C.3", "Y57A10C.4", "Y57A10C.6", "Y57A10C.8", "Y57A10C.7", "Y57A10C.9", "Y57A10C.10")
      # F57G9.3, Y57A10C.1, Y57A10C.11, sre-27, sre-26, daf-22, Y57A10C.8, dct-12, Y57A10C.9, Y57A10C.10

df_gene_del <- gene_ref_flat %>%
  dplyr::filter(gene %in% gene_del) %>%
  dplyr::mutate(gene_name = ifelse(gene=="Y57A10C.3", "sre-27",
                                   ifelse(gene=="Y57A10C.4", "sre-26",
                                          ifelse(gene=="Y57A10C.6", "daf-22", 
                                                 ifelse(gene=="Y57A10C.7", "dct-12", gene)))))

plot_deletion <- ggplot(df_gene_del)+
  #geom_segment(aes(x = txend, xend = txstart+(.2*(txstart-txend)), y = 0, yend = 0), color = "black", arrow = arrow(length = unit(0.1, "inches")))+
  geom_rect(aes(xmin =  txstart, xmax = txend, ymin = -1 , ymax = 0), fill = 'blue', color = "black", 
            alpha = 0.7)+
  geom_rect(aes(xmin=12411041, xmax=12440052, ymin = -2.3, ymax = -1.8), fill='red', alpha = 0.7) +
  geom_text(aes(x=(txstart+txend)/2, y=1, label=gene_name), fontface = "italic", size=2.8) +
  theme_bw() +
  theme(axis.title.x.top = element_text(size=10, color = "black"), 
        axis.text.x.top = element_text(size=9, color = "black"), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  scale_x_continuous(position = "top", expand=c(0.015,0.015), breaks = c(12410000, 12415000, 12420000, 12425000,12430000,12435000,12440000)) +
  scale_y_continuous(expand=c(0.2,0.2))+
  coord_cartesian(xlim=c(12409000, 12442000)) +
  theme(plot.margin = margin(t=0.1, l=0.2, r=0.1, b=0.1, unit="in"),
        panel.grid = element_blank()) +
  labs(x= "Chr II: 12409000-12442000")

plot_deletion

#ggsave(plot_deletion, file = "Plots/plot_deletion.pdf", width = 7.5, height = 1.2)
  


