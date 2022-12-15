library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(PopGenome)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

strains <- colnames(read_tsv(file="Processed_Data/WI_328_reference_list.tsv"))

### Tajima's D

tajimas_d_temp <- function (vcf_path, chromosome = "II", interval_start = 11021073, 
                            interval_end = 12008179, window_size = 300, slide_distance = 100, 
                            samples = sort(colnames(snps[, 5:ncol(snps)])), outgroup = "XZ1516", 
                            site_of_interest = 11875145) 
{
  gen <- PopGenome::readVCF(vcf_path, numcols = 10000, tid = chromosome, 
                            frompos = interval_start, topos = interval_end, samplenames = samples, 
                            approx = T)
  gen1 <- PopGenome::set.populations(gen, list(samples), diploid = TRUE)
  gen2 <- PopGenome::set.outgroup(gen1, outgroup, diploid = TRUE)
  s_gen <- PopGenome::sliding.window.transform(gen2, width = window_size, 
                                               jump = slide_distance, whole.data = FALSE)
  test <- data.frame(snps = 1:length(as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][, 
                                                                                           ])))), position = as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][, 
                                                                                                                                                                   ]))))
  slides <- list()
  for (i in 1:length(s_gen@SLIDE.POS)) {
    slides[[i]] <- data.frame(window = i, snps = s_gen@SLIDE.POS[[i]])
  }
  slides <- dplyr::rbind_all(slides) %>% dplyr::left_join(data.frame(test), 
                                                          ., by = "snps")
  s_gen_stats <- PopGenome::neutrality.stats(s_gen)
  td <- data.frame(window = 1:length(s_gen@SLIDE.POS), Td = s_gen_stats@Tajima.D) %>% 
    dplyr::left_join(slides, ., by = "window") %>% dplyr::group_by(window) %>% 
    dplyr::mutate(swindow = min(as.numeric(position)), ewindow = max(as.numeric(position))) %>% 
    dplyr::mutate(midwindow = (swindow + ewindow)/2) %>% 
    dplyr::rename(Td = pop.1) %>% dplyr::distinct(Td, window, 
                                                  .keep_all = T)
  tajimas_d_plot <- ggplot2::ggplot(td) + ggplot2::geom_point(size = 0.8, alpha = 0.8, aes(x = midwindow/1e+06, y = Td)) + ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, 
                                                       color = "black"), axis.text.y = ggplot2::element_text(size = 10, 
                                                                                                             color = "black"), axis.title.x = ggplot2::element_text(size = 12, 
                                                                                                                                                                    color = "black", vjust = -0.3), axis.title.y = ggplot2::element_text(size = 12, 
                                                                                                                                                                                                                                         color = "black"), strip.text.x = ggplot2::element_text(size = 10, 
                                                                                                                                                                                                                                                                                                color = "black"), strip.text.y = ggplot2::element_text(size = 10, color = "black"), plot.title = ggplot2::element_text(size = 10, color = "black",                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              vjust = 1), legend.position = "none") + ggplot2::labs(x = "Genomic Position (Mb)", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      y = "Tajima's D")
  if (!is.null(site_of_interest)) {
    tajimas_d_plot <- tajimas_d_plot + ggplot2::geom_vline(ggplot2::aes(xintercept = site_of_interest/1e+06), 
                                                           color = "red", alpha = 0.7, size = 2)
  }
  return(list(td, tajimas_d_plot))
}

### daf-22 Tajima's D - window size : 50 snps, mean window size in bp = xxx, slide_distance = 1 snp

tajimaD_II <- tajimas_d_temp(vcf_path = "Processed_Data/WI.HARD-FILTERED_II_12400000..12500000.vcf.gz",
                            chromosome = "II", interval_start = 12400000, interval_end = 12450000, 
                            site_of_interest = NULL, slide_distance = 10, window_size = 100, samples = strains)

tajimaD_II[[2]] +
   geom_rect(data=df_gene_del, aes(xmin =  txstart/1e6, xmax = txend/1e6, ymin = -4 , ymax = -3), fill = 'blue', color = "black", 
            alpha = 0.7)+
   geom_text(data=df_gene_del, aes(x=(txstart+txend)/2e6, y=-4.5, label=gene_name), fontface = "italic", size=2.8) +
  theme_bw() +
  coord_cartesian(xlim=c(12415000/1e6, 12435000/1e6)) +
  theme(plot.margin = margin(t=0.1, l=0.2, r=0.1, b=0.1, unit="in"),
        panel.grid = element_blank())
