library(tidyverse)
library(RColorBrewer)
source("Scripts/Scripts_Alfonso/plot_motif_positions.R")


motif_number <- 1
p1 <- plot_motif_positions(motif_number)

motif_number <- 2
p2 <- plot_motif_positions(motif_number)

motif_number <- 3
p3 <- plot_motif_positions(motif_number)

motif_number <- 4
p4 <- plot_motif_positions(motif_number)

motif_number <- 5
p5 <- plot_motif_positions(motif_number)

motif_number <- 6
p6 <- plot_motif_positions(motif_number)

motif_number <- 7
p7 <- plot_motif_positions(motif_number)

motif_number <- 8   
p8 <- plot_motif_positions(motif_number)

motif_number <- 9
p9 <- plot_motif_positions(motif_number)

motif_number <- 10
p10 <- plot_motif_positions(motif_number)

motif_number <- 11
p11 <- plot_motif_positions(motif_number)

motif_number <- 12
p12 <- plot_motif_positions(motif_number)

motif_number <- 13
p13 <- plot_motif_positions(motif_number)

motif_number <- 14
p14 <- plot_motif_positions(motif_number)

motif_number <- 15
p15 <- plot_motif_positions(motif_number)

motif_number <- 16
p16 <- plot_motif_positions(motif_number)

motif_number <- 17
p17 <- plot_motif_positions(motif_number)


library(patchwork)

(p1 & ylab("2 species")| p2 & ylab("3 species") | p3 )/ (p4 & ylab("4 species")| p5 | p6 | p7)/ (p8  & ylab("5 species")| p9 | p10 | p11 | p12) / (p13 & ylab("5 species")| p14 | p15 | p16 | p17)+
 plot_annotation(tag_levels = '1') & theme(axis.title.y = element_text(color="black", size=14, face="bold"))

