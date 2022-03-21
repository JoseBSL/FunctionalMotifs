library(tidyverse)
library(RColorBrewer)
source("Scripts/trial_motif_positions.R")

motif_number <- 1
p1 <- plot_motif_positions(motif_number) +
geom_point(aes(color = as.factor(label)),size=10,color = "grey90")+
geom_text(aes(x, y , label = label),size=3)


motif_number <- 2
p2 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "mediumaquamarine")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 3
p3 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "mediumaquamarine")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 4
p4 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "mediumaquamarine")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 5
p5 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "lightcoral")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 6
p6 <- plot_motif_positions(motif_number) +
geom_point(aes(color = as.factor(label)),size=10,color = "goldenrod2")+
geom_text(aes(x, y , label = label),size=3)


motif_number <- 7
p7 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "mediumaquamarine")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 8   
p8 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "mediumaquamarine")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 9
p9 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "lightcoral")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 10
p10 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "lightcoral")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 11
p11 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "skyblue")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 12
p12 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "goldenrod2")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 13
p13 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "lightcoral")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 14
p14 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "lightcoral")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 15
p15 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "skyblue")+
geom_text(aes(x, y , label = label),size=3)


motif_number <- 16
p16 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "goldenrod2")+
geom_text(aes(x, y , label = label),size=3)

motif_number <- 17
p17 <- plot_motif_positions(motif_number)+
geom_point(aes(color = as.factor(label)),size=10,color = "mediumaquamarine")+
geom_text(aes(x, y , label = label),size=3) 




library(patchwork)

p <- (p1 & ylab("2 species")| p2 & ylab("3 species") | p3 ) /
(p4 & ylab("4 species")| p5 | p6 | p7) / 
(p8  & ylab("5 species")| p9 | p10 | p11 | p12) /
(p13 & ylab("5 species")| p14 | p15 | p16 | p17) &
theme(axis.title.y = element_text(color="black", size=14, face="bold"))



#Create custom legend
a <- c(1,2,3,4)
b <- c(1,2,3,4)
c <- c("a","b", "c", "d")
abc <- data.frame(a,b,c)
p1 <- ggplot(abc, aes(x=a,y=b, color=c))+geom_point() + theme_bw() + 
  scale_colour_manual(name = 'Path length \n classification', 
                      values =c(a= 'goldenrod2',b='mediumaquamarine',
                                c='skyblue', d='lightcoral'), 
                      labels = c(a='Strong \n mean path length = 1.38', 
                                 b='Medium-strong \n mean path length = 1.48',
                                 c='Medium-weak \n mean path length = 1.60',
                                 d='Weak \n mean path length = 1.85'))+
theme(legend.position="bottom")

legend <- cowplot::get_legend(p1+ theme(legend.position = "bottom"))

cowplot::plot_grid(p, legend, ncol = 1, rel_heights = c(6, 1, 0.2))





