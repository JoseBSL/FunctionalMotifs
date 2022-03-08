library(ggplot2)
library(scatterpie)
library(RColorBrewer)
library(tidyverse)


direct_indirect <- read_csv("Data/Data_processing/Motifs_connections/motif_interactions.csv")
colnames(direct_indirect)
direct_indirect <- direct_indirect %>% select(motif_id,position,direct_interactions, indirect_interactions,ind_same_type_over_directed_int)

plant_means_reordered <- read_csv("Data/Csv/plant_abs_freq_means.csv")
pollinator_means_reordered <- read_csv("Data/Csv/pollinator_abs_freq_means.csv")
#Filter out in order to have just 6 FG'S
pollinator_means_reordered <- pollinator_means_reordered %>% filter(!Node_FG %in% c("Birds", "Lizards", "Other_insects"))

#Check pol levels
levels(factor(pollinator_means_reordered$Node_FG))

#Dplyr merge
joined_tibble <- left_join(direct_indirect, pollinator_means_reordered)

#Bees
bees <- joined_tibble %>% filter(Node_FG=="Bee")
bees$ind_same_type_over_directed_int <- round(bees$ind_same_type_over_directed_int, digits = 1)
#Convert to factor
bees$position <- as.factor(bees$position)
bees$ind_same_type_over_directed_int <- as.factor(bees$ind_same_type_over_directed_int)
bee_plot <- ggplot(bees,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~ind_same_type_over_directed_int,scale='free_x', space = "free_x") + ggtitle("Bees")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Frequency")

#Coleoptera
cole <- joined_tibble %>% filter(Node_FG=="Coleoptera")
cole$ind_same_type_over_directed_int <- round(cole$ind_same_type_over_directed_int, digits = 1)
#Convert to factor
cole$position <- as.factor(cole$position)
cole$ind_same_type_over_directed_int <- as.factor(cole$ind_same_type_over_directed_int)
cole_plot <- ggplot(cole,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~ind_same_type_over_directed_int,scale='free_x', space = "free_x") + ggtitle("Coleoptera")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Frequency")

#Lepidoptera
lepi <- joined_tibble %>% filter(Node_FG=="Lepidoptera")
lepi$ind_same_type_over_directed_int <- round(lepi$ind_same_type_over_directed_int, digits = 1)
#Convert to factor
lepi$position <- as.factor(lepi$position)
lepi$ind_same_type_over_directed_int <- as.factor(lepi$ind_same_type_over_directed_int)
lepi_plot <- ggplot(lepi,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~ind_same_type_over_directed_int,scale='free_x', space = "free_x") + ggtitle("Lepidoptera") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Frequency")

#Non-bee-Hymenoptera
non_bee <- joined_tibble %>% filter(Node_FG=="Non-bee-Hymenoptera")
non_bee$ind_same_type_over_directed_int <- round(non_bee$ind_same_type_over_directed_int, digits = 2)
#Convert to factor
non_bee$position <- as.factor(non_bee$position)
non_bee$ind_same_type_over_directed_int <- as.factor(non_bee$ind_same_type_over_directed_int)
non_bee_plot <- ggplot(non_bee,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~ind_same_type_over_directed_int,scale='free_x', space = "free_x") + ggtitle("Non-bee-Hymenoptera") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Frequency")

#Non-syrphids-diptera
non_syr <- joined_tibble %>% filter(Node_FG=="Non-syrphids-diptera")
non_syr$ind_same_type_over_directed_int <- round(non_syr$ind_same_type_over_directed_int, digits = 1)
#Convert to factor
non_syr$position <- as.factor(non_syr$position)
non_syr$ind_same_type_over_directed_int <- as.factor(non_syr$ind_same_type_over_directed_int)
non_syr_plot <- ggplot(non_syr,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~ind_same_type_over_directed_int,scale='free_x', space = "free_x") + ggtitle("Non-syrphids-diptera")+
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
  ylab("Frequency")

#Syrphids
syr <- joined_tibble %>% filter(Node_FG=="Syrphids")
syr$ind_same_type_over_directed_int <- round(syr$ind_same_type_over_directed_int, digits = 1)
#Convert to factor
syr$position <- as.factor(syr$position)
syr$ind_same_type_over_directed_int <- as.factor(syr$ind_same_type_over_directed_int)
syr_plot <- ggplot(syr,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~ind_same_type_over_directed_int,scale='free_x', space = "free_x") + ggtitle("Syrphids") 
  

library(patchwork)
bee_plot / cole_plot / lepi_plot / non_bee_plot / non_syr_plot / syr_plot
