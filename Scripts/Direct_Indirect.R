library(ggplot2)
library(scatterpie)
library(RColorBrewer)
library(tidyverse)


direct_indirect <- read_csv("Data/Data_processing/Motifs_connections/motif_interactions.csv")
colnames(direct_indirect)
direct_indirect <- direct_indirect %>% select(motif_id,position,direct_interactions, indirect_interactions,percentage_indirect_interactions)

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
bees$percentage_indirect_interactions <- round(bees$percentage_indirect_interactions, digits = 2)
#Convert to factor
bees$position <- as.factor(bees$position)
bees$percentage_indirect_interactions <- as.factor(bees$percentage_indirect_interactions)
bee_plot <- ggplot(bees,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~percentage_indirect_interactions,scale='free_x', space = "free_x") + ggtitle("Bees")

#Coleoptera
cole <- joined_tibble %>% filter(Node_FG=="Coleoptera")
cole$percentage_indirect_interactions <- round(cole$percentage_indirect_interactions, digits = 2)
#Convert to factor
cole$position <- as.factor(cole$position)
cole$percentage_indirect_interactions <- as.factor(cole$percentage_indirect_interactions)
cole_plot <- ggplot(cole,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~percentage_indirect_interactions,scale='free_x', space = "free_x") + ggtitle("Coleoptera")

#Lepidoptera
lepi <- joined_tibble %>% filter(Node_FG=="Lepidoptera")
lepi$percentage_indirect_interactions <- round(lepi$percentage_indirect_interactions, digits = 2)
#Convert to factor
lepi$position <- as.factor(lepi$position)
lepi$percentage_indirect_interactions <- as.factor(lepi$percentage_indirect_interactions)
lepi_plot <- ggplot(lepi,aes(position,mean_natural_units))+
  geom_col(position = "dodge")+
  facet_grid(.~percentage_indirect_interactions,scale='free_x', space = "free_x") + ggtitle("Lepidoptera")


library(patchwork)
bee_plot / cole_plot / lepi_plot
