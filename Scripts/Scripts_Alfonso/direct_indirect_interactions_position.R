
library(tidyverse)

motifs_raw_data <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

plant_interactions <- motifs_raw_data %>% group_by(motif_id,nodes,plant) %>% count() %>% 
  rename(direct_interactions = n, position = plant) %>% 
  mutate(type = "plant",
         indirect_interactions = (nodes - 1) - direct_interactions,
         percentage_indirect_interactions = indirect_interactions/(nodes-1))

pollinator_interactions <- motifs_raw_data %>% group_by(motif_id,nodes,pollinator) %>% count() %>% 
  rename(direct_interactions = n, position = pollinator) %>% 
  mutate(type = "pollinator",
         indirect_interactions = (nodes - 1) - direct_interactions,
         percentage_indirect_interactions = indirect_interactions/(nodes-1))

motif_interactions_aux <- bind_rows(pollinator_interactions,plant_interactions)

# Remove letters from positions
motif_interactions_aux$position <- as.numeric(gsub("[^0-9.-]", "", motif_interactions_aux$position))
motif_interactions <- motif_interactions_aux %>% distinct ()

# Facets per percentage of indirect interactions:
motif_interactions$percentage_indirect_interactions %>% unique() %>% sort()

write_csv(motif_interactions, "Data/Data_processing/Motifs_connections/motif_interactions.csv")
