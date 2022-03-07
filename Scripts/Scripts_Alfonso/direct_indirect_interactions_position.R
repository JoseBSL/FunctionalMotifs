
library(tidyverse)

motifs_raw_data <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")


plant_indirect_int_shame_group <- motifs_raw_data %>% select(-pollinator) %>% 
  distinct () %>% 
  group_by(motif_id,nodes) %>% count() %>% 
  rename(indirect_same_type = n) %>% mutate(indirect_same_type = indirect_same_type-1)

plant_interactions_aux <- motifs_raw_data %>% group_by(motif_id,nodes,plant) %>% count() %>% 
  rename(direct_interactions = n, position = plant) %>% 
  mutate(type = "plant",
         indirect_interactions = (nodes - 1) - direct_interactions,
         percentage_indirect_interactions = indirect_interactions/(nodes-1))

plant_interactions <- plant_interactions_aux %>% left_join(plant_indirect_int_shame_group,
                                                           by = c("motif_id","nodes")) %>%
  mutate(ind_same_type_over_directed_int=indirect_same_type/direct_interactions)


pollinator_indirect_int_shame_group <- motifs_raw_data %>% select(-plant) %>% 
  distinct () %>% 
  group_by(motif_id,nodes) %>% count() %>% 
  rename(indirect_same_type = n) %>% mutate(indirect_same_type = indirect_same_type-1)

pollinator_interactions_aux <- motifs_raw_data %>% group_by(motif_id,nodes,pollinator) %>% count() %>% 
  rename(direct_interactions = n, position = pollinator) %>% 
  mutate(type = "pollinator",
         indirect_interactions = (nodes - 1) - direct_interactions,
         percentage_indirect_interactions = indirect_interactions/(nodes-1))

pollinator_interactions <- pollinator_interactions_aux %>% 
  left_join(pollinator_indirect_int_shame_group,by = c("motif_id","nodes")) %>%
  mutate(ind_same_type_over_directed_int=indirect_same_type/direct_interactions)



motif_interactions_aux <- bind_rows(pollinator_interactions,plant_interactions)

# Remove letters from positions
motif_interactions_aux$position <- as.numeric(gsub("[^0-9.-]", "", motif_interactions_aux$position))
motif_interactions <- motif_interactions_aux %>% distinct ()

# Facets per percentage of indirect interactions:
motif_interactions$ind_same_type_over_directed_int %>% unique() %>% sort()

write_csv(motif_interactions, "Data/Data_processing/Motifs_connections/motif_interactions.csv")
