
library(tidyverse)
library(lme4)
library(nlme)
library(RColorBrewer)

# Estimate pollinator richness per network

int.threshold <- 1

networks <- read_csv("Data/Csv/data_for_motifs_analysis_1.csv") %>% select(-X1) %>%
  filter(Interaction >= int.threshold)

pollinator_richness <- networks %>% select(Network_id,Pollinator_species,Pollinator_functional_group) %>%
  unique() %>% group_by(Network_id,Pollinator_functional_group) %>% count() %>% 
  rename(pollinator_richness = n, Node_FG = Pollinator_functional_group)


# Extract information about the positions
file_pos <- "Data/Csv/Motifs_positions_and_null_models/GF_positions_frequency_percentile.csv"
all_position_percentiles <- read_csv(file_pos)

# Create different dataframes for plant and pollinators, and remove data for postions np1 and np2,
# and leave only those positions that can be occupied either by plants or poll. 
list_plants_FG <- as.character(1:10)

pollinator_position_percentiles_filtered <- all_position_percentiles %>% 
  filter(!Node_FG %in% list_plants_FG, !position %in% c("np1","np2"),
         !is.na(expected_freq_its_GF))

pollinator_position_richness <- pollinator_position_percentiles_filtered %>%
  left_join(pollinator_richness, by = c("Network_id","Node_FG")) %>%
  mutate(log10_ob_frec = log10(observed_freq), log10_richness = log10(pollinator_richness))


# Plot positions

ggplot(pollinator_position_richness %>% filter(! Node_FG %in% c("Birds","Lizards","Other_insects")), 
       aes(x = pollinator_richness, y = observed_freq, color = Node_FG))+
  geom_point(alpha=0.3)+
  geom_smooth(method = "lm", formula = y~x)+
  facet_wrap(~position,scales = "free")+
  scale_color_brewer(palette = "Set1")+
  labs(x="Richness",y="Observed frequencies", color="Group")+
  theme_bw()

unique(pollinator_position_richness$Node_FG )
