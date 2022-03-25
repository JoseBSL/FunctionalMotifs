
library(tidyverse)
library(lme4)
library(nlme)
library(RColorBrewer)
library(stringi)
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

###########################
# ESTIMATING SPECIFICITY
###########################
# Reference: 10.1371/journal.pone.0114674

motif_pattern_connections <- read_csv("Data/Data_processing/Motifs_connections/motif_pattern_connections.csv")

number_plant <- motif_pattern_connections %>% group_by(motif_id,plant) %>% 
  count() %>% group_by(motif_id) %>% count() %>% rename(number_plants = n)

number_pollinator <- motif_pattern_connections %>% group_by(motif_id,pollinator) %>% 
  count() %>% group_by(motif_id) %>% count() %>% rename(number_pollinators = n)

specificity_pollinator_aux <- motif_pattern_connections %>% 
  group_by(motif_id,pollinator) %>% count()

specificity_pollinator_aux$position <- stri_extract_first_regex(specificity_pollinator_aux$pollinator,
                                                                "[0-9]+") %>% as.numeric()

specificity_pollinator_aux2 <- specificity_pollinator_aux %>% 
  left_join(number_plant, by = "motif_id")

specificity_pollinator <- specificity_pollinator_aux2 %>% ungroup %>%
  mutate(s=(number_plants-n)/(number_plants-1)) %>%
  select(position,s) %>% unique() %>% mutate(position=paste0("np",position))

# We replace NaN of some specialist positions by one
specificity_pollinator$s[is.nan(specificity_pollinator$s)] <- 1.0


pollinator_position_richness_specifity <- pollinator_position_richness %>%
  left_join(specificity_pollinator, by = c("position"))

# Plot positions
specificity_pollinator_panels <- unique(pollinator_position_richness_specifity[,c('s','position')])

ggplot(data = pollinator_position_richness_specifity %>% filter(! Node_FG %in% 
                                                                  c("Birds","Lizards","Other_insects")))+
  geom_point(aes(x = pollinator_richness, y = observed_freq, color = Node_FG), alpha=0.3)+
  geom_smooth(aes(x = pollinator_richness, y = observed_freq, color = Node_FG),method = "lm", formula = y~x)+
  geom_rect(data = specificity_pollinator_panels, aes(fill = as.factor(s)), xmin = -Inf, xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.2)+
  facet_wrap(~position,scales = "free")+
  scale_color_brewer(palette = "Set1")+
  labs(x="Richness",y="Observed frequencies", color="Group", fill = "Specificity")

unique(pollinator_position_richness$Node_FG )
